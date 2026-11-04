"""Bounded, deterministic calibration from strong MS1 peaks, separately by polarity.

This estimates repeatability and chromatographic timing, not absolute mass
accuracy. Noise settings are never tuned. Unsupported estimates retain the
configured fallback values and carry explicit reasons in provenance.
"""

from collections import Counter
from contextlib import closing
from dataclasses import asdict, dataclass, replace
import heapq
import time

import numpy as np
from numba import njit

from .centroid import centroid_profile
from .calibration import measure_representative, reference_measurements
from .config import Config, TUNABLE_FIELDS
from .features import TraceTracker
from .mass_tolerance import MassRegion
from .readers import MS2Scan, iter_scans


POLICY_VERSION = 2


def config_for(polarity, config, configs=None):
    if not configs:
        return config
    if polarity in configs:
        return configs[polarity]
    if polarity == "unknown" and len(configs) == 1:
        return next(iter(configs.values()))
    return config


@dataclass(frozen=True, slots=True)
class Representative:
    mz: float
    rt: float
    start: float
    end: float
    height: float
    fwhm: float
    width: float
    mass_error_da: float
    mass_error_ppm: float
    observations: int
    gaussian_similarity: float = 1.0


@dataclass
class TuningResult:
    configs: dict
    diagnostics: dict


@njit(cache=True, nogil=True)
def competing_peaks(mz, intensity, indices, ppm, da):
    flags = np.zeros(len(indices))
    for k, i in enumerate(indices):
        tolerance = max(da, mz[i] * ppm * 1e-6)
        left = np.searchsorted(mz, mz[i] - tolerance)
        right = np.searchsorted(mz, mz[i] + tolerance, side="right")
        for j in range(left, right):
            if i != j and intensity[j] >= intensity[i] * 0.2:
                flags[k] = 1.0
                break
    return flags


def representative(trace, config):
    """Measure a complete, isolated peak for the calibration candidate pool."""
    values, reason = measure_representative(
        trace,
        float(config.tuning_search_da),
        float(config.tuning_search_ppm),
        float(config.tuning_min_gaussian_similarity),
    )
    return (
        (Representative(*reference_measurements(trace)), None)
        if values is not None
        else (None, reason)
    )


def _inliers(values, floor=1e-9):
    """MAD screen with a multiplicative floor; zero spread cannot reject everything."""
    values = np.asarray(values, dtype=float)
    logs = np.log(np.maximum(values, floor))
    center = np.median(logs)
    scale = 1.4826 * np.median(np.abs(logs - center))
    return np.abs(logs - center) <= max(4.5 * scale, np.log(4.0))


def select_representatives(candidates, maximum):
    """Prefer strong peaks within m/z/RT strata, avoiding a single dominant region."""
    if not candidates:
        return []
    masses = np.array([p.mz for p in candidates])
    times = np.array([p.rt for p in candidates])
    mz_edges = np.quantile(masses, [0.25, 0.5, 0.75])
    rt_edges = np.quantile(times, [0.25, 0.5, 0.75])
    cells = {}
    for peak in sorted(candidates, key=lambda p: (-p.height, p.rt, p.mz)):
        key = (int(np.searchsorted(mz_edges, peak.mz)), int(np.searchsorted(rt_edges, peak.rt)))
        cells.setdefault(key, []).append(peak)
    selected = []
    depth = 0
    while len(selected) < maximum:
        added = False
        for key in sorted(cells):
            if depth < len(cells[key]):
                selected.append(cells[key][depth])
                added = True
                if len(selected) == maximum:
                    break
        if not added:
            break
        depth += 1
    return selected


def _screen_representatives(peaks, config):
    selected = select_representatives(peaks, config.tuning_max_peaks)
    before = len(selected)
    if selected:
        keep = _inliers([p.fwhm for p in selected]) & _inliers(
            [p.mass_error_ppm for p in selected], floor=0.5
        )
        selected = [p for p, use in zip(selected, keep) if use]
    return selected, before


def _mass_estimate(errors, config, minimum=(0.0005, 5.0)):
    """One Q95 absolute residual per peak; mean + k * sample SD across peaks."""
    values = np.array([[error, error / mass * 1e6] for error, mass in errors])
    mean, sd = values.mean(axis=0), values.std(axis=0, ddof=1)
    estimate = np.minimum(
        np.maximum(minimum, mean + config.tuning_mass_sigma * sd),
        [config.tuning_search_da, config.tuning_search_ppm],
    )
    return tuple(float(v) for v in estimate), {
        "count": len(errors),
        "mean_da_ppm": mean.tolist(),
        "sample_sd_da_ppm": sd.tolist(),
        "sigma_multiplier": config.tuning_mass_sigma,
    }


def _precursor_errors(peaks, scans, config, polarity):
    """Match against the union of selected peaks, so regional edges hide no competitor."""
    residuals = [[] for _ in peaks]
    for scan in scans:
        if (
            scan.polarity != polarity
            or scan.precursor_count != 1
            or scan.precursor_mz is None
            or scan.rt is None
            or scan.precursor_mz_source != "selected_ion"
        ):
            continue
        tolerance = max(
            config.tuning_search_da, scan.precursor_mz * config.tuning_search_ppm * 1e-6
        )
        hits = [
            i
            for i, p in enumerate(peaks)
            if p.start <= scan.rt <= p.end and abs(p.mz - scan.precursor_mz) <= tolerance
        ]
        if len(hits) == 1:
            index = hits[0]
            residuals[index].append(abs(scan.precursor_mz - peaks[index].mz))
    return {p: float(np.quantile(v, 0.95)) for p, v in zip(peaks, residuals) if v}


def _screen_precursor_errors(peaks, errors):
    values = [(errors[p], p.mz) for p in peaks if p in errors]
    if not values:
        return []
    keep = _inliers([e / m * 1e6 for e, m in values], floor=0.5)
    return [value for value, use in zip(values, keep) if use]


def derive_settings(peaks, intervals, ms2, config, polarity, *, mass_range=None):
    """Independent MS1 and precursor estimates, each requiring >=20 representatives."""
    defaults = Config()
    locked = set(config.tuning_locked) | {
        key for key in TUNABLE_FIELDS if getattr(config, key) != getattr(defaults, key)
    }
    # Legacy scan-unit overrides also protect the corresponding physical setting.
    if config.smoothing_sigma != defaults.smoothing_sigma:
        locked.add("smoothing_seconds")
    if config.peak_distance != defaults.peak_distance:
        locked.add("peak_distance_seconds")
    info = {"status": "insufficient_representatives", "protected_parameters": sorted(locked)}
    updates = {}
    intervals = np.asarray(intervals, dtype=float)
    intervals = intervals[np.isfinite(intervals) & (intervals > 0)]
    if len(intervals):
        intervals = intervals[_inliers(intervals)]
    dt = float(np.median(intervals)) if len(intervals) else None
    info["ms1_interval_seconds"] = dt
    if dt is not None:
        # With insufficient peaks, only translate legacy scan units to seconds.
        updates.update(
            smoothing_seconds=config.smoothing_sigma * dt,
            peak_distance_seconds=config.peak_distance * dt,
        )
    qualified = [p for p in peaks if p.gaussian_similarity >= config.tuning_min_gaussian_similarity]
    selected, before = _screen_representatives(qualified, config)
    info["representatives_before_outlier_screen"] = before
    info["representative_count"] = len(selected)
    info["representatives"] = [asdict(p) for p in selected]
    info["ms2_status"] = "insufficient_representatives"
    if len(selected) >= config.tuning_min_peaks and dt is not None and polarity != "unknown":
        info["status"] = "tuned"
        # Each representative gets one vote, regardless of intensity/scan count.
        (da, ppm), info["ms1_error_statistics"] = _mass_estimate(
            [(p.mass_error_da, p.mz) for p in selected], config
        )
        widths = np.array([p.width for p in selected])
        fwhm = np.array([p.fwhm for p in selected])
        typical = float(np.median(fwhm))
        fast = float(np.quantile(fwhm, 0.1))
        cycle_high = float(np.quantile(intervals, 0.95))
        updates.update(
            mz_tolerance_da=da,
            mz_tolerance_ppm=ppm,
            gapfill_mz_tolerance_da=da,
            gapfill_mz_tolerance_ppm=ppm,
            smoothing_seconds=max(0.5 * dt, 0.08 * fast),
            peak_distance_seconds=max(dt, 0.5 * fast),
            max_gap_seconds=cycle_high * (config.max_gap_scans + 1) * 1.25,
            precursor_max_age_seconds=cycle_high * 2,
            gapfill_max_gap_seconds=cycle_high * (config.max_gap_scans + 1) * 1.25,
            ms2_rt_tolerance=max(cycle_high, 0.15 * typical),
            ms2_gapfill_rt_window=min(
                config.tuning_block_seconds, max(3 * cycle_high, float(np.quantile(widths, 0.9)))
            ),
            ms2_seed_rt_tolerance=max(dt, 0.5 * fast),
            group_rt_tolerance=max(dt, 0.25 * typical),
        )
        info["fwhm_seconds_p10_p50_p90"] = np.quantile(fwhm, [0.1, 0.5, 0.9]).tolist()
        info["mass_range"] = [min(p.mz for p in selected), max(p.mz for p in selected)]
        info["rt_range"] = [min(p.rt for p in selected), max(p.rt for p in selected)]
    elif polarity == "unknown":
        info["status"] = "unknown_polarity"
    regions = _derive_mass_regions(
        qualified, selected, ms2, config, polarity, mass_range, updates, info, locked
    )
    info["estimated_parameters"] = dict(updates)
    applied = {key: value for key, value in updates.items() if key not in locked}
    info["applied_parameters"] = applied
    return replace(config, **applied, mass_regions=regions), info


def _derive_mass_regions(peaks, selected, ms2, config, polarity, mass_range, updates, info, locked):
    """Regional fits have independent evidence; sparse regions use a file-wide envelope."""
    edges = []
    if mass_range is not None:
        low, high = mass_range
        if np.isfinite(low) and np.isfinite(high) and 0 <= low < high:
            edges = np.linspace(low, high, config.tuning_mz_regions + 1).tolist()
    pools, records = [], []
    for index, (low, high) in enumerate(zip(edges, edges[1:])):
        candidates = [
            p
            for p in peaks
            if low <= p.mz
            and (p.mz < high or index == config.tuning_mz_regions - 1 and p.mz == high)
        ]
        regional, before = _screen_representatives(candidates, config)
        pools.append(regional)
        records.append(
            {
                "lower": low,
                "upper": high,
                "candidate_count": len(candidates),
                "representatives_before_outlier_screen": before,
                "representative_count": len(regional),
            }
        )
    # Distinct peaks only. A peak chosen both globally and locally still gets one vote.
    union = sorted(set(selected).union(*(set(p) for p in pools)), key=lambda p: (p.mz, p.rt))
    errors = _precursor_errors(union, ms2, config, polarity) if info["status"] == "tuned" else {}
    global_errors = _screen_precursor_errors(selected, errors)
    info["ms2_representative_count"] = len(global_errors)
    ms1 = (
        updates.get("mz_tolerance_da", config.mz_tolerance_da),
        updates.get("mz_tolerance_ppm", config.mz_tolerance_ppm),
    )
    if len(global_errors) >= config.tuning_min_peaks:
        pair, info["ms2_error_statistics"] = _mass_estimate(global_errors, config, ms1)
        updates["ms2_mz_tolerance_da"], updates["ms2_mz_tolerance_ppm"] = pair
        info["ms2_status"] = "tuned"
    estimates = []
    for regional, record in zip(pools, records):
        values = {}
        eligible = info["status"] == "tuned" and len(regional) >= config.tuning_min_peaks
        record["ms1_status"] = "tuned" if eligible else "file_wide_fallback"
        if eligible:
            pair, record["ms1_error_statistics"] = _mass_estimate(
                [(p.mass_error_da, p.mz) for p in regional], config
            )
            values.update(
                mz_tolerance_da=pair[0],
                mz_tolerance_ppm=pair[1],
                gapfill_mz_tolerance_da=pair[0],
                gapfill_mz_tolerance_ppm=pair[1],
            )
        regional_errors = _screen_precursor_errors(regional, errors)
        record["ms2_representative_count"] = len(regional_errors)
        record["ms2_status"] = "file_wide_fallback"
        if eligible and len(regional_errors) >= config.tuning_min_peaks:
            pair, record["ms2_error_statistics"] = _mass_estimate(
                regional_errors, config, (values["mz_tolerance_da"], values["mz_tolerance_ppm"])
            )
            values.update(ms2_mz_tolerance_da=pair[0], ms2_mz_tolerance_ppm=pair[1])
            record["ms2_status"] = "tuned"
        estimates.append(values)
    # Unsupported regions must not inherit a narrower tolerance than any supported
    # region. If the global fit failed, also retain the configured fallback.
    mass_fields = (
        "mz_tolerance_da",
        "mz_tolerance_ppm",
        "ms2_mz_tolerance_da",
        "ms2_mz_tolerance_ppm",
        "gapfill_mz_tolerance_da",
        "gapfill_mz_tolerance_ppm",
    )
    for key in mass_fields:
        regional_values = [values[key] for values in estimates if key in values]
        if regional_values:
            updates[key] = max(updates.get(key, getattr(config, key)), *regional_values)
    fallback = {key: updates.get(key, getattr(config, key)) for key in mass_fields}
    models = []
    indices = {p: index for index, p in enumerate(union)}
    for regional, record, estimate in zip(pools, records, estimates):
        effective = {**fallback, **estimate}
        # Sparse MS2 evidence cannot reduce the local MS1-derived minimum.
        if record["ms2_status"] != "tuned":
            for unit in ("da", "ppm"):
                effective["ms2_mz_tolerance_" + unit] = max(
                    effective["ms2_mz_tolerance_" + unit], effective["mz_tolerance_" + unit]
                )
        record["estimated_parameters"] = dict(effective)
        for key in locked & set(mass_fields):
            effective[key] = getattr(config, key)
        models.append(MassRegion(record["lower"], record["upper"], **effective))
        record["representative_indices"] = [indices[p] for p in regional]
        record["effective_parameters"] = effective
    info["regional_representatives"] = [asdict(p) for p in union]
    info["regions"] = records
    info["regional_status"] = (
        "configured_regions"
        if config.mass_regions
        else "enabled"
        if records
        else "missing_mass_range"
    )
    return config.mass_regions or tuple(models)


def calibrate_file(source, config, *, scan_source=None):
    if not config.auto_tune:
        return TuningResult({}, {"enabled": False, "policy_version": POLICY_VERSION})
    started = time.perf_counter()
    pilot_config = replace(
        config,
        min_scans=7,
        max_gap_scans=1,
        group_ions=False,
        mz_tolerance_da=config.tuning_search_da,
        mz_tolerance_ppm=config.tuning_search_ppm,
        mass_regions=(),
    )
    states = {}
    ms2 = []
    previous = {}
    sequence = 0

    def consume(state, traces):
        nonlocal sequence
        for trace in traces:
            peak, reason = representative(trace, config)
            if peak is None:
                state["rejected"][reason] += 1
                continue
            state["candidates"] += 1
            sequence += 1
            item = (peak.height, sequence, peak)
            if len(state["heap"]) < config.tuning_max_peaks * config.tuning_mz_regions * 10:
                heapq.heappush(state["heap"], item)
            elif item[:2] > state["heap"][0][:2]:
                heapq.heapreplace(state["heap"], item)

    with closing(
        scan_source.scans()
        if scan_source is not None
        else iter_scans(
            source, assume_centroid=config.assume_centroid, include_ms2=config.associate_ms2
        )
    ) as scans:
        for scan in scans:
            if isinstance(scan, MS2Scan):
                # Calibration needs precursor metadata only. The bounded source
                # keeps its own replay copy; never pin full fragment arrays in
                # this file-wide list after the source cache overflows.
                scan.mz = scan.intensity = None
                if scan.polarity == "unknown" and scan.precursor_ref:
                    parents = [s for s in previous.values() if s.native_id == scan.precursor_ref]
                    if len(parents) == 1:
                        scan.polarity = parents[0].polarity
                ms2.append(scan)
                continue
            previous[scan.polarity] = scan
            if scan.polarity not in states:
                states[scan.polarity] = {
                    "tracker": TraceTracker(pilot_config),
                    "start": scan.rt,
                    "block_scans": 0,
                    "last_rt": None,
                    "intervals": [],
                    "interval_count": 0,
                    "rng": np.random.default_rng(
                        {"positive": 1, "negative": 2, "unknown": 3}[scan.polarity]
                    ),
                    "heap": [],
                    "rejected": Counter(),
                    "candidates": 0,
                    "scans": 0,
                    "observed_range": [np.inf, -np.inf],
                    "acquisition_range": [np.inf, -np.inf],
                    "range_metadata_scans": 0,
                    "invalid_range_scans": 0,
                }
            state = states[scan.polarity]
            if len(scan.mz):
                state["observed_range"][0] = min(state["observed_range"][0], float(scan.mz[0]))
                state["observed_range"][1] = max(state["observed_range"][1], float(scan.mz[-1]))
            if scan.acquisition_range is not None:
                state["range_metadata_scans"] += 1
                state["acquisition_range"][0] = min(
                    state["acquisition_range"][0], scan.acquisition_range[0]
                )
                state["acquisition_range"][1] = max(
                    state["acquisition_range"][1], scan.acquisition_range[1]
                )
            state["invalid_range_scans"] += int(scan.invalid_acquisition_range)
            if state["last_rt"] is not None:
                delta = scan.rt - state["last_rt"]
                if delta <= 0:
                    raise ValueError(
                        "MS1 retention times must increase strictly within each polarity"
                    )
                state["interval_count"] += 1
                if len(state["intervals"]) < 4096:
                    state["intervals"].append(delta)
                else:
                    position = int(state["rng"].integers(state["interval_count"]))
                    if position < 4096:
                        state["intervals"][position] = delta
            state["last_rt"] = scan.rt
            state["scans"] += 1
            tracker = state["tracker"]
            if (
                scan.rt - state["start"] >= config.tuning_block_seconds
                or state["block_scans"] >= config.tuning_block_scans
            ):
                consume(state, tracker.flush())
                state["start"], state["block_scans"] = scan.rt, 0
            state["block_scans"] += 1
            mz, y = scan.mz, scan.intensity
            if not scan.centroided:
                if config.profile_mode == "reject":
                    raise ValueError("Profile MS1 spectrum; profile_mode is 'reject'")
                zeros = np.zeros(len(mz))
                mz, y, _ = centroid_profile(
                    mz, y, zeros, zeros, config.profile_min_points, float(config.profile_gap_factor)
                )
            if len(y) > config.tuning_peak_limit:
                cutoff = np.partition(y, len(y) - config.tuning_peak_limit)[
                    len(y) - config.tuning_peak_limit
                ]
                high = np.flatnonzero(y > cutoff)
                tied = np.flatnonzero(y == cutoff)[: config.tuning_peak_limit - len(high)]
                indices = np.sort(np.r_[high, tied])
            else:
                indices = np.arange(len(y))
            flags = competing_peaks(
                mz, y, indices, config.tuning_search_ppm, config.tuning_search_da
            )
            consume(state, tracker.add(mz[indices], y[indices], flags, scan.rt))
    configs, diagnostics = {}, {}
    for polarity, state in states.items():
        consume(state, state["tracker"].flush())
        observed, acquired = state["observed_range"], state["acquisition_range"]
        bounds = (min(observed[0], acquired[0]), max(observed[1], acquired[1]))
        range_source = (
            "acquisition_metadata" if state["range_metadata_scans"] else "observed_ms1_extent"
        )
        if state["range_metadata_scans"] and (bounds[0] < acquired[0] or bounds[1] > acquired[1]):
            range_source = "acquisition_metadata_with_observed_extension"
        configs[polarity], info = derive_settings(
            [item[2] for item in state["heap"]],
            state["intervals"],
            ms2,
            config,
            polarity,
            mass_range=bounds,
        )
        info.update(
            ms1_scans=state["scans"],
            candidates=state["candidates"],
            rejected_candidates=dict(sorted(state["rejected"].items())),
            acquired_mz_range=list(bounds) if all(np.isfinite(bounds)) else None,
            acquired_mz_range_source=range_source,
            range_metadata_scans=state["range_metadata_scans"],
            invalid_range_metadata_scans=state["invalid_range_scans"],
            retained_candidates=len(state["heap"]),
        )
        diagnostics[polarity] = info
    return TuningResult(
        configs,
        {
            "enabled": True,
            "policy_version": POLICY_VERSION,
            "elapsed_seconds": time.perf_counter() - started,
            "minimum_representatives": config.tuning_min_peaks,
            "mass_regions": config.tuning_mz_regions,
            "minimum_gaussian_similarity": config.tuning_min_gaussian_similarity,
            "mass_sigma": config.tuning_mass_sigma,
            "search_da": config.tuning_search_da,
            "search_ppm": config.tuning_search_ppm,
            "polarities": diagnostics,
            "effective_configs": {p: c.to_dict() for p, c in configs.items()},
        },
    )
