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
from scipy.ndimage import gaussian_filter1d
from scipy.signal import find_peaks

from .centroid import centroid_profile
from .config import Config, TUNABLE_FIELDS
from .features import TraceTracker
from .readers import MS2Scan, iter_scans


POLICY_VERSION = 1


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


def _crossings(t, y, apex, level):
    left = right = apex
    while left > 0 and y[left] > level:
        left -= 1
    while right < len(y) - 1 and y[right] > level:
        right += 1
    if y[left] > level or y[right] > level or left == apex or right == apex:
        return None
    a = t[left] + (t[left + 1] - t[left]) * (level - y[left]) / (y[left + 1] - y[left])
    b = t[right - 1] + (t[right] - t[right - 1]) * (y[right - 1] - level) / (
        y[right - 1] - y[right]
    )
    return float(a), float(b)


def representative(trace, config):
    """Require a resolved, complete rise/fall and a stable central mass trace."""
    if len(trace) < 7:
        return None, "short_trace"
    mz, y, t, competing = trace.T
    if not np.all(np.isfinite(trace)) or np.max(y) <= 0:
        return None, "invalid_trace"
    smooth = gaussian_filter1d(y, 0.6)
    baseline = float(np.quantile(y, 0.1))
    apex = int(np.argmax(smooth))
    amplitude = float(smooth[apex] - baseline)
    if amplitude <= 0.5 * float(y.max()):
        return None, "flat_or_low_contrast"
    peaks, _ = find_peaks(smooth, prominence=amplitude * 0.3)
    if len(peaks) != 1 or peaks[0] != apex:
        return None, "multiple_or_edge_peaks"
    if np.count_nonzero(y == y.max()) >= max(3, len(y) // 3):
        return None, "plateau"
    central = y >= baseline + amplitude * 0.5
    if np.count_nonzero(central) < 5:
        return None, "undersampled_peak"
    if np.any(competing[central] > 0):
        return None, "competing_mass"
    half = _crossings(t, smooth, apex, baseline + amplitude * 0.5)
    outer = _crossings(t, smooth, apex, baseline + amplitude * 0.1)
    if half is None or outer is None:
        return None, "truncated_peak"
    if np.sum(np.abs(np.diff(smooth))) > 4 * amplitude:
        return None, "rough_trace"
    center = float(np.median(mz[central]))
    error = float(np.quantile(np.abs(mz[central] - center), 0.95))
    search = max(config.tuning_search_da, center * config.tuning_search_ppm * 1e-6)
    if error > search * 0.4:
        return None, "mass_near_search_boundary"
    return Representative(
        center,
        float(t[apex]),
        half[0],
        half[1],
        float(y.max()),
        half[1] - half[0],
        outer[1] - outer[0],
        error,
        error / center * 1e6,
        int(np.count_nonzero(central)),
    ), None


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


def derive_settings(peaks, intervals, ms2, config, polarity):
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
    selected = select_representatives(peaks, config.tuning_max_peaks)
    info["representatives_before_outlier_screen"] = len(selected)
    if selected:
        keep = _inliers([p.fwhm for p in selected]) & _inliers(
            [p.mass_error_ppm for p in selected], floor=0.5
        )
        selected = [p for p, use in zip(selected, keep) if use]
    info["representative_count"] = len(selected)
    info["representatives"] = [asdict(p) for p in selected]
    info["ms2_status"] = "insufficient_representatives"
    if len(selected) >= config.tuning_min_peaks and dt is not None and polarity != "unknown":
        info["status"] = "tuned"
        # Each representative gets one vote, regardless of intensity/scan count.
        da = max(0.0005, 2 * float(np.quantile([p.mass_error_da for p in selected], 0.95)))
        ppm = max(5.0, 2 * float(np.quantile([p.mass_error_ppm for p in selected], 0.95)))
        da, ppm = min(da, config.tuning_search_da), min(ppm, config.tuning_search_ppm)
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
        residuals = [[] for _ in selected]
        for scan in ms2:
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
                for i, peak in enumerate(selected)
                if peak.start <= scan.rt <= peak.end
                and abs(peak.mz - scan.precursor_mz) <= tolerance
            ]
            if len(hits) == 1:
                i = hits[0]
                residuals[i].append(abs(scan.precursor_mz - selected[i].mz))
        errors = [
            (float(np.quantile(v, 0.95)), selected[i].mz) for i, v in enumerate(residuals) if v
        ]
        if errors:
            keep = _inliers([error / mass * 1e6 for error, mass in errors], floor=0.5)
            errors = [e for e, use in zip(errors, keep) if use]
        info["ms2_representative_count"] = len(errors)
        if len(errors) >= config.tuning_min_peaks:
            updates["ms2_mz_tolerance_da"] = min(
                config.tuning_search_da,
                max(da, 1.5 * float(np.quantile([e for e, _ in errors], 0.95))),
            )
            updates["ms2_mz_tolerance_ppm"] = min(
                config.tuning_search_ppm,
                max(ppm, 1.5 * float(np.quantile([e / m * 1e6 for e, m in errors], 0.95))),
            )
            info["ms2_status"] = "tuned"
    elif polarity == "unknown":
        info["status"] = "unknown_polarity"
    info["estimated_parameters"] = dict(updates)
    applied = {key: value for key, value in updates.items() if key not in locked}
    info["applied_parameters"] = applied
    return replace(config, **applied), info


def calibrate_file(source, config):
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
            if len(state["heap"]) < config.tuning_max_peaks * 10:
                heapq.heappush(state["heap"], item)
            elif item[:2] > state["heap"][0][:2]:
                heapq.heapreplace(state["heap"], item)

    with closing(
        iter_scans(source, assume_centroid=config.assume_centroid, include_ms2=config.associate_ms2)
    ) as scans:
        for scan in scans:
            if isinstance(scan, MS2Scan):
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
                }
            state = states[scan.polarity]
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
        configs[polarity], info = derive_settings(
            [item[2] for item in state["heap"]], state["intervals"], ms2, config, polarity
        )
        info.update(
            ms1_scans=state["scans"],
            candidates=state["candidates"],
            rejected_candidates=dict(sorted(state["rejected"].items())),
        )
        diagnostics[polarity] = info
    return TuningResult(
        configs,
        {
            "enabled": True,
            "policy_version": POLICY_VERSION,
            "elapsed_seconds": time.perf_counter() - started,
            "minimum_representatives": config.tuning_min_peaks,
            "search_da": config.tuning_search_da,
            "search_ppm": config.tuning_search_ppm,
            "polarities": diagnostics,
            "effective_configs": {p: c.to_dict() for p, c in configs.items()},
        },
    )
