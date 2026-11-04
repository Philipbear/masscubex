"""MS2-seeded features quantified from a second, streaming pass over raw MS1.

These are forced local EIC measurements, not new noise-qualified detections.
Only unmatched usable precursors are queried. A bounded decoded-scan cache can
provide this pass; larger files are streamed again from disk.
"""

from bisect import bisect_left, bisect_right, insort
from collections import Counter
from contextlib import closing
from dataclasses import dataclass, field
import heapq

import numpy as np
from numba import njit

from .centroid import centroid_profile
from .features import Feature
from .gapfill_refinement import _radius, refine_measurement, merge_refined
from .ms2 import MS2Assignment, apply_ms2_assignments, ineligible_reason
from .readers import iter_ms1
from .tuning import config_for
from .noise_regions import regional_thresholds


@njit(cache=True, nogil=True)
def eic_intensities(mz, intensity, targets, radii):
    """Take the largest raw centroid/profile sample inside each target interval."""
    result = np.zeros(len(targets))
    for i in range(len(targets)):
        tolerance = radii[i]
        left = np.searchsorted(mz, targets[i] - tolerance)
        right = np.searchsorted(mz, targets[i] + tolerance, side="right")
        for j in range(left, right):
            result[i] = max(result[i], intensity[j])
    return result


@dataclass(slots=True)
class _Measurement:
    feature_index: int
    polarity: str | None
    scans: int = 0
    previous_rt: float | None = None
    previous_intensity: float = 0.0
    area: float = 0.0
    trace: list = field(default_factory=list)
    refine: bool = False
    samples: list = field(default_factory=list)
    sampled_values: int = 0
    sampling_limit: bool = False
    retained_values: int = 0
    finalized: bool = False
    noise_references: list = field(default_factory=list)


def _impute_missing_height(feature, measurement, config):
    """Use a noise-derived placeholder only when MS1 is absent, never as an area."""
    if feature.ms1_quant_status not in {"no_ms1_signal", "no_local_ms1_signal"}:
        return
    feature.gapfill_status = feature.ms1_quant_status
    feature.peak_height = feature.peak_area = feature.intensity = None
    # Up to three closest compatible MS1 scans in the query window. Unsupported
    # regions (threshold zero) provide no basis for a positive imputation.
    thresholds = [row[2] for row in measurement.noise_references if row[2] > 0]
    if not thresholds:
        return
    threshold = float(np.median(thresholds))
    feature.apex_noise_threshold = threshold
    if config.gapfill_noise_fraction == 0:
        return
    feature.peak_height = config.gapfill_noise_fraction * threshold
    feature.intensity = feature.peak_height if config.quant_method == "peak_height" else None
    feature.ms1_quant_status = "imputed_noise"
    # No measured peak, area, trace or scan count is manufactured. Keep the
    # gapfill_status reason and MS2 seed metadata so the placeholder is auditable.


def seed_ms2_features(
    features, records, assignments, source, config, stats, configs=None, *, scan_source=None
):
    """Keep existing assignments and measure unmatched precursors on raw MS1.

    Each seed is anchored to its first MS2 m/z and time. Later scans can reuse it
    only inside that fixed time window and the normal MS2 mass tolerance; the
    window never expands through a chain of scans. No polarity is guessed in a
    mixed-polarity run. IDs are feature-list positions plus one, as in the TSV.
    """
    stats.update(
        ms2_matched_detected=stats["ms2_assigned"],
        ms2_initially_unassigned=stats["ms2_unassigned"],
        ms2_initial_unassigned_reasons=stats["ms2_unassigned_reasons"].copy(),
        ms2_seeded_scans=0,
        ms2_seeded_features=0,
        ms2_gapfill_passes=0,
        ms2_gapfill_sample_points=0,
        ms2_gapfill_peak_retained_values=0,
        ms2_gapfill_status_counts={},
    )
    if not config.ms2_seeded_features:
        return assignments
    polarities = stats["ms1_polarities"]
    ms2_polarities = {s.polarity for s in records} - {"unknown"}
    measurements = []
    seed_index = []  # Sorted (m/z, feature index); only added features are indexed.
    first_feature = len(features)
    for i in sorted(range(len(records)), key=lambda i: records[i].spectrum_index):
        scan = records[i]
        if assignments[i].feature_index is not None or ineligible_reason(scan, config):
            continue
        polarity = scan.polarity
        if polarity == "unknown" and len(polarities) == 1:
            polarity = polarities[0]
        effective = config_for(polarity, config, configs)
        refined = (
            stats.get("parameter_tuning", {}).get("polarities", {}).get(polarity, {}).get("status")
            == "tuned"
        )
        extraction_polarity = polarity
        if polarity == "unknown" and len(polarities) > 1:
            extraction_polarity = None
        elif polarities == ["unknown"]:
            extraction_polarity = "unknown" if len(ms2_polarities) <= 1 else None
        tolerance = max(
            effective.ms2_mz_tolerance_da, scan.precursor_mz * effective.ms2_mz_tolerance_ppm * 1e-6
        )
        left = bisect_left(seed_index, (scan.precursor_mz - tolerance, -1))
        right = bisect_right(seed_index, (scan.precursor_mz + tolerance, len(features)))
        candidates = [
            j
            for _, j in seed_index[left:right]
            if not refined
            and features[j].polarity == polarity
            and abs(features[j].seed_rt_seconds - scan.rt) <= effective.ms2_seed_rt_tolerance
        ]
        if candidates:
            best = min(
                candidates,
                key=lambda j: (
                    abs(features[j].mz - scan.precursor_mz),
                    abs(features[j].rt_apex_seconds - scan.rt),
                    j,
                ),
            )
        else:
            best = len(features)
            features.append(
                Feature(
                    mz=scan.precursor_mz,
                    rt_apex_seconds=scan.rt,
                    rt_start_seconds=max(0.0, scan.rt - effective.ms2_gapfill_rt_window),
                    rt_end_seconds=scan.rt + effective.ms2_gapfill_rt_window,
                    peak_height=None,
                    peak_area=None,
                    scan_count=0,
                    apex_noise_threshold=None,
                    apex_threshold_ratio=None,
                    polarity=polarity,
                    intensity=None,
                    quant_method=config.quant_method,
                    feature_origin="ms2_seeded",
                    seed_precursor_mz=scan.precursor_mz,
                    seed_rt_seconds=scan.rt,
                    gapfill_status="fixed_window",
                    ms1_quant_status=(
                        "ambiguous_polarity" if extraction_polarity is None else "no_ms1_scans"
                    ),
                    trace=np.empty((0, 2)) if config.group_ions else None,
                )
            )
            measurements.append(_Measurement(best, extraction_polarity, refine=refined))
            insort(seed_index, (scan.precursor_mz, best))
        assignments[i] = MS2Assignment(best, "seeded")
        stats["ms2_seeded_scans"] += 1

    # One bounded target schedule per polarity permits interleaved acquisition
    # times between polarities. Existing MS1 validation guarantees order within it.
    schedules = {}
    for polarity in polarities:
        ids = [i for i, m in enumerate(measurements) if m.polarity == polarity]
        if ids:
            ids.sort(key=lambda i: features[measurements[i].feature_index].rt_start_seconds)
            schedules[polarity] = [ids, 0, {}, []]
    targets = np.array([features[m.feature_index].mz for m in measurements])
    retained_values = 0

    def finish_refined(measurement):
        nonlocal retained_values
        if not measurement.refine or measurement.finalized:
            return
        feature = features[measurement.feature_index]
        effective = config_for(feature.polarity, config, configs)
        refine_measurement(feature, measurement, effective)
        _impute_missing_height(feature, measurement, effective)
        retained_values -= measurement.retained_values
        measurement.retained_values = 0
        measurement.samples.clear()
        measurement.finalized = True

    if schedules:
        stats["ms2_gapfill_passes"] = 1
        with closing(
            scan_source.ms1()
            if scan_source is not None
            else iter_ms1(source, assume_centroid=config.assume_centroid)
        ) as scans:
            for scan in scans:
                if scan.polarity not in schedules:
                    continue
                schedule = schedules[scan.polarity]
                ids, cursor, active, ends = schedule
                while cursor < len(ids):
                    i = ids[cursor]
                    feature = features[measurements[i].feature_index]
                    if feature.rt_start_seconds > scan.rt:
                        break
                    active[i] = None
                    heapq.heappush(ends, (feature.rt_end_seconds, i))
                    cursor += 1
                schedule[1] = cursor
                while ends and ends[0][0] < scan.rt:
                    expired = heapq.heappop(ends)[1]
                    active.pop(expired, None)
                    finish_refined(measurements[expired])
                if not active:
                    continue
                query_ids = np.fromiter(active, dtype=np.int64, count=len(active))
                effective = config_for(scan.polarity, config, configs)
                raw_mz, raw_y = scan.mz, scan.intensity
                if not scan.centroided and any(measurements[i].refine for i in query_ids):
                    zeros = np.zeros(len(raw_mz))
                    raw_mz, raw_y, _ = centroid_profile(
                        raw_mz,
                        raw_y,
                        zeros,
                        zeros,
                        config.profile_min_points,
                        float(config.profile_gap_factor),
                    )
                values = eic_intensities(
                    scan.mz,
                    scan.intensity,
                    targets[query_ids],
                    np.array([_radius(targets[i], effective) for i in query_ids]),
                )
                thresholds = regional_thresholds(
                    scan.mz, scan.intensity, targets[query_ids], effective
                )
                for i, value, threshold in zip(query_ids, values, thresholds):
                    measurement = measurements[i]
                    feature = features[measurement.feature_index]
                    measurement.scans += 1
                    insort(
                        measurement.noise_references,
                        (abs(scan.rt - feature.seed_rt_seconds), scan.rt, float(threshold)),
                    )
                    del measurement.noise_references[3:]
                    if measurement.refine:
                        if measurement.sampling_limit:
                            continue
                        radius = _radius(targets[i], effective)
                        a, b = (
                            np.searchsorted(raw_mz, targets[i] - radius),
                            np.searchsorted(raw_mz, targets[i] + radius, side="right"),
                        )
                        measurement.sampled_values += b - a
                        if (
                            b - a > 64
                            or measurement.sampled_values > 100000
                            or len(measurement.samples) >= 10000
                        ):
                            measurement.sampling_limit = True
                            retained_values -= measurement.retained_values
                            measurement.retained_values = 0
                            measurement.samples.clear()
                        else:
                            measurement.samples.append(
                                (scan.rt, raw_mz[a:b].copy(), raw_y[a:b].copy())
                            )
                            measurement.retained_values += b - a
                            retained_values += b - a
                            stats["ms2_gapfill_peak_retained_values"] = max(
                                stats["ms2_gapfill_peak_retained_values"], int(retained_values)
                            )
                        continue
                    if measurement.previous_rt is not None:
                        delta = scan.rt - measurement.previous_rt
                        if delta <= effective.gapfill_max_gap_seconds:
                            measurement.area += (
                                delta * (measurement.previous_intensity + value) * 0.5
                            )
                    measurement.previous_rt = scan.rt
                    measurement.previous_intensity = float(value)
                    if value > 0:
                        feature.scan_count += 1
                        if value > (feature.peak_height or 0.0):
                            feature.peak_height = float(value)
                            feature.rt_apex_seconds = scan.rt
                    if config.group_ions:
                        measurement.trace.append((float(value), scan.rt))
    for measurement in measurements:
        feature = features[measurement.feature_index]
        effective = config_for(feature.polarity, config, configs)
        if measurement.refine:
            finish_refined(measurement)
            continue
        feature.gapfill_mz_radius_da = _radius(feature.mz, effective)
        if measurement.scans:
            if not np.isfinite(measurement.area):
                raise ValueError("MS2-seeded MS1 peak area overflow")
            feature.peak_height = feature.peak_height or 0.0
            feature.peak_area = measurement.area
            feature.intensity = (
                feature.peak_height if config.quant_method == "peak_height" else feature.peak_area
            )
            feature.ms1_quant_status = "gap_filled" if feature.scan_count else "no_ms1_signal"
            if config.group_ions:
                feature.trace = np.asarray(measurement.trace)
            _impute_missing_height(feature, measurement, effective)
    assignments = merge_refined(features, assignments, first_feature, config, configs)
    stats["ms2_seeded_features"] = len(features) - first_feature
    stats["ms2_gapfill_sample_points"] = sum(m.scans for m in measurements)
    stats["ms2_gapfill_status_counts"] = dict(
        sorted(Counter(f.ms1_quant_status for f in features[first_feature:]).items())
    )
    apply_ms2_assignments(features, records, assignments, config, stats)
    return assignments
