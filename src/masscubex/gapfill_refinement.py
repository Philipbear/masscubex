"""Refine a broad MS2 query into one narrow, local, unfiltered MS1 measurement."""

from bisect import bisect_left, bisect_right, insort

import numpy as np
from numba import njit
from scipy.ndimage import gaussian_filter1d
from scipy.signal import find_peaks

from .ms2 import MS2Assignment
from .tuning import config_for


def _radius(mass, config):
    return max(config.gapfill_mz_tolerance_da, mass * config.gapfill_mz_tolerance_ppm * 1e-6)


@njit(cache=True, nogil=True)
def _eic(masses, intensities, offsets, center, radius):
    """Query packed raw scan slices without Python dispatch per observation."""
    y = np.zeros(len(offsets) - 1)
    for index in range(len(y)):
        start, stop = offsets[index], offsets[index + 1]
        mz = masses[start:stop]
        left = start + np.searchsorted(mz, center - radius)
        right = start + np.searchsorted(mz, center + radius, side="right")
        for point in range(left, right):
            y[index] = max(y[index], intensities[point])
    return y


def _local_bounds(t, y, seed_rt, config):
    """Smooth in physical time for boundaries only; quantification stays raw."""
    if len(t) < 3 or y.max() <= 0:
        return 0, len(t) - 1
    step = float(np.median(np.diff(t)))
    count = min(4 * len(t), max(2, int(np.ceil((t[-1] - t[0]) / step)) + 1))
    grid = np.linspace(t[0], t[-1], count)
    step = grid[1] - grid[0]
    smooth = np.interp(grid, t, y)
    sigma = config.smoothing_seconds / step if config.smoothing_seconds else config.smoothing_sigma
    smooth = gaussian_filter1d(smooth, max(sigma, 1e-6))
    distance = (
        max(1, int(np.ceil(config.peak_distance_seconds / step)))
        if config.peak_distance_seconds
        else config.peak_distance
    )
    peaks, _ = find_peaks(
        np.r_[0.0, smooth, 0.0], prominence=smooth.max() * config.peak_prominence, distance=distance
    )
    peaks = peaks - 1
    peaks = peaks[(peaks >= 0) & (peaks < count)]
    if not len(peaks):
        return 0, len(t) - 1
    valleys = (
        [0]
        + [int(a + np.argmin(smooth[a : b + 1])) for a, b in zip(peaks[:-1], peaks[1:])]
        + [count - 1]
    )
    choices = [i for i in range(len(peaks)) if grid[valleys[i]] <= seed_rt <= grid[valleys[i + 1]]]
    if not choices:
        choices = list(range(len(peaks)))
    chosen = min(choices, key=lambda i: (abs(grid[peaks[i]] - seed_rt), -smooth[peaks[i]]))
    left, right, apex = valleys[chosen], valleys[chosen + 1], peaks[chosen]
    baseline = float(np.quantile(smooth[left : right + 1], 0.1))
    level = baseline + 0.05 * (smooth[apex] - baseline)
    # Flat traces remain explicit forced measurements; no invented baseline subtraction.
    if smooth[apex] > baseline:
        a = b = apex
        while a > left and smooth[a] > level:
            a -= 1
        while b < right and smooth[b] > level:
            b += 1
        left, right = a, b
    a = max(0, int(np.searchsorted(t, grid[left], side="right")) - 1)
    b = min(len(t) - 1, int(np.searchsorted(t, grid[right])))
    return a, b


def refine_measurement(feature, measurement, config):
    if measurement.sampling_limit:
        feature.ms1_quant_status = feature.gapfill_status = "sampling_limit"
        return
    samples = measurement.samples
    if not samples:
        return
    t = np.array([sample[0] for sample in samples])
    masses = np.concatenate([s[1] for s in samples])
    if not len(masses):
        feature.peak_height = feature.peak_area = feature.intensity = 0.0
        feature.ms1_quant_status = feature.gapfill_status = "no_ms1_signal"
        return
    intensities = np.concatenate([sample[2] for sample in samples])
    offsets = np.empty(len(samples) + 1, dtype=np.int64)
    offsets[0] = 0
    np.cumsum([len(sample[1]) for sample in samples], out=offsets[1:])
    ordered_masses = np.sort(masses)
    # Bound each cluster's full diameter: no chain of close masses can span the query.
    centers = []
    start = 0
    while start < len(masses):
        end = int(
            np.searchsorted(
                ordered_masses,
                ordered_masses[start] + 2 * _radius(ordered_masses[start], config),
                side="right",
            )
        )
        centers.append(float(np.median(ordered_masses[start:end])))
        start = end
    dt = float(np.median(np.diff(t))) if len(t) > 1 else 0.0
    local = np.abs(t - feature.seed_rt_seconds) <= max(config.ms2_rt_tolerance, 2 * dt)
    candidates = []
    search = max(
        config.tuning_search_da,
        feature.seed_precursor_mz * config.tuning_search_ppm * 1e-6,
        _radius(feature.seed_precursor_mz, config),
    )
    for center in centers:
        if abs(center - feature.seed_precursor_mz) > search:
            continue
        y = _eic(masses, intensities, offsets, center, _radius(center, config))
        height = float(np.max(y[local], initial=0.0))
        if height > 0:
            candidates.append((height, center, y))
    if not candidates:
        # An ion elsewhere in the broad time query is not evidence for this MS2.
        feature.ms1_quant_status = feature.gapfill_status = "no_local_ms1_signal"
        return
    candidates.sort(key=lambda item: (-item[0], abs(item[1] - feature.seed_precursor_mz)))
    # Prefer candidates inside the independently learned precursor tolerance.
    precursor_radius = max(
        config.ms2_mz_tolerance_da, feature.seed_precursor_mz * config.ms2_mz_tolerance_ppm * 1e-6
    )
    close = [c for c in candidates if abs(c[1] - feature.seed_precursor_mz) <= precursor_radius]
    candidates = close or candidates
    if len(candidates) > 1 and candidates[0][0] < 3 * candidates[1][0]:
        feature.ms1_quant_status = feature.gapfill_status = "ambiguous_mass"
        return
    _, center, y = candidates[0]
    feature.mz = center
    feature.gapfill_mz_radius_da = _radius(center, config)
    a, b = _local_bounds(t, y, feature.seed_rt_seconds, config)
    if (
        not t[a] - config.ms2_rt_tolerance
        <= feature.seed_rt_seconds
        <= t[b] + config.ms2_rt_tolerance
    ):
        feature.ms1_quant_status = feature.gapfill_status = "no_local_peak"
        return
    t, y = t[a : b + 1], y[a : b + 1]
    deltas = np.diff(t)
    area = float(
        np.sum(
            deltas[deltas <= config.gapfill_max_gap_seconds]
            * ((y[1:] + y[:-1]) * 0.5)[deltas <= config.gapfill_max_gap_seconds]
        )
    )
    if not np.isfinite(area):
        raise ValueError("MS2-seeded MS1 peak area overflow")
    apex = int(np.argmax(y))
    feature.rt_start_seconds, feature.rt_end_seconds = float(t[0]), float(t[-1])
    feature.rt_apex_seconds = float(t[apex])
    feature.peak_height, feature.peak_area = float(y[apex]), area
    feature.scan_count = int(np.count_nonzero(y > 0))
    feature.intensity = feature.peak_height if config.quant_method == "peak_height" else area
    feature.ms1_quant_status = "gap_filled" if feature.scan_count else "no_ms1_signal"
    feature.gapfill_status = "peak_bounds" if feature.scan_count > 1 else "single_observation"
    if config.group_ions:
        feature.trace = np.column_stack((y, t))


def merge_refined(features, assignments, first, config, configs):
    """Deduplicate only seeds that recover the same polarity, mass and MS1 apex."""
    retained = features[:first]
    mapping = {i: i for i in range(first)}
    index = []
    for old_index, feature in enumerate(features[first:], first):
        effective = config_for(feature.polarity, config, configs)
        radius = _radius(feature.mz, effective)
        matches = []
        if feature.gapfill_status == "peak_bounds" and feature.polarity != "unknown":
            left = bisect_left(index, (feature.mz - radius, -1))
            right = bisect_right(index, (feature.mz + radius, len(retained)))
            for _, j in index[left:right]:
                previous = retained[j]
                if (
                    previous.polarity == feature.polarity
                    and abs(previous.rt_apex_seconds - feature.rt_apex_seconds) <= 1e-6
                    and previous.rt_start_seconds - effective.ms2_rt_tolerance
                    <= feature.seed_rt_seconds
                    <= previous.rt_end_seconds + effective.ms2_rt_tolerance
                    and feature.rt_start_seconds - effective.ms2_rt_tolerance
                    <= previous.seed_rt_seconds
                    <= feature.rt_end_seconds + effective.ms2_rt_tolerance
                ):
                    matches.append(j)
        if matches:
            mapping[old_index] = min(matches)
        else:
            mapping[old_index] = len(retained)
            if feature.gapfill_status == "peak_bounds":
                insort(index, (feature.mz, len(retained)))
            retained.append(feature)
    features[:] = retained
    return [
        MS2Assignment(mapping[a.feature_index], a.reason) if a.feature_index is not None else a
        for a in assignments
    ]
