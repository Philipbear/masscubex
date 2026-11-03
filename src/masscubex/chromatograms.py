"""Compiled chromatogram segmentation and float64 peak measurements.

The helpers implement the same Gaussian reflection, plateau, distance and
prominence rules as the SciPy reference path. No fast-math approximations are
enabled. Equal-height competing peaks use SciPy's own tie-breaking path.
"""

from typing import NamedTuple

import numpy as np
from numba import njit


class PeakSummary(NamedTuple):
    first: int
    last: int
    mz: float
    apex_time: float
    start_time: float
    end_time: float
    height: float
    area: float
    observations: int
    noise_threshold: float


@njit(cache=True)
def median_spacing(times):
    intervals = np.empty(len(times) - 1)
    for index in range(len(intervals)):
        intervals[index] = times[index + 1] - times[index]
    return np.median(intervals)


@njit(cache=True)
def uniform_grid(times):
    """Match the endpoint-preserving NumPy linspace used by the reference."""
    interval = median_spacing(times)
    duration = times[-1] - times[0]
    size = min(len(times) * 4, max(2, int(np.ceil(duration / interval - 1e-9)) + 1))
    grid = np.empty(size)
    step = duration / (size - 1)
    for index in range(size):
        grid[index] = index * step + times[0]
    grid[-1] = times[-1]
    return grid


@njit(cache=True, inline="always")
def reflected_index(index, length):
    """Half-sample reflection, including kernels wider than the signal."""
    folded = index % (2 * length)
    return folded if folded < length else 2 * length - 1 - folded


@njit(cache=True)
def gaussian_smooth(values, sigma, pad=False):
    """Order-zero Gaussian with reflected boundaries and optional zero ends."""
    radius = int(4.0 * sigma + 0.5)
    offsets = np.arange(-radius, radius + 1)
    weights = np.exp((-0.5 / (sigma * sigma)) * offsets**2)
    weights /= weights.sum()
    signal = values
    if pad:
        signal = np.empty(len(values) + 2)
        signal[0] = signal[-1] = 0.0
        signal[1:-1] = values
    result = np.empty(len(signal))
    for index in range(len(signal)):
        value = signal[index] * weights[radius]
        # Symmetric pairs, outside-in, as in SciPy's order-zero correlation.
        for offset in range(radius, 0, -1):
            left = reflected_index(index - offset, len(signal))
            right = reflected_index(index + offset, len(signal))
            value += (signal[left] + signal[right]) * weights[radius - offset]
        result[index] = value
    return result


@njit(cache=True)
def local_maxima(values):
    positions = np.empty(len(values) // 2, dtype=np.int64)
    count = 0
    index = 1
    while index < len(values) - 1:
        if values[index] > values[index - 1]:
            end = index
            while end + 1 < len(values) and values[end + 1] == values[index]:
                end += 1
            if end + 1 < len(values) and values[end] > values[end + 1]:
                positions[count] = (index + end) // 2
                count += 1
            index = end
        index += 1
    return positions[:count]


@njit(cache=True)
def select_peaks(signal, minimum_distance, minimum_prominence):
    """Apply distance before prominence; flag ambiguous distance-priority ties."""
    peaks = local_maxima(signal)
    keep = np.ones(len(peaks), dtype=np.bool_)
    if minimum_distance > 1:
        priorities = signal[peaks]
        ordered = np.argsort(priorities)
        # NumPy's unstable tie order depends on its version/SIMD implementation.
        # Also defer near ties: Gaussian arithmetic can differ by a few ulps.
        for index in range(1, len(ordered)):
            higher, lower = priorities[ordered[index]], priorities[ordered[index - 1]]
            if higher - lower <= 8 * np.finfo(np.float64).eps * max(abs(higher), abs(lower)):
                return peaks, True
        for index in range(len(ordered) - 1, -1, -1):
            candidate = ordered[index]
            if not keep[candidate]:
                continue
            neighbor = candidate - 1
            while neighbor >= 0 and peaks[candidate] - peaks[neighbor] < minimum_distance:
                keep[neighbor] = False
                neighbor -= 1
            neighbor = candidate + 1
            while neighbor < len(peaks) and peaks[neighbor] - peaks[candidate] < minimum_distance:
                keep[neighbor] = False
                neighbor += 1
    for candidate in range(len(peaks)):
        if not keep[candidate]:
            continue
        peak = peaks[candidate]
        height = signal[peak]
        left_min = right_min = height
        index = peak
        while index >= 0 and signal[index] <= height:
            left_min = min(left_min, signal[index])
            index -= 1
        index = peak
        while index < len(signal) and signal[index] <= height:
            right_min = min(right_min, signal[index])
            index += 1
        keep[candidate] = height - max(left_min, right_min) >= minimum_prominence
    return peaks[keep], False


@njit(cache=True)
def valley_boundaries(times, grid, smooth, peaks, physical_time):
    """Map first valley minima back to original observations, without overlap."""
    cuts = [0]
    for index in range(len(peaks) - 1):
        left, right = peaks[index] - 1, peaks[index + 1] - 1
        valley = left + np.argmin(smooth[left + 1 : right + 2])
        if physical_time:
            nearest = min(np.searchsorted(times, grid[valley]), len(times) - 1)
            if nearest and abs(times[nearest - 1] - grid[valley]) <= abs(
                times[nearest] - grid[valley]
            ):
                nearest -= 1
            valley = nearest
        if cuts[-1] < valley < len(times) - 1:
            cuts.append(valley)
    cuts.append(len(times) - 1)
    return cuts


@njit(cache=True, nogil=True)
def summarize_chromatogram(
    trace,
    minimum_scans,
    segment,
    smoothing_sigma,
    smoothing_seconds,
    peak_distance,
    peak_distance_seconds,
    relative_prominence,
):
    """Return named summaries and whether SciPy must resolve a tied peak order."""
    first = 0
    while first < len(trace) and trace[first, 1] <= 0:
        first += 1
    last = len(trace) - 1
    while last >= first and trace[last, 1] <= 0:
        last -= 1
    summaries = []
    if last - first + 1 < minimum_scans:
        return summaries, False
    values = trace[first : last + 1, 1]
    times = trace[first : last + 1, 2]
    cuts = [0, len(values) - 1]
    if segment:
        physical = smoothing_seconds > 0 or peak_distance_seconds > 0
        sigma, distance = smoothing_sigma, peak_distance
        grid, signal = times, values
        if physical and len(times) > 1:
            grid = uniform_grid(times)
            step = grid[1] - grid[0]
            signal = np.interp(grid, times, values)
            if smoothing_seconds > 0:
                sigma = smoothing_seconds / step
            if peak_distance_seconds > 0:
                distance = max(1, int(np.ceil(peak_distance_seconds / step - 1e-9)))
        smooth = gaussian_smooth(signal, sigma, pad=True)
        peaks, tied = select_peaks(smooth, distance, smooth.max() * relative_prominence)
        if tied:
            return summaries, True
        cuts = valley_boundaries(times, grid, smooth, peaks, physical)
    for left, right in zip(cuts[:-1], cuts[1:]):
        start, stop = first + left, first + right
        observations = 0
        total_intensity = area = 0.0
        apex = start
        # Quantify original observations in two strided loops. Avoid allocating
        # masses, intensities, time differences and trapezoids for every peak.
        for index in range(start, stop + 1):
            intensity = trace[index, 1]
            total_intensity += intensity
            observations += intensity > 0
            if intensity > trace[apex, 1]:
                apex = index
            if index > start:
                interval = trace[index, 2] - trace[index - 1, 2]
                area += interval * ((trace[index - 1, 1] + intensity) * 0.5)
        if observations < minimum_scans:
            continue
        mz = 0.0
        for index in range(start, stop + 1):
            mz += trace[index, 0] * (trace[index, 1] / total_intensity)
        summaries.append(
            PeakSummary(
                start,
                stop,
                mz,
                trace[apex, 2],
                trace[start, 2],
                trace[stop, 2],
                trace[apex, 1],
                area,
                observations,
                trace[apex, 3],
            )
        )
    return summaries, False
