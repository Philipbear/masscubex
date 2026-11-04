"""SciPy-equivalent primitives for calibration and exact peak selection.

Extraction uses the installed SciPy smoothing kernel and falls back to SciPy
peak selection for ambiguous distance-priority ties.
"""

import numpy as np
from numba import njit


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
    if len(signal) == 0:
        return np.empty(0)
    # Reflect each position once. Doing this inside the convolution would need
    # two integer divisions per pair, for every output point and kernel offset.
    extended = np.empty(len(signal) + 2 * radius)
    for index in range(len(extended)):
        extended[index] = signal[reflected_index(index - radius, len(signal))]
    result = np.empty(len(signal))
    for index in range(len(signal)):
        value = signal[index] * weights[radius]
        # Symmetric pairs, outside-in, as in SciPy's order-zero correlation.
        for offset in range(radius, 0, -1):
            left, right = index + radius - offset, index + radius + offset
            value += (extended[left] + extended[right]) * weights[radius - offset]
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
