"""Readable SciPy compatibility path for ambiguous peak-distance ties.

This is the unchanged 0.6 segmentation/quantification path. It also provides
an independent numerical oracle for the compiled chromatogram implementation.
"""

from functools import lru_cache
import numpy as np
from numba import njit
from scipy.ndimage import correlate1d
from scipy.signal import find_peaks
from .config import Config
from .features import Feature


@lru_cache(maxsize=512)
def _gaussian_weights(sigma):
    radius = int(4.0 * sigma + 0.5)
    positions = np.arange(-radius, radius + 1)
    weights = np.exp((-0.5 / (sigma * sigma)) * positions**2)
    weights /= weights.sum()
    weights.flags.writeable = False
    return weights


def _smooth_padded(signal, sigma):
    """Same zero pad and order-0 reflect Gaussian as scipy.gaussian_filter1d.

    Cache only the small kernel, never a signal or a peak result. A bounded
    cache avoids rebuilding identical kernels for every chromatographic trace.
    """
    padded = np.empty(len(signal) + 2, dtype=np.float64)
    padded[0] = padded[-1] = 0.0
    padded[1:-1] = signal
    return correlate1d(padded, _gaussian_weights(float(sigma)), mode="reflect")


@njit(cache=True, nogil=True)
def _has_multiple_maxima(values):
    """Count local maxima with SciPy's plateau definition, stopping after two.

    With zero/one maximum, prominence and distance cannot produce a split, so
    the expensive peak-property calculation is unnecessary. Multiple maxima
    always use the original SciPy prominence/distance algorithm.
    """
    peaks = 0
    i = 1
    while i < len(values) - 1:
        if values[i] > values[i - 1]:
            end = i
            while end < len(values) - 1 and values[end] == values[end + 1]:
                end += 1
            if end < len(values) - 1 and values[end] > values[end + 1]:
                peaks += 1
                if peaks > 1:
                    return True
            i = end
        i += 1
    return False


@njit(cache=True, nogil=True)
def _median_spacing(rt):
    # Trace columns are strided; Numba's np.diff requires a contiguous input.
    differences = np.empty(len(rt) - 1, dtype=np.float64)
    for i in range(len(differences)):
        differences[i] = rt[i + 1] - rt[i]
    return np.median(differences)


def summarize_reference_trace(trace, config: Config, polarity="unknown"):
    """Split at chromatographic valleys; quantify original, un-smoothed signals.

    Adjacent segments share a valley point, but no trapezoid interval is counted
    twice. Area is intensity * seconds, and missing MS1 scans contribute zeros.
    """
    valid = np.flatnonzero(trace[:, 1] > 0)
    if len(valid) < config.min_scans:
        return
    trace = trace[valid[0] : valid[-1] + 1]
    intensity = trace[:, 1]
    cuts = [0]
    if config.segment:
        rt = trace[:, 2]
        physical = config.smoothing_seconds > 0 or config.peak_distance_seconds > 0
        if physical and len(rt) > 1:
            step = float(_median_spacing(rt))
            size = min(len(rt) * 4, max(2, int(np.ceil((rt[-1] - rt[0]) / step - 1e-9)) + 1))
            grid = np.linspace(rt[0], rt[-1], size)
            step = float(grid[1] - grid[0])
            signal = np.interp(grid, rt, intensity)
            sigma = (
                config.smoothing_seconds / step
                if config.smoothing_seconds > 0
                else config.smoothing_sigma
            )
            distance = (
                max(1, int(np.ceil(config.peak_distance_seconds / step - 1e-9)))
                if config.peak_distance_seconds > 0
                else config.peak_distance
            )
        else:
            grid, signal, sigma, distance = (
                rt,
                intensity,
                config.smoothing_sigma,
                config.peak_distance,
            )
        smooth = _smooth_padded(signal, sigma)
        peaks = np.empty(0, dtype=np.int64)
        if _has_multiple_maxima(smooth):
            peaks, _ = find_peaks(
                smooth, prominence=smooth.max() * config.peak_prominence, distance=distance
            )
        peaks = peaks - 1
        peaks = peaks[(peaks >= 0) & (peaks < len(signal))]
        for left, right in zip(peaks[:-1], peaks[1:]):
            valley = left + int(np.argmin(smooth[left + 1 : right + 2]))
            if physical:
                index = int(np.searchsorted(rt, grid[valley]))
                index = min(index, len(rt) - 1)
                if index and abs(rt[index - 1] - grid[valley]) <= abs(rt[index] - grid[valley]):
                    index -= 1
                valley = index
            if cuts[-1] < valley < len(intensity) - 1:
                cuts.append(valley)
    cuts.append(len(trace) - 1)
    for left, right in zip(cuts[:-1], cuts[1:]):
        part = trace[left : right + 1]
        y = part[:, 1]
        count = int(np.count_nonzero(y))
        if count < config.min_scans:
            continue
        apex = int(np.argmax(y))
        rt = part[:, 2]
        height = float(y[apex])
        threshold = float(part[apex, 3])
        area = float(np.dot(np.diff(rt), (y[:-1] + y[1:]) * 0.5))
        yield Feature(
            mz=float(np.dot(part[:, 0], y / y.sum())),
            rt_apex_seconds=float(rt[apex]),
            rt_start_seconds=float(rt[0]),
            rt_end_seconds=float(rt[-1]),
            peak_height=height,
            peak_area=area,
            scan_count=count,
            apex_noise_threshold=threshold,
            apex_threshold_ratio=height / threshold if threshold > 0 else None,
            polarity=polarity,
            intensity=height if config.quant_method == "peak_height" else area,
            quant_method=config.quant_method,
            trace=part[:, 1:3].copy() if config.group_ions else None,
        )
