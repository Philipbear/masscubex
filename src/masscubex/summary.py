"""MassCube trace rules using intervals instead of temporary feature objects.

The scientific functions and NumPy arithmetic are unchanged. Intervals avoid
deep copies, repeated parent measurements and peak-shape arrays that extraction
never consumes. Source functions remain pinned in upstream.py.
"""

from functools import lru_cache

import numpy as np
from numba import njit
from scipy.ndimage import gaussian_filter1d

try:
    from scipy.ndimage import _filters, _nd_image
except ImportError:  # Future SciPy versions may reorganize private helpers.
    _filters = _nd_image = None

from . import upstream
from .chromatograms import select_peaks
from .measurements import gaussian_similarity, integrated_area, noise_score_at_apex


@lru_cache(maxsize=128)
def gaussian_weights(sigma):
    """Same NumPy weight construction as scipy.ndimage's order-zero kernel."""
    radius = int(4.0 * sigma + 0.5)
    weights = np.ascontiguousarray(_filters._gaussian_kernel1d(sigma, 0, radius)[::-1])
    weights.flags.writeable = False
    return weights


def smooth_padded(values, weights):
    """Call the same SciPy C kernel, retaining platform-specific accumulation."""
    signal = np.empty(len(values) + 2, dtype=np.float64)
    signal[0] = signal[-1] = 0.0
    signal[1:-1] = values
    result = np.empty(len(signal))
    _nd_image.correlate1d(signal, weights, 0, result, 2, 0.0, 0)
    return result


@njit(cache=True, nogil=True)
def _select_raw_peaks(values, smooth, distance, prominence, float32_values):
    """Select smoothed peaks and apply the unchanged raw-intensity median gate."""
    peaks, ambiguous = select_peaks(smooth, distance, smooth.max() * prominence)
    if ambiguous or not len(peaks):
        return peaks, ambiguous
    # NumPy median propagates a NaN anywhere in the raw signal. Its comparison
    # then rejects every peak, including those away from that observation.
    if np.isnan(values).any():
        return np.empty(0, dtype=np.int64), False
    middle = len(values) // 2
    partitioned = np.partition(values, middle)
    median = np.float64(partitioned[middle])
    if len(values) % 2 == 0:
        lower = partitioned[:middle].max()
        if float32_values:
            # NumPy's float32 mean adds before dividing, including overflow.
            median = np.float64(np.float32(np.float32(lower + partitioned[middle]) / np.float32(2)))
        else:
            median = (np.float64(lower) + np.float64(partitioned[middle])) / 2.0
    selected = np.empty(len(peaks), dtype=np.int64)
    count = 0
    for peak in peaks:
        raw_peak = peak - 1
        if values[raw_peak] > median:
            selected[count] = raw_peak
            count += 1
    return selected[:count], False


def maxima(values, config):
    if config.peak_distance > np.iinfo(np.int64).max:
        return upstream.gaussian_filter_and_prominence_method(
            values, config.smoothing_sigma, config.peak_prominence, config.peak_distance
        )
    try:
        smooth = smooth_padded(values, gaussian_weights(config.smoothing_sigma))
    except (AttributeError, TypeError):
        # Private SciPy API changes must affect speed only, never the result.
        smooth = gaussian_filter1d(np.r_[0.0, values.astype(float), 0.0], config.smoothing_sigma)
    peaks, ambiguous = _select_raw_peaks(
        values, smooth, config.peak_distance, config.peak_prominence, values.dtype == np.float32
    )
    if ambiguous:
        # Preserve the installed NumPy/SciPy unstable priority order for ties.
        return upstream.gaussian_filter_and_prominence_method(
            values, config.smoothing_sigma, config.peak_prominence, config.peak_distance
        )
    return peaks


@njit(cache=True, nogil=True)
def _span_state(y, left, right):
    """Trim zeros, count observations and select the first maximum in one pass."""
    while left < right and y[left] == 0:
        left += 1
    while right > left and y[right - 1] == 0:
        right -= 1
    apex, count = left, 0
    if left < right:
        height = y[left]
        for i in range(left, right):
            value = y[i]
            if value > 0:
                count += 1
            # NumPy argmax selects the first NaN as well as the first tie.
            if value > height or (np.isnan(value) and not np.isnan(height)):
                apex, height = i, value
    return left, right, count, apex


def summarize_intervals(trace, config, polarity, feature_type, status):
    signals = np.asarray(trace[:, :2], dtype=np.float32)
    y = signals[:, 1]
    left, right, count, apex = _span_state(y, 0, len(y))
    if left == right:
        from .features import _summarize_trace_reference

        yield from _summarize_trace_reference(trace, config, polarity)
        return
    minutes = (trace[:, 2] / 60.0).tolist()
    parent_mass = float(signals[apex, 0])
    # left, right, was split, terminal, cached noise score, count, first apex.
    intervals = [(left, right, False, False, None, count, apex)]
    noise_limit = config.peak_noise_score_limit if config.peak_shape_guard else np.inf
    if config.segment:
        for _ in range(config.segmentation_passes):
            following = []
            for left, right, split, terminal, noise, count, apex in intervals:
                if terminal:
                    following.append((left, right, split, terminal, noise, count, apex))
                    continue
                values = y[left:right]
                noise = noise_score_at_apex(values, apex - left)
                if count < 5 or noise > noise_limit:
                    following.append((left, right, split, True, noise, count, apex))
                    continue
                peaks = maxima(values, config)
                if len(peaks) < 2:
                    following.append((left, right, split, True, noise, count, apex))
                    continue
                valleys = [left]
                valleys.extend(
                    left + int(a + np.argmin(values[a:b])) for a, b in zip(peaks[:-1], peaks[1:])
                )
                valleys.append(right - 1)
                for a, b in zip(valleys[:-1], valleys[1:]):
                    a, b, observed, top = _span_state(y, a, b + 1)
                    following.append((a, b, True, False, None, observed, top))
            intervals = following
            if all(item[3] for item in intervals):
                break
    for left, right, split, _, noise, count, apex in intervals:
        values = y[left:right]
        if count < config.min_scans:
            continue
        mass = float(signals[apex, 0])
        height = float(int(y[apex]))
        times = minutes[left:right]
        area = integrated_area(times, values)
        if noise is None:
            noise = noise_score_at_apex(values, apex - left)
        threshold = float(trace[apex, 3])
        yield feature_type(
            mz=mass,
            rt_apex_seconds=float(minutes[apex] * 60.0),
            rt_start_seconds=float(times[0] * 60.0),
            rt_end_seconds=float(times[-1] * 60.0),
            peak_height=height,
            peak_area=area,
            scan_count=count,
            apex_noise_threshold=threshold,
            apex_threshold_ratio=height / threshold if threshold > 0 else None,
            polarity=polarity,
            intensity=height if config.quant_method == "peak_height" else area,
            quant_method=config.quant_method,
            gaussian_similarity=float(gaussian_similarity(times, values, apex - left)),
            segmentation_status=status(
                config.segment, split, config.peak_shape_guard and noise > noise_limit
            ),
            trace=np.column_stack((values, trace[left:right, 2])) if config.group_ions else None,
            _parent_mz=parent_mass,
        )
