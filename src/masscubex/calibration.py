"""Compiled quality screening for robust calibration representatives.

Selection thresholds and per-polarity estimation policy live in tuning.py.
This module only measures one original trace and explains a rejection.
"""

from typing import NamedTuple
import numpy as np
from numba import njit
from scipy.ndimage import gaussian_filter1d
from .chromatograms import gaussian_smooth
from .peak_shapes import gaussian_similarity


class RepresentativeValues(NamedTuple):
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
    gaussian_similarity: float


def reference_measurements(trace):
    """Preserve NumPy/SciPy rounding for the small pool of accepted candidates.

    Calibration values set later mass/time boundaries. Re-measure only accepted
    peaks with the original arithmetic so tiny kernel rounding differences do
    not propagate into every seeded query's window or learned mass radius.
    """
    mz, intensity, times, _ = trace.T
    smooth = gaussian_filter1d(intensity, 0.6)
    baseline = float(np.quantile(intensity, 0.1))
    apex = int(np.argmax(smooth))
    amplitude = float(smooth[apex] - baseline)
    crossing = getattr(_crossings, "py_func", _crossings)
    half = crossing(times, smooth, apex, baseline + amplitude * 0.5)
    outer = crossing(times, smooth, apex, baseline + amplitude * 0.1)
    if np.isnan(outer[0]):
        outer = crossing(times, intensity, int(np.argmax(intensity)), baseline + amplitude * 0.1)
    central = (intensity >= baseline + amplitude * 0.5) & (times >= half[0]) & (times <= half[1])
    center = float(np.median(mz[central]))
    error = float(np.quantile(np.abs(mz[central] - center), 0.95))
    return RepresentativeValues(
        center,
        float(times[apex]),
        half[0],
        half[1],
        float(intensity.max()),
        half[1] - half[0],
        outer[1] - outer[0],
        error,
        error / center * 1e6,
        int(np.count_nonzero(central)),
        float(gaussian_similarity(times, np.maximum(intensity - baseline, 0.0))),
    )


@njit(cache=True)
def _crossings(t, y, apex, level):
    """Interpolate rise/fall times; NaNs mark an incomplete crossing pair."""
    left = right = apex
    while left > 0 and y[left] > level:
        left -= 1
    while right < len(y) - 1 and y[right] > level:
        right += 1
    if y[left] > level or y[right] > level or left == apex or right == apex:
        return np.nan, np.nan
    a = t[left] + (t[left + 1] - t[left]) * (level - y[left]) / (y[left + 1] - y[left])
    b = t[right - 1] + (t[right] - t[right - 1]) * (y[right - 1] - level) / (
        y[right - 1] - y[right]
    )
    return float(a), float(b)


@njit(cache=True, nogil=True)
def measure_representative(trace, search_da, search_ppm, minimum_similarity=0.8):
    """Require a resolved, complete rise/fall and a stable central mass trace."""
    if len(trace) < 5 or np.count_nonzero(trace[:, 1] > 0) < 5:
        return None, "short_trace"
    mz, y, t, competing = trace[:, 0], trace[:, 1], trace[:, 2], trace[:, 3]
    if not np.all(np.isfinite(trace)) or np.max(y) <= 0:
        return None, "invalid_trace"
    smooth = gaussian_smooth(y, 0.6)
    baseline = float(np.quantile(y, 0.1))
    apex = int(np.argmax(smooth))
    amplitude = float(smooth[apex] - baseline)
    if amplitude <= 0:
        return None, "flat_trace"
    if np.count_nonzero(y == y.max()) >= max(3, len(y) // 3):
        return None, "plateau"
    observed_flank = (y > 0) & (y <= baseline + amplitude * 0.5)
    if not np.any(observed_flank & (t < t[apex])) or not np.any(observed_flank & (t > t[apex])):
        return None, "truncated_peak"
    half = _crossings(t, smooth, apex, baseline + amplitude * 0.5)
    outer = _crossings(t, smooth, apex, baseline + amplitude * 0.1)
    # Smoothing can lift both endpoints of a five-scan peak above 10% height.
    # Real raw flanks can still establish its outer width; never invent padding.
    if np.isnan(outer[0]):
        outer = _crossings(t, y, int(np.argmax(y)), baseline + amplitude * 0.1)
    if np.isnan(half[0]) or np.isnan(outer[0]):
        return None, "truncated_peak"
    # Only the main apex's central basin supplies mass measurements. A secondary
    # peak elsewhere in a trace cannot inflate its mass error or sample count.
    central = (y >= baseline + amplitude * 0.5) & (t >= half[0]) & (t <= half[1])
    if np.count_nonzero(central) < 3:
        return None, "undersampled_peak"
    if np.any(competing[central] > 0):
        return None, "competing_mass"
    center = float(np.median(mz[central]))
    error = float(np.quantile(np.abs(mz[central] - center), 0.95))
    search = max(search_da, center * search_ppm * 1e-6)
    if error > search * 0.4:
        return None, "mass_near_search_boundary"
    similarity = gaussian_similarity(t, np.maximum(y - baseline, 0.0))
    if similarity < minimum_similarity:
        return None, "poor_gaussian_shape"
    return RepresentativeValues(
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
        similarity,
    ), None
