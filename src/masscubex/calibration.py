"""Compiled quality screening for robust calibration representatives.

Selection thresholds and per-polarity estimation policy live in tuning.py.
This module only measures one original trace and explains a rejection.
"""

from typing import NamedTuple
import numpy as np
from numba import njit
from .chromatograms import gaussian_smooth, select_peaks


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
def measure_representative(trace, search_da, search_ppm):
    """Require a resolved, complete rise/fall and a stable central mass trace."""
    if len(trace) < 7:
        return None, "short_trace"
    mz, y, t, competing = trace[:, 0], trace[:, 1], trace[:, 2], trace[:, 3]
    if not np.all(np.isfinite(trace)) or np.max(y) <= 0:
        return None, "invalid_trace"
    smooth = gaussian_smooth(y, 0.6)
    baseline = float(np.quantile(y, 0.1))
    apex = int(np.argmax(smooth))
    amplitude = float(smooth[apex] - baseline)
    if amplitude <= 0.5 * float(y.max()):
        return None, "flat_or_low_contrast"
    peaks, _ = select_peaks(smooth, 1, amplitude * 0.3)
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
    if np.isnan(half[0]) or np.isnan(outer[0]):
        return None, "truncated_peak"
    if np.sum(np.abs(np.diff(smooth))) > 4 * amplitude:
        return None, "rough_trace"
    center = float(np.median(mz[central]))
    error = float(np.quantile(np.abs(mz[central] - center), 0.95))
    search = max(search_da, center * search_ppm * 1e-6)
    if error > search * 0.4:
        return None, "mass_near_search_boundary"
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
    ), None
