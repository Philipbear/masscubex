"""v0.6 SciPy calibration oracle for differential numerical tests."""

import numpy as np
from scipy.ndimage import gaussian_filter1d
from scipy.signal import find_peaks
from masscubex.tuning import Representative


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
