"""Shape evidence for chromatographic splits, separate from MS1 mass-axis noise.

An unconvincing split is suppressed, never used to discard a mass trace. All
measurements still use the original observations, including missing-scan zeros.
"""

import numpy as np
from numba import njit


@njit(cache=True)
def trace_noise_score(values):
    """MassCube's total-variation score; >3 triggers scrutiny, not rejection."""
    if len(values) < 5 or values.max() <= 0:
        return 0.0
    height = values.max()
    previous = values[0] / height if values[0] >= 0.05 * height else 0.0
    variation = previous
    for index in range(1, len(values)):
        current = values[index] / height if values[index] >= 0.05 * height else 0.0
        variation += abs(current - previous)
        previous = current
    return max(0.0, (variation + previous) * 0.5 - 1)


@njit(cache=True)
def gaussian_similarity(times, values, asymmetric=False):
    """MassCube's apex/FWHM Gaussian Pearson score, clamped to [0, 1].

    The optional two-sided template is used only by the split guard, allowing
    fronting/tailing peaks. Exported scores always use the symmetric template.
    This is a diagnostic, not a fitted area or a chemical identification.
    """
    if len(values) < 5:
        return 0.0
    apex = np.argmax(values)
    height = values[apex]
    if height <= 0:
        return 0.0
    left, right = apex, apex
    while left > 0 and values[left] >= height * 0.5:
        left -= 1
    while right < len(values) - 1 and values[right] >= height * 0.5:
        right += 1
    sigma = max((times[right] - times[left]) / 2.355, 1e-6)
    observed = values / height
    left_sigma = right_sigma = sigma
    if asymmetric:
        # At half height, |t - apex| / sigma = sqrt(2 log(2)).
        left_sigma = max((times[apex] - times[left]) / np.sqrt(2 * np.log(2)), 1e-6)
        right_sigma = max((times[right] - times[apex]) / np.sqrt(2 * np.log(2)), 1e-6)
    model = np.empty(len(times))
    for index in range(len(times)):
        width = left_sigma if times[index] < times[apex] else right_sigma
        model[index] = np.exp(-0.5 * ((times[index] - times[apex]) / width) ** 2)
    mean_y, mean_g = np.mean(observed), np.mean(model)
    covariance = variance_y = variance_g = 0.0
    for index in range(len(values)):
        y = observed[index] - mean_y
        g = model[index] - mean_g
        covariance += y * g
        variance_y += y * y
        variance_g += g * g
    denominator = np.sqrt(variance_y * variance_g)
    if denominator == 0:
        return 0.0
    return min(1.0, max(0.0, covariance / denominator))


@njit(cache=True)
def temporal_noise(times, values):
    """Robust per-observation noise from local linear interpolation residuals.

    Only triples of consecutive, observed (positive) MS1 signals contribute.
    A normalized residual removes a locally linear trend and accounts for
    unequal scan spacing. MAD / Phi^-1(0.75) estimates Gaussian noise sigma.
    Fewer than five residuals provide no reliable estimate: use shape alone.
    """
    residuals = np.empty(max(0, len(values) - 2))
    count = 0
    for index in range(1, len(values) - 1):
        if min(values[index - 1], values[index], values[index + 1]) <= 0:
            continue
        weight = (times[index] - times[index - 1]) / (times[index + 1] - times[index - 1])
        prediction = (1 - weight) * values[index - 1] + weight * values[index + 1]
        residuals[count] = (values[index] - prediction) / np.sqrt(1 + (1 - weight) ** 2 + weight**2)
        count += 1
    if count < 5:
        return 0.0
    residuals = residuals[:count]
    return np.median(np.abs(residuals - np.median(residuals))) / 0.6744897501960817


@njit(cache=True)
def interpolate_missing(times, values):
    """Bridge the tracker's short missing-scan runs for split detection only."""
    observed = values > 0
    return np.interp(times, times[observed], values[observed])


@njit(cache=True)
def filter_shape_peaks(times, values, grid, smooth, peaks, minimum_similarity, minimum_snr):
    """Require a supported rise and fall before accepting a split candidate.

    Evaluate each candidate inside its original adjacent-valley basin. Rejecting
    a neighbor must not enlarge that basin and make a fluctuation look better.
    The zero padding used by peak finding supplies no baseline evidence.
    """
    if len(peaks) < 2:
        return peaks
    noise = temporal_noise(times, values)
    observed = values > 0
    observed_times, observed_values = times[observed], values[observed]
    valleys = np.empty(len(peaks) + 1, dtype=np.int64)
    valleys[0], valleys[-1] = 0, len(grid) - 1
    for index in range(len(peaks) - 1):
        left, right = peaks[index], peaks[index + 1]
        valleys[index + 1] = left + np.argmin(smooth[left : right + 1]) - 1
    selected = np.empty(len(peaks), dtype=np.int64)
    count = 0
    for index in range(len(peaks)):
        peak = peaks[index]
        left, right = valleys[index], valleys[index + 1]
        baseline = max(smooth[left + 1], smooth[right + 1])
        amplitude = smooth[peak] - baseline
        if amplitude <= 0 or amplitude < minimum_snr * noise:
            continue
        first = np.searchsorted(observed_times, grid[left])
        stop = np.searchsorted(observed_times, grid[right], side="right")
        if stop - first < 5:
            continue
        rt = observed_times[first:stop]
        y = np.maximum(observed_values[first:stop] - baseline, 0.0)
        # Interpolated/smoothed points cannot manufacture sampling support.
        # Half the RAW excess height excludes a lone spike on noisy background.
        support_level = 0.5 * y.max()
        support = 0
        rising = falling = False
        for sample in range(len(y)):
            if y[sample] > support_level:
                support += 1
                rising |= rt[sample] < grid[peak - 1]
                falling |= rt[sample] > grid[peak - 1]
        if support < 3 or not rising or not falling:
            continue
        score = gaussian_similarity(rt, y)
        if score < minimum_similarity and gaussian_similarity(rt, y, True) < minimum_similarity:
            continue
        selected[count] = peak
        count += 1
    return selected[:count]
