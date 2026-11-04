"""Shape evidence for chromatographic splits, separate from MS1 mass-axis noise.

An unconvincing split is suppressed, never used to discard a mass trace. All
measurements still use the original observations, including missing-scan zeros.
"""

import numpy as np
from numba import njit


@njit(cache=True)
def trace_noise_score(values):
    """MassCube's total-variation score; >3 prevents splitting."""
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
def gaussian_similarity(times, values):
    """MassCube's apex/FWHM Gaussian Pearson score, clamped to [0, 1].

    This diagnostic is also used for mass-calibration anchor selection.
    It is not a splitting gate or a fitted quantitative area.
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
    model = np.empty(len(times))
    for index in range(len(times)):
        model[index] = np.exp(-0.5 * ((times[index] - times[apex]) / sigma) ** 2)
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
