"""Moving MS1 background from low-intensity neighbors along physical m/z.

The default uses the weakest fraction with a minimum count. A sliding rank
tree supplies low-tail counts/moments without sorting every window.
The previous clipped-moments estimator remains selectable for comparison.
"""

import numpy as np
from numba import njit

from .config import Config


@njit(cache=True, inline="always")
def _compensated_add(total, correction, value):
    updated = total + value
    if abs(total) >= abs(value):
        correction += (total - updated) + value
    else:
        correction += (value - updated) + total
    return updated, correction


@njit(cache=True, nogil=True)
def _local_moments(mz, intensity, valid, half_window, min_neighbors):
    n = mz.size
    mean = np.zeros(n)
    std = np.zeros(n)
    supported = np.zeros(n, dtype=np.bool_)
    left = 0
    right = 0
    total = 0.0
    squares = 0.0
    total_correction = 0.0
    squares_correction = 0.0
    count = 0
    for i in range(n):
        while right < n and mz[right] <= mz[i] + half_window:
            if valid[right]:
                total, total_correction = _compensated_add(
                    total, total_correction, intensity[right]
                )
                squares, squares_correction = _compensated_add(
                    squares, squares_correction, intensity[right] * intensity[right]
                )
                count += 1
            right += 1
        while left < right and mz[left] < mz[i] - half_window:
            if valid[left]:
                total, total_correction = _compensated_add(
                    total, total_correction, -intensity[left]
                )
                squares, squares_correction = _compensated_add(
                    squares, squares_correction, -intensity[left] * intensity[left]
                )
                count -= 1
            left += 1
        # Leave the candidate out of its own noise estimate.
        c = count - int(valid[i])
        if c >= min_neighbors:
            # Subtract before applying the compensation so an isolated huge
            # peak does not erase the low background after its window moves on.
            s = (total - (intensity[i] if valid[i] else 0.0)) + total_correction
            q = (squares - (intensity[i] * intensity[i] if valid[i] else 0.0)) + squares_correction
            mean[i] = max(0.0, s / c)
            std[i] = np.sqrt(max(0.0, q / c - mean[i] * mean[i]))
            supported[i] = True
    return mean, std, supported


@njit(cache=True, nogil=True)
def _estimate(mz, intensity, window, sigma, clip_sigma, iterations, minimum, floor):
    valid = np.ones(mz.size, dtype=np.bool_)
    for _ in range(iterations):
        mean, std, supported = _local_moments(mz, intensity, valid, window / 2, minimum)
        # Exclusions are monotonic. Unsupported windows keep their candidates.
        valid &= (~supported) | (intensity <= mean + clip_sigma * std)
    mean, std, supported = _local_moments(mz, intensity, valid, window / 2, minimum)
    threshold = np.full(mz.size, floor)
    for i in range(mz.size):
        if supported[i]:
            threshold[i] = max(floor, mean[i] + sigma * std[i])
    return mean, std, threshold


@njit(cache=True, inline="always")
def _merge_moments(n_a, mean_a, m2_a, n_b, mean_b, m2_b):
    """Combine centered moments without subtracting two large raw squares."""
    if n_a == 0:
        return n_b, mean_b, m2_b
    if n_b == 0:
        return n_a, mean_a, m2_a
    n = n_a + n_b
    delta = mean_b - mean_a
    mean = mean_a + delta * (n_b / n)
    m2 = m2_a + m2_b + delta * delta * (n_a * (n_b / n))
    return n, mean, m2


@njit(cache=True, inline="always")
def _tree_update(leaf, value, active, counts, means, m2):
    """Set one rank's membership, then rebuild centered moments to the root."""
    counts[leaf] = int(active)
    means[leaf] = value if active else 0.0
    m2[leaf] = 0.0
    node = leaf // 2
    while node:
        a, b = 2 * node, 2 * node + 1
        counts[node], means[node], m2[node] = _merge_moments(
            counts[a], means[a], m2[a], counts[b], means[b], m2[b]
        )
        node //= 2


@njit(cache=True, inline="always")
def _lowest_moments(k, counts, means, m2):
    """Select a prefix containing exactly k active ranks in O(log N)."""
    node, remaining, count = 1, k, 0
    mean = total_m2 = 0.0
    while remaining:
        if counts[node] <= remaining:
            count, mean, total_m2 = _merge_moments(
                count, mean, total_m2, counts[node], means[node], m2[node]
            )
            remaining -= counts[node]
            break
        left = 2 * node
        if counts[left] <= remaining:
            count, mean, total_m2 = _merge_moments(
                count, mean, total_m2, counts[left], means[left], m2[left]
            )
            remaining -= counts[left]
            node = left + 1
        else:
            node = left
    return mean, np.sqrt(max(0.0, total_m2 / k))


@njit(cache=True, nogil=True)
def _estimate_low_tail(mz, intensity, window, sigma, fraction, minimum, floor):
    n = len(mz)
    mean = np.zeros(n)
    std = np.zeros(n)
    threshold = np.full(n, floor)
    eligible = np.flatnonzero(intensity > 0)
    if len(eligible) < minimum:
        return mean, std, threshold
    order = eligible[np.argsort(intensity[eligible], kind="mergesort")]
    span = intensity[order[-1]] - intensity[order[0]]
    if not np.isfinite(span * span * len(eligible)):
        raise ValueError("MS1 intensity range overflows noise-estimation moments")
    base = 1
    while base < len(order):
        base *= 2
    ranks = np.zeros(n, dtype=np.int64)
    for i in range(len(order)):
        ranks[order[i]] = base + i
    counts = np.zeros(2 * base, dtype=np.int64)
    means = np.zeros(2 * base)
    m2 = np.zeros(2 * base)
    left = right = active = 0
    for i in range(n):
        while right < n and mz[right] <= mz[i] + window / 2:
            if ranks[right]:
                _tree_update(ranks[right], intensity[right], True, counts, means, m2)
                active += 1
            right += 1
        while left < right and mz[left] < mz[i] - window / 2:
            if ranks[left]:
                _tree_update(ranks[left], intensity[left], False, counts, means, m2)
                active -= 1
            left += 1
        # Positive profile samples / centroids are references; zero padding is
        # not. A candidate is excluded from its own background sample.
        count = active - int(ranks[i] > 0)
        if count < minimum:
            continue
        if ranks[i]:
            _tree_update(ranks[i], intensity[i], False, counts, means, m2)
        k = max(minimum, int(np.ceil(fraction * count)))
        mean[i], deviation = _lowest_moments(k, counts, means, m2)
        std[i] = deviation
        threshold[i] = max(floor, mean[i] + sigma * deviation)
        if ranks[i]:
            _tree_update(ranks[i], intensity[i], True, counts, means, m2)
    return mean, std, threshold


def estimate_noise(mz, intensity, config: Config | None = None):
    """Return baseline, standard deviation, threshold for sorted MS1 samples.

    Sparse windows use ``min_intensity`` alone: without neighboring background
    samples, a local noise floor cannot be identified reliably. Arrays are not
    modified. This is an empirical filter, not a false discovery probability.
    """
    config = config or Config()
    mz = np.ascontiguousarray(mz, dtype=np.float64)
    intensity = np.ascontiguousarray(intensity, dtype=np.float64)
    if mz.ndim != 1 or intensity.ndim != 1 or mz.shape != intensity.shape:
        raise ValueError("m/z and intensity must be equal-length one-dimensional arrays")
    if (
        not np.all(np.isfinite(mz))
        or not np.all(np.isfinite(intensity))
        or np.any(mz <= 0)
        or np.any(intensity < 0)
        or np.any(mz[1:] < mz[:-1])
    ):
        raise ValueError("Noise estimation requires sorted finite m/z and nonnegative intensity")
    if config.noise_mode == "none":
        return np.zeros(mz.size), np.zeros(mz.size), np.full(mz.size, config.min_intensity)
    if config.noise_estimator == "low_tail":
        return _estimate_low_tail(
            mz,
            intensity,
            float(config.noise_window_da),
            float(config.noise_sigma),
            float(config.noise_fraction),
            config.noise_min_neighbors,
            float(config.min_intensity),
        )
    return _estimate(
        mz,
        intensity,
        float(config.noise_window_da),
        float(config.noise_sigma),
        float(config.noise_clip_sigma),
        config.noise_iterations,
        config.noise_min_neighbors,
        float(config.min_intensity),
    )
