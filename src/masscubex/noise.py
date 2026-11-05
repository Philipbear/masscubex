"""MS1 background from low-intensity samples in physical m/z regions.

The default uses one estimate per fixed bin. The historical sliding rank tree
and clipped-moments estimator remain selectable for explicit comparisons.
"""

import numpy as np
from numba import njit

from .config import Config
from .array_validation import normalized_arrays
from .noise_regions import estimate_bins


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
def _tree_update(leaf, value, active, counts, means, m2, dirty):
    """Update membership immediately; defer moments until a query needs them.

    Most sliding-window changes affect intensities above the selected low tail.
    Their counts are needed for rank selection, but their moments are not.
    """
    change = int(active) - counts[leaf]
    counts[leaf] = int(active)
    means[leaf] = value if active else 0.0
    m2[leaf] = 0.0
    node = leaf // 2
    while node:
        counts[node] += change
        dirty[node] = True
        node //= 2


@njit(cache=True, inline="always")
def _refresh_moments(root, counts, means, m2, dirty):
    """Refresh a requested subtree in postorder without recursion or a stack.

    Clean subtrees are already correct. Empty subtrees need no descent. The
    merge order is identical to the eager tree, including float64 rounding.
    """
    node = root
    while True:
        if dirty[node]:
            left, right = 2 * node, 2 * node + 1
            if counts[node] == 0:
                means[node] = m2[node] = 0.0
            elif dirty[left]:
                node = left
                continue
            elif dirty[right]:
                node = right
                continue
            else:
                _, means[node], m2[node] = _merge_moments(
                    counts[left],
                    means[left],
                    m2[left],
                    counts[right],
                    means[right],
                    m2[right],
                )
            dirty[node] = False
        if node == root:
            return
        node //= 2


@njit(cache=True, inline="always")
def _lowest_moments(k, counts, means, m2, dirty):
    """Return moments and the last covered rank of the lowest k active leaves."""
    node, remaining, count = 1, k, 0
    last_subtree = 1
    mean = total_m2 = 0.0
    while remaining:
        if counts[node] <= remaining:
            _refresh_moments(node, counts, means, m2, dirty)
            count, mean, total_m2 = _merge_moments(
                count, mean, total_m2, counts[node], means[node], m2[node]
            )
            remaining -= counts[node]
            last_subtree = node
            break
        left = 2 * node
        if counts[left] <= remaining:
            _refresh_moments(left, counts, means, m2, dirty)
            count, mean, total_m2 = _merge_moments(
                count, mean, total_m2, counts[left], means[left], m2[left]
            )
            remaining -= counts[left]
            last_subtree = left
            node = left + 1
        else:
            node = left
    # Every active rank through this subtree belongs to the selected prefix.
    # Inactive leaves at its right edge do not affect candidate membership.
    while last_subtree < len(counts) // 2:
        last_subtree = 2 * last_subtree + 1
    return mean, np.sqrt(max(0.0, total_m2 / k)), last_subtree


@njit(cache=True, inline="always")
def _prefix_guard(k, counts):
    """Find the first unselected rank, or the tree's end if none exists.

    This extra rank fixes the query's subtree partition. Guarding only the
    selected ranks would preserve their values but could change the order of
    floating-point merges when a brighter neighboring rank leaves the window.
    """
    if counts[1] <= k:
        return len(counts) - 1
    node, remaining = 1, k + 1
    while node < len(counts) // 2:
        left = 2 * node
        if remaining <= counts[left]:
            node = left
        else:
            remaining -= counts[left]
            node = left + 1
    return node


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
    dirty = np.zeros(2 * base, dtype=np.bool_)
    left = right = active = 0
    cached_count, cached_last_rank, cached_guard = -1, 0, 0
    background = deviation = 0.0
    for i in range(n):
        while right < n and mz[right] <= mz[i] + window / 2:
            if ranks[right]:
                if ranks[right] <= cached_guard:
                    cached_count = -1
                _tree_update(ranks[right], intensity[right], True, counts, means, m2, dirty)
                active += 1
            right += 1
        while left < right and mz[left] < mz[i] - window / 2:
            if ranks[left]:
                if ranks[left] <= cached_guard:
                    cached_count = -1
                _tree_update(ranks[left], intensity[left], False, counts, means, m2, dirty)
                active -= 1
            left += 1
        # Positive profile samples / centroids are references; zero padding is
        # not. A candidate is excluded from its own background sample.
        count = active - int(ranks[i] > 0)
        if count < minimum:
            continue
        k = max(minimum, int(np.ceil(fraction * count)))
        # Reuse only when the size, selected ranks and first unselected rank
        # are unchanged. Higher ranks cannot alter the query or its moments.
        # The cached prefix includes the candidate; exclusion below is local
        # to this output point and never overwrites the cached values.
        if k != cached_count:
            background, deviation, cached_last_rank = _lowest_moments(k, counts, means, m2, dirty)
            cached_count = k
            cached_guard = _prefix_guard(k, counts)
        candidate_background, candidate_deviation = background, deviation
        # Most candidates exceed the low tail and cannot contribute to their
        # own estimate. Only rebuild the tree when the candidate is in that tail.
        excluded = 0 < ranks[i] <= cached_last_rank
        if excluded:
            _tree_update(ranks[i], intensity[i], False, counts, means, m2, dirty)
            candidate_background, candidate_deviation, _ = _lowest_moments(
                k, counts, means, m2, dirty
            )
        mean[i] = candidate_background
        std[i] = candidate_deviation
        threshold[i] = max(floor, candidate_background + sigma * candidate_deviation)
        if excluded:
            _tree_update(ranks[i], intensity[i], True, counts, means, m2, dirty)
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
    if not normalized_arrays(mz, intensity, False, False):
        raise ValueError("Noise estimation requires sorted finite m/z and nonnegative intensity")
    if config.noise_mode == "none":
        return np.zeros(mz.size), np.zeros(mz.size), np.full(mz.size, config.min_intensity)
    if config.noise_estimator == "consistency":
        raise ValueError(
            "Consistency noise requires neighboring scans; use process_file/extract CLI"
        )
    if config.noise_layout == "bins":
        return estimate_bins(
            mz,
            intensity,
            float(config.noise_window_da),
            float(config.noise_fraction),
            config.noise_min_neighbors,
            float(config.min_intensity),
            config.noise_estimator == "clipped",
            config.noise_iterations,
            float(config.noise_clip_sigma),
        )
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
