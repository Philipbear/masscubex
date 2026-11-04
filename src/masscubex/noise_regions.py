"""One stable MS1 noise estimate per fixed m/z region, including empty queries."""

import numpy as np
from numba import njit


@njit(cache=True)
def _moments(values):
    # Center before summation: a large baseline must not erase a small spread.
    origin = values[0]
    shifted = values - origin
    center = np.mean(shifted)
    variance = np.mean((shifted - center) ** 2)
    mean, std = origin + center, np.sqrt(variance)
    if not np.isfinite(mean + 3 * std):
        raise ValueError("MS1 intensity range overflows noise-estimation moments")
    return mean, std


@njit(cache=True)
def _region_moments(intensity, fraction, minimum, clipped, iterations, clip_sigma):
    values = intensity[intensity > 0]
    if len(values) < minimum:
        return 0.0, 0.0, False
    if not clipped:
        count = max(minimum, int(np.ceil(fraction * len(values))))
        selected = np.partition(values, count - 1)[:count]
        mean, std = _moments(selected)
    else:
        for _ in range(iterations):
            mean, std = _moments(values)
            values = values[values <= mean + clip_sigma * std]
            if len(values) < minimum:
                return 0.0, 0.0, False
        mean, std = _moments(values)
    return mean, std, True


@njit(cache=True, nogil=True)
def estimate_bins(mz, intensity, width, fraction, minimum, floor, clipped, iterations, clip_sigma):
    """Half-open bins [k*width, (k+1)*width), anchored at m/z zero.

    All points in a bin share an estimate, including its low-tail points. There
    is no per-candidate exclusion. Zero padding is never a background sample.
    """
    mean, std = np.zeros(len(mz)), np.zeros(len(mz))
    threshold = np.full(len(mz), floor)
    start = 0
    while start < len(mz):
        key = np.floor(mz[start] / width)
        end = start + 1
        while end < len(mz) and np.floor(mz[end] / width) == key:
            end += 1
        mu, sd, supported = _region_moments(
            intensity[start:end], fraction, minimum, clipped, iterations, clip_sigma
        )
        if supported:
            mean[start:end], std[start:end] = mu, sd
            threshold[start:end] = max(floor, mu + 3 * sd)
        start = end
    return mean, std, threshold


@njit(cache=True, nogil=True)
def _query_thresholds(
    mz, intensity, targets, width, fraction, minimum, floor, bins, clipped, iterations, clip_sigma
):
    result = np.full(len(targets), floor)
    previous_key, previous_value = np.nan, floor
    # Reuse a bin's estimate when several simultaneous MS2 targets share it.
    for i in np.argsort(targets):
        target = targets[i]
        key = np.floor(target / width) if bins else target
        if key == previous_key:
            result[i] = previous_value
            continue
        lower = key * width if bins else target - width / 2
        upper = (key + 1) * width if bins else target + width / 2
        start = np.searchsorted(mz, lower)
        end = np.searchsorted(mz, upper) if bins else np.searchsorted(mz, upper, side="right")
        mu, sd, supported = _region_moments(
            intensity[start:end], fraction, minimum, clipped, iterations, clip_sigma
        )
        if supported:
            result[i] = max(floor, mu + 3 * sd)
        previous_key, previous_value = key, result[i]
    return result


def regional_thresholds(mz, intensity, targets, config):
    """Evaluate noise where an MS2 target has no centroid, without adding fake peaks.

    Inputs are already normalized raw MS1 arrays. Profile inputs remain profile
    samples, exactly as in the first-pass noise estimate before centroiding.
    """
    if config.noise_mode == "none":
        return np.full(len(targets), config.min_intensity)
    return _query_thresholds(
        mz,
        intensity,
        targets,
        float(config.noise_window_da),
        float(config.noise_fraction),
        config.noise_min_neighbors,
        float(config.min_intensity),
        config.noise_layout == "bins",
        config.noise_estimator == "clipped",
        config.noise_iterations,
        float(config.noise_clip_sigma),
    )
