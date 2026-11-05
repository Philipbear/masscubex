"""Experimental profile-specific detection thresholds; never edit raw samples.

Estimate baseline and electronic variation on dense, regular profile segments.
Zero-padded/censored segments cannot identify that variation. There use only
the existing physical mass-peak detection guards, followed by the common filter,
rather than apply a second, unsupported amplitude cutoff to exported peaks.
"""

from collections import Counter
import math

import numpy as np
from scipy.ndimage import median_filter

from .centroid import centroid_profile
from .noise import estimate_noise
from .upstream_adapter import merge_centroids


def profile_thresholds(mz, intensity, config, *, median_width=31, minimum_samples=64):
    baseline, threshold = np.zeros(len(mz)), np.zeros(len(mz))
    status = Counter()
    keys = np.floor(mz / 100)
    for key in np.unique(keys):
        indices = np.flatnonzero(keys == key)
        if len(indices) < minimum_samples:
            status["sparse_bin_fallback_points"] += len(indices)
            continue
        local_mz = mz[indices].astype(float)
        spacing = np.diff(local_mz)
        ordinary_spacing = float(np.median(spacing))
        boundaries = np.r_[0, np.flatnonzero(spacing > 4 * ordinary_spacing) + 1, len(indices)]
        lengths = np.diff(boundaries)
        short_points = int(lengths[lengths < minimum_samples].sum())
        if short_points:
            status["short_segment_fallback_points"] += short_points
        # Sparse exports have hundreds of thousands of short padded segments.
        # Their lengths identify abstention without allocating/casting samples.
        for position in np.flatnonzero(lengths >= minimum_samples):
            start, stop = boundaries[position : position + 2]
            region = indices[start:stop]
            values = intensity[region].astype(float)
            if np.mean(values > 0) < 0.95:
                status["censored_segment_fallback_points"] += len(region)
                continue
            second_difference = np.diff(values, n=2)
            center = np.median(second_difference)
            # For independent electronic variation Var(Δ²y) = 6σ². The
            # Gaussian MAD constant applies to the *difference residuals*,
            # not to all peak heights in the m/z region.
            difference_sigma = float(np.median(abs(second_difference - center))) / (
                0.6744897501960817 * math.sqrt(6)
            )
            local_baseline = median_filter(values, size=median_width, mode="nearest")
            residual = values - local_baseline
            residual_center = float(np.median(residual))
            residual_sigma = float(np.median(abs(residual - residual_center))) / 0.6744897501960817
            # Difference variance alone underestimates smooth/colored noise.
            # Keep the larger residual estimate as a conservative check.
            sigma = max(difference_sigma, residual_sigma)
            if not np.isfinite(sigma) or sigma <= 0:
                status["unidentified_variation_fallback_points"] += len(region)
                continue
            baseline[region] = local_baseline + residual_center
            threshold[region] = baseline[region] + 3 * sigma
            status["profile_specific_points"] += len(region)
    return baseline, threshold, dict(status)


def detect_profile(scan, config, mode):
    mz, intensity = scan.mz.astype(np.float32), scan.intensity.astype(np.float32)
    if mode == "post_only":
        baseline = np.zeros(len(mz))
        threshold = np.zeros(len(mz))
        status = {"permissive_detection_points": len(mz)}
    elif mode == "profile_specific":
        baseline, threshold, status = profile_thresholds(mz, intensity, config)
    elif mode == "existing":
        baseline, _, threshold = estimate_noise(mz, intensity, config)
        status = {"existing_profile_points": len(mz)}
    else:
        raise ValueError(mode)
    if mode == "profile_specific":
        # Seed at μ + 3σ, but follow the physical mass-peak shoulders down to
        # μ + σ. Requiring three samples above the *apex* threshold clips weak
        # real peaks; this retains the same apex significance and raw height.
        support = baseline + (threshold - baseline) / 3
        detected = centroid_profile(
            mz,
            intensity,
            baseline,
            support,
            config.profile_min_points,
            config.profile_gap_factor,
            apex_threshold=threshold,
        )
    else:
        detected = centroid_profile(
            mz, intensity, baseline, threshold, config.profile_min_points, config.profile_gap_factor
        )
    return (*merge_centroids(*detected, config.centroid_mz_tolerance_da), status)
