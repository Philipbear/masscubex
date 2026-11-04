"""Adapters around pinned MassCube scientific functions; no shared mutation."""

from functools import lru_cache
from dataclasses import replace
import gzip
from types import FunctionType

import numpy as np

from . import upstream
from .centroid import centroid_profile
from .noise import estimate_noise


class ExtractionFeature(upstream.Feature):
    """Upstream feature with a cheap clone for its immediate subset operation.

    This private extraction state has no MS2 records or annotations attached.
    The pinned splitter deep-copies a parent and immediately calls ``subset``:
    that replaces both sequence lists and the MS2 list, then ``finalize`` copies
    the signal array and rebuilds peak_shape before modifying the child. Until
    then these buffers are only read. Copying the attributes is therefore enough;
    recursively copying the entire parent once per child is redundant.

    Do not use this class for general mutable upstream feature state. Its copy
    contract is specific to the pinned segmentation/subset/finalize functions.
    """

    def __deepcopy__(self, memo):
        child = type(self).__new__(type(self))
        child.__dict__ = self.__dict__.copy()
        memo[id(self)] = child
        return child


def resolve_ms2_filter(source, config):
    """Match upstream's first-300-lines Orbitrap keyword/default selection."""
    if config.ms2_min_intensity is not None:
        return config
    opener = gzip.open if str(source).lower().endswith(".gz") else open
    with opener(source, "rb") as stream:
        text = (
            b"".join(stream.readline() for _ in range(300)).decode("utf-8", errors="ignore").lower()
        )
    return replace(config, ms2_min_intensity=10000.0 if "orbitrap" in text else 500.0)


@lru_cache(maxsize=128)
def segmenter(sigma, prominence, distance):
    """Bind explicit parameter overrides without changing upstream function bodies."""
    if (sigma, prominence, distance) == (0.6, 0.03, 3):
        return upstream.segment_feature

    def maxima(values):
        return upstream.gaussian_filter_and_prominence_method(values, sigma, prominence, distance)

    return FunctionType(
        upstream.segment_feature.__code__,
        {**upstream.segment_feature.__globals__, "gaussian_filter_and_prominence_method": maxima},
        argdefs=upstream.segment_feature.__defaults__,
    )


def merge_centroids(mz, intensity, threshold, tolerance):
    """Upstream float32 adjacent-centroid merging; thresholds are diagnostics."""
    signals = np.column_stack((mz, intensity)).astype(np.float32)
    if tolerance <= 0 or len(signals) <= 1:
        return signals[:, 0], signals[:, 1], threshold
    order = signals[:, 0].argsort()
    starts = np.r_[0, np.flatnonzero(np.diff(signals[order, 0]) >= tolerance) + 1]
    merged = upstream.centroid_signals(signals, tolerance)
    return merged[:, 0], merged[:, 1], np.maximum.reduceat(threshold[order], starts)


def preprocess_ms1(scan, config):
    # Match the source's float32 scan representation before scientific decisions.
    mz = scan.mz.astype(np.float32)
    intensity = scan.intensity.astype(np.float32)
    baseline, _, threshold = estimate_noise(mz, intensity, config)
    if scan.centroided:
        keep = (intensity > threshold) & (mz > 0) & (mz < 100000)
        mz, intensity, threshold = mz[keep], intensity[keep], threshold[keep]
    else:
        if config.profile_mode != "centroid":
            raise ValueError(
                f"Profile MS1 spectrum {scan.native_id}; set --profile-mode centroid to process it"
            )
        mz, intensity, threshold = centroid_profile(
            mz, intensity, baseline, threshold, config.profile_min_points, config.profile_gap_factor
        )
    return merge_centroids(mz, intensity, threshold, config.centroid_mz_tolerance_da)


def preprocess_ms2(scan, config):
    """Original absolute/relative MS2 filter and adjacent-centroid merge."""
    if scan.mz is None:
        return
    mz, intensity = scan.mz.astype(np.float32), scan.intensity.astype(np.float32)
    absolute = 500.0 if config.ms2_min_intensity is None else config.ms2_min_intensity
    cutoff = max(absolute, float(intensity.max(initial=0)) * config.ms2_relative_intensity)
    keep = (mz > 0) & (mz < 100000) & (intensity > cutoff)
    mz, intensity = mz[keep], intensity[keep]
    scan.mz, scan.intensity, _ = merge_centroids(
        mz, intensity, np.zeros(len(mz)), config.centroid_mz_tolerance_da
    )
    # find_best_ms2 uses NumPy's float32 sum, not the separate float64 TIC.
    scan.fragment_intensity = float(np.sum(scan.intensity))
