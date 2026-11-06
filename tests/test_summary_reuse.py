"""Exact parity when reusing finalized traces and avoiding redundant parent copies."""

from types import FunctionType

import numpy as np
import pytest

from masscubex import Config
from masscubex.features import summarize_trace
from masscubex.upstream_adapter import ExtractionFeature
import upstream_masscube as upstream


def reference_segments(trace, config):
    """Unoptimized upstream workflow: recursive deep copies and repeated passes."""
    parent = upstream.Feature()
    parent.signals = trace[:, :2]
    parent.rt_seq = (trace[:, 2] / 60).tolist()
    parent.scan_idx_seq = list(range(len(trace)))
    parent.finalize(pa=False, g_score=False, a_score=False)
    parent_mass = float(parent.mz)

    def maxima(values):
        return upstream.gaussian_filter_and_prominence_method(
            values, config.smoothing_sigma, config.peak_prominence, config.peak_distance
        )

    split = FunctionType(
        upstream.segment_feature.__code__,
        {**upstream.segment_feature.__globals__, "gaussian_filter_and_prominence_method": maxima},
        argdefs=upstream.segment_feature.__defaults__,
    )
    segments = [parent]
    if config.segment:
        for _ in range(config.segmentation_passes):
            segments = [
                child
                for peak in segments
                for child in split(
                    peak,
                    noise_tol=config.peak_noise_score_limit if config.peak_shape_guard else np.inf,
                )
            ]
    for peak in segments:
        peak.finalize()
    return parent_mass, segments


@pytest.mark.parametrize(
    "config",
    [
        Config(min_scans=1),
        Config(min_scans=1, segment=False),
        Config(min_scans=1, peak_shape_guard=False),
        Config(min_scans=1, segmentation_passes=3),
        Config(min_scans=1, smoothing_sigma=1.2, peak_distance=5, peak_prominence=0.02),
        Config(min_scans=1, group_ions=False, quant_method="peak_area"),
    ],
)
def test_reused_summary_matches_original_full_precision(config):
    rng = np.random.default_rng(841)
    for n in [1, 2, 4, 5, 8, 20, 100, 500]:
        for shape in range(5):
            times = np.cumsum(rng.uniform(0.1, 2, n))
            intensity = rng.lognormal(4, 3, n)
            intensity[rng.random(n) < 0.2] = 0
            intensity[n // 2] = max(intensity.max(), 1)
            if shape < 2:
                x = np.arange(n)
                intensity = 1e5 * np.exp(-0.5 * ((x - n * 0.35) / (n * 0.12)) ** 2)
                intensity += 1e4 * np.exp(-0.5 * ((x - n * 0.65) / (n * 0.08)) ** 2)
            if shape == 1 and n > 5:
                intensity[:2] = intensity[-2:] = 0
            trace = np.column_stack(
                (100 + rng.uniform(-0.003, 0.003, n), intensity, times, rng.uniform(1, 5, n))
            )
            saved = trace.copy()
            parent_mass, expected = reference_segments(trace, config)
            actual = list(summarize_trace(trace, config, "positive"))
            assert len(actual) == len(expected)
            for feature, peak in zip(actual, expected):
                threshold = float(trace[peak.scan_idx, 3])
                values = {
                    "mz": float(peak.mz),
                    "rt_apex_seconds": peak.rt * 60,
                    "rt_start_seconds": peak.rt_seq[0] * 60,
                    "rt_end_seconds": peak.rt_seq[-1] * 60,
                    "peak_height": float(peak.peak_height),
                    "peak_area": float(peak.peak_area),
                    "scan_count": int(peak.length),
                    "apex_noise_threshold": threshold,
                    "apex_threshold_ratio": peak.peak_height / threshold,
                    "gaussian_similarity": float(peak.gaussian_similarity),
                    "_parent_mz": parent_mass,
                    "intensity": float(getattr(peak, config.quant_method)),
                }
                assert {key: getattr(feature, key) for key in values} == values
                guarded = config.peak_shape_guard and peak.noise_score > config.peak_noise_score_limit
                status = ("guarded_" if guarded else "") + (
                    "split" if peak.is_segmented else "unsplit"
                )
                assert feature.segmentation_status == (status if config.segment else "disabled")
                if config.group_ions:
                    np.testing.assert_array_equal(feature.trace[:, 0], peak.signals[:, 1])
                    np.testing.assert_array_equal(feature.trace[:, 1], trace[peak.scan_idx_seq, 2])
                else:
                    assert feature.trace is None
            np.testing.assert_array_equal(trace, saved)


def test_split_children_own_all_mutable_buffers_used_after_subsetting():
    parent = ExtractionFeature()
    parent.rt_seq = list(np.arange(100) / 60)
    parent.scan_idx_seq = list(range(100))
    x = np.arange(100)
    intensity = np.exp(-0.5 * ((x - 20) / 3) ** 2) * 1000
    intensity += np.exp(-0.5 * ((x - 70) / 3) ** 2) * 500
    parent.signals = np.column_stack((np.full(100, 100), intensity))
    parent.finalize()
    children = upstream.segment_feature(parent)
    assert len(children) == 2
    for child in children:
        assert child is not parent
        for name in ("rt_seq", "scan_idx_seq", "ms2_seq"):
            assert getattr(child, name) is not getattr(parent, name)
        for name in ("signals", "peak_shape"):
            assert not np.shares_memory(getattr(child, name), getattr(parent, name))
    for name in ("rt_seq", "scan_idx_seq", "ms2_seq"):
        assert getattr(children[0], name) is not getattr(children[1], name)
    for name in ("signals", "peak_shape"):
        assert not np.shares_memory(getattr(children[0], name), getattr(children[1], name))


@pytest.mark.parametrize("dtype", [np.float32, np.float64])
@pytest.mark.parametrize("minimum", [1, 2, 4])
def test_short_span_arithmetic_matches_upstream_at_precision_boundaries(dtype, minimum):
    settings = Config(min_scans=minimum)
    cases = [
        [1.0],
        [1.0, 1.0],  # First tied maximum remains the apex.
        [np.nextafter(np.float32(1), np.float32(0)), 1, 1.00000006, 1.1],
        [1.0, 0.0, 2.0, 0.5],  # Interior missing scans count toward the span.
        [1e-38, 2e-38, 1e-38],
        [1e30, 1.1e30, 1e30],
    ]
    for values in cases:
        y = np.array([0, 0, *values, *([0] * 10)], dtype=dtype)
        if np.count_nonzero(y > 0) < minimum:
            continue
        n = len(y)
        trace = np.column_stack(
            (
                100 + np.arange(n) * 0.000003814697265625,
                y,
                119.731 + np.arange(n) * 0.17357,
                np.full(n, 0.2),
            )
        ).astype(dtype)
        _, expected = reference_segments(trace, settings)
        (peak,) = expected
        (actual,) = summarize_trace(trace, settings, "negative")
        assert actual.mz == float(peak.mz)
        assert actual.peak_height == float(peak.peak_height)
        assert actual.peak_area == float(peak.peak_area)
        assert actual.scan_count == int(peak.length)
        assert actual.rt_apex_seconds == peak.rt * 60
        assert actual.rt_start_seconds == peak.rt_seq[0] * 60
        assert actual.rt_end_seconds == peak.rt_seq[-1] * 60
        assert actual.gaussian_similarity == peak.gaussian_similarity == 0
        assert actual.segmentation_status == "unsplit"
        np.testing.assert_array_equal(actual.trace[:, 0], peak.signals[:, 1])
        np.testing.assert_array_equal(actual.trace[:, 1], trace[peak.scan_idx_seq, 2])


def test_sparse_long_span_keeps_upstream_scores_and_splitting_path():
    # Three observed ions spanning seven acquisitions are not a short span.
    y = np.array([0, 2, 0, 0, 10, 0, 0, 2, 0, 0], dtype=float)
    trace = np.column_stack((np.full(len(y), 100), y, np.arange(len(y)), np.ones(len(y))))
    config = Config(min_scans=3)
    _, expected = reference_segments(trace, config)
    (actual,) = summarize_trace(trace, config)
    assert actual.gaussian_similarity == expected[0].gaussian_similarity > 0
    assert actual.scan_count == 3
