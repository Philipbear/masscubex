"""Chromatographic evidence: suppress background splits without deleting traces."""

from dataclasses import replace

import numpy as np
import pytest

from masscubex import Config
from masscubex.features import summarize_trace
from masscubex.peak_shapes import gaussian_similarity, trace_noise_score
import upstream_masscube as upstream
from masscubex.reference_chromatograms import reference_gaussian_similarity


def trace(times, values):
    return np.column_stack((np.full(len(values), 100.0), values, times, np.ones(len(values))))


def upstream_peak(times, values):
    peak = upstream.Feature()
    peak.signals = np.column_stack((np.full(len(values), 100), values))
    peak.rt_seq = (times / 60).tolist()
    peak.scan_idx_seq = list(range(len(times)))
    peak.finalize()
    return peak


@pytest.mark.parametrize("seed", range(10))
@pytest.mark.parametrize("missing", [False, True])
def test_persistent_background_is_one_trace_with_original_quantification(seed, missing):
    rng = np.random.default_rng(seed)
    times = np.cumsum(rng.uniform(0.2, 0.5, 2000))
    values = 10000 + rng.normal(0, 1500, len(times))
    if missing:
        values[17:-1:23] = 0
        values[18:-1:23] = 0
    source = trace(times, values)
    config = Config()
    legacy = list(summarize_trace(source, replace(config, peak_shape_guard=False)))
    (feature,) = summarize_trace(source, config)
    assert len(legacy) > 40
    assert feature.segmentation_status == "guarded_unsplit"
    assert feature.gaussian_similarity < 0.3
    expected = upstream_peak(times, values)
    assert feature.rt_start_seconds == expected.rt_seq[0] * 60
    assert feature.rt_end_seconds == expected.rt_seq[-1] * 60
    assert feature.peak_height == expected.peak_height
    assert feature.peak_area == expected.peak_area
    assert feature.scan_count == np.count_nonzero(values)
    np.testing.assert_array_equal(feature.trace[:, 0], values.astype(np.float32))
    np.testing.assert_array_equal(feature.trace[:, 1], times)


@pytest.mark.parametrize("spacing", [0.1, 0.3, 0.8])
@pytest.mark.parametrize("tail_factor", [1, 4, 10])
def test_noisy_trace_is_never_split_even_when_it_contains_real_peaks(spacing, tail_factor):
    times = np.arange(0, 650, spacing)
    centers = np.array([50, 170, 290, 410, 530])
    values = sum(
        100 * np.exp(-0.5 * ((times - center) / np.where(times < center, 2, 2 * tail_factor)) ** 2)
        for center in centers
    )
    assert trace_noise_score(values) > 3
    config = Config()
    (feature,) = summarize_trace(trace(times, values), config)
    assert feature.segmentation_status == "guarded_unsplit"
    assert feature.peak_area == upstream_peak(times, values).peak_area


@pytest.mark.parametrize("kind", ["broad", "overlapping", "weak_on_background", "tailing"])
def test_plausible_two_peaks_are_preserved(kind):
    rng = np.random.default_rng(42)
    times = np.arange(0, 200, 0.2)
    centers = [60, 90] if kind == "overlapping" else [50, 150]
    width = 10 if kind in {"broad", "overlapping"} else 2
    values = sum(
        200
        * np.exp(
            -0.5
            * (
                (times - center)
                / np.where((times > center) & (kind == "tailing"), 10 * width, width)
            )
            ** 2
        )
        for center in centers
    )
    if kind == "weak_on_background":
        values += 1000 + rng.normal(0, 10, len(times))
    fs = list(summarize_trace(trace(times, values), Config()))
    assert len(fs) == (1 if trace_noise_score(values) > 3 else 2)
    if len(fs) == 2:
        np.testing.assert_allclose([f.rt_apex_seconds for f in fs], centers, atol=1)


def test_isolated_spikes_do_not_make_chromatographic_splits():
    rng = np.random.default_rng(290)
    times = np.arange(2000) * 0.3
    values = 1000 + rng.normal(0, 100, len(times))
    values[[200, 600, 1000, 1400, 1800]] = 3000
    (feature,) = summarize_trace(trace(times, values), Config())
    assert feature.peak_height == 3000
    assert feature.segmentation_status == "guarded_unsplit"


def test_gap_bridge_is_only_for_split_detection():
    times = np.arange(2000) * 0.3
    values = np.full(2000, 1000.0)
    values[10:-10:20] = 0
    values[11:-10:20] = 0
    (feature,) = summarize_trace(trace(times, values), Config())
    assert feature.peak_area == upstream_peak(times, values).peak_area
    assert feature.scan_count == np.count_nonzero(values)
    np.testing.assert_array_equal(feature.trace[:, 0], values)


@pytest.mark.parametrize("scale", [1e-12, 1.0, 1e12])
def test_gaussian_similarity_is_scale_invariant_and_matches_reference(scale):
    times = np.linspace(0, 30, 301)
    values = scale * np.exp(-0.5 * ((times - 14) / 2) ** 2)
    result = gaussian_similarity(times, values)
    assert result > 0.999
    assert result == pytest.approx(reference_gaussian_similarity(times, values), abs=1e-14)
    assert gaussian_similarity(times, np.full(len(times), scale)) == 0
    assert gaussian_similarity(times[:4], values[:4]) == 0


def test_shape_guard_settings_are_validated_and_round_trip():
    assert Config.from_dict(Config().to_dict()) == Config()
    for key in ["peak_noise_score_limit"]:
        with pytest.raises(ValueError):
            Config(**{key: -1})
        with pytest.raises(ValueError):
            Config(**{key: float("nan")})
    with pytest.raises(ValueError):
        Config(peak_shape_guard="yes")


@pytest.mark.parametrize("seed", range(20))
def test_split_boundaries_match_pinned_upstream(seed):
    rng = np.random.default_rng(seed)
    t = np.arange(0, 100, 0.3)
    y = np.full(len(t), 1.0)
    for center, amplitude, width in zip(
        rng.uniform(10, 90, 4), rng.uniform(50, 200, 4), rng.uniform(1, 5, 4)
    ):
        y += amplitude * np.exp(-0.5 * ((t - center) / width) ** 2)
    if seed % 3 == 0:
        y += rng.uniform(0, 20, len(t))
    y = y.astype(np.float32).astype(np.float64)
    original = upstream.Feature()
    original.signals = np.column_stack((np.full(len(y), 100.0), y)).astype(np.float32)
    original.rt_seq = t / 60
    original.scan_idx_seq = np.arange(len(t))
    original.finalize(pa=False, ta=False, g_score=False, a_score=False)
    segments = [original]
    for _ in range(2):
        segments = [part for f in segments for part in upstream.segment_feature(f)]
    expected = [(f.rt_seq[0] * 60, f.rt_seq[-1] * 60) for f in segments]
    result = list(summarize_trace(trace(t, y), Config(auto_tune=False)))
    assert [(f.rt_start_seconds, f.rt_end_seconds) for f in result] == expected
