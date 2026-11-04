"""Chromatographic evidence: suppress background splits without deleting traces."""

from dataclasses import replace

import numpy as np
import pytest

from masscubex import Config
from masscubex.features import summarize_trace
from masscubex.peak_shapes import gaussian_similarity, temporal_noise, trace_noise_score
from masscubex.reference_chromatograms import reference_gaussian_similarity


def trace(times, values):
    return np.column_stack((np.full(len(values), 100.0), values, times, np.ones(len(values))))


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
    config = Config(smoothing_seconds=0.3, peak_distance_seconds=1.5)
    legacy = list(summarize_trace(source, replace(config, peak_shape_guard=False)))
    (feature,) = summarize_trace(source, config)
    assert len(legacy) > 40
    assert feature.segmentation_status == "guarded_unsplit"
    assert feature.gaussian_similarity < 0.3
    assert feature.rt_start_seconds == times[0]
    assert feature.rt_end_seconds == times[-1]
    assert feature.peak_height == values.max()
    assert feature.peak_area == pytest.approx(
        np.dot(np.diff(times), (values[:-1] + values[1:]) / 2)
    )
    assert feature.scan_count == np.count_nonzero(values)
    np.testing.assert_array_equal(feature.trace, source[:, 1:3])


@pytest.mark.parametrize("spacing", [0.1, 0.3, 0.8])
@pytest.mark.parametrize("tail_factor", [1, 4, 10])
def test_many_real_peaks_survive_whole_trace_noise_score_guard(spacing, tail_factor):
    times = np.arange(0, 650, spacing)
    centers = np.array([50, 170, 290, 410, 530])
    values = sum(
        100 * np.exp(-0.5 * ((times - center) / np.where(times < center, 2, 2 * tail_factor)) ** 2)
        for center in centers
    )
    assert trace_noise_score(values) > 3
    config = Config(smoothing_seconds=0.4, peak_distance_seconds=3)
    features = list(summarize_trace(trace(times, values), config))
    assert len(features) == len(centers)
    np.testing.assert_allclose([f.rt_apex_seconds for f in features], centers, atol=spacing)
    assert sum(f.peak_area for f in features) == pytest.approx(np.trapz(values, times))


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
    fs = list(summarize_trace(trace(times, values), Config(smoothing_seconds=0.4)))
    assert len(fs) == 2
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
    assert feature.peak_area == pytest.approx(np.trapz(values, times))
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


def test_temporal_noise_estimate_handles_uneven_spacing_and_linear_trends():
    rng = np.random.default_rng(251)
    times = np.cumsum(rng.uniform(0.1, 0.5, 20000))
    values = 10000 + 0.5 * times + rng.normal(0, 100, len(times))
    assert temporal_noise(times, values) == pytest.approx(100, rel=0.04)
    values[10:-10:20] = 0
    assert temporal_noise(times, values) == pytest.approx(100, rel=0.04)
    assert temporal_noise(times[:5], values[:5]) == 0


def test_shape_guard_settings_are_validated_and_round_trip():
    assert Config.from_dict(Config().to_dict()) == Config()
    for key in ["peak_min_gaussian_similarity", "peak_min_snr", "peak_noise_score_limit"]:
        with pytest.raises(ValueError):
            Config(**{key: -1})
        with pytest.raises(ValueError):
            Config(**{key: float("nan")})
    with pytest.raises(ValueError):
        Config(peak_min_gaussian_similarity=1.1)
    with pytest.raises(ValueError):
        Config(peak_shape_guard="yes")
