from dataclasses import replace

import numpy as np
import pytest

from masscubex.config import Config
from masscubex.noise import estimate_noise


def test_legacy_clipped_background_tracks_mz_and_rejects_uniform_noise():
    rng = np.random.default_rng(76)
    mz = np.linspace(100, 300, 4001)
    y = np.where(mz < 200, 10, 1000) * rng.uniform(0.8, 1.2, mz.size)
    y[500] = 100  # A fixed 500/30000 cutoff would lose this true peak.
    y[3500] = 10000
    _, _, threshold = estimate_noise(
        mz, y, Config(noise_estimator="clipped", noise_window_da=25, noise_min_neighbors=8)
    )
    assert 10 < threshold[500] < 30
    assert 1000 < threshold[3500] < 3000
    assert np.count_nonzero(y > threshold) < 20
    assert y[500] > threshold[500]
    assert y[3500] > threshold[3500]


def test_strong_peak_does_not_mask_nearby_weak_peak():
    mz = np.linspace(100, 125, 501)
    y = np.full(mz.size, 10.0)
    y[250] = 1e9
    y[252] = 100
    _, _, threshold = estimate_noise(mz, y)
    assert threshold[252] < y[252]
    assert threshold[250] < y[250]
    assert np.count_nonzero(y > threshold) == 2


def test_intensity_scale_invariance_and_no_mutation():
    rng = np.random.default_rng(9)
    mz = np.linspace(100, 200, 1001)
    intensity = rng.lognormal(3, 0.2, mz.size)
    original = intensity.copy()
    base = estimate_noise(mz, intensity)[2]
    scaled = estimate_noise(mz, intensity * 1e5)[2]
    np.testing.assert_allclose(scaled, base * 1e5, rtol=1e-10)
    np.testing.assert_array_equal(intensity, original)


def test_physical_window_does_not_bridge_empty_mz_regions():
    mz = np.r_[np.linspace(100, 105, 100), np.linspace(200, 205, 100)]
    y = np.r_[np.full(100, 10.0), np.full(100, 10000.0)]
    threshold = estimate_noise(mz, y)[2]
    np.testing.assert_allclose(threshold[:100], 10)
    np.testing.assert_allclose(threshold[100:], 10000)


def test_isolated_huge_peak_cannot_corrupt_later_windows():
    mz = np.r_[100.0, np.linspace(200, 300, 1001)]
    y = np.r_[1e15, np.tile([9.0, 10.0, 11.0], 334)[:1001]]
    full = estimate_noise(mz, y)[2]
    without_isolated_peak = estimate_noise(mz[1:], y[1:])[2]
    np.testing.assert_allclose(full[1:], without_isolated_peak, rtol=1e-12)


def test_sparse_spectra_empty_and_floor():
    cfg = Config(min_intensity=5)
    np.testing.assert_array_equal(estimate_noise([100, 110], [10, 10000], cfg)[2], [5, 5])
    assert estimate_noise([], [], cfg)[2].size == 0
    dense = np.arange(100, 200, 0.1)
    np.testing.assert_array_equal(estimate_noise(dense, np.ones(len(dense)), cfg)[2], 5)


@pytest.mark.parametrize(
    "kwargs",
    [
        {"noise_window_da": 0},
        {"noise_sigma": np.nan},
        {"noise_sigma": 2},
        {"noise_sigma": 4},
        {"noise_iterations": 1.5},
        {"min_scans": True},
        {"min_scans": 1},
        {"mz_tolerance_da": 0, "mz_tolerance_ppm": 0},
        {"segment": "yes"},
        {"peak_prominence": 2},
        {"noise_mode": "instrument"},
        {"noise_estimator": "mean"},
        {"noise_fraction": 0},
        {"noise_fraction": 1.01},
        {"noise_fraction": float("nan")},
        {"noise_fraction": True},
        {"noise_min_neighbors": 1},
    ],
)
def test_invalid_config_rejected(kwargs):
    with pytest.raises(ValueError):
        Config(**kwargs)


def test_noise_disabled_preserves_explicit_floor():
    cfg = replace(Config(), noise_mode="none", min_intensity=7)
    np.testing.assert_array_equal(estimate_noise([100, 101], [10, 100], cfg)[2], [7, 7])


def test_legacy_clipped_chromatographic_truth_across_three_noise_scales():
    from masscubex.pipeline import extract_features
    from masscubex.readers import Scan

    rng = np.random.default_rng(321)
    mz = np.linspace(100, 400, 3001)
    background = np.where(mz < 200, 10.0, np.where(mz < 300, 100.0, 1000.0))
    scans = []
    for i in range(81):
        y = background * rng.uniform(0.9, 1.1, len(mz))
        for position, amplitude, center in [(500, 80, 20), (1500, 800, 40), (2500, 8000, 60)]:
            y[position] += amplitude * np.exp(-0.5 * ((i - center) / 5) ** 2)
        scans.append(Scan(float(i), mz, y))
    features = sorted(
        extract_features(
            iter(scans),
            Config(noise_estimator="clipped", noise_window_da=25, noise_min_neighbors=8),
        ),
        key=lambda feature: feature.mz,
    )
    assert len(features) == 3
    np.testing.assert_allclose([f.mz for f in features], [150, 250, 350], atol=1e-9)
    np.testing.assert_allclose([f.rt_apex_seconds for f in features], [20, 40, 60], atol=1)
    # A fixed 30,000 threshold erases all three, including the high-background peak.
    assert list(extract_features(iter(scans), Config(noise_mode="none", min_intensity=30000))) == []


def brute_force_low_tail(mz, intensity, config):
    """Independent small-array reference: sort each physical neighborhood."""
    mean, std = np.zeros(len(mz)), np.zeros(len(mz))
    cutoff = np.full(len(mz), config.min_intensity)
    for i, mass in enumerate(mz):
        neighbors = (abs(mz - mass) <= config.noise_window_da / 2) & (intensity > 0)
        neighbors[i] = False
        values = np.sort(intensity[neighbors])
        if len(values) < config.noise_min_neighbors:
            continue
        count = max(config.noise_min_neighbors, int(np.ceil(config.noise_fraction * len(values))))
        subset = values[:count]
        mean[i], std[i] = np.mean(subset), np.std(subset, ddof=0)
        cutoff[i] = max(config.min_intensity, mean[i] + config.noise_sigma * std[i])
    return mean, std, cutoff


@pytest.mark.parametrize("fraction", [0.01, 0.05, 0.10, 1.0])
@pytest.mark.parametrize("window", [4.0, 50.0, 1000.0])
def test_low_tail_agrees_with_independent_sorted_window_reference(fraction, window):
    rng = np.random.default_rng(183)
    mz = 100 + np.cumsum(rng.uniform(0, 0.7, 300))
    intensity = np.round(rng.lognormal(4, 2, len(mz)))  # Ties and zero values.
    intensity[::9] = 0
    cfg = Config(noise_fraction=fraction, noise_window_da=window, min_intensity=3.0)
    expected = brute_force_low_tail(mz, intensity, cfg)
    actual = estimate_noise(mz, intensity, cfg)
    np.testing.assert_allclose(actual, expected, rtol=1e-11, atol=1e-10)


def test_low_tail_minimum_expands_10_percent_to_five_and_excludes_candidate():
    mz = np.linspace(100, 110, 41)
    intensity = np.arange(1, 42, dtype=float)
    mean, std, cutoff = estimate_noise(mz, intensity)
    # Forty neighbors -> ceil(10% * 40) = 4, expanded to 5.
    # Candidate 1 is excluded, so its five background points are 2,...,6.
    assert mean[0] == np.mean(np.arange(2, 7))
    assert std[0] == np.std(np.arange(2, 7))
    assert cutoff[0] == mean[0] + 3 * std[0]
    assert mean[-1] == np.mean(np.arange(1, 6))


def test_low_tail_uses_ceiling_when_fraction_exceeds_minimum():
    mz = np.linspace(100, 110, 52)
    intensity = np.arange(1, 53, dtype=float)
    mean, std, _ = estimate_noise(mz, intensity)
    # ceil(51 * .10) = 6, not floor(...)=5 / the minimum of 5.
    assert mean[-1] == 3.5
    assert std[-1] == np.std(np.arange(1, 7))


def test_profile_zero_padding_is_not_a_background_reference():
    mz = np.linspace(100, 110, 101)
    intensity = np.zeros(101)
    intensity[::10] = np.arange(1, 12)
    actual = estimate_noise(mz, intensity)
    expected = brute_force_low_tail(mz, intensity, Config())
    np.testing.assert_allclose(actual, expected)
    assert np.all(actual[0] > 0)
    assert actual[0][1] == 3.0  # A zero sample excludes no positive neighbor.


def test_exact_minimum_and_inclusive_physical_window_boundary():
    mz = np.array([100, 110, 120, 125, 140, 150], dtype=float)
    intensity = np.arange(1, 7, dtype=float)
    cfg = Config()
    result = estimate_noise(mz, intensity, cfg)
    expected = brute_force_low_tail(mz, intensity, cfg)
    np.testing.assert_allclose(result, expected)
    assert result[0][3] == pytest.approx(3.4)  # Five others, including both ±25 Da endpoints.
    assert np.count_nonzero(result[2]) == 1
    assert np.all(estimate_noise(mz, np.zeros(6), cfg)[2] == 0)
    assert np.all(estimate_noise(mz[:-1], intensity[:-1], cfg)[2] == 0)  # Only four neighbors.


def test_high_baseline_in_later_region_does_not_erase_small_standard_deviation():
    rng = np.random.default_rng(44)
    mz = np.r_[np.linspace(100, 105, 100), np.linspace(200, 205, 100)]
    intensity = np.r_[np.ones(100), 1e10 + rng.uniform(-5, 5, 100)]
    actual = estimate_noise(mz, intensity)
    expected = brute_force_low_tail(mz, intensity, Config())
    np.testing.assert_allclose(actual[0], expected[0], rtol=1e-15)
    np.testing.assert_allclose(actual[1], expected[1], atol=2e-6)
    np.testing.assert_allclose(actual[2], expected[2], rtol=1e-15)


def test_low_tail_does_not_depend_on_legacy_clipping_controls():
    rng = np.random.default_rng(91)
    mz = np.linspace(100, 150, 201)
    intensity = rng.lognormal(3, 1, len(mz))
    np.testing.assert_array_equal(
        estimate_noise(mz, intensity),
        estimate_noise(mz, intensity, Config(noise_clip_sigma=10, noise_iterations=0)),
    )


def test_low_tail_sigma_is_not_a_full_distribution_noise_significance_test():
    # This known limitation is intentional: truncation narrows the estimated
    # sigma and cannot reject most uniformly distributed electronic background.
    rng = np.random.default_rng(76)
    mz = np.linspace(100, 125, 1001)
    intensity = rng.uniform(8, 12, len(mz))
    cutoff = estimate_noise(mz, intensity)[2]
    assert np.count_nonzero(intensity > cutoff) > 0.8 * len(mz)


@pytest.mark.parametrize("fraction", [0.1, 0.5, 1.0])
def test_lazy_noise_moments_follow_repeated_empty_and_refilled_rank_subtrees(fraction):
    # Alternating disjoint intensity ranges empty and refill different parts
    # of the rank tree as the m/z window moves across each acquisition region.
    rng = np.random.default_rng(204)
    mz = np.concatenate([100 * region + np.arange(64) * 0.25 for region in range(1, 7)])
    intensity = np.concatenate(
        [rng.integers(1, 100, 64) * scale for scale in (1.0, 1e8, 1.0, 1e8, 1.0, 1e8)]
    )
    intensity[::7] = 0
    config = Config(noise_fraction=fraction, noise_window_da=8.0)
    actual = estimate_noise(mz, intensity, config)
    expected = brute_force_low_tail(mz, intensity, config)
    np.testing.assert_allclose(actual, expected, rtol=1e-12, atol=1e-12)
