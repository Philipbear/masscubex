import numpy as np
import pytest

from masscubex.centroid import centroid_profile
from masscubex.config import Config
from masscubex.noise import estimate_noise
from masscubex.pipeline import extract_features
from masscubex.readers import Scan


def detect(mz, y, config=None):
    baseline, _, threshold = estimate_noise(mz, y, config)
    return centroid_profile(mz, y, baseline, threshold, 3, 4.0)


@pytest.mark.parametrize("center", [100.00013, 500.123456, 1000.99871])
def test_gaussian_mass_accuracy_on_nonuniform_sampling(center):
    rng = np.random.default_rng(7)
    offsets = np.linspace(-0.05, 0.05, 501)
    offsets[1:-1] += rng.uniform(-0.00002, 0.00002, 499)
    mz = center + offsets
    y = 10 + 10000 * np.exp(-0.5 * (offsets / 0.002) ** 2)
    masses, heights, _ = detect(mz, y)
    assert len(masses) == 1
    assert abs(masses[0] - center) / center * 1e6 < 0.1
    assert heights[0] > 9900


def test_profile_resolves_neighbors_and_rejects_single_bin_spike():
    mz = np.linspace(100, 100.2, 1001)
    y = 5 + 1000 * np.exp(-0.5 * ((mz - 100.05) / 0.002) ** 2)
    y += 500 * np.exp(-0.5 * ((mz - 100.10) / 0.002) ** 2)
    y[850] = 1e5
    masses, _, _ = detect(mz, y)
    np.testing.assert_allclose(masses, [100.05, 100.10], atol=2e-5)


def test_flat_top_empty_and_sparse_profiles():
    mz = np.arange(100, 107, dtype=float)
    y = np.array([0, 1, 5, 5, 5, 1, 0.0])
    masses, heights, _ = detect(mz, y, Config(noise_mode="none"))
    np.testing.assert_array_equal(masses, [103])
    np.testing.assert_array_equal(heights, [5])
    assert len(detect(np.array([]), np.array([]))[0]) == 0
    assert len(detect(np.array([100.0]), np.array([100.0]))[0]) == 0


def test_centroid_does_not_bridge_large_sampling_gap():
    mz = np.array([100, 100.001, 100.002, 110, 110.001])
    y = np.array([0, 10, 20, 10, 0.0])
    assert len(detect(mz, y, Config(noise_mode="none"))[0]) == 0


def test_hysteresis_keeps_weak_peak_shoulders_with_same_apex_significance():
    mz = np.arange(100, 105, dtype=float)
    y = np.array([0.0, 4.0, 9.0, 4.0, 0.0])
    high = np.full(5, 8.0)
    low = np.full(5, 3.0)
    baseline = np.zeros(5)
    assert not len(centroid_profile(mz, y, baseline, high)[0])
    masses, heights, thresholds = centroid_profile(mz, y, baseline, low, apex_threshold=high)
    np.testing.assert_array_equal(masses, [102.0])
    np.testing.assert_array_equal(heights, [9.0])
    np.testing.assert_array_equal(thresholds, [8.0])
    y[2] = 7.0
    assert not len(centroid_profile(mz, y, baseline, low, apex_threshold=high)[0])


def test_profile_to_chromatographic_feature():
    mz = np.linspace(99.95, 100.05, 501)
    scans = [
        Scan(
            float(i),
            mz,
            10 + amplitude * np.exp(-0.5 * ((mz - 100) / 0.002) ** 2),
            centroided=False,
        )
        for i, amplitude in enumerate([100, 200, 400, 200, 100])
    ]
    stats = {}
    (feature,) = extract_features(scans, Config(), stats)
    assert feature.mz == pytest.approx(100, abs=1e-8)
    assert feature.rt_apex_seconds == 2
    assert feature.peak_height == pytest.approx(410)
    assert stats["profile_scans"] == 5
    with pytest.raises(ValueError, match="Profile MS1"):
        list(extract_features(scans, Config(profile_mode="reject")))
