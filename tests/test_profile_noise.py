"""Profile-specific boundaries with known signals, padding and sample gaps."""

import numpy as np

from masscubex.config import Config
from masscubex.readers import Scan
from masscubex.profile_noise import detect_profile, profile_thresholds


def test_dense_profile_scale_and_raw_peak_height_are_preserved():
    rng = np.random.default_rng(44)
    mz = np.linspace(100, 101, 5001).astype(np.float32)
    y = (10 + rng.normal(0, 2, len(mz)) + 100 * np.exp(-0.5 * ((mz - 100.5) / 0.001) ** 2)).astype(
        np.float32
    )
    original = y.copy()
    baseline, threshold, status = profile_thresholds(mz, y, Config())
    assert status["profile_specific_points"] == len(y)
    assert 5 < np.median(threshold - baseline) < 7
    masses, heights, _, _ = detect_profile(
        Scan(0, mz, y, centroided=False), Config(), "profile_specific"
    )
    k = int(np.argmin(abs(masses - 100.5)))
    assert abs(masses[k] - 100.5) < 0.0002
    assert heights[k] == np.max(y[abs(mz - 100.5) < 0.003])
    np.testing.assert_array_equal(y, original)


def test_zero_padded_profile_falls_back_without_inventing_zero_noise():
    mz = np.linspace(100, 101, 5001).astype(np.float32)
    y = np.zeros(len(mz), dtype=np.float32)
    peak = abs(mz - 100.5) < 0.004
    y[peak] = 100 * np.exp(-0.5 * ((mz[peak] - 100.5) / 0.001) ** 2)
    scan = Scan(0, mz, y, centroided=False)
    original = detect_profile(scan, Config(), "post_only")
    candidate = detect_profile(scan, Config(), "profile_specific")
    assert candidate[3]["censored_segment_fallback_points"] == len(y)
    for before, after in zip(original[:3], candidate[:3]):
        np.testing.assert_array_equal(before, after)


def test_disconnected_short_segments_do_not_identify_profile_noise():
    mz = np.r_[np.linspace(100, 100.001, 20), np.linspace(150, 150.001, 60)].astype(np.float32)
    y = np.full(len(mz), 10, dtype=np.float32)
    _, _, status = profile_thresholds(mz, y, Config())
    assert status.get("profile_specific_points", 0) == 0
    assert status["short_segment_fallback_points"] == len(mz)
