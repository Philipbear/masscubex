from dataclasses import asdict

import numpy as np
import pytest
from scipy.ndimage import gaussian_filter1d
from scipy.signal import find_peaks

from masscubex import Config
from masscubex import features


@pytest.mark.parametrize("seed", range(20))
def test_fast_summarization_is_identical_to_full_scipy_path(seed, monkeypatch):
    rng = np.random.default_rng(seed)
    n = int(rng.integers(8, 200))
    rt = np.cumsum(rng.uniform(0.1, 2, n))
    y = rng.integers(0, 30, n).astype(float)
    trace = np.column_stack((rng.normal(100, 0.0001, n), y, rt, np.ones(n)))
    config = Config(min_scans=2, group_ions=False, smoothing_seconds=0.7, peak_distance_seconds=1.0)
    fast = [asdict(f) for f in features.summarize_trace(trace, config)]
    monkeypatch.setattr(features, "_smooth_padded", lambda y, s: gaussian_filter1d(np.pad(y, 1), s))
    monkeypatch.setattr(features, "_has_multiple_maxima", lambda y: True)
    monkeypatch.setattr(features, "_median_spacing", lambda t: np.median(np.diff(t)))
    reference = [asdict(f) for f in features.summarize_trace(trace, config)]
    assert fast == reference


@pytest.mark.parametrize("sigma", [0.05, 0.5, 1, 1.25, 4, 20])
def test_cached_smoothing_and_plateaus_match_scipy(sigma):
    for y in ([0, 1, 1, 0], [10, 1, 10], [1, 1, 1], [0], [0, 0], [1, 2, 3, 2, 1]):
        expected = gaussian_filter1d(np.pad(np.asarray(y, dtype=float), 1), sigma)
        smooth = features._smooth_padded(y, sigma)
        np.testing.assert_array_equal(smooth, expected)
        assert features._has_multiple_maxima(smooth) == (len(find_peaks(smooth)[0]) > 1)
