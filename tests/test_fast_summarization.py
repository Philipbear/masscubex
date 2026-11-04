"""Calibration primitives; full feature parity lives in test_upstream_defaults."""

from itertools import product
import numpy as np
import pytest
from scipy.ndimage import gaussian_filter1d
from scipy.signal import find_peaks
from masscubex.chromatograms import gaussian_smooth, local_maxima, select_peaks


@pytest.mark.parametrize("sigma", [0.05, 0.5, 1, 1.25, 4, 20])
def test_compiled_gaussian_matches_scipy(sigma):
    for y in ([], [0, 1, 1, 0], [10, 1, 10], [1, 1, 1], [0], [0, 0], [1, 2, 3, 2, 1]):
        values = np.asarray(y, dtype=float)
        expected = gaussian_filter1d(np.pad(values, 1), sigma)
        smooth = gaussian_smooth(values, sigma, pad=True)
        np.testing.assert_allclose(smooth, expected, rtol=5e-15, atol=1e-15)
    assert gaussian_smooth(np.empty(0), sigma).size == 0


def test_plateaus_prominence_and_distance_order_exhaustively():
    for values in product(range(3), repeat=6):
        signal = np.asarray(values, dtype=float)
        np.testing.assert_array_equal(local_maxima(signal), find_peaks(signal)[0])
        for distance in (1, 3, 5):
            actual, tied = select_peaks(signal, distance, 1.0)
            if not tied:
                np.testing.assert_array_equal(
                    actual, find_peaks(signal, distance=distance, prominence=1.0)[0]
                )
