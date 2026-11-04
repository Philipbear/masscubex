"""Exact library-kernel and interval boundary checks for the extraction adapter."""

from itertools import product

import numpy as np
import pytest
from scipy.ndimage import gaussian_filter1d

from masscubex import Config, upstream
from masscubex.summary import _span_state, gaussian_weights, maxima, smooth_padded


@pytest.mark.parametrize("sigma", [0.05, 0.5, 0.6, 1.0, 1.25, 4.0, 20.0])
def test_cached_kernel_is_bitwise_equal_to_public_scipy(sigma):
    for n in (1, 2, 5, 8, 20, 100, 2000):
        values = np.random.default_rng(n).lognormal(8, 3, n).astype(np.float32)
        actual = smooth_padded(values, gaussian_weights(sigma))
        padded = np.zeros(n + 2)
        padded[1:-1] = values
        np.testing.assert_array_equal(actual, gaussian_filter1d(padded, sigma))


@pytest.mark.parametrize("distance", [1, 3, 5])
def test_peak_ties_plateaus_and_prominence_match_original(distance):
    config = Config(peak_distance=distance)
    for values in product(range(3), repeat=6):
        values = np.asarray(values, dtype=np.float32)
        expected = upstream.gaussian_filter_and_prominence_method(values, distance=distance)
        np.testing.assert_array_equal(maxima(values, config), expected)


@pytest.mark.parametrize("values", [[0, 0], [0, 3, 0, 3, 0], [0, 1, np.nan, 2, np.nan, 0]])
def test_span_state_keeps_first_apex_and_interior_zeros(values):
    y = np.array(values, dtype=np.float32)
    left, right, count, apex = _span_state(y, 0, len(y))
    nonzero = np.flatnonzero(y)
    if not len(nonzero):
        assert left == right and count == 0
    else:
        assert (left, right) == (nonzero[0], nonzero[-1] + 1)
        assert count == np.count_nonzero(y[left:right] > 0)
        assert apex == left + np.argmax(y[left:right])


def test_private_scipy_api_fallback_preserves_exact_result(monkeypatch):
    from masscubex import summary

    values = np.array([1, 5, 1, 1, 10, 1], dtype=np.float32)
    monkeypatch.setattr(summary, "_nd_image", None)
    np.testing.assert_array_equal(
        maxima(values, Config()), upstream.gaussian_filter_and_prominence_method(values)
    )


@pytest.mark.parametrize("dtype", [np.float32, np.float64])
def test_native_raw_median_gate_matches_numpy_for_random_and_strided_traces(dtype):
    rng = np.random.default_rng(321)
    for length in [1, 2, 3, 4, 5, 8, 17, 64, 129]:
        for _ in range(40):
            values = rng.lognormal(5, 3, (length, 2)).astype(dtype)[:, 1]
            config = Config(peak_distance=int(rng.integers(1, 8)))
            np.testing.assert_array_equal(
                maxima(values, config), upstream.gaussian_filter_and_prominence_method(
                    values, distance=config.peak_distance
                )
            )


def test_native_median_retains_float32_extremes_and_nan_behavior():
    import warnings

    rng = np.random.default_rng(543)
    cases = [
        np.array([2e38, 3e38, 2e38, 3e38, 2e38, 3e38], dtype=np.float32),
        np.array([1, 4, np.nan, 1, 8, 1], dtype=np.float32),
        np.array([0, 4, 0, 0, 8, 0], dtype=np.float32) * np.float32(1e-40),
    ]
    for _ in range(200):
        # Exercise finite, nonnegative values across the float32 exponent range.
        bits = rng.integers(0, 0x7F800000, 12, dtype=np.uint32)
        cases.append(bits.view(np.float32))
    with warnings.catch_warnings(), np.errstate(over="ignore", invalid="ignore"):
        warnings.simplefilter("ignore", RuntimeWarning)
        for values in cases:
            np.testing.assert_array_equal(
                maxima(values, Config()), upstream.gaussian_filter_and_prominence_method(values)
            )
