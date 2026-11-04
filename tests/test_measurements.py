"""Specialized measurements must retain unrounded reference results."""

import numpy as np
import pytest

from masscubex import upstream
from masscubex.measurements import (
    gaussian_similarity, integrated_area, noise_score_at_apex, pearson_correlation,
)


@pytest.mark.parametrize("dtype", [np.float32, np.float64])
def test_area_retains_numpy_reduction_and_integer_rounding(dtype):
    rng = np.random.default_rng(651)
    for length in [1, 2, 3, 4, 5, 7, 8, 63, 128, 129, 1024]:
        for _ in range(100):
            values = rng.lognormal(7, 5, length).astype(dtype)
            minutes = (np.cumsum(rng.uniform(0.01, 2, length)) / 60).tolist()
            expected = float(int(np.trapz(y=values, x=minutes) * 60))
            assert integrated_area(minutes, values).hex() == expected.hex()


@pytest.mark.parametrize("dtype", [np.float32, np.float64])
def test_trace_noise_matches_upstream_bits_and_preserves_input(dtype):
    rng = np.random.default_rng(230)
    for length in [1, 2, 4, 5, 6, 63, 128, 129, 1024]:
        for case in range(100):
            values = rng.lognormal(4, 6, length).astype(dtype)
            if case % 7 == 0:
                values[:] = values[0]
            if case % 11 == 0:
                values[rng.random(length) > 0.5] = 0
            before = values.copy()
            expected = float(upstream.calculate_noise_score(values))
            actual = float(noise_score_at_apex(values, int(np.argmax(values))))
            assert actual.hex() == expected.hex()
            np.testing.assert_array_equal(values, before)


@pytest.mark.parametrize("dtype", [np.float32, np.float64])
def test_gaussian_similarity_matches_upstream_bits(dtype):
    rng = np.random.default_rng(931)
    for length in [1, 2, 4, 5, 6, 15, 64, 129, 1024]:
        for case in range(80):
            times = np.cumsum(rng.uniform(0.001, 0.01, length)).tolist()
            values = rng.lognormal(5, 3, length).astype(dtype)
            if case % 7 == 0:
                values[:] = values[0]
            if case % 11 == 0:
                values[rng.random(length) > 0.5] = 0
            expected = float(upstream.calculate_gaussian_similarity(times, values))
            actual = float(gaussian_similarity(times, values, int(np.argmax(values))))
            assert actual.hex() == expected.hex()


@pytest.mark.parametrize("dtype", [np.float32, np.float64])
def test_pearson_matches_numpy_bits_including_constant_vectors(dtype):
    rng = np.random.default_rng(32)
    for length in [3, 4, 5, 12, 63, 128, 1024]:
        for case in range(100):
            x = rng.lognormal(0, 5, length).astype(dtype)
            y = rng.lognormal(0, 5, length).astype(dtype)
            if case % 7 == 0:
                x[:] = x[0]
            if case % 11 == 0:
                y[:] = y[0]
            with np.errstate(invalid="ignore", divide="ignore"):
                expected = float(np.corrcoef(x, y)[0, 1])
                actual = float(pearson_correlation(x, y))
            assert (np.isnan(actual) and np.isnan(expected)) or actual.hex() == expected.hex()
