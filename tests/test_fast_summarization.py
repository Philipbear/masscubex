"""Differential checks against the independent SciPy implementation."""

from dataclasses import fields
from itertools import product

import numpy as np
import pytest
from scipy.ndimage import gaussian_filter1d
from scipy.signal import find_peaks

from masscubex import Config
from masscubex.features import summarize_trace
from masscubex.chromatograms import gaussian_smooth, local_maxima, select_peaks
from masscubex.reference_chromatograms import summarize_reference_trace


def assert_features_equal(actual, expected):
    assert len(actual) == len(expected)
    for fast, reference in zip(actual, expected):
        for field in fields(fast):
            a, b = getattr(fast, field.name), getattr(reference, field.name)
            if isinstance(b, np.ndarray):
                np.testing.assert_array_equal(a, b)
            elif (
                field.name in {"mz", "peak_area", "intensity", "gaussian_similarity"}
                and b is not None
            ):
                # Compiled sums/dots can differ in their final float64 bits.
                assert a == pytest.approx(b, rel=3e-14, abs=1e-14)
            else:
                assert a == b, field.name


@pytest.mark.parametrize("seed", range(30))
@pytest.mark.parametrize(
    "settings",
    [
        {},
        {"smoothing_seconds": 0.7, "peak_distance_seconds": 1.0},
        {"smoothing_sigma": 0.05, "peak_distance": 1},
        {"smoothing_sigma": 4.0, "peak_distance": 10},
        {"segment": False, "quant_method": "peak_area"},
        {"peak_shape_guard": False},
    ],
)
def test_compiled_summarization_matches_scipy(seed, settings):
    rng = np.random.default_rng(seed)
    n = int(rng.integers(2, 200))
    rt = np.cumsum(rng.uniform(0.1, 2, n))
    y = rng.integers(0, 30, n).astype(float)
    y[: seed % 4] = 0
    y[-(seed % 3 + 1) :] = 0
    trace = np.column_stack((rng.normal(100, 0.0001, n), y, rt, np.ones(n)))
    config = Config(min_scans=2, **settings)
    assert_features_equal(
        list(summarize_trace(trace, config)), list(summarize_reference_trace(trace, config))
    )


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


def test_equal_priority_uses_reference_tie_order():
    y = np.array([0, 2, 0, 2, 0, 2, 0], dtype=float)
    assert select_peaks(y, 3, 0.0)[1]
    trace = np.column_stack((np.full(len(y), 100.0), y, np.arange(len(y)), np.ones(len(y))))
    config = Config(min_scans=2, smoothing_sigma=0.05, peak_distance=3)
    assert_features_equal(
        list(summarize_trace(trace, config)), list(summarize_reference_trace(trace, config))
    )


@pytest.mark.parametrize("layout", ["strided", "fortran", "readonly"])
@pytest.mark.parametrize("group", [False, True])
def test_summary_accepts_noncontiguous_trace_views(layout, group):
    values = np.tile([1.0, 2.0, 1.0], 10)
    trace = np.column_stack(
        (np.full(len(values), 100.0), values, np.arange(len(values)), np.ones(len(values)))
    )
    trace = trace[::2] if layout == "strided" else np.asfortranarray(trace)
    if layout == "readonly":
        trace = np.ascontiguousarray(trace)
        trace.flags.writeable = False
    config = Config(min_scans=2, smoothing_sigma=0.05, peak_distance=3, group_ions=group)
    assert_features_equal(
        list(summarize_trace(trace, config)), list(summarize_reference_trace(trace, config))
    )
