"""Deterministic randomized checks against independent numerical references."""

import math

import numpy as np
import pytest

from masscubex import Config
from masscubex.noise import estimate_noise
from masscubex.pipeline import extract_features
from masscubex.readers import Scan
from test_noise import brute_force_low_tail


@pytest.mark.parametrize("seed", range(24))
def test_randomized_low_tail_matches_sorted_reference(seed):
    rng = np.random.default_rng(seed)
    n = [0, 1, 4, 5, 6, 49, 50, 51, 80, 256][seed % 10]
    scale = [1e-12, 1.0, 1e12][seed % 3]
    mz = 100 + np.cumsum(rng.uniform(0.001, 3, n))
    intensity = np.round(rng.lognormal(3, 2, n)) * scale
    intensity[::7] = 0
    config = Config(
        noise_window_da=[1.0, 50.0, 1000.0][seed % 3],
        noise_fraction=[0.01, 0.05, 0.10, 0.5, 1.0][seed % 5],
        noise_min_neighbors=[2, 5, 8][seed % 3],
        min_intensity=scale,
    )
    expected = brute_force_low_tail(mz, intensity, config)
    actual = estimate_noise(mz, intensity, config)
    np.testing.assert_allclose(actual, expected, rtol=1e-10, atol=1e-10 * scale)
    assert np.all(np.isfinite(actual))
    assert np.all(actual[2] >= scale)


def oracle_segments(t, y, masses, max_scans, max_seconds):
    """Known isolated-ion truth: partition observation indices, then integrate zeros."""
    observed = np.flatnonzero(y > 0)
    segments = []
    for index in observed:
        if (
            not segments
            or index - segments[-1][-1] - 1 > max_scans
            or t[index] - t[segments[-1][-1]] > max_seconds
        ):
            segments.append([index])
        else:
            segments[-1].append(index)
    result = []
    for indices in segments:
        if len(indices) < 2:
            continue
        first, last = indices[0], indices[-1]
        values = y[first : last + 1]
        area = math.fsum((t[i + 1] - t[i]) * (y[i] + y[i + 1]) / 2 for i in range(first, last))
        center = np.average(masses[indices], weights=y[indices])
        result.append((first, last, max(values), area, center, len(indices)))
    return result


@pytest.mark.parametrize("seed", range(20))
def test_randomized_isolated_traces_preserve_mass_height_area_and_gaps(seed):
    rng = np.random.default_rng(seed + 500)
    t = np.cumsum(rng.uniform(0.1, 1.2, 60))
    t[30:] += 20
    scale = [1e-10, 1.0, 1e10][seed % 3]
    y = rng.lognormal(4, 1, (len(t), 6)) * scale
    y[rng.random(y.shape) < 0.3] = 0
    masses = 100 + 10 * np.arange(6) + rng.uniform(-0.0002, 0.0002, y.shape)
    config = Config(
        auto_tune=False,
        noise_mode="none",
        segment=False,
        group_ions=False,
        min_scans=2,
        max_gap_scans=2,
        max_gap_seconds=3.2,
    )
    scans = [
        Scan(rt, mass[values > 0], values[values > 0], "positive")
        for rt, mass, values in zip(t, masses, y)
    ]
    found = list(extract_features(scans, config))
    for ion in range(6):
        expected = oracle_segments(t, y[:, ion], masses[:, ion], 2, 3.2)
        actual = sorted(
            [f for f in found if abs(f.mz - (100 + 10 * ion)) < 0.01],
            key=lambda f: f.rt_start_seconds,
        )
        assert len(actual) == len(expected)
        for f, (first, last, height, area, center, count) in zip(actual, expected):
            assert (f.rt_start_seconds, f.rt_end_seconds, f.scan_count) == (
                t[first],
                t[last],
                count,
            )
            np.testing.assert_allclose([f.peak_height, f.peak_area], [height, area], rtol=1e-12)
            assert f.mz == pytest.approx(center, abs=1e-11)
            assert f.rt_start_seconds <= f.rt_apex_seconds <= f.rt_end_seconds
