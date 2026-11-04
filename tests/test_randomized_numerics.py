"""Deterministic randomized checks against independent numerical references."""

import numpy as np
import pytest

from masscubex import Config
from masscubex.noise import estimate_noise
from masscubex.pipeline import extract_features
from masscubex.readers import Scan
from test_noise import brute_force_low_tail
from test_summary_reuse import reference_segments
import upstream_masscube as original


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
        trace = np.column_stack((
            masses[first:last + 1], y[first:last + 1],
            t[first:last + 1], np.zeros(last - first + 1),
        ))
        # Independent pinned finalization supplies float32 apex mass and
        # integer height/area; the partition above remains a hand-built oracle.
        _, (peak,) = reference_segments(trace, Config(segment=False))
        result.append((first, last, peak))
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
    normalized_masses = masses.astype(np.float32)
    normalized_intensity = y.astype(np.float32)
    for i, (mass, values) in enumerate(zip(masses, y)):
        keep = values > 0
        # In a multi-ion scan, upstream applies its float32 weighted-centroid
        # arithmetic even to isolated groups. A whole singleton scan is exempt.
        signals = original.centroid_signals(
            np.column_stack((mass[keep], values[keep])).astype(np.float32)
        )
        normalized_masses[i, keep] = signals[:, 0]
        normalized_intensity[i, keep] = signals[:, 1]
    for ion in range(6):
        expected = oracle_segments(t, normalized_intensity[:, ion], normalized_masses[:, ion], 2, 3.2)
        actual = sorted(
            [f for f in found if abs(f.mz - (100 + 10 * ion)) < 0.01],
            key=lambda f: f.rt_start_seconds,
        )
        assert len(actual) == len(expected)
        for f, (first, last, peak) in zip(actual, expected):
            assert (f.rt_start_seconds, f.rt_end_seconds, f.scan_count) == (
                peak.rt_seq[0] * 60,
                peak.rt_seq[-1] * 60,
                peak.length,
            )
            assert f.peak_height == peak.peak_height
            assert f.peak_area == peak.peak_area
            assert f.mz == float(peak.mz)
            assert f.rt_apex_seconds == peak.rt * 60
            assert f.rt_start_seconds <= f.rt_apex_seconds <= f.rt_end_seconds
