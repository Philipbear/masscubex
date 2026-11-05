"""Exact search boundaries, evidence isolation and disk-backed extraction."""

from dataclasses import replace

import numpy as np
import pytest

from masscubex.config import Config
from masscubex.consistency_cache import ConsistencyCache
from masscubex.consistency_kernels import histogram_counts
from masscubex.consistency_noise import (
    Parameters,
    evaluate,
    intensity_groups,
    mass_recurrence,
    precursor_protection,
)
from masscubex.mass_tolerance import MassRegion, mass_radius
from masscubex.ms2 import MS2Collector
from masscubex.noise import estimate_noise
from masscubex.readers import MS2Scan, Scan
from masscubex.upstream_adapter import merge_centroids


@pytest.mark.parametrize("dtype", [np.float32, np.float64])
@pytest.mark.parametrize("unique", [False, True])
def test_recurrence_at_region_jumps_matches_independent_strict_windows(dtype, unique):
    config = Config(
        mass_regions=(
            MassRegion(0.0, 200.0, 0.0001, 0.0, 0.01, 0.0, 0.0001, 0.0),
            MassRegion(200.0, 400.0, 0.01, 0.0, 0.01, 0.0, 0.01, 0.0),
            MassRegion(400.0, None, 0.0001, 0.0, 0.01, 0.0, 0.0001, 0.0),
        )
    )
    mz = np.array(
        [199.9998, 199.9999, 200.0, 200.0001, 399.9999, 400.0, 400.0001, 450.0], dtype=dtype
    )
    target = Scan(0.0, mz, np.ones(len(mz)))
    rng = np.random.default_rng(13)
    others = []
    for t in range(4):
        masses = np.sort(
            np.r_[mz.astype(float) + rng.normal(0, 0.0002, len(mz)), 200.005, 399.995]
        ).astype(dtype)
        others.append(Scan(t + 1.0, masses, np.ones(len(masses))))
    expected = []
    for mass in mz.astype(float):
        radius = mass_radius(mass, config)
        count = 0
        for other in others:
            values = other.mz.astype(float)
            nearest = np.min(abs(values - mass)) < radius
            members = np.sum((values > mass - radius) & (values < mass + radius))
            count += bool(nearest and (not unique or members == 1))
        expected.append(count >= 2)
    np.testing.assert_array_equal(mass_recurrence(target, others, config, 2, unique), expected)


@pytest.mark.parametrize("dtype", [np.float32, np.float64])
def test_intensity_partition_boundaries_and_float32_spans_match_reference(dtype):
    rng = np.random.default_rng(9)
    mz = np.r_[rng.uniform(100, 900, 500), np.linspace(1.111111, 51.11111, 20)].astype(dtype)
    values = np.r_[rng.choice([14.0, 50.0, 50.05005, 70.0, 90.0], 500), np.full(20, 200.0)].astype(
        dtype
    )
    parameters = Parameters()
    for published in (False, True):
        order = np.argsort(values, kind="stable")
        sorted_values = values[order]
        expected = np.zeros(len(mz), dtype=bool)
        groups = []
        stop = len(mz)
        minimum = 4 if published else max(8, int(np.ceil(np.log2(len(mz)))))
        while stop:
            high = float(sorted_values[stop - 1])
            low = high * 0.999
            start = np.searchsorted(sorted_values, low)
            members = order[start:stop]
            if len(members) >= minimum and (
                published
                or (np.ptp(mz[members]) >= 50 and len(np.unique(np.floor(mz[members] / 10))) >= 3)
            ):
                expected[members] = True
                groups.append((low, high, len(members)))
            stop = start
        mask, actual = intensity_groups(mz, values, parameters, published_rule=published)
        np.testing.assert_array_equal(mask, expected)
        assert actual == sorted(groups)


def test_histogram_preserves_explicit_edges_and_final_closed_bin():
    edges = np.array([0.0, 0.1, 0.7, 1.0, 2.0])
    values = np.r_[edges, np.nextafter(edges, -np.inf), np.nextafter(edges, np.inf)]
    np.testing.assert_array_equal(histogram_counts(values, edges), np.histogram(values, edges)[0])


def example_scans():
    rng = np.random.default_rng(4)
    for t in range(7):
        mz = np.r_[rng.uniform(100, 900, 100), 950.0]
        order = np.argsort(mz)
        yield Scan(
            float(t),
            mz[order].astype("f4"),
            np.full(len(mz), 50.0, dtype="f4"),
            "positive",
            f"scan={t + 1}",
        )


def test_disk_workspace_matches_in_memory_and_preserves_raw_acquisition_metadata(tmp_path):
    config = Config(noise_estimator="consistency")
    frames = list(example_scans())
    raw_record = MS2Scan(
        0.1,
        "scan=20",
        20,
        20,
        "unknown",
        precursor_mz=950.0,
        precursor_ref="scan=1",
        precursor_count=1,
    )
    records = [raw_record]
    evidence = precursor_protection(frames, records, {"positive": config})
    reference, _, reference_status = evaluate(frames, {"positive": config}, ms2_evidence=evidence)
    collector = MS2Collector(config)
    stream = [frames[0], replace(raw_record), *frames[1:]]
    with ConsistencyCache(tmp_path) as cache:
        cache.build(
            collector.ms1_stream(stream), collector.noise_records, config, {"positive": config}
        )
        assert collector.records[0].polarity == "positive"
        assert collector.noise_records[0].polarity == "unknown"
        assert cache.stats["component_status_counts"] == reference_status
        for i, (old, new) in enumerate(zip(reference, cache.decisions)):
            for name in ("coherent", "candidate_noise", "threshold", "ms2_supported"):
                np.testing.assert_array_equal(getattr(old, name), getattr(new, name))
            assert old.components == new.components
            keep = old.keep()
            expected = merge_centroids(
                frames[i].mz[keep],
                frames[i].intensity[keep],
                old.threshold[keep],
                config.centroid_mz_tolerance_da,
            )
            for a, b in zip(cache.preprocess(cache.frames[i], config), expected):
                np.testing.assert_array_equal(a, b)
    assert not list(tmp_path.glob(".masscubex-noise-*"))


def test_noise_workspace_cleans_up_after_profile_rejection(tmp_path):
    scan = Scan(0.0, np.array([100.0, 101.0, 102.0]), np.array([1.0, 2.0, 1.0]), centroided=False)
    with pytest.raises(ValueError, match="Profile"):
        with ConsistencyCache(tmp_path) as cache:
            cache.build(
                [scan], [], Config(noise_estimator="consistency", profile_mode="reject"), {}
            )
    assert not list(tmp_path.glob(".masscubex-noise-*"))


def test_consistency_requires_multi_scan_context_and_fixed_bins():
    with pytest.raises(ValueError, match="neighboring scans"):
        estimate_noise(
            np.array([100.0, 101.0]), np.array([10.0, 20.0]), Config(noise_estimator="consistency")
        )
    with pytest.raises(ValueError, match="fixed m/z bins"):
        Config(noise_estimator="consistency", noise_layout="sliding")
    config = Config(noise_estimator="consistency", noise_mode="none")
    _, _, threshold = estimate_noise(np.array([100.0, 101.0]), np.array([10.0, 20.0]), config)
    np.testing.assert_array_equal(threshold, [0.0, 0.0])
