"""Regional pilot selection must preserve weak regions and complete mass bounds."""

from dataclasses import replace

import numpy as np
import pytest

from conftest import mzml_spectrum, write_mzml
from masscubex import Config
from masscubex.pilot_sampling import strongest_per_region
from masscubex.readers import MS2Scan
from masscubex.scan_cache import ScanSource
from masscubex.tuning import calibrate_file, competing_peaks
from test_ms2 import precursor_mzml


@pytest.mark.parametrize("seed", range(12))
def test_partial_selection_matches_independent_full_sort(seed):
    rng = np.random.default_rng(seed)
    masses = np.sort(np.r_[rng.uniform(1, 1200, 2000), 200, 500])
    values = rng.integers(1, 15, len(masses)).astype(float)
    breaks = np.array([200.0, 500.0])
    selected, counts = strongest_per_region(masses, values, breaks, 128)
    region = np.searchsorted(breaks, masses, side="right")
    expected = []
    for index in range(3):
        candidates = [i for i in range(len(masses)) if region[i] == index]
        expected.extend(sorted(candidates, key=lambda i: (-values[i], i))[:128])
    np.testing.assert_array_equal(selected, sorted(expected))
    np.testing.assert_array_equal(counts, np.bincount(region[selected], minlength=3))


def test_weak_regions_ties_boundaries_and_unused_capacity():
    masses = np.r_[np.arange(1, 190), 200, 201, 500, 501]
    values = np.r_[np.full(189, 1e9), np.ones(4)]
    selected, counts = strongest_per_region(masses, values, [200, 500], 128)
    np.testing.assert_array_equal(counts, [128, 2, 2])
    np.testing.assert_array_equal(selected, np.r_[np.arange(128), np.arange(189, 193)])
    selected, counts = strongest_per_region(np.array([]), np.array([]), [200, 500], 128)
    assert selected.dtype.kind == "i" and selected.size == 0
    np.testing.assert_array_equal(counts, [0, 0, 0])


def test_competitor_across_boundary_is_checked_even_if_not_selected():
    masses = np.r_[199.99, 200.01, np.arange(201, 333)]
    values = np.r_[100.0, 25.0, np.full(132, 1000.0)]
    selected, _ = strongest_per_region(masses, values, [200, 500], 128)
    assert 0 in selected and 1 not in selected
    flags = competing_peaks(masses, values, selected, 100.0, 0.05)
    assert flags[0] == 1  # The unselected 200.01 ion still invalidates 199.99.


@pytest.mark.parametrize("budget", [0, 0.001, 64])
def test_equal_width_sampling_uses_complete_range_separately_by_polarity(tmp_path, budget):
    spectra = [
        mzml_spectrum(0, 0, [100, 200], [100, 10]),
        mzml_spectrum(1, 0.1, [50, 100], [10, 100], polarity="negative"),
        mzml_spectrum(2, 1, [800, 1000], [100, 10]),
        mzml_spectrum(3, 1.1, [400, 650], [10, 100], polarity="negative"),
        mzml_spectrum(4, 1.2, [50, 60], [5, 10], level=2, precursor=precursor_mzml()),
    ]
    path = write_mzml(tmp_path / "expanding.mzML", spectra)
    config = Config(tuning_mz_regions=3, scan_cache_mb=budget)
    source = ScanSource(path, config)
    result = calibrate_file(path, config, scan_source=source)
    expected = {
        "positive": [(100, 400), (400, 700), (700, 1000)],
        "negative": [(50, 250), (250, 450), (450, 650)],
    }
    for polarity, edges in expected.items():
        info = result.diagnostics["polarities"][polarity]
        assert [(r["lower"], r["upper"]) for r in info["regions"]] == edges
        assert [(r["lower"], r["upper"]) for r in info["pilot_sampling_regions"]] == edges
        assert sum(r["selected_centroids"] for r in info["pilot_sampling_regions"]) == 4
    assert result.diagnostics["range_discovery_pass"]
    assert source.stats["reader_passes"] == (1 if budget == 64 else 2)
    replayed_ms2 = [scan for scan in source.scans() if isinstance(scan, MS2Scan)]
    assert len(replayed_ms2) == 1
    np.testing.assert_array_equal(replayed_ms2[0].mz, [50, 60])
    np.testing.assert_array_equal(replayed_ms2[0].intensity, [5, 10])


def test_fixed_layout_does_not_add_a_range_discovery_pass(tmp_path):
    path = write_mzml(tmp_path / "fixed.mzML", [mzml_spectrum(0, 0, [100, 800], [50, 10])])
    config = Config()
    source = ScanSource(path, config)
    result = calibrate_file(path, config, scan_source=source)
    assert not result.diagnostics["range_discovery_pass"]
    assert source.stats["reader_passes"] == 1 and source.stats["replay_passes"] == 0
    assert result.diagnostics["minimum_gaussian_similarity"] == 0.8
    assert result.diagnostics["minimum_representatives"] == 20
    single = calibrate_file(path, replace(config, tuning_mz_breaks=()))
    assert len(single.diagnostics["polarities"]["positive"]["pilot_sampling_regions"]) == 1


def test_candidate_pool_cannot_be_filled_by_another_region(tmp_path):
    masses = np.r_[20 + 0.2 * np.arange(700), 600.0, 610.0]
    height = np.r_[1e8 + np.arange(700), 1000, 900]
    spectra = [
        mzml_spectrum(i, float(i), masses, height * np.exp(-0.5 * ((i - 20) / 3) ** 2))
        for i in range(41)
    ]
    path = write_mzml(tmp_path / "crowded.mzML", spectra)
    config = Config(tuning_peak_limit=1024, tuning_max_peaks=20)
    result = calibrate_file(path, config)
    info = result.diagnostics["polarities"]["positive"]
    pools = info["pilot_sampling_regions"]
    assert [r["qualified_candidates"] for r in pools] == [700, 0, 2]
    assert [r["retained_candidates"] for r in pools] == [200, 0, 2]
    high = [r for r in info["regional_representatives"] if r["mz"] >= 500]
    assert len(high) == 2 and all(r["gaussian_similarity"] >= 0.9 for r in high)
