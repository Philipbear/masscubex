"""Pinned upstream core behavior and the default MassCubeX extensions."""

import ast
import hashlib
import json
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pytest

from masscubex import Config, iter_feature_rows, process_file
from masscubex.features import Feature
from masscubex.ms2 import associate_ms2
from masscubex.pipeline import extract_features
from masscubex.readers import Scan
from masscubex import upstream as runtime_source
from conftest import mzml_spectrum, write_mzml
from test_ms2 import ms2, precursor_mzml
import upstream_masscube as original


def upstream_config(**kwargs):
    """Disable added stages for a strict upstream algorithm comparison."""
    return Config(auto_tune=False, profile_mode="reject", ms2_seeded_features=False, **kwargs)


def test_scientific_source_bodies_are_unchanged():
    source = Path(runtime_source.__file__).read_text()
    expected = json.loads(Path(__file__).with_name("upstream_function_hashes.json").read_text())
    for node in ast.parse(source).body:
        if isinstance(node, (ast.ClassDef, ast.FunctionDef)):
            body = ast.get_source_segment(source, node)
            assert hashlib.sha256(body.encode()).hexdigest() == expected[node.name]["sha256"]


def reference_features(scans):
    data = SimpleNamespace(
        scans=[
            SimpleNamespace(
                time=s.rt / 60,
                signals=original.centroid_signals(
                    np.column_stack((s.mz, s.intensity)).astype(np.float32)
                ),
            )
            for s in scans
        ],
        ms1_idx_arr=np.arange(len(scans)),
        params=SimpleNamespace(mz_tol_ms1=0.01, feature_gap_tol=10),
    )
    peaks = original.detect_features(data)
    for _ in range(2):
        peaks = [part for peak in peaks for part in original.segment_feature(peak)]
    for peak in peaks:
        peak.finalize()
    return sorted(peaks, key=lambda peak: peak.mz)


@pytest.mark.parametrize("seed", range(12))
def test_compiled_tracking_and_final_measurements_match_original(seed):
    rng = np.random.default_rng(seed)
    scans = []
    centers = np.array([100, 100.017, 150, 200, 300], dtype=np.float32)
    for i in range(45):
        mz = centers + rng.uniform(-0.009, 0.009, len(centers)).astype(np.float32)
        intensity = rng.integers(1, 10, len(mz)).astype(np.float32)
        keep = rng.random(len(mz)) > 0.4
        order = np.argsort(mz[keep])
        scans.append(Scan(i * 1.3, mz[keep][order], intensity[keep][order], "positive"))
    actual = list(extract_features(scans, upstream_config(noise_mode="none")))
    actual.sort(key=lambda f: (f.mz, f._parent_mz))
    expected = reference_features(scans)
    assert len(actual) == len(expected)
    for ours, old in zip(actual, expected):
        assert ours.mz == old.mz
        assert ours.peak_height == old.peak_height
        assert ours.peak_area == old.peak_area
        assert ours.scan_count == old.length
        assert ours.gaussian_similarity == old.gaussian_similarity
        np.testing.assert_array_equal(ours.trace[:, 0], old.signals[:, 1])
        np.testing.assert_allclose(ours.trace[:, 1] / 60, old.rt_seq, rtol=0, atol=1e-15)


def test_short_traces_and_gap_rule_are_upstream_defaults():
    cfg = Config()
    assert cfg.min_scans == 1 and cfg.max_gap_scans == 10 and cfg.max_gap_seconds == 0
    assert cfg.mz_tolerance_da == 0.01 and cfg.mz_tolerance_ppm == 0
    assert cfg.auto_tune and cfg.ms2_seeded_features and cfg.profile_mode == "centroid"
    scans = [Scan(i * 100.0, np.array([100.0]), np.array([10.0]), "positive") for i in range(2)]
    (feature,) = extract_features(scans, Config(noise_mode="none"))
    assert feature.scan_count == 2  # No extra elapsed-time cutoff.
    (single,) = extract_features(scans[:1], Config(noise_mode="none"))
    assert single.scan_count == 1 and single.peak_area == 0


def test_ms2_highest_height_and_first_spectrum_tie_match_original():
    peaks = [
        Feature(mass, 5, 0, 10, height, 0, 1, 0, None, "positive")
        for mass, height in [(100, 50), (100.05, 100)]
    ]
    records = [
        ms2(1, rt=8, mz=100.001, fragment_intensity=100),
        ms2(2, rt=5, mz=100.001, fragment_intensity=100),
    ]
    for record in records:
        record.precursor_count = 1
    stats = {}
    associate_ms2(peaks, records, upstream_config(), stats)
    assert peaks[0].ms2_count == 0 and peaks[1].ms2_count == 2
    assert peaks[1].best_ms2_scan_number == 1  # First acquisition wins an equal sum.
    reference = SimpleNamespace(
        features=[
            SimpleNamespace(
                mz=f.mz, rt_seq=[0, 10 / 60], peak_edges=None, peak_height=f.peak_height, ms2_seq=[]
            )
            for f in peaks
        ],
        feature_mz_arr=np.array([f.mz for f in peaks]),
        ms2_idx_arr=np.arange(2),
        scans=[SimpleNamespace(time=r.rt / 60, precursor_mz=r.precursor_mz) for r in records],
    )
    original.allocate_ms2_to_features(reference)
    assert [len(f.ms2_seq) for f in reference.features] == [f.ms2_count for f in peaks]


def test_extensions_are_enabled_by_default_and_tables_agree(tmp_path):
    scans = [mzml_spectrum(i, i, [100], [v]) for i, v in enumerate([1, 5, 10, 5, 1])]
    scans.append(mzml_spectrum(99, 4.1, [70], [1000], level=2, precursor=precursor_mzml(200)))
    source = write_mzml(tmp_path / "extensions.mzML", scans)
    a, b = tmp_path / "upstream.tsv", tmp_path / "default.parquet"
    process_file(source, a, config=upstream_config())
    process_file(source, b)
    baseline, extended = list(iter_feature_rows(a)), list(iter_feature_rows(b))
    assert len(baseline) == 1 and len(extended) == 2
    assert extended[-1]["feature_origin"] == "ms2_seeded"
    meta = json.loads(Path(str(b) + ".json").read_text())
    assert meta["parameter_tuning"]["enabled"]
    assert meta["ms2_unassigned"] == 0
    paired = tmp_path / "default.tsv"
    process_file(source, paired)
    assert extended == list(iter_feature_rows(paired))
    baseline_meta = json.loads(Path(str(a) + ".json").read_text())
    assert not baseline_meta["parameter_tuning"]["enabled"]
    assert baseline_meta["ms2_unassigned"] == 1


def test_profile_centroiding_is_enabled_by_default_and_can_be_disabled(tmp_path):
    mz = np.linspace(99.99, 100.01, 101)
    intensity = 1000 * np.exp(-0.5 * ((mz - 100) / 0.002) ** 2)
    source = write_mzml(
        tmp_path / "profile.mzML", [mzml_spectrum(1, 1, mz, intensity, mode="profile")]
    )
    with pytest.raises(ValueError, match="Profile MS1"):
        process_file(source, tmp_path / "rejected.parquet", config=Config(profile_mode="reject"))
    result = process_file(source, tmp_path / "profile.parquet")
    assert result.features > 0
