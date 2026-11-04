"""Boundary cases found during the pre-publication audit."""

import csv
import json
from pathlib import Path

import numpy as np
import pytest

from masscubex import Config, process_file
from conftest import mzml_spectrum, write_mzml
from test_ms2 import precursor_mzml


def test_output_rounding_does_not_change_raw_inclusive_mass_boundary(tmp_path):
    center = 100.123456789123
    radius = 0.002
    source = write_mzml(
        tmp_path / "boundary.mzML",
        [
            mzml_spectrum(1, 1, [center - radius], [10]),
            mzml_spectrum(2, 1.1, [50], [1], level=2, precursor=precursor_mzml(center)),
        ],
    )
    output = tmp_path / "out.tsv"
    process_file(source, output, config=Config(auto_tune=False, min_intensity=100))
    with output.open() as stream:
        (row,) = csv.DictReader(stream, delimiter="\t")
    # Rounding this seed before measuring would incorrectly exclude the raw
    # peak. Preserve the measurement; round only the exported coordinate.
    assert row["mz"] == row["seed_precursor_mz"] == "100.12346"
    assert center - radius < float(row["mz"]) - float(row["gapfill_mz_radius_da"])
    assert float(row["peak_height"]) == 10
    assert json.loads(row["ms2_scan_numbers"]) == [2]


@pytest.mark.parametrize("level", [1, 2])
def test_nonfinite_quantification_cannot_publish_success(tmp_path, level):
    spectra = [mzml_spectrum(i, i, [100], [1e308 if level == 1 else 10]) for i in range(6)]
    if level == 2:
        spectra.append(
            mzml_spectrum(100, 3.1, [50, 60], [1e308, 1e308], level=2, precursor=precursor_mzml())
        )
    source = write_mzml(tmp_path / "overflow.mzML", spectra)
    output = tmp_path / "out.tsv"
    with np.errstate(over="ignore", invalid="ignore"):
        with pytest.raises(ValueError, match="finite|overflow"):
            process_file(source, output, config=Config(auto_tune=False, noise_mode="none"))
    assert not output.exists()
    assert not Path(str(output) + ".json").exists()
    assert not Path(str(output) + ".lock").exists()
    assert not list(tmp_path.glob("*.partial"))


def test_ambiguous_ms1_polarity_never_supplies_precursor_intensity(tmp_path):
    spectra = []
    for i, polarity in enumerate(["positive", "negative"]):
        spectra.append(mzml_spectrum(i * 2, i, [100], [1000], polarity=None))
        spectra.append(
            mzml_spectrum(
                i * 2 + 1,
                i + 0.1,
                [50],
                [1],
                level=2,
                polarity=polarity,
                precursor=precursor_mzml(intensity=None, ref=f"scan={i * 2}"),
            )
        )
    source = write_mzml(tmp_path / "unknown.mzML", spectra)
    output = tmp_path / "out.tsv"
    process_file(source, output, config=Config(best_ms2="precursor_intensity"))
    with output.open() as stream:
        seeds = [
            r for r in csv.DictReader(stream, delimiter="\t") if r["feature_origin"] == "ms2_seeded"
        ]
    assert len(seeds) == 2
    for row in seeds:
        assert row["ms1_quant_status"] == "ambiguous_polarity"
        assert row["best_ms2_precursor_intensity"] == ""
        assert "best_ms2_selection_method" not in row
    meta = json.loads(Path(str(output) + ".json").read_text())
    assert meta["ms2_assigned"] == 2
    assert meta["best_ms2_apex_fallback_feature_ids"] == [int(row["feature_id"]) for row in seeds]


@pytest.mark.parametrize("dense", [False, True])
def test_refined_queries_release_finished_slices_and_enforce_collection_limit(
    tmp_path, monkeypatch, dense
):
    from masscubex import gapfill, pipeline
    from masscubex.tuning import TuningResult

    config = Config(min_intensity=100, ms2_gapfill_rt_window=3, group_ions=False)
    monkeypatch.setattr(
        pipeline,
        "calibrate_file",
        lambda *args, **kwargs: TuningResult(
            {"positive": config}, {"polarities": {"positive": {"status": "tuned"}}}
        ),
    )
    spectra = []
    for t in range(81):
        mass = np.linspace(99.96, 100.04, 81) if dense else np.array([100.0])
        y = np.full(len(mass), 1 + 4 * np.exp(-0.5 * ((t % 20 - 10) / 2) ** 2))
        spectra.append(mzml_spectrum(t, t, mass, y))
        if t in (10, 30, 50, 70):
            spectra.append(
                mzml_spectrum(
                    100 + t, t + 0.1, [50], [1], level=2, precursor=precursor_mzml(ref=f"scan={t}")
                )
            )
    stage, finished, measurements = {"scans": 0}, [], []
    from masscubex.scan_cache import ScanSource

    original_read, original_refine = ScanSource.ms1, gapfill.refine_measurement

    def read(*args, **kwargs):
        for scan in original_read(*args, **kwargs):
            stage["scans"] += 1
            yield scan

    def refine(f, m, c):
        finished.append(stage["scans"])
        measurements.append(m)
        original_refine(f, m, c)

    monkeypatch.setattr(ScanSource, "ms1", read)
    monkeypatch.setattr(gapfill, "refine_measurement", refine)
    output = tmp_path / "out.tsv"
    process_file(write_mzml(tmp_path / "queries.mzML", spectra), output, config=config)
    with output.open() as stream:
        rows = list(csv.DictReader(stream, delimiter="\t"))
    assert len(rows) == 4 and all(r["ms2_count"] == "1" for r in rows)
    assert finished[0] < 20 and len(finished) == 4
    assert all(not m.samples and m.finalized for m in measurements)
    metadata = json.loads(Path(str(output) + ".json").read_text())
    assert metadata["ms2_gapfill_peak_retained_values"] <= 7
    if dense:
        assert all(
            r["ms1_quant_status"] == "sampling_limit" and r["peak_height"] == "" for r in rows
        )
    else:
        assert all(r["ms1_quant_status"] == "gap_filled" for r in rows)
