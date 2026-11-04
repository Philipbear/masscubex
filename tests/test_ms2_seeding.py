"""MS2-to-feature completeness and independently calculated forced MS1 measurements."""

import csv
import json
from dataclasses import replace
from pathlib import Path

import pytest

from masscubex import Config, process_file
from conftest import mzml_spectrum, mzxml_scan, write_mzml, write_mzxml
from test_ms2 import precursor_mzml


def outputs(path, output, config=None):
    process_file(path, output, config=config)
    with output.open() as stream:
        rows = list(csv.DictReader(stream, delimiter="\t"))
    metadata = json.loads(Path(str(output) + ".json").read_text())
    # All successful-file MS2 scans must occur once in the feature table or in
    # the explicit unresolved list. Numeric scan numbers alone are not unique.
    indices = [i for row in rows for i in json.loads(row["ms2_spectrum_indices"])]
    indices += [r["spectrum_index"] for r in metadata["unassigned_ms2"]]
    assert len(indices) == len(set(indices)) == metadata["ms2_scans"]
    for row in rows:
        if int(row["ms2_count"]):
            assert int(row["best_ms2_spectrum_index"]) in json.loads(row["ms2_spectrum_indices"])
    return rows, metadata


@pytest.mark.parametrize("format_", ["mzml", "mzxml"])
@pytest.mark.parametrize("profile", [False, True])
def test_weak_ms1_gapfill_all_ids_widths_and_quantification(tmp_path, format_, profile):
    # All MS1 samples are below the detection floor. The local raw EIC has a
    # known height of 4 and trapezoidal area of 15 intensity*seconds.
    scan_factory = mzml_spectrum if format_ == "mzml" else mzxml_scan
    spectra = []
    for i, height in enumerate([1, 2, 3, 4, 3, 2, 1]):
        options = (
            {"mode": "profile" if profile else "centroid"}
            if format_ == "mzml"
            else {"centroid": "0" if profile else "1"}
        )
        mz = [99.999, 100, 100.001] if profile else [100]
        y = [height / 2, height, height / 2] if profile else [height]
        spectra.append(scan_factory(i * 10 + 1, i, mz, y, **options))
        if i in (2, 3):
            precursor = (
                precursor_mzml(width=25 if i == 2 else 400)
                if format_ == "mzml"
                else f'<precursorMz windowWideness="{25 if i == 2 else 400}">100</precursorMz>'
            )
            spectra.append(
                scan_factory(100 + i, i + 0.1, [50], [i * 10], level=2, precursor=precursor)
            )
    writer = write_mzml if format_ == "mzml" else write_mzxml
    path = writer(tmp_path / ("run.mzML.gz" if format_ == "mzml" else "run.mzXML"), spectra)
    config = Config(noise_mode="none", min_intensity=100, quant_method="peak_area")
    (row,), metadata = outputs(path, tmp_path / "out.tsv", config)
    assert row["feature_origin"] == "ms2_seeded"
    assert row["ms1_quant_status"] == "gap_filled"
    assert float(row["peak_height"]) == 4
    assert float(row["peak_area"]) == 15
    assert row["scan_count"] == "7"
    assert row["apex_noise_threshold"] == ""
    assert json.loads(row["ms2_scan_numbers"]) == [102, 103]
    assert row["best_ms2_scan_number"] == "103"
    assert metadata["ms2_matched_detected"] == 0
    assert metadata["ms2_seeded_features"] == 1
    assert metadata["ms2_assigned"] == 2 and metadata["ms2_unassigned"] == 0
    assert metadata["ms2_skipped_wide_isolation"] == 0
    assert metadata["ms2_gapfill_passes"] == 1


def test_short_ms1_trace_is_preserved_when_below_min_scan_count(tmp_path):
    spectra = [mzml_spectrum(1, 0, [100], [10]), mzml_spectrum(2, 1, [100], [20])]
    spectra.append(mzml_spectrum(3, 1.1, [50], [5], level=2, precursor=precursor_mzml()))
    path = write_mzml(tmp_path / "short.mzML", spectra)
    (row,), _ = outputs(path, tmp_path / "out.tsv")
    assert row["feature_origin"] == "ms2_seeded" and row["scan_count"] == "2"
    assert float(row["peak_height"]) == 20 and float(row["peak_area"]) == 15


def test_absent_ms1_without_supported_noise_and_unavailable_ms1_are_blank(tmp_path):
    spectra = [mzml_spectrum(i + 1, i, [200], [10]) for i in range(3)]
    spectra += [
        mzml_spectrum(10, 1, [50], [5], level=2, precursor=precursor_mzml()),
        mzml_spectrum(11, 30, [50], [5], level=2, precursor=precursor_mzml()),
    ]
    path = write_mzml(tmp_path / "absent.mzML", spectra)
    rows, metadata = outputs(path, tmp_path / "out.tsv", Config(group_ions=True))
    assert len(rows) == 2
    assert rows[0]["ms1_quant_status"] == "no_ms1_signal"
    assert rows[0]["peak_height"] == rows[0]["peak_area"] == ""
    assert rows[1]["ms1_quant_status"] == "no_ms1_scans"
    assert rows[1]["peak_height"] == rows[1]["peak_area"] == ""
    assert metadata["ms2_assigned"] == 2


def test_mixed_polarities_never_share_ms1_evidence(tmp_path):
    spectra = []
    for i in range(3):
        spectra += [
            mzml_spectrum(i * 10, i, [100], [10], polarity="positive"),
            mzml_spectrum(i * 10 + 1, i + 0.01, [100], [1000], polarity="negative"),
        ]
    spectra += [
        mzml_spectrum(
            100, 1.1, [50], [5], level=2, polarity="positive", precursor=precursor_mzml()
        ),
        mzml_spectrum(
            101, 1.2, [50], [5], level=2, polarity="negative", precursor=precursor_mzml()
        ),
        mzml_spectrum(102, 1.3, [50], [5], level=2, polarity=None, precursor=precursor_mzml()),
    ]
    path = write_mzml(tmp_path / "polarity.mzML", spectra)
    rows, metadata = outputs(path, tmp_path / "out.tsv")
    assert len(rows) == 3
    assert [r["peak_height"] for r in rows] == ["10.00", "1000.00", ""]
    assert rows[2]["ms1_quant_status"] == "ambiguous_polarity"
    assert metadata["ms2_assigned"] == 3


def test_unusable_precursors_are_explicit_and_seeding_can_be_disabled(tmp_path):
    precursor = precursor_mzml()
    multiple = precursor.replace("</precursorList>", precursor.split("<precursorList>")[1])
    spectra = [mzml_spectrum(1, 0, [200], [1])]
    spectra += [
        mzml_spectrum(10, 0.1, [50], [5], level=2, precursor=precursor),
        mzml_spectrum(11, 0.2, [50], [5], level=2),
        mzml_spectrum(12, 0.3, [50], [5], level=2, precursor=multiple),
        mzml_spectrum(13, 0.4, [50], [5], level=2, precursor=precursor).replace(
            '<cvParam accession="MS:1000016" value="0.4" unitAccession="UO:0000010"/>', ""
        ),
    ]
    path = write_mzml(tmp_path / "metadata.mzML", spectra)
    rows, metadata = outputs(path, tmp_path / "out.tsv")
    assert len(rows) == 1
    assert metadata["ms2_assigned"] == 1
    assert [r["reason"] for r in metadata["unassigned_ms2"]] == [
        "missing_precursor_mz",
        "multiple_precursors",
        "missing_retention_time",
    ]
    rows, metadata = outputs(path, tmp_path / "off.tsv", Config(ms2_seeded_features=False))
    assert not rows and metadata["ms2_unassigned"] == 4
    assert metadata["ms2_gapfill_passes"] == 0


def test_detected_feature_assignment_is_unchanged_and_has_no_second_pass(tmp_path):
    spectra = [mzml_spectrum(i, i, [100], [v]) for i, v in enumerate([1, 2, 3, 4, 3, 2, 1])]
    spectra.append(mzml_spectrum(100, 3.1, [50], [10], level=2, precursor=precursor_mzml(width=25)))
    path = write_mzml(tmp_path / "detected.mzML", spectra)
    (row,), metadata = outputs(path, tmp_path / "on.tsv")
    (original,), _ = outputs(path, tmp_path / "off.tsv", Config(ms2_seeded_features=False))
    assert row == original and row["feature_origin"] == "ms1_detected"
    assert metadata["ms2_gapfill_passes"] == 0


def test_old_config_width_can_be_overridden_explicitly(tmp_path):
    spectra = [mzml_spectrum(1, 0, [100], [10])]
    spectra.append(mzml_spectrum(2, 0.1, [50], [5], level=2, precursor=precursor_mzml(width=25)))
    path = write_mzml(tmp_path / "wide.mzML", spectra)
    old_config = Config(ms2_max_isolation_width=4)
    rows, metadata = outputs(path, tmp_path / "guard.tsv", old_config)
    assert not rows and metadata["ms2_unassigned_reasons"] == {"isolation_width_limit": 1}
    rows, metadata = outputs(
        path, tmp_path / "default.tsv", replace(old_config, ms2_max_isolation_width=0)
    )
    assert len(rows) == 1 and metadata["ms2_assigned"] == 1


def test_repeated_seeds_are_bounded_and_native_numbers_are_not_used_as_keys(tmp_path):
    spectra = [mzml_spectrum(i + 1, i, [100], [1]) for i in range(21)]
    for i, rt in enumerate([1.0, 2.0, 6.0, 7.0, 20.0]):
        # Different native IDs can contain the same numeric scan number.
        spectrum = mzml_spectrum(100, rt, [50], [i + 1], level=2, precursor=precursor_mzml())
        spectra.append(spectrum.replace('id="scan=100"', f'id="controllerNumber={i} scan=100"'))
    # Include a native ID with no numeric scan number as well.
    spectra.append(
        mzml_spectrum(101, 2.1, [50], [99], level=2, precursor=precursor_mzml()).replace(
            'id="scan=101"',
            'id="cycle=1 experiment=2"',
        )
    )
    path = write_mzml(tmp_path / "repeated.mzML", spectra)
    rows, metadata = outputs(path, tmp_path / "out.tsv", Config(min_intensity=10))
    assert len(rows) == 3  # A seed at 1 s covers [0, 6]; its 6 s member cannot expand it.
    assert [int(r["ms2_count"]) for r in rows] == [4, 1, 1]
    assert json.loads(rows[0]["ms2_scan_numbers"]) == [100, 100, 100, None]
    assert rows[0]["best_ms2_scan_number"] == ""
    assert rows[0]["best_ms2_native_id"] == "cycle=1 experiment=2"
    assert metadata["ms2_assigned"] == 6


def test_gapfill_uses_actual_times_and_does_not_integrate_across_long_acquisition_gaps(tmp_path):
    spectra = [mzml_spectrum(1, 0, [100], [2]), mzml_spectrum(2, 60, [100], [2])]
    spectra.append(mzml_spectrum(3, 30, [50], [5], level=2, precursor=precursor_mzml()))
    path = write_mzml(tmp_path / "gaps.mzML", spectra)
    (row,), _ = outputs(path, tmp_path / "out.tsv", Config(ms2_gapfill_rt_window=30))
    assert row["peak_height"] == "2.00" and row["peak_area"] == "0.00"


def test_gapfill_window_controls_actual_measurement_and_resume(tmp_path):
    spectra = [mzml_spectrum(i, i, [100], [v]) for i, v in enumerate([1, 2, 3, 4, 3, 2, 1])]
    spectra.append(mzml_spectrum(100, 3, [50], [5], level=2, precursor=precursor_mzml()))
    path = write_mzml(tmp_path / "window.mzML", spectra)
    output = tmp_path / "out.tsv"
    config = Config(min_intensity=10, ms2_gapfill_rt_window=1)
    (row,), _ = outputs(path, output, config)
    assert float(row["peak_area"]) == 7
    assert process_file(path, output, config=config, resume=True).status == "skipped"
    assert (
        process_file(
            path, output, config=replace(config, ms2_gapfill_rt_window=5), resume=True
        ).status
        == "ok"
    )
    with output.open() as stream:
        (row,) = csv.DictReader(stream, delimiter="\t")
    assert float(row["peak_area"]) == 15


def test_gapfill_failure_does_not_publish_partial_outputs(tmp_path, monkeypatch):

    path = write_mzml(
        tmp_path / "failure.mzML",
        [
            mzml_spectrum(1, 0, [100], [2]),
            mzml_spectrum(2, 0.1, [50], [5], level=2, precursor=precursor_mzml()),
        ],
    )
    output = tmp_path / "out.tsv"
    outputs(path, output)
    before = output.read_bytes()
    metadata = Path(str(output) + ".json")
    before_json = metadata.read_bytes()

    def broken(*args, **kwargs):
        raise ValueError("second-pass read failed")
        yield  # A streaming reader, not a pre-materialized list.

    from masscubex.scan_cache import ScanSource

    monkeypatch.setattr(ScanSource, "ms1", broken)
    with pytest.raises(ValueError, match="second-pass"):
        process_file(path, output, overwrite=True)
    assert output.read_bytes() == before and metadata.read_bytes() == before_json
    assert not list(tmp_path.glob("*.partial")) and not list(tmp_path.glob("*.lock"))


def test_float32_mass_boundary_and_duplicate_centroids_use_normalized_raw_ms1(tmp_path):
    precursor = "<precursorMz>327.192180998545</precursorMz>"
    # Float32 subtraction would incorrectly place this peak outside 0.01 Da.
    # The reader converts to float64 and sums exact duplicate m/z entries first.
    spectra = [mzxml_scan(1, 0, [327.2021789550781] * 2, [10, 20], precision=32)]
    spectra.append(mzxml_scan(2, 0.1, [50], [5], level=2, precursor=precursor))
    path = write_mzxml(tmp_path / "precision.mzXML", spectra)
    (row,), _ = outputs(
        path, tmp_path / "out.tsv", Config(mz_tolerance_da=0.01, gapfill_mz_tolerance_da=0.01)
    )
    assert row["peak_height"] == "30.00" and row["peak_area"] == "0.00"


@pytest.mark.parametrize("executor", ["thread", "process"])
@pytest.mark.parametrize("output_format", ["parquet", "tsv"])
def test_seeded_outputs_are_identical_across_file_workers(tmp_path, executor, output_format):
    from masscubex.cli import run_batch

    spectra = [mzml_spectrum(1, 0, [100], [2]), mzml_spectrum(2, 1, [100], [4])]
    spectra.append(mzml_spectrum(3, 1.1, [50], [5], level=2, precursor=precursor_mzml(width=25)))
    a = write_mzml(tmp_path / "a.mzML", spectra)
    b = write_mzml(tmp_path / "b.mzML", spectra)
    reference = tmp_path / ("reference." + output_format)
    process_file(a, reference)
    results = list(
        run_batch(
            [a, b],
            tmp_path / "parallel",
            Config(),
            workers=2,
            executor=executor,
            output_format=output_format,
        )
    )
    assert len(results) == 2 and all(r["status"] == "ok" for r in results)
    for result in results:
        assert Path(result["output"]).read_bytes() == reference.read_bytes()
        metadata = json.loads(Path(result["output"] + ".json").read_text())
        assert metadata["ms2_seeded_features"] == metadata["ms2_assigned"] == 1


@pytest.mark.parametrize("polarity", ["positive", None])
def test_ms2_only_files_preserve_features_without_inventing_ms1_quantification(tmp_path, polarity):
    spectra = [
        mzml_spectrum(
            i,
            i,
            [50],
            [i + 1],
            level=2,
            polarity=polarity,
            precursor=precursor_mzml(),
        )
        for i in range(3)
    ]
    path = write_mzml(tmp_path / "only_ms2.mzML", spectra)
    output = tmp_path / "out.tsv"
    (row,), metadata = outputs(path, output, Config(group_ions=True))
    assert row["ms1_quant_status"] == "no_ms1_scans"
    assert row["peak_height"] == row["peak_area"] == ""
    assert json.loads(row["ms2_scan_numbers"]) == [0, 1, 2]
    assert metadata["ms1_scans"] == 0 and metadata["ms2_assigned"] == 3
    assert metadata["ms2_gapfill_passes"] == 0
    assert (
        process_file(path, output, config=Config(group_ions=True), resume=True).status == "skipped"
    )
