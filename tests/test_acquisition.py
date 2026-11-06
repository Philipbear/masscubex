"""Acquisition evidence, conservative exclusions and batch continuation."""

import csv
import gzip
import json
from pathlib import Path
import subprocess
import sys

import pytest

from masscubex.acquisition import DIA_TERMS, inspect_acquisition
from masscubex.config import Config
from masscubex.pipeline import process_file
from masscubex.output import iter_feature_rows
from masscubex.reporting import BatchReport

from conftest import mzml_spectrum, mzxml_scan, write_mzml, write_mzxml


def cv(accession, value="", name=""):
    return f'<cvParam accession="{accession}" name="{name}" value="{value}"/>'


def survey(index=1, polarity="positive"):
    spectrum = mzml_spectrum(index, index, [100], [500], polarity=polarity)
    limits = "<scanWindowList><scanWindow>" + cv("MS:1000501", 100)
    limits += cv("MS:1000500", 1100) + "</scanWindow></scanWindowList>"
    return spectrum.replace("</scan>", limits + "</scan>")


def window(index, center, width, *, polarity="positive"):
    precursor = "<precursorList><precursor><isolationWindow>"
    precursor += cv("MS:1000827", center)
    precursor += cv("MS:1000828", width / 2) + cv("MS:1000829", width / 2)
    precursor += "</isolationWindow></precursor></precursorList>"
    return mzml_spectrum(index, index, [50], [30], level=2, precursor=precursor, polarity=polarity)


def cycles(centers, width, *, repeats=3, polarity="positive"):
    scans = []
    for _ in range(repeats):
        scans.append(survey(len(scans) + 1, polarity))
        scans.extend(
            window(len(scans) + i + 1, mass, width, polarity=polarity)
            for i, mass in enumerate(centers)
        )
    return scans


@pytest.mark.parametrize("accession", DIA_TERMS)
def test_declared_dia_cv_terms_skip_before_decoding(tmp_path, accession):
    path = write_mzml(
        tmp_path / "declared.mzML",
        ["<spectrum><binary>not base64</binary></spectrum>"],
        groups="<fileDescription><fileContent>"
        + cv(accession, name="DIA")
        + "</fileContent></fileDescription>",
    )
    info = inspect_acquisition(path)
    assert info.is_dia and info.confidence == "high"
    assert not info.inspection_complete and not info.binary_arrays_decoded
    assert accession in info.evidence[0]


@pytest.mark.parametrize(
    "mode",
    [
        "SWATH",
        "DIA",
        "AIF",
        "MS^E",
        "MSE",
        "MSᴱ",
        "MS(E)",
        "HDMSE",
        "HDMSᴱ",
        "bbCID",
        "SONAR",
        "diaPASEF",
    ],
)
def test_structured_vendor_acquisition_declarations(tmp_path, mode):
    path = write_mzml(
        tmp_path / "mode.mzML",
        [survey()],
        groups=f'<fileDescription><fileContent><userParam name="acquisition mode" value="{mode}"/>'
        "</fileContent></fileDescription>",
    )
    assert inspect_acquisition(path).is_dia


def test_positive_dia_boolean_flag(tmp_path):
    path = write_mzml(
        tmp_path / "flag.mzML",
        [survey()],
        groups='<fileDescription><fileContent><userParam name="DIA" value="true"/>'
        "</fileContent></fileDescription>",
    )
    assert inspect_acquisition(path).is_dia


@pytest.mark.parametrize(
    "field,value",
    [
        ("acquisition mode", "DDA"),
        ("acquisition mode", "not DIA"),
        ("acquisition mode", "not data-independent acquisition"),
        ("DIA", "false"),
        ("acquisition mode", "off"),
        ("sample description", "DIA reference"),
        ("software name", "DIA-NN"),
    ],
)
def test_unrelated_or_disabled_mode_mentions_are_not_dia(tmp_path, field, value):
    path = write_mzml(
        tmp_path / "SWATH_AIF_DIA.mzML",
        [survey()],
        groups=f'<fileDescription><fileContent><userParam name="{field}" value="{value}"/>'
        "</fileContent></fileDescription>",
    )
    assert not inspect_acquisition(path).is_dia


def test_only_referenced_parameter_groups_label_acquisition(tmp_path):
    groups = '<referenceableParamGroupList><referenceableParamGroup id="unused">'
    groups += cv("MS:1003227", name="full mass range")
    groups += "</referenceableParamGroup></referenceableParamGroupList>"
    path = write_mzml(tmp_path / "unused.mzML", [survey()], groups=groups)
    assert not inspect_acquisition(path).is_dia
    path = write_mzml(
        tmp_path / "used.mzML",
        [survey()],
        groups=groups + '<fileDescription><fileContent><referenceableParamGroupRef ref="unused"/>'
        "</fileContent></fileDescription>",
    )
    assert inspect_acquisition(path).is_dia


@pytest.mark.parametrize("compressed", [False, True])
@pytest.mark.parametrize(
    "width,centers",
    [
        (25, list(range(112, 388, 25))),
        (2, list(range(101, 221, 2))),
        (30, list(range(110, 391, 25))),
    ],
)
def test_repeating_swath_including_narrow_and_overlapping_windows(
    tmp_path, compressed, width, centers
):
    path = tmp_path / ("swath.mzML.gz" if compressed else "swath.mzML")
    info = inspect_acquisition(write_mzml(path, cycles(centers, width)))
    assert info.is_dia and "windowed" in info.dia_families


def test_variable_width_swath_and_ms2_only_export(tmp_path):
    windows = [(100, 120), (120, 140), (140, 170), (170, 210), (210, 260), (260, 340)]
    scans = [
        window(j * len(windows) + i, (a + b) / 2, b - a)
        for j in range(3)
        for i, (a, b) in enumerate(windows)
    ]
    assert inspect_acquisition(write_mzml(tmp_path / "ms2-only.mzML", scans)).is_dia


def test_alternating_staggered_swath_cycles(tmp_path):
    scans = []
    for _ in range(3):
        for shift in (0, 12.5):
            scans.append(survey(len(scans) + 1))
            scans += [
                window(len(scans) + i + 1, mass + shift, 25)
                for i, mass in enumerate(range(112, 263, 25))
            ]
    assert inspect_acquisition(write_mzml(tmp_path / "staggered.mzML", scans)).is_dia


def test_selected_ion_center_when_isolation_target_is_missing(tmp_path):
    scans = cycles(list(range(112, 388, 25)), 25)
    converted = []
    for spectrum in scans:
        for mass in range(112, 388, 25):
            if cv("MS:1000827", mass) in spectrum:
                spectrum = spectrum.replace(cv("MS:1000827", mass), "")
                spectrum = spectrum.replace(
                    "</isolationWindow>",
                    "</isolationWindow><selectedIonList><selectedIon>"
                    + cv("MS:1000744", mass)
                    + "</selectedIon></selectedIonList>",
                )
                break
        converted.append(spectrum)
    scans = converted
    assert inspect_acquisition(write_mzml(tmp_path / "selected.mzML", scans)).is_dia


def test_polarity_switching_cycles_are_evaluated_separately(tmp_path):
    scans = []
    for _ in range(3):
        for polarity in ("positive", "negative"):
            scans += [survey(len(scans) + 1, polarity)]
            scans += [
                window(len(scans) + i + 1, mass, 25, polarity=polarity)
                for i, mass in enumerate(range(112, 388, 25))
            ]
    assert inspect_acquisition(write_mzml(tmp_path / "switching.mzML", scans)).is_dia


def test_late_mixed_acquisition_is_detected(tmp_path):
    # More spectra than a short header sample: do not assume the whole file is DDA.
    scans = [survey(i) for i in range(1, 1002)] + cycles(list(range(112, 388, 25)), 25)
    path = write_mzml(
        tmp_path / "mixed.mzML",
        scans,
        groups="<fileDescription><fileContent>"
        + cv("MS:1003221", name="DDA")
        + "</fileContent></fileDescription>",
    )
    info = inspect_acquisition(path)
    assert info.mode == "mixed_dia" and info.spectra_inspected > 1000


def test_prm_sparse_repeats_and_one_wide_dda_scan_are_not_dia(tmp_path):
    scans = cycles([120, 220, 350, 480, 650, 900], 1)
    scans.append(window(1000, 600, 1000))
    assert not inspect_acquisition(write_mzml(tmp_path / "prm.mzML", scans)).is_dia


def test_missing_windows_do_not_create_false_contiguous_cycles(tmp_path):
    scans = []
    for _ in range(3):
        scans += [survey(len(scans) + 1)]
        scans += [
            window(len(scans) + i + 1, mass, 25) for i, mass in enumerate(range(112, 388, 25))
        ]
        scans.append(mzml_spectrum(len(scans) + 1, len(scans) + 1, [], [], level=2))
    assert not inspect_acquisition(write_mzml(tmp_path / "missing.mzML", scans)).is_dia


def test_repeated_full_range_isolation_is_all_ion(tmp_path):
    scans = [survey()] + [window(i, 600, 1000) for i in range(2, 5)]
    info = inspect_acquisition(write_mzml(tmp_path / "waters.mzML", scans))
    assert info.is_dia and "all_ion" in info.dia_families
    assert "acquired MS1 range" in info.evidence[0]


def test_few_broad_swath_tiles_need_full_survey_coverage(tmp_path):
    full = cycles([350, 850], 500)
    assert inspect_acquisition(write_mzml(tmp_path / "broad-swath.mzML", full)).is_dia
    restricted = cycles([125, 175], 50)
    assert not inspect_acquisition(write_mzml(tmp_path / "broad-prm.mzML", restricted)).is_dia


def test_one_of_multiple_ms1_ranges_is_not_the_full_acquisition(tmp_path):
    first = survey().replace(cv("MS:1000500", 1100), cv("MS:1000500", 400))
    extra = "<scanWindow>" + cv("MS:1000501", 700) + cv("MS:1000500", 1100) + "</scanWindow>"
    first = first.replace("</scanWindowList>", extra + "</scanWindowList>")
    scans = [first] + [window(i, 900, 400) for i in range(2, 5)]
    assert not inspect_acquisition(write_mzml(tmp_path / "multiple-ranges.mzML", scans)).is_dia


def test_wide_isolation_without_a_survey_range_remains_unknown(tmp_path):
    path = write_mzml(tmp_path / "unknown.mzML", [window(i, 600, 1000) for i in range(3)])
    assert inspect_acquisition(path).mode == "unknown"


@pytest.mark.parametrize("kind", ["mzml", "mzxml"])
def test_thermo_all_ion_filter_and_mzxml_scan_type(tmp_path, kind):
    if kind == "mzml":
        spectrum = mzml_spectrum(2, 2, [], [], level=2)
        spectrum = spectrum.replace(
            "<scan>", "<scan>" + cv("MS:1000512", "FTMS + p NSI Full ms2 0.00@hcd30.00 [50-1200]")
        )
        path = write_mzml(tmp_path / "aif.mzML", [spectrum])
    else:
        scan = mzxml_scan(2, 2, [], [], level=2).replace(
            'msLevel="2"', 'msLevel="2" scanType="AIF"'
        )
        path = write_mzxml(tmp_path / "aif.mzXML", [scan])
    assert inspect_acquisition(path).is_dia


def test_nested_mzxml_swath_parent_order(tmp_path):
    scans = []
    for cycle in range(3):
        nested = "".join(
            mzxml_scan(
                100 * cycle + i + 1,
                cycle + 0.1 * i,
                [],
                [],
                level=2,
                precursor=f'<precursorMz windowWideness="25">{mass}</precursorMz>',
            )
            for i, mass in enumerate(range(112, 388, 25))
        )
        scans.append(mzxml_scan(100 * cycle, cycle, [], [], nested=nested))
    assert inspect_acquisition(write_mzxml(tmp_path / "nested.mzXML", scans)).is_dia


def test_ms3_is_not_called_ms1_only_and_bad_xml_is_an_error(tmp_path):
    path = write_mzml(tmp_path / "ms3.mzML", [survey(), mzml_spectrum(2, 2, [], [], level=3)])
    assert inspect_acquisition(path).mode == "unknown"
    path.write_bytes(path.read_bytes()[:-25])
    with pytest.raises(ValueError, match="Acquisition metadata"):
        inspect_acquisition(path)


def test_vendor_identity_from_generic_model_value(tmp_path):
    header = "<instrumentConfigurationList><instrumentConfiguration>"
    header += cv("MS:1000031", "Agilent 6545 Q-TOF", "instrument model")
    header += "</instrumentConfiguration></instrumentConfigurationList>"
    path = write_mzml(tmp_path / "agilent.mzML", [survey()], groups=header)
    assert inspect_acquisition(path).vendor == "Agilent"


def test_chromatogram_only_is_reported_before_calibration(tmp_path, monkeypatch):
    path = write_mzml(
        tmp_path / "mrm.mzML",
        [],
        groups='<run><chromatogramList><chromatogram id="SRM"/></chromatogramList></run>',
    )
    info = inspect_acquisition(path)
    assert info.mode == "chromatogram_only" and info.chromatograms_inspected == 1
    monkeypatch.setattr(
        "masscubex.pipeline.calibrate_file", lambda *a, **k: pytest.fail("calibrated")
    )
    with pytest.raises(ValueError, match="Chromatogram-only export"):
        process_file(path, tmp_path / "mrm.parquet")


@pytest.mark.parametrize("compressed,suffix", [(True, ".mzML"), (False, ".mzML.gz")])
def test_compression_is_detected_from_content_consistently(tmp_path, compressed, suffix):
    from masscubex.readers import iter_ms1
    from masscubex.upstream_adapter import resolve_ms2_filter

    original = write_mzml(
        tmp_path / "original.mzML",
        [survey()],
        groups="<instrumentConfigurationList><instrumentConfiguration>"
        + cv("MS:1002536", name="Orbitrap Fusion")
        + "</instrumentConfiguration></instrumentConfigurationList>",
    )
    path = tmp_path / ("renamed" + suffix)
    data = original.read_bytes()
    path.write_bytes(gzip.compress(data) if compressed else data)
    assert inspect_acquisition(path).vendor == "Thermo"
    assert len(list(iter_ms1(path))) == 1
    assert resolve_ms2_filter(path, Config()).ms2_min_intensity == 10000


@pytest.mark.parametrize(
    "native_id,expected",
    [
        ("scanId=119107", 119107),
        ("scanId=0", 0),
        ("sample=1 period=1 cycle=200 experiment=2", None),
        ("index=119107", None),
        ("scanId=1 scanId=2", None),
    ],
)
def test_vendor_scan_numbers_are_not_confused_with_indices(tmp_path, native_id, expected):
    from masscubex.readers import iter_scans

    spectrum = mzml_spectrum(1, 1, [], [], level=2).replace('id="scan=1"', f'id="{native_id}"')
    path = write_mzml(tmp_path / "vendor.mzML", [spectrum])
    (scan,) = iter_scans(path)
    assert scan.scan_number == expected and scan.native_id == native_id


def test_dia_exclusion_does_not_calibrate_or_create_a_table(tmp_path, monkeypatch):
    path = write_mzml(tmp_path / "dia.mzML", [survey()], groups=cv("MS:1003227", name="AIF"))
    monkeypatch.setattr(
        "masscubex.pipeline.calibrate_file", lambda *a, **k: pytest.fail("DIA calibrated")
    )
    output = tmp_path / "dia.parquet"
    result = process_file(path, output, config=Config(dia_policy="skip"))
    assert result.status == "skipped" and result.skip_reason == "dia_acquisition"
    assert result.output is None and not output.exists()
    assert not list(tmp_path.glob("*.lock")) and not list(tmp_path.glob("*.partial"))
    # Old feature tables must never be silently deleted by the exclusion policy.
    output.write_bytes(b"old output")
    excluded = process_file(path, output, config=Config(dia_policy="skip"))
    assert excluded.skip_reason == "dia_acquisition"
    assert output.read_bytes() == b"old output"
    with pytest.raises(ValueError, match="DIA acquisition"):
        process_file(path, output, config=Config(dia_policy="error"))


def test_dia_report_distinguishes_policy_exclusion_from_resume(tmp_path):
    result = {
        "input": "/input/dia.mzML",
        "status": "skipped",
        "features": 0,
        "skip_reason": "dia_acquisition",
        "acquisition": {"mode": "dia", "vendor": "Waters", "evidence": ["all-ion isolation"]},
    }
    with BatchReport(tmp_path / "reports", Config(), {}) as report:
        report.record(result)
        report.record(
            {
                "input": "/input/done.mzML",
                "status": "skipped",
                "features": 1,
                "skip_reason": "already_complete",
            }
        )
    summary = json.loads(report.summary_path.read_text())
    assert summary["skip_reasons"] == {"dia_acquisition": 1, "already_complete": 1}
    assert summary["acquisition_modes"] == {"dia": 1}
    with (report.directory / "skipped.tsv").open() as stream:
        rows = list(csv.DictReader(stream, delimiter="\t"))
    assert len(rows) == 1 and rows[0]["reason"] == "dia_acquisition"


def test_default_dia_processing_and_resume_provenance(tmp_path):
    path = write_mzml(tmp_path / "dia.mzML", [survey()], groups=cv("MS:1003227", name="AIF"))
    config = Config()
    assert config.dia_policy == "process"
    result = process_file(path, tmp_path / "features.tsv", config=config)
    assert result.status == "ok" and result.acquisition["mode"] == "dia"
    resumed = process_file(path, tmp_path / "features.tsv", config=config, resume=True)
    assert resumed.skip_reason == "already_complete" and resumed.acquisition == result.acquisition
    metadata = json.loads((tmp_path / "features.tsv.json").read_text())
    assert metadata["config"]["associate_ms2"] is True
    assert metadata["effective_config"]["associate_ms2"] is False
    assert metadata["effective_config"]["ms2_seeded_features"] is False
    assert metadata["effective_config"]["group_isf"] is False
    assert metadata["effective_config"]["group_ions"] is True
    assert metadata["acquisition_policy"]["applied"] == "dia_ms1_only"
    assert metadata["algorithm_policy"]["ms2_seeding_extension"] is False
    excluded = process_file(
        path, tmp_path / "features.tsv", config=Config(dia_policy="skip"), resume=True
    )
    assert excluded.skip_reason == "dia_acquisition"


@pytest.mark.parametrize("mixed", [False, True])
def test_dia_extracts_same_ms1_features_without_fragment_links(tmp_path, monkeypatch, mixed):
    ms1_scans = [
        mzml_spectrum(index, index, [100.0, 121.981942], [height, height / 2])
        for index, height in enumerate([1000, 1500, 2000, 1500, 1000], start=1)
    ]
    # A strong fragment/precursor must not add seeded features or ISF labels.
    precursor = (
        '<precursorList><precursor><selectedIonList><selectedIon>'
        + cv("MS:1000744", 300)
        + '</selectedIon></selectedIonList></precursor></precursorList>'
    )
    scans = ms1_scans.copy()
    scans.insert(1, mzml_spectrum(6, 1.5, [100], [100000], level=2, precursor=precursor))
    declaration = (cv("MS:1003221", name="DDA") if mixed else "") + cv(
        "MS:1003227", name="AIF"
    )
    dia = write_mzml(tmp_path / "dia.mzML", scans, groups=declaration)
    # Keep baseline settings identical to the automatically selected DIA policy.
    baseline = write_mzml(tmp_path / "ms1.mzML", ms1_scans)
    expected = process_file(
        baseline,
        tmp_path / "baseline.parquet",
        config=Config(associate_ms2=False, ms2_seeded_features=False, group_isf=False),
    )
    monkeypatch.setattr(
        "masscubex.pipeline.associate_ms2", lambda *a, **k: pytest.fail("DIA fragments associated")
    )
    monkeypatch.setattr(
        "masscubex.gapfill.seed_ms2_features", lambda *a, **k: pytest.fail("DIA seeds created")
    )
    result = process_file(dia, tmp_path / "dia.parquet")
    assert result.acquisition["mode"] == ("mixed_dia" if mixed else "dia")
    assert list(iter_feature_rows(result.output)) == list(iter_feature_rows(expected.output))
    rows = list(iter_feature_rows(result.output))
    assert len(rows) == 2
    assert all(not row["ms2_scan_numbers"] for row in rows)
    assert all(row["best_ms2_scan_number"] is None for row in rows)
    assert any(row["adduct"] for row in rows)
    metadata = json.loads((tmp_path / "dia.parquet.json").read_text())
    assert metadata["acquisition_policy"]["scope"] == "entire_file"
    assert metadata["parameter_tuning"]["enabled"] is True


def test_dia_ms2_only_does_not_create_precursor_seeded_rows(tmp_path):
    path = write_mzml(
        tmp_path / "dia_only.mzML", [window(1, 200, 100)], groups=cv("MS:1003227", name="AIF")
    )
    with pytest.raises(ValueError, match="DIA MS1-only extraction requires MS1"):
        process_file(path, tmp_path / "features.parquet")
    assert not (tmp_path / "features.parquet").exists()
    assert not list(tmp_path.glob("*.partial"))


def test_cli_default_processes_dia_ms1_features(tmp_path):
    path = write_mzml(tmp_path / "dia.mzML", [survey()], groups=cv("MS:1003227", name="AIF"))
    output = tmp_path / "output"
    run = subprocess.run(
        [sys.executable, "-m", "masscubex", "extract", str(path), "-o", str(output)],
        capture_output=True,
        text=True,
    )
    assert run.returncode == 0, run.stderr
    result = json.loads(run.stdout)
    assert result["status"] == "ok"
    assert result["acquisition"]["mode"] == "dia"
    metadata = json.loads(Path(result["output"] + ".json").read_text())
    assert metadata["acquisition_policy"]["applied"] == "dia_ms1_only"


@pytest.mark.parametrize("value", ["invalid", None, True, [], {}])
def test_dia_policy_validation(value):
    with pytest.raises(ValueError, match="dia_policy"):
        Config(dia_policy=value)


@pytest.mark.parametrize("workers", [1, 2])
def test_cli_skips_dia_reports_corruption_and_processes_later_files(tmp_path, workers):
    dia = write_mzml(tmp_path / "dia.mzML", [survey()], groups=cv("MS:1003227", name="AIF"))
    broken = tmp_path / "broken.mzML"
    broken.write_text("<mzML><run>")
    dda = write_mzml(tmp_path / "dda.mzML", [survey()])
    output = tmp_path / "output"
    run = subprocess.run(
        [
            sys.executable,
            "-m",
            "masscubex",
            "extract",
            str(dia),
            str(broken),
            str(dda),
            "-o",
            str(output),
            "--workers",
            str(workers),
            "--dia-policy",
            "skip",
        ],
        capture_output=True,
        text=True,
    )
    assert run.returncode == 1, run.stderr
    rows = [json.loads(line) for line in run.stdout.splitlines()]
    assert {row["status"] for row in rows} == {"skipped", "error", "ok"}
    excluded = next(row for row in rows if row["status"] == "skipped")
    assert excluded["skip_reason"] == "dia_acquisition" and excluded["output"] is None
    (summary_path,) = output.glob("reports/*/summary.json")
    summary = json.loads(summary_path.read_text())
    assert summary["counts"] == {"total": 3, "succeeded": 1, "skipped": 1, "failed": 1}
    assert summary["skip_reasons"] == {"dia_acquisition": 1}
    assert len(list(output.rglob("*.parquet"))) == 1
    assert not list(output.rglob("*.lock")) and not list(output.rglob("*.partial"))
