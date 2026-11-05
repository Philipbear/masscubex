"""Co-exported optical spectra must never become mass features or mask bad MS."""

import json
from pathlib import Path

import pytest

from masscubex import iter_feature_rows, process_file
from masscubex.acquisition import inspect_acquisition
from masscubex.readers import iter_scans
from masscubex.spectrum_types import NON_MS_SPECTRUM_TERMS
from conftest import mzml_spectrum, write_mzml


def _cv(accession, **attrs):
    extra = " ".join(f'{key}="{value}"' for key, value in attrs.items())
    return f'<cvParam accession="{accession}" {extra}/>'


def _optical(index=1, *, accession="MS:1000804", referenced=False):
    designation = (
        '<referenceableParamGroupRef ref="optical"/>' if referenced else _cv(accession)
    )
    # Deliberately unusable as mass data: no MS level, centroid flag or time,
    # with invalid base64 arrays that must not be decoded by an MS reader.
    return (
        f'<spectrum id="optical={index}" defaultArrayLength="261">{designation}'
        + _cv("MS:1000619", value="200", unitAccession="UO:0000018")
        + _cv("MS:1000618", value="460", unitAccession="UO:0000018")
        + '<binaryDataArrayList><binaryDataArray>'
        + _cv("MS:1000617", unitAccession="UO:0000018")
        + '<binary>not base64!</binary></binaryDataArray></binaryDataArrayList></spectrum>'
    )


def _group(accession="MS:1000804"):
    return (
        '<referenceableParamGroupList><referenceableParamGroup id="optical">'
        + _cv(accession)
        + '</referenceableParamGroup></referenceableParamGroupList>'
    )


@pytest.mark.parametrize("accession", sorted(NON_MS_SPECTRUM_TERMS))
@pytest.mark.parametrize("referenced", [False, True])
def test_mixed_ms_and_optical_stream_preserves_ms_arrays_and_source_indices(
    tmp_path, accession, referenced
):
    optical = _optical(accession=accession, referenced=referenced)
    path = write_mzml(
        tmp_path / "mixed.mzML",
        [
            mzml_spectrum(1, 0, [100], [1000]),
            optical,
            mzml_spectrum(2, 0.5, [50], [500], level=2),
            optical.replace('optical=1', 'optical=2'),
            mzml_spectrum(3, 1, [100], [2000]),
        ],
        groups=_group(accession),
    )
    scans = list(iter_scans(path, decode_ms2=True))
    assert [scan.native_id for scan in scans] == ["scan=1", "scan=2", "scan=3"]
    assert [scan.spectrum_index for scan in scans] == [0, 2, 4]
    assert [scan.mz.tolist() for scan in scans] == [[100], [50], [100]]
    info = inspect_acquisition(path)
    assert info.spectra_inspected == 5
    assert info.non_ms_spectra_inspected == 2
    assert info.ms_levels == {"1": 2, "2": 1}
    assert info.inspection_complete and not info.is_dia


def test_optical_entries_do_not_change_ms1_only_classification(tmp_path):
    path = write_mzml(
        tmp_path / "ms1_uv.mzML", [_optical(), mzml_spectrum(1, 0, [100], [1000])]
    )
    info = inspect_acquisition(path)
    assert info.mode == "ms1_only"
    assert info.ms_levels == {"1": 1}
    assert info.spectra_inspected == 2 and info.non_ms_spectra_inspected == 1


@pytest.mark.parametrize("type_first", [False, True])
def test_explicit_optical_type_precedes_mass_metadata_and_array_handling(tmp_path, type_first):
    spectrum = mzml_spectrum(1, 0, [200, 460], [1000, 2000]).replace(
        'accession="MS:1000511" value="1"', 'accession="MS:1000511" value="not-an-ms-level"'
    )
    designation = _cv("MS:1000804")
    if type_first:
        spectrum = spectrum.replace('defaultArrayLength="2">', 'defaultArrayLength="2">' + designation)
    else:
        spectrum = spectrum.replace("</spectrum>", designation + "</spectrum>")
    path = write_mzml(tmp_path / "optical_mass_fields.mzML", [spectrum])
    assert list(iter_scans(path, decode_ms2=True)) == []
    info = inspect_acquisition(path)
    assert info.mode == "non_ms_only" and info.ms_levels == {}
    assert info.non_ms_spectra_inspected == 1


def test_optical_only_file_reports_no_mass_spectra_before_calibration(tmp_path, monkeypatch):
    path = write_mzml(tmp_path / "optical.mzML", [_optical()])
    info = inspect_acquisition(path)
    assert info.mode == "non_ms_only" and info.ms_levels == {}
    monkeypatch.setattr(
        "masscubex.pipeline.calibrate_file", lambda *a, **k: pytest.fail("optical calibrated")
    )
    output = tmp_path / "features.parquet"
    with pytest.raises(ValueError, match="Non-MS-only export contains no mass spectra"):
        process_file(path, output)
    assert not output.exists()


@pytest.mark.parametrize(
    "extra", ["none", "header", "unreferenced_group", "nested", "wavelength", "name"]
)
def test_missing_ms_level_is_still_an_error_without_spectrum_type(tmp_path, extra):
    spectrum = mzml_spectrum(1, 0, [100], [1000]).replace(
        '<cvParam accession="MS:1000511" value="1" />', ""
    )
    assert 'accession="MS:1000511"' not in spectrum
    groups = ""
    if extra == "header":
        groups = '<fileDescription><fileContent>' + _cv("MS:1000804")
        groups += '</fileContent></fileDescription>'
    elif extra == "unreferenced_group":
        groups = _group()
    elif extra == "nested":
        spectrum = spectrum.replace("</scan>", _cv("MS:1000804") + "</scan>")
    elif extra == "wavelength":
        spectrum = spectrum.replace(
            "</spectrum>",
            _cv("MS:1000619", value="200", unitAccession="UO:0000018") + "</spectrum>",
        )
    elif extra == "name":
        spectrum = spectrum.replace(
            "</spectrum>",
            _cv("MS:9999999", name="electromagnetic radiation spectrum") + "</spectrum>",
        )
    path = write_mzml(tmp_path / "missing_level.mzML", [spectrum], groups=groups)
    info = inspect_acquisition(path)
    assert info.non_ms_spectra_inspected == 0
    assert info.ms_levels == {"0": 1}
    with pytest.raises(ValueError, match="mzML spectrum has no MS level"):
        list(iter_scans(path, decode_ms2=True))


def test_optical_content_does_not_change_mass_feature_table(sample_file, tmp_path):
    baseline = process_file(sample_file, tmp_path / "baseline.parquet")
    sample_file.write_text(
        sample_file.read_text()
        .replace("<spectrumList>", "<spectrumList>" + _optical())
        .replace("</spectrumList>", _optical(2, referenced=True) + "</spectrumList>")
        .replace("<run>", _group() + "<run>")
    )
    result = process_file(sample_file, tmp_path / "mixed.parquet")
    assert result.ms1_scans == baseline.ms1_scans
    assert list(iter_feature_rows(result.output)) == list(iter_feature_rows(baseline.output))
    metadata = json.loads(Path(result.output + ".json").read_text())
    assert metadata["acquisition"]["non_ms_spectra_inspected"] == 2
    assert metadata["acquisition"]["ms_levels"] == {"1": 9}
