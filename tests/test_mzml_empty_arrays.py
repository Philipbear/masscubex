"""Empty-export compatibility must not hide damaged nonempty binary arrays."""

import base64
import zlib

from lxml import etree
import numpy as np
import pytest

from masscubex import iter_feature_rows, process_file
from masscubex.readers import iter_scans
from conftest import mzml_spectrum, write_mzml


def _empty_payload_spectrum(*, level=1, precision=64, compress=True):
    element = etree.fromstring(
        mzml_spectrum(23, 1, [], [], level=level, precision=precision, compress=compress).encode()
    )
    for array in element.findall("binaryDataArrayList/binaryDataArray"):
        array.set("encodedLength", "0")
        array.find("binary").text = "\n  \t"
    return element


@pytest.mark.parametrize("level", [1, 2])
@pytest.mark.parametrize("precision", [32, 64])
@pytest.mark.parametrize("compress", [False, True])
def test_declared_empty_arrays_can_omit_binary_payload(tmp_path, level, precision, compress):
    element = _empty_payload_spectrum(level=level, precision=precision, compress=compress)
    path = write_mzml(tmp_path / "empty.mzML", [etree.tostring(element).decode()])
    (scan,) = iter_scans(path, decode_ms2=True)
    assert scan.native_id == "scan=23"
    assert scan.mz.size == scan.intensity.size == 0
    if level == 2:
        assert scan.fragment_intensity == 0


def test_empty_ms2_payload_can_be_summed_without_retaining_arrays(tmp_path):
    element = _empty_payload_spectrum(level=2)
    path = write_mzml(tmp_path / "empty.mzML", [etree.tostring(element).decode()])
    (scan,) = iter_scans(path, sum_ms2=True)
    assert scan.fragment_intensity == 0
    assert scan.mz is scan.intensity is None


@pytest.mark.parametrize("expected", [0, 1])
def test_nonempty_truncated_zlib_payload_is_never_treated_as_empty(tmp_path, expected):
    element = _empty_payload_spectrum()
    element.set("defaultArrayLength", str(expected))
    raw = zlib.compress(np.array([100.0], dtype="<f8").tobytes())[:-1]
    for array in element.findall("binaryDataArrayList/binaryDataArray"):
        encoded = base64.b64encode(raw).decode()
        array.set("encodedLength", str(len(encoded)))
        array.find("binary").text = encoded
    path = write_mzml(tmp_path / "truncated.mzML", [etree.tostring(element).decode()])
    with pytest.raises(ValueError, match="incomplete or truncated stream"):
        list(iter_scans(path))


def test_empty_payload_cannot_replace_a_declared_nonempty_array(tmp_path):
    element = _empty_payload_spectrum()
    element.set("defaultArrayLength", "1")
    path = write_mzml(tmp_path / "missing.mzML", [etree.tostring(element).decode()])
    with pytest.raises(ValueError, match="incomplete or truncated stream"):
        list(iter_scans(path))


def test_nonempty_array_cannot_hide_behind_zero_default_length(tmp_path):
    spectrum = mzml_spectrum(1, 0, [100], [20]).replace(
        'defaultArrayLength="1"', 'defaultArrayLength="0"'
    )
    path = write_mzml(tmp_path / "wrong_length.mzML", [spectrum])
    with pytest.raises(ValueError, match="length 1 != defaultArrayLength 0"):
        list(iter_scans(path))


@pytest.mark.parametrize("damage", ["encoded_length", "precision", "compression", "binary"])
def test_empty_payload_still_requires_consistent_encoding_metadata(tmp_path, damage):
    element = _empty_payload_spectrum()
    array = element.find("binaryDataArrayList/binaryDataArray")
    if damage == "encoded_length":
        array.set("encodedLength", "12")
    elif damage == "binary":
        array.remove(array.find("binary"))
    else:
        accession = "MS:1000523" if damage == "precision" else "MS:1000574"
        parameter = next(cv for cv in array.findall("cvParam") if cv.get("accession") == accession)
        array.remove(parameter)
    path = write_mzml(tmp_path / "bad_metadata.mzML", [etree.tostring(element).decode()])
    with pytest.raises(ValueError):
        list(iter_scans(path))


def test_empty_ms2_export_does_not_interrupt_default_feature_extraction(sample_file, tmp_path):
    baseline = process_file(sample_file, tmp_path / "baseline.parquet")
    empty_ms2 = etree.tostring(_empty_payload_spectrum(level=2)).decode()
    sample_file.write_text(
        sample_file.read_text().replace("</spectrumList>", empty_ms2 + "</spectrumList>")
    )
    result = process_file(sample_file, tmp_path / "with_empty_ms2.parquet")
    assert result.status == "ok" and result.ms1_scans == baseline.ms1_scans
    assert list(iter_feature_rows(result.output)) == list(iter_feature_rows(baseline.output))
