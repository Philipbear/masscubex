import numpy as np
import pytest
import base64
import zlib
from lxml import etree

from masscubex.readers import _duration, iter_ms1
from conftest import mzml_spectrum, mzxml_scan, write_mzml, write_mzxml


@pytest.mark.parametrize("precision", [32, 64])
@pytest.mark.parametrize("compress", [False, True])
@pytest.mark.parametrize("gz", [False, True])
def test_mzml_encodings_and_units(tmp_path, precision, compress, gz):
    path = tmp_path / ("data.mzML.gz" if gz else "data.mzML")
    write_mzml(
        path,
        [
            mzml_spectrum(
                1,
                1.5,
                [100.125, 200.5],
                [123, 456],
                unit="UO:0000031",
                precision=precision,
                compress=compress,
            )
        ],
        indexed=True,
    )
    (scan,) = iter_ms1(path)
    assert scan.rt == 90
    assert scan.mz.dtype == np.float64
    np.testing.assert_array_equal(scan.mz, [100.125, 200.5])
    np.testing.assert_array_equal(scan.intensity, [123, 456])


@pytest.mark.parametrize("precision", [32, 64])
@pytest.mark.parametrize("endian", ["<", ">"])
@pytest.mark.parametrize("compress", [False, True])
def test_mzxml_binary_and_nested_ms2(tmp_path, precision, endian, compress):
    nested = '<scan num="2" msLevel="2" retentionTime="PT1.1S" peaksCount="1"><peaks>invalid base64!</peaks></scan>'
    path = write_mzxml(
        tmp_path / "data.mzXML.gz",
        [
            mzxml_scan(
                1,
                1,
                [100.125],
                [123],
                precision=precision,
                endian=endian,
                compress=compress,
                nested=nested,
            ),
            mzxml_scan(3, 2, [100.125], [321]),
        ],
    )
    scans = list(iter_ms1(path))
    assert [s.rt for s in scans] == [1, 2]
    np.testing.assert_array_equal(scans[0].mz, [100.125])


def test_mzml_skips_ms2_without_decoding(tmp_path):
    bad_ms2 = '<spectrum id="2"><cvParam accession="MS:1000511" value="2"/><binaryDataArrayList><binaryDataArray><binary>invalid</binary></binaryDataArray></binaryDataArrayList></spectrum>'
    path = write_mzml(tmp_path / "test.mzML", [bad_ms2, mzml_spectrum(3, 2, [100], [10])])
    assert len(list(iter_ms1(path))) == 1


def test_referenceable_params_and_sanitizing(tmp_path):
    group = (
        '<referenceableParamGroupList><referenceableParamGroup id="spectrum_params">'
        '<cvParam accession="MS:1000511" value="1"/><cvParam accession="MS:1000127"/>'
        "</referenceableParamGroup></referenceableParamGroupList>"
    )
    path = write_mzml(
        tmp_path / "data.mzML",
        [
            mzml_spectrum(1, 0, [300, 100, 200, 100], [0, 5, 8, 7], refs=True),
        ],
        groups=group,
    )
    (scan,) = iter_ms1(path)
    # Raw reading preserves duplicate centroids; upstream preprocessing merges them.
    np.testing.assert_array_equal(scan.mz, [100, 100, 200])
    np.testing.assert_array_equal(scan.intensity, [5, 7, 8])
    from masscubex import Config
    from masscubex.upstream_adapter import preprocess_ms1

    masses, intensities, _ = preprocess_ms1(scan, Config(noise_mode="none"))
    np.testing.assert_array_equal(masses, [100, 200])
    np.testing.assert_array_equal(intensities, [12, 8])


def test_global_centroid_flag(tmp_path):
    path = write_mzxml(
        tmp_path / "test.mzXML",
        [mzxml_scan(1, 0, [100], [20], centroid=None)],
        processing='<dataProcessing centroided="1"/>',
    )
    assert len(list(iter_ms1(path))) == 1


@pytest.mark.parametrize("kind", ["mzML", "mzXML"])
def test_profile_metadata_preserved_even_when_assuming_centroid(tmp_path, kind):
    path = tmp_path / f"profile.{kind}"
    if kind == "mzML":
        write_mzml(path, [mzml_spectrum(1, 0, [100], [20], mode="profile")])
    else:
        write_mzxml(path, [mzxml_scan(1, 0, [100], [20], centroid="0")])
    (scan,) = iter_ms1(path, assume_centroid=True)
    assert scan.centroided is False


def test_missing_mode_needs_explicit_assumption(tmp_path):
    path = write_mzml(tmp_path / "test.mzML", [mzml_spectrum(1, 0, [100], [1], mode=None)])
    with pytest.raises(ValueError, match="Unknown centroid"):
        list(iter_ms1(path))
    assert len(list(iter_ms1(path, assume_centroid=True))) == 1


@pytest.mark.parametrize(
    "change,match",
    [
        (lambda s: s.replace('defaultArrayLength="1"', 'defaultArrayLength="2"'), "length"),
        (lambda s: s.replace("UO:0000010", "invalid"), "time unit"),
        (lambda s: s.replace("MS:1000574", "MS:9999999"), "compression"),
        (
            lambda s: s.replace('value="0" unitAccession', 'value="nan" unitAccession'),
            "retention time",
        ),
    ],
)
def test_corrupt_metadata_fails(tmp_path, change, match):
    path = write_mzml(tmp_path / "bad.mzML", [change(mzml_spectrum(1, 0, [100], [20]))])
    with pytest.raises(ValueError, match=match):
        list(iter_ms1(path))


@pytest.mark.parametrize("mz,intensity", [([np.nan], [10]), ([100], [-1]), ([100], [np.inf])])
def test_invalid_numeric_arrays_fail(tmp_path, mz, intensity):
    path = write_mzml(tmp_path / "bad.mzML", [mzml_spectrum(1, 0, mz, intensity)])
    with pytest.raises(ValueError):
        list(iter_ms1(path))


def test_empty_scan_preserved(tmp_path):
    path = write_mzml(tmp_path / "empty.mzML", [mzml_spectrum(1, 0, [], [])])
    (scan,) = iter_ms1(path)
    assert len(scan.mz) == 0


def test_mz_precision_preserved(tmp_path):
    value = 1000.123456789123
    path = write_mzml(tmp_path / "accurate.mzML", [mzml_spectrum(1, 0, [value], [10])])
    assert next(iter_ms1(path)).mz[0] == value


@pytest.mark.parametrize("encoding", ["linear", "pic", "slof"])
@pytest.mark.parametrize("compressed", [False, True])
def test_numpress_encodings(tmp_path, encoding, compressed):
    pynumpress = pytest.importorskip("pynumpress")
    values = np.array([100.0, 200.0, 300.0])
    encode = getattr(pynumpress, "encode_" + encoding)
    if encoding == "pic":
        raw = encode(values).tobytes()
    else:
        fixed = getattr(pynumpress, "optimal_" + encoding + "_fixed_point")(values)
        raw = encode(values, fixed).tobytes()
    if compressed:
        raw = zlib.compress(raw)
    term = {
        "linear": "MS:1002746" if compressed else "MS:1002312",
        "pic": "MS:1002747" if compressed else "MS:1002313",
        "slof": "MS:1002748" if compressed else "MS:1002314",
    }[encoding]
    element = etree.fromstring(mzml_spectrum(1, 0, values, values).encode())
    for array in element.findall("binaryDataArrayList/binaryDataArray"):
        for cv in array.findall("cvParam"):
            if cv.get("accession") == "MS:1000574":
                cv.set("accession", term)
        array.find("binary").text = base64.b64encode(raw).decode()
    path = write_mzml(tmp_path / "numpress.mzML", [etree.tostring(element).decode()])
    (scan,) = iter_ms1(path)
    np.testing.assert_allclose(scan.mz, values, rtol=1e-4)
    np.testing.assert_allclose(scan.intensity, values, rtol=1e-4)


def test_profile_zeros_and_negative_baseline_preserved_as_boundaries(tmp_path):
    path = write_mzml(
        tmp_path / "profile.mzML",
        [mzml_spectrum(1, 0, [99, 100, 101], [-5, 10, 0], mode="profile")],
    )
    (scan,) = iter_ms1(path)
    np.testing.assert_array_equal(scan.mz, [99, 100, 101])
    np.testing.assert_array_equal(scan.intensity, [0, 10, 0])


def test_empty_mzxml_peaks_and_duration_units(tmp_path):
    path = write_mzxml(
        tmp_path / "empty.mzXML",
        [
            '<scan num="1" msLevel="1" peaksCount="0" centroided="1" retentionTime="PT1H2M3.5S"><peaks/></scan>',
        ],
    )
    (scan,) = iter_ms1(path)
    assert scan.rt == 3723.5
    assert scan.mz.size == 0
    assert _duration("P1D") == 86400
    for invalid in ("PT", "P", "P1DT", "PT-1S"):
        with pytest.raises(ValueError):
            _duration(invalid)


@pytest.mark.parametrize("kind", ["mzML", "mzXML"])
def test_streaming_reader_ignores_stale_offset_indexes(tmp_path, kind):
    path = tmp_path / ("indexed." + kind)
    offsets = '<index name="scan">' + '<offset id="1">invalid-offset</offset>' * 1000 + "</index>"
    if kind == "mzML":
        write_mzml(path, [mzml_spectrum(1, 0, [100], [50])], indexed=True)
        path.write_text(
            path.read_text().replace(
                "</indexedmzML>", "<indexList>" + offsets + "</indexList></indexedmzML>"
            )
        )
    else:
        write_mzxml(path, [mzxml_scan(1, 0, [100], [50])])
        path.write_text(path.read_text().replace("</mzXML>", offsets + "</mzXML>"))
    (scan,) = iter_ms1(path)
    np.testing.assert_array_equal(scan.mz, [100])
