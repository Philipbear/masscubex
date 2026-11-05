from dataclasses import replace
import csv
import json

import numpy as np
import pytest

from masscubex import Config, process_file
from masscubex.features import Feature
from masscubex.ms2 import MS2Collector, associate_ms2
from masscubex.readers import MS2Scan, Scan, iter_scans
from conftest import mzml_spectrum, mzxml_scan, write_mzml, write_mzxml
from test_summary_reuse import reference_segments


def precursor_mzml(mass=100, intensity=200, charge=1, width=2, ref="scan=10"):
    ion = f'<cvParam accession="MS:1000744" value="{mass}"/>' if mass else ""
    if intensity is not None:
        ion += f'<cvParam accession="MS:1000042" value="{intensity}"/>'
    if charge is not None:
        ion += f'<cvParam accession="MS:1000041" value="{charge}"/>'
    window = '<cvParam accession="MS:1000827" value="100"/>'
    if width is not None:
        window += (
            f'<cvParam accession="MS:1000828" value="{width / 2}"/>'
            f'<cvParam accession="MS:1000829" value="{width / 2}"/>'
        )
    return (
        f'<precursorList><precursor spectrumRef="{ref}">'
        f"<isolationWindow>{window}</isolationWindow>"
        f"<selectedIonList><selectedIon>{ion}</selectedIon></selectedIonList>"
        "</precursor></precursorList>"
    )


def feature(mz=100, polarity="positive", start=0, end=6, apex=3):
    return Feature(mz, apex, start, end, 100, 330, 7, 0, None, polarity)


def ms2(index, rt=3, mz=100, intensity=100, polarity="positive", **kwargs):
    return MS2Scan(
        rt,
        f"scan={index}",
        index,
        index,
        polarity,
        precursor_mz=mz,
        precursor_intensity=intensity,
        precursor_count=1,
        intensity_source="metadata" if intensity else "unavailable",
        **kwargs,
    )


def test_mzml_precursors_native_ids_indices_and_lazy_fragment_decoding(tmp_path):
    first = mzml_spectrum(10, 2, [100], [500])
    second = mzml_spectrum(
        20, 2.1, [50, 80], [20, 70], level=2, precursor=precursor_mzml()
    ).replace('id="scan=20"', 'id="controllerType=0 controllerNumber=1 scan=20"')
    third = mzml_spectrum(
        30, 2.2, [40], [1], level=2, precursor=precursor_mzml(intensity=None)
    ).replace('id="scan=30"', 'id="sample=1 period=1 cycle=8 experiment=2"')
    path = write_mzml(tmp_path / "ids.mzML.gz", [first, second, third], indexed=True)
    scans = list(iter_scans(path))
    assert [s.spectrum_index for s in scans] == [0, 1, 2]
    a, b = scans[1:]
    assert a.scan_number == 20
    assert a.precursor_mz == 100
    assert a.precursor_intensity == 200
    assert a.precursor_charge == 1
    assert a.precursor_ref == "scan=10"
    assert a.isolation_width == 2
    assert a.mz is None
    assert b.scan_number is None  # Do not invent scan numbers from array indices.
    assert b.precursor_intensity is None
    decoded = list(iter_scans(path, decode_ms2=True))
    np.testing.assert_array_equal(decoded[1].mz, [50, 80])
    np.testing.assert_array_equal(decoded[1].intensity, [20, 70])


@pytest.mark.parametrize("mode", [None, "profile"])
def test_noncentroid_ms2_metadata_associable_but_not_fragment_evidence(tmp_path, mode):
    spectrum = mzml_spectrum(5, 1, [50], [10], level=2, mode=mode, precursor=precursor_mzml())
    path = write_mzml(tmp_path / "ms2.mzML", [spectrum])
    (record,) = iter_scans(path, decode_ms2=True)
    assert record.mz is None
    assert record.precursor_mz == 100


def test_isolation_target_only_requires_known_width_and_multiprecursors_stay_ambiguous(tmp_path):
    one = mzml_spectrum(1, 1, [], [], level=2, precursor=precursor_mzml(mass=None))
    two = mzml_spectrum(2, 2, [], [], level=2, precursor=precursor_mzml(mass=None, width=None))
    multiple = precursor_mzml().replace(
        "</selectedIonList>",
        '<selectedIon><cvParam accession="MS:1000744" value="150"/></selectedIon></selectedIonList>',
    )
    three = mzml_spectrum(3, 3, [], [], level=2, precursor=multiple)
    path = write_mzml(tmp_path / "precursors.mzML", [one, two, three])
    a, b, c = iter_scans(path)
    assert (a.precursor_mz, a.precursor_mz_source) == (100, "isolation_target")
    assert b.precursor_mz is None
    assert c.precursor_count == 2 and c.precursor_mz is None


def test_nested_mzxml_keeps_parent_before_ms2_and_preserves_precursor_metadata(tmp_path):
    precursor = '<precursorMz precursorScanNum="10" precursorIntensity="900" precursorCharge="-1" windowWideness="3">100</precursorMz>'
    child = mzxml_scan(15, 2.1, [50], [80], level=2, polarity="-", precursor=precursor)
    parent = mzxml_scan(10, 2, [100], [500], polarity="-", nested=child)
    path = write_mzxml(tmp_path / "nested.mzXML", [parent])
    a, b = iter_scans(path, decode_ms2=True)
    assert isinstance(a, Scan) and isinstance(b, MS2Scan)
    assert (a.spectrum_index, b.spectrum_index) == (0, 1)
    assert (b.scan_number, b.native_id, b.precursor_ref) == (15, "15", "10")
    assert (b.precursor_intensity, b.precursor_charge, b.isolation_width) == (900, -1, 3)
    assert b.polarity == "negative"
    np.testing.assert_array_equal(b.mz, [50])


def test_best_ms2_is_precursor_intensity_or_configurable_apex_and_all_ids_survive():
    f = feature()
    records = [
        ms2(19, rt=1, intensity=900),
        ms2(42, rt=3, intensity=100),
        ms2(85, rt=4, intensity=None),
    ]
    records[0].intensity = np.array([1.0])
    records[1].intensity = np.array([1e9])  # Fragment TIC must not affect the best scan.
    stats = {}
    associate_ms2([f], records, Config(best_ms2="precursor_intensity"), stats)
    assert f.ms2_scan_numbers == (19, 42, 85)
    assert f.ms2_native_ids == ("scan=19", "scan=42", "scan=85")
    assert f.ms2_spectrum_indices == (19, 42, 85)
    assert f.best_ms2_scan_number == 19
    assert f.best_ms2_precursor_intensity == 900
    assert f.best_ms2_precursor_intensity_source == "metadata"
    associate_ms2([f], records, Config(best_ms2="apex"), {})
    assert f.best_ms2_scan_number == 42
    assert stats["ms2_assigned"] == 3


def test_missing_fragment_intensities_keep_first_acquisition():
    f = feature()
    records = [
        ms2(30, rt=4, intensity=None),
        ms2(20, rt=2, intensity=None),
        ms2(10, rt=0, intensity=None),
    ]
    associate_ms2([f], records, Config(), {})
    assert f.ms2_scan_numbers == (10, 20, 30)
    assert f.best_ms2_scan_number == 10
    assert f.best_ms2_precursor_intensity is None
    assert f.best_ms2_precursor_intensity_source == "unavailable"
    assert f.best_ms2_selection_method == "apex_fallback"
    associate_ms2([f], [], Config(), {})
    assert f.ms2_count == 0 and f.best_ms2_scan_number is None
    assert f.best_ms2_score is None


@pytest.mark.parametrize("number", [None, 19, 42])
def test_best_scan_identifier_fallback_only_for_missing_or_duplicate_numbers(number):
    f = feature()
    records = [ms2(19, rt=1), ms2(42, rt=3)]
    records[1].scan_number = number
    records[1].native_id = "cycle=2 experiment=1"
    stats = {}
    associate_ms2([f], records, Config(best_ms2="apex"), stats)
    assert f.best_ms2_spectrum_index == 42
    assert stats["best_ms2_identifier_fallbacks"] == (
        [{"feature_id": 1, "native_id": "cycle=2 experiment=1", "spectrum_index": 42}]
        if number in (None, 19)
        else []
    )


@pytest.mark.parametrize("format_", ["mzml", "mzxml"])
def test_default_best_ms2_uses_actual_fragment_sum_before_pruning(tmp_path, format_):
    if format_ == "mzml":
        scans = [
            mzml_spectrum(10, 2, [30], [100], level=2, precursor=precursor_mzml(intensity=900)),
            mzml_spectrum(
                20, 3, [30, 40, 50], [80, 70, 60], level=2, precursor=precursor_mzml(intensity=100)
            ),
        ]
        # A vendor TIC may refer to pre-export samples; compute from the exported arrays.
        scans[0] = scans[0].replace(
            "<scanList>", '<cvParam accession="MS:1000285" value="999999"/><scanList>'
        )
        path = write_mzml(tmp_path / "tic.mzML", scans)
    else:
        scans = [
            mzxml_scan(
                10,
                2,
                [30],
                [100],
                level=2,
                precursor='<precursorMz precursorIntensity="900">100</precursorMz>',
            ),
            mzxml_scan(
                20,
                3,
                [30, 40, 50],
                [80, 70, 60],
                level=2,
                precursor='<precursorMz precursorIntensity="100">100</precursorMz>',
            ),
        ]
        scans[0] = scans[0].replace('msLevel="2"', 'msLevel="2" totIonCurrent="999999"')
        path = write_mzxml(tmp_path / "tic.mzXML", scans)
    records = list(iter_scans(path, sum_ms2=True))
    assert [r.fragment_intensity for r in records] == [100, 210]
    assert all(r.mz is None for r in records)  # Scalar-only requests discard arrays.
    # Toy intensities are below the instrument default; isolate selection and
    # evidence pruning with an explicit zero MS2 floor.
    config = Config(isf_max_peaks=1, ms2_min_intensity=0)
    collector = MS2Collector(config)
    list(collector.ms1_stream(iter_scans(path, decode_ms2=True, sum_ms2=True)))
    assert all(len(r.mz) == 1 for r in collector.records)
    f = feature()
    associate_ms2([f], collector.records, config, {})
    assert f.best_ms2_scan_number == 20
    assert f.best_ms2_score == f.best_ms2_fragment_intensity == 210
    assert f.best_ms2_selection_method == "fragment_intensity"
    assert f.best_ms2_precursor_intensity == 100


def test_profile_ms2_sum_request_preserves_decoded_samples_for_preprocessing(tmp_path):
    path = write_mzml(
        tmp_path / "profile.mzML",
        [
            mzml_spectrum(
                1,
                1,
                [49.99, 50, 50.01],
                [-1, 10, 20],
                level=2,
                mode="profile",
                precursor=precursor_mzml(),
            )
        ],
    )
    (record,) = iter_scans(path, decode_ms2=True, sum_ms2=True)
    assert record.fragment_intensity == 30
    np.testing.assert_array_equal(record.mz, [49.99, 50, 50.01])
    np.testing.assert_array_equal(record.intensity, [0, 10, 20])


def test_preceding_ms1_fallback_checks_parent_polarity_time_and_mass():
    parent = Scan(2, np.array([99.999, 100.006]), np.array([500.0, 1000.0]), "positive", "scan=10")
    records = [
        ms2(11, rt=2.1, intensity=None, precursor_ref="scan=10"),
        ms2(12, rt=2.2, intensity=2000),
        ms2(13, rt=2.3, intensity=0, precursor_ref="scan=9"),
        ms2(14, rt=2.4, intensity=None, polarity="negative"),
        ms2(15, rt=50, intensity=None),
        ms2(16, rt=2.5, mz=200, intensity=None),
    ]
    collector = MS2Collector(Config())
    assert list(collector.ms1_stream([parent, *records])) == [parent]
    assert records[0].precursor_intensity == 500  # Closest centroid, not highest nearby ion.
    assert records[0].intensity_source == "preceding_ms1"
    assert records[1].precursor_intensity == 2000
    assert records[1].intensity_source == "metadata"
    assert all(r.precursor_intensity is None for r in records[2:])


def test_profile_precursor_fallback_uses_apex_not_sum():
    parent = Scan(
        2, np.array([99.999, 100.0, 100.001]), np.array([50.0, 100.0, 60.0]), "positive", "p", False
    )
    record = ms2(2, rt=2.1, intensity=None)
    list(MS2Collector(Config()).ms1_stream([parent, record]))
    assert record.precursor_intensity == 100


def test_assignment_handles_split_peaks_explicit_width_guard_and_usable_precursors():
    features = [feature(start=0, end=4, apex=2), feature(start=4, end=9, apex=7)]
    records = [
        ms2(1, rt=1),
        ms2(2, rt=8),
        ms2(3, rt=4),
        ms2(4, rt=3, isolation_width=25),
        replace(ms2(5), precursor_count=2),
        ms2(6, rt=None),
        ms2(7, mz=None),
        ms2(8, rt=30),
        ms2(9, polarity="negative"),
    ]
    stats = {}
    associate_ms2(features, records, Config(ms2_max_isolation_width=4), stats)
    assert features[0].ms2_scan_numbers == (1, 3, 5)
    assert features[1].ms2_scan_numbers == (2,)
    assert stats["ms2_assigned"] == 4 and stats["ms2_unassigned"] == 5
    assert stats["ms2_skipped_wide_isolation"] == 1
    assert stats["ms2_skipped_multiprecursor"] == 0
    assert stats["ms2_missing_metadata"] == 2


@pytest.mark.parametrize("width", [None, 2, 4, 4.01, 25, 400])
def test_default_has_no_isolation_width_limit(width):
    f = feature()
    stats = {}
    associate_ms2([f], [ms2(9, isolation_width=width)], Config(), stats)
    assert f.ms2_scan_numbers == (9,)
    assert stats["ms2_skipped_wide_isolation"] == 0


def test_every_unassigned_record_has_a_reason_and_original_identity():
    records = [
        ms2(10, mz=None),
        ms2(11, rt=None),
        replace(ms2(12), precursor_count=2, precursor_mz=None),
        ms2(13, isolation_width=25),
        ms2(14, mz=200),
        ms2(15, rt=100),
        ms2(16, polarity="negative"),
        ms2(17),
    ]
    stats = {}
    assignments = associate_ms2([feature()], records, Config(ms2_max_isolation_width=4), stats)
    expected = [
        "missing_precursor_mz",
        "missing_retention_time",
        "missing_precursor_mz",
        "isolation_width_limit",
        "no_feature_mz",
        "no_feature_rt",
        "no_feature_polarity",
        "matched",
    ]
    assert [a.reason for a in assignments] == expected
    assert [r["spectrum_index"] for r in stats["unassigned_ms2"]] == list(range(10, 17))
    assert stats["ms2_assigned"] + stats["ms2_unassigned"] == len(records)
    assert sum(stats["ms2_unassigned_reasons"].values()) == 7


def test_unknown_ms2_polarity_does_not_cross_assign_switching_runs():
    features = [feature(), feature(polarity="negative")]
    stats = {}
    associate_ms2(features, [ms2(1, polarity="unknown")], Config(), stats)
    assert stats["ms2_ambiguous_polarity"] == 1
    assert not any(f.ms2_count for f in features)


def test_end_to_end_ms2_and_quantification_tsv_resume(tmp_path):
    spectra = []
    for i, height in enumerate([10, 40, 70, 100, 70, 40, 10]):
        spectra.append(mzml_spectrum(i * 10, i, [100], [height]))
        if i in (1, 2, 3):
            precursor = precursor_mzml(intensity={1: 900, 2: None, 3: 100}[i], ref=f"scan={i * 10}")
            spectrum = mzml_spectrum(
                100 + i, i + 0.1, [50], [1e5 * i], level=2, precursor=precursor
            )
            if i == 2:
                spectrum = spectrum.replace('id="scan=102"', 'id="cycle=2 experiment=1"')
            spectra.append(spectrum)
    path = write_mzml(tmp_path / "run.mzML", spectra)
    output = tmp_path / "run.tsv"
    config = Config(
        quant_method="peak_area", best_ms2="precursor_intensity", include_peak_area=True
    )
    process_file(path, output, config=config)
    with output.open() as stream:
        (row,) = csv.DictReader(stream, delimiter="\t")
    assert "peak_height" not in row
    trace = np.column_stack(
        (
            np.full(7, 100),
            [10, 40, 70, 100, 70, 40, 10],
            np.arange(7),
            np.zeros(7),
        )
    )
    _, (reference,) = reference_segments(trace, Config(segment=False))
    assert float(row["peak_area"]) == reference.peak_area
    assert "intensity" not in row and "quant_method" not in row
    assert json.loads(row["ms2_scan_numbers"]) == [101, None, 103]
    assert "ms2_spectrum_indices" not in row and "ms2_native_ids" not in row
    metadata = json.loads((output.parent / (output.name + ".json")).read_text())
    assert metadata["ms2_identifier_fallbacks"] == [
        {
            "feature_id": 1,
            "ms2_list_index": 1,
            "scan_number": None,
            "native_id": "cycle=2 experiment=1",
            "spectrum_index": 4,
        }
    ]
    assert row["best_ms2_scan_number"] == "101"
    assert "best_ms2_precursor_intensity" not in row
    assert "best_ms2_precursor_intensity_source" not in row
    assert "trace" not in row and "_ms2_records" not in row
    assert process_file(path, output, config=config, resume=True).status == "skipped"
    process_file(
        path,
        output,
        config=replace(config, quant_method="peak_height", associate_ms2=False),
        resume=True,
    )
    with output.open() as stream:
        (row,) = csv.DictReader(stream, delimiter="\t")
    assert row["peak_height"] == "100.0" and row["peak_area"] == f"{reference.peak_area:.1f}"
    assert json.loads(row["ms2_scan_numbers"]) == [] and row["best_ms2_scan_number"] == ""
