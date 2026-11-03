from dataclasses import replace
import csv
import json

import numpy as np
import pytest

from masscubex import Config, process_file
from masscubex.adducts import Adduct, DEFAULT_ADDUCT_IDS, catalog, parse_selection
from masscubex.features import Feature
from masscubex.grouping import group_ions, trace_correlation
from masscubex.ms2 import associate_ms2
from masscubex.readers import MS2Scan
from conftest import mzml_spectrum, write_mzml


def feature(mz, *, polarity="positive", y=None, times=None):
    y = np.array([1, 3, 8, 12, 8, 3, 1] if y is None else y, dtype=float)
    times = np.arange(len(y), dtype=float) if times is None else np.array(times, dtype=float)
    return Feature(
        mz,
        float(times[np.argmax(y)]),
        float(times[0]),
        float(times[-1]),
        float(y.max()),
        100,
        int(np.count_nonzero(y)),
        0,
        None,
        polarity,
        trace=np.column_stack((y, times)),
    )


def test_user_selection_and_ionic_mass_formula_for_both_polarities_and_charges():
    assert DEFAULT_ADDUCT_IDS == tuple(range(1, 16)) + (22, 23)
    ions = catalog()
    assert ions[1].mass_shift == pytest.approx(1.0072764523, abs=1e-10)
    assert ions[5].mass_shift == pytest.approx(-1.0072764523, abs=1e-10)
    assert ions[6].mass_shift == pytest.approx(34.9694012619, abs=1e-10)
    assert ions[8].mass_shift == pytest.approx(44.9982028513, abs=1e-10)
    for ion in ions.values():
        assert ion.mz(300) > 0
        assert ion.neutral_mass(ion.mz(300)) == pytest.approx(300, abs=1e-10)
    assert ions[18].mz(300) == pytest.approx(148.9927235477, abs=1e-10)
    assert ions[1].mz(300) - ions[22].mz(300) == pytest.approx(18.010564684, abs=1e-9)
    assert ions[1].mz(300) - ions[23].mz(300) == pytest.approx(36.021129368, abs=1e-9)
    assert ions[22].kind == ions[23].kind == "neutral_loss"


def test_selection_validation_and_custom_ions_roundtrip():
    custom = Adduct(100, "[M+Li]+", 1, 1, 7.015454857)
    assert parse_selection("1-3,5,[M+Li]+", (custom,)) == (1, 2, 3, 5, 100)
    config = Config(adduct_ids=(1, 100), custom_adducts=(custom,))
    assert Config.from_dict(json.loads(json.dumps(config.to_dict()))) == config
    for value in ("22", "999", "5-1", "", "1,unknown"):
        with pytest.raises(ValueError):
            parse_selection(value)
    with pytest.raises(ValueError, match="Duplicate"):
        Config(custom_adducts=(replace(custom, id=1),))
    with pytest.raises(ValueError, match="compatible parent"):
        Config(adduct_ids=(5, 22))
    with pytest.raises(ValueError, match="Unknown configuration"):
        Config.from_dict({"typo": 1})


@pytest.mark.parametrize(
    "ids,polarity", [((1, 2), "positive"), ((5, 8), "negative"), ((1, 12), "positive")]
)
def test_adduct_hypotheses_require_coeluting_correlated_traces(ids, polarity):
    ions = catalog()
    a, b = [ions[i] for i in ids]
    features = [feature(a.mz(300), polarity=polarity), feature(b.mz(300), polarity=polarity)]
    stats = {}
    group_ions(features, [], Config(group_ions=True, adduct_ids=ids), stats)
    assert stats["adduct_relations"] == 1
    left, right = [f.ion_relations[0] for f in features]
    assert left["ion"] == a.name and right["ion"] == b.name
    assert (left["partner_feature_id"], right["partner_feature_id"]) == (2, 1)
    assert left["neutral_mass"] == pytest.approx(300)
    assert abs(left["mass_error_ppm"]) < 1e-6
    assert left["shared_scans"] == 7 and left["correlation"] == pytest.approx(1)
    assert left["status"] == "candidate"
    assert all(f.trace is None for f in features)


@pytest.mark.parametrize(
    "case",
    [
        "constant",
        "short",
        "disjoint",
        "opposite",
        "anticorrelated",
        "shifted_mass",
        "low_overlap",
        "apex_far",
    ],
)
def test_insufficient_grouping_evidence_never_assumed_perfect(case):
    a, b = catalog()[1], catalog()[2]
    first = feature(a.mz(300))
    second = feature(b.mz(300))
    if case == "constant":
        second = feature(b.mz(300), y=[5] * 7)
    elif case == "short":
        first = feature(a.mz(300), y=[1, 5, 1])
        second = feature(b.mz(300), y=[2, 10, 2])
    elif case == "disjoint":
        second = feature(b.mz(300), times=np.arange(7) + 0.5)
    elif case == "opposite":
        second.polarity = "negative"
    elif case == "anticorrelated":
        second = feature(b.mz(300), y=[12, 10, 5, 1, 5, 10, 12])
    elif case == "shifted_mass":
        second.mz += 0.02
    elif case == "low_overlap":
        second = feature(b.mz(300), y=[1, 3, 8, 12, 8, 3, 1] + [1] * 13)
    elif case == "apex_far":
        second.rt_apex_seconds = 10
    stats = {}
    group_ions([first, second], [], Config(group_ions=True, adduct_ids=(1, 2)), stats)
    assert stats["adduct_relations"] == 0
    assert not first.ion_relations and not second.ion_relations


def test_correlation_pairs_acquisition_times_and_excludes_missing_scans():
    a = np.array([[1.0, 0], [0, 1], [3, 2], [5, 4]])
    b = np.array([[2.0, 0], [500, 1], [6, 2], [9, 3], [10, 4]])
    count, correlation = trace_correlation(a, b)
    assert count == 3 and correlation == pytest.approx(1)
    assert trace_correlation(a[:2], b)[1] == 0


def test_water_losses_distinguished_from_adducts_and_ms2_supported_isf():
    ions = catalog()
    parent = ions[1].mz(300)
    water = ions[22].mz(300)
    two_water = ions[23].mz(300)
    features = [feature(parent), feature(water), feature(two_water), feature(210)]
    records = [
        MS2Scan(
            2,
            "scan=20",
            20,
            20,
            "positive",
            parent,
            5000,
            1,
            precursor_count=1,
            centroided=True,
            mz=np.array([30.0]),
            intensity=np.array([100.0]),
        ),
        MS2Scan(
            3,
            "scan=30",
            30,
            30,
            "positive",
            parent,
            100,
            1,
            precursor_count=1,
            centroided=True,
            mz=np.array([210.0, water]),
            intensity=np.array([100.0, 50.0]),
        ),
    ]
    config = Config(group_ions=True, adduct_ids=(1, 22, 23), best_ms2="precursor_intensity")
    stats = {}
    associate_ms2(features, records, config, stats)
    assert features[0].best_ms2_scan_number == 20
    group_ions(features, records, config, stats)
    assert stats["adduct_relations"] == 0
    assert stats["neutral_loss_relations"] == 3  # Also the successive 18 Da loss.
    assert stats["isf_relations"] == 1
    links = {link["partner_feature_id"]: link for link in features[0].ion_relations}
    assert links[2]["relation"] == links[3]["relation"] == "neutral_loss"
    assert links[2]["ms2_supported"] is True
    assert links[2]["evidence_ms2_scan_number"] == 30  # Any MS2, not only best MS2.
    assert links[3]["ms2_supported"] is False
    assert links[4]["relation"] == "in_source_fragment"
    assert links[4]["ms2_supported"] is True
    assert links[2]["fragment_ion"] == "[M+H-H2O]+"
    assert links[3]["fragment_ion"] == "[M+H-2H2O]+"
    assert all(link["status"] == "candidate" for link in links.values())


def test_precursor_charge_restricts_adduct_and_fragment_interpretations():
    ions = catalog()
    parent = ions[1].mz(300)
    features = [feature(parent), feature(ions[2].mz(300)), feature(210)]
    record = MS2Scan(
        3,
        "scan=1",
        1,
        1,
        "positive",
        parent,
        100,
        2,
        precursor_count=1,
        centroided=True,
        mz=np.array([210.0]),
        intensity=np.array([100.0]),
    )
    config = Config(group_ions=True, adduct_ids=(1, 2))
    associate_ms2(features, [record], config, {})
    stats = {}
    group_ions(features, [record], config, stats)
    assert not any(f.ion_relations for f in features)
    assert stats["adduct_relations"] == stats["isf_relations"] == 0


def test_optional_grouping_end_to_end_and_disabled_isf(tmp_path):
    ions = catalog()
    masses = sorted([ions[1].mz(300), ions[2].mz(300), ions[22].mz(300)])
    spectra = [
        mzml_spectrum(i, i, masses, [y, 2 * y, 3 * y]) for i, y in enumerate([1, 3, 8, 12, 8, 3, 1])
    ]
    path = write_mzml(tmp_path / "ions.mzML", spectra)
    config = Config(group_ions=True, adduct_ids=(1, 2, 22))
    output = tmp_path / "ions.tsv"
    process_file(path, output, config=config)
    with output.open() as stream:
        rows = list(csv.DictReader(stream, delimiter="\t"))
    links = [link for row in rows for link in json.loads(row["ion_relations"])]
    assert {link["relation"] for link in links} == {"adduct", "neutral_loss"}
    assert len(links) == 4  # Each relationship appears on both features.
    assert process_file(path, output, config=config, resume=True).status == "skipped"
    config = replace(config, group_isf=False, associate_ms2=False)
    process_file(path, output, config=config, resume=True)
    with output.open() as stream:
        rows = list(csv.DictReader(stream, delimiter="\t"))
    assert {link["relation"] for row in rows for link in json.loads(row["ion_relations"])} == {
        "adduct"
    }
