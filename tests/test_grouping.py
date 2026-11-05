from dataclasses import replace
import json

import numpy as np
import pytest

from masscubex import Config, process_file
from masscubex.adducts import Adduct, DEFAULT_ADDUCT_IDS, catalog, parse_selection
from masscubex.features import Feature
from masscubex.grouping import group_ions, trace_correlation
from masscubex.grouping_spectra import GroupingSpectra
from masscubex import iter_feature_rows
from types import SimpleNamespace
import upstream_masscube as upstream
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


def raw_scans(features):
    """Independent synthetic MS1 source, including signals outside partner bounds."""
    scans = {}
    for f in features:
        for intensity, time in f.trace:
            masses = scans.setdefault((f.polarity, float(time)), {})
            masses[f.mz] = max(masses.get(f.mz, 0), intensity)
    return [
        (polarity, time, np.array(sorted(masses.items()), dtype=float))
        for (polarity, time), masses in sorted(scans.items())
    ]


def run_grouping(features, records=(), config=None, scans=None):
    config = config or Config()
    stats = {}
    with GroupingSpectra() as spectra:
        for polarity, time, points in scans or raw_scans(features):
            spectra.append(time, polarity, points[:, 0], points[:, 1])
        group_ions(features, records, config, stats, spectra=spectra)
    return stats


def reference_grouping(features, records, config, monkeypatch):
    """Run the original function with independent NumPy EIC queries."""
    ions = catalog()
    for polarity, attribute in [("positive", "POS_ADDUCTS"), ("negative", "NEG_ADDUCTS")]:
        dictionary = {
            rule.name: SimpleNamespace(
                mass_shift=rule.mass_shift,
                charge=abs(rule.charge),
                mol_multiplier=rule.molecules,
                considered=rule.id in config.adduct_ids,
            )
            for rule in ions.values()
            if rule.polarity == polarity
        }
        monkeypatch.setattr(upstream, attribute, dictionary)
    scans = raw_scans(features)
    times = np.array([s[1] / 60 for s in scans])
    original = []
    for index, f in enumerate(features):
        g = upstream.Feature()
        g.id, g.mz, g.rt, g.peak_height = index, f.mz, f.rt_apex_seconds / 60, f.peak_height
        g.signals = np.column_stack((np.full(len(f.trace), f.mz), f.trace[:, 0]))
        g.rt_seq = f.trace[:, 1] / 60
        g.scan_idx = int(np.argmin(abs(times - g.rt)))
        best = next((r for r in records if r.spectrum_index == f.best_ms2_spectrum_index), None)
        g.ms2 = (
            SimpleNamespace(signals=np.column_stack((best.mz, best.intensity)))
            if best is not None and best.mz is not None
            else None
        )
        original.append(g)

    def get_eic_data(mass, rt, rt_range, mz_tol=0.005):
        selected = np.flatnonzero((times >= rt_range[0]) & (times <= rt_range[1]))
        y = []
        for scan in selected:
            points = scans[scan][2]
            keep = (points[:, 0] >= mass - mz_tol) & (points[:, 0] <= mass + mz_tol)
            y.append(np.max(points[keep, 1], initial=0.0))
        return times[selected], np.column_stack((np.full(len(y), mass), y)), selected

    data = SimpleNamespace(
        features=original,
        metadata=SimpleNamespace(ion_mode=features[0].polarity),
        scans=[SimpleNamespace(signals=s[2]) for s in scans],
        get_eic_data=get_eic_data,
        params=SimpleNamespace(
            rt_tol_feature_grouping=config.group_rt_tolerance / 60,
            mz_tol_feature_grouping=config.group_mz_tolerance_da,
            scan_scan_cor_tol=config.group_min_correlation,
            isotope_rel_int_limit=config.isotope_rel_int_limit,
        ),
    )
    upstream.group_features_single_file(data)
    return sorted(data.features, key=lambda f: f.id)


@pytest.mark.parametrize("polarity,base_id,partner_id", [("positive", 1, 2), ("negative", 5, 8)])
@pytest.mark.parametrize(
    "case", ["pair", "strong_partner", "isotope", "child_isotope", "fragment", "bad_shape", "short"]
)
def test_grouping_assignments_match_upstream(polarity, base_id, partner_id, case, monkeypatch):
    ions = catalog()
    root = ions[base_id].mz(300)
    second = ions[partner_id].mz(300)
    y = np.array([1, 3, 8, 12, 8, 3, 1], dtype=float)
    if case == "short":
        y = np.array([2.0, 1.0])
    features = [
        feature(root, polarity=polarity, y=y * 10),
        feature(second, polarity=polarity, y=y * 3),
    ]
    if case == "bad_shape":
        features[1] = feature(second, polarity=polarity, y=np.full(len(y), 5.0))
    if case == "strong_partner":
        features[1] = feature(second, polarity=polarity, y=y * 20)
    if case in {"isotope", "child_isotope"}:
        center = root if case == "isotope" else second
        features.extend(
            [
                feature(center + offset, polarity=polarity, y=y)
                for offset in [1.003355, 2.00671, 3.010065]
            ]
        )
    records = []
    config = Config(adduct_ids=(base_id, partner_id))
    if case == "fragment":
        features.append(feature(210, polarity=polarity, y=y * 2))
        records = [
            MS2Scan(
                3.0,
                "scan=1",
                1,
                1,
                polarity,
                root,
                1000,
                1,
                precursor_count=1,
                centroided=True,
                mz=np.array([210.0]),
                intensity=np.array([100.0]),
            )
        ]
        associate_ms2(features, records, config, {})
    expected = reference_grouping(features, records, config, monkeypatch)
    run_grouping(features, records, config)
    for ours, old in zip(features, expected):
        assert (
            ours.group_id,
            ours.adduct,
            ours.is_isotope,
            ours.isotope_state,
            ours.is_in_source_fragment,
        ) == (
            old.feature_group_id,
            old.adduct_type,
            old.is_isotope,
            old.isotope_state,
            old.is_in_source_fragment,
        )


def test_monomer_pair_is_not_duplicated_as_dimer_pair():
    ions = catalog()
    y = np.array([1, 3, 8, 12, 8, 3, 1])
    features = [feature(ions[1].mz(300), y=y * 10), feature(ions[2].mz(300), y=y)]
    stats = run_grouping(features, config=Config(adduct_ids=(1, 2, 12, 13)))
    assert [f.adduct for f in features] == ["[M+H]+", "[M+Na]+"]
    assert stats["adduct_relations"] == 1 and features[0].group_id == features[1].group_id
    assert all("2M" not in str(f.ion_relations) for f in features)


def test_real_dimer_at_correct_mass_is_allowed_but_anchor_remains_monomer():
    ions = catalog()
    features = [
        feature(ions[1].mz(300), y=np.arange(1.0, 8.0) * 10),
        feature(ions[12].mz(300), y=np.arange(1.0, 8.0)),
    ]
    run_grouping(features, config=Config(adduct_ids=(1, 12)))
    assert [f.adduct for f in features] == ["[M+H]+", "[2M+H]+"]
    assert features[0].group_id == features[1].group_id


def test_only_best_ms2_supplies_fragments_and_water_losses_are_ion_forms():
    ions = catalog()
    parent = ions[1].mz(300)
    y = np.array([1, 3, 8, 12, 8, 3, 1])
    features = [
        feature(parent, y=y * 10),
        feature(ions[22].mz(300), y=y * 2),
        feature(210, y=y),
        feature(220, y=y),
    ]
    records = [
        MS2Scan(
            3.0,
            f"scan={i}",
            i,
            i,
            "positive",
            parent,
            100,
            1,
            precursor_count=1,
            centroided=True,
            mz=np.array([mass]),
            intensity=np.array([height]),
        )
        for i, mass, height in [(1, 210.0, 100.0), (2, 220.0, 50.0)]
    ]
    config = Config(adduct_ids=(1, 22, 23))
    associate_ms2(features, records, config, {})
    assert features[0].best_ms2_scan_number == 1
    run_grouping(features, records, config)
    assert features[1].adduct == "[M+H-H2O]+" and not features[1].is_in_source_fragment
    assert features[2].is_in_source_fragment and not features[3].is_in_source_fragment
    assert features[2].adduct == "[M+H]+"
    assert features[3].group_id != features[0].group_id


def test_upstream_later_isf_assignment_replaces_prior_link(monkeypatch):
    ions = catalog()
    y = np.array([1, 3, 8, 12, 8, 3, 1])
    parent, sodium = ions[1].mz(300), ions[2].mz(300)
    features = [feature(parent, y=y * 10), feature(sodium, y=y)]
    record = MS2Scan(
        3.0,
        "scan=1",
        1,
        1,
        "positive",
        parent,
        100,
        1,
        precursor_count=1,
        mz=np.array([sodium]),
        intensity=np.array([100.0]),
    )
    associate_ms2(features, [record], Config(), {})
    config = Config(adduct_ids=(1, 2))
    expected = reference_grouping(features, [record], config, monkeypatch)
    stats = run_grouping(features, [record], config)
    assert features[1].adduct == expected[1].adduct_type == "[M+H]+"
    assert features[1].is_in_source_fragment == expected[1].is_in_source_fragment is True
    assert stats["isf_relations"] == 1 and stats["adduct_relations"] == 0
    assert len(features[0].ion_relations) == len(features[1].ion_relations) == 1


def test_correlation_retains_upstream_sparse_and_constant_rules():
    a = np.column_stack(([1.0, 2.0], [0.0, 1.0]))
    assert trace_correlation(a, a)[1] == 1
    constant = np.column_stack((np.ones(5), np.arange(5.0)))
    assert np.isnan(trace_correlation(constant, constant)[1])


def test_grouping_output_and_disabled_flags(tmp_path):
    ions = catalog()
    masses = [ions[1].mz(300), ions[1].mz(300) + 1.003355, ions[2].mz(300)]
    source = write_mzml(
        tmp_path / "ions.mzML",
        [
            mzml_spectrum(i, i, masses, [10 * y, y, 3 * y])
            for i, y in enumerate([1, 3, 8, 12, 8, 3, 1])
        ],
    )
    config = Config(auto_tune=False, adduct_ids=(1, 2), noise_mode="none")
    tsv, parquet = tmp_path / "ions.tsv", tmp_path / "ions.parquet"
    process_file(source, tsv, config=config)
    process_file(source, parquet, config=config)
    rows = list(iter_feature_rows(tsv))
    assert rows == list(iter_feature_rows(parquet))
    assert len(rows) == 3 and len({r["group_id"] for r in rows}) == 1
    assert sum(r["isotope_state"] != "M+0" for r in rows) == 1
    assert {r["adduct"] for r in rows} == {"[M+H]+", "[M+Na]+"}
    assert all("seed_precursor_mz" not in r and "seed_rt_seconds" not in r for r in rows)
    process_file(source, tsv, config=replace(config, group_ions=False), overwrite=True)
    assert all(
        r["adduct"] is None and r["group_id"] is None and r["isotope_state"] is None
        for r in iter_feature_rows(tsv)
    )


def test_file_backed_eic_uses_source_outside_partner_trace_and_cleans_up(tmp_path):
    with GroupingSpectra(tmp_path) as spectra:
        for t in range(7):
            spectra.append(
                float(t),
                "positive",
                np.array([100.0, 100.003, 200.0]),
                np.array([1.0, 3.0 + t, 50.0]),
            )
            spectra.append(float(t) + 0.1, "negative", np.array([100.0]), np.array([1000.0]))
        np.testing.assert_array_equal(
            spectra.eic(100, np.arange(7.0), "positive", 0.005), np.arange(3.0, 10.0)
        )
        assert spectra.eic(100, np.arange(7.0), "negative", 0.005).sum() == 0
        handle = spectra.file
    assert handle.closed and not list(tmp_path.iterdir())
