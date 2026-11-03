"""Known-truth interleaved LCMS fixtures, including unequal acquisition cadence."""

import json

import numpy as np
import pytest

from masscubex import Config, process_file
from masscubex.readers import iter_scans
from masscubex.tuning import calibrate_file
from conftest import mzml_spectrum, mzxml_scan, write_mzml, write_mzxml
from test_ms2 import precursor_mzml
from test_ms2_seeding import outputs


def switching_file(
    tmp_path,
    format_="mzml",
    profile=False,
    weak_seeds=False,
    *,
    precision=64,
    compress=True,
    gz=None,
):
    masses = np.arange(24) * 20 + 100.0
    apices = np.arange(24) * 4 + 30.0
    events = []
    for polarity, dt, offset, sigma, ppm, height in [
        ("positive", 0.5, 0, 3, 3, 1e5),
        ("negative", 1.0, 0.25, 6, 12, 2e5),
    ]:
        for t in np.arange(0, 170, dt) + offset:
            mz = masses * (1 + ppm * 1e-6 * np.sin(t * 1.3))
            y = 1 + height * np.exp(-0.5 * ((t - apices) / sigma) ** 2)
            if weak_seeds and polarity == "positive":
                weak = 4 * np.exp(-0.5 * ((t - 40) / 2) ** 2) + 8 * np.exp(
                    -0.5 * ((t - 52) / 2) ** 2
                )
                mz, y = np.r_[mz, 800.0, 800.04], np.r_[y, weak, 1000 * weak]
            if profile:
                delta = np.arange(-4, 5) * 0.001
                mz = (mz[:, None] + delta).ravel()
                y = (y[:, None] * np.exp(-0.5 * (delta / 0.001) ** 2)).ravel()
            events.append((t, polarity, 1, mz, y, None))
        for i, (mass, t) in enumerate(zip(masses, apices)):
            events.append(
                (
                    t + offset + 0.1,
                    polarity,
                    2,
                    [50],
                    [100 + i],
                    mass * (1 + (7 if polarity == "positive" else 20) * 1e-6),
                )
            )
        if weak_seeds and polarity == "positive":
            for t in (39.6, 40.1, 40.6, 51.6, 52.1, 52.6):
                events.append((t, polarity, 2, [50], [t], 800.0001))
    events.sort(key=lambda e: (e[0], e[2], e[1]))
    spectra = []
    parents = {}
    for i, (t, polarity, level, mz, y, precursor) in enumerate(events, 1):
        if format_ == "mzml":
            spectra.append(
                mzml_spectrum(
                    i,
                    t,
                    mz,
                    y,
                    level=level,
                    precision=precision,
                    compress=compress,
                    polarity=polarity,
                    mode="profile" if profile and level == 1 else "centroid",
                    precursor=precursor_mzml(precursor, ref=f"scan={parents.get(polarity, 0)}")
                    if level == 2
                    else "",
                )
            )
        else:
            spectra.append(
                mzxml_scan(
                    i,
                    t,
                    mz,
                    y,
                    level=level,
                    precision=precision,
                    compress=compress,
                    polarity="+" if polarity == "positive" else "-",
                    centroid="0" if profile and level == 1 else "1",
                    precursor=f'<precursorMz precursorScanNum="{parents.get(polarity, 0)}">{precursor}</precursorMz>'
                    if level == 2
                    else "",
                )
            )
        if level == 1:
            parents[polarity] = i
    writer = write_mzml if format_ == "mzml" else write_mzxml
    gz = format_ == "mzml" if gz is None else gz
    suffix = ("mzML" if format_ == "mzml" else "mzXML") + (".gz" if gz else "")
    return writer(tmp_path / ("switch." + suffix), spectra)


@pytest.mark.parametrize("format_", ["mzml", "mzxml"])
@pytest.mark.parametrize("profile", [False, True])
@pytest.mark.parametrize("precision", [32, 64])
@pytest.mark.parametrize("compress", [False, True])
@pytest.mark.parametrize("gz", [False, True])
def test_separate_tuning_quantification_ms2_and_groups(
    tmp_path, format_, profile, precision, compress, gz
):
    source = switching_file(
        tmp_path, format_, profile, precision=precision, compress=compress, gz=gz
    )
    config = Config(noise_mode="none", tuning_block_scans=2048, tuning_block_seconds=300)
    rows, meta = outputs(source, tmp_path / "out.tsv", config)
    tune = meta["parameter_tuning"]["polarities"]
    pos, neg = tune["positive"], tune["negative"]
    assert pos["status"] == neg["status"] == "tuned"
    assert pos["ms2_status"] == neg["ms2_status"] == "tuned"
    assert pos["representative_count"] >= 20 and neg["representative_count"] >= 20
    assert pos["ms1_interval_seconds"] == 0.5 and neg["ms1_interval_seconds"] == 1
    assert (
        neg["applied_parameters"]["mz_tolerance_ppm"]
        > 2 * pos["applied_parameters"]["mz_tolerance_ppm"]
    )
    assert neg["fwhm_seconds_p10_p50_p90"][1] == pytest.approx(
        2 * pos["fwhm_seconds_p10_p50_p90"][1], rel=0.03
    )
    assert len(rows) == 48 and meta["ms2_assigned"] == 48 and meta["ms2_seeded_features"] == 0
    scans = {s.native_id: s for s in iter_scans(source) if hasattr(s, "precursor_mz")}
    for row in rows:
        polarity = row["polarity"]
        assert row["feature_origin"] == "ms1_detected"
        assert all(
            scans[native].polarity == polarity for native in json.loads(row["ms2_native_ids"])
        )
        assert int(row["ms2_count"]) == 1
        if polarity == "positive":
            assert float(row["peak_height"]) == pytest.approx(100001, rel=1e-6)
            assert float(row["peak_area"]) == pytest.approx(
                100000 * 3 * np.sqrt(2 * np.pi) + 169.5, rel=1e-5
            )
        else:
            assert float(row["peak_height"]) == pytest.approx(
                1 + 200000 * np.exp(-0.5 * (0.25 / 6) ** 2), rel=1e-6
            )
            assert float(row["peak_area"]) == pytest.approx(
                200000 * 6 * np.sqrt(2 * np.pi) + 169, rel=1e-5
            )
        for relation in json.loads(row["ion_relations"]):
            assert rows[relation["partner_feature_id"] - 1]["polarity"] == polarity
    assert (
        process_file(source, tmp_path / "out.tsv", config=config, resume=True).status == "skipped"
    )


def test_unknown_ms2_uses_explicit_parent_only_and_never_a_single_candidate_guess(tmp_path):
    spectra = []
    for i, value in enumerate([1, 3, 8, 12, 8, 3, 1]):
        spectra += [
            mzml_spectrum(10 * i + 1, i, [100], [value]),
            mzml_spectrum(10 * i + 2, i + 0.1, [200], [value * 10], polarity="negative"),
        ]
        if i == 3:
            spectra += [
                mzml_spectrum(
                    90,
                    3.2,
                    [50],
                    [1],
                    level=2,
                    polarity=None,
                    precursor=precursor_mzml(ref="scan=31"),
                ),
                mzml_spectrum(
                    91, 3.3, [50], [2], level=2, polarity=None, precursor=precursor_mzml(ref="")
                ),
            ]
    path = write_mzml(tmp_path / "unknown.mzML", spectra)
    rows, meta = outputs(path, tmp_path / "out.tsv")
    assigned = {native: row for row in rows for native in json.loads(row["ms2_native_ids"])}
    assert assigned["scan=90"]["feature_origin"] == "ms1_detected"
    assert assigned["scan=90"]["polarity"] == "positive"
    assert assigned["scan=91"]["feature_origin"] == "ms2_seeded"
    assert assigned["scan=91"]["ms1_quant_status"] == "ambiguous_polarity"
    assert assigned["scan=91"]["peak_height"] == ""
    assert meta["ms2_polarity_inferences"][0]["native_id"] == "scan=90"


def test_ms1_without_polarity_cannot_be_shared_by_known_opposite_ms2(tmp_path):
    spectra = [mzml_spectrum(i, i, [100], [10], polarity=None) for i in range(7)]
    spectra += [
        mzml_spectrum(
            100 + i, 3.1 + i * 0.1, [50], [1], level=2, polarity=p, precursor=precursor_mzml(ref="")
        )
        for i, p in enumerate(["positive", "negative"])
    ]
    rows, meta = outputs(write_mzml(tmp_path / "missing.mzML", spectra), tmp_path / "out.tsv")
    seeds = [r for r in rows if r["feature_origin"] == "ms2_seeded"]
    assert len(seeds) == 2 and all(r["ms1_quant_status"] == "ambiguous_polarity" for r in seeds)
    assert all(r["peak_height"] == "" for r in seeds)
    assert meta["ms2_assigned"] == 2


def test_conflicting_polarity_metadata_fails_explicitly(tmp_path):
    text = mzml_spectrum(1, 0, [100], [1]).replace(
        "<scanList>", '<cvParam accession="MS:1000129"/><scanList>'
    )
    path = write_mzml(tmp_path / "conflict.mzML", [text])
    with pytest.raises(ValueError, match="polarity"):
        calibrate_file(path, Config())


@pytest.mark.parametrize("profile", [False, True])
def test_tuned_pipeline_preserves_below_noise_seeds_and_distinct_local_peaks(tmp_path, profile):
    source = switching_file(tmp_path, profile=profile, weak_seeds=True)
    rows, meta = outputs(
        source,
        tmp_path / "out.tsv",
        Config(
            min_intensity=100, noise_mode="none", tuning_block_scans=2048, tuning_block_seconds=300
        ),
    )
    seeds = [r for r in rows if r["feature_origin"] == "ms2_seeded"]
    assert meta["parameter_tuning"]["polarities"]["positive"]["status"] == "tuned"
    assert len(seeds) == 2 and meta["ms2_seeded_scans"] == 6
    assert [int(r["ms2_count"]) for r in seeds] == [3, 3]
    assert [float(r["rt_apex_seconds"]) for r in seeds] == [40, 52]
    for row, height in zip(seeds, [4, 8]):
        assert row["gapfill_status"] == "peak_bounds"
        assert float(row["peak_height"]) == pytest.approx(height, rel=1e-5)
        assert float(row["peak_area"]) == pytest.approx(height * 2 * np.sqrt(2 * np.pi), rel=0.05)
        assert float(row["mz"]) == pytest.approx(800, abs=1e-7)
        assert float(row["seed_precursor_mz"]) == 800.0001


def test_grouping_on_by_default_stays_within_polarity(tmp_path):
    from masscubex.adducts import selected_adducts, DEFAULT_ADDUCT_IDS

    ions = {r.id: r for r in selected_adducts(DEFAULT_ADDUCT_IDS, ())}
    spectra = []
    for i, value in enumerate([1, 3, 8, 12, 8, 3, 1]):
        spectra += [
            mzml_spectrum(2 * i, i, [ions[1].mz(300), ions[2].mz(300)], [value, 2 * value]),
            mzml_spectrum(
                2 * i + 1,
                i + 0.1,
                [ions[5].mz(300), ions[6].mz(300)],
                [value, 2 * value],
                polarity="negative",
            ),
        ]
    rows, meta = outputs(write_mzml(tmp_path / "ions.mzML", spectra), tmp_path / "out.tsv")
    assert meta["adduct_relations"] >= 2
    for row in rows:
        links = json.loads(row["ion_relations"])
        assert links
        for link in links:
            assert rows[link["partner_feature_id"] - 1]["polarity"] == row["polarity"]
