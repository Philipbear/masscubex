"""Signal preservation and background-model boundary checks."""

import numpy as np
import pytest

from masscubex.config import Config
from masscubex.readers import MS2Scan, Scan
from masscubex.consistency_noise import (
    Parameters,
    evaluate,
    fit_component,
    intensity_groups,
    mass_recurrence,
    precursor_protection,
)


def frame(mz, intensity, rt=0, polarity="positive", centroided=True):
    order = np.argsort(mz)
    return Scan(
        rt,
        np.array(mz, dtype=np.float32)[order],
        np.array(intensity, dtype=np.float32)[order],
        polarity=polarity,
        centroided=centroided,
    )


def test_distinct_weak_ion_survives_repeated_band():
    mz = np.r_[np.linspace(100, 900, 30), 950]
    intensity = np.r_[np.full(30, 50.0), 30.123]
    scan = frame(mz, intensity)
    decisions, _, _ = evaluate([scan], {"positive": Config()})
    assert decisions[0].keep()[-1]
    assert not decisions[0].keep()[:-1].any()


def test_floor_atom_does_not_collapse_broader_background_spread():
    rng = np.random.default_rng(42)
    mz = rng.uniform(100, 200, 500)
    values = np.r_[np.full(200, 14.0), np.round(rng.normal(65, 12, 300) / 2) * 2]
    atom = values == 14
    component, status = fit_component(mz, values, atom, Parameters())
    assert status == "density_mode"
    assert component.sd > 5
    assert 60 < component.threshold < 200
    assert component.threshold == component.mean + 3 * component.sd


def test_coherence_rescues_quantized_true_ion_without_protecting_random_grass():
    rng = np.random.default_rng(12)
    scans = []
    true_positions = []
    for rt in range(7):
        noise = rng.uniform(100, 900, 100)
        noise = noise[abs(noise - 500) > 0.1]
        scan = frame(np.r_[noise, 500.0], np.full(len(noise) + 1, 50.0), rt=rt)
        true_positions.append(int(np.searchsorted(scan.mz, 500)))
        scans.append(scan)
    original = [(s.mz.copy(), s.intensity.copy()) for s in scans]
    decisions, _, _ = evaluate(scans, {"positive": Config()})
    for i, position in enumerate(true_positions):
        assert not decisions[i].keep(rescue=False)[position]
        assert decisions[i].keep()[position]
        assert np.array_equal(scans[i].mz, original[i][0])
        assert np.array_equal(scans[i].intensity, original[i][1])
    retained = sum(int(d.keep().sum()) for d in decisions)
    assert retained < sum(len(s.mz) for s in scans) * 0.05


def test_polarity_streams_cannot_supply_each_others_recurrence():
    mz = np.linspace(100, 900, 100)
    scans = [
        frame(mz, np.full(100, 50.0), rt=i, polarity="positive" if i == 0 else "negative")
        for i in range(5)
    ]
    decisions, _, _ = evaluate(scans, {"positive": Config(), "negative": Config()})
    assert not decisions[0].coherent.any()
    assert not decisions[0].keep().any()
    assert decisions[2].keep().all()


def test_a_local_envelope_cannot_supply_global_electronic_evidence():
    scan = frame(np.linspace(200, 205, 40), np.full(40, 100.0))
    repeated, _ = intensity_groups(scan.mz, scan.intensity, Parameters())
    assert not repeated.any()
    decisions, _, _ = evaluate([scan], {"positive": Config()})
    assert decisions[0].keep().all()


def test_competing_masses_do_not_count_as_unambiguous_recurrence():
    target = frame([200], [50])
    competitor = frame([199.995, 200.005], [50, 50])
    assert not mass_recurrence(target, [competitor, competitor], Config(), 2)[0]


def test_profile_samples_must_be_mass_detected_first():
    with pytest.raises(ValueError, match="mass-detected"):
        evaluate([frame([100, 101, 102], [10, 20, 10], centroided=False)], {"positive": Config()})


def test_raw_ms2_evidence_protects_an_otherwise_unresolved_short_ion():
    scan = frame(np.r_[np.linspace(100, 900, 100), 950.0], np.full(101, 50.0))
    scan.native_id = "scan=1"
    record = MS2Scan(
        0.1,
        "scan=2",
        1,
        2,
        "positive",
        precursor_mz=950.0,
        precursor_ref="scan=1",
        precursor_count=1,
    )
    evidence = precursor_protection([scan], [record], {"positive": Config()})
    decisions, _, _ = evaluate([scan], {"positive": Config()}, ms2_evidence=evidence)
    assert not decisions[0].coherent[-1]
    assert not decisions[0].keep(rescue=False)[-1]
    assert decisions[0].keep()[-1]
    assert decisions[0].keep().sum() == 1


def test_ms2_protection_obeys_polarity_time_and_ordinary_mass_radius():
    config = Config(mz_tolerance_da=0.001, mz_tolerance_ppm=0)
    scans = [
        frame([500.0, 500.002], [50.0, 50.0], rt=0, polarity="negative"),
        frame([500.0, 500.002], [50.0, 50.0], rt=1),
        frame([500.0, 500.002], [50.0, 50.0], rt=30),
    ]
    record = MS2Scan(1.1, "scan=4", 3, 4, "positive", precursor_mz=500.0, precursor_count=1)
    evidence = precursor_protection(scans, [record], {"positive": config, "negative": config})
    assert not evidence[0].any()  # Opposite polarity cannot contribute.
    assert evidence[1].tolist() == [True, False]  # No enlarged precursor window.
    assert not evidence[2].any()  # A distant survey cannot contribute.


def test_conflicting_ms2_polarity_and_ambiguous_nearest_mass_are_not_protection():
    scan = frame([499.5, 500.5], [50.0, 50.0])
    scan.native_id = "scan=1"
    config = Config(mz_tolerance_da=1.0, mz_tolerance_ppm=0)
    conflicting = MS2Scan(
        0.1,
        "scan=2",
        1,
        2,
        "negative",
        precursor_mz=499.5,
        precursor_ref="scan=1",
        precursor_count=1,
    )
    tied = MS2Scan(
        0.1,
        "scan=3",
        2,
        3,
        "positive",
        precursor_mz=500.0,
        precursor_ref="scan=1",
        precursor_count=1,
    )
    evidence = precursor_protection([scan], [conflicting, tied], {"positive": config})
    assert not evidence[0].any()


def test_repeated_intensity_above_regional_ceiling_is_not_a_veto():
    rng = np.random.default_rng(91)
    scan = frame(
        np.r_[rng.uniform(100, 900, 1500), np.linspace(110, 890, 20)],
        np.r_[rng.choice([14.0, 20.0, 30.0, 40.0, 50.0, 60.0, 70.0], 1500), np.full(20, 250.0)],
    )
    decisions, _, _ = evaluate([scan], {"positive": Config()})
    decision = decisions[0]
    high = scan.intensity == 250.0
    assert intensity_groups(scan.mz, scan.intensity, Parameters())[0][high].all()
    assert (scan.intensity[high] > decision.threshold[high]).all()
    assert decision.keep(rescue=False)[high].all()
    assert not (decision.candidate_noise & (scan.intensity > decision.threshold)).any()


def test_intermittent_physical_ion_is_not_lost_between_two_scan_pools():
    rng = np.random.default_rng(133)
    scans = []
    for rt in range(13):
        mz = rng.uniform(100, 900, 100)
        mz = mz[abs(mz - 950) > 0.1]
        if rt in (0, 4, 8, 12):
            mz = np.r_[mz, 950.0]
        scans.append(frame(mz, np.full(len(mz), 50.0), rt=rt))
    decisions, _, _ = evaluate(scans, {"positive": Config()})
    for index in (0, 4, 8, 12):
        assert decisions[index].coherent[-1]
        assert decisions[index].keep()[-1]
