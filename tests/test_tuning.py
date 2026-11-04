from dataclasses import replace

import numpy as np
import pytest

from masscubex import Config
from masscubex.cli import main
from masscubex.tuning import Representative, derive_settings, representative
from test_ms2 import ms2


def peaks(count=24):
    return [
        Representative(
            100 + 20 * i,
            30 + 4 * i,
            27 + 4 * i,
            33 + 4 * i,
            1e5 + i,
            6,
            12,
            (100 + 20 * i) * 3e-6,
            3,
            13,
        )
        for i in range(count)
    ]


def test_independent_minimum_evidence_and_outliers():
    base = Config()
    anchors = peaks()
    # One extreme high-intensity and mass-error outlier must not broaden every EIC.
    anchors += [
        replace(
            anchors[0], mz=800, rt=140, height=1e12, fwhm=80, mass_error_da=0.04, mass_error_ppm=50
        )
    ]
    spectra = [
        ms2(i, rt=p.rt, mz=p.mz * (1 + 14e-6), precursor_mz_source="selected_ion")
        for i, p in enumerate(anchors[:24])
    ]
    tuned, info = derive_settings(anchors, [0.5] * 100 + [60], spectra, base, "positive")
    assert info["representative_count"] == 24
    assert info["status"] == info["ms2_status"] == "tuned"
    assert tuned.mz_tolerance_ppm == pytest.approx(5)
    assert tuned.ms2_mz_tolerance_ppm == pytest.approx(14)
    assert tuned.precursor_max_age_seconds == 1
    assert tuned.noise_sigma == 3 and tuned.noise_fraction == 0.1
    _, few = derive_settings(peaks(19), [0.5] * 100, spectra, base, "positive")
    assert few["status"] == "insufficient_representatives"
    partial, info = derive_settings(peaks(), [0.5] * 100, spectra[:19], base, "positive")
    assert info["status"] == "tuned" and info["ms2_status"] != "tuned"
    assert partial.ms2_mz_tolerance_da == base.ms2_mz_tolerance_da
    assert partial.ms2_mz_tolerance_ppm == base.ms2_mz_tolerance_ppm


def test_quantized_masses_do_not_reject_every_nonzero_error():
    anchors = [
        replace(p, mass_error_da=0, mass_error_ppm=0)
        if i < 20
        else replace(p, mass_error_da=p.mz * 1e-6, mass_error_ppm=1)
        for i, p in enumerate(peaks())
    ]
    _, info = derive_settings(anchors, [1] * 20, [], Config(), "positive")
    assert info["representative_count"] == 24


def test_overrides_and_unknown_polarity_are_protected():
    base = Config(mz_tolerance_ppm=50, tuning_locked=("ms2_rt_tolerance",))
    tuned, info = derive_settings(peaks(), [0.5] * 20, [], base, "positive")
    assert tuned.mz_tolerance_ppm == 50 and tuned.ms2_rt_tolerance == 2
    assert "mz_tolerance_ppm" in info["protected_parameters"]
    unknown, info = derive_settings(peaks(), [0.5] * 20, [], Config(), "unknown")
    assert info["status"] == "unknown_polarity"
    assert unknown.mz_tolerance_ppm == 10


@pytest.mark.parametrize("defect", ["flat", "edge", "competitor", "undersampled", "plateau"])
def test_strong_bad_traces_are_not_representatives(defect):
    t = np.linspace(0, 20, 81)
    y = 1e6 * np.exp(-0.5 * ((t - 10) / 2) ** 2)
    competing = np.zeros(len(t))
    if defect == "flat":
        y[:] = 1e9
    elif defect == "edge":
        y = 1e6 * np.exp(-0.5 * (t / 2) ** 2)
    elif defect == "competitor":
        competing[40] = 1
    elif defect == "undersampled":
        y = 1e6 * np.exp(-0.5 * ((t - 10) / 0.1) ** 2)
    else:
        y = np.minimum(y, 3e5)
    peak, reason = representative(
        np.column_stack((np.full(len(t), 100), y, t, competing)), Config()
    )
    assert peak is None and reason


def test_valid_representative_and_deterministic_learning():
    t = np.linspace(0, 20, 81)
    y = 1e6 * np.exp(-0.5 * ((t - 10) / 2) ** 2)
    peak, reason = representative(
        np.column_stack((100 + 0.0002 * np.sin(t), y, t, np.zeros(len(t)))), Config()
    )
    assert reason is None and peak.fwhm == pytest.approx(4.71, abs=0.1)
    a, ai = derive_settings(peaks(120), [0.5] * 100, [], Config(), "positive")
    b, bi = derive_settings(list(reversed(peaks(120))), [0.5] * 100, [], Config(), "positive")
    assert a == b and ai == bi


def test_five_real_scans_support_calibration_but_missing_scans_do_not():
    t = np.arange(5.0)
    y = 1e6 * np.exp(-0.5 * (t - 2) ** 2)
    trace = np.column_stack((150 + 1e-5 * t, y, t, np.zeros(5)))
    peak, reason = representative(trace, Config())
    assert reason is None and peak.observations == 3
    assert peak.width > peak.fwhm > 0
    trace[0, 1] = 0  # Four measured scans plus a tracker zero are not five scans.
    assert representative(trace, Config()) == (None, "short_trace")


def test_missing_scan_cannot_manufacture_a_peak_flank():
    t = np.arange(7.0)
    y = np.array([1, 10, 40, 70, 90, 100, 0.0]) * 1e4
    trace = np.column_stack((np.full(7, 150.0), y, t, np.zeros(7)))
    assert representative(trace, Config()) == (None, "truncated_peak")


def test_background_and_minor_secondary_peak_do_not_veto_good_main_peak():
    t = np.linspace(0, 40, 161)
    main = 1e6 * np.exp(-0.5 * ((t - 10) / 2) ** 2)
    second = 4e5 * np.exp(-0.5 * ((t - 25) / 1) ** 2)
    y = 1e7 + main + second
    mz = np.where(t < 20, 150.0, 150.005)
    trace = np.column_stack((mz, y, t, np.zeros(len(t))))
    peak, reason = representative(trace, Config())
    assert reason is None and peak.gaussian_similarity >= 0.8
    assert peak.mz == 150 and peak.mass_error_da == 0
    trace[np.argmax(main), 3] = 1
    assert representative(trace, Config()) == (None, "competing_mass")


def test_mass_outliers_are_removed_before_the_representative_quota():
    from masscubex.tuning import _screen_representatives

    good = peaks(20)
    bad = [replace(p, height=1e12, mass_error_da=p.mz * 50e-6, mass_error_ppm=50) for p in peaks(8)]
    selected, before = _screen_representatives(good + bad, Config(tuning_max_peaks=20))
    assert before == 28 and set(selected) == set(good)


def test_cli_defaults_and_explicit_default_value_lock(capsys):
    import json

    main(["extract", "--print-config", "--mz-ppm", "10", "--no-ms2"])
    config = json.loads(capsys.readouterr().out)
    assert config["auto_tune"] and config["group_ions"]
    assert config["tuning_min_peaks"] == 20
    assert config["tuning_locked"] == ["mz_tolerance_ppm"]


@pytest.mark.parametrize("seed", range(30))
def test_compiled_calibration_matches_independent_scipy_measurements(seed):
    from dataclasses import asdict
    from reference_calibration import representative as reference

    rng = np.random.default_rng(seed)
    t = np.cumsum(rng.uniform(0.2, 0.5, 100))
    for noise in (0.0, 0.01, 0.2, 1.0):
        y = 1e5 * (np.exp(-0.5 * ((t - np.median(t)) / 3) ** 2) + noise * rng.random(len(t)))
        trace = np.column_stack((100 + rng.normal(0, 0.0001, len(t)), y, t, np.zeros(len(t))))
        a, a_reason = representative(trace, Config())
        b, b_reason = reference(trace, Config())
        assert a_reason == b_reason
        if b is not None:
            actual, expected = asdict(a), asdict(b)
            assert actual.pop("gaussian_similarity") == pytest.approx(
                expected.pop("gaussian_similarity")
            )
            assert actual == expected
