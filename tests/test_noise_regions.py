"""Independent numerical and end-to-end checks of binned noise and imputation."""

from dataclasses import replace

import numpy as np
import pytest

from masscubex import Config
from masscubex.noise import estimate_noise
from masscubex.noise_regions import regional_thresholds
from masscubex.gapfill_refinement import _radius
from conftest import mzml_spectrum, write_mzml
from test_ms2 import precursor_mzml
from test_ms2_seeding import outputs


def test_one_threshold_per_half_open_bin_including_boundary_and_zero_padding():
    mz = np.array([90, 91, 92, 93, 94, 99.999, 100, 101, 102, 103, 104, 105, 199.999, 200.0])
    intensity = np.array([1, 2, 3, 4, 5, 100, 11, 12, 13, 14, 15, 1000, 0, 1e8], dtype=float)
    mean, sd, threshold = estimate_noise(mz, intensity)
    for a, b, background in ((0, 6, np.arange(1, 6)), (6, 13, np.arange(11, 16))):
        np.testing.assert_allclose(mean[a:b], background.mean())
        np.testing.assert_allclose(sd[a:b], background.std(ddof=0))
        np.testing.assert_allclose(threshold[a:b], background.mean() + 3 * background.std(ddof=0))
    assert threshold[-1] == 0  # One positive sample is insufficient; no neighboring-bin borrowing.
    targets = np.array([250, 199, 99.999, 100, 110, 0], dtype=float)
    queried = regional_thresholds(mz, intensity, targets, Config())
    np.testing.assert_allclose(queried, threshold[[13, 6, 0, 6, 6, 0]])


@pytest.mark.parametrize("estimator", ["low_tail", "clipped"])
def test_query_noise_agrees_with_full_spectrum_bins_and_explicit_floor(estimator):
    rng = np.random.default_rng(109)
    mz = np.sort(rng.uniform(20, 500, 1000))
    intensity = rng.lognormal(4, 3, len(mz))
    config = Config(noise_estimator=estimator, min_intensity=2.0)
    expected = estimate_noise(mz, intensity, config)[2]
    np.testing.assert_array_equal(
        regional_thresholds(mz, intensity, mz[::-1], config), expected[::-1]
    )
    assert regional_thresholds(mz, intensity, np.array([600.0]), config)[0] == 2


@pytest.mark.parametrize(
    "invalid",
    [
        dict(noise_layout="anything"),
        dict(gapfill_noise_fraction=-1),
        dict(gapfill_noise_fraction=1.01),
    ],
)
def test_new_settings_validation(invalid):
    with pytest.raises(ValueError):
        Config(**invalid)


@pytest.mark.parametrize("tuned", [False, True])
@pytest.mark.parametrize("profile", [False, True])
@pytest.mark.parametrize("quant", ["peak_height", "peak_area"])
def test_missing_signal_imputes_height_only_from_compatible_nearest_scans(
    tmp_path, monkeypatch, tuned, profile, quant
):
    from masscubex import pipeline
    from masscubex.tuning import TuningResult

    config = Config(
        auto_tune=False, quant_method=quant, min_scans=20, include_apex_noise_threshold=True
    )
    if tuned:
        monkeypatch.setattr(
            pipeline,
            "calibrate_file",
            lambda source, cfg, **k: TuningResult(
                {"positive": cfg, "negative": cfg},
                {"polarities": {"positive": {"status": "tuned"}, "negative": {"status": "tuned"}}},
            ),
        )
    spectra = []
    for i, level in enumerate([2, 10, 20, 30, 500]):
        for j, polarity in enumerate(("positive", "negative")):
            spectra.append(
                mzml_spectrum(
                    10 * i + j,
                    i + j * 0.01,
                    [110, 120, 130, 140, 150, 160],
                    [level * (1 if j == 0 else 100)] * 6,
                    polarity=polarity,
                    mode="profile" if profile else "centroid",
                )
            )
    spectra.append(
        mzml_spectrum(99, 2.05, [50], [1], level=2, precursor=precursor_mzml(), polarity="positive")
    )
    source = write_mzml(tmp_path / "missing.mzML", spectra)
    (row,), metadata = outputs(source, tmp_path / "output.tsv", config)
    # Positive scan times 1,2,3 are nearest; median thresholds 10,20,30 -> 20.
    assert row["feature_origin"] == "ms2_seeded"
    assert row["ms1_quant_status"] == "imputed_noise"
    if quant == "peak_height":
        assert row["peak_height"] == "2.0"
    else:
        assert "peak_height" not in row
    assert row["peak_area"] == "" and "scan_count" not in row
    assert row["apex_noise_threshold"] == "20.0"
    assert row["ion_relations"] == "[]"
    assert metadata["ms2_assigned"] == 1
    (disabled,), _ = outputs(
        source, tmp_path / "disabled.tsv", replace(config, gapfill_noise_fraction=0)
    )
    assert disabled["ms1_quant_status"] == "no_ms1_signal"
    assert disabled["peak_area"] == ""
    if quant == "peak_height":
        assert disabled["peak_height"] == ""
    else:
        assert "peak_height" not in disabled


def test_tolerance_cap_is_minimum_of_physical_radii_not_separate_components():
    config = Config(
        mz_tolerance_da=0.002,
        mz_tolerance_ppm=1,
        gapfill_mz_tolerance_da=0.001,
        gapfill_mz_tolerance_ppm=10,
    )
    assert _radius(100, config) == 0.001
    assert _radius(500, config) == 0.002


def test_explicit_broad_gapfill_override_is_capped_by_ms1_tolerance(tmp_path):
    spectra = [mzml_spectrum(i, i, [100.008], [1000]) for i in range(3)]
    spectra.append(mzml_spectrum(10, 1, [50], [5], level=2, precursor=precursor_mzml()))
    source = write_mzml(tmp_path / "neighbor.mzML", spectra)
    config = Config(
        auto_tune=False,
        mz_tolerance_da=0.002,
        gapfill_mz_tolerance_da=0.05,
        min_intensity=2000,
    )
    (row,), _ = outputs(source, tmp_path / "out.tsv", config)
    assert "gapfill_mz_radius_da" not in row
    assert row["peak_height"] == "200.0" and row["ms1_quant_status"] == "imputed_noise"
    assert row["peak_area"] == "" and "scan_count" not in row
