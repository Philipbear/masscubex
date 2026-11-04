import numpy as np
import pytest

from masscubex import Config
from masscubex.gapfill import _Measurement
from masscubex.gapfill_refinement import merge_refined, refine_measurement
from masscubex.ms2 import MS2Assignment
from test_ms2 import feature


def measure(seed=30.0, second=0.0, close_competitor=False):
    f = feature(mz=100.0001, start=10, end=60)
    f.feature_origin = "ms2_seeded"
    f.seed_precursor_mz, f.seed_rt_seconds = f.mz, seed
    f.peak_height = f.peak_area = f.intensity = None
    f.scan_count = 0
    f.ms1_quant_status = "no_ms1_scans"
    f.trace = np.empty((0, 2))
    m = _Measurement(0, "positive", refine=True)
    for t in np.arange(10.0, 60.25, 0.25):
        y = 4 * np.exp(-0.5 * ((t - 30) / 2) ** 2) + second * np.exp(-0.5 * ((t - 42) / 2) ** 2)
        mass = 100.008 if close_competitor else 100.03
        # A much stronger adjacent ion must not replace the selected precursor.
        m.samples.append(
            (t, np.array([100.0, mass]), np.array([y, y * (1.2 if close_competitor else 100)]))
        )
    m.scans = len(m.samples)
    return f, m


def test_raw_weak_measurement_uses_narrow_mass_and_local_boundaries():
    f, m = measure(second=8)
    config = Config(min_intensity=100, smoothing_sigma=1.2, peak_distance=8)
    refine_measurement(f, m, config)
    assert f.ms1_quant_status == "gap_filled" and f.gapfill_status == "peak_bounds"
    assert f.mz == f.seed_precursor_mz == 100.0001 and f.peak_height == pytest.approx(4, rel=1e-6)
    assert f.rt_start_seconds > 10 and f.rt_end_seconds < 36.1
    assert f.peak_area == pytest.approx(4 * 2 * np.sqrt(2 * np.pi), rel=0.03)
    assert f.gapfill_mz_radius_da == 0.002
    # Independently sum the actual original-time trapezoids inside reported bounds.
    selected = [(t, y[0]) for t, _, y in m.samples if f.rt_start_seconds <= t <= f.rt_end_seconds]
    t, y = np.asarray(selected).T
    assert f.peak_area == pytest.approx(np.sum(np.diff(t) * (y[:-1] + y[1:]) / 2))


def test_ion_outside_ms1_tolerance_cannot_compete_for_gapfill():
    f, m = measure(close_competitor=True)
    refine_measurement(f, m, Config())
    assert f.ms1_quant_status == "gap_filled"
    assert f.peak_height == pytest.approx(4)
    assert f.mz == f.seed_precursor_mz


def test_repeated_queries_merge_by_recovered_apex_not_broad_window():
    features = []
    for rt in (29.5, 30.5, 42):
        f, m = measure(seed=rt, second=8)
        refine_measurement(f, m, Config(smoothing_sigma=1.2, peak_distance=8))
        features.append(f)
    assignments = merge_refined(
        features, [MS2Assignment(i, "seeded") for i in range(3)], 0, Config(), {}
    )
    assert len(features) == 2
    assert [a.feature_index for a in assignments] == [0, 0, 1]
    assert [f.rt_apex_seconds for f in features] == [30, 42]


def test_limits_and_no_signal_are_explicit():
    f, m = measure()
    m.sampling_limit = True
    refine_measurement(f, m, Config())
    assert f.ms1_quant_status == "sampling_limit" and f.peak_height is None
    f, m = measure()
    m.samples = [(t, np.array([]), np.array([])) for t, _, _ in m.samples]
    refine_measurement(f, m, Config())
    assert f.ms1_quant_status == "no_ms1_signal" and f.peak_height is None


def test_upstream_scan_based_segmentation_recovers_resolved_peaks_across_scan_rates():
    from masscubex.features import summarize_trace

    apices = []
    for dt in (0.2, 0.4, 0.8):
        t = np.arange(0, 30.001, dt)
        y = 100 * np.exp(-0.5 * ((t - 8) / 1.5) ** 2) + 150 * np.exp(-0.5 * ((t - 20) / 1.5) ** 2)
        trace = np.column_stack((np.full(len(t), 100.0), y, t, np.zeros(len(t))))
        found = list(summarize_trace(trace, Config(smoothing_sigma=0.6, peak_distance=3)))
        assert len(found) == 2
        apices.append([f.rt_apex_seconds for f in found])
    np.testing.assert_allclose(apices, [[8, 20]] * 3, atol=0.4)


def test_broad_search_cannot_borrow_a_different_ms1_ion(tmp_path, monkeypatch):
    from masscubex.tuning import TuningResult
    from masscubex import pipeline
    from conftest import mzml_spectrum, write_mzml
    from test_ms2 import precursor_mzml
    from test_ms2_seeding import outputs

    config = Config(min_intensity=100, smoothing_sigma=1.2)
    # Isolate the target-query schedule from the independently tested calibration pass.
    monkeypatch.setattr(
        pipeline,
        "calibrate_file",
        lambda *args, **kwargs: TuningResult(
            {"positive": config}, {"polarities": {"positive": {"status": "tuned"}}}
        ),
    )
    y = [1, 2, 3, 4, 3, 2, 1]
    spectra = [
        mzml_spectrum(i, i, [100.0001 if i % 2 == 0 else 99.9995], [v]) for i, v in enumerate(y)
    ]
    spectra.append(
        mzml_spectrum(100, 3.1, [50], [1], level=2, precursor=precursor_mzml(100.0498, ref=""))
    )
    source = write_mzml(tmp_path / "edge.mzML", spectra)
    (row,), _ = outputs(source, tmp_path / "out.tsv", config)
    assert "gapfill_status" not in row
    assert row["ms1_quant_status"] == "imputed_noise"
    assert float(row["peak_height"]) == 10  # 0.1 * explicit noise floor, not the neighboring ion.
    assert row["peak_area"] == ""
    assert row["scan_count"] == "0"
    assert row["ms2_count"] == "1"
