"""Regional calibration must survive sparse evidence, boundaries and overrides."""

from dataclasses import replace
import json

import numpy as np
import pytest

from masscubex import Config
from masscubex.features import TraceTracker
from masscubex.gapfill_refinement import _radius
from masscubex.mass_tolerance import MassRegion, mass_radius, ms1_radius, ms1_region_array
from masscubex.readers import iter_ms1
from masscubex.tuning import calibrate_file, derive_settings, representative
from conftest import mzml_spectrum, mzxml_scan, write_mzml, write_mzxml
from test_ms2 import ms2, feature
from test_tuning import peaks


def regional_peaks():
    result = []
    for region in range(5):
        for i, peak in enumerate(peaks()):
            mass = 110 + 200 * region + 7 * i
            error_ppm = (region + 1) * 3 + i % 4
            result.append(
                replace(
                    peak,
                    mz=mass,
                    mass_error_da=mass * error_ppm * 1e-6,
                    mass_error_ppm=error_ppm,
                    gaussian_similarity=0.98,
                )
            )
    return result


def model(lower, upper, da, ppm=0):
    return MassRegion(lower, upper, da, ppm, da * 2, ppm * 2, da * 3, ppm * 3)


def test_five_independent_regions_and_mean_plus_sample_sd():
    anchors = regional_peaks()
    scans = [
        ms2(i, rt=p.rt, mz=p.mz * (1 + 30e-6), precursor_mz_source="selected_ion")
        for i, p in enumerate(anchors)
    ]
    config, info = derive_settings(
        anchors, [0.5] * 100, scans, Config(tuning_mz_regions=5), "positive", mass_range=(100, 1100)
    )
    assert len(config.mass_regions) == 5
    for index, region in enumerate(config.mass_regions):
        selected = anchors[index * 24 : (index + 1) * 24]
        da = np.array([p.mass_error_da for p in selected])
        ppm = np.array([p.mass_error_ppm for p in selected])
        assert region.lower == 100 + index * 200
        assert region.upper == 300 + index * 200
        assert region.mz_tolerance_da == pytest.approx(da.mean() + 3 * da.std(ddof=1))
        assert region.mz_tolerance_ppm == pytest.approx(ppm.mean() + 3 * ppm.std(ddof=1))
        assert region.ms2_mz_tolerance_ppm == pytest.approx(30)
        assert info["regions"][index]["representative_count"] == 24
        assert info["regions"][index]["ms2_status"] == "tuned"
    # Relative errors transfer across regions; absolute errors stay local.
    assert config.mz_tolerance_da == 0.0005
    assert config.mz_tolerance_ppm >= max(r.mz_tolerance_ppm for r in config.mass_regions)
    assert Config.from_dict(json.loads(json.dumps(config.to_dict()))) == config


def test_sparse_region_falls_back_without_relaxing_quality_or_minimum_count():
    anchors = regional_peaks()[: 24 + 19]
    # Numerous huge, badly shaped peaks cannot rescue the sparse region.
    anchors += [replace(p, mz=p.mz + 600, height=1e12, gaussian_similarity=0.79) for p in peaks(40)]
    config, info = derive_settings(
        anchors, [1] * 30, [], Config(tuning_mz_regions=5), "positive", mass_range=(100, 1100)
    )
    assert info["regions"][0]["ms1_status"] == "tuned"
    assert info["regions"][1]["representative_count"] == 19
    for region, details in zip(config.mass_regions[1:], info["regions"][1:]):
        assert details["ms1_status"] == "file_wide_fallback"
        assert region.mz_tolerance_da == config.mz_tolerance_da
        assert region.mz_tolerance_ppm == config.mz_tolerance_ppm
        assert details["ms2_status"] == "file_wide_fallback"
        assert region.ms2_mz_tolerance_da == Config().ms2_mz_tolerance_da
        assert region.ms2_mz_tolerance_ppm == max(
            Config().ms2_mz_tolerance_ppm, region.mz_tolerance_ppm,
        )


def test_high_mass_da_error_is_not_copied_into_low_mass_fallback():
    anchors = [
        replace(p, mz=250 + i, mass_error_da=(250 + i) * 2e-6, mass_error_ppm=2)
        for i, p in enumerate(peaks())
    ]
    anchors += [
        replace(p, mz=1000 + i, mass_error_da=(1000 + i) * 10e-6, mass_error_ppm=10)
        for i, p in enumerate(peaks())
    ]
    config, info = derive_settings(anchors, [1] * 30, [], Config(), "positive")
    low, _, high = config.mass_regions
    assert info["regions"][0]["ms1_status"] == "file_wide_fallback"
    assert high.mz_tolerance_da > 0.01
    statistics = info["ms1_error_statistics"]
    expected_ppm = np.sqrt(2) * (
        statistics["mean_da_ppm"][1] + 3 * statistics["sample_sd_da_ppm"][1]
    )
    assert low.mz_tolerance_da == 0.0005 and low.mz_tolerance_ppm == expected_ppm
    assert mass_radius(150, config) == pytest.approx(150 * expected_ppm * 1e-6)
    assert mass_radius(150, config) < 0.004
    assert mass_radius(150, config, "gapfill") == mass_radius(150, config)


@pytest.mark.parametrize("polarity", ["positive", "negative", "unknown"])
def test_overrides_are_applied_to_every_region_and_gapfill_is_capped(polarity):
    base = Config(
        mz_tolerance_da=0.001,
        mz_tolerance_ppm=10,
        ms2_mz_tolerance_da=0.02,
        gapfill_mz_tolerance_da=0.04,
        tuning_mz_regions=5,
        tuning_locked=("mz_tolerance_ppm",),
    )
    config, info = derive_settings(
        regional_peaks(), [0.5] * 40, [], base, polarity, mass_range=(100, 1100)
    )
    for region in config.mass_regions:
        assert region.mz_tolerance_da == 0.001 and region.mz_tolerance_ppm == 10
        assert region.ms2_mz_tolerance_da == 0.02
        assert region.gapfill_mz_tolerance_da == 0.04
        for mass in (region.lower, (region.lower + region.upper) / 2, region.upper):
            assert _radius(mass, config) <= mass_radius(mass, config)
    if polarity == "unknown":
        assert info["status"] == "unknown_polarity"
        assert all(r["ms1_status"] == "file_wide_fallback" for r in info["regions"])


def test_lookup_edges_and_real_tracker_match_across_region_boundary():
    config = Config(
        min_scans=2, segment=False, mz_tolerance_da=0.002, mz_tolerance_ppm=10,
        mass_regions=(model(100, 200, 0.001), model(200, 300, 0.01)),
    )
    array = ms1_region_array(config)
    for mass, expected in [
        (99, 0.002),
        (100, 0.001),
        (199.9999, 0.001),
        (200, 0.01),
        (300, 0.01),
        (301, 0.00301),
    ]:
        assert mass_radius(mass, config) == pytest.approx(expected)
        assert ms1_radius(
            mass, config.mz_tolerance_da, config.mz_tolerance_ppm, array
        ) == pytest.approx(expected)
        assert _radius(mass, config) == mass_radius(mass, config)
    tracker = TraceTracker(config)
    tracker.add(np.array([199.9995]), np.array([100.0]), np.array([0.0]), 1.0)
    tracker.add(np.array([200.0004]), np.array([100.0]), np.array([0.0]), 2.0)
    assert len(tracker.flush()) == 1


def test_five_sigma_only_widens_and_one_region_is_supported():
    three, _ = derive_settings(
        regional_peaks(),
        [1] * 30,
        [],
        Config(tuning_mz_regions=5),
        "positive",
        mass_range=(100, 1100),
    )
    five, _ = derive_settings(
        regional_peaks(),
        [1] * 30,
        [],
        Config(tuning_mass_sigma=5, tuning_mz_regions=5),
        "positive",
        mass_range=(100, 1100),
    )
    for a, b in zip(three.mass_regions, five.mass_regions):
        assert b.mz_tolerance_da > a.mz_tolerance_da
        assert b.mz_tolerance_ppm > a.mz_tolerance_ppm
    single, _ = derive_settings(
        peaks(), [1] * 30, [], Config(tuning_mz_regions=1), "positive", mass_range=(50, 1000)
    )
    assert len(single.mass_regions) == 1


def test_regional_ms2_assignment_and_precursor_intensity_lookup():
    from masscubex.ms2 import MS2Collector, associate_ms2
    from masscubex.readers import Scan

    config = Config(mass_regions=(model(100, 200, 0.001), model(200, 300, 0.01)))
    features = [feature(150.004), feature(250.015)]
    records = [ms2(1, mz=150), ms2(2, mz=250)]
    assignments = associate_ms2(features, records, config, {})
    assert assignments[0].feature_index is None  # global 0.01 Da would incorrectly match
    assert assignments[1].feature_index == 1  # local 0.02 Da is needed
    collector = MS2Collector(config)
    parent = Scan(2, np.array([150.004, 250.015]), np.array([100.0, 200.0]), "positive")
    list(
        collector.ms1_stream(
            [parent, ms2(1, mz=150, intensity=None), ms2(2, mz=250, intensity=None)]
        )
    )
    assert collector.records[0].precursor_intensity is None
    assert collector.records[1].precursor_intensity == 200


def test_repeated_ms2_does_not_satisfy_distinct_anchor_minimum():
    anchors = regional_peaks()
    spectra = [
        ms2(i, rt=anchors[0].rt, mz=anchors[0].mz, precursor_mz_source="selected_ion")
        for i in range(100)
    ]
    config, info = derive_settings(
        anchors, [1] * 30, spectra, Config(), "positive", mass_range=(100, 1100)
    )
    assert info["regions"][0]["ms2_representative_count"] == 1
    assert all(r["ms2_status"] == "file_wide_fallback" for r in info["regions"])
    assert all(
        r.ms2_mz_tolerance_ppm == max(Config().ms2_mz_tolerance_ppm, r.mz_tolerance_ppm)
        for r in config.mass_regions
    )
    assert all(r.ms2_mz_tolerance_da == Config().ms2_mz_tolerance_da for r in config.mass_regions)


def test_range_windows_union_and_observed_extension(tmp_path):
    spectrum = mzml_spectrum(0, 1, [100, 900], [100, 10])
    windows = "".join(
        f'<scanWindow><cvParam accession="MS:1000501" value="{low}"/>'
        f'<cvParam accession="MS:1000500" value="{high}"/></scanWindow>'
        for low, high in [(50, 400), (500, 800)]
    )
    spectrum = spectrum.replace(
        "</scan>", "<scanWindowList>" + windows + "</scanWindowList></scan>"
    )
    path = write_mzml(tmp_path / "windows.mzML", [spectrum])
    assert next(iter_ms1(path)).acquisition_range == (50, 800)
    info = calibrate_file(path, Config()).diagnostics["polarities"]["positive"]
    assert info["acquired_mz_range"] == [50, 900]
    assert info["acquired_mz_range_source"] == "acquisition_metadata_with_observed_extension"


def test_gaussian_gate_uses_raw_baseline_corrected_shape():
    t = np.linspace(0, 40, 161)
    # A smooth, complete but strongly tailing peak passes old single-peak checks.
    y = 1e6 * np.exp(-0.5 * ((t - 10) / np.where(t < 10, 1.0, 6.0)) ** 2) + 100
    trace = np.column_stack((np.full(len(t), 200.0), y, t, np.zeros(len(t))))
    loose, reason = representative(trace, Config(tuning_min_gaussian_similarity=0))
    assert reason is None and loose.gaussian_similarity < 0.9
    strict, reason = representative(trace, Config(tuning_min_gaussian_similarity=0.9))
    assert strict is None and reason == "poor_gaussian_shape"


@pytest.mark.parametrize("format_", ["mzml", "mzxml"])
@pytest.mark.parametrize("bounds", [(50, 1500), (900, 100), ("nan", 1000), (None, None)])
def test_range_metadata_and_observed_extent(tmp_path, format_, bounds):
    if format_ == "mzml":
        scan = mzml_spectrum(0, 1, [100, 800], [100, 10])
        if bounds[0] is not None:
            window = (
                '<scanWindowList><scanWindow><cvParam accession="MS:1000501" '
                f'value="{bounds[0]}"/><cvParam accession="MS:1000500" '
                f'value="{bounds[1]}"/></scanWindow></scanWindowList>'
            )
            scan = scan.replace("</scan>", window + "</scan>")
        path = write_mzml(tmp_path / "range.mzML", [scan])
    else:
        scan = mzxml_scan(0, 1, [100, 800], [100, 10])
        if bounds[0] is not None:
            scan = scan.replace("<scan ", f'<scan startMz="{bounds[0]}" endMz="{bounds[1]}" ')
        path = write_mzxml(tmp_path / "range.mzXML", [scan])
    record = next(iter_ms1(path))
    assert record.acquisition_range == ((50, 1500) if bounds == (50, 1500) else None)
    info = calibrate_file(path, Config()).diagnostics["polarities"]["positive"]
    expected = [50, 1500] if bounds == (50, 1500) else [100, 800]
    assert info["acquired_mz_range"] == expected
    assert len(info["regions"]) == 3
    assert info["invalid_range_metadata_scans"] == int(bounds not in [(50, 1500), (None, None)])


@pytest.mark.parametrize(
    "kwargs",
    [
        dict(tuning_mz_regions=0),
        dict(tuning_mz_regions=101),
        dict(tuning_mass_sigma=0),
        dict(tuning_mass_sigma=float("nan")),
        dict(tuning_min_gaussian_similarity=1.1),
        dict(mass_regions=[]),
        dict(mass_regions=(model(100, 200, 0.001), model(250, 300, 0.001))),
    ],
)
def test_invalid_configuration(kwargs):
    with pytest.raises(ValueError):
        Config(**kwargs)
