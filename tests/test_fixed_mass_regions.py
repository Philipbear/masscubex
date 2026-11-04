"""Fixed calibration bins must retain exact boundaries and an unbounded tail."""

from dataclasses import replace
import json

import numpy as np
import pytest

from masscubex import Config
from masscubex.cli import main
from masscubex.mass_tolerance import MassRegion, mass_radius, ms1_radius, ms1_region_array
from masscubex.tuning import derive_settings
from test_tuning import peaks


def test_fixed_default_ranges_are_independent_of_observed_extent():
    for span in (None, (100, 1500), (70, 350), (650, 2000)):
        config, info = derive_settings(
            peaks(), [0.5] * 50, [], Config(), "positive", mass_range=span
        )
        assert [(r.lower, r.upper) for r in config.mass_regions] == [
            (0, 200),
            (200, 500),
            (500, None),
        ]
        assert info["regional_layout"] == "fixed_breaks"
        encoded = json.dumps(config.to_dict(), allow_nan=False)
        assert Config.from_dict(json.loads(encoded)) == config


def test_fixed_regions_learn_independently_with_known_scatter():
    anchors = []
    for index in range(4):
        for i, p in enumerate(peaks()):
            mz = 20 + 200 * index + 7 * i
            ppm = 6 + index * 3 + i % 3
            anchors.append(replace(p, mz=mz, mass_error_da=mz * ppm * 1e-6, mass_error_ppm=ppm))
    config, info = derive_settings(
        anchors,
        [1] * 50,
        [],
        Config(tuning_mz_breaks=(200, 400, 600)),
        "positive",
        mass_range=(20, 1000),
    )
    assert all(r["ms1_status"] == "tuned" for r in info["regions"])
    assert [r["representative_count"] for r in info["regions"]] == [24] * 4
    for index, region in enumerate(config.mass_regions):
        errors = np.array([p.mass_error_ppm for p in anchors[index * 24 : (index + 1) * 24]])
        assert region.mz_tolerance_ppm == pytest.approx(errors.mean() + 3 * errors.std(ddof=1))


def test_half_open_boundaries_and_unbounded_tail_match_in_both_lookups():
    regions = tuple(
        MassRegion(lower, upper, da, 0, da, 0, da, 0)
        for (lower, upper), da in zip(
            [(0, 200), (200, 400), (400, 600), (600, None)], [0.001, 0.002, 0.003, 0.004]
        )
    )
    config = Config(mass_regions=regions)
    packed = ms1_region_array(config)
    for mass, expected in [
        (199.999, 0.001),
        (200, 0.002),
        (400, 0.003),
        (600, 0.004),
        (1e6, 0.004),
    ]:
        assert mass_radius(mass, config) == expected
        assert mass_radius(mass, config, "ms2") == expected
        assert ms1_radius(mass, 0.01, 20.0, packed) == expected


def test_fixed_boundary_customization_and_legacy_equal_width_cli(capsys):
    main(["extract", "--print-config", "--tuning-mz-breaks", "150", "300", "900"])
    values = json.loads(capsys.readouterr().out)
    assert values["tuning_mz_breaks"] == [150, 300, 900] and values["tuning_mz_regions"] is None
    main(["extract", "--print-config", "--tuning-mz-regions", "5"])
    values = json.loads(capsys.readouterr().out)
    assert values["tuning_mz_regions"] == 5
    with pytest.raises(SystemExit):
        main(["extract", "--print-config", "--tuning-mz-regions", "5", "--tuning-mz-breaks", "100"])


def test_empty_breaks_means_one_global_region():
    config, _ = derive_settings(peaks(), [1] * 20, [], Config(tuning_mz_breaks=()), "positive")
    assert len(config.mass_regions) == 1 and config.mass_regions[0].upper is None


@pytest.mark.parametrize(
    "breaks",
    [
        (0, 200),
        (200, 200),
        (400, 200),
        (200, float("inf")),
        (float("nan"),),
        (True, 200),
        [200, 400],
    ],
)
def test_invalid_breaks(breaks):
    with pytest.raises(ValueError):
        Config(tuning_mz_breaks=breaks)


def test_unbounded_region_cannot_precede_another_region():
    a = MassRegion(0, None, 0.001, 5, 0.001, 5, 0.001, 5)
    b = replace(a, lower=200)
    with pytest.raises(ValueError):
        Config(mass_regions=(a, b))
