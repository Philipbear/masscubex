"""Sparse query discovery must preserve the scalar NumPy membership rules."""

import numpy as np
import pytest

from masscubex import Config
from masscubex.adducts import catalog, selected_adducts
from masscubex.grouping import isotope_targets
from masscubex.grouping_queries import anchor_snapshot, static_query_index, target_tables


@pytest.mark.parametrize("ppm", [0.0, 20.0])
def test_static_index_matches_all_scalar_candidates_and_priorities(ppm):
    config = Config(group_mz_tolerance_ppm=ppm)
    rng = np.random.default_rng(122)
    masses = np.r_[rng.uniform(100, 130, 64), np.linspace(100, 100.02, 30)].astype(np.float32)
    polarity = rng.choice(np.array([-1, 1], dtype=np.int8), len(masses))
    observed = rng.random(len(masses)) > 0.1
    apices = rng.uniform(0, 0.1, len(masses))
    priority = rng.permutation(len(masses)).astype(np.int64)
    rules = selected_adducts(config.adduct_ids)
    defaults = {"positive": catalog()[1], "negative": catalog()[5]}
    by_polarity = {
        p: [r for r in rules if r.polarity == p and r.name != defaults[p].name] for p in defaults
    }
    tables = target_tables(masses, by_polarity, defaults, config, set(defaults))
    if int(np.__version__.split(".")[0]) != 1:
        assert tables == {}  # The general scalar path handles NumPy 2.x.
        return
    ids = np.argsort(masses, kind="stable")
    for p, table in tables.items():
        code = 1 if p == "positive" else -1
        result, offsets = static_query_index(
            masses[ids].astype(float),
            ids,
            masses,
            table,
            polarity,
            observed,
            apices,
            code,
            0.05,
            priority,
            config.group_mz_tolerance_da,
            ppm,
        )
        assert result.dtype == offsets.dtype == np.dtype("int32")
        for anchor, mass in enumerate(masses):
            neutral = defaults[p].neutral_mass(mass)
            targets = [t for _, t in isotope_targets(mass)] + [
                float(r.mz(neutral)) for r in by_polarity[p]
            ]
            np.testing.assert_array_equal(table[anchor], targets)
            for q, target in enumerate(targets):
                key = anchor * table.shape[1] + q
                expected = np.empty(0, dtype=np.int64)
                if polarity[anchor] == code:
                    tolerance = max(config.group_mz_tolerance_da, target * ppm * 1e-6)
                    mask = (
                        (np.abs(masses - target) < tolerance)
                        & (polarity == code)
                        & observed
                        & (np.abs(apices - apices[anchor]) < 0.05)
                    )
                    expected = np.flatnonzero(mask)
                    expected = expected[np.argsort(priority[expected])]
                np.testing.assert_array_equal(result[offsets[key] : offsets[key + 1]], expected)
        empty, no_offsets = static_query_index(
            masses[ids].astype(float),
            ids,
            masses,
            table,
            polarity,
            observed,
            apices,
            code,
            0.05,
            priority,
            0.01,
            ppm,
            budget=0,
        )
        assert not len(empty) and not len(no_offsets)


def test_anchor_snapshot_keeps_priority_and_original_eligibility():
    packed = np.array([3, 1, 0, 2, 3, 2, 1], dtype=np.int32)
    offsets = np.array([0, 0, 0, 0, 3, 5, 7], dtype=np.int32)
    assigned = np.array([4, -1, -1, -1])
    result, boundaries = anchor_snapshot(packed, offsets, assigned, 1, 3)
    np.testing.assert_array_equal(result, [3, 1, 2, 3, 2, 1])
    np.testing.assert_array_equal(boundaries, [0, 2, 4, 6])
    # A child assigned during the anchor remains in its later snapshot queries.
    assigned[3] = 5
    np.testing.assert_array_equal(result[boundaries[1]:boundaries[2]], [2, 3])


def test_sparse_storage_growth_respects_exact_capacity_and_budget():
    n = 128
    masses = np.full(n, 100.0, dtype=np.float32)
    ids = np.arange(n, dtype=np.int64)
    targets = np.full((n, 3), 100.0)
    arguments = (
        masses.astype(float), ids, masses, targets, np.ones(n, dtype=np.int8),
        np.ones(n, dtype=bool), np.zeros(n), 1, 0.05, ids[::-1].copy(), 0.01, 0.0,
    )
    entries = n * n * 3
    result, offsets = static_query_index(*arguments, budget=entries)
    np.testing.assert_array_equal(result.reshape(-1, n), np.tile(ids[::-1], (n * 3, 1)))
    np.testing.assert_array_equal(offsets, np.arange(n * 3 + 1) * n)
    result, offsets = static_query_index(*arguments, budget=entries - 1)
    assert not len(result) and not len(offsets)


def test_prepared_eic_queries_preserve_float32_mass_boundaries(tmp_path):
    from masscubex.grouping_spectra import GroupingSpectra

    mass = np.float32(300.0)
    neighbors = np.array([
        np.nextafter(mass, np.float32(0)), mass,
        np.nextafter(mass, np.float32(np.inf)),
    ])
    intensities = np.array([100.0, 200.0, 300.0], dtype=np.float32)
    times = np.array([3.0, 1.0, 2.0, 4.0, 2.0])
    with GroupingSpectra(tmp_path) as reference, GroupingSpectra(tmp_path, dtype="float32") as compact:
        for spectrum in [reference, compact]:
            spectrum.append(1.0, "positive", neighbors, intensities)
            spectrum.append(2.0, "positive", neighbors, intensities * 2)
            spectrum.append(3.0, "positive", neighbors, intensities * 3)
        for radius in [0.0, float(neighbors[2] - mass), 0.005]:
            query = compact.prepare_query(times, "positive", radius)
            for target in [float(mass), float(neighbors[2]), (float(mass) + float(neighbors[2])) / 2]:
                np.testing.assert_array_equal(
                    query(target), reference.eic(target, times, "positive", radius)
                )
        assert compact.file.seek(0, 2) * 2 == reference.file.seek(0, 2)
        np.testing.assert_array_equal(compact.prepare_query(times, "negative", 0.01)(300.0), 0)
