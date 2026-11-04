"""Numerical equivalence checks for grouping's allocation shortcuts."""

import numpy as np
import pytest

from masscubex.grouping import (
    _eligible_candidates,
    _matching_candidates,
    _positive_correlation,
    trace_correlation,
)


@pytest.mark.parametrize("length", [1, 2, 3, 7, 64, 1001])
def test_aligned_correlation_is_bitwise_equal_to_time_join(length):
    rng = np.random.default_rng(length)
    times = np.cumsum(rng.uniform(0.1, 2, size=length))
    x = rng.lognormal(8, 3, size=length)
    for y in (x, rng.lognormal(8, 3, size=length), np.ones(length)):
        # Real traces contain missing observations; positivity filtering must
        # stay identical before the unchanged NumPy Pearson calculation.
        y = y.copy()
        y[::4] = 0
        expected = trace_correlation(np.column_stack((x, times)), np.column_stack((y, times)))
        actual = _positive_correlation(x, y)
        assert actual[0] == expected[0]
        if np.isnan(expected[1]):
            assert np.isnan(actual[1])
        else:
            assert actual[1] == expected[1]


@pytest.mark.parametrize("snapshot", [False, True])
@pytest.mark.parametrize("radius", [0.0, 0.05, 0.1, 1.0 / 3])
def test_compiled_candidate_filter_preserves_boundary_snapshot_and_priority(snapshot, radius):
    # Include both signs of the exact boundary and their immediate neighbors.
    # The Python subtraction below, rather than decimal intuition, is the rule.
    apices = np.array(
        [0.0, radius, -radius, np.nextafter(radius, 0), np.nextafter(radius, np.inf)] * 6,
        dtype=np.float64,
    )
    count = len(apices)
    ids = np.random.default_rng(3).permutation(count).astype(np.int64)
    assigned_epoch = np.tile(np.array([-1, 5, 4, 6, -1], dtype=np.int64), 6)
    polarities = np.tile(np.array([1, 1, 1, -1, 0, 1], dtype=np.int8), 5)
    observed = np.array([i % 7 != 0 for i in range(count)])
    # Priorities are unique but arbitrary. Equal RT distances must therefore
    # retain this order when the caller selects only its nearest RT candidate.
    priority = np.random.default_rng(7).permutation(count).astype(np.int64)
    epoch, anchor = 5, 0
    expected = sorted(
        [
            int(j)
            for j in ids
            if (assigned_epoch[j] < 0 or (snapshot and assigned_epoch[j] == epoch))
            and polarities[j] == polarities[anchor]
            and observed[j]
            and abs(apices[j] - apices[anchor]) < radius
        ],
        key=lambda j: priority[j],
    )
    actual = _eligible_candidates(
        ids, assigned_epoch, epoch, snapshot, polarities, observed, apices, anchor, radius, priority
    ).tolist()
    assert actual == expected
    if actual:

        def distance(j):
            return abs(apices[j] - apices[anchor])

        assert min(actual, key=distance) == min(expected, key=distance)
    assert (
        _eligible_candidates(
            ids[:0],
            assigned_epoch,
            epoch,
            snapshot,
            polarities,
            observed,
            apices,
            anchor,
            radius,
            priority,
        ).size
        == 0
    )


@pytest.mark.parametrize("target_type", [float, np.float32, np.float64])
@pytest.mark.parametrize("tolerance_type", [float, np.float32, np.float64])
def test_compiled_mass_query_matches_numpy_scalar_promotion(target_type, tolerance_type):
    for center, width in [
        (1e-35, 1e-38),
        (0.1, 0.01),
        (100.1234567, 0.01),
        (500.0, 0.0005),
        (10000.0, 0.1),
        (1e30, 1e25),
        (1e38, 1e37),
        (1e39, 1e39),
    ]:
        with np.errstate(over="ignore"):
            target, tolerance = target_type(center), tolerance_type(width)
            bounds = np.array([center - width, center, center + width], dtype=np.float32)
        if not np.isfinite(target) or not np.isfinite(tolerance):
            # These values retain the NumPy fallback, not the compiled query.
            continue
        masses = np.unique(
            np.r_[
                bounds,
                np.nextafter(bounds, -np.inf, dtype=np.float32),
                np.nextafter(bounds, np.inf, dtype=np.float32),
                np.float32(0),
            ]
        )
        masses = masses[np.isfinite(masses)]
        sorted_ids = np.arange(len(masses), dtype=np.int64)
        sorted_masses = masses.astype(np.float64)
        padding = 2 * np.finfo(np.float32).eps * abs(target)
        lower, upper = target - tolerance - padding, target + tolerance + padding
        first, stop = sorted_masses.searchsorted(lower), sorted_masses.searchsorted(upper, "right")
        candidates = sorted_ids[first:stop]
        with np.errstate(over="ignore", invalid="ignore"):
            expected = candidates[np.abs(masses[candidates] - target) < tolerance]
        # Use reverse priority to also check sorting after mass filtering.
        priority = np.arange(len(masses), dtype=np.int64)[::-1].copy()
        expected = expected[np.argsort(priority[expected])]
        subtract_dtype = np.result_type(masses, target)
        compare_dtype = np.result_type(np.empty(0, dtype=subtract_dtype), tolerance)
        actual = _matching_candidates(
            sorted_masses,
            sorted_ids,
            masses,
            float(lower),
            float(upper),
            float(target),
            float(tolerance),
            subtract_dtype == np.float32,
            compare_dtype == np.float32,
            np.full(len(masses), -1, dtype=np.int64),
            1,
            True,
            np.ones(len(masses), dtype=np.int8),
            np.ones(len(masses), dtype=bool),
            np.zeros(len(masses)),
            0,
            0.05,
            priority,
        )
        np.testing.assert_array_equal(actual, expected)


@pytest.mark.parametrize(
    "subtract_float32,compare_float32", [(True, True), (True, False), (False, False)]
)
def test_mass_query_reproduces_each_supported_promotion_combination(
    subtract_float32, compare_float32
):
    # Exercise strong float64 promotion even on NumPy 1.x, where scalar values
    # in this range ordinarily permit the float32 array operation.
    target = 100.123456789
    tolerance = 0.0100000000001
    masses = np.linspace(100.1134, 100.1335, 2001).astype(np.float32)
    dtype = np.float32 if subtract_float32 else np.float64
    distances = np.abs(masses.astype(dtype) - dtype(target))
    comparison_dtype = np.float32 if compare_float32 else np.float64
    expected = np.flatnonzero(distances.astype(comparison_dtype) < comparison_dtype(tolerance))
    ids = np.arange(len(masses), dtype=np.int64)
    actual = _matching_candidates(
        masses.astype(np.float64),
        ids,
        masses,
        100.0,
        101.0,
        target,
        tolerance,
        subtract_float32,
        compare_float32,
        np.full(len(masses), -1, dtype=np.int64),
        1,
        True,
        np.ones(len(masses), dtype=np.int8),
        np.ones(len(masses), dtype=bool),
        np.zeros(len(masses)),
        0,
        0.05,
        ids,
    )
    np.testing.assert_array_equal(actual, expected)
