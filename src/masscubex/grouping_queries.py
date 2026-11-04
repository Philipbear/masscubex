"""Batch the independent mass queries in one greedy grouping anchor.

Only candidate discovery is batched. Assignments, correlations, overwrites and
child-isotope checks retain their original order in grouping.py.
"""

import numpy as np
from numba import njit


@njit(cache=True, nogil=True)
def anchor_snapshot(packed, offsets, assigned, anchor, width):
    """Filter one anchor's cached queries against its initial assignment mask.

    The greedy algorithm keeps this mask fixed while processing an anchor.
    Taking it once also preserves eligibility when a later ion rule overwrites
    an assignment made by an earlier rule for the same anchor.
    """
    first = anchor * width
    start, stop = offsets[first], offsets[first + width]
    selected = np.empty(stop - start, dtype=packed.dtype)
    boundaries = np.empty(width + 1, dtype=np.int64)
    boundaries[0] = 0
    count = 0
    for query in range(width):
        for i in range(offsets[first + query], offsets[first + query + 1]):
            candidate = packed[i]
            if assigned[candidate] < 0:
                selected[count] = candidate
                count += 1
        boundaries[query + 1] = count
    return selected[:count], boundaries


@njit(cache=True, nogil=True)
def static_query_index(
    sorted_masses,
    sorted_ids,
    masses,
    targets,
    polarities,
    observed,
    apices,
    polarity,
    rt_radius,
    priority,
    da,
    ppm,
    budget=8_000_000,
):
    """Build a bounded CSR index of immutable mass/RT candidates.

    Group assignment is deliberately absent: the caller checks its current
    assignment generation when consuming a row. This caches discovery without
    changing any greedy decisions or their order.
    """
    rows, width = targets.shape
    # The bounded index fits signed 32-bit IDs and offsets. Use the general
    # query path when an exceptional caller requests a larger representation.
    if rows * width + 1 > budget or rows > 2**31 - 1 or budget > 2**31 - 1:
        return np.empty(0, dtype=np.int32), np.empty(0, dtype=np.int32)
    offsets = np.zeros(rows * width + 1, dtype=np.int32)
    result = np.empty(min(max(rows, 1), budget, 65536), dtype=np.int32)
    count = 0
    # Grow bounded packed storage rather than repeating every binary search in
    # a count/fill pair of passes. Capacity, including unused entries, remains
    # within the same budget. No dense rows-by-feature candidate matrix exists.
    for anchor in range(rows):
        for query in range(width):
            key = anchor * width + query
            if polarities[anchor] != polarity:
                offsets[key + 1] = count
                continue
            target = targets[anchor, query]
            tolerance = max(da, target * ppm * 1e-6)
            padding = 2.0 * 1.1920928955078125e-7 * abs(target)
            first = np.searchsorted(sorted_masses, target - tolerance - padding)
            stop = np.searchsorted(sorted_masses, target + tolerance + padding, side="right")
            target32, limit32 = np.float32(target), np.float32(tolerance)
            start = count
            for j in range(first, stop):
                candidate = sorted_ids[j]
                if (
                    abs(np.float32(masses[candidate] - target32)) < limit32
                    and polarities[candidate] == polarity
                    and observed[candidate]
                    and abs(apices[candidate] - apices[anchor]) < rt_radius
                ):
                    if count == len(result):
                        if count == budget:
                            return np.empty(0, dtype=np.int32), np.empty(0, dtype=np.int32)
                        grown = np.empty(min(budget, len(result) * 2), dtype=np.int32)
                        grown[:count] = result
                        result = grown
                    result[count] = candidate
                    count += 1
            offsets[key + 1] = count
            if count - start > 1:
                selected = result[start:count].copy()
                result[start:count] = selected[np.argsort(priority[selected])]
    return result[:count], offsets


def target_tables(masses, rules, defaults, config, present):
    """Precompute scalar-equivalent targets on the supported NumPy 1.x path.

    NumPy 2.x scalar promotion differs; it retains the general query path.
    Bound memory and the scalar-promotion range rather than changing semantics
    for unusual catalogs or extreme externally supplied values.
    """
    if (
        int(np.__version__.split(".")[0]) != 1
        or type(config.group_mz_tolerance_da) is not float
        or config.group_mz_tolerance_da >= 1e12
        or not np.isfinite(masses).all()
        or np.max(np.abs(masses), initial=0) >= 1e20
        or sum(3 + len(rules[p]) for p in present if p in rules) * len(masses) * 8 > 64 * 1024**2
    ):
        return {}
    values = masses.astype(np.float64)
    tables = {}
    for polarity in present:
        if polarity not in rules:
            continue
        base = defaults[polarity]
        neutral = (abs(base.charge) * values - base.mass_shift) / base.molecules
        targets = np.empty((len(masses), 3 + len(rules[polarity])))
        for i in range(3):
            targets[:, i] = values + 1 + i * 1.003355
        for i, rule in enumerate(rules[polarity], 3):
            targets[:, i] = (rule.molecules * neutral + rule.mass_shift) / abs(rule.charge)
        tables[polarity] = targets
    return tables


@njit(cache=True, nogil=True)
def snapshot_queries(
    sorted_masses,
    sorted_ids,
    masses,
    targets,
    limits,
    lower,
    upper,
    assigned,
    polarities,
    observed,
    apices,
    anchor,
    rt_radius,
    priority,
):
    """Discover all initial candidates against the same eligibility snapshot."""
    starts = np.searchsorted(sorted_masses, lower)
    stops = np.searchsorted(sorted_masses, upper, side="right")
    capacity = np.sum(stops - starts)
    # Avoid a large temporary for an unusually broad configured query.
    if capacity > max(65536, len(masses)):
        return np.empty(0, dtype=np.int64), np.empty(0, dtype=np.int64)
    result = np.empty(capacity, dtype=np.int64)
    offsets = np.empty(len(targets) + 1, dtype=np.int64)
    count = 0
    offsets[0] = 0
    for query in range(len(targets)):
        start = count
        target, limit = np.float32(targets[query]), np.float32(limits[query])
        for i in range(starts[query], stops[query]):
            candidate = sorted_ids[i]
            distance = abs(np.float32(masses[candidate] - target))
            if (
                distance < limit
                and assigned[candidate] < 0
                and polarities[candidate] == polarities[anchor]
                and observed[candidate]
                and abs(apices[candidate] - apices[anchor]) < rt_radius
            ):
                result[count] = candidate
                count += 1
        if count - start > 1:
            selected = result[start:count].copy()
            result[start:count] = selected[np.argsort(priority[selected])]
        offsets[query + 1] = count
    return result[:count], offsets
