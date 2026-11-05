"""Compiled searches for the signal-preserving MS1 consistency estimator.

No fast-math or approximate windows: strict bounds, ties and float64 mass
differences reproduce the reference NumPy decisions on float32 scan arrays.
Background moments remain in the NumPy implementation to preserve its reduction
order and threshold rounding.
"""

import numpy as np
from numba import njit

from .mass_tolerance_kernels import ms1_radius


@njit(cache=True, nogil=True)
def resolve_radii(mz, da, ppm, regions):
    radii = np.empty(len(mz), dtype=np.float64)
    for i in range(len(mz)):
        radii[i] = ms1_radius(np.float64(mz[i]), da, ppm, regions)
    return radii


@njit(cache=True, nogil=True)
def add_recurrences(targets, others, radii, counts, minimum, require_unique):
    """Merge sorted scans, preserving all three reference search boundaries.

    Bounds can move backward at a regional tolerance discontinuity. Advancing
    and rewinding each cursor handles that without assuming monotone radii.
    Completed targets need no further searches: only count >=minimum is used.
    """
    if not len(others):
        return
    nearest = lower = upper = 0
    for i in range(len(targets)):
        if counts[i] >= minimum:
            continue
        target = np.float64(targets[i])
        radius = radii[i]
        while nearest < len(others) and others[nearest] < target:
            nearest += 1
        right = min(nearest, len(others) - 1)
        left = max(nearest - 1, 0)
        distance = min(
            abs(np.float64(others[right]) - target), abs(np.float64(others[left]) - target)
        )
        if distance >= radius:
            continue
        if require_unique:
            low_value, high_value = target - radius, target + radius
            while lower > 0 and others[lower - 1] > low_value:
                lower -= 1
            while lower < len(others) and others[lower] <= low_value:
                lower += 1
            while upper > 0 and others[upper - 1] >= high_value:
                upper -= 1
            while upper < len(others) and others[upper] < high_value:
                upper += 1
            if upper - lower != 1:
                continue
        counts[i] += 1


@njit(cache=True, nogil=True)
def has_distinct(values, minimum):
    """Answer a small support query without sorting a complete population."""
    if minimum <= 0:
        return True
    distinct = np.empty(minimum, dtype=values.dtype)
    size = 0
    for value in values:
        found = False
        for j in range(size):
            if value == distinct[j]:
                found = True
                break
        if not found:
            distinct[size] = value
            size += 1
            if size == minimum:
                return True
    return False


@njit(cache=True, nogil=True)
def group_intensities(
    mz, values, order, mass_cells, width, minimum, span, minimum_cells, published_only=False
):
    """Fixed-diameter groups on the exact stable NumPy intensity order."""
    rejected = np.zeros(len(values), dtype=np.bool_)
    paper = np.zeros(len(values), dtype=np.bool_)
    lows = np.empty(len(values))
    highs = np.empty(len(values))
    sizes = np.empty(len(values), dtype=np.int64)
    cells = np.empty(max(1, minimum_cells), dtype=np.float64)
    groups = 0
    stop = len(values)
    while stop:
        high = np.float64(values[order[stop - 1]])
        low = high * (1 - width)
        # Search the sorted order without casting/copying the full spectrum.
        left, right = 0, stop
        while left < right:
            middle = (left + right) // 2
            if values[order[middle]] < low:
                left = middle + 1
            else:
                right = middle
        start = left
        count = stop - start
        enough = count >= minimum
        if count >= 4:
            for j in range(start, stop):
                paper[order[j]] = True
        if enough and not published_only:
            low_mass = high_mass = mz[order[start]]
            cell_count = 0
            for j in range(start, stop):
                mass = mz[order[j]]
                low_mass, high_mass = min(low_mass, mass), max(high_mass, mass)
                cell = mass_cells[order[j]]
                if cell_count < minimum_cells:
                    found = False
                    for k in range(cell_count):
                        if cells[k] == cell:
                            found = True
                            break
                    if not found:
                        cells[cell_count] = cell
                        cell_count += 1
            enough = high_mass - low_mass >= span and cell_count >= minimum_cells
        if enough:
            for j in range(start, stop):
                rejected[order[j]] = True
            lows[groups], highs[groups], sizes[groups] = low, high, count
            groups += 1
        stop = start
    return rejected, paper, lows[:groups], highs[:groups], sizes[:groups]


@njit(cache=True, nogil=True)
def closest_precursor(mz, mass, radius):
    if not len(mz):
        return -1
    after = np.searchsorted(mz, mass)
    if after == 0:
        best = 0
    elif after == len(mz):
        best = after - 1
    else:
        before_distance = abs(np.float64(mz[after - 1]) - mass)
        after_distance = abs(np.float64(mz[after]) - mass)
        if before_distance == after_distance:
            return -1
        best = after - 1 if before_distance < after_distance else after
    return best if abs(np.float64(mz[best]) - mass) < radius else -1


@njit(cache=True, nogil=True)
def histogram_counts(values, edges):
    """Exact NumPy explicit-edge bins, including the closed final edge."""
    counts = np.zeros(len(edges) - 1, dtype=np.int64)
    for value in values:
        if value < edges[0] or value > edges[-1]:
            continue
        index = np.searchsorted(edges, value, side="right") - 1
        if index == len(counts):
            index -= 1
        counts[index] += 1
    return counts
