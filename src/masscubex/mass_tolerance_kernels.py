"""Numerical mass-window kernel; imported only by extraction workers."""

from numba import njit


@njit(cache=True, inline="always")
def ms1_radius(mass, da, ppm, regions):
    for index in range(len(regions)):
        if regions[index, 0] <= mass and (
            mass < regions[index, 1] or index == len(regions) - 1 and mass == regions[index, 1]
        ):
            da, ppm = regions[index, 2], regions[index, 3]
            break
    return max(da, mass * ppm * 1e-6)
