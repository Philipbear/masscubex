"""Learned mass windows, shared by tracing, precursor matching and gap filling.

Regions are half-open [lower, upper), except that a finite last upper bound is
included. A null last upper bound means unbounded (500+ in the default layout).
A query can match across a region boundary. Outside a finite layout, callers
use the file-wide fallback.
"""

from dataclasses import dataclass, fields
import math

import numpy as np
from numba import njit


@dataclass(frozen=True, slots=True)
class MassRegion:
    lower: float
    upper: float | None
    mz_tolerance_da: float
    mz_tolerance_ppm: float
    ms2_mz_tolerance_da: float
    ms2_mz_tolerance_ppm: float
    gapfill_mz_tolerance_da: float
    gapfill_mz_tolerance_ppm: float

    def __post_init__(self):
        for field in fields(self):
            value = getattr(self, field.name)
            if field.name == "upper" and value is None:
                continue
            if (
                isinstance(value, bool)
                or not isinstance(value, (int, float))
                or not math.isfinite(value)
                or value < 0
            ):
                raise ValueError(f"Invalid mass region {field.name}: {value}")
        if self.upper is not None and self.upper <= self.lower:
            raise ValueError("Mass region upper bound must exceed its lower bound")
        for prefix in ("", "ms2_", "gapfill_"):
            if not (
                getattr(self, prefix + "mz_tolerance_da")
                or getattr(self, prefix + "mz_tolerance_ppm")
            ):
                raise ValueError("Mass region requires a positive mass tolerance")


def mass_radius(mass, config, kind="ms1"):
    """Resolve one window at its target mass; ion grouping has separate settings."""
    prefix = {"ms1": "", "ms2": "ms2_", "gapfill": "gapfill_"}[kind]
    settings = config
    for index, region in enumerate(config.mass_regions):
        if region.lower <= mass and (
            region.upper is None
            or mass < region.upper
            or index == len(config.mass_regions) - 1
            and mass == region.upper
        ):
            settings = region
            break
    return max(
        getattr(settings, prefix + "mz_tolerance_da"),
        mass * getattr(settings, prefix + "mz_tolerance_ppm") * 1e-6,
    )


def ms1_region_array(config):
    """Pack once per tracker, avoiding Python calls in the centroid hot loop."""
    return np.array(
        [
            (r.lower, np.inf if r.upper is None else r.upper, r.mz_tolerance_da, r.mz_tolerance_ppm)
            for r in config.mass_regions
        ],
        dtype=np.float64,
    ).reshape(-1, 4)


@njit(cache=True, inline="always")
def ms1_radius(mass, da, ppm, regions):
    for index in range(len(regions)):
        if regions[index, 0] <= mass and (
            mass < regions[index, 1] or index == len(regions) - 1 and mass == regions[index, 1]
        ):
            da, ppm = regions[index, 2], regions[index, 3]
            break
    return max(da, mass * ppm * 1e-6)
