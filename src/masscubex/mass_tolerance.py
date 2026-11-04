"""Learned mass windows, shared by tracing, precursor matching and gap filling.

Regions are half-open [lower, upper), except that a finite last upper bound is
included. A null last upper bound means unbounded (500+ in the default layout).
A query can match across a region boundary. Outside a finite layout, callers
use the file-wide fallback.
"""

from dataclasses import dataclass, fields
import math

def calibration_edges(config, mass_range=None):
    """Use one layout for pilot sampling and the subsequent regional mass fits."""
    if config.tuning_mz_regions is None:
        return [0.0, *config.tuning_mz_breaks, None]
    if mass_range is not None:
        import numpy as np

        low, high = mass_range
        if np.isfinite(low) and np.isfinite(high) and 0 <= low < high:
            return np.linspace(low, high, config.tuning_mz_regions + 1).tolist()
    return []


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
    import numpy as np

    return np.array(
        [
            (r.lower, np.inf if r.upper is None else r.upper, r.mz_tolerance_da, r.mz_tolerance_ppm)
            for r in config.mass_regions
        ],
        dtype=np.float64,
    ).reshape(-1, 4)


def __getattr__(name):
    # Config/CLI discovery use only lightweight tolerance definitions. Keep the
    # numerical kernel available at its existing import path without loading
    # NumPy, Numba and LLVM into the batch coordinator.
    if name == "ms1_radius":
        from .mass_tolerance_kernels import ms1_radius

        return ms1_radius
    raise AttributeError(f"module {__name__!r} has no attribute {name!r}")
