"""Bounded mass-stratified sampling for the calibration pilot.

Fixed boundaries need no range-discovery pass. Optional equal-width layouts
need the complete per-polarity range first, so a later scan cannot move the
sampling boundaries used for earlier scans. ScanSource can replay that pass.
"""

from contextlib import closing
from dataclasses import dataclass
import math

import numpy as np

from .readers import MS2Scan, iter_scans


@dataclass
class AcquisitionRange:
    """Track acquisition metadata and raw observed extents before any filtering."""

    observed_low: float = math.inf
    observed_high: float = -math.inf
    acquired_low: float = math.inf
    acquired_high: float = -math.inf
    metadata_scans: int = 0
    invalid_scans: int = 0

    def add(self, scan):
        if len(scan.mz):
            self.observed_low = min(self.observed_low, float(scan.mz[0]))
            self.observed_high = max(self.observed_high, float(scan.mz[-1]))
        if scan.acquisition_range is not None:
            self.metadata_scans += 1
            self.acquired_low = min(self.acquired_low, scan.acquisition_range[0])
            self.acquired_high = max(self.acquired_high, scan.acquisition_range[1])
        self.invalid_scans += int(scan.invalid_acquisition_range)

    @property
    def bounds(self):
        return min(self.observed_low, self.acquired_low), max(
            self.observed_high, self.acquired_high
        )

    def diagnostics(self):
        low, high = self.bounds
        source = "acquisition_metadata" if self.metadata_scans else "observed_ms1_extent"
        if self.metadata_scans and (low < self.acquired_low or high > self.acquired_high):
            source = "acquisition_metadata_with_observed_extension"
        return {
            "acquired_mz_range": [low, high]
            if math.isfinite(low) and math.isfinite(high)
            else None,
            "acquired_mz_range_source": source,
            "range_metadata_scans": self.metadata_scans,
            "invalid_range_metadata_scans": self.invalid_scans,
        }


def discover_ranges(source, config, scan_source=None):
    """Read all MS1 extents by polarity; do not retain arrays or alter MS2 records."""
    ranges = {}
    with closing(
        scan_source.scans()
        if scan_source is not None
        else iter_scans(source, assume_centroid=config.assume_centroid, include_ms2=False)
    ) as scans:
        for scan in scans:
            if not isinstance(scan, MS2Scan):
                ranges.setdefault(scan.polarity, AcquisitionRange()).add(scan)
    return ranges


def strongest_per_region(mz, intensity, breaks, maximum):
    """Select up to maximum per half-open region, returning sorted original indices.

    Use partial selection instead of sorting every intensity. Equal-intensity
    cutoff ties prefer lower original mass-order indices, including when every
    point ties. Exact boundary masses enter the higher region. Sparse regions
    keep all their points; unused capacity is not borrowed by another region.
    """
    cuts = np.r_[0, np.searchsorted(mz, breaks, side="left"), len(mz)]
    selected = []
    counts = np.zeros(len(cuts) - 1, dtype=np.int64)
    for region, (start, stop) in enumerate(zip(cuts, cuts[1:])):
        values = intensity[start:stop]
        if len(values) <= maximum:
            indices = np.arange(len(values))
        else:
            cutoff = np.partition(values, len(values) - maximum)[len(values) - maximum]
            higher = np.flatnonzero(values > cutoff)
            tied = np.flatnonzero(values == cutoff)[: maximum - len(higher)]
            indices = np.sort(np.r_[higher, tied])
        selected.append(start + indices)
        counts[region] = len(indices)
    return np.concatenate(selected), counts
