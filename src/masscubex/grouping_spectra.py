"""File-backed processed MS1 spectra for upstream grouping EIC queries.

Keep only scan offsets/times in memory. The anonymous temporary file is closed
on every exit; partner queries can extend beyond a detected trace's boundaries.
"""

import tempfile
from functools import partial

import numpy as np
from numba import njit


@njit(cache=True, nogil=True)
def extract_eic(points, starts, stops, scan_ids, mass, radius):
    values = np.zeros(len(scan_ids))
    for i, scan in enumerate(scan_ids):
        if scan < 0:
            continue
        first, last = starts[scan], stops[scan]
        masses = points[first:last, 0]
        left = first + np.searchsorted(masses, mass - radius)
        right = first + np.searchsorted(masses, mass + radius, side="right")
        for j in range(left, right):
            values[i] = max(values[i], points[j, 1])
    return values


class GroupingSpectra:
    def __init__(self, directory=None, *, dtype=np.float64):
        self.file = tempfile.TemporaryFile(dir=directory)
        self.dtype = np.dtype(dtype)
        self.scans = {}
        self.count = 0
        self.points = None

    def append(self, time, polarity, masses, intensities):
        if self.points is not None:
            raise RuntimeError("Cannot append after grouping spectra are finalized")
        self.scans.setdefault(polarity, []).append((time, self.count, self.count + len(masses)))
        self.file.write(np.column_stack((masses, intensities)).astype(self.dtype, copy=False).tobytes())
        self.count += len(masses)

    def finish(self):
        if self.points is None:
            self.file.flush()
            self.points = (
                np.memmap(self.file, dtype=self.dtype, mode="r", shape=(self.count, 2))
                if self.count
                else np.empty((0, 2))
            )
            self.scans = {
                polarity: (
                    np.array([s[0] for s in rows]),
                    np.array([s[1] for s in rows], dtype=np.int64),
                    np.array([s[2] for s in rows], dtype=np.int64),
                )
                for polarity, rows in self.scans.items()
            }

    def prepare_query(self, times, polarity, radius):
        """Join an anchor's time grid once, then query any partner mass."""
        self.finish()
        if polarity not in self.scans:
            return lambda mass: np.zeros(len(times))
        acquired, starts, stops = self.scans[polarity]
        indices = np.searchsorted(acquired, times)
        valid = indices < len(acquired)
        valid[valid] &= acquired[indices[valid]] == times[valid]
        indices[~valid] = -1
        return partial(extract_eic, self.points, starts, stops, indices, radius=radius)

    def eic(self, mass, times, polarity, radius):
        return self.prepare_query(times, polarity, radius)(mass)

    def close(self):
        if isinstance(self.points, np.memmap):
            self.points._mmap.close()
        self.points = None
        self.file.close()

    def __enter__(self):
        return self

    def __exit__(self, *exc):
        self.close()
