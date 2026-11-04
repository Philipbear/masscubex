"""Bounded replay of decoded scans within a single file extraction.

Small files need one XML/decompression pass for calibration, extraction and raw
MS1 gap filling. On overflow we discard the cache and continue streaming. Large
profile files therefore cannot cause an unbounded whole-file allocation.
"""

from contextlib import closing
from dataclasses import fields, replace
import sys

import numpy as np

from .readers import MS2Scan, iter_scans


def _retained_size(scan):
    """Conservative charge for arrays, scan fields and the cache-list reference."""
    size = sys.getsizeof(scan) + 16
    for field in fields(scan):
        value = getattr(scan, field.name)
        size += sys.getsizeof(value)
        if isinstance(value, np.ndarray):
            # A view may keep a larger decoding buffer alive. Charge its owner;
            # shared owners can be counted twice, deliberately conservatively.
            owner = value
            while isinstance(owner.base, np.ndarray):
                owner = owner.base
            size += owner.nbytes
    return size


class ScanSource:
    """One reader configuration; replay MS2 metadata independently for each pass.

    Consumers treat raw MS1 arrays as read-only. MS2 consumers replace metadata
    and pruned arrays, so every pass gets a fresh shallow copy of each MS2 record.
    A partial or failed first pass never becomes a complete replayable source.
    """

    def __init__(self, path, config):
        self.path = path
        self.options = dict(
            assume_centroid=config.assume_centroid,
            include_ms2=config.associate_ms2,
            decode_ms2=config.associate_ms2,
            sum_ms2=config.best_ms2 == "fragment_intensity",
        )
        self.budget = int(config.scan_cache_mb * 1024**2)
        self._scans = []
        self._complete = False
        self._attempted = False
        self.stats = dict(budget_bytes=self.budget, peak_bytes=0, reader_passes=0, replay_passes=0)

    def scans(self, *, ms1_only=False):
        if self._complete:
            self.stats["replay_passes"] += 1
            for scan in self._scans:
                if not ms1_only or not isinstance(scan, MS2Scan):
                    yield replace(scan) if isinstance(scan, MS2Scan) else scan
            return
        retaining = not ms1_only and not self._attempted and self.budget > 0
        self._attempted = True
        self.stats["reader_passes"] += 1
        retained = 0
        try:
            options = {**self.options, "include_ms2": False} if ms1_only else self.options
            with closing(iter_scans(self.path, **options)) as scans:
                for scan in scans:
                    if retaining:
                        size = _retained_size(scan)
                        if retained + size <= self.budget:
                            self._scans.append(replace(scan) if isinstance(scan, MS2Scan) else scan)
                            retained += size
                            self.stats["peak_bytes"] = retained
                        else:
                            self._scans.clear()
                            retaining = False
                    yield scan
            self._complete = retaining
        finally:
            if not self._complete:
                self._scans.clear()

    def ms1(self):
        yield from self.scans(ms1_only=True)

    def release(self):
        """Release decoded arrays once all numerical passes have finished."""
        self._scans.clear()
        self._complete = False
