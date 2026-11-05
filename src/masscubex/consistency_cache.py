"""Per-file disk-backed centroids and decisions for consistency extraction.

Raw profile samples are mass-detected one scan at a time. Persistent centroid
arrays and point decisions use temporary memory maps; background-training
copies stay bounded by the estimator's short temporal pool. Independent files
never share mutable preprocessing state.
"""

from collections import Counter
from dataclasses import asdict, replace
from pathlib import Path
from tempfile import TemporaryDirectory
import time

import numpy as np

from . import consistency_noise as MODEL
from . import profile_noise as PROFILE
from .upstream_adapter import merge_centroids
from .tuning import config_for


class ConsistencyCache:
    def __init__(self, directory):
        self._temporary = TemporaryDirectory(prefix=".masscubex-noise-", dir=directory)
        self.directory = Path(self._temporary.name)
        self.frames = []
        self.decisions = []
        self._maps = []
        self._offsets = [0]
        self._profile = []
        self._indices = {}
        self.raw_points = 0
        self.disk_bytes = 0
        self.stats = {}
        self.closed = False

    def __enter__(self):
        return self

    def _map(self, name, dtype):
        size = self._offsets[-1]
        dtype = np.dtype(dtype)
        if not size:
            return np.empty(0, dtype=dtype)
        path = self.directory / (name + ".bin")
        mode = "r" if path.exists() else "w+"
        array = np.memmap(path, mode=mode, dtype=dtype, shape=(size,))
        self._maps.append(array)
        self.disk_bytes += array.nbytes
        return array

    def zeros(self, name, dtype):
        # New extended file bytes are zero-filled by the filesystem; do not
        # dirty every mapped page just to initialize a full-file work array.
        if (self.directory / (name + ".bin")).exists():
            raise RuntimeError(f"Noise workspace already exists: {name}")
        array = self._map(name, dtype)
        return [array[a:b] for a, b in zip(self._offsets[:-1], self._offsets[1:])]

    def build(self, scans, records, config, configs):
        started = time.perf_counter()
        profile_status = Counter()
        empty = np.empty(0, dtype=np.float32)
        with (
            (self.directory / "mz.bin").open("wb") as mz_file,
            (self.directory / "intensity.bin").open("wb") as intensity_file,
        ):
            for scan in scans:
                self.raw_points += len(scan.mz)
                was_profile = not scan.centroided
                if was_profile:
                    if config.profile_mode != "centroid":
                        raise ValueError(
                            f"Profile MS1 spectrum {scan.native_id}; profile_mode is 'reject'"
                        )
                    mz, intensity, _, status = PROFILE.detect_profile(
                        scan, config, "profile_specific"
                    )
                    profile_status.update(status)
                else:
                    mz, intensity = scan.mz.astype(np.float32), scan.intensity.astype(np.float32)
                mz_file.write(np.ascontiguousarray(mz, dtype=np.float32).tobytes())
                intensity_file.write(np.ascontiguousarray(intensity, dtype=np.float32).tobytes())
                self._offsets.append(self._offsets[-1] + len(mz))
                self.frames.append(replace(scan, mz=empty, intensity=empty, centroided=True))
                self._profile.append(was_profile)
        mz = self._map("mz", np.float32)
        intensity = self._map("intensity", np.float32)
        for i, frame in enumerate(self.frames):
            a, b = self._offsets[i : i + 2]
            frame.mz, frame.intensity = mz[a:b], intensity[a:b]
            self._indices[id(frame)] = i
        prepared = time.perf_counter()
        parameters = MODEL.Parameters(bin_width_da=float(config.noise_window_da))
        effective = {f.polarity: config_for(f.polarity, config, configs) for f in self.frames}
        # Reference hooks exist only for independent benchmark processes. The
        # installed implementation supplies disk-backed workspace support.
        if MODEL.__name__ == "masscubex.consistency_noise":
            evidence = MODEL.precursor_protection(
                self.frames, records, effective, parameters, storage=self
            )
            self.decisions, _, statuses = MODEL.evaluate(
                self.frames,
                effective,
                parameters,
                ms2_evidence=evidence,
                storage=self,
                include_published=False,
            )
        else:
            evidence = MODEL.precursor_protection(self.frames, records, effective, parameters)
            self.decisions, _, statuses = MODEL.evaluate(
                self.frames,
                effective,
                parameters,
                ms2_evidence=evidence,
            )
        self.stats = {
            "parameters": asdict(parameters),
            "profile_status": dict(profile_status),
            "component_status_counts": statuses,
            "centroids_before_filter": self._offsets[-1],
            "temporary_disk_bytes": self.disk_bytes,
            "preparation_seconds": prepared - started,
            "estimator_seconds": time.perf_counter() - prepared,
        }
        return self

    def preprocess(self, scan, config):
        if self.closed:
            raise RuntimeError("Consistency cache is closed")
        index = self._indices[id(scan)]
        decision = self.decisions[index]
        keep = decision.keep()
        threshold = decision.threshold
        if config.min_intensity:
            keep &= scan.intensity > config.min_intensity
            threshold = np.maximum(threshold, config.min_intensity)
        mz, intensity, threshold = scan.mz[keep], scan.intensity[keep], threshold[keep]
        if self._profile[index]:
            return mz, intensity, threshold
        return merge_centroids(mz, intensity, threshold, config.centroid_mz_tolerance_da)

    def apply_stats(self, stats):
        stats["input_points"] = self.raw_points
        stats["profile_scans"] = sum(self._profile)
        stats["consistency_noise"] = self.stats

    def close(self):
        if self.closed:
            return
        self.decisions.clear()
        self.frames.clear()
        self._indices.clear()
        for array in self._maps:
            array._mmap.close()
        self._maps.clear()
        self._temporary.cleanup()
        self.closed = True

    def __exit__(self, *exc):
        self.close()
