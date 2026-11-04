"""Explain missing calibration evidence without relaxing production defaults.

Compare global and region-stratified pilot sampling on the saved centroid-mzXML
cohort. Use the same tracing/quality rules, and retain pre-Gaussian candidates
so shape thresholds and evidence minima can be compared without rerunning data.
This is an exploratory calibration audit, not a feature-table benchmark.
"""

from collections import Counter
from dataclasses import asdict, replace
import argparse
import json
from pathlib import Path
import time

import numpy as np

from masscubex import Config
from masscubex.features import TraceTracker, MEAN_MZ, OBSERVATION_COUNT
from masscubex.pipeline import implementation_id
from masscubex.readers import iter_ms1
from masscubex.tuning import representative, competing_peaks, _screen_representatives


BREAKS = np.array([200.0, 400.0, 600.0])
MODES = {
    "global128": (False, 128),
    "regional32": (True, 32),
    "global512": (False, 512),
    "regional128": (True, 128),
}


def strongest(values, maximum):
    if len(values) <= maximum:
        return np.arange(len(values))
    cutoff = np.partition(values, len(values) - maximum)[len(values) - maximum]
    higher = np.flatnonzero(values > cutoff)
    tied = np.flatnonzero(values == cutoff)[: maximum - len(higher)]
    return np.sort(np.r_[higher, tied])


def indices_for(masses, values, regional, maximum):
    if not regional:
        return strongest(values, maximum)
    cuts = [0, *np.searchsorted(masses, BREAKS, side="left"), len(masses)]
    return np.concatenate(
        [start + strongest(values[start:stop], maximum) for start, stop in zip(cuts, cuts[1:])]
    )


class AuditedTracker(TraceTracker):
    def __init__(self, config):
        super().__init__(config)
        self.rejections = [Counter() for _ in range(4)]

    def _finish(self, finished_tracks):
        short = finished_tracks[self.track_layout[finished_tracks, OBSERVATION_COUNT] < 7]
        counts = np.bincount(
            np.searchsorted(BREAKS, self.track_values[short, MEAN_MZ], side="right"), minlength=4
        )
        for index, count in enumerate(counts):
            self.rejections[index]["fewer_than_7_observations"] += int(count)
        return super()._finish(finished_tracks)


def summarize(candidates, config):
    result = {}
    for shape in (0.0, 0.8, 0.9):
        qualified = [p for p in candidates if p.gaussian_similarity >= shape]
        selected, _ = _screen_representatives(qualified, config)
        bins = []
        for index in range(4):
            region = [p for p in qualified if np.searchsorted(BREAKS, p.mz, side="right") == index]
            accepted, before = _screen_representatives(region, config)
            bins.append(
                {
                    "candidate_count": len(region),
                    "before_outliers": before,
                    "accepted": len(accepted),
                }
            )
        result[str(shape)] = {"global_accepted": len(selected), "regions": bins}
    return result


def audit_file(path):
    config = Config(tuning_min_gaussian_similarity=0)
    pilot = replace(
        config,
        min_scans=7,
        max_gap_scans=1,
        group_ions=False,
        mz_tolerance_da=config.tuning_search_da,
        mz_tolerance_ppm=config.tuning_search_ppm,
        mass_regions=(),
    )
    trackers = {name: AuditedTracker(pilot) for name in MODES}
    candidates = {name: [] for name in MODES}
    selected_counts = {name: np.zeros(4, dtype=np.int64) for name in MODES}
    raw_counts = np.zeros(4, dtype=np.int64)
    raw_intensity = np.zeros(4)
    start, block_scans, total_scans = None, 0, 0

    def consume(name, traces):
        for trace in traces:
            peak, reason = representative(trace, config)
            if peak is None:
                mass = np.median(trace[trace[:, 1] > 0, 0])
                index = int(np.searchsorted(BREAKS, mass, side="right"))
                trackers[name].rejections[index][reason] += 1
            else:
                candidates[name].append(peak)

    started = time.perf_counter()
    for scan in iter_ms1(path):
        assert scan.centroided and scan.polarity == "positive", "Audit is specific to this cohort"
        if start is None:
            start = scan.rt
        if (
            scan.rt - start >= config.tuning_block_seconds
            or block_scans >= config.tuning_block_scans
        ):
            for name, tracker in trackers.items():
                consume(name, tracker.flush())
            start, block_scans = scan.rt, 0
        block_scans += 1
        total_scans += 1
        bins = np.searchsorted(BREAKS, scan.mz, side="right")
        raw_counts += np.bincount(bins, minlength=4)
        raw_intensity += np.bincount(bins, weights=scan.intensity, minlength=4)
        for name, (regional, maximum) in MODES.items():
            selected = indices_for(scan.mz, scan.intensity, regional, maximum)
            selected_counts[name] += np.bincount(bins[selected], minlength=4)
            flags = competing_peaks(
                scan.mz, scan.intensity, selected, config.tuning_search_ppm, config.tuning_search_da
            )
            consume(
                name,
                trackers[name].add(scan.mz[selected], scan.intensity[selected], flags, scan.rt),
            )
    for name, tracker in trackers.items():
        consume(name, tracker.flush())
    return {
        "file": path.name,
        "scans": total_scans,
        "raw_centroids": raw_counts.tolist(),
        "raw_intensity": raw_intensity.tolist(),
        "audit_seconds": time.perf_counter() - started,
        "modes": {
            name: {
                "selected_centroids": selected_counts[name].tolist(),
                "rejections": [dict(c) for c in trackers[name].rejections],
                "summary": summarize(candidates[name], config),
                "candidates": [asdict(p) for p in candidates[name]],
            }
            for name in MODES
        },
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--limit", type=int)
    parser.add_argument("--workers", type=int, default=4)
    parser.add_argument(
        "--output",
        type=Path,
        default=Path("validation/calibration-coverage-v0.13/pilot-audit.json"),
    )
    args = parser.parse_args()
    root = Path("validation/MSV000086989")
    manifest = json.loads((root / "download_manifest.json").read_text())[: args.limit]
    from concurrent.futures import ProcessPoolExecutor, as_completed

    rows = []
    with ProcessPoolExecutor(max_workers=args.workers) as pool:
        futures = [pool.submit(audit_file, root / r["name"]) for r in manifest]
        for future in as_completed(futures):
            rows.append(future.result())
            if len(rows) % 10 == 0 or len(rows) == len(manifest):
                print("Audited files:", len(rows), flush=True)
    report = {
        "implementation_sha256": implementation_id()[0],
        "region_breaks": BREAKS.tolist(),
        "modes": MODES,
        "rows": sorted(rows, key=lambda r: r["file"]),
    }
    args.output.parent.mkdir(exist_ok=True, parents=True)
    args.output.write_text(json.dumps(report, indent=2, allow_nan=False) + "\n")
    for name in MODES:
        for shape in ("0.8", "0.9"):
            print(
                name,
                shape,
                {
                    n: [
                        sum(
                            r["modes"][name]["summary"][shape]["global_accepted"] >= 20
                            and r["modes"][name]["summary"][shape]["regions"][i]["accepted"] >= n
                            for r in rows
                        )
                        for i in range(4)
                    ]
                    for n in (10, 20)
                },
                flush=True,
            )


if __name__ == "__main__":
    main()
