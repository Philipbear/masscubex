"""Inspect archived v0.9 features against independently decoded, unfiltered MS1.

This diagnostic does not run or change feature detection. Archived v0.8 tables
contain the unrounded coordinates of the identical v0.9 detections. Original
MassCube tables come from audit_upstream_features.py at its pinned commit.
"""

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
import pyarrow.parquet as pq
from numba import njit
from pyteomics import mzxml


COLUMNS = [
    "mz",
    "rt_apex_seconds",
    "rt_start_seconds",
    "rt_end_seconds",
    "peak_height",
    "scan_count",
]


def raw_ms1(source, cache):
    """Cache independent Pyteomics decoding; verify the source bytes on reuse."""
    digest = hashlib.sha256(source.read_bytes()).hexdigest()
    path = cache / (source.name + ".npz")
    if path.exists():
        data = np.load(path)
        if str(data["sha256"]) == digest:
            return tuple(data[key] for key in ("times", "masses", "intensities", "offsets"))
    times, masses, intensities, offsets = [], [], [], [0]
    with mzxml.MzXML(str(source), use_index=False) as reader:
        for spectrum in reader:
            if spectrum["msLevel"] != 1:
                continue
            rt = spectrum["retentionTime"]
            if str(rt.unit_info) != "minute":
                raise ValueError("Unexpected Pyteomics RT unit")
            mz = np.asarray(spectrum["m/z array"], dtype=np.float64)
            y = np.asarray(spectrum["intensity array"], dtype=np.float64)
            if np.any(np.diff(mz) < 0):
                raise ValueError("Raw centroids must be sorted")
            times.append(float(rt) * 60)
            masses.append(mz)
            intensities.append(y)
            offsets.append(offsets[-1] + len(mz))
    arrays = np.array(times), np.concatenate(masses), np.concatenate(intensities), np.array(offsets)
    cache.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(
        path, sha256=digest, **dict(zip(("times", "masses", "intensities", "offsets"), arrays))
    )
    return arrays


@njit(cache=True)
def extract_eics(masses, intensities, offsets, targets, radii):
    """Maximum centroid per scan, matching the measurement convention, without filtering."""
    result = np.zeros((len(targets), len(offsets) - 1))
    for i in range(len(offsets) - 1):
        start, stop = offsets[i : i + 2]
        mz = masses[start:stop]
        for j in range(len(targets)):
            left = start + np.searchsorted(mz, targets[j] - radii[j])
            right = start + np.searchsorted(mz, targets[j] + radii[j], side="right")
            for k in range(left, right):
                result[j, i] = max(result[j, i], intensities[k])
    return result


def tables(root, name):
    current = pq.read_table(
        root / "optimization-v0.8/final-parquet" / (name + ".parquet"),
        columns=COLUMNS + ["feature_id", "feature_origin"],
    ).to_pandas()
    original = pq.read_table(
        root / "output-v0.9/upstream-orbitrap" / (name + ".parquet"), columns=COLUMNS
    ).to_pandas()
    floor = pq.read_table(
        root / "output-v0.9/floor-50000" / (name + ".parquet"),
        columns=COLUMNS + ["feature_id", "feature_origin"],
    ).to_pandas()
    metadata = json.loads((root / "output-v0.9/parquet" / (name + ".parquet.json")).read_text())
    config = metadata["parameter_tuning"]["effective_configs"]["positive"]
    return (
        current[current.feature_origin == "ms1_detected"],
        original,
        floor[floor.feature_origin == "ms1_detected"],
        config,
    )


def screen(root, name, cache):
    current, original, floor, config = tables(root, name)
    t, masses, y, offsets = raw_ms1(root / "MSV000086989" / name, cache)
    width = current.rt_end_seconds - current.rt_start_seconds
    candidates = current[
        (current.peak_height.between(10000, 48000))
        & (current.scan_count >= 10)
        & width.between(2, 40)
        & current.rt_apex_seconds.between(20, 580)
    ]
    targets = candidates.mz.to_numpy()
    radii = np.maximum(config["mz_tolerance_da"], targets * config["mz_tolerance_ppm"] * 1e-6)
    traces = extract_eics(masses, y, offsets, targets, radii)
    weak = []
    for row, trace in zip(candidates.itertuples(), traces):
        start, end = row.rt_start_seconds, row.rt_end_seconds
        inside = (t >= start) & (t <= end)
        flank = ((t >= start - 15) & (t < start - 2)) | ((t > end + 2) & (t <= end + 15))
        if not flank.any() or not inside.any():
            continue
        height = float(trace[inside].max())
        baseline = float(np.quantile(trace[flank], 0.75))
        ratio = height / max(baseline, 1)
        original_count = int(
            (
                (abs(original.mz - row.mz) <= 0.01)
                & (original.rt_apex_seconds >= start - 2)
                & (original.rt_apex_seconds <= end + 2)
            ).sum()
        )
        if ratio >= 4 and height <= 50000 and original_count == 0:
            weak.append(
                dict(
                    feature_id=row.feature_id,
                    mz=row.mz,
                    rt=row.rt_apex_seconds,
                    start=start,
                    end=end,
                    height=height,
                    flank_q75=baseline,
                    peak_to_flank=ratio,
                    scans=row.scan_count,
                )
            )
    weak.sort(key=lambda r: (-r["scans"], -r["peak_to_flank"], r["feature_id"]))
    split = []
    for row in original[original.scan_count >= 20].itertuples():
        a, b = row.rt_start_seconds, row.rt_end_seconds
        if b - a < 20:
            continue
        match = floor[
            (abs(floor.mz - row.mz) < max(0.001, row.mz * 5e-6))
            & floor.rt_apex_seconds.between(a, b)
        ]
        if len(match) >= 4:
            split.append(
                dict(
                    original_row=int(row.Index + 1),
                    mz=row.mz,
                    rt=row.rt_apex_seconds,
                    start=a,
                    end=b,
                    original_scans=row.scan_count,
                    floor_x_count=len(match),
                    height=row.peak_height,
                )
            )
    split.sort(key=lambda r: (-r["floor_x_count"], r["original_row"]))
    return dict(
        file=name,
        current_rows=len(current),
        original_rows=len(original),
        original_at_least_5=int((original.scan_count >= 5).sum()),
        weak_candidates=len(weak),
        weak=weak[:20],
        split_candidates=len(split),
        split=split[:20],
    )


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--validation", type=Path, default=Path("validation"))
    parser.add_argument("--output", type=Path, default=Path("validation/eic-review"))
    parser.add_argument("--files", nargs="*")
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=True)
    manifest = json.loads((args.validation / "MSV000086989/download_manifest.json").read_text())
    names = args.files or [r["name"] for r in manifest if r["group"] in ("blank", "QC")]
    results = []
    for name in names:
        result = screen(args.validation, name, args.output / "raw")
        results.append(result)
        print(json.dumps(result), flush=True)
        (args.output / "screen.json").write_text(json.dumps(results, indent=2) + "\n")


if __name__ == "__main__":
    main()
