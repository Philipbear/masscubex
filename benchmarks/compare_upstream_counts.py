"""Compare real-file feature counts and isolate the absolute MS1 floor effect."""

import argparse
from collections import Counter, defaultdict
from dataclasses import replace
import hashlib
import json
from pathlib import Path
import subprocess

import pyarrow.parquet as pq

from masscubex import Config, __version__, process_file
from masscubex.pipeline import implementation_id


UPSTREAM_COMMIT = "356748ede691862b9a87f67812db5e728a703b9d"


def table_counts(path, origin_column=True):
    columns = ["peak_height", "scan_count", "ms2_count"]
    if origin_column:
        columns += ["feature_origin", "ms1_quant_status"]
    data = pq.read_table(path, columns=columns).to_pydict()
    counts = Counter()
    for i, height in enumerate(data["peak_height"]):
        count = data["scan_count"][i]
        ms2 = data["ms2_count"][i]
        origin = data["feature_origin"][i] if origin_column else "ms1_detected"
        counts["features"] += 1
        counts[origin] += 1
        counts["at_least_5_scans"] += count >= 5
        counts["upstream_alignment_count_rule"] += count > 5 or ms2 > 0
        counts["with_ms2"] += ms2 > 0
        if height is not None:
            counts["height_below_50000"] += height < 50000
            counts["height_at_least_50000"] += height >= 50000
            counts["height_at_least_1000000"] += height >= 1e6
        if origin_column:
            counts["ambiguous_mass"] += data["ms1_quant_status"][i] == "ambiguous_mass"
    return dict(counts)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, required=True)
    parser.add_argument("--upstream", type=Path, required=True)
    parser.add_argument("--upstream-source", type=Path, required=True)
    parser.add_argument("--candidate", type=Path, required=True)
    parser.add_argument("--floor-output", type=Path, required=True)
    parser.add_argument("--unsegmented-output", type=Path, required=True)
    parser.add_argument("--report", type=Path, required=True)
    args = parser.parse_args()
    # Verify that every Python source in the audit checkout is the pinned
    # original, not a modified local algorithm bearing the same name.
    source_hashes = {}
    for path in sorted((args.upstream_source / "src" / "masscube").rglob("*.py")):
        relative = path.relative_to(args.upstream_source).as_posix()
        expected = subprocess.check_output(["git", "show", f"{UPSTREAM_COMMIT}:{relative}"])
        assert path.read_bytes() == expected, relative
        source_hashes[relative] = hashlib.sha256(expected).hexdigest()
    manifest = json.loads((args.input / "download_manifest.json").read_text())
    config = replace(Config(), noise_mode="none", min_intensity=50000)
    results = []
    totals = defaultdict(Counter)
    groups = defaultdict(lambda: defaultdict(Counter))
    for entry in manifest:
        name = entry["name"]
        source = args.input / name
        assert hashlib.sha256(source.read_bytes()).hexdigest() == entry["sha256"]
        original = args.upstream / (name + ".parquet")
        candidate = args.candidate / (name + ".parquet")
        floor = args.floor_output / (name + ".parquet")
        unsegmented = args.unsegmented_output / (name + ".parquet")
        process_file(source, floor, config=config, resume=True)
        process_file(source, unsegmented, config=replace(config, segment=False), resume=True)
        old_meta = json.loads(Path(str(original) + ".json").read_text())
        new_meta = json.loads(Path(str(candidate) + ".json").read_text())
        floor_meta = json.loads(Path(str(floor) + ".json").read_text())
        assert old_meta["input_sha256"] == entry["sha256"]
        # Calibration deliberately uses raw representative peaks, so changing
        # the detection floor should not change any learned tolerance/window.
        assert (
            floor_meta["parameter_tuning"]["polarities"]
            == new_meta["parameter_tuning"]["polarities"]
        )
        row = {
            "file": name,
            "group": entry["group"],
            "input_sha256": entry["sha256"],
            "original": table_counts(original, False),
            "masscubex": table_counts(candidate),
            "masscubex_floor_50000": table_counts(floor),
            "masscubex_floor_50000_unsplit": table_counts(unsegmented),
        }
        row["original"]["unsegmented_traces_at_least_5_scans"] = old_meta[
            "unsegmented_traces_at_least_5_scans"
        ]
        row["original"]["retained_ms1_points"] = old_meta["retained_ms1_points"]
        row["masscubex"]["retained_ms1_points"] = new_meta["retained_centroids"]
        row["masscubex_floor_50000"]["retained_ms1_points"] = floor_meta["retained_centroids"]
        for method in (
            "original",
            "masscubex",
            "masscubex_floor_50000",
            "masscubex_floor_50000_unsplit",
        ):
            totals[method].update(row[method])
            groups[entry["group"]][method].update(row[method])
        results.append(row)
        print(json.dumps(row), flush=True)
    report = {
        "dataset": "MSV000086989",
        "files": len(results),
        "version": __version__,
        "implementation_sha256": implementation_id()[0],
        "upstream_commit": UPSTREAM_COMMIT,
        "upstream_version": "1.2.22",
        "upstream_source_hashes": source_hashes,
        "upstream_parameters": {
            "instrument": "orbitrap",
            "ms1_abs_int_tol": 50000,
            "ms2_abs_int_tol": 10000,
            "centroid_mz_tol": 0.002,
            "mz_tol_ms1": 0.01,
            "feature_gap_tol": 10,
            "segmentation_passes": 2,
        },
        "masscubex_parameters": Config().to_dict(),
        "floor_parameters": config.to_dict(),
        "totals": totals,
        "by_group": groups,
        "rows": results,
        "limitations": [
            "Original pinned MassCube does not accept mzXML. Pyteomics records were translated to its mzML scan-extractor metadata schema; arrays and minute RTs are passed unchanged. Original preprocessing, detection, segmentation and finalization were not modified.",
            "Original single-file output retains short traces. At least five nonzero scans is a common diagnostic, not its single-file filter. Its downstream alignment count rule is >5 scans OR associated MS2; other alignment filters also exist.",
            "MassCubeX changes matching, segmentation and calibration as well as noise filtering. The fixed-floor comparison isolates the detection-floor change within MassCubeX, not every inter-software difference.",
            "Counts on this positive-mode centroid Orbitrap cohort cannot establish false-positive rates or generalize to every repository/instrument.",
        ],
    }
    args.report.write_text(json.dumps(report, indent=2, allow_nan=False) + "\n")
    print(json.dumps(totals))


if __name__ == "__main__":
    main()
