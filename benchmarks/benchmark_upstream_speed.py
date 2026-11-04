"""Five-file, identical-workload benchmark against pinned MassCube functions.

Both engines receive the same adaptive MS1 preprocessing, fixed mass windows,
MS2 thresholds and ion catalog. Extensions absent upstream are disabled only in
this comparison. Normal defaults are measured separately against saved tables.
The upstream timer excludes serialization; MassCubeX includes its full Parquet
transaction, making the speed ratio conservative for MassCubeX.
"""

import argparse
import ast
import gc
import hashlib
import json
from pathlib import Path
import platform
import pickle
import subprocess
import time

import verify_upstream_defaults as verify
from masscubex import Config, iter_feature_rows, process_file
from masscubex.features import warmup
from masscubex.pipeline import implementation_id


COMMIT = "356748ede691862b9a87f67812db5e728a703b9d"


def verify_sources():
    """Check the actual benchmark functions against the pinned Git objects."""
    snapshot = Path(verify.reference.__file__).read_text()
    nodes = {
        n.name: n
        for n in ast.parse(snapshot).body
        if isinstance(n, (ast.FunctionDef, ast.ClassDef))
    }
    files = {
        "src/masscube/feature_detection.py": [
            "Feature",
            "detect_features",
            "segment_feature",
            "gaussian_filter_and_prominence_method",
        ],
        "src/masscube/feature_grouping.py": [
            "group_features_single_file",
            "enumerate_source_ions",
            "scan_to_scan_cor_intensity",
        ],
    }
    hashes = {}
    for path, names in files.items():
        source = subprocess.check_output(["git", "show", f"{COMMIT}:{path}"]).decode()
        original = {
            n.name: n
            for n in ast.parse(source).body
            if isinstance(n, (ast.FunctionDef, ast.ClassDef))
        }
        for name in names:
            text = ast.get_source_segment(source, original[name])
            assert ast.get_source_segment(snapshot, nodes[name]) == text, name
            hashes[name] = hashlib.sha256(text.encode()).hexdigest()
    return hashes


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("inputs", type=Path, nargs=5)
    parser.add_argument("--output-directory", type=Path, required=True)
    parser.add_argument("--report", type=Path, required=True)
    parser.add_argument("--mode", choices=["matched", "orbitrap"], default="matched")
    args = parser.parse_args()
    config = Config(
        auto_tune=False,
        profile_mode="reject",
        ms2_seeded_features=False,
        ms2_min_intensity=500 if args.mode == "matched" else 10000,
        noise_mode="adaptive" if args.mode == "matched" else "none",
        min_intensity=0 if args.mode == "matched" else 50000,
    )
    args.output_directory.mkdir(parents=True, exist_ok=True)
    args.report.parent.mkdir(parents=True, exist_ok=True)
    report = {
        "upstream_commit": COMMIT,
        "upstream_version": "1.2.22",
        "source_hashes": verify_sources(),
        "implementation_sha256": implementation_id()[0],
        "dependencies": implementation_id()[1],
        "config": config.to_dict(),
        "python": platform.python_version(),
        "platform": platform.platform(),
        "mode": args.mode,
        "workers": 1,
        "rows": [],
        "notes": [
            "Same five input files and identical scientific parameters for both engines.",
            "Matched mode feeds both engines the same adaptive MS1 threshold; stock upstream does not implement it.",
            "Auto-tuning, profile centroiding and MS2 seeding disabled because upstream lacks these extensions.",
            "Selected adduct catalog and grouping enabled for both; upstream defaults normally disable single-file grouping.",
            "Pyteomics metadata adapter supplies mzXML arrays; pinned upstream only accepts mzML.",
            "Original timer includes decode/preprocessing, detection, splitting, finalization, MS2 and grouping; excludes serialization.",
            "MassCubeX timer includes the full process_file Parquet transaction, including serialization and provenance.",
            "Output parity verification and JIT/full-path warmup are outside timers. Single warm run per engine/file.",
        ],
    }
    warmup()
    process_file(
        args.inputs[0], args.output_directory / "warmup.parquet", config=config, overwrite=True
    )
    for source in args.inputs:
        gc.collect()
        print(json.dumps({"file": source.name, "stage": "original_masscube"}), flush=True)
        start = time.perf_counter()
        expected = verify.original_run(source, config, grouping=True)
        upstream_seconds = time.perf_counter() - start
        print(json.dumps({"file": source.name, "upstream_seconds": upstream_seconds}), flush=True)
        # Do not pin upstream's large object graph during the candidate timer.
        # This trusted, self-generated snapshot is used only for later parity.
        saved = args.output_directory / (source.name + ".reference.pkl")
        with saved.open("wb") as stream:
            pickle.dump(expected, stream, protocol=5)
        del expected
        gc.collect()
        output = args.output_directory / (source.name + ".parquet")
        result = process_file(source, output, config=config, overwrite=True)
        actual = list(iter_feature_rows(output))
        with saved.open("rb") as stream:
            expected = pickle.load(stream)
        saved.unlink()
        errors = verify.compare(expected, actual, True)
        row = {
            "file": source.name,
            "input_sha256": hashlib.sha256(source.read_bytes()).hexdigest(),
            "original_features": len(expected),
            "masscubex_features": len(actual),
            "upstream_processing_seconds": upstream_seconds,
            "masscubex_end_to_end_seconds": result.elapsed_seconds,
            "conservative_speedup": upstream_seconds / result.elapsed_seconds,
            "mismatched_rows": len(errors),
            "first_mismatches": errors[:3],
        }
        report["rows"].append(row)
        report["aggregate_speedup"] = sum(
            r["upstream_processing_seconds"] for r in report["rows"]
        ) / sum(r["masscubex_end_to_end_seconds"] for r in report["rows"])
        args.report.write_text(json.dumps(report, indent=2) + "\n")
        print(json.dumps(row), flush=True)
        if errors:
            raise AssertionError(f"Upstream parity failed: {source.name}")
        del expected, actual
        gc.collect()


if __name__ == "__main__":
    main()
