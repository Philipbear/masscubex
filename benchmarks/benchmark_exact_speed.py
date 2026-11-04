"""Compare extraction implementations on specified files, without algorithm changes.

Pass --package-path to select a saved source package. Reuse the same inputs for
both runs; --compare points to the baseline output directory. Warmup/profiling
is excluded from measured repeats. Parquet bytes must match, including every
measurement, label, relation and MS/MS reference.
"""

import argparse
import cProfile
import hashlib
import json
from pathlib import Path
import platform
import pstats
import statistics
import sys


def sha256(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def scientific_metadata(value):
    """Exclude only performance measurements and the changed source fingerprint."""
    if isinstance(value, dict):
        return {
            k: scientific_metadata(v)
            for k, v in value.items()
            if k not in {"elapsed_seconds", "timings_seconds", "implementation_sha256"}
        }
    if isinstance(value, list):
        return [scientific_metadata(v) for v in value]
    return value


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("inputs", type=Path, nargs="+")
    parser.add_argument("--package-path", type=Path)
    parser.add_argument("--output-directory", type=Path, required=True)
    parser.add_argument("--report", type=Path, required=True)
    parser.add_argument("--compare", type=Path)
    parser.add_argument("--config", type=Path, help="Optional shared scientific configuration JSON")
    parser.add_argument("--repeats", type=int, default=2)
    parser.add_argument("--profile-first", action="store_true")
    args = parser.parse_args()
    if args.repeats < 1:
        parser.error("--repeats must be positive")
    if args.package_path:
        sys.path.insert(0, str(args.package_path.resolve()))

    from masscubex import Config, process_file
    from masscubex.features import warmup
    from masscubex.pipeline import implementation_id

    config = Config.from_dict(json.loads(args.config.read_text())) if args.config else Config()

    args.output_directory.mkdir(parents=True, exist_ok=True)
    args.report.parent.mkdir(parents=True, exist_ok=True)
    tasks = [(p.resolve(), args.output_directory / (p.name + ".parquet")) for p in args.inputs]
    report = {
        "implementation_sha256": implementation_id()[0],
        "dependencies": implementation_id()[1],
        "python": platform.python_version(),
        "platform": platform.platform(),
        "config": config.to_dict(),
        "format": "parquet",
        "workers": 1,
        "repeats": args.repeats,
        "inputs": [{"name": p.name, "sha256": sha256(p)} for p, _ in tasks],
        "runs": [],
        "notes": [
            f"Exactly {len(tasks)} specified inputs; "
            + ("explicit shared configuration." if args.config else "default settings."),
            "One worker, warm filesystem/JIT; no profiler in measured repeats.",
            "Table bytes compared exactly; timings and implementation provenance necessarily differ.",
        ],
    }
    warmup()
    if args.profile_first:
        profiler = cProfile.Profile()
        profiler.runcall(process_file, *tasks[0], config=config, overwrite=True)
        profile_path = args.report.with_suffix(".prof")
        profiler.dump_stats(profile_path)
        with args.report.with_suffix(".profile.txt").open("w") as stream:
            pstats.Stats(profiler, stream=stream).strip_dirs().sort_stats("cumtime").print_stats(35)
    else:
        process_file(*tasks[0], config=config, overwrite=True)
    print("Warmup complete; starting measured repeats", flush=True)
    for repeat in range(args.repeats):
        for source, output in tasks:
            result = process_file(source, output, config=config, overwrite=True)
            metadata = json.loads(Path(str(output) + ".json").read_text())
            row = {
                "name": source.name,
                "repeat": repeat + 1,
                "seconds": result.elapsed_seconds,
                "stages": metadata["timings_seconds"],
                "features": result.features,
                "ms2_assigned": metadata["ms2_assigned"],
                "output_sha256": metadata["output_sha256"],
                "output_bytes": metadata["output_bytes"],
            }
            if args.compare:
                baseline = args.compare / output.name
                row["baseline_sha256"] = sha256(baseline)
                row["exact_table_match"] = row["output_sha256"] == row["baseline_sha256"]
                if not row["exact_table_match"]:
                    raise AssertionError(f"Feature table changed: {source.name}")
                old_metadata = json.loads(Path(str(baseline) + ".json").read_text())
                row["scientific_metadata_match"] = scientific_metadata(
                    old_metadata
                ) == scientific_metadata(metadata)
                if not row["scientific_metadata_match"]:
                    raise AssertionError(f"Scientific metadata changed: {source.name}")
            previous = [r for r in report["runs"] if r["name"] == source.name]
            if previous and previous[0]["output_sha256"] != row["output_sha256"]:
                raise AssertionError(f"Nondeterministic output: {source.name}")
            report["runs"].append(row)
            args.report.write_text(json.dumps(report, indent=2) + "\n")
            print(json.dumps(row), flush=True)
    report["medians"] = []
    for source, _ in tasks:
        runs = [r for r in report["runs"] if r["name"] == source.name]
        report["medians"].append(
            {
                "name": source.name,
                "seconds": statistics.median(r["seconds"] for r in runs),
                "stages": {
                    k: statistics.median(r["stages"][k] for r in runs) for k in runs[0]["stages"]
                },
            }
        )
    args.report.write_text(json.dumps(report, indent=2) + "\n")


if __name__ == "__main__":
    main()
