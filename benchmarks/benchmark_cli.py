"""Measure the complete extraction command, including reporting and startup.

Use the same manifest/order as benchmark_extraction.py. Run this separately
from all other measurements and audits; the optional warmup file is excluded.
"""

import argparse
from collections import Counter
import json
from pathlib import Path
import statistics
import subprocess
import sys
import time

from masscubex import __version__
from masscubex.pipeline import implementation_id


def command(arguments):
    start = time.perf_counter()
    result = subprocess.run(
        [sys.executable, "-m", "masscubex", "extract", *arguments],
        capture_output=True,
        text=True,
    )
    elapsed = time.perf_counter() - start
    if result.returncode:
        raise RuntimeError(result.stderr)
    rows = [json.loads(line) for line in result.stdout.splitlines()]
    if any(row["status"] != "ok" for row in rows):
        raise RuntimeError("Expected fresh successful outputs for every benchmark input")
    return {"wall_seconds": elapsed, "files": rows}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--directory", type=Path, default=Path("validation/MSV000086989"))
    parser.add_argument("--output-directory", type=Path, required=True)
    parser.add_argument("--report", type=Path, required=True)
    parser.add_argument("--workers", type=int, default=4)
    parser.add_argument("--repeats", type=int, default=3)
    args = parser.parse_args()
    records = json.loads((args.directory / "download_manifest.json").read_text())
    sources = [(args.directory / row["name"]).resolve() for row in records]
    args.output_directory.mkdir(parents=True, exist_ok=True)
    manifest = args.output_directory / "benchmark-inputs.txt"
    manifest.write_text("".join(str(path) + "\n" for path in sources))
    shared = ["-o", str(args.output_directory), "--format", "parquet", "--overwrite"]
    command([str(sources[0]), *shared])
    rounds = []
    for repeat in range(args.repeats):
        row = command(["--manifest", str(manifest), "--workers", str(args.workers), *shared])
        if len(row["files"]) != len(sources):
            raise RuntimeError("Benchmark did not process the complete manifest")
        rounds.append(row)
        print(json.dumps({"repeat": repeat + 1, "seconds": row["wall_seconds"]}), flush=True)
    report = {
        "version": __version__,
        "implementation_sha256": implementation_id()[0],
        "workers": args.workers,
        "executor": "process",
        "format": "parquet",
        "input_files": len(sources),
        "median_wall_seconds": statistics.median(row["wall_seconds"] for row in rounds),
        "rounds": rounds,
        "notes": [
            "Full CLI command, including import/startup, bounded scheduling and batch reports.",
            "Same input manifest order and scientific defaults as the processing API benchmark.",
            "Warmup file excluded; fresh outputs on every pass; no reference audit during timing.",
        ],
    }
    stages = Counter()
    for row in rounds[-1]["files"]:
        metadata = json.loads(Path(row["output"] + ".json").read_text())
        stages.update(metadata["timings_seconds"])
    report["stage_seconds_last_repeat"] = dict(stages)
    args.report.parent.mkdir(parents=True, exist_ok=True)
    args.report.write_text(json.dumps(report, indent=2) + "\n")


if __name__ == "__main__":
    main()
