"""Reproducible warm-cache extraction benchmark; run each implementation separately.

Use PYTHONPATH to select an archived baseline. Raw inputs and outputs stay under
validation/. Timings exclude downloads, independent audits and profiler overhead.
"""

import argparse
from concurrent.futures import ProcessPoolExecutor
from dataclasses import asdict
from collections import Counter
import json
import multiprocessing
from pathlib import Path
import platform
import statistics
import time

from masscubex import Config, __version__, process_file
from masscubex.features import warmup
from masscubex.pipeline import implementation_id


def job(arguments):
    source, output = arguments
    return asdict(process_file(source, output, overwrite=True))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--directory", type=Path, default=Path("validation/MSV000086989"))
    parser.add_argument("--output-directory", type=Path, required=True)
    parser.add_argument("--report", type=Path, required=True)
    parser.add_argument("--format", choices=["tsv", "parquet"], default="tsv")
    parser.add_argument("--repeats", type=int, default=3)
    parser.add_argument("--workers", type=int, default=1)
    parser.add_argument("--limit", type=int)
    args = parser.parse_args()
    rows = json.loads((args.directory / "download_manifest.json").read_text())[: args.limit]
    args.output_directory.mkdir(parents=True, exist_ok=True)
    tasks = [
        (args.directory / row["name"], args.output_directory / (row["name"] + "." + args.format))
        for row in rows
    ]
    warmup()
    # Warm the whole processing path, including optional kernels and libraries.
    job(tasks[0])
    rounds = []
    for repeat in range(args.repeats):
        start = time.perf_counter()
        if args.workers == 1:
            outcomes = [job(task) for task in tasks]
        else:
            with ProcessPoolExecutor(
                max_workers=args.workers,
                initializer=warmup,
                mp_context=multiprocessing.get_context("spawn"),
            ) as pool:
                outcomes = list(pool.map(job, tasks))
        elapsed = time.perf_counter() - start
        rounds.append({"wall_seconds": elapsed, "files": outcomes})
        print(
            json.dumps({"repeat": repeat + 1, "files": len(tasks), "seconds": elapsed}), flush=True
        )
    report = {
        "version": __version__,
        "implementation_sha256": implementation_id()[0],
        "dependencies": implementation_id()[1],
        "python": platform.python_version(),
        "platform": platform.platform(),
        "config": Config().to_dict(),
        "format": args.format,
        "workers": args.workers,
        "repeats": args.repeats,
        "input_files": len(tasks),
        "input_bytes": sum(source.stat().st_size for source, _ in tasks),
        "output_bytes": sum(output.stat().st_size for _, output in tasks),
        "sidecar_bytes": sum(Path(str(output) + ".json").stat().st_size for _, output in tasks),
        "median_wall_seconds": statistics.median(r["wall_seconds"] for r in rounds),
        "rounds": rounds,
        "limitations": [
            "One positive-mode centroid mzXML cohort; warm filesystem cache and JIT.",
            "Includes extraction and output publication; excludes reference audits.",
            "Multiple workers include worker startup; no concurrent benchmarks.",
        ],
    }
    stages = Counter()
    for _, output in tasks:
        metadata = json.loads(Path(str(output) + ".json").read_text())
        stages.update(metadata.get("timings_seconds", {}))
    report["stage_seconds_last_repeat"] = dict(stages)
    if args.workers == 1:
        try:
            import resource

            rss = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
            report["process_peak_rss_mib"] = rss / (
                1024**2 if platform.system() == "Darwin" else 1024
            )
        except ImportError:
            pass
    args.report.parent.mkdir(parents=True, exist_ok=True)
    args.report.write_text(json.dumps(report, indent=2) + "\n")


if __name__ == "__main__":
    main()
