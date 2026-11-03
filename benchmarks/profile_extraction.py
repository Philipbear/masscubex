"""Time coarse extraction stages on one file from each validation sample group.

Timers are inclusive within a stage. Extraction substeps are a breakdown of
feature_extraction, not additional elapsed time. Run separately from benchmarks.
"""

import argparse
from collections import defaultdict
from contextlib import contextmanager
import json
from pathlib import Path
import time

from masscubex import __version__, pipeline, gapfill, grouping
from masscubex.features import warmup
from masscubex.pipeline import implementation_id


@contextmanager
def instrument():
    timers = defaultdict(float)
    originals = []

    def patch(module, name, label, generator=False):
        original = getattr(module, name)
        originals.append((module, name, original))

        def call(*args, **kwargs):
            started = time.perf_counter()
            try:
                return original(*args, **kwargs)
            finally:
                timers[label] += time.perf_counter() - started

        def iterate(*args, **kwargs):
            iterator = iter(original(*args, **kwargs))
            try:
                while True:
                    started = time.perf_counter()
                    try:
                        item = next(iterator)
                    except StopIteration:
                        timers[label] += time.perf_counter() - started
                        return
                    timers[label] += time.perf_counter() - started
                    yield item
            finally:
                if hasattr(iterator, "close"):
                    iterator.close()

        setattr(module, name, iterate if generator else call)

    patch(pipeline, "calibrate_file", "calibration")
    patch(pipeline, "extract_features", "feature_extraction", True)
    patch(pipeline, "associate_ms2", "ms2_association")
    patch(gapfill, "seed_ms2_features", "ms2_gapfill")
    patch(grouping, "group_ions", "ion_grouping")
    if hasattr(pipeline, "iter_scans"):
        patch(pipeline, "iter_scans", "extraction.reading", True)
    else:
        # The same bounded source supplies calibration, extraction and gapfill.
        patch(pipeline.ScanSource, "scans", "shared.scan_reading_and_replay", True)
    patch(pipeline, "estimate_noise", "extraction.noise")
    patch(pipeline, "summarize_trace", "extraction.peak_summarization", True)
    # This class is shared with calibration, so its timer overlaps both stages.
    patch(pipeline.TraceTracker, "add", "shared.trace_tracking")
    try:
        yield timers
    finally:
        for module, name, original in reversed(originals):
            setattr(module, name, original)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--directory", type=Path, default=Path("validation/MSV000086989"))
    parser.add_argument("--output-directory", type=Path, required=True)
    parser.add_argument("--report", type=Path, required=True)
    parser.add_argument("--format", default="tsv", choices=["tsv", "parquet"])
    parser.add_argument("--repeats", type=int, default=3)
    args = parser.parse_args()
    manifest = json.loads((args.directory / "download_manifest.json").read_text())
    selected = {}
    for row in manifest:
        selected.setdefault(row["group"], row)
    args.output_directory.mkdir(parents=True, exist_ok=True)
    warmup()
    rows = []
    for group, row in selected.items():
        source = args.directory / row["name"]
        out = args.output_directory / (row["name"] + "." + args.format)
        pipeline.process_file(source, out, overwrite=True)
        for repeat in range(args.repeats):
            with instrument() as timers:
                result = pipeline.process_file(source, out, overwrite=True)
            stages = {
                k: v for k, v in timers.items() if not k.startswith(("extraction.", "shared."))
            }
            stages["output_and_provenance"] = result.elapsed_seconds - sum(stages.values())
            rows.append(
                {
                    "name": row["name"],
                    "group": group,
                    "repeat": repeat + 1,
                    "seconds": result.elapsed_seconds,
                    "stages": stages,
                    "extraction_substeps": {
                        k.split(".")[1]: v for k, v in timers.items() if k.startswith("extraction.")
                    },
                    "shared_substeps": {
                        k.split(".")[1]: v for k, v in timers.items() if k.startswith("shared.")
                    },
                }
            )
    report = {
        "version": __version__,
        "implementation_sha256": implementation_id()[0],
        "dependencies": implementation_id()[1],
        "format": args.format,
        "rows": rows,
        "notes": [
            "Lightweight timers, not cProfile; one file per sample group.",
            "Extraction includes reader, noise, tracking and peak summarization.",
            "Substeps are inclusive and must not be added to top-level stages.",
            "Shared trace tracking includes calibration and extraction; it is not exclusive to either stage.",
        ],
    }
    args.report.write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps({"version": __version__, "profiled_runs": len(rows)}))


if __name__ == "__main__":
    main()
