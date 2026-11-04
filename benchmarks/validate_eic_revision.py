"""Validate the v0.10 scientific change separately from unchanged MS1 tracing.

Run after validate_ms2_coverage.py. This checks default Parquet/TSV equality
on all 50 files, and replays the v0.9 noise setting on every blank/QC to prove
that its detected-feature measurements are unchanged by the gap-fill fix.
"""

from collections import Counter, defaultdict
from dataclasses import replace
import json
from itertools import zip_longest
from pathlib import Path

import pyarrow.parquet as pq

from masscubex import Config, __version__, iter_feature_rows, process_file
from masscubex.features import warmup
from masscubex.pipeline import implementation_id


ROOT = Path("validation")
OUT = ROOT / "eic-review/v0.10"
CORE = [
    "feature_id",
    "mz",
    "rt_apex_seconds",
    "rt_start_seconds",
    "rt_end_seconds",
    "peak_height",
    "peak_area",
    "scan_count",
    "apex_noise_threshold",
    "polarity",
    "ms2_count",
    "ms2_scan_numbers",
    "ms2_native_ids",
    "ms2_spectrum_indices",
    "best_ms2_scan_number",
    "best_ms2_native_id",
    "best_ms2_spectrum_index",
]


def main():
    manifest = json.loads((ROOT / "MSV000086989/download_manifest.json").read_text())
    coverage = json.loads(
        Path(f"benchmarks/results/MSV000086989-ms2-v{__version__}.json").read_text()
    )
    assert coverage["files"] == 50 and coverage["implementation_sha256"] == implementation_id()[0]
    OUT.mkdir(parents=True, exist_ok=True)
    warmup()
    # Warm read-only signatures and output paths before timing the serial sweep.
    process_file(
        ROOT / "MSV000086989" / manifest[0]["name"], OUT / "warmup.parquet", overwrite=True
    )
    rows = []
    for record in manifest:
        source = ROOT / "MSV000086989" / record["name"]
        output = OUT / "parquet" / (source.name + ".parquet")
        result = process_file(source, output, overwrite=True)
        metadata = json.loads(Path(str(output) + ".json").read_text())
        tsv = ROOT / "MSV000086989" / f"v{__version__}-ms2" / (source.name + ".tsv")
        tsv_metadata = json.loads(Path(str(tsv) + ".json").read_text())
        assert (
            metadata["implementation_sha256"]
            == tsv_metadata["implementation_sha256"]
            == implementation_id()[0]
        )
        counts = Counter()
        for a, b in zip_longest(iter_feature_rows(output), iter_feature_rows(tsv)):
            assert a == b, (source.name, a and a["feature_id"])
            counts[a["feature_origin"]] += 1
            counts[a["ms1_quant_status"]] += 1
        rows.append(
            dict(
                file=source.name,
                group=record["group"],
                input_sha256=record["sha256"],
                features=result.features,
                counts=dict(counts),
                elapsed_seconds=result.elapsed_seconds,
                timings_seconds=metadata["timings_seconds"],
                parquet_bytes=output.stat().st_size,
                tsv_bytes=tsv.stat().st_size,
            )
        )
        if len(rows) % 10 == 0:
            print("Parquet/TSV files verified:", len(rows), flush=True)
    unchanged = []
    old_noise = replace(Config(), noise_layout="sliding", noise_window_da=50)
    for record in manifest:
        if record["group"] not in ("blank", "QC"):
            continue
        source = ROOT / "MSV000086989" / record["name"]
        output = OUT / "old-noise-strict-gapfill" / (source.name + ".parquet")
        process_file(source, output, config=old_noise, overwrite=True)
        before = pq.read_table(
            ROOT / "output-v0.9/parquet" / (source.name + ".parquet"),
            columns=CORE + ["feature_origin"],
        ).to_pylist()
        after = pq.read_table(output, columns=CORE + ["feature_origin"]).to_pylist()
        before = [r for r in before if r["feature_origin"] == "ms1_detected"]
        after = [r for r in after if r["feature_origin"] == "ms1_detected"]
        assert before == after, source.name
        unchanged.append(dict(file=source.name, detected_features_unchanged=len(before)))
    group_counts = defaultdict(Counter)
    for row in rows:
        group_counts[row["group"]].update(row["counts"])
    report = dict(
        version=__version__,
        implementation_sha256=implementation_id()[0],
        config=Config().to_dict(),
        files=50,
        all_parquet_tsv_rows_equal=True,
        group_counts=dict(group_counts),
        serial_warm_processing_seconds=sum(r["elapsed_seconds"] for r in rows),
        timing_note="One warm serial sweep on local input; independent audit and comparisons excluded. Not a repeated cross-version benchmark.",
        unchanged_ms1_with_old_noise=unchanged,
        rows=rows,
    )
    Path("benchmarks/results/eic-revision-v0.10.0.json").write_text(
        json.dumps(report, indent=2) + "\n"
    )
    print(
        json.dumps(
            {
                "files": len(rows),
                "rows": sum(r["features"] for r in rows),
                "unchanged_blank_qc_detections": sum(
                    r["detected_features_unchanged"] for r in unchanged
                ),
            }
        )
    )


if __name__ == "__main__":
    main()
