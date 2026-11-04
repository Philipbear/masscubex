"""Audit the shape-guard revision against archived v0.10 tables and raw-MS2 checks.

Run validate_ms2_coverage.py first. Time file processing serially, with no other
benchmark/test processes running. Input data and large outputs stay in validation/.
"""

from collections import Counter, defaultdict
from dataclasses import replace
from itertools import zip_longest
import json
from pathlib import Path
import statistics

import pyarrow.parquet as pq

from masscubex import Config, __version__, iter_feature_rows, process_file
from masscubex.features import warmup
from masscubex.pipeline import implementation_id

ROOT = Path("validation")
OUT = ROOT / "peak-shape-v0.11"


def run(record, guarded):
    source = ROOT / "MSV000086989" / record["name"]
    path = OUT / ("parquet" if guarded else "guard-off") / (source.name + ".parquet")
    result = process_file(
        source, path, config=replace(Config(), peak_shape_guard=guarded), overwrite=True
    )
    metadata = json.loads(Path(str(path) + ".json").read_text())
    return result, path, metadata


def main():
    manifest = json.loads((ROOT / "MSV000086989/download_manifest.json").read_text())
    coverage = json.loads(
        Path(f"benchmarks/results/MSV000086989-ms2-v{__version__}.json").read_text()
    )
    assert coverage["files"] == 50
    assert coverage["implementation_sha256"] == implementation_id()[0]
    warmup()
    run(manifest[0], True)
    rows = []
    groups = defaultdict(Counter)
    for record in manifest:
        result, path, metadata = run(record, True)
        tsv = ROOT / "MSV000086989" / f"v{__version__}-ms2" / (record["name"] + ".tsv")
        tsv_meta = json.loads(Path(str(tsv) + ".json").read_text())
        assert (
            metadata["implementation_sha256"]
            == tsv_meta["implementation_sha256"]
            == implementation_id()[0]
        )
        counts = Counter()
        for a, b in zip_longest(iter_feature_rows(path), iter_feature_rows(tsv)):
            assert a == b, (record["name"], a and a["feature_id"])
            counts[a["feature_origin"]] += 1
            counts[a["segmentation_status"]] += 1
            if a["feature_origin"] == "ms2_seeded":
                assert a["gaussian_similarity"] is None
            else:
                assert 0 <= a["gaussian_similarity"] <= 1
        before_path = ROOT / "eic-review/v0.10/parquet" / (record["name"] + ".parquet")
        before = pq.read_table(before_path, columns=["feature_origin"]).column(0).to_pylist()
        prior = Counter(before)
        original_path = ROOT / "output-v0.9/upstream-orbitrap" / (record["name"] + ".parquet")
        original = pq.read_table(original_path, columns=["scan_count"]).column(0).to_numpy()
        groups[record["group"]].update(
            {
                "files": 1,
                "before_detected": prior["ms1_detected"],
                "after_detected": counts["ms1_detected"],
                "before_seeded": prior["ms2_seeded"],
                "after_seeded": counts["ms2_seeded"],
                "original_at_least_5_scans": int((original >= 5).sum()),
            }
        )
        rows.append(
            dict(
                file=record["name"],
                group=record["group"],
                input_sha256=record["sha256"],
                features=result.features,
                counts=dict(counts),
                before_counts=dict(prior),
                elapsed_seconds=result.elapsed_seconds,
                timings_seconds=metadata["timings_seconds"],
                parquet_bytes=path.stat().st_size,
                before_parquet_bytes=before_path.stat().st_size,
                tsv_bytes=tsv.stat().st_size,
            )
        )
        if len(rows) % 10 == 0:
            print("Parquet/TSV files verified:", len(rows), flush=True)
    controls = [r for r in manifest if r["group"] in {"blank", "QC"}]
    equality = []
    for record in controls:
        _, path, _ = run(record, False)
        old = ROOT / "eic-review/v0.10/parquet" / (record["name"] + ".parquet")
        count = 0
        for a, b in zip_longest(iter_feature_rows(old), iter_feature_rows(path)):
            assert a is not None and b is not None
            assert a == {key: b[key] for key in a}, (record["name"], a["feature_id"])
            count += 1
        equality.append(dict(file=record["name"], all_previous_columns_equal=True, rows=count))
    # Independent serial timing, alternating order; equality/audit work excluded.
    timing = []
    for repeat in range(3):
        for guarded in [False, True] if repeat % 2 == 0 else [True, False]:
            elapsed = 0.0
            stages = Counter()
            for record in controls:
                result, _, metadata = run(record, guarded)
                elapsed += result.elapsed_seconds
                stages.update(metadata["timings_seconds"])
            timing.append(
                dict(repeat=repeat + 1, guard=guarded, seconds=elapsed, stages_seconds=dict(stages))
            )
            print("Timed blank/QC pass:", timing[-1], flush=True)
    median_on = statistics.median(r["seconds"] for r in timing if r["guard"])
    median_off = statistics.median(r["seconds"] for r in timing if not r["guard"])
    report = dict(
        version=__version__,
        implementation_sha256=implementation_id()[0],
        config=Config().to_dict(),
        files=len(rows),
        all_parquet_tsv_rows_equal=True,
        group_counts=dict(groups),
        rows=rows,
        guard_off_v010_equality=equality,
        serial_50_file_processing_seconds=sum(r["elapsed_seconds"] for r in rows),
        serial_50_file_timing_note="One warm default pass; validation and comparisons excluded; not repeated.",
        blank_qc_timings=timing,
        blank_qc_median_seconds_guard_on=median_on,
        blank_qc_median_seconds_guard_off=median_off,
        blank_qc_timing_note="Three warm serial passes of all 10 blanks/QCs per setting, alternating order, same v0.11 code; guard-off includes the new shape reporting. Scientific outputs differ.",
        limitations=[
            "One positive-mode centroid mzXML cohort; MS2/provenance correctness is not chemical peak ground truth.",
            "A shape guard can merge poorly sampled, weak or irregular real peaks in noisy traces.",
            "Original MassCube uses a 50,000 intensity floor and broader mass tolerance; counts are not an equal-sensitivity comparison.",
        ],
    )
    Path("benchmarks/results/peak-shapes-v0.11.0.json").write_text(
        json.dumps(report, indent=2) + "\n"
    )
    print(
        json.dumps(
            dict(
                groups=report["group_counts"],
                rows=sum(r["features"] for r in rows),
                seconds=report["serial_50_file_processing_seconds"],
                timing_ratio=median_off / median_on,
            )
        )
    )


if __name__ == "__main__":
    main()
