"""Compare 3/5-SD regional calibration with archived v0.11 feature counts.

Run the independent validate_ms2_coverage.py audit first. This serial benchmark
verifies typed Parquet/TSV equality and records per-region evidence/fallbacks.
Raw files and complete sidecars remain under validation/.
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
OUT = ROOT / "regional-tuning-v0.12"


def run(record, sigma):
    source = ROOT / "MSV000086989" / record["name"]
    target = OUT / ("parquet" if sigma == 3 else "sigma5") / (record["name"] + ".parquet")
    result = process_file(
        source, target, config=replace(Config(), tuning_mass_sigma=sigma), overwrite=True
    )
    assert result.status == "ok", result
    meta = json.loads(Path(str(target) + ".json").read_text())
    return result, target, meta


def main():
    manifest = json.loads((ROOT / "MSV000086989/download_manifest.json").read_text())
    audit = json.loads(Path(f"benchmarks/results/MSV000086989-ms2-v{__version__}.json").read_text())
    assert audit["files"] == 50 and audit["implementation_sha256"] == implementation_id()[0]
    warmup()
    run(manifest[0], 3)
    rows, groups = [], defaultdict(Counter)
    for record in manifest:
        result, path, meta = run(record, 3)
        tsv = ROOT / "MSV000086989" / f"v{__version__}-ms2" / (record["name"] + ".tsv")
        counts = Counter()
        for actual, reference in zip_longest(iter_feature_rows(path), iter_feature_rows(tsv)):
            assert actual == reference, (record["name"], actual and actual["feature_id"])
            counts[actual["feature_origin"]] += 1
        old = ROOT / "peak-shape-v0.11/parquet" / (record["name"] + ".parquet")
        prior = Counter(pq.read_table(old, columns=["feature_origin"]).column(0).to_pylist())
        five_result, five_path, five_meta = run(record, 5)
        five_counts = Counter(
            pq.read_table(five_path, columns=["feature_origin"]).column(0).to_pylist()
        )
        tune = meta["parameter_tuning"]["polarities"]["positive"]
        five_tune = five_meta["parameter_tuning"]["polarities"]["positive"]
        region_details = []
        for a, b in zip(tune["regions"], five_tune["regions"]):
            region_details.append(
                {
                    k: a[k]
                    for k in (
                        "lower",
                        "upper",
                        "representative_count",
                        "ms2_representative_count",
                        "ms1_status",
                        "ms2_status",
                        "effective_parameters",
                    )
                }
            )
            region_details[-1]["sigma5_effective_parameters"] = b["effective_parameters"]
            for key, value in a["effective_parameters"].items():
                assert b["effective_parameters"][key] >= value
        groups[record["group"]].update(
            {
                "files": 1,
                "v0.11_detected": prior["ms1_detected"],
                "sigma3_detected": counts["ms1_detected"],
                "sigma5_detected": five_counts["ms1_detected"],
                "v0.11_seeded": prior["ms2_seeded"],
                "sigma3_seeded": counts["ms2_seeded"],
                "sigma5_seeded": five_counts["ms2_seeded"],
            }
        )
        rows.append(
            {
                "file": record["name"],
                "group": record["group"],
                "input_sha256": record["sha256"],
                "features": result.features,
                "sigma5_features": five_result.features,
                "counts": dict(counts),
                "sigma5_counts": dict(five_counts),
                "before_counts": dict(prior),
                "elapsed_seconds": result.elapsed_seconds,
                "timings_seconds": meta["timings_seconds"],
                "parquet_bytes": path.stat().st_size,
                "global_status": tune["status"],
                "global_ms2_status": tune["ms2_status"],
                "representative_count": tune["representative_count"],
                "ms2_representative_count": tune["ms2_representative_count"],
                "range": tune["acquired_mz_range"],
                "range_source": tune["acquired_mz_range_source"],
                "shape_rejections": tune["rejected_candidates"].get("poor_gaussian_shape", 0),
                "regions": region_details,
            }
        )
        if len(rows) % 10 == 0:
            print("Regional calibration / table equality verified:", len(rows), flush=True)
    controls = [r for r in manifest if r["group"] in {"blank", "QC"}]
    repeats = []
    for _ in range(3):
        samples = [run(record, 3)[2] for record in controls]
        repeats.append(
            {
                "elapsed_seconds": sum(s["elapsed_seconds"] for s in samples),
                "timings_seconds": {
                    key: sum(s["timings_seconds"][key] for s in samples)
                    for key in samples[0]["timings_seconds"]
                },
            }
        )
    report = {
        "version": __version__,
        "implementation_sha256": implementation_id()[0],
        "dependencies": implementation_id()[1],
        "config": Config().to_dict(),
        "files": len(rows),
        "groups": dict(groups),
        "rows": rows,
        "source_audit": f"MSV000086989-ms2-v{__version__}.json",
        "ms2_audit_totals": audit["totals"],
        "single_pass_50_files_seconds": sum(r["elapsed_seconds"] for r in rows),
        "blank_qc_three_serial_repeats": repeats,
        "blank_qc_median_seconds": statistics.median(r["elapsed_seconds"] for r in repeats),
        "limitations": [
            "One positive-mode centroid mzXML dataset; shape/tolerance behavior is not a chemical-identity validation.",
            "Sparse high-m/z calibration uses explicit fallbacks; five regions do not imply five supported fits.",
            "v0.11 feature counts are archived controls, not a simultaneous timing comparison.",
            "Mean plus sample SD is applied to screened per-peak Q95 absolute errors; no normal-distribution coverage claim.",
        ],
    }
    target = Path(f"benchmarks/results/regional-tuning-v{__version__}.json")
    target.write_text(json.dumps(report, indent=2, allow_nan=False) + "\n")
    print(
        json.dumps(
            {
                "groups": report["groups"],
                "50_files_seconds": report["single_pass_50_files_seconds"],
                "10_file_median_seconds": report["blank_qc_median_seconds"],
            }
        ),
        flush=True,
    )


if __name__ == "__main__":
    main()
