"""Paired global-128 versus regional-128 calibration and full-extraction benchmark.

Run after the independent MS2 audit and tests finish; no other benchmark should
run concurrently. The baseline source is an unchanged v0.13 package copy.
Alternate version order across three serial 50-file passes, with warmup excluded.
Use --analyze-only to rebuild the report from completed pass measurements.
"""

import argparse
from collections import Counter, defaultdict
from itertools import zip_longest
import json
import os
from pathlib import Path
import platform
import statistics
import subprocess
import sys


ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "validation/regional-sampling-v0.14"
DATA = ROOT / "validation/MSV000086989"


def measure_pass(variant, repeat):
    # Import only in the worker after its isolated PYTHONPATH has been selected.
    from masscubex import Config, __version__, process_file
    from masscubex.features import warmup
    from masscubex.pipeline import implementation_id

    manifest = json.loads((DATA / "download_manifest.json").read_text())
    config = Config()
    assert config.tuning_min_gaussian_similarity == 0.9 and config.tuning_min_peaks == 20
    assert config.tuning_peak_limit == 128 and config.tuning_mz_breaks == (200, 500)
    warmup()
    # Exercise both a blank and a QC, including calibrated gap-fill paths.
    for name in ("Blank_00.mzXML", "QC_Mix_01.mzXML"):
        process_file(DATA / name, OUT / variant / "warmup.parquet", overwrite=True)
    rows = []
    for record in manifest:
        target = OUT / variant / (record["name"] + ".parquet")
        result = process_file(DATA / record["name"], target, overwrite=True)
        assert result.status == "ok", result
        meta = json.loads(Path(str(target) + ".json").read_text())
        tune = meta["parameter_tuning"]["polarities"]["positive"]
        assert meta["ms2_assigned"] == meta["ms2_scans"] and meta["ms2_unassigned"] == 0
        rows.append(
            {
                "file": record["name"],
                "group": record["group"],
                "input_sha256": record["sha256"],
                "features": result.features,
                "seeded": meta["ms2_seeded_features"],
                "ms2_assigned": meta["ms2_assigned"],
                "output_sha256": meta["output_sha256"],
                "elapsed_seconds": result.elapsed_seconds,
                "timings_seconds": meta["timings_seconds"],
                "global_status": tune["status"],
                "global_ms2_status": tune["ms2_status"],
                "global_representatives": tune["representative_count"],
                "regions": [
                    {
                        k: r[k]
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
                    for r in tune["regions"]
                ],
                "pilot_sampling_regions": tune.get("pilot_sampling_regions"),
            }
        )
    report = {
        "variant": variant,
        "repeat": repeat,
        "version": __version__,
        "implementation_sha256": implementation_id()[0],
        "dependencies": implementation_id()[1],
        "config": config.to_dict(),
        "rows": rows,
        "elapsed_seconds": sum(r["elapsed_seconds"] for r in rows),
        "timings_seconds": dict(sum((Counter(r["timings_seconds"]) for r in rows), Counter())),
    }
    (OUT / f"{variant}-pass-{repeat}.json").write_text(json.dumps(report, indent=2) + "\n")
    print(f"{variant}, pass {repeat}: {report['elapsed_seconds']:.3f} s", flush=True)


def analyze(repeats):
    import numpy as np
    import pyarrow.parquet as pq
    from masscubex import __version__, iter_feature_rows
    from masscubex.pipeline import implementation_id

    audit = json.loads(
        (ROOT / f"benchmarks/results/MSV000086989-ms2-v{__version__}.json").read_text()
    )
    assert audit["implementation_sha256"] == implementation_id()[0]
    identity = json.loads((OUT / "baseline/identity.json").read_text())
    expected = {"global": identity["implementation_sha256"], "regional": implementation_id()[0]}
    reports, summary = {}, {}
    for variant in ("global", "regional"):
        passes = [
            json.loads((OUT / f"{variant}-pass-{i}.json").read_text()) for i in range(repeats)
        ]
        assert all(p["implementation_sha256"] == expected[variant] for p in passes)
        first = passes[0]
        for p in passes[1:]:
            assert [(r["file"], r["output_sha256"]) for r in p["rows"]] == [
                (r["file"], r["output_sha256"]) for r in first["rows"]
            ]
        totals, groups = Counter(), defaultdict(Counter)
        region_counts, reference_counts = [], np.zeros(3, dtype=int)
        supported_counts = np.zeros(3, dtype=int)
        for row in first["rows"]:
            name = row["file"]
            table = pq.read_table(
                ROOT / "validation/calibration-coverage-v0.13/three" / (name + ".parquet"),
                columns=["mz", "feature_origin"],
            )
            detected = np.array(table["feature_origin"].to_pylist()) == "ms1_detected"
            masses = table["mz"].to_numpy()[detected]
            counts = np.bincount(np.searchsorted([200, 500], masses, side="right"), minlength=3)
            supported = np.array([r["ms1_status"] == "tuned" for r in row["regions"]])
            reference_counts += counts
            supported_counts += counts * supported
            region_counts.append([r["representative_count"] for r in row["regions"]])
            row["reference_features_by_region"] = counts.tolist()
            row["reference_in_supported_region"] = int((counts * supported).sum())
            totals.update(
                files=1,
                features=row["features"],
                detected=row["features"] - row["seeded"],
                seeded=row["seeded"],
                ms2_assigned=row["ms2_assigned"],
            )
            groups[row["group"]].update(
                files=1, detected=row["features"] - row["seeded"], seeded=row["seeded"]
            )
            if variant == "global":
                old_meta = json.loads(
                    (
                        ROOT
                        / "validation/calibration-coverage-v0.13/three"
                        / (name + ".parquet.json")
                    ).read_text()
                )
                assert old_meta["output_sha256"] == row["output_sha256"], name
            else:
                tsv = DATA / f"v{__version__}-ms2" / (name + ".tsv")
                for actual, reference in zip_longest(
                    iter_feature_rows(OUT / variant / (name + ".parquet")), iter_feature_rows(tsv)
                ):
                    assert actual == reference, (name, actual and actual["feature_id"])
        summary[variant] = {
            "version": first["version"],
            "implementation_sha256": expected[variant],
            "totals": dict(totals),
            "groups": dict(groups),
            "global_fit_files": sum(r["global_status"] == "tuned" for r in first["rows"]),
            "global_precursor_fit_files": sum(
                r["global_ms2_status"] == "tuned" for r in first["rows"]
            ),
            "regional_ms1_fit_files": [
                sum(r["regions"][i]["ms1_status"] == "tuned" for r in first["rows"])
                for i in range(3)
            ],
            "regional_ms2_fit_files": [
                sum(r["regions"][i]["ms2_status"] == "tuned" for r in first["rows"])
                for i in range(3)
            ],
            "median_regional_representatives": np.median(region_counts, axis=0).tolist(),
            "reference_detected_by_region": reference_counts.tolist(),
            "reference_in_supported_regions": supported_counts.tolist(),
            "reference_coverage_fraction": float(supported_counts.sum() / reference_counts.sum()),
            "pass_seconds": [p["elapsed_seconds"] for p in passes],
            "median_seconds": statistics.median(p["elapsed_seconds"] for p in passes),
            "stage_median_seconds": {
                key: statistics.median(p["timings_seconds"][key] for p in passes)
                for key in first["timings_seconds"]
            },
            "rows": first["rows"],
        }
        reports[variant] = passes
    assert reports["global"][0]["config"] == reports["regional"][0]["config"]
    report = {
        "version": __version__,
        "implementation_sha256": implementation_id()[0],
        "dependencies": implementation_id()[1],
        "config": reports["regional"][0]["config"],
        "baseline_commit": identity["commit"],
        "platform": platform.platform(),
        "python": platform.python_version(),
        "repeats": repeats,
        "order": [
            (["global", "regional"] if i % 2 == 0 else ["regional", "global"])
            for i in range(repeats)
        ],
        "summary": summary,
        "passes": {
            key: [
                {k: p[k] for k in ("repeat", "elapsed_seconds", "timings_seconds")} for p in value
            ]
            for key, value in reports.items()
        },
        "ms2_audit_totals": audit["totals"],
        "limitations": [
            "Single positive centroid-mzXML cohort; synthetic tests cover other formats and polarity switching.",
            "Warm serial file processing; excludes startup, JIT warmup, downloads and independent audits.",
            "Coverage uses the same archived v0.13 detected feature masses for both variants; it does not establish chemical identity or absolute mass accuracy.",
            "More calibration anchors also change learned time settings; feature-count differences cannot be attributed solely to mass tolerances.",
        ],
    }
    (ROOT / f"benchmarks/results/regional-sampling-v{__version__}.json").write_text(
        json.dumps(report, indent=2, allow_nan=False) + "\n"
    )
    print(
        json.dumps(
            {
                variant: {k: v for k, v in details.items() if k != "rows"}
                for variant, details in summary.items()
            },
            indent=2,
        ),
        flush=True,
    )


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--worker", choices=("global", "regional"))
    parser.add_argument("--pass-index", type=int, default=0)
    parser.add_argument("--repeats", type=int, default=3)
    parser.add_argument("--analyze-only", action="store_true")
    args = parser.parse_args()
    if args.worker:
        measure_pass(args.worker, args.pass_index)
        return
    assert args.repeats >= 3, "Use at least three paired passes"
    if not args.analyze_only:
        for index in range(args.repeats):
            for variant in ["global", "regional"] if index % 2 == 0 else ["regional", "global"]:
                source = OUT / "baseline/src" if variant == "global" else ROOT / "src"
                env = {**os.environ, "PYTHONPATH": str(source)}
                subprocess.run(
                    [sys.executable, __file__, "--worker", variant, "--pass-index", str(index)],
                    env=env,
                    cwd=ROOT,
                    check=True,
                )
    analyze(args.repeats)


if __name__ == "__main__":
    main()
