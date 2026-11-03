"""Compare every feature and provenance field across the baseline and both formats.

Run after benchmark_extraction.py. This audit is deliberately outside the timed
benchmark and streams Parquet/TSV rows instead of accumulating a dataset table.
"""

import argparse
from collections import Counter
import hashlib
from itertools import zip_longest
import json
from pathlib import Path

from masscubex import __version__, iter_feature_rows
from masscubex.pipeline import implementation_id


def digest(path):
    result = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            result.update(block)
    return result.hexdigest()


def scientific_metadata(path):
    data = json.loads(Path(str(path) + ".json").read_text())
    for key in (
        "schema_version",
        "version",
        "implementation_sha256",
        "dependencies",
        "elapsed_seconds",
        "timings_seconds",
        "output_format",
        "output_bytes",
        "output_sha256",
        "tsv_bytes",
        "tsv_sha256",
    ):
        data.pop(key, None)
    data["parameter_tuning"].pop("elapsed_seconds", None)
    return data


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--directory", type=Path, default=Path("validation/MSV000086989"))
    parser.add_argument(
        "--baseline", type=Path, default=Path("validation/performance/baseline-tsv")
    )
    parser.add_argument("--tsv", type=Path, default=Path("validation/performance/optimized-tsv"))
    parser.add_argument(
        "--parquet", type=Path, default=Path("validation/performance/optimized-parquet")
    )
    parser.add_argument(
        "--report",
        type=Path,
        default=Path(f"benchmarks/results/extraction-equivalence-v{__version__}.json"),
    )
    args = parser.parse_args()
    manifest = json.loads((args.directory / "download_manifest.json").read_text())
    totals = Counter()
    rows = []
    for entry in manifest:
        name = entry["name"]
        old = args.baseline / (name + ".tsv")
        tsv = args.tsv / (name + ".tsv")
        parquet = args.parquet / (name + ".parquet")
        old_sha, new_sha = digest(old), digest(tsv)
        assert old_sha == new_sha, f"TSV changed: {name}"
        metadata = scientific_metadata(old)
        assert metadata == scientific_metadata(tsv) == scientific_metadata(parquet), name
        counts = Counter()
        for number, (a, b) in enumerate(
            zip_longest(iter_feature_rows(tsv), iter_feature_rows(parquet)), 1
        ):
            assert a == b, f"Parquet round trip changed {name} feature {number}"
            counts["features"] += 1
            counts[a["feature_origin"]] += 1
            counts["assigned_ms2"] += a["ms2_count"]
            counts["ion_relation_rows"] += len(a["ion_relations"])
        assert counts["features"] == metadata["features"]
        assert counts["assigned_ms2"] == metadata["ms2_assigned"]
        totals.update(counts)
        rows.append(
            {
                "name": name,
                **counts,
                "tsv_sha256": new_sha,
                "parquet_sha256": digest(parquet),
                "input_sha256": entry["sha256"],
            }
        )
    report = {
        "version": __version__,
        "implementation_sha256": implementation_id()[0],
        "files": len(rows),
        "totals": dict(totals),
        "all_optimized_tsvs_byte_identical_to_v0_5_1": True,
        "all_parquet_values_exactly_equal_to_tsv": True,
        "scientific_sidecar_metadata_identical": True,
        "scope": "Every scalar, list element, relationship and row order; scientific metadata excludes version, dependency, formatting and timing provenance.",
        "rows": rows,
    }
    args.report.write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps({k: v for k, v in report.items() if k != "rows"}, indent=2))


if __name__ == "__main__":
    main()
