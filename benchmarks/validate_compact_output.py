"""Compare schema-6 output to the independently retained full-precision baseline."""

import argparse
from collections import Counter
from itertools import zip_longest
import json
from pathlib import Path

from masscubex import __version__, iter_feature_rows
from masscubex.pipeline import implementation_id


REMOVED = {
    "intensity",
    "quant_method",
    "best_ms2_selection_method",
    "best_ms2_score",
    "apex_threshold_ratio",
}
MASS_FIELDS = {"mz", "seed_precursor_mz", "gapfill_mz_radius_da", "neutral_mass", "fragment_mz"}


def without_operational_fields(value, ignored):
    if isinstance(value, dict):
        return {
            k: without_operational_fields(v, ignored) for k, v in value.items() if k not in ignored
        }
    if isinstance(value, list):
        return [without_operational_fields(item, ignored) for item in value]
    return value


def expected_row(row):
    """Independent specification: all remaining scalar floats use 2 decimals,
    except masses (5) and evidence correlation (6). IDs/lists stay exact.
    """
    if isinstance(row, list):
        return [expected_row(item) for item in row]
    if not isinstance(row, dict):
        return row
    return {
        name: round(value, 5 if name in MASS_FIELDS else 6 if name == "correlation" else 2)
        if isinstance(value, float)
        else expected_row(value)
        for name, value in row.items()
        if name not in REMOVED
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--reference", type=Path, required=True)
    parser.add_argument("--candidate", type=Path, required=True)
    parser.add_argument("--tsv", type=Path, required=True)
    parser.add_argument("--reference-tsv", type=Path, required=True)
    parser.add_argument("--report", type=Path, required=True)
    args = parser.parse_args()
    names = {p.name for p in args.reference.glob("*.parquet")}
    assert names and names == {p.name for p in args.candidate.glob("*.parquet")}
    totals = Counter()
    rows = []
    for name in sorted(names):
        old, new = args.reference / name, args.candidate / name
        tsv = args.tsv / name.removesuffix(".parquet")
        tsv = tsv.with_name(tsv.name + ".tsv")
        old_tsv = args.reference_tsv / tsv.name
        counts = Counter()
        fallback_ids = []
        for baseline, parquet, text in zip_longest(
            iter_feature_rows(old), iter_feature_rows(new), iter_feature_rows(tsv)
        ):
            assert baseline is not None and parquet is not None and text is not None, name
            assert expected_row(baseline) == parquet == text, (name, baseline["feature_id"])
            if baseline["best_ms2_selection_method"] == "apex_fallback":
                fallback_ids.append(baseline["feature_id"])
            counts["features"] += 1
            counts[parquet["feature_origin"]] += 1
            counts["associated_ms2"] += parquet["ms2_count"]
            counts["ion_relationship_entries"] += len(parquet["ion_relations"])
        old_meta = json.loads(Path(str(old) + ".json").read_text())
        for output in (new, tsv):
            meta = json.loads(Path(str(output) + ".json").read_text())
            assert meta["best_ms2_apex_fallback_feature_ids"] == fallback_ids
            # Every scientific diagnostic and effective parameter must match.
            ignored = {
                "schema_version",
                "version",
                "implementation_sha256",
                "elapsed_seconds",
                "timings_seconds",
                "output_format",
                "output_bytes",
                "output_sha256",
                "tsv_bytes",
                "tsv_sha256",
                "output_precision",
                "best_ms2_apex_fallback_feature_ids",
            }
            assert without_operational_fields(meta, ignored) == without_operational_fields(
                old_meta, ignored
            ), name
        sizes = {
            "parquet_before_bytes": old.stat().st_size,
            "parquet_after_bytes": new.stat().st_size,
            "tsv_before_bytes": old_tsv.stat().st_size,
            "tsv_after_bytes": tsv.stat().st_size,
        }
        totals.update(counts)
        totals.update(sizes)
        rows.append({"file": name.removesuffix(".parquet"), **counts, **sizes})
        print(json.dumps({"file": name, "verified_rows": counts["features"]}), flush=True)
    report = {
        "version": __version__,
        "implementation_sha256": implementation_id()[0],
        "files": len(rows),
        "totals": dict(totals),
        "rows": rows,
        "criteria": "Every row equals the old full-precision row after only specified rounding and five column removals; all scientific sidecar fields unchanged.",
    }
    args.report.write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps(report["totals"]))


if __name__ == "__main__":
    main()
