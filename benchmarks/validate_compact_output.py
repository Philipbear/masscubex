"""Historical schema-9 compact-output audit against a full-precision baseline.

Schema 10 changes the splitting/grouping algorithms, so old feature rows are
not a valid numerical baseline. Use fresh paired TSV/Parquet runs for schema 10.
"""

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
    "best_ms2_native_id",
    "best_ms2_spectrum_index",
    "gapfill_status",
    "gapfill_mz_radius_da",
    "best_ms2_fragment_intensity",
    "best_ms2_precursor_intensity_source",
    "best_ms2_precursor_intensity",
    "seed_precursor_mz",
    "seed_rt_seconds",
}
INTENSITY_FIELDS = {
    "peak_height",
    "peak_area",
    "apex_noise_threshold",
    "fragment_relative_intensity",
}
CORRELATION_FIELDS = {"correlation", "gaussian_similarity"}
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
    """Independent precision policy: masses 5, intensities 1, correlations 4,
    and other floats 2 decimals. IDs/lists stay exact.
    """
    if isinstance(row, list):
        return [expected_row(item) for item in row]
    if not isinstance(row, dict):
        return row
    return {
        name: round(
            value,
            5
            if name in MASS_FIELDS
            else 4
            if name in CORRELATION_FIELDS
            else 1
            if name in INTENSITY_FIELDS
            else 2,
        )
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
        old_meta = json.loads(Path(str(old) + ".json").read_text())
        new_meta = json.loads(Path(str(new) + ".json").read_text())
        if new_meta["schema_version"] != 9:
            raise ValueError("This historical baseline comparison requires schema 9 outputs")
        counts = Counter()
        fallback_ids = []
        identifier_fallbacks = []
        for baseline, parquet, text in zip_longest(
            iter_feature_rows(old), iter_feature_rows(new), iter_feature_rows(tsv)
        ):
            assert baseline is not None and parquet is not None and text is not None, name
            expected = expected_row(baseline)
            # Independently derive the new labels from unchanged relationship
            # evidence, without calling the production grouping summarizer.
            forms = set()
            fragment = False
            for link in baseline["ion_relations"]:
                form = (
                    link.get("ion")
                    if link["relation"] == "adduct"
                    else link.get(link["role"] + "_ion")
                )
                if form is not None:
                    forms.add(form)
                fragment |= link.get("role") == "fragment"
            expected["adduct"] = sorted(forms)
            checked = (
                old_meta["config"]["group_ions"]
                and old_meta["config"]["group_isf"]
                and baseline["polarity"] in {"positive", "negative"}
            )
            expected["is_in_source_fragment"] = fragment if checked else None
            assert expected == parquet == text, (name, baseline["feature_id"])
            best = baseline["best_ms2_scan_number"]
            if baseline["ms2_count"] and (
                best is None or baseline["ms2_scan_numbers"].count(best) > 1
            ):
                identifier_fallbacks.append(
                    {
                        "feature_id": baseline["feature_id"],
                        "native_id": baseline["best_ms2_native_id"],
                        "spectrum_index": baseline["best_ms2_spectrum_index"],
                    }
                )
            if baseline["best_ms2_selection_method"] == "apex_fallback":
                fallback_ids.append(baseline["feature_id"])
            counts["features"] += 1
            counts[parquet["feature_origin"]] += 1
            counts["associated_ms2"] += parquet["ms2_count"]
            counts["ion_relationship_entries"] += len(parquet["ion_relations"])
        for output in (new, tsv):
            meta = json.loads(Path(str(output) + ".json").read_text())
            assert meta["best_ms2_apex_fallback_feature_ids"] == fallback_ids
            assert meta["best_ms2_identifier_fallbacks"] == identifier_fallbacks
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
                "best_ms2_identifier_fallbacks",
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
        "criteria": "Every row equals the full-precision baseline after specified rounding, column removals and labels derived from unchanged ion relationships; scientific diagnostics and settings are unchanged.",
    }
    args.report.write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps(report["totals"]))


if __name__ == "__main__":
    main()
