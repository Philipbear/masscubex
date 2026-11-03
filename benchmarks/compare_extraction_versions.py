"""Compare full feature tables and scientific sidecars across numerical engines.

Identifiers, order, bounds, heights and scan mappings must match exactly.
Allow only explicitly listed float64 rounding differences, never missing rows.
"""

import argparse
from collections import Counter, defaultdict
from itertools import zip_longest
import json
from pathlib import Path

from masscubex import __version__
from masscubex.output import iter_feature_rows
from masscubex.pipeline import implementation_id


ROUNDING_FIELDS = {
    "mz",
    "peak_area",
    "intensity",
    "apex_noise_threshold",
    "apex_threshold_ratio",
    "neutral_mass",
}
OPERATIONAL_FIELDS = {
    "version",
    "implementation_sha256",
    "elapsed_seconds",
    "timings_seconds",
    "output_sha256",
    "output_bytes",
    "scan_cache",
}


def scientific_metadata(value):
    if isinstance(value, dict):
        return {
            k: scientific_metadata(v)
            for k, v in value.items()
            if k not in OPERATIONAL_FIELDS and k != "scan_cache_mb"
        }
    if isinstance(value, list):
        return [scientific_metadata(v) for v in value]
    return value


class Comparison:
    def __init__(self):
        self.differences = Counter()
        self.absolute = defaultdict(float)
        self.relative = defaultdict(float)
        self.failures = []
        self.failure_count = 0

    def check(self, actual, expected, context, key="row"):
        if actual == expected:
            return
        if (
            isinstance(actual, dict)
            and isinstance(expected, dict)
            and actual.keys() == expected.keys()
        ):
            for field in actual:
                self.check(actual[field], expected[field], context, key + "." + field)
            return
        if (
            isinstance(actual, (list, tuple))
            and isinstance(expected, (list, tuple))
            and len(actual) == len(expected)
        ):
            for a, b in zip(actual, expected):
                self.check(a, b, context, key)
            return
        if isinstance(actual, float) and isinstance(expected, float):
            error = abs(actual - expected)
            self.differences[key] += 1
            self.absolute[key] = max(self.absolute[key], error)
            self.relative[key] = max(self.relative[key], error / max(abs(expected), 1e-300))
            field = key.rsplit(".", 1)[-1]
            limit = (
                2e-7
                if field == "mass_error_ppm"
                else 3e-14 * max(abs(expected), 1.0)
                if field in ROUNDING_FIELDS or key.startswith("metadata.")
                else 0
            )
            if error <= limit:
                return
        self.failure_count += 1
        if len(self.failures) < 30:
            self.failures.append(
                {"context": context, "field": key, "actual": actual, "expected": expected}
            )


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--reference", type=Path, required=True)
    parser.add_argument("--candidate", type=Path, required=True)
    parser.add_argument("--report", type=Path, required=True)
    args = parser.parse_args()
    reference_files = {p.name for p in args.reference.glob("*.parquet")}
    candidate_files = {p.name for p in args.candidate.glob("*.parquet")}
    if not reference_files or reference_files != candidate_files:
        raise ValueError("Directories must contain the same nonempty set of Parquet tables")
    comparison = Comparison()
    totals = Counter()
    for name in sorted(reference_files):
        source, output = args.reference / name, args.candidate / name
        for index, (actual, expected) in enumerate(
            zip_longest(iter_feature_rows(output), iter_feature_rows(source)), 1
        ):
            comparison.check(actual, expected, f"{name}:{index}")
            if actual is not None:
                totals["features"] += 1
                totals[actual["feature_origin"]] += 1
                totals["associated_ms2"] += actual["ms2_count"]
                totals["ion_relationship_entries"] += len(actual["ion_relations"])
        old_meta = scientific_metadata(json.loads(Path(str(source) + ".json").read_text()))
        new_meta = scientific_metadata(json.loads(Path(str(output) + ".json").read_text()))
        comparison.check(new_meta, old_meta, name, "metadata")
    report = {
        "version": __version__,
        "implementation_sha256": implementation_id()[0],
        "reference": str(args.reference),
        "candidate": str(args.candidate),
        "files": len(reference_files),
        **totals,
        "failure_count": comparison.failure_count,
        "failures": comparison.failures,
        "differing_float_values": dict(comparison.differences),
        "max_absolute_difference": dict(comparison.absolute),
        "max_relative_difference": dict(comparison.relative),
        "criteria": {
            "measurement_rtol": 3e-14,
            "measurement_atol_near_zero": 3e-14,
            "mass_error_ppm_atol": 2e-7,
            "identity_boundaries_heights_scan_lists": "exact",
        },
    }
    args.report.write_text(json.dumps(report, indent=2, allow_nan=False) + "\n")
    print(
        json.dumps({"files": len(reference_files), **totals, "failures": comparison.failure_count})
    )
    if comparison.failure_count:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
