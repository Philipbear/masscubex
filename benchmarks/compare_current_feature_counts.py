"""Fresh six-file count comparison against pinned, unmodified MassCube.

Run with upstream-only dependencies on PYTHONPATH. Original source and all
feature tables are temporary; retain only counts, settings and source hashes.
The mzXML metadata adapter is documented in audit_upstream_features.py.
"""

import argparse
from collections import Counter, defaultdict
from dataclasses import replace
import hashlib
import json
from importlib.metadata import version
from pathlib import Path
import subprocess
import sys
from tempfile import TemporaryDirectory
import traceback

import pyarrow.parquet as pq

from audit_upstream_features import audit_file
from masscubex import Config, process_file
from masscubex.features import warmup
from masscubex.pipeline import implementation_id


UPSTREAM_COMMIT = "356748ede691862b9a87f67812db5e728a703b9d"
DEFAULT_FILES = (
    "NIST_POS_Samp_01-01.mzXML",
    "NIST_POS_Samp_14-03.mzXML",
    "Blank_00.mzXML",
    "Blank_Post_14-03.mzXML",
    "QC_Mix_01.mzXML",
    "QC_Mix_04.mzXML",
)


def copy_upstream(destination):
    """Load byte-exact original sources from the existing pinned Git object."""
    names = subprocess.check_output(
        ["git", "ls-tree", "-r", "--name-only", UPSTREAM_COMMIT, "--", "src/masscube"],
        text=True,
    ).splitlines()
    hashes = {}
    for name in names:
        content = subprocess.check_output(["git", "show", f"{UPSTREAM_COMMIT}:{name}"])
        target = destination / name
        target.parent.mkdir(parents=True, exist_ok=True)
        target.write_bytes(content)
        assert target.read_bytes() == content
        hashes[name] = hashlib.sha256(content).hexdigest()
    return hashes


def count_table(path, *, original=False):
    columns = ["scan_count", "peak_height", "ms2_count"]
    if not original:
        columns += ["feature_origin"]
    table = pq.read_table(path, columns=columns).to_pydict()
    result = Counter(total=0, detected=0, seeded=0, detected_at_least_5_scans=0)
    for i, scans in enumerate(table["scan_count"]):
        detected = original or table["feature_origin"][i] == "ms1_detected"
        height = table["peak_height"][i]
        result["total"] += 1
        result["detected" if detected else "seeded"] += 1
        result["detected_at_least_5_scans"] += detected and scans >= 5
        result["features_with_ms2"] += table["ms2_count"][i] > 0
        result["assigned_ms2"] += table["ms2_count"][i]
        result["height_below_50000"] += height is not None and height < 50000
    return dict(result)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, default=Path("validation/MSV000086989"))
    parser.add_argument("--report", type=Path, required=True)
    parser.add_argument("--files", nargs="+", default=DEFAULT_FILES)
    args = parser.parse_args()
    manifest = {
        r["name"]: r for r in json.loads((args.input / "download_manifest.json").read_text())
    }
    config = Config()
    fixed_floor = replace(config, noise_mode="none", min_intensity=50000)
    report = {
        "dataset": "MSV000086989",
        "selection": "Files specified before running; default selection is two samples, two blanks and two QCs.",
        "selected_files": list(args.files),
        "implementation_sha256": implementation_id()[0],
        "dependencies": implementation_id()[1],
        "upstream_only_dependencies": {
            "IsoSpecPy": version("IsoSpecPy"),
            "pyteomics": version("pyteomics"),
        },
        "upstream_commit": UPSTREAM_COMMIT,
        "upstream_version": "1.2.22",
        "masscubex_config": config.to_dict(),
        "floor_control_config": fixed_floor.to_dict(),
        "upstream_settings": {
            "instrument": "orbitrap",
            "ms1_absolute_threshold": 50000,
            "ms1_mz_tolerance_da": 0.01,
            "segmentation_passes": 2,
            "single_file_minimum_observations": None,
        },
        "rows": [],
        "errors": [],
        "limitations": [
            "The pinned upstream rejects mzXML: unchanged arrays and minute RTs pass through a metadata-only adapter to its mzML extractor.",
            "Upstream counting stops after original finalization/MS2 assignment; grouping adds labels without changing the row set.",
            "The >=5-observation count is a common diagnostic, not upstream's single-file output filter.",
            "The fixed 50000 floor is a diagnostic control only; no production settings are changed.",
            "Counts do not establish chemical identity or false-positive rates; this is one positive-mode centroid Orbitrap dataset.",
        ],
    }
    totals = defaultdict(Counter)
    args.report.parent.mkdir(parents=True, exist_ok=True)
    warmup()
    with TemporaryDirectory(prefix="masscubex-count-audit-") as temporary:
        root = Path(temporary)
        upstream = root / "upstream"
        report["upstream_source_hashes"] = copy_upstream(upstream)
        sys.path.insert(0, str(upstream / "src"))
        for name in args.files:
            try:
                entry = manifest[name]
                source = args.input / name
                assert hashlib.sha256(source.read_bytes()).hexdigest() == entry["sha256"]
                original_path = root / (name + ".upstream.parquet")
                original_meta = audit_file(source, original_path, "orbitrap")
                default_path = root / (name + ".current.parquet")
                floor_path = root / (name + ".floor.parquet")
                process_file(source, default_path, config=config)
                process_file(source, floor_path, config=fixed_floor)
                default_meta = json.loads(Path(str(default_path) + ".json").read_text())
                floor_meta = json.loads(Path(str(floor_path) + ".json").read_text())
                default_tuning = default_meta["parameter_tuning"]
                floor_tuning = floor_meta["parameter_tuning"]
                assert default_tuning["polarities"] == floor_tuning["polarities"]
                for polarity, effective in default_tuning["effective_configs"].items():
                    unchanged = {
                        k: v
                        for k, v in effective.items()
                        if k not in {"noise_mode", "min_intensity"}
                    }
                    control = floor_tuning["effective_configs"][polarity]
                    assert unchanged == {
                        k: v for k, v in control.items() if k not in {"noise_mode", "min_intensity"}
                    }
                row = {
                    "file": name,
                    "group": entry["group"],
                    "input_sha256": entry["sha256"],
                    "original": count_table(original_path, original=True),
                    "masscubex": count_table(default_path),
                    "masscubex_floor_50000": count_table(floor_path),
                    "original_diagnostics": original_meta,
                    "retained_centroids": default_meta["retained_centroids"],
                    "floor_retained_centroids": floor_meta["retained_centroids"],
                    "masscubex_stage_seconds": default_meta["timings_seconds"],
                }
                for method in ("original", "masscubex", "masscubex_floor_50000"):
                    totals[method].update(row[method])
                report["rows"].append(row)
                print(json.dumps(row), flush=True)
            except Exception as error:
                traceback.print_exc()
                failure = {"file": name, "error": f"{type(error).__name__}: {error}"}
                report["errors"].append(failure)
                print(json.dumps(failure), flush=True)
            report["totals"] = {key: dict(value) for key, value in totals.items()}
            args.report.write_text(json.dumps(report, indent=2, allow_nan=False) + "\n")
    if report["errors"]:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
