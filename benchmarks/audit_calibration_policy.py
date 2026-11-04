"""Local policy-5 calibration comparison; not a release-validation workflow.

Use the existing policy-4 sidecars as the scientific baseline. Recalibrate the
same 50 files and smoke-check extraction on a blank, QC and the affected sample.
The one-pass elapsed time is diagnostic only, not a paired speed benchmark.
"""

import json
from pathlib import Path

import numpy as np

from masscubex import Config, process_file
from masscubex.pipeline import implementation_id
from masscubex.tuning import calibrate_file


ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "validation/MSV000086989"
BASELINE = ROOT / "validation/regional-sampling-v0.14/regional"
OUT = ROOT / "validation/calibration-policy5"


def main():
    OUT.mkdir(exist_ok=True)
    manifest = json.loads((DATA / "download_manifest.json").read_text())
    config = Config()
    calibrate_file(DATA / "Blank_00.mzXML", config)  # Exclude import/JIT warmup.
    rows = []
    for number, record in enumerate(manifest, 1):
        old = json.loads((BASELINE / (record["name"] + ".parquet.json")).read_text())
        assert old["parameter_tuning"]["policy_version"] == 4
        result = calibrate_file(DATA / record["name"], config)
        rows.append(
            {
                "file": record["name"],
                "group": record["group"],
                "input_sha256": record["sha256"],
                "before": old["parameter_tuning"],
                "after": result.diagnostics,
            }
        )
        if number % 10 == 0:
            print(f"Calibrated {number}/{len(manifest)}", flush=True)
    summary = {}
    for mode in ("before", "after"):
        info = [row[mode]["polarities"]["positive"] for row in rows]
        summary[mode] = {
            "independent_ms1_fits": [
                sum(d["regions"][i]["ms1_status"] == "tuned" for d in info) for i in range(3)
            ],
            "median_anchors": [
                float(np.median([d["regions"][i]["representative_count"] for d in info]))
                for i in range(3)
            ],
            "min_anchors": [
                min(d["regions"][i]["representative_count"] for d in info) for i in range(3)
            ],
        }
    smoke = []
    for name in ("Blank_00.mzXML", "QC_Mix_01.mzXML", "NIST_POS_Samp_14-03.mzXML"):
        target = OUT / (name + ".parquet")
        result = process_file(DATA / name, target, config=config, overwrite=True)
        assert result.status == "ok", result
        meta = json.loads(Path(str(target) + ".json").read_text())
        old = json.loads((BASELINE / (name + ".parquet.json")).read_text())
        assert meta["ms2_assigned"] == meta["ms2_scans"] == old["ms2_scans"]
        assert meta["ms2_unassigned"] == 0
        assert meta["parameter_tuning"]["polarities"] == next(
            row["after"]["polarities"] for row in rows if row["file"] == name
        )
        smoke.append(
            {
                "file": name,
                "before_features": old["features"],
                "after_features": meta["features"],
                "before_seeded": old["ms2_seeded_features"],
                "after_seeded": meta["ms2_seeded_features"],
                "ms2_assigned": meta["ms2_assigned"],
            }
        )
    report = {
        "implementation_sha256": implementation_id()[0],
        "scope": "50 calibration comparisons and 3 extraction smoke checks; not release validation",
        "config": config.to_dict(),
        "summary": summary,
        "smoke": smoke,
        "files": rows,
    }
    (OUT / "calibration.json").write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps({"summary": summary, "smoke": smoke}, indent=2))


if __name__ == "__main__":
    main()
