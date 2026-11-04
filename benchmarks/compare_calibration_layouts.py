"""Compare the requested three/four mass layouts without relaxing quality rules.

Use one common set of archived v0.12 detected features to measure how much data
each layout places in supported regions. Cross-validation holds out each strong
candidate once, testing its Q95 error summary (not individual raw observations).
The audit of alternative pilot samplers is exploratory and not a runtime change.
"""

from collections import Counter
from dataclasses import replace
from itertools import zip_longest
import json
from pathlib import Path

import numpy as np
import pyarrow.parquet as pq

from masscubex import Config, __version__, process_file, iter_feature_rows
from masscubex.features import warmup
from masscubex.mass_tolerance import mass_radius
from masscubex.pipeline import implementation_id
from masscubex.tuning import Representative, _screen_representatives, derive_settings


ROOT = Path("validation")
OUT = ROOT / "calibration-coverage-v0.13"
LAYOUTS = {"three": (200.0, 500.0), "four": (200.0, 400.0, 600.0)}


def threshold_counts(candidates, breaks):
    result = {}
    for shape in (0.8, 0.9):
        qualified = [p for p in candidates if p.gaussian_similarity >= shape]
        global_peaks, _ = _screen_representatives(qualified, Config())
        counts = []
        for index in range(len(breaks) + 1):
            selected, _ = _screen_representatives(
                [p for p in qualified if np.searchsorted(breaks, p.mz, side="right") == index],
                Config(),
            )
            counts.append(len(selected))
        result[str(shape)] = {
            "global_accepted": len(global_peaks),
            "regional_accepted": counts,
            "supported": {
                str(n): [len(global_peaks) >= 20 and count >= n for count in counts]
                for n in (10, 20)
            },
        }
    return result


def holdout(candidates, breaks):
    peaks = sorted(
        [p for p in candidates if p.gaussian_similarity >= 0.9], key=lambda p: (p.rt, p.mz)
    )
    totals = Counter()
    for fold in range(5):
        train = [p for i, p in enumerate(peaks) if i % 5 != fold]
        test = [p for i, p in enumerate(peaks) if i % 5 == fold]
        # Cadence is used only as an availability gate for this mass-only check.
        config, info = derive_settings(
            train, [1.0] * 100, [], Config(tuning_mz_breaks=breaks), "positive"
        )
        for peak in test:
            index = int(np.searchsorted(breaks, peak.mz, side="right"))
            local = info["regions"][index]["ms1_status"] == "tuned"
            covered = peak.mass_error_da <= mass_radius(peak.mz, config)
            totals.update(
                held_out=1,
                covered=int(covered),
                local_model=int(local),
                locally_covered=int(local and covered),
            )
    return dict(totals)


def main():
    manifest = json.loads((ROOT / "MSV000086989/download_manifest.json").read_text())
    pilot = json.loads((OUT / "pilot-audit.json").read_text())
    candidates_by_file = {row["file"]: row for row in pilot["rows"]}
    raw_audit = json.loads(
        Path(f"benchmarks/results/MSV000086989-ms2-v{__version__}.json").read_text()
    )
    assert raw_audit["implementation_sha256"] == implementation_id()[0]
    warmup()
    process_file(
        ROOT / "MSV000086989" / manifest[0]["name"], OUT / "warmup.parquet", overwrite=True
    )
    rows, totals = [], {name: Counter() for name in LAYOUTS}
    for record in manifest:
        name = record["name"]
        pilot_row = candidates_by_file[name]
        candidates = [Representative(**p) for p in pilot_row["modes"]["global128"]["candidates"]]
        # The scientific comparison uses identical baseline feature masses.
        table = pq.read_table(
            ROOT / "regional-tuning-v0.12/parquet" / (name + ".parquet"),
            columns=["mz", "feature_origin"],
        )
        masses = table["mz"].to_numpy()[
            np.array(table["feature_origin"].to_pylist()) == "ms1_detected"
        ]
        row = {
            "file": name,
            "group": record["group"],
            "input_sha256": record["sha256"],
            "layouts": {},
        }
        for layout, breaks in LAYOUTS.items():
            target = OUT / layout / (name + ".parquet")
            config = replace(Config(), tuning_mz_breaks=breaks)
            result = process_file(
                ROOT / "MSV000086989" / name, target, config=config, overwrite=True
            )
            assert result.status == "ok"
            meta = json.loads(Path(str(target) + ".json").read_text())
            info = meta["parameter_tuning"]["polarities"]["positive"]
            assert meta["ms2_assigned"] == meta["ms2_scans"] and meta["ms2_unassigned"] == 0
            if layout == "three":
                tsv = ROOT / "MSV000086989" / f"v{__version__}-ms2" / (name + ".tsv")
                for a, b in zip_longest(iter_feature_rows(target), iter_feature_rows(tsv)):
                    assert a == b, (name, a and a["feature_id"])
            thresholds = threshold_counts(candidates, breaks)
            supported = [r["ms1_status"] == "tuned" for r in info["regions"]]
            assert supported == thresholds["0.9"]["supported"]["20"]
            local = np.array(supported)[np.searchsorted(breaks, masses, side="right")]
            measured = Counter(
                pq.read_table(target, columns=["feature_origin"]).column(0).to_pylist()
            )
            validation = holdout(candidates, breaks)
            totals[layout].update(
                files=1,
                features=result.features,
                ms2_assigned=meta["ms2_assigned"],
                reference_detected=len(masses),
                reference_in_supported_region=int(local.sum()),
                fitted_regions=sum(supported),
                elapsed_seconds=result.elapsed_seconds,
                **validation,
            )
            row["layouts"][layout] = {
                "counts": dict(measured),
                "features": result.features,
                "elapsed_seconds": result.elapsed_seconds,
                "timings_seconds": meta["timings_seconds"],
                "representative_counts": [r["representative_count"] for r in info["regions"]],
                "local_fits": supported,
                "global_status": info["status"],
                "reference_in_supported_region": int(local.sum()),
                "reference_detected": len(masses),
                "threshold_comparison": thresholds,
                "holdout": validation,
            }
        rows.append(row)
        if len(rows) % 10 == 0:
            print("Layout comparison files:", len(rows), flush=True)
    sampling = {}
    for mode in pilot["modes"]:
        variants = {}
        for shape in ("0.8", "0.9"):
            variants[shape] = {
                str(n): [
                    sum(
                        r["modes"][mode]["summary"][shape]["global_accepted"] >= 20
                        and r["modes"][mode]["summary"][shape]["regions"][i]["accepted"] >= n
                        for r in pilot["rows"]
                    )
                    for i in range(4)
                ]
                for n in (10, 20)
            }
        sampling[mode] = {
            "four_bin_fit_counts": variants,
            "selected_centroids": np.sum(
                [r["modes"][mode]["selected_centroids"] for r in pilot["rows"]], axis=0
            ).tolist(),
            "high_mass_rejections": dict(
                sum((Counter(r["modes"][mode]["rejections"][3]) for r in pilot["rows"]), Counter())
            ),
        }
    summary = {}
    for layout, breaks in LAYOUTS.items():
        summary[layout] = {
            "totals": dict(totals[layout]),
            "breaks": list(breaks),
            "local_fit_files": [
                sum(r["layouts"][layout]["local_fits"][i] for r in rows)
                for i in range(len(breaks) + 1)
            ],
            "threshold_variants": {
                shape: {
                    str(n): [
                        sum(
                            r["layouts"][layout]["threshold_comparison"][shape]["supported"][
                                str(n)
                            ][i]
                            for r in rows
                        )
                        for i in range(len(breaks) + 1)
                    ]
                    for n in (10, 20)
                }
                for shape in ("0.8", "0.9")
            },
        }
    report = {
        "version": __version__,
        "implementation_sha256": implementation_id()[0],
        "config": Config().to_dict(),
        "pilot_audit_implementation_sha256": pilot["implementation_sha256"],
        "layouts": summary,
        "sampling_experiments": sampling,
        "raw_four_bin_centroids": np.sum(
            [r["raw_centroids"] for r in pilot["rows"]], axis=0
        ).tolist(),
        "raw_four_bin_intensity": np.sum(
            [r["raw_intensity"] for r in pilot["rows"]], axis=0
        ).tolist(),
        "rows": rows,
        "limitations": [
            "A supported broader region covers more feature masses; this alone is not proof of more accurate matching.",
            "Holdout tests Q95 absolute-error summaries of strong Gaussian-qualified pilot peaks, not weak/raw-ion or absolute mass accuracy.",
            "Alternative 0.8 shape / 10-peak counts are projected evidence availability; runtime defaults remain 0.9 and 20.",
            "Regional sampling experiments use four bins and are not applied in production.",
            "Single positive-mode centroid mzXML dataset; numerical floors dominate many estimated windows.",
        ],
    }
    output = Path(f"benchmarks/results/calibration-layouts-v{__version__}.json")
    output.write_text(json.dumps(report, indent=2, allow_nan=False) + "\n")
    print(json.dumps(summary), flush=True)


if __name__ == "__main__":
    main()
