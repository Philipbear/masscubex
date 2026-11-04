"""Independently audit every MS2 mapping and every seeded MS1 measurement.

Raw files/output TSVs remain under validation/. Only the compact report is
published. Pyteomics is used as an independent decoder, not as a runtime reader.
"""

import argparse
from concurrent.futures import ProcessPoolExecutor, as_completed
from collections import Counter
import csv
import hashlib
import json
from pathlib import Path
import platform
import time

import numpy as np

from masscubex import Config, __version__
from validate_repository import process_for_audit, reference_check


def audit_seeded_quantification(path, output, config, full_precision_rows=None):
    from pyteomics import mzxml

    with output.open() as stream:
        rows = (
            full_precision_rows
            if full_precision_rows is not None
            else list(csv.DictReader(stream, delimiter="\t"))
        )
    metadata = json.loads(Path(str(output) + ".json").read_text())
    if metadata.get("schema_version", 0) >= 6 and full_precision_rows is None:
        raise ValueError("Rounded tables need unrounded audit coordinates; use process_for_audit")
    resolved = metadata.get("parameter_tuning", {}).get("effective_configs", {})
    seeded = [r for r in rows if r["feature_origin"] == "ms2_seeded"]
    masses = np.array([float(r["mz"]) for r in seeded])
    start = np.array([float(r["rt_start_seconds"]) for r in seeded])
    end = np.array([float(r["rt_end_seconds"]) for r in seeded])
    traces = [[] for _ in seeded]
    noise_references = [[] for _ in seeded]
    source_ms2 = set()
    known_wide = 0
    missing_width = 0
    with mzxml.MzXML(str(path), use_index=False) as reader:
        for index, scan in enumerate(reader):
            if scan["msLevel"] == 2:
                source_ms2.add(index)
                precursors = scan.get("precursorMz", [])
                width = precursors[0].get("windowWideness") if len(precursors) == 1 else None
                known_wide += width is not None and float(width) > 4
                missing_width += width is None
                continue
            if scan["msLevel"] != 1:
                continue
            rt = float(scan["retentionTime"]) * 60
            polarity = {"+": "positive", "-": "negative"}.get(scan.get("polarity"), "unknown")
            # Evaluate m/z boundaries in float64, as the runtime does. Subtracting
            # a Python float directly from a float32 export rounds the target and
            # can exclude a real peak near the 0.01 Da boundary.
            raw_mz = scan["m/z array"].astype(np.float64)
            raw_y = scan["intensity array"].astype(np.float64)
            unique_mz, inverse = np.unique(raw_mz, return_inverse=True)
            summed_y = np.bincount(inverse, weights=raw_y)
            for i in np.flatnonzero((start <= rt + 1e-9) & (rt <= end + 1e-9)):
                if seeded[i]["polarity"] != polarity:
                    continue
                if not seeded[i]["peak_height"]:
                    continue
                effective = resolved.get(polarity, config.to_dict())
                tolerance = (
                    float(seeded[i]["gapfill_mz_radius_da"])
                    if seeded[i].get("gapfill_mz_radius_da")
                    else max(
                        effective["ms2_mz_tolerance_da"],
                        masses[i] * effective["ms2_mz_tolerance_ppm"] * 1e-6,
                    )
                )
                in_window = (unique_mz >= masses[i] - tolerance) & (
                    unique_mz <= masses[i] + tolerance
                )
                value = float(summed_y[in_window].max()) if in_window.any() else 0.0
                traces[i].append((rt, value))
                if seeded[i]["ms1_quant_status"] == "imputed_noise":
                    # Independent, sorted NumPy reference; do not call the runtime noise kernel.
                    assert (
                        effective["noise_layout"] == "bins"
                        and effective["noise_estimator"] == "low_tail"
                    )
                    width = effective["noise_window_da"]
                    bin_id = np.floor(masses[i] / width)
                    values = np.sort(
                        summed_y[(np.floor(unique_mz / width) == bin_id) & (summed_y > 0)]
                    )
                    level = effective["min_intensity"]
                    if (
                        effective["noise_mode"] != "none"
                        and len(values) >= effective["noise_min_neighbors"]
                    ):
                        count = max(
                            effective["noise_min_neighbors"],
                            int(np.ceil(effective["noise_fraction"] * len(values))),
                        )
                        selected = values[:count]
                        level = max(level, float(selected.mean() + 3 * selected.std(ddof=0)))
                    noise_references[i].append(
                        (abs(rt - float(seeded[i]["seed_rt_seconds"])), rt, level)
                    )
    seen = set()
    for row in rows:
        for index in json.loads(row["ms2_spectrum_indices"]):
            assert index not in seen
            seen.add(index)
    unresolved = {r["spectrum_index"] for r in metadata["unassigned_ms2"]}
    assert not seen & unresolved
    assert seen | unresolved == source_ms2
    for row, trace, references in zip(seeded, traces, noise_references):
        if not trace:
            assert row["peak_height"] == row["peak_area"] == ""
            continue
        t, y = np.asarray(trace).T
        effective = resolved.get(row["polarity"], config.to_dict())
        if "noise_layout" in effective:
            assert float(row["mz"]) == float(row["seed_precursor_mz"])
            assert (
                float(row["gapfill_mz_radius_da"])
                <= max(
                    effective["mz_tolerance_da"],
                    float(row["mz"]) * effective["mz_tolerance_ppm"] * 1e-6,
                )
                + 1e-14
            )
        if row["ms1_quant_status"] == "imputed_noise":
            chosen = [r[2] for r in sorted(references)[:3] if r[2] > 0]
            assert chosen and row["peak_area"] == "" and int(row["scan_count"]) == 0
            threshold = float(np.median(chosen))
            np.testing.assert_allclose(float(row["apex_noise_threshold"]), threshold, rtol=1e-12)
            np.testing.assert_allclose(
                float(row["peak_height"]),
                effective["gapfill_noise_fraction"] * threshold,
                rtol=1e-12,
            )
            if row["gapfill_status"] == "no_local_ms1_signal":
                dt = np.median(np.diff(t)) if len(t) > 1 else 0
                local = abs(t - float(row["seed_rt_seconds"])) <= max(
                    effective["ms2_rt_tolerance"], 2 * dt
                )
                assert np.max(y[local], initial=0.0) == 0
            else:
                assert y.max() == 0
            continue
        dt = np.diff(t)
        expected_area = np.sum(
            np.where(
                dt <= effective.get("gapfill_max_gap_seconds", config.max_gap_seconds),
                dt * (y[:-1] + y[1:]) / 2,
                0,
            )
        )
        np.testing.assert_allclose(float(row["peak_height"]), y.max(), rtol=1e-10)
        np.testing.assert_allclose(
            float(row["peak_area"]),
            expected_area,
            rtol=1e-9,
            atol=1e-8,
            err_msg=f"{path.name} feature {row['feature_id']} {row['mz']} {row['gapfill_status']}",
        )
        assert int(row["scan_count"]) == np.count_nonzero(y)
        assert row["ms1_quant_status"] == ("gap_filled" if y.max() else "no_ms1_signal")
    return {
        "all_source_ms2_accounted_for": True,
        "unique_feature_per_assigned_ms2": True,
        "seeded_quantification_independently_verified": len(seeded),
        "seeded_numeric_measurements": sum(r["ms1_quant_status"] == "gap_filled" for r in seeded),
        "seeded_imputed_heights": sum(r["ms1_quant_status"] == "imputed_noise" for r in seeded),
        "known_isolation_width_above_4": known_wide,
        "missing_isolation_width": missing_width,
    }


def check_file(arguments):
    directory, row = arguments
    source = directory / row["name"]
    assert hashlib.sha256(source.read_bytes()).hexdigest() == row["sha256"]
    output = directory / f"v{__version__}-ms2" / (source.name + ".tsv")
    result, audit_rows = process_for_audit(source, output, overwrite=True)
    metadata = json.loads(Path(str(output) + ".json").read_text())
    original = reference_check(source, output, audit_rows)
    seeded = audit_seeded_quantification(source, output, Config(), audit_rows)
    return {
        "name": source.name,
        "input_sha256": row["sha256"],
        "group": row["group"],
        "features": result.features,
        "extraction_seconds": result.elapsed_seconds,
        "parameter_tuning": metadata["parameter_tuning"],
        "ion_relations": {
            k: metadata[k] for k in ("adduct_relations", "neutral_loss_relations", "isf_relations")
        },
        **{key: value for key, value in metadata.items() if key.startswith("ms2_")},
        "tsv_sha256": metadata["tsv_sha256"],
        "reference": {**original, **seeded},
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--directory", type=Path, default=Path("validation/MSV000086989"))
    parser.add_argument("--workers", type=int, default=4)
    parser.add_argument("--limit", type=int)
    parser.add_argument(
        "--output",
        type=Path,
        default=Path(f"benchmarks/results/MSV000086989-ms2-v{__version__}.json"),
    )
    args = parser.parse_args()
    manifest = json.loads((args.directory / "download_manifest.json").read_text())[: args.limit]
    started = time.perf_counter()
    rows = []
    with ProcessPoolExecutor(max_workers=args.workers) as pool:
        futures = [pool.submit(check_file, (args.directory, row)) for row in manifest]
        for future in as_completed(futures):
            rows.append(future.result())
            if len(rows) % 10 == 0 or len(rows) == len(manifest):
                print(json.dumps({"verified_files": len(rows), "files": len(manifest)}), flush=True)
    rows.sort(key=lambda row: row["name"])
    from masscubex.pipeline import implementation_id

    report = {
        "dataset": "MSV000086989",
        "version": __version__,
        "config": Config().to_dict(),
        "implementation_sha256": implementation_id()[0],
        "dependencies": implementation_id()[1],
        "python": platform.python_version(),
        "platform": platform.platform(),
        "files": len(rows),
        "wall_seconds_including_reference_checks": time.perf_counter() - started,
        "totals": {
            key: sum(row[key] for row in rows)
            for key in (
                "features",
                "ms2_scans",
                "ms2_matched_detected",
                "ms2_seeded_features",
                "ms2_seeded_scans",
                "ms2_assigned",
                "ms2_unassigned",
                "ms2_gapfill_sample_points",
            )
        },
        "initial_unassigned_reasons": dict(
            sum((Counter(r["ms2_initial_unassigned_reasons"]) for r in rows), Counter())
        ),
        "gapfill_status_counts": dict(
            sum((Counter(r["ms2_gapfill_status_counts"]) for r in rows), Counter())
        ),
        "rows": rows,
        "limitations": [
            "Raw-data interval/quantification audits use captured unrounded writer inputs; compact table rounding is checked separately.",
            "One centroid-mzXML dataset; this is bookkeeping/quantification validation, not proof of biological precursor identity.",
            "Profile mzML/mzXML, missing metadata and broad isolation windows are covered by synthetic tests.",
            "No claim of DIA classification or DIA fragment deconvolution.",
        ],
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(report, indent=2, allow_nan=False) + "\n")
    print(json.dumps(report["totals"], sort_keys=True))


if __name__ == "__main__":
    main()
