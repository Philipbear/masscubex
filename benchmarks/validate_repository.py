"""Validate and time the downloaded MSV000086989 sample, excluding downloads.

Run after download_massive.py. Optional --reference-check requires Pyteomics.
Raw data and TSVs stay in validation/; a compact public report is exported.
"""

from dataclasses import replace
from pathlib import Path
import argparse
import csv
import hashlib
import json
import platform
import resource
import sys
import time
from unittest.mock import patch

import numpy as np

from masscubex import Config, __version__, process_file
from masscubex.cli import run_batch
from masscubex.features import warmup
from masscubex.readers import MS2Scan, iter_scans
from masscubex.output import FEATURE_COLUMNS, FLOAT_COLUMNS, NESTED_COLUMNS, write_features


def reference_mass_radius(mass, settings, prefix=""):
    """Independent lookup using serialized bounds, without runtime lookup helpers."""
    regions = settings.get("mass_regions", [])
    hits = [r for r in regions if r["lower"] <= mass < r["upper"]]
    if regions and mass == regions[-1]["upper"]:
        hits = [regions[-1]]
    assert len(hits) <= 1
    selected = hits[0] if hits else settings
    return max(
        selected[prefix + "mz_tolerance_da"], mass * selected[prefix + "mz_tolerance_ppm"] * 1e-6
    )


def process_for_audit(path, output, **options):
    """Keep unrounded rows in this audit process, never in production outputs.

    Rounded coordinates cannot reproduce exact EIC interval membership. Capture
    the writer's input so independent raw-data audits keep their original tight
    tolerances. This wrapper is for serial/process workers, not shared threads.
    """
    rows = []

    def capture(destination, features, output_format):
        features = list(features)
        for index, feature in enumerate(features, 1):
            row = {"feature_id": str(index)}
            for name in FEATURE_COLUMNS[1:]:
                value = getattr(feature, name)
                row[name] = (
                    ""
                    if value is None
                    else format(value, ".17g")
                    if name in FLOAT_COLUMNS
                    else json.dumps(value)
                    if name in NESTED_COLUMNS
                    else str(value)
                )
            rows.append(row)
        return write_features(destination, features, output_format)

    with patch("masscubex.pipeline.write_features", capture):
        result = process_file(path, output, **options)
    return result, rows


def reference_check(path, output, full_precision_rows=None):
    """Check MS1 arrays and every MS2 ID/score against independent Pyteomics decoding.

    Audit all TSV associations and representative choices against those reference
    spectra. This verifies bookkeeping, not biological identity or DDA purity.
    """
    from pyteomics import mzxml

    ours = iter_scans(path, sum_ms2=True)
    checked = 0
    ms2 = {}
    with mzxml.MzXML(str(path), use_index=False) as reader:
        for index, record in enumerate(reader):
            if record["msLevel"] not in (1, 2):
                continue
            scan = next(ours)
            assert scan.spectrum_index == index
            assert scan.native_id == str(record["num"])
            np.testing.assert_allclose(scan.rt, float(record["retentionTime"]) * 60, atol=1e-10)
            if isinstance(scan, MS2Scan):
                score = float(np.sum(np.maximum(record["intensity array"], 0), dtype=np.float64))
                assert scan.fragment_intensity == score
                assert scan.scan_number == int(record["num"])
                precursors = record.get("precursorMz", [])
                assert scan.precursor_count == len(precursors)
                precursor = precursors[0] if len(precursors) == 1 else {}
                if precursor:
                    assert scan.precursor_mz == precursor["precursorMz"]
                    assert scan.precursor_intensity == precursor.get("precursorIntensity")
                    assert scan.precursor_ref == precursor.get("precursorScanNum")
                ms2[scan.native_id] = {
                    "score": score,
                    "rt": scan.rt,
                    "index": index,
                    "mz": precursor.get("precursorMz"),
                    "polarity": record.get("polarity"),
                }
            elif checked < 5:
                keep = record["intensity array"] > 0
                np.testing.assert_array_equal(scan.mz, record["m/z array"][keep])
                np.testing.assert_array_equal(scan.intensity, record["intensity array"][keep])
                checked += 1
    assert next(ours, None) is None
    ours.close()
    assigned = set()
    features_with_ms2 = 0
    metadata = json.loads(Path(str(output) + ".json").read_text())
    if metadata.get("schema_version", 0) >= 6 and full_precision_rows is None:
        raise ValueError("Rounded tables need unrounded audit coordinates; use process_for_audit")
    configured = metadata.get("config", Config().to_dict())
    resolved = metadata.get("parameter_tuning", {}).get("effective_configs", {})
    with output.open() as stream:
        for feature in (
            full_precision_rows
            if full_precision_rows is not None
            else csv.DictReader(stream, delimiter="\t")
        ):
            native_ids = json.loads(feature["ms2_native_ids"])
            numbers = json.loads(feature["ms2_scan_numbers"])
            indices = json.loads(feature["ms2_spectrum_indices"])
            assert len(native_ids) == len(numbers) == len(indices) == int(feature["ms2_count"])
            if not native_ids:
                assert not feature["best_ms2_native_id"]
                continue
            features_with_ms2 += 1
            apex = float(feature["rt_apex_seconds"])
            effective = resolved.get(feature["polarity"], configured)
            for native_id, number, index in zip(native_ids, numbers, indices):
                assert native_id not in assigned
                assigned.add(native_id)
                record = ms2[native_id]
                assert number == int(native_id) and index == record["index"]
                tolerance = reference_mass_radius(record["mz"], effective, "ms2_")
                if feature.get("seed_precursor_mz"):
                    # Refined seeds can recover an ion outside the narrow matching
                    # tolerance, but never outside the initial broad search.
                    if (
                        "noise_layout" not in effective
                        and feature.get("gapfill_status") != "fixed_window"
                    ):
                        tolerance = max(
                            tolerance,
                            effective["tuning_search_da"],
                            record["mz"] * effective["tuning_search_ppm"] * 1e-6,
                        )
                    assert abs(record["mz"] - float(feature["mz"])) <= tolerance + 1e-8
                    if feature.get("gapfill_status") == "peak_bounds":
                        assert (
                            float(feature["rt_start_seconds"])
                            - effective["ms2_rt_tolerance"]
                            - 1e-8
                            <= record["rt"]
                            <= float(feature["rt_end_seconds"])
                            + effective["ms2_rt_tolerance"]
                            + 1e-8
                        )
                    else:
                        assert (
                            abs(record["rt"] - float(feature["seed_rt_seconds"]))
                            <= effective["ms2_seed_rt_tolerance"] + 1e-8
                        )
                else:
                    assert abs(record["mz"] - float(feature["mz"])) <= tolerance + 1e-8
                    assert (
                        float(feature["rt_start_seconds"]) - effective["ms2_rt_tolerance"] - 1e-8
                        <= record["rt"]
                        <= float(feature["rt_end_seconds"]) + effective["ms2_rt_tolerance"] + 1e-8
                    )
                assert feature["polarity"] == {"+": "positive", "-": "negative"}[record["polarity"]]
            expected = max(
                native_ids,
                key=lambda id: (ms2[id]["score"], -abs(ms2[id]["rt"] - apex), -ms2[id]["index"]),
            )
            assert feature["best_ms2_native_id"] == expected
            assert int(feature["best_ms2_scan_number"]) == int(expected)
            assert int(feature["best_ms2_spectrum_index"]) == ms2[expected]["index"]
            assert configured["best_ms2"] == "fragment_intensity"
            np.testing.assert_allclose(
                float(feature["best_ms2_fragment_intensity"]), ms2[expected]["score"], rtol=1e-10
            )
    return {
        "ms1_spectra": checked,
        "ms2_spectra": len(ms2),
        "assigned_ms2": len(assigned),
        "features_with_ms2": features_with_ms2,
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--directory", type=Path, default=Path("validation/MSV000086989"))
    parser.add_argument(
        "--output", type=Path, default=Path(f"benchmarks/results/MSV000086989-v{__version__}.json")
    )
    parser.add_argument("--reference-check", action="store_true")
    args = parser.parse_args()
    downloads = json.loads((args.directory / "download_manifest.json").read_text())
    paths = [args.directory / row["name"] for row in downloads]
    config = Config()
    result_dir = args.directory / f"v{__version__}"
    for row, path in zip(downloads, paths):
        if hashlib.sha256(path.read_bytes()).hexdigest() != row["sha256"]:
            raise ValueError(f"Input checksum mismatch: {path}")
    warmup()
    rows = []
    checked = []
    started = time.perf_counter()
    for download, path in zip(downloads, paths):
        output = result_dir / "features" / (path.name + ".tsv")
        if args.reference_check:
            result, audit_rows = process_for_audit(path, output, overwrite=True)
            checked.append(reference_check(path, output, audit_rows))
        else:
            result = process_file(path, output, overwrite=True)
        metadata = json.loads(Path(result.output + ".json").read_text())
        row = {
            **download,
            **{
                k: metadata[k]
                for k in (
                    "features",
                    "ms1_scans",
                    "profile_scans",
                    "input_points",
                    "retained_centroids",
                    "peak_active_traces",
                    "peak_live_nodes",
                    "tsv_sha256",
                    "ms2_scans",
                    "ms2_assigned",
                    "ms2_skipped_wide_isolation",
                    "ms2_skipped_multiprecursor",
                    "ms2_missing_metadata",
                )
            },
            "seconds": result.elapsed_seconds,
        }
        rows.append(row)
        print(
            json.dumps({"name": path.name, "features": row["features"], "seconds": row["seconds"]}),
            flush=True,
        )
    sequential_wall = time.perf_counter() - started
    parallel = {}
    expected = {r["name"]: r["tsv_sha256"] for r in rows}
    for executor in ("process", "thread"):
        started = time.perf_counter()
        results = list(
            run_batch(
                paths,
                result_dir / (executor + "_features"),
                config,
                workers=4,
                executor=executor,
                overwrite=True,
                output_format="tsv",
            )
        )
        elapsed = time.perf_counter() - started
        if not all(r["status"] == "ok" for r in results):
            raise RuntimeError(results)
        for result in results:
            checksum = hashlib.sha256(Path(result["output"]).read_bytes()).hexdigest()
            assert checksum == expected[Path(result["input"]).name]
        parallel[executor] = {"workers": 4, "wall_seconds": elapsed, "identical_tsvs": True}
    started = time.perf_counter()
    resumed = list(
        run_batch(
            paths,
            result_dir / "process_features",
            config,
            workers=1,
            resume=True,
            output_format="tsv",
        )
    )
    assert all(r["status"] == "skipped" for r in resumed)
    resume_seconds = time.perf_counter() - started
    comparisons = []
    for group in ("biological", "blank", "QC"):
        row = next(r for r in rows if r["group"] == group)
        for floor in (500, 30000):
            result = process_file(
                args.directory / row["name"],
                result_dir / "fixed_noise" / (f"{floor}_" + row["name"] + ".tsv"),
                config=replace(config, noise_mode="none", min_intensity=floor),
                overwrite=True,
            )
            metadata = json.loads(Path(result.output + ".json").read_text())
            comparisons.append(
                {
                    "name": row["name"],
                    "group": group,
                    "fixed_floor": floor,
                    "features": result.features,
                    "retained_centroids": metadata["retained_centroids"],
                    "adaptive_features": row["features"],
                    "adaptive_retained_centroids": row["retained_centroids"],
                }
            )
    total_points = sum(r["input_points"] for r in rows)
    total_retained = sum(r["retained_centroids"] for r in rows)
    report = {
        "dataset": "MSV000086989",
        "dataset_page": "https://massive.ucsd.edu/ProteoSAFe/dataset.jsp?task=2f73277b4e034948acebfdf1edab17ed",
        "version": __version__,
        "python": platform.python_version(),
        "platform": platform.platform(),
        "implementation_sha256": metadata["implementation_sha256"],
        "dependencies": metadata["dependencies"],
        "config": config.to_dict(),
        "files": len(rows),
        "bytes": sum(r["bytes"] for r in rows),
        "ms1_scans": sum(r["ms1_scans"] for r in rows),
        "profile_scans": sum(r["profile_scans"] for r in rows),
        "features": sum(r["features"] for r in rows),
        "ms2_scans": sum(r["ms2_scans"] for r in rows),
        "ms2_assigned": sum(r["ms2_assigned"] for r in rows),
        "input_points": total_points,
        "retained_centroids": total_retained,
        "removed_fraction": 1 - total_retained / total_points,
        "sequential_wall_seconds": sequential_wall,
        "median_file_seconds": float(np.median([r["seconds"] for r in rows])),
        "parent_peak_rss_mb": resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
        / (1e6 if sys.platform == "darwin" else 1e3),
        "parallel": parallel,
        "resume_all_skipped": True,
        "resume_seconds": resume_seconds,
        "reference_validation": {
            key: sum(r[key] for r in checked)
            for key in ("ms1_spectra", "ms2_spectra", "assigned_ms2", "features_with_ms2")
        },
        "fixed_threshold_comparisons": comparisons,
        "limitations": [
            "Selected public mzXML files are centroid exports; profile mass accuracy is tested separately with synthetic truth.",
            "Feature counts and filtering fractions are not identification accuracy or measured false discovery rates.",
            "Timing excludes downloads and initial parent-process JIT/imports; parallel timings include worker startup.",
            "Parent RSS includes Python, scientific libraries, the JIT runtime and validation checks; it excludes child processes.",
        ],
        "results": rows,
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(report, indent=2) + "\n")
    print(f"Saved {args.output}", flush=True)


if __name__ == "__main__":
    main()
