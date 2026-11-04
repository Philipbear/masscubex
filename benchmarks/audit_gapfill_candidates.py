"""Reconstruct the historical v0.9 3x rule using independent raw MS1 decoding."""

import json
from pathlib import Path

import numpy as np
import pyarrow.parquet as pq

from audit_eics import raw_ms1


def audit(name, root=Path("validation")):
    table = pq.read_table(
        root / "optimization-v0.8/final-parquet" / (name + ".parquet"),
        columns=[
            "feature_id",
            "seed_precursor_mz",
            "seed_rt_seconds",
            "rt_start_seconds",
            "rt_end_seconds",
            "ms1_quant_status",
        ],
    ).to_pandas()
    config = json.loads((root / "output-v0.9/parquet" / (name + ".parquet.json")).read_text())[
        "parameter_tuning"
    ]["effective_configs"]["positive"]
    times, masses, intensity, offsets = raw_ms1(
        root / "MSV000086989" / name, root / "eic-review/raw"
    )
    results = []

    def radius(mass):
        return max(
            config["gapfill_mz_tolerance_da"], mass * config["gapfill_mz_tolerance_ppm"] * 1e-6
        )

    for feature in table[table.ms1_quant_status == "ambiguous_mass"].itertuples():
        seed = feature.seed_precursor_mz
        scan_ids = np.flatnonzero(
            (times >= feature.rt_start_seconds) & (times <= feature.rt_end_seconds)
        )
        search = max(
            config["tuning_search_da"], seed * config["tuning_search_ppm"] * 1e-6, radius(seed)
        )
        retained_radius = search + 2 * radius(seed + search)
        ordered = []
        for i in scan_ids:
            mz = masses[offsets[i] : offsets[i + 1]]
            ordered.extend(mz[(mz >= seed - retained_radius) & (mz <= seed + retained_radius)])
        ordered = np.sort(ordered)
        centers, start = [], 0
        while start < len(ordered):
            end = np.searchsorted(
                ordered, ordered[start] + 2 * radius(ordered[start]), side="right"
            )
            centers.append(float(np.median(ordered[start:end])))
            start = end
        dt = np.median(np.diff(times[scan_ids])) if len(scan_ids) > 1 else 0
        local_ids = scan_ids[
            abs(times[scan_ids] - feature.seed_rt_seconds)
            <= max(config["ms2_rt_tolerance"], 2 * dt)
        ]
        candidates = []
        for center in centers:
            if abs(center - seed) > search:
                continue
            best_height, best_point = 0, None
            for i in local_ids:
                mz = masses[offsets[i] : offsets[i + 1]]
                a = offsets[i] + np.searchsorted(mz, center - radius(center))
                b = offsets[i] + np.searchsorted(mz, center + radius(center), side="right")
                if b > a:
                    point = a + int(np.argmax(intensity[a:b]))
                    if intensity[point] > best_height:
                        best_height, best_point = float(intensity[point]), int(point)
            if best_height:
                candidates.append(
                    dict(
                        center=center,
                        height=best_height,
                        point_index=best_point,
                        raw_centroid_mz=float(masses[best_point]),
                    )
                )
        candidates.sort(key=lambda c: (-c["height"], abs(c["center"] - seed)))
        precursor_radius = max(
            config["ms2_mz_tolerance_da"], seed * config["ms2_mz_tolerance_ppm"] * 1e-6
        )
        close = [c for c in candidates if abs(c["center"] - seed) <= precursor_radius]
        eligible = close or candidates
        if len(eligible) < 2 or not eligible[0]["height"] < 3 * eligible[1]["height"]:
            raise AssertionError(
                f"Historical ambiguity not reproduced: {name} #{feature.feature_id}"
            )
        results.append(
            dict(
                file=name,
                feature_id=feature.feature_id,
                seed_mz=seed,
                seed_rt_seconds=feature.seed_rt_seconds,
                ms1_radius_da=radius(seed),
                search_radius_da=search,
                top_two=eligible[:2],
                height_ratio=eligible[0]["height"] / eligible[1]["height"],
                same_raw_centroid=eligible[0]["point_index"] == eligible[1]["point_index"],
            )
        )
    return results


if __name__ == "__main__":
    manifest = json.loads(Path("validation/MSV000086989/download_manifest.json").read_text())
    results = []
    for row in manifest:
        if row["group"] in ("blank", "QC"):
            found = audit(row["name"])
            results.extend(found)
            print(row["name"], len(found), sum(r["same_raw_centroid"] for r in found), flush=True)
    output = dict(
        ambiguous_seeds=len(results),
        same_raw_centroid=sum(r["same_raw_centroid"] for r in results),
        scope="All six blanks and four QC files; archived v0.9 algorithm",
        examples=results,
    )
    Path("benchmarks/results/gapfill-ambiguity-v0.9.json").write_text(
        json.dumps(output, indent=2) + "\n"
    )
