"""Compare every measured row/MS2 assignment against pinned upstream functions.

Run fixed-floor controls and adaptive-noise runs on the same six public files.
The reference uses independently decoded Pyteomics arrays and the unmodified
function-body snapshot in tests/upstream_masscube.py. Temporary tables are removed.
"""

from pathlib import Path
from types import MethodType, SimpleNamespace
from tempfile import TemporaryDirectory
import json
import sys

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "tests"))
import upstream_masscube as reference

from audit_upstream_features import mzml_records
from compare_current_feature_counts import DEFAULT_FILES
from masscubex import Config, iter_feature_rows, process_file
from masscubex.adducts import catalog
from masscubex.noise import estimate_noise
from masscubex.pipeline import implementation_id


def original_run(source, config, *, grouping=False):
    scans, ms1, ms2 = [], [], []
    for index, record in enumerate(mzml_records(source)):
        signals = np.column_stack((record["m/z array"], record["intensity array"])).astype(
            np.float32
        )
        time = float(record["scanList"]["scan"][0]["scan start time"])
        level = record["ms level"]
        precursor = record["precursorList"]["precursor"]
        mass = (
            precursor[0]["selectedIonList"]["selectedIon"][0]["selected ion m/z"]
            if precursor
            else None
        )
        if level == 1:
            threshold = estimate_noise(signals[:, 0], signals[:, 1], config)[2]
            keep = (signals[:, 1] > threshold) & (signals[:, 0] > 0) & (signals[:, 0] < 100000)
        else:
            cutoff = max(
                config.ms2_min_intensity,
                np.max(signals[:, 1], initial=0) * config.ms2_relative_intensity,
            )
            keep = (signals[:, 1] > cutoff) & (signals[:, 0] > 0) & (signals[:, 0] < 100000)
        signals = reference.centroid_signals(signals[keep])
        scans.append(
            SimpleNamespace(signals=signals, time=time, precursor_mz=mass, raw_file_id=index)
        )
        (ms1 if level == 1 else ms2).append(index)
    data = SimpleNamespace(
        scans=scans,
        ms1_idx_arr=np.array(ms1),
        ms2_idx_arr=np.array(ms2),
        ms1_time_arr=np.array([scans[i].time for i in ms1]),
        metadata=SimpleNamespace(ion_mode="positive"),
        params=SimpleNamespace(
            mz_tol_ms1=0.01,
            feature_gap_tol=10,
            mz_tol_feature_grouping=0.01,
            rt_tol_feature_grouping=0.05,
            isotope_rel_int_limit=1.5,
            scan_scan_cor_tol=0.9,
        ),
    )
    features = reference.detect_features(data)
    for _ in range(2):
        features = [part for peak in features for part in reference.segment_feature(peak)]
    for feature in features:
        feature.finalize()
        reference._set_feature_peak_edges(data, feature)
    data.features = sorted(features, key=lambda f: f.mz)
    data.feature_mz_arr = np.array([f.mz for f in data.features])
    for index, feature in enumerate(data.features, 1):
        feature.id = index
    reference.allocate_ms2_to_features(data)
    for feature in data.features:
        feature.ms2 = reference.find_best_ms2(feature.ms2_seq)
    if grouping:
        reference.POS_ADDUCTS = {
            a.name: SimpleNamespace(
                mass_shift=a.mass_shift,
                charge=abs(a.charge),
                mol_multiplier=a.molecules,
                considered=a.id in config.adduct_ids,
            )
            for a in catalog().values()
            if a.polarity == "positive"
        }
        reference._EMPTY_F32 = np.empty(0, dtype=np.float32)
        reference._EMPTY_I32 = np.empty(0, dtype=np.int32)
        reference._EMPTY_SIG = np.empty((0, 2), dtype=np.float32)
        data.get_eic_data = MethodType(reference.get_eic_data, data)
        reference.group_features_single_file(data)
    return sorted(data.features, key=lambda f: f.id)


def compare(expected, actual, grouping):
    errors = []
    if len(expected) != len(actual):
        return [{"counts": [len(expected), len(actual)]}]
    for old, new in zip(expected, actual):
        values = dict(
            mz=round(float(old.mz), 5),
            rt_apex_seconds=round(old.rt * 60, 2),
            rt_start_seconds=round(old.rt_seq[0] * 60, 2),
            rt_end_seconds=round(old.rt_seq[-1] * 60, 2),
            peak_height=float(old.peak_height),
            peak_area=float(old.peak_area),
            scan_count=int(old.length),
            gaussian_similarity=round(float(old.gaussian_similarity), 4),
            ms2_spectrum_indices=[s.raw_file_id for s in old.ms2_seq],
        )
        if grouping:
            values.update(
                group_id=old.feature_group_id,
                adduct=old.adduct_type,
                is_isotope=old.is_isotope,
                isotope_state=old.isotope_state,
                is_in_source_fragment=old.is_in_source_fragment,
            )
        changed = {k: [v, new[k]] for k, v in values.items() if v != new[k]}
        if old.ms2 is not None:
            # The best scan number and corresponding source index must refer
            # to the same member of this row's complete provenance lists.
            pos = (
                new["ms2_spectrum_indices"].index(old.ms2.raw_file_id)
                if old.ms2.raw_file_id in new["ms2_spectrum_indices"]
                else None
            )
            if pos is None or new["ms2_scan_numbers"][pos] != new["best_ms2_scan_number"]:
                changed["best_ms2"] = [old.ms2.raw_file_id, new["best_ms2_scan_number"]]
        if changed:
            errors.append({"feature_id": new["feature_id"], "differences": changed})
    return errors


def main():
    root = Path("validation/MSV000086989")
    report_path = Path("benchmarks/results/upstream-restoration-2026-11-04.json")
    report = {"implementation_sha256": implementation_id()[0], "rows": [], "errors": []}
    with TemporaryDirectory(prefix="masscubex-upstream-parity-") as directory:
        for mode in ("fixed", "adaptive"):
            config = Config(
                auto_tune=False,
                profile_mode="reject",
                ms2_seeded_features=False,
                noise_mode="none" if mode == "fixed" else "adaptive",
                min_intensity=50000 if mode == "fixed" else 0,
                ms2_min_intensity=10000 if mode == "fixed" else 500,
            )
            for name in DEFAULT_FILES:
                grouping = name == "NIST_POS_Samp_14-03.mzXML"
                expected = original_run(root / name, config, grouping=grouping)
                output = Path(directory) / (mode + name + ".parquet")
                result = process_file(root / name, output, config=config)
                actual = list(iter_feature_rows(output))
                errors = compare(expected, actual, grouping)
                row = dict(
                    file=name,
                    mode=mode,
                    features=len(actual),
                    original_features=len(expected),
                    at_least_5_scans=sum(f["scan_count"] >= 5 for f in actual),
                    assigned_ms2=sum(f["ms2_count"] for f in actual),
                    grouping_compared=grouping,
                    mismatched_rows=len(errors),
                    first_mismatches=errors[:5],
                    seconds=result.elapsed_seconds,
                )
                report["rows"].append(row)
                print(json.dumps(row), flush=True)
                report_path.write_text(json.dumps(report, indent=2) + "\n")
                if errors:
                    raise AssertionError(f"Upstream parity failed: {mode} {name}")


if __name__ == "__main__":
    main()
