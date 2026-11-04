"""Count features using an untouched, pinned original MassCube checkout.

The pinned upstream no longer reads mzXML. A metadata-only adapter feeds
Pyteomics mzXML records into its existing mzML scan extractor; preprocessing,
tracking, two-pass segmentation and finalization are all upstream code.
Install upstream-only dependencies separately and put them on PYTHONPATH.
"""

import argparse
import hashlib
import json
from pathlib import Path
import sys
import time

import numpy as np
import pyarrow as pa
import pyarrow.parquet as pq
from pyteomics import mzxml


def mzml_records(source):
    """Translate scan metadata only; pass original arrays and minute RTs through."""
    with mzxml.MzXML(str(source), use_index=False) as reader:
        for spectrum in reader:
            rt = spectrum["retentionTime"]
            if str(rt.unit_info) != "minute":
                raise ValueError(f"Unexpected Pyteomics mzXML RT unit: {rt.unit_info}")
            record = {
                "ms level": spectrum["msLevel"],
                "scanList": {"scan": [{"scan start time": rt}]},
                "m/z array": spectrum["m/z array"],
                "intensity array": spectrum["intensity array"],
            }
            precursors = []
            for precursor in spectrum.get("precursorMz", []):
                entry = {
                    "selectedIonList": {
                        "selectedIon": [{"selected ion m/z": precursor["precursorMz"]}]
                    }
                }
                if "windowWideness" in precursor:
                    half = float(precursor["windowWideness"]) / 2
                    entry["isolationWindow"] = {
                        "isolation window lower offset": half,
                        "isolation window upper offset": half,
                    }
                precursors.append(entry)
            record["precursorList"] = {"precursor": precursors}
            yield record


def audit_file(source, output, instrument):
    from masscube.params import Params
    from masscube.raw_data_utils.core import MSData

    start = time.perf_counter()
    data = MSData()
    data.params = Params()
    data.params.set_default(instrument)
    data.update_metadata({"file_path": str(source), "ms_type": instrument, "ion_mode": "positive"})
    data.extract_scans_mzml(mzml_records(source))
    timings = {"read_preprocess": time.perf_counter() - start}
    retained = sum(len(data.scans[int(i)].signals) for i in data.ms1_idx_arr)
    started = time.perf_counter()
    data.detect_features()
    traces = len(data.features)
    long_traces = sum(int(f.length >= 5) for f in data.features)
    timings["detect"] = time.perf_counter() - started
    started = time.perf_counter()
    data.segment_features(iteration=2)
    timings["segment"] = time.perf_counter() - started
    started = time.perf_counter()
    data.finalize_features()
    timings["finalize"] = time.perf_counter() - started
    # This is exactly the row set used by upstream output_single_file().
    rows = {
        "mz": [float(f.mz) for f in data.features],
        "rt_apex_seconds": [float(f.rt) * 60 for f in data.features],
        "rt_start_seconds": [float(f.rt_seq[0]) * 60 for f in data.features],
        "rt_end_seconds": [float(f.rt_seq[-1]) * 60 for f in data.features],
        "peak_height": [float(f.peak_height) for f in data.features],
        "peak_area": [float(f.peak_area) for f in data.features],
        "scan_count": [int(f.length) for f in data.features],
        "ms2_count": [len(f.ms2_seq) for f in data.features],
    }
    pq.write_table(pa.table(rows), output, compression="zstd")
    result = {
        "file": source.name,
        "input_sha256": hashlib.sha256(source.read_bytes()).hexdigest(),
        "instrument": instrument,
        "ms1_absolute_threshold": data.params.ms1_abs_int_tol,
        "ms1_scans": len(data.ms1_idx_arr),
        "retained_ms1_points": retained,
        "unsegmented_traces": traces,
        "unsegmented_traces_at_least_5_scans": long_traces,
        "features": len(data.features),
        "features_at_least_5_scans": sum(int(f.length >= 5) for f in data.features),
        "features_with_ms2": sum(bool(f.ms2_seq) for f in data.features),
        "median_peak_height": float(np.median(rows["peak_height"])),
        "timings_seconds": timings,
        "elapsed_seconds": time.perf_counter() - start,
    }
    Path(str(output) + ".json").write_text(json.dumps(result, indent=2) + "\n")
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--upstream", type=Path, required=True)
    parser.add_argument("--input", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--instrument", choices=["orbitrap", "qtof"], default="orbitrap")
    parser.add_argument("--limit", type=int)
    parser.add_argument("--resume", action="store_true")
    args = parser.parse_args()
    sys.path.insert(0, str(args.upstream.resolve() / "src"))
    args.output.mkdir(parents=True, exist_ok=True)
    paths = sorted(args.input.glob("*.mzXML"))[: args.limit]
    results = []
    for path in paths:
        output = args.output / (path.name + ".parquet")
        metadata = Path(str(output) + ".json")
        if args.resume and output.exists() and metadata.exists():
            result = json.loads(metadata.read_text())
            if result["instrument"] != args.instrument:
                raise ValueError("Resume instrument mismatch")
        else:
            try:
                result = audit_file(path, output, args.instrument)
            except Exception as error:
                result = {"file": path.name, "error": f"{type(error).__name__}: {error}"}
        results.append(result)
        print(json.dumps(result), flush=True)
        (args.output / "summary.json").write_text(json.dumps(results, indent=2) + "\n")
    if any("error" in row for row in results):
        raise SystemExit(1)


if __name__ == "__main__":
    main()
