"""Reproducible extraction-engine comparison with the cloned upstream commit.

Run: python benchmarks/benchmark_engine.py --scans 500 --peaks 500
Both engines see identical centroid arrays; parsing/noise/segmentation are
excluded. This is an engine benchmark, NOT an end-to-end speedup claim.
"""

from pathlib import Path
from types import ModuleType, SimpleNamespace
import argparse
import importlib.util
import json
import platform
import subprocess
import sys
import tempfile
import time

import numpy as np

from masscubex import Config
from masscubex.features import warmup
from masscubex.pipeline import extract_features
from masscubex.readers import Scan


UPSTREAM = "356748ede691862b9a87f67812db5e728a703b9d"


def load_upstream(directory):
    package = ModuleType("_masscubex_upstream")
    package.__path__ = [str(directory)]
    sys.modules[package.__name__] = package
    for name in ("feature_evaluation", "feature_detection"):
        source = subprocess.check_output(
            ["git", "show", f"{UPSTREAM}:src/masscube/{name}.py"],
            cwd=Path(__file__).resolve().parents[1],
        )
        path = directory / (name + ".py")
        path.write_bytes(source)
        spec = importlib.util.spec_from_file_location(package.__name__ + "." + name, path)
        module = importlib.util.module_from_spec(spec)
        sys.modules[spec.name] = module
        spec.loader.exec_module(module)
    return module.detect_features


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--scans", type=int, default=500)
    parser.add_argument("--peaks", type=int, default=500)
    parser.add_argument("--repeats", type=int, default=3)
    args = parser.parse_args()
    rng = np.random.default_rng(1729)
    masses = np.linspace(100, 1000, args.peaks)
    scans = []
    for i in range(args.scans):
        mz = masses + rng.normal(0, 0.0001, args.peaks)
        intensity = (
            100 + 1000 * np.exp(-0.5 * ((i - args.scans / 2) / (args.scans / 8)) ** 2)
        ) * np.ones(args.peaks)
        scans.append(Scan(float(i), mz, intensity))
    old = SimpleNamespace(
        scans=[
            SimpleNamespace(time=s.rt / 60, signals=np.column_stack((s.mz, s.intensity)))
            for s in scans
        ],
        ms1_idx_arr=np.arange(args.scans),
        params=SimpleNamespace(mz_tol_ms1=0.01, feature_gap_tol=10),
    )
    config = Config(
        noise_mode="none",
        segment=False,
        min_scans=2,
        max_gap_scans=10,
        mz_tolerance_da=0.01,
        mz_tolerance_ppm=0,
    )
    warmup()
    timings = {"upstream_seconds": [], "masscubex_seconds": []}
    with tempfile.TemporaryDirectory() as temporary:
        baseline = load_upstream(Path(temporary))
        for _ in range(args.repeats):
            start = time.perf_counter()
            old_features = baseline(old)
            timings["upstream_seconds"].append(time.perf_counter() - start)
            start = time.perf_counter()
            new_features = list(extract_features(iter(scans), config))
            timings["masscubex_seconds"].append(time.perf_counter() - start)
            assert len(old_features) == len(new_features) == args.peaks
    print(
        json.dumps(
            {
                "benchmark": "centroid trace detection and summary, no XML/noise/segmentation",
                "upstream_commit": UPSTREAM,
                "python": platform.python_version(),
                "platform": platform.platform(),
                "scans": args.scans,
                "peaks_per_scan": args.peaks,
                **timings,
                "median_speedup": float(
                    np.median(timings["upstream_seconds"]) / np.median(timings["masscubex_seconds"])
                ),
                "feature_counts_match": True,
            },
            indent=2,
        )
    )


if __name__ == "__main__":
    main()
