"""Time the actual MS1 estimators after JIT warmup, outside full extraction.

Synthetic sorted scans exercise sparse through dense grids and narrow/broad
windows. This is a kernel cost comparison, not a noise-accuracy benchmark.
"""

from pathlib import Path
from statistics import median
import argparse
import json
import platform
import time

import numpy as np

from masscubex import Config, __version__
from masscubex.noise import estimate_noise
from masscubex.pipeline import implementation_id


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--repeats", type=int, default=5)
    parser.add_argument(
        "--output", type=Path, default=Path(f"benchmarks/results/noise-v{__version__}.json")
    )
    args = parser.parse_args()
    if args.repeats < 1:
        parser.error("repeats must be positive")
    for name in ("low_tail", "clipped"):
        estimate_noise(np.linspace(100, 110, 20), np.ones(20), Config(noise_estimator=name))
    rng = np.random.default_rng(20261103)
    rows = []
    for n in (1000, 10000, 100000):
        mz = np.linspace(100, 1100, n)
        intensity = rng.lognormal(5, 0.7, n)
        intensity[::101] *= 1000
        for width in (25, 500):
            for name in ("low_tail", "clipped"):
                config = Config(noise_estimator=name, noise_window_da=width)
                times = []
                for _ in range(args.repeats):
                    start = time.perf_counter()
                    cutoff = estimate_noise(mz, intensity, config)[2]
                    times.append(time.perf_counter() - start)
                rows.append(
                    {
                        "points": n,
                        "window_da": width,
                        "estimator": name,
                        "median_seconds": median(times),
                        "seconds": times,
                        "retained_points": int((intensity > cutoff).sum()),
                    }
                )
    report = {
        "version": __version__,
        "python": platform.python_version(),
        "platform": platform.platform(),
        "implementation_sha256": implementation_id()[0],
        "seed": 20261103,
        "repeats": args.repeats,
        "results": rows,
        "limitations": [
            "Synthetic positive lognormal intensities with strong peaks every 101st point; not instrument noise truth.",
            "After JIT/import warmup; includes estimate_noise input checks, sort, temporary allocation and kernel.",
            "Excludes XML decoding, centroiding, feature extraction, MS2 handling and output.",
            "Different estimators retain different subsets; this is cost, not an equivalence/speedup comparison.",
        ],
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(report, indent=2) + "\n")
    for row in rows:
        print(
            f"{row['points']:>6} points / {row['window_da']:>3} Da / {row['estimator']:<8}: {1000 * row['median_seconds']:.3f} ms"
        )


if __name__ == "__main__":
    main()
