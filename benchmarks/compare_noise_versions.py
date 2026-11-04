"""Compare every MS1 noise estimate with an archived implementation.

The baseline argument is an archived src/masscubex/noise.py using the same
Config API. Both engines receive the same decoded arrays and default settings.
This is a correctness audit, not an end-to-end throughput benchmark.
"""

import argparse
import hashlib
import importlib.util
import json
from pathlib import Path
import sys
from tempfile import TemporaryDirectory

import numpy as np

from masscubex import Config, __version__
from masscubex.noise import estimate_noise
from masscubex.pipeline import implementation_id
from masscubex.readers import iter_scans


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--baseline", type=Path, required=True)
    parser.add_argument("--directory", type=Path, default=Path("validation/MSV000086989"))
    parser.add_argument("--report", type=Path, required=True)
    args = parser.parse_args()
    # Give the aliased module its own cache location. Otherwise its compiled
    # cache can break later canonical imports from the archived baseline.
    with TemporaryDirectory(prefix="masscubex-noise-audit-") as directory:
        source = Path(directory) / "noise.py"
        source.write_bytes(args.baseline.read_bytes())
        spec = importlib.util.spec_from_file_location("masscubex._baseline_noise", source)
        baseline = importlib.util.module_from_spec(spec)
        sys.modules[spec.name] = baseline
        spec.loader.exec_module(baseline)
        try:
            compare(args, baseline)
        finally:
            sys.modules.pop(spec.name, None)


def compare(args, baseline):
    config = Config()
    files = json.loads((args.directory / "download_manifest.json").read_text())
    rows = []
    for row in files:
        scans = points = differences = changed_decisions = 0
        for scan in iter_scans(args.directory / row["name"], include_ms2=False):
            expected = baseline.estimate_noise(scan.mz, scan.intensity, config)
            actual = estimate_noise(scan.mz, scan.intensity, config)
            scans += 1
            points += len(scan.mz)
            differences += sum(np.count_nonzero(a != b) for a, b in zip(actual, expected))
            changed_decisions += np.count_nonzero(
                (scan.intensity > actual[2]) != (scan.intensity > expected[2])
            )
        rows.append(
            dict(
                name=row["name"],
                scans=scans,
                points=points,
                differing_values=int(differences),
                changed_decisions=int(changed_decisions),
            )
        )
    report = {
        "version": __version__,
        "implementation_sha256": implementation_id()[0],
        "baseline_source_sha256": hashlib.sha256(args.baseline.read_bytes()).hexdigest(),
        "config": config.to_dict(),
        "files": len(rows),
        "ms1_scans": sum(row["scans"] for row in rows),
        "ms1_points": sum(row["points"] for row in rows),
        "differing_values": sum(row["differing_values"] for row in rows),
        "changed_decisions": sum(row["changed_decisions"] for row in rows),
        "rows": rows,
        "notes": [
            "Exact comparison of mean, population standard deviation and threshold.",
            "All MS1 points, including points removed before feature detection.",
        ],
    }
    args.report.parent.mkdir(parents=True, exist_ok=True)
    args.report.write_text(json.dumps(report, indent=2) + "\n")
    print(
        json.dumps({key: value for key, value in report.items() if key not in {"rows", "config"}})
    )
    if report["differing_values"]:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
