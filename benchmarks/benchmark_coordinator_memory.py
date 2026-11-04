"""Measure the coordinator's own RSS before/after process-pool initialization."""

import argparse
import json
from pathlib import Path
import sys
import tempfile


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--package-path", type=Path, required=True)
    parser.add_argument("--report", type=Path, required=True)
    args = parser.parse_args()
    sys.path.insert(0, str(args.package_path.absolute()))
    import psutil

    process = psutil.Process()
    before = process.memory_info().rss
    from masscubex import Config
    from masscubex.cli import run_batch

    imported = process.memory_info().rss
    with tempfile.TemporaryDirectory(prefix="masscubex-parent-probe-") as directory:
        list(run_batch(iter(()), Path(directory), Config(), workers=2))
    report = {
        "before_import_bytes": before, "after_cli_import_bytes": imported,
        "after_pool_startup_bytes": process.memory_info().rss,
        "numerical_modules_in_parent": [
            name for name in ["numpy", "numba", "scipy", "pyarrow"] if name in sys.modules
        ],
    }
    args.report.write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps(report))


if __name__ == "__main__":
    main()
