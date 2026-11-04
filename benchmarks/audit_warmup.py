"""Find production kernel signatures missing from the explicit worker warmup."""

import argparse
import json
from pathlib import Path
import sys
import tempfile


def signatures():
    from numba.core.registry import CPUDispatcher

    result = {}
    for name, module in list(sys.modules.items()):
        if not name.startswith("masscubex.") or module is None:
            continue
        for value in vars(module).values():
            if isinstance(value, CPUDispatcher):
                function = value.py_func
                key = function.__module__ + "." + function.__name__
                result[key] = {str(signature) for signature in value.signatures}
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("input", type=Path)
    parser.add_argument("--report", type=Path, required=True)
    args = parser.parse_args()
    from masscubex import process_file
    from masscubex.features import warmup
    from masscubex.pipeline import implementation_id

    warmup()
    before = signatures()
    with tempfile.TemporaryDirectory(prefix="masscubex-warmup-audit-") as temporary:
        process_file(args.input, Path(temporary) / "features.parquet", overwrite=True)
    after = signatures()
    added = {key: sorted(value - before.get(key, set())) for key, value in after.items()}
    added = {key: value for key, value in added.items() if value}
    report = {"implementation_sha256": implementation_id()[0],
              "input": str(args.input.absolute()), "additional_signatures": added}
    args.report.write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
