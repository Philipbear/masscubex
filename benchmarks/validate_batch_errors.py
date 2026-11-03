"""Run the 50 public files with five deliberately bad inputs through the CLI.

Requires the downloaded corpus and the versioned validate_repository.py report.
Only validation/corruption_check fixtures are created; public inputs are immutable.
"""

from pathlib import Path
import argparse
import csv
import gzip
import hashlib
import json
import subprocess
import sys
import time

from masscubex import __version__


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--executor", choices=["process", "thread"], default="process")
    parser.add_argument("--format", choices=["parquet", "tsv"], default="parquet")
    parser.add_argument(
        "--reference-tables",
        type=Path,
        help="Current-version tables to compare byte-for-byte after a separate scientific audit",
    )
    parser.add_argument(
        "--reference",
        type=Path,
        default=Path(f"benchmarks/results/MSV000086989-ms2-v{__version__}.json"),
    )
    args = parser.parse_args()
    root = Path(__file__).resolve().parents[1]
    source_dir = root / "validation/MSV000086989"
    directory = (
        root / "validation/corruption_check" / f"v{__version__}-{args.executor}-{args.format}"
    )
    directory.mkdir(parents=True, exist_ok=True)
    reference = json.loads(args.reference.read_text())
    good = [source_dir / row["name"] for row in reference["rows"]]
    cases = {
        "truncated.mzML": b'<mzML xmlns="http://psi.hupo.org/ms/mzml"><run>',
        "truncated.mzXML": b'<mzXML><msRun><scan num="1" msLevel="1"',
        "bad_binary.mzXML": b'<mzXML><msRun><scan num="1" msLevel="1" peaksCount="1" centroided="1" polarity="+" retentionTime="PT1S"><peaks precision="64">!!!notbase64!!!</peaks></scan></msRun></mzXML>',
        "truncated.mzML.gz": gzip.compress(
            b'<mzML xmlns="http://psi.hupo.org/ms/mzml"><run><spectrumList/></run></mzML>'
        )[:-8],
    }
    bad = []
    for name, data in cases.items():
        path = directory / name
        path.write_bytes(data)
        bad.append(path)
    missing = directory / "deliberately_missing.mzML"
    assert not missing.exists(), (
        "Missing-file fixture unexpectedly exists; choose a clean validation directory"
    )
    bad.append(missing)
    paths = []
    for index, path in enumerate(good):
        if index % 10 == 0:
            paths.append(bad[index // 10])
        paths.append(path)
    manifest = directory / "files.txt"
    manifest.write_text("".join(str(path) + "\n" for path in paths))
    started = time.perf_counter()
    run = subprocess.run(
        [
            sys.executable,
            "-m",
            "masscubex",
            "extract",
            "--manifest",
            str(manifest),
            "-o",
            str(directory / "features"),
            "--workers",
            "4",
            "--executor",
            args.executor,
            "--format",
            args.format,
            "--overwrite",
        ],
        capture_output=True,
        text=True,
    )
    elapsed = time.perf_counter() - started
    assert run.returncode == 1, run.stderr
    outcomes = [json.loads(line) for line in run.stdout.splitlines()]
    assert len(outcomes) == 55
    assert {Path(row["input"]) for row in outcomes if row["status"] == "error"} == set(bad)
    assert {Path(row["input"]) for row in outcomes if row["status"] == "ok"} == set(good)
    expected = {row["name"]: row["tsv_sha256"] for row in reference["rows"]}
    for row in outcomes:
        if row["status"] == "ok":
            if args.reference_tables is not None:
                ref = args.reference_tables / (Path(row["input"]).name + "." + args.format)
                assert (
                    hashlib.sha256(Path(row["output"]).read_bytes()).digest()
                    == hashlib.sha256(ref.read_bytes()).digest()
                )
            elif args.format == "tsv":
                assert (
                    hashlib.sha256(Path(row["output"]).read_bytes()).hexdigest()
                    == expected[Path(row["input"]).name]
                )
            else:
                from masscubex import iter_feature_rows

                ref = (
                    source_dir / f"v{reference['version']}-ms2" / (Path(row["input"]).name + ".tsv")
                )
                assert list(iter_feature_rows(row["output"])) == list(iter_feature_rows(ref))
    line = next(line for line in run.stderr.splitlines() if line.startswith("Batch report: "))
    summary_path = Path(line.removeprefix("Batch report: "))
    summary = json.loads(summary_path.read_text())
    assert summary["counts"] == {"total": 55, "succeeded": 50, "skipped": 0, "failed": 5}
    assert summary["status"] == "completed_with_errors" and summary["fatal_error"] is None
    with (summary_path.parent / "errors.tsv").open(newline="") as stream:
        errors = list(csv.DictReader(stream, delimiter="\t"))
    assert {Path(row["input"]) for row in errors} == set(bad)
    assert not list((directory / "features").rglob("*.lock"))
    assert not list((directory / "features").rglob("*.partial"))
    resumed = subprocess.run(
        [
            sys.executable,
            "-m",
            "masscubex",
            "extract",
            "--manifest",
            str(manifest),
            "-o",
            str(directory / "features"),
            "--workers",
            "4",
            "--executor",
            args.executor,
            "--format",
            args.format,
            "--resume",
        ],
        capture_output=True,
        text=True,
    )
    assert resumed.returncode == 1, resumed.stderr
    resumed_outcomes = [json.loads(line) for line in resumed.stdout.splitlines()]
    assert len(resumed_outcomes) == 55
    assert sum(row["status"] == "skipped" for row in resumed_outcomes) == 50
    assert {Path(row["input"]) for row in resumed_outcomes if row["status"] == "error"} == set(bad)
    report = {
        "version": __version__,
        "counts": summary["counts"],
        "exit_code": run.returncode,
        "workers": 4,
        "executor": args.executor,
        "output_format": args.format,
        "reference_report": str(args.reference),
        "reference_tables": str(args.reference_tables) if args.reference_tables else None,
        "resume_skipped_all_50_completed_files": True,
        "wall_seconds_including_cli_startup": elapsed,
        "all_valid_tables_match_reference": True,
        "error_types": summary["error_types"],
        "error_files": [
            {"name": Path(row["input"]).name, "error_type": row["error_type"]} for row in errors
        ],
        "run_report": str(summary_path.relative_to(root)),
        "temporary_outputs_cleaned": True,
        "scope": "50 real centroid mzXML files plus five corrupt/missing fixtures. Error continuation/retry is also tested in sequential and thread modes by the test suite.",
    }
    destination = (
        root / f"benchmarks/results/batch-errors-{args.executor}-{args.format}-v{__version__}.json"
    )
    destination.write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
