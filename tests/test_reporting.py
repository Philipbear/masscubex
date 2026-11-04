from concurrent.futures import Future
from concurrent.futures.process import BrokenProcessPool
from pathlib import Path
import csv
import gzip
import json
import multiprocessing
import os
import subprocess
import sys
import time

import pytest

from masscubex import Config
from masscubex import cli
from masscubex.reporting import BatchReport, InputFailure, error_result, retry_sources
from conftest import mzml_spectrum, mzxml_scan, write_mzml, write_mzxml


def _crash_test_job(source, output_dir, config, resume, overwrite, output_format="parquet"):
    # Exercise an actual lost child, without touching raw inputs or table files.
    assert multiprocessing.current_process().name != "MainProcess"
    if Path(source).name == "crash.mzML":
        os._exit(19)
    # Hold the other initially submitted jobs until the parent observes the
    # broken pool. Otherwise a faster-starting worker can consume the whole
    # input stream before the crashing worker starts, leaving no undispatched
    # inputs with which to exercise pool replacement.
    deadline = time.monotonic() + 30
    while not (Path(output_dir) / "crash_reported").exists():
        if time.monotonic() > deadline:
            raise TimeoutError("Parent did not report the injected worker crash")
        time.sleep(0.01)
    return {"input": str(source), "status": "ok", "features": 0}


def read_errors(path):
    with path.open(encoding="utf-8", newline="") as stream:
        return list(csv.DictReader(stream, delimiter="\t"))


@pytest.mark.parametrize("workers,executor", [(1, "process"), (2, "process"), (2, "thread")])
def test_corrupt_files_are_reported_and_later_files_complete_then_retry_only_failures(
    sample_file, tmp_path, workers, executor
):
    malformed = tmp_path / "truncated.mzML"
    malformed.write_bytes(sample_file.read_bytes()[:-20])
    binary = write_mzxml(tmp_path / "binary.mzXML", [mzxml_scan(1, 1, [100], [10])])
    binary.write_text(binary.read_text().replace("</peaks>", "%%%BAD%%%</peaks>"))
    compressed = tmp_path / "truncated.mzML.gz"
    compressed.write_bytes(gzip.compress(sample_file.read_bytes())[:-10])
    nonfinite = write_mzml(
        tmp_path / "nonfinite.mzML", [mzml_spectrum(1, 1, [100], [float("nan")])]
    )
    missing = tmp_path / "missing.mzML"
    good_xml = write_mzxml(
        tmp_path / "good.mzXML",
        [
            mzxml_scan(i, i * 2, [100], [y])
            for i, y in enumerate([10, 20, 40, 80, 100, 80, 40, 20, 10])
        ],
    )
    last = tmp_path / "last.mzML"
    last.write_bytes(sample_file.read_bytes())
    sources = [sample_file, malformed, binary, compressed, nonfinite, missing, good_xml, last]
    manifest = tmp_path / "inputs.txt"
    manifest.write_text("\n".join(str(path) for path in sources) + "\n")
    output = tmp_path / "features"
    report_root = tmp_path / "batch_reports"
    common = [
        sys.executable,
        "-m",
        "masscubex",
        "extract",
        "-o",
        str(output),
        "--workers",
        str(workers),
        "--executor",
        executor,
        "--report-dir",
        str(report_root),
    ]
    run = subprocess.run(common + ["--manifest", str(manifest)], capture_output=True, text=True)
    assert run.returncode == 1, run.stderr  # Failure is signalled after all files finish.
    results = [json.loads(line) for line in run.stdout.splitlines()]
    assert len(results) == 8
    assert {Path(r["input"]) for r in results if r["status"] == "ok"} == {
        sample_file,
        good_xml,
        last,
    }
    failures = {malformed, binary, compressed, nonfinite, missing}
    assert {Path(r["input"]) for r in results if r["status"] == "error"} == failures
    for source in failures:
        assert not cli.output_for(source, output).exists()
    assert not list(output.rglob("*.lock"))
    assert not list(output.rglob("*.partial"))
    (summary_path,) = report_root.glob("*/summary.json")
    summary = json.loads(summary_path.read_text())
    assert summary["status"] == "completed_with_errors"
    assert summary["counts"] == {"total": 8, "succeeded": 3, "skipped": 0, "failed": 5}
    assert summary["fatal_error"] is None
    assert summary["features_written"] == 3
    errors = summary_path.parent / "errors.tsv"
    error_rows = read_errors(errors)
    assert {Path(row["input"]) for row in error_rows} == failures
    assert all(
        row["error_type"] and row["error_message"] and row["stage"] == "processing"
        for row in error_rows
    )
    assert [
        json.loads(line) for line in (summary_path.parent / "files.jsonl").read_text().splitlines()
    ] == results
    assert str(summary_path) in run.stderr
    before = errors.read_bytes()
    successful_times = {
        Path(r["output"]): Path(r["output"]).stat().st_mtime_ns
        for r in results
        if r["status"] == "ok"
    }
    for source in (malformed, nonfinite, missing):
        source.write_bytes(sample_file.read_bytes())
    binary.write_bytes(good_xml.read_bytes())
    compressed.write_bytes(gzip.compress(sample_file.read_bytes()))
    run = subprocess.run(
        common + ["--retry-errors", str(errors), "--resume"], capture_output=True, text=True
    )
    assert run.returncode == 0, run.stderr
    retried = [json.loads(line) for line in run.stdout.splitlines()]
    assert len(retried) == 5 and all(r["status"] == "ok" for r in retried)
    assert {Path(r["input"]) for r in retried} == failures
    assert all(path.stat().st_mtime_ns == before for path, before in successful_times.items())
    assert errors.read_bytes() == before  # Previous reports are never overwritten.
    new_summary = next(path for path in report_root.glob("*/summary.json") if path != summary_path)
    assert json.loads(new_summary.read_text())["counts"] == {
        "total": 5,
        "succeeded": 5,
        "skipped": 0,
        "failed": 0,
    }
    assert read_errors(new_summary.parent / "errors.tsv") == []


def test_report_checkpoints_and_preserves_counts_on_interrupt(tmp_path):
    with pytest.raises(KeyboardInterrupt):
        with BatchReport(tmp_path, Config(), {}, checkpoint_seconds=0) as report:
            report.record(
                error_result(tmp_path / "bad\tline\nfile.mzML", ValueError("bad\tvalue\nnext"))
            )
            running = json.loads(report.summary_path.read_text())
            assert running["status"] == "running" and running["counts"]["failed"] == 1
            assert running["finished_at_utc"] is None
            raise KeyboardInterrupt()
    final = json.loads(report.summary_path.read_text())
    assert final["status"] == "interrupted" and final["counts"]["failed"] == 1
    assert final["fatal_error"]["type"] == "KeyboardInterrupt"
    assert final["finished_at_utc"] is not None
    rows = read_errors(report.directory / "errors.tsv")
    assert rows[0]["error_message"] == "bad\tvalue\nnext"
    assert list(retry_sources(report.directory / "errors.tsv")) == [
        (tmp_path / "bad\tline\nfile.mzML", "file")
    ]


def test_discovery_failure_does_not_hide_other_directories(sample_file, tmp_path, monkeypatch):
    blocked = tmp_path / "blocked"
    blocked.mkdir()
    (blocked / "hidden.mzML").write_bytes(sample_file.read_bytes())
    scandir = os.scandir

    def fail_one(path):
        if Path(path) == blocked:
            raise PermissionError("fixture: inaccessible directory")
        return scandir(path)

    monkeypatch.setattr(cli.os, "scandir", fail_one)
    inputs = list(cli.iter_inputs([tmp_path]))
    assert sample_file in inputs
    (failure,) = [p for p in inputs if isinstance(p, InputFailure)]
    assert failure.result["input"] == str(blocked)
    assert failure.result["input_kind"] == "directory"
    results = list(cli.run_batch(iter(inputs), tmp_path / "out", Config()))
    assert {r["status"] for r in results} == {"ok", "error"}


def test_manifest_discovery_failure_is_reported_and_can_be_retried(sample_file, tmp_path, capsys):
    manifest = tmp_path / "not_downloaded.txt"
    output = tmp_path / "out"
    assert cli.main(["extract", "--manifest", str(manifest), "-o", str(output)]) == 1
    capsys.readouterr()
    (errors,) = output.glob("reports/*/errors.tsv")
    assert read_errors(errors)[0]["input_kind"] == "manifest"
    manifest.write_text(sample_file.name + "\n")
    assert cli.main(["extract", "--retry-errors", str(errors), "-o", str(output)]) == 0
    assert json.loads(capsys.readouterr().out)["status"] == "ok"


def test_invalid_retry_report_and_empty_directory_leave_aborted_summary(tmp_path, capsys):
    invalid = tmp_path / "bad_report.tsv"
    invalid.write_text("wrong\tcolumns\n")
    output = tmp_path / "out"
    assert cli.main(["extract", "--retry-errors", str(invalid), "-o", str(output)]) == 1
    capsys.readouterr()
    (summary,) = output.glob("reports/*/summary.json")
    state = json.loads(summary.read_text())
    assert state["status"] == "aborted"
    assert state["fatal_error"]["type"] == "ValueError"
    assert state["counts"]["total"] == 0


def test_broken_worker_pool_is_replaced_without_losing_undispatched_inputs(
    sample_file, tmp_path, monkeypatch
):
    instances = []

    class Pool:
        def __init__(self, **kwargs):
            self.first = not instances
            self.submissions = 0
            instances.append(self)

        def submit(self, fn, *args):
            if fn is cli._worker_ready:
                future = Future()
                future.set_result(True)
                return future
            if self.first and self.submissions:
                raise BrokenProcessPool("fixture: worker exited")
            self.submissions += 1
            future = Future()
            if self.first:
                future.set_exception(BrokenProcessPool("fixture: uncertain in-flight result"))
            else:
                future.set_result(fn(*args))
            return future

        def shutdown(self, **kwargs):
            pass

    monkeypatch.setattr(cli, "ProcessPoolExecutor", Pool)
    files = []
    for i in range(7):
        path = tmp_path / f"{i}.mzML"
        path.write_bytes(sample_file.read_bytes())
        files.append(path)
    results = list(cli.run_batch(iter(files), tmp_path / "out", Config(), workers=2))
    assert len(results) == 7
    assert len(instances) == 2
    assert results[0]["status"] == "error" and results[0]["stage"] == "worker"
    assert all(r["status"] == "ok" for r in results[1:])
    assert {Path(r["input"]) for r in results} == set(files)


def test_failed_numerical_initializer_does_not_consume_or_mislabel_inputs(tmp_path, monkeypatch):
    consumed, closed = [], []

    class Pool:
        def __init__(self, **kwargs):
            pass

        def submit(self, fn, *args):
            assert fn is cli._worker_ready
            future = Future()
            future.set_exception(BrokenProcessPool("fixture: numerical initializer failed"))
            return future

        def shutdown(self, **kwargs):
            closed.append(True)

    def paths():
        consumed.append(True)
        yield tmp_path / "input.mzML"

    monkeypatch.setattr(cli, "ProcessPoolExecutor", Pool)
    with pytest.raises(BrokenProcessPool, match="numerical initializer"):
        list(cli.run_batch(paths(), tmp_path / "out", Config(), workers=2))
    assert not consumed
    assert closed == [True]


@pytest.mark.parametrize("count,expected", [(6, [2, 2, 2]), (7, [2, 2, 2, 1])])
def test_worker_recycling_drains_every_job_before_replacing_pool(tmp_path, monkeypatch, count, expected):
    pools = []

    class Pool:
        def __init__(self, **kwargs):
            self.jobs = []
            self.closed = False
            pools.append(self)

        def submit(self, fn, *args):
            assert not self.closed
            future = Future()
            if fn is cli._worker_ready:
                future.set_result(True)
            else:
                self.jobs.append(args[0])
                future.set_result({"input": str(args[0]), "status": "ok"})
            return future

        def shutdown(self, **kwargs):
            self.closed = True

    monkeypatch.setattr(cli, "ProcessPoolExecutor", Pool)
    paths = [tmp_path / f"{i}.mzML" for i in range(count)]
    rows = list(cli.run_batch(paths, tmp_path / "out", Config(), workers=2, worker_recycle_files=1))
    assert len(rows) == count
    assert {r["input"] for r in rows} == {str(p) for p in paths}
    assert [len(p.jobs) for p in pools] == expected
    assert all(p.closed for p in pools)


def test_resume_hits_do_not_trigger_worker_recycling(tmp_path, monkeypatch):
    pools = []

    class Pool:
        def __init__(self, **kwargs):
            pools.append(self)

        def submit(self, fn, *args):
            future = Future()
            future.set_result(True if fn is cli._worker_ready else {"status": "skipped"})
            return future

        def shutdown(self, **kwargs):
            pass

    monkeypatch.setattr(cli, "ProcessPoolExecutor", Pool)
    paths = [tmp_path / f"{i}.mzML" for i in range(20)]
    rows = list(cli.run_batch(paths, tmp_path / "out", Config(), workers=2, worker_recycle_files=1))
    assert len(rows) == 20 and len(pools) == 1


def test_abrupt_child_exit_does_not_prevent_remaining_files(tmp_path, monkeypatch):
    monkeypatch.setattr(cli, "_job", _crash_test_job)
    paths = [tmp_path / "crash.mzML"] + [tmp_path / f"good_{i}.mzML" for i in range(10)]
    output = tmp_path / "out"
    output.mkdir()
    results = []
    for result in cli.run_batch(paths, output, Config(), workers=2):
        results.append(result)
        if result["status"] == "error":
            (output / "crash_reported").touch()
    assert len(results) == len(paths)
    assert {r["input"] for r in results} == {str(p) for p in paths}
    failed = [r for r in results if r["status"] == "error"]
    assert str(paths[0]) in {r["input"] for r in failed}
    assert all(r["stage"] == "worker" and r["error_type"] == "BrokenProcessPool" for r in failed)
    assert len(failed) <= 4
    # Up to four submitted tasks can be affected; undispatched inputs survive.
    assert all(r["status"] == "ok" for r in results if r["input"] in {str(p) for p in paths[-3:]})
