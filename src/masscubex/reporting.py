"""Bounded-memory, persistent batch outcomes and corruption reports."""

from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
import csv
import json
import os
import time
import uuid

from . import __version__


def error_result(source, error, *, stage="processing", input_kind="file"):
    return {
        "input": os.path.abspath(os.fspath(source)),
        "input_kind": input_kind,
        "status": "error",
        "stage": stage,
        "error_type": type(error).__name__,
        "error_message": str(error),
        "error": f"{type(error).__name__}: {error}",
    }


@dataclass(frozen=True)
class InputFailure:
    """A discovery failure that travels through the scheduler like an outcome."""

    result: dict


class BatchReport:
    """Stream outcomes to a unique run directory; atomically checkpoint counts.

    Report memory is constant in the number of inputs. Logs are flushed per
    outcome and fsynced at checkpoints/close, not fsynced for every successful
    file. An interruption cannot erase previous run reports.
    """

    ERROR_COLUMNS = ("input", "input_kind", "stage", "error_type", "error_message")
    SKIP_COLUMNS = ("input", "reason", "acquisition_mode", "vendor", "evidence")

    def __init__(self, root, config, options, *, checkpoint_seconds=30.0):
        self.started = datetime.now(timezone.utc)
        self.start_clock = time.monotonic()
        self.last_checkpoint = self.start_clock
        self.checkpoint_seconds = checkpoint_seconds
        self.run_id = self.started.strftime("%Y%m%dT%H%M%SZ") + "-" + uuid.uuid4().hex[:12]
        self.directory = Path(root).absolute() / self.run_id
        self.summary_path = self.directory / "summary.json"
        self.options = options
        self.config = config.to_dict()
        self.counts = {"total": 0, "succeeded": 0, "skipped": 0, "failed": 0}
        self.error_types = {}
        self.error_stages = {}
        self.skip_reasons = {}
        self.acquisition_modes = {}
        self.features_written = 0
        self.records = self.errors = self.skips = None

    def __enter__(self):
        self.directory.mkdir(parents=True, exist_ok=False)
        try:
            self.records = (self.directory / "files.jsonl").open("x", encoding="utf-8")
            self.errors = (self.directory / "errors.tsv").open("x", encoding="utf-8", newline="")
            self.writer = csv.DictWriter(
                self.errors, self.ERROR_COLUMNS, delimiter="\t", lineterminator="\n"
            )
            self.writer.writeheader()
            self.skips = (self.directory / "skipped.tsv").open("x", encoding="utf-8", newline="")
            self.skip_writer = csv.DictWriter(
                self.skips, self.SKIP_COLUMNS, delimiter="\t", lineterminator="\n"
            )
            self.skip_writer.writeheader()
            self.checkpoint("running")
        except BaseException:
            self._close()
            raise
        return self

    def record(self, result):
        status = result["status"]
        key = {"ok": "succeeded", "skipped": "skipped", "error": "failed"}[status]
        self.records.write(json.dumps(result, allow_nan=False) + "\n")
        self.records.flush()
        if status == "error":
            self.writer.writerow({key: result.get(key, "") for key in self.ERROR_COLUMNS})
            self.errors.flush()
            error_type = result.get("error_type", "UnknownError")
            stage = result.get("stage", "processing")
            self.error_types[error_type] = self.error_types.get(error_type, 0) + 1
            self.error_stages[stage] = self.error_stages.get(stage, 0) + 1
        elif status == "ok":
            self.features_written += result["features"]
        elif status == "skipped":
            reason = result.get("skip_reason") or "already_complete"
            self.skip_reasons[reason] = self.skip_reasons.get(reason, 0) + 1
            # Resume hits already have a complete feature table and a JSONL
            # outcome. Keep this file focused on excluded inputs.
            if reason != "already_complete":
                acquisition = result.get("acquisition") or {}
                self.skip_writer.writerow(
                    {
                        "input": result["input"],
                        "reason": reason,
                        "acquisition_mode": acquisition.get("mode", "unknown"),
                        "vendor": acquisition.get("vendor", "unknown"),
                        "evidence": "; ".join(acquisition.get("evidence", [])),
                    }
                )
                self.skips.flush()
        acquisition = result.get("acquisition")
        if acquisition:
            mode = acquisition["mode"]
            self.acquisition_modes[mode] = self.acquisition_modes.get(mode, 0) + 1
        self.counts["total"] += 1
        self.counts[key] += 1
        if time.monotonic() - self.last_checkpoint >= self.checkpoint_seconds:
            self.checkpoint("running")

    def checkpoint(self, status, fatal_error=None):
        now = datetime.now(timezone.utc).isoformat()
        for stream in (self.records, self.errors, self.skips):
            stream.flush()
            os.fsync(stream.fileno())
        summary = {
            "report_schema_version": 2,
            "software": "masscubex",
            "version": __version__,
            "run_id": self.run_id,
            "status": status,
            "started_at_utc": self.started.isoformat(),
            "updated_at_utc": now,
            "finished_at_utc": None if status == "running" else now,
            "elapsed_seconds": time.monotonic() - self.start_clock,
            "counts": self.counts,
            "features_written": self.features_written,
            "error_types": self.error_types,
            "error_stages": self.error_stages,
            "skip_reasons": self.skip_reasons,
            "acquisition_modes": self.acquisition_modes,
            "fatal_error": fatal_error,
            "working_directory": os.getcwd(),
            "options": self.options,
            "config": self.config,
            "reports": {
                "all_outcomes": "files.jsonl",
                "errors": "errors.tsv",
                "skipped": "skipped.tsv",
            },
        }
        temporary = self.directory / "summary.json.partial"
        try:
            with temporary.open("w", encoding="utf-8") as stream:
                json.dump(summary, stream, indent=2, allow_nan=False)
                stream.write("\n")
                stream.flush()
                os.fsync(stream.fileno())
            os.replace(temporary, self.summary_path)
        finally:
            temporary.unlink(missing_ok=True)
        self.last_checkpoint = time.monotonic()

    def _close(self):
        for stream in (self.records, self.errors, self.skips):
            if stream is not None:
                stream.close()

    def __exit__(self, error_type, error, traceback):
        try:
            if error_type is None:
                status = "completed_with_errors" if self.counts["failed"] else "completed"
                self.checkpoint(status)
            else:
                status = "interrupted" if isinstance(error, KeyboardInterrupt) else "aborted"
                self.checkpoint(status, {"type": error_type.__name__, "message": str(error)})
        finally:
            self._close()


def retry_sources(path):
    """Read an error TSV lazily, including quoted tabs/newlines in filenames."""
    path = Path(path).absolute()
    with path.open(encoding="utf-8", newline="") as stream:
        rows = csv.DictReader(stream, delimiter="\t", strict=True)
        if not rows.fieldnames or not {"input", "input_kind"} <= set(rows.fieldnames):
            raise ValueError(f"Not a MassCubeX error report: {path}")
        for row in rows:
            if not row.get("input"):
                raise ValueError(f"Missing input path in error report: {path}, row {rows.line_num}")
            source = Path(row["input"])
            yield (source if source.is_absolute() else path.parent / source), row["input_kind"]
