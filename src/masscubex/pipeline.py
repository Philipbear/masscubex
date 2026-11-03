"""Streaming extraction and crash-safe, independently resumable file outputs."""

from dataclasses import dataclass
from contextlib import closing
from datetime import datetime, timezone
from functools import lru_cache
from importlib.metadata import PackageNotFoundError, version
from pathlib import Path
import hashlib
import json
import os
import socket
import tempfile
import time

from . import __version__
from .config import Config
from .centroid import centroid_profile
from .features import TraceTracker, summarize_trace
from .output import FEATURE_COLUMNS, resolve_output, write_features
from .noise import estimate_noise
from .readers import file_format, iter_scans
from .ms2 import MS2Collector, associate_ms2
from .tuning import calibrate_file, config_for


SCHEMA_VERSION = 5
TSV_COLUMNS = FEATURE_COLUMNS  # Compatibility for table readers and existing audits.


@lru_cache(maxsize=1)
def implementation_id():
    """Invalidate resume after editable-code or numerical-dependency changes."""
    digest = hashlib.sha256()
    for name in (
        "config.py",
        "readers.py",
        "noise.py",
        "centroid.py",
        "features.py",
        "pipeline.py",
        "ms2.py",
        "gapfill.py",
        "gapfill_refinement.py",
        "tuning.py",
        "adducts.py",
        "grouping.py",
        "output.py",
    ):
        digest.update(Path(__file__).with_name(name).read_bytes())
    dependencies = {name: version(name) for name in ("numpy", "numba", "scipy", "lxml", "pyarrow")}
    # Optional decoding can affect rounding, so it belongs in provenance too.
    try:
        dependencies["pynumpress"] = version("pynumpress")
    except PackageNotFoundError:
        pass
    digest.update(json.dumps(dependencies, sort_keys=True).encode())
    return digest.hexdigest(), dependencies


@dataclass(frozen=True)
class ProcessingResult:
    input: str
    output: str
    status: str
    features: int
    ms1_scans: int
    elapsed_seconds: float


def extract_features(scans, config: Config, stats=None, configs=None):
    """Yield features from a scan iterator without accumulating completed traces."""
    stats = {} if stats is None else stats
    stats.update(ms1_scans=0, input_points=0, retained_centroids=0, profile_scans=0)
    trackers = {}
    for scan in scans:
        effective = config_for(scan.polarity, config, configs)
        tracker = trackers.get(scan.polarity)
        if tracker is None:
            tracker = trackers[scan.polarity] = TraceTracker(effective)
        stats["ms1_scans"] += 1
        stats["input_points"] += len(scan.mz)
        baseline, _, threshold = estimate_noise(scan.mz, scan.intensity, config)
        if scan.centroided:
            keep = scan.intensity > threshold
            mz, intensity, threshold = scan.mz[keep], scan.intensity[keep], threshold[keep]
        else:
            if config.profile_mode == "reject":
                raise ValueError(f"Profile MS1 spectrum {scan.native_id}; profile_mode is 'reject'")
            stats["profile_scans"] += 1
            mz, intensity, threshold = centroid_profile(
                scan.mz,
                scan.intensity,
                baseline,
                threshold,
                config.profile_min_points,
                float(config.profile_gap_factor),
            )
        stats["retained_centroids"] += len(mz)
        for trace in tracker.add(mz, intensity, threshold, scan.rt):
            yield from summarize_trace(trace, effective, scan.polarity)
    for polarity, tracker in trackers.items():
        for trace in tracker.flush():
            yield from summarize_trace(trace, tracker.config, polarity)
    stats["peak_active_traces"] = sum(t.peak_active_traces for t in trackers.values())
    stats["peak_live_nodes"] = sum(t.peak_live_nodes for t in trackers.values())
    stats["ms1_polarities"] = sorted(trackers)


def _signature(path):
    stat = path.stat()
    return {"path": str(path), "size": stat.st_size, "mtime_ns": stat.st_mtime_ns}


def _sha256(path):
    digest = hashlib.sha256()
    with open(path, "rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _timed_extraction(features, timings):
    iterator = iter(features)
    timings["feature_extraction"] = 0.0
    while True:
        started = time.perf_counter()
        try:
            feature = next(iterator)
        except StopIteration:
            timings["feature_extraction"] += time.perf_counter() - started
            return
        timings["feature_extraction"] += time.perf_counter() - started
        yield feature


def _resume_metadata(output, metadata_path, signature, config, output_format):
    try:
        metadata = json.loads(metadata_path.read_text())
        if (
            metadata["schema_version"] == SCHEMA_VERSION
            and metadata["version"] == __version__
            and metadata["implementation_sha256"] == implementation_id()[0]
            and metadata["input"] == signature
            and metadata["config"] == config.to_dict()
            and type(metadata["features"]) is int
            and metadata["features"] >= 0
            and type(metadata["ms1_scans"]) is int
            and metadata["ms1_scans"] >= 0
            and (
                metadata["ms1_scans"] > 0
                or (
                    config.associate_ms2 and config.ms2_seeded_features and metadata["features"] > 0
                )
            )
            and metadata["output_format"] == output_format
            and metadata["output_bytes"] == output.stat().st_size
            and metadata["output_sha256"] == _sha256(output)
        ):
            return metadata
    except (OSError, ValueError, KeyError, TypeError):
        pass
    return None


def process_file(
    input_path,
    output_path=None,
    *,
    config: Config | None = None,
    resume=False,
    overwrite=False,
    output_format=None,
):
    """Write one Parquet (default) or TSV table plus a JSON completion sidecar.

    An explicit output path selects its format by extension. ``output_format``
    may select the default name, or validate an explicitly supplied extension.

    ``resume=True`` skips only a matching, checksum-verified result and replaces
    stale results. No input-content hash is read: input identity is absolute path,
    byte size and nanosecond mtime. Keep inputs immutable during repository jobs.
    Each output has an exclusive lock. Worker exceptions propagate to callers.
    """
    start = time.perf_counter()
    config = config or Config()
    source = Path(input_path).resolve(strict=True)
    file_format(source)
    if not source.is_file():
        raise ValueError(f"Input is not a regular file: {source}")
    output, output_format = resolve_output(source, output_path, output_format)
    output.parent.mkdir(parents=True, exist_ok=True)
    metadata_path = Path(str(output) + ".json")
    lock = Path(str(output) + ".lock")
    try:
        fd = os.open(lock, os.O_CREAT | os.O_EXCL | os.O_WRONLY, 0o600)
    except FileExistsError as error:
        raise FileExistsError(
            f"Output is locked: {lock}. Remove a stale lock only after confirming its worker stopped."
        ) from error
    temporary = []
    try:
        with os.fdopen(fd, "w") as stream:
            json.dump(
                {
                    "pid": os.getpid(),
                    "host": socket.gethostname(),
                    "created": datetime.now(timezone.utc).isoformat(),
                },
                stream,
            )
        signature = _signature(source)
        if resume:
            metadata = _resume_metadata(output, metadata_path, signature, config, output_format)
            if metadata is not None:
                return ProcessingResult(
                    str(source),
                    str(output),
                    "skipped",
                    metadata["features"],
                    metadata["ms1_scans"],
                    time.perf_counter() - start,
                )
        if (output.exists() or metadata_path.exists()) and not (overwrite or resume):
            raise FileExistsError(f"Output exists: {output}; use resume=True or overwrite=True")
        started = time.perf_counter()
        tuning = calibrate_file(source, config)
        timings = {"calibration": time.perf_counter() - started}
        stats = {"parameter_tuning": tuning.diagnostics, "timings_seconds": timings}
        with tempfile.NamedTemporaryFile(
            mode="wb",
            dir=output.parent,
            prefix="." + output.name + ".",
            suffix=".partial",
            delete=False,
        ) as stream:
            temp_table = Path(stream.name)
            temporary.append(temp_table)
        with closing(
            iter_scans(
                source,
                assume_centroid=config.assume_centroid,
                include_ms2=config.associate_ms2,
                decode_ms2=config.group_ions and config.group_isf,
                sum_ms2=config.best_ms2 == "fragment_intensity",
            )
        ) as scans:
            collector = MS2Collector(config, tuning.configs)
            stream_ms1 = collector.ms1_stream(scans) if config.associate_ms2 else scans
            features = _timed_extraction(
                extract_features(stream_ms1, config, stats, tuning.configs), timings
            )
            if config.associate_ms2 or config.group_ions:
                # Compact summaries resolve overlaps globally; full MS1 spectra
                # remain streamed. The no-MS2/no-grouping path writes lazily.
                features = list(features)
                if config.associate_ms2:
                    from .gapfill import seed_ms2_features

                    started = time.perf_counter()
                    assignments = associate_ms2(
                        features, collector.records, config, stats, tuning.configs
                    )
                    timings["ms2_association"] = time.perf_counter() - started
                    started = time.perf_counter()
                    seed_ms2_features(
                        features,
                        collector.records,
                        assignments,
                        source,
                        config,
                        stats,
                        tuning.configs,
                    )
                    timings["ms2_gapfill"] = time.perf_counter() - started
                if config.group_ions:
                    from .grouping import group_ions

                    started = time.perf_counter()
                    group_ions(features, collector.records, config, stats, tuning.configs)
                    timings["ion_grouping"] = time.perf_counter() - started
            extraction_before_write = timings.get("feature_extraction", 0.0)
            started = time.perf_counter()
            count = write_features(temp_table, features, output_format)
            timings["table_write"] = (
                time.perf_counter()
                - started
                - timings["feature_extraction"]
                + extraction_before_write
            )
        if not stats["ms1_scans"] and not count:
            raise ValueError(f"No MS1 scans or usable MS2 precursor features found in {source}")
        if signature != _signature(source):
            raise RuntimeError(f"Input changed during extraction: {source}")
        elapsed = time.perf_counter() - start
        metadata = {
            "schema_version": SCHEMA_VERSION,
            "software": "masscubex",
            "version": __version__,
            "implementation_sha256": implementation_id()[0],
            "dependencies": implementation_id()[1],
            "input": signature,
            "config": config.to_dict(),
            "features": count,
            **stats,
            "elapsed_seconds": elapsed,
            "output_format": output_format,
            "output_bytes": temp_table.stat().st_size,
            "output_sha256": _sha256(temp_table),
        }
        if output_format == "tsv":
            metadata.update(
                tsv_bytes=metadata["output_bytes"], tsv_sha256=metadata["output_sha256"]
            )
        with tempfile.NamedTemporaryFile(
            mode="w",
            encoding="utf-8",
            dir=output.parent,
            prefix="." + output.name + ".",
            suffix=".partial",
            delete=False,
        ) as stream:
            temp_json = Path(stream.name)
            temporary.append(temp_json)
            json.dump(metadata, stream, indent=2, sort_keys=True, allow_nan=False)
            stream.write("\n")
            stream.flush()
            os.fsync(stream.fileno())
        os.replace(temp_table, output)
        # JSON is the completion marker. A crash between renames cannot pass resume validation.
        os.replace(temp_json, metadata_path)
        return ProcessingResult(
            str(source), str(output), "ok", count, stats["ms1_scans"], time.perf_counter() - start
        )
    finally:
        for path in temporary:
            path.unlink(missing_ok=True)
        lock.unlink(missing_ok=True)
