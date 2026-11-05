"""Streaming extraction and crash-safe, independently resumable file outputs."""

from dataclasses import dataclass, replace
from bisect import bisect_left
from contextlib import closing, ExitStack
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
from .acquisition import inspect_acquisition
from .config import Config
from .features import TraceTracker
from .output import (
    FEATURE_COLUMNS,
    OUTPUT_PRECISION,
    TABLE_SCHEMA_VERSION,
    feature_columns,
    resolve_output,
    write_features,
)
from .upstream_adapter import preprocess_ms1, resolve_ms2_filter
from .readers import file_format
from .scan_cache import ScanSource
from .ms2 import MS2Collector, associate_ms2
from .tuning import calibrate_file, config_for


SCHEMA_VERSION = TABLE_SCHEMA_VERSION
TSV_COLUMNS = FEATURE_COLUMNS  # Compatibility for table readers and existing audits.


@lru_cache(maxsize=1)
def implementation_id():
    """Invalidate resume after editable-code or numerical-dependency changes."""
    digest = hashlib.sha256()
    for name in (
        "config.py",
        "upstream.py",
        "upstream_adapter.py",
        "mass_tolerance.py",
        "mass_tolerance_kernels.py",
        "formats.py",
        "spectrum_types.py",
        "acquisition.py",
        "readers.py",
        "array_validation.py",
        "scan_cache.py",
        "noise.py",
        "noise_regions.py",
        "consistency_noise.py",
        "consistency_kernels.py",
        "consistency_cache.py",
        "profile_noise.py",
        "centroid.py",
        "features.py",
        "summary.py",
        "measurements.py",
        "chromatograms.py",
        "peak_shapes.py",
        "reference_chromatograms.py",
        "pipeline.py",
        "ms2.py",
        "gapfill.py",
        "gapfill_refinement.py",
        "tuning.py",
        "pilot_sampling.py",
        "calibration.py",
        "adducts.py",
        "grouping.py",
        "grouping_queries.py",
        "grouping_spectra.py",
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
    output: str | None
    status: str
    features: int
    ms1_scans: int
    elapsed_seconds: float
    skip_reason: str | None = None
    acquisition: dict | None = None


def extract_features(
    scans, config: Config, stats=None, configs=None, spectra=None, *, preprocessor=None
):
    """Yield features from a scan iterator without accumulating completed traces."""
    stats = {} if stats is None else stats
    preprocessor = preprocessor or preprocess_ms1
    stats.update(ms1_scans=0, input_points=0, retained_centroids=0, profile_scans=0)
    trackers, acquisition_times = {}, {}

    def with_edges(features, polarity):
        # MassCube extends MS2 association bounds to the adjacent acquired MS1
        # scans, independently of the trimmed trace's exported start/end.
        times = acquisition_times[polarity]
        for feature in features:
            # Upstream minutes -> exported seconds can change the last bit.
            # Locate the original acquired endpoints before adding neighbors.
            # Two scalar bisections avoid allocating several NumPy arrays per
            # feature. Preserve the strict nearest-time comparison and ties.
            left = _nearest_acquisition(times, feature.rt_start_seconds)
            right = _nearest_acquisition(times, feature.rt_end_seconds)
            feature._peak_edges = (
                times[max(0, left - 1)] / 60.0,
                times[min(len(times) - 1, right + 1)] / 60.0,
            )
            yield feature

    for scan in scans:
        if not 0 <= scan.rt <= 600000:
            continue
        effective = config_for(scan.polarity, config, configs)
        tracker = trackers.get(scan.polarity)
        if tracker is None:
            tracker = trackers[scan.polarity] = TraceTracker(effective)
            acquisition_times[scan.polarity] = []
        acquisition_times[scan.polarity].append(scan.rt)
        stats["ms1_scans"] += 1
        stats["input_points"] += len(scan.mz)
        if not scan.centroided:
            stats["profile_scans"] += 1
        mz, intensity, threshold = preprocessor(scan, config)
        stats["retained_centroids"] += len(mz)
        if spectra is not None:
            spectra.append(scan.rt, scan.polarity, mz, intensity)
        yield from with_edges(
            tracker.add_features(mz, intensity, threshold, scan.rt, scan.polarity), scan.polarity
        )
    for polarity, tracker in trackers.items():
        yield from with_edges(tracker.flush_features(polarity), polarity)
    stats["peak_active_traces"] = sum(t.peak_active_traces for t in trackers.values())
    stats["peak_live_nodes"] = sum(t.peak_live_nodes for t in trackers.values())
    stats["ms1_polarities"] = sorted(trackers)


def _nearest_acquisition(times, target):
    """Locate an original acquisition after the minutes/seconds roundtrip."""
    after = min(bisect_left(times, target), len(times) - 1)
    before = max(after - 1, 0)
    return before if abs(times[before] - target) < abs(times[after] - target) else after


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
    grouping_spectra = None
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
                    skip_reason="already_complete",
                    acquisition=metadata.get("acquisition"),
                )
        started = time.perf_counter()
        acquisition = inspect_acquisition(source)
        inspection_seconds = time.perf_counter() - started
        if signature != _signature(source):
            raise RuntimeError(f"Input changed during acquisition inspection: {source}")
        if acquisition.mode == "chromatogram_only":
            raise ValueError(
                f"Chromatogram-only export has no mass spectra: {source}; "
                "MRM/SRM chromatograms cannot be used for untargeted MS1 feature extraction"
            )
        if acquisition.mode == "non_ms_only":
            raise ValueError(f"Non-MS-only export contains no mass spectra: {source}")
        if acquisition.is_dia:
            if config.dia_policy == "skip":
                return ProcessingResult(
                    str(source),
                    None,
                    "skipped",
                    0,
                    0,
                    time.perf_counter() - start,
                    skip_reason="dia_acquisition",
                    acquisition=acquisition.to_dict(),
                )
            if config.dia_policy == "error":
                raise ValueError(f"DIA acquisition: {source}: {'; '.join(acquisition.evidence)}")
        if (output.exists() or metadata_path.exists()) and not (overwrite or resume):
            raise FileExistsError(f"Output exists: {output}; use resume=True or overwrite=True")
        started = time.perf_counter()
        # Broad-window fragments do not identify individual MS1 precursors.
        # Apply one conservative policy to the entire file, including mixed
        # DDA/DIA acquisitions, while retaining MS1 isotope/adduct grouping.
        acquisition_config = (
            replace(config, associate_ms2=False, ms2_seeded_features=False, group_isf=False)
            if acquisition.is_dia
            else config
        )
        runtime_config = resolve_ms2_filter(source, acquisition_config)
        scan_source = ScanSource(source, runtime_config)
        tuning = calibrate_file(source, runtime_config, scan_source=scan_source)
        timings = {
            "acquisition_inspection": inspection_seconds,
            "calibration": time.perf_counter() - started,
        }
        stats = {
            "acquisition": acquisition.to_dict(),
            "parameter_tuning": tuning.diagnostics,
            "timings_seconds": timings,
            "scan_cache": scan_source.stats,
            "ms2_processing_enabled": runtime_config.associate_ms2,
            "ms2_min_intensity_used": (
                runtime_config.ms2_min_intensity if runtime_config.associate_ms2 else None
            ),
        }
        with tempfile.NamedTemporaryFile(
            mode="wb",
            dir=output.parent,
            prefix="." + output.name + ".",
            suffix=".partial",
            delete=False,
        ) as stream:
            temp_table = Path(stream.name)
            temporary.append(temp_table)
        with closing(scan_source.scans()) as scans, ExitStack() as preprocessing:
            if runtime_config.group_ions:
                from .grouping_spectra import GroupingSpectra

                # Preprocessing already produces exact upstream float32 values.
                # Keep that representation on disk instead of doubling its size.
                grouping_spectra = GroupingSpectra(output.parent, dtype="float32")
            collector = MS2Collector(runtime_config, tuning.configs)
            stream_ms1 = collector.ms1_stream(scans) if runtime_config.associate_ms2 else scans
            noise_cache = None
            if (
                runtime_config.noise_mode == "adaptive"
                and runtime_config.noise_estimator == "consistency"
            ):
                from .consistency_cache import ConsistencyCache

                started = time.perf_counter()
                noise_cache = preprocessing.enter_context(ConsistencyCache(output.parent))
                noise_cache.build(
                    stream_ms1,
                    collector.noise_records if collector.noise_records is not None else [],
                    runtime_config,
                    tuning.configs,
                )
                timings["adaptive_noise"] = time.perf_counter() - started
                stream_ms1 = noise_cache.frames
            features = _timed_extraction(
                extract_features(
                    stream_ms1,
                    runtime_config,
                    stats,
                    tuning.configs,
                    grouping_spectra,
                    preprocessor=noise_cache.preprocess if noise_cache is not None else None,
                ),
                timings,
            )
            features = list(features)
            if noise_cache is not None:
                noise_cache.apply_stats(stats)
                preprocessing.close()
            # Upstream sorts parent traces by apex m/z before segmentation,
            # then sorts the resulting features by apex m/z again.
            features.sort(key=lambda f: (f.mz, f._parent_mz))
            if runtime_config.associate_ms2 or runtime_config.group_ions:
                # Compact summaries resolve overlaps globally; full MS1 spectra
                # remain streamed. The no-MS2/no-grouping path writes lazily.
                if runtime_config.associate_ms2:
                    from .gapfill import seed_ms2_features

                    started = time.perf_counter()
                    assignments = associate_ms2(
                        features, collector.records, runtime_config, stats, tuning.configs
                    )
                    timings["ms2_association"] = time.perf_counter() - started
                    started = time.perf_counter()
                    seed_ms2_features(
                        features,
                        collector.records,
                        assignments,
                        source,
                        runtime_config,
                        stats,
                        tuning.configs,
                        scan_source=scan_source,
                    )
                    timings["ms2_gapfill"] = time.perf_counter() - started
                scan_source.release()
                if runtime_config.group_ions:
                    from .grouping import group_ions

                    started = time.perf_counter()
                    group_ions(
                        features, collector.records, runtime_config, stats, spectra=grouping_spectra
                    )
                    timings["ion_grouping"] = time.perf_counter() - started
            extraction_before_write = timings.get("feature_extraction", 0.0)
            started = time.perf_counter()
            count = write_features(temp_table, features, output_format, config=runtime_config)
            timings["table_write"] = (
                time.perf_counter()
                - started
                - timings["feature_extraction"]
                + extraction_before_write
            )
        if acquisition.is_dia and not stats["ms1_scans"]:
            raise ValueError(f"DIA MS1-only extraction requires MS1 scans: {source}")
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
            "effective_config": runtime_config.to_dict(),
            "acquisition_policy": {
                "version": 1,
                "requested_dia_policy": config.dia_policy,
                "applied": "dia_ms1_only" if acquisition.is_dia else "standard",
                "scope": "entire_file",
                "disabled_operations": (
                    ["ms2_association", "ms2_seeded_features", "ms2_in_source_fragments"]
                    if acquisition.is_dia
                    else []
                ),
                "dia_deconvolution": False,
            },
            "algorithm_policy": {
                "upstream_commit": "356748ede691862b9a87f67812db5e728a703b9d",
                "upstream_version": "1.2.22",
                "core_algorithms": "upstream_with_configurable_extensions",
                "adaptive_ms1_noise": runtime_config.noise_mode == "adaptive",
                "auto_tuning_extension": runtime_config.auto_tune,
                "profile_centroiding_extension": runtime_config.profile_mode == "centroid",
                "ms2_seeding_extension": runtime_config.ms2_seeded_features,
                "scan_precision": "float32_as_upstream",
            },
            "features": count,
            **stats,
            "elapsed_seconds": elapsed,
            "output_format": output_format,
            "output_columns": feature_columns(runtime_config),
            "output_precision": OUTPUT_PRECISION,
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
            str(source),
            str(output),
            "ok",
            count,
            stats["ms1_scans"],
            time.perf_counter() - start,
            acquisition=acquisition.to_dict(),
        )
    finally:
        if grouping_spectra is not None:
            grouping_spectra.close()
        for path in temporary:
            path.unlink(missing_ok=True)
        lock.unlink(missing_ok=True)
