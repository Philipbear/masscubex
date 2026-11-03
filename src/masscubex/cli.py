"""Bounded file scheduling and a streaming JSON-lines job report."""

from concurrent.futures import FIRST_COMPLETED, ProcessPoolExecutor, ThreadPoolExecutor, wait
from concurrent.futures.process import BrokenProcessPool
from contextlib import closing
from dataclasses import asdict
from pathlib import Path
import argparse
import csv
import hashlib
import json
import multiprocessing
import os
import sys

from . import __version__
from .config import Config, TUNABLE_FIELDS
from .readers import file_format
from .reporting import BatchReport, InputFailure, error_result, retry_sources


def _supported(path):
    try:
        file_format(path)
        return True
    except ValueError:
        return False


def _walk_files(directory):
    # os.walk materializes every filename in a directory. Keep scandir iterators
    # instead so even a single directory containing millions of files is lazy.
    stack = []
    try:
        try:
            stack.append((directory, os.scandir(directory)))
        except OSError as error:
            yield InputFailure(
                error_result(directory, error, stage="discovery", input_kind="directory")
            )
            return
        while stack:
            try:
                entry = next(stack[-1][1])
            except StopIteration:
                stack.pop()[1].close()
                continue
            except OSError as error:
                failed, iterator = stack.pop()
                iterator.close()
                yield InputFailure(
                    error_result(failed, error, stage="discovery", input_kind="directory")
                )
                continue
            try:
                if entry.is_dir(follow_symlinks=False):
                    try:
                        stack.append((entry.path, os.scandir(entry.path)))
                    except OSError as error:
                        yield InputFailure(
                            error_result(
                                entry.path, error, stage="discovery", input_kind="directory"
                            )
                        )
                elif entry.is_file() and _supported(entry.name):
                    yield Path(entry.path)
            except OSError as error:
                yield InputFailure(error_result(entry.path, error, stage="discovery"))
    finally:
        for _, iterator in stack:
            iterator.close()


def iter_inputs(paths, manifest=None, retry_errors=None):
    """Lazy directory/manifest traversal: never collect a million input paths."""
    for raw in paths:
        path = Path(raw)
        try:
            if path.is_dir():
                yield from _walk_files(path)
            else:
                yield path
        except OSError as error:
            yield InputFailure(error_result(path, error, stage="discovery"))
    if manifest:
        manifest = Path(manifest).absolute()
        try:
            with manifest.open(encoding="utf-8") as stream:
                for line in stream:
                    line = line.strip()
                    if line and not line.startswith("#"):
                        path = Path(line)
                        yield path if path.is_absolute() else manifest.parent / path
        except (OSError, UnicodeError) as error:
            yield InputFailure(
                error_result(manifest, error, stage="discovery", input_kind="manifest")
            )
    if retry_errors:
        for source, kind in retry_sources(retry_errors):
            if kind == "manifest":
                yield from iter_inputs([], manifest=source)
            else:
                yield from iter_inputs([source])


def output_for(source, output_dir, output_format="parquet"):
    """128-bit path IDs avoid collisions; a prefix shards large output folders."""
    path = Path(source).resolve()
    digest = hashlib.sha256(os.fsencode(path)).hexdigest()[:32]
    # Keep long input basenames from exceeding filesystem filename limits.
    name = path.name.encode("utf-8")[:100].decode("utf-8", errors="ignore")
    if output_format not in {"parquet", "tsv"}:
        raise ValueError("Output format must be parquet or tsv")
    return Path(output_dir) / digest[:2] / f"{name}.{digest}.features.{output_format}"


def _job(source, output_dir, config, resume, overwrite, output_format="parquet"):
    from .pipeline import process_file

    try:
        return asdict(
            process_file(
                source,
                output_for(source, output_dir, output_format),
                config=config,
                resume=resume,
                overwrite=overwrite,
            )
        )
    except Exception as error:
        return error_result(source, error)


def run_batch(
    paths,
    output_dir,
    config,
    *,
    workers=1,
    executor="process",
    resume=False,
    overwrite=False,
    output_format="parquet",
):
    """Yield each outcome as soon as ready, with at most 2 * workers futures."""
    if type(workers) is not int or workers < 1:
        raise ValueError("workers must be a positive integer")
    if executor not in {"process", "thread"}:
        raise ValueError("executor must be process or thread")
    if output_format not in {"parquet", "tsv"}:
        raise ValueError("Output format must be parquet or tsv")
    if workers == 1:
        for path in paths:
            yield (
                path.result
                if isinstance(path, InputFailure)
                else _job(path, output_dir, config, resume, overwrite, output_format)
            )
        return
    from .features import warmup

    # Validate compilation in the parent so a broken environment is not treated
    # as a sequence of corrupt files. Child workers load these cached kernels.
    warmup()

    def new_pool():
        return (
            ThreadPoolExecutor(max_workers=workers)
            if executor == "thread"
            else ProcessPoolExecutor(
                max_workers=workers,
                initializer=warmup,
                mp_context=multiprocessing.get_context("spawn"),
            )
        )

    pool = new_pool()
    iterator = iter(paths)
    pending = {}
    try:
        exhausted = False
        while pending or not exhausted:
            while not exhausted and len(pending) < 2 * workers:
                try:
                    path = next(iterator)
                except StopIteration:
                    exhausted = True
                    break
                if isinstance(path, InputFailure):
                    yield path.result
                    continue
                try:
                    future = pool.submit(
                        _job, path, output_dir, config, resume, overwrite, output_format
                    )
                except BrokenProcessPool:
                    # This input was not submitted. Drain the failed pool first,
                    # then dispatch it into a fresh one without losing its path.
                    for future, affected in pending.items():
                        try:
                            yield future.result()
                        except Exception as error:
                            yield error_result(affected, error, stage="worker")
                    pending.clear()
                    pool.shutdown(wait=True, cancel_futures=True)
                    pool = new_pool()
                    future = pool.submit(
                        _job, path, output_dir, config, resume, overwrite, output_format
                    )
                pending[future] = path
            if pending:
                done, _ = wait(pending, return_when=FIRST_COMPLETED)
                for future in done:
                    path = pending.pop(future)
                    try:
                        yield future.result()
                    except Exception as error:
                        yield error_result(path, error, stage="worker")
    finally:
        pool.shutdown(wait=True, cancel_futures=True)


# Every processing setting is available through JSON; routine controls also have flags.
_SETTING_FLAGS = (
    ("tuning_min_peaks", "--tuning-min-peaks", int),
    ("tuning_max_peaks", "--tuning-max-peaks", int),
    ("tuning_search_ppm", "--tuning-search-ppm", float),
    ("tuning_search_da", "--tuning-search-da", float),
    ("tuning_peak_limit", "--tuning-peak-limit", int),
    ("tuning_block_scans", "--tuning-block-scans", int),
    ("tuning_block_seconds", "--tuning-block-seconds", float),
    ("mz_tolerance_ppm", "--mz-ppm", float),
    ("mz_tolerance_da", "--mz-da", float),
    ("quant_method", "--quant-method", str),
    ("noise_mode", "--noise", str),
    ("noise_estimator", "--noise-estimator", str),
    ("noise_fraction", "--noise-fraction", float),
    ("noise_window_da", "--noise-window-da", float),
    ("scan_cache_mb", "--scan-cache-mb", float),
    ("noise_sigma", "--noise-sigma", float),
    ("noise_clip_sigma", "--noise-clip-sigma", float),
    ("noise_iterations", "--noise-iterations", int),
    ("noise_min_neighbors", "--noise-min-neighbors", int),
    ("min_intensity", "--min-intensity", float),
    ("min_scans", "--min-scans", int),
    ("max_gap_scans", "--max-gap-scans", int),
    ("max_gap_seconds", "--max-gap-seconds", float),
    ("profile_mode", "--profile-mode", str),
    ("profile_min_points", "--profile-min-points", int),
    ("profile_gap_factor", "--profile-gap-factor", float),
    ("smoothing_sigma", "--smoothing-sigma", float),
    ("smoothing_seconds", "--smoothing-seconds", float),
    ("peak_distance_seconds", "--peak-distance-seconds", float),
    ("peak_prominence", "--peak-prominence", float),
    ("peak_distance", "--peak-distance", int),
    ("ms2_mz_tolerance_ppm", "--ms2-mz-ppm", float),
    ("ms2_mz_tolerance_da", "--ms2-mz-da", float),
    ("ms2_rt_tolerance", "--ms2-rt-tolerance", float),
    ("ms2_max_isolation_width", "--ms2-max-isolation-width", float),
    ("ms2_gapfill_rt_window", "--ms2-gapfill-rt-window", float),
    ("ms2_seed_rt_tolerance", "--ms2-seed-rt-tolerance", float),
    ("gapfill_mz_tolerance_ppm", "--gapfill-mz-ppm", float),
    ("gapfill_mz_tolerance_da", "--gapfill-mz-da", float),
    ("precursor_max_age_seconds", "--precursor-max-age-seconds", float),
    ("gapfill_max_gap_seconds", "--gapfill-max-gap-seconds", float),
    ("best_ms2", "--best-ms2", str),
    ("group_mz_tolerance_ppm", "--group-mz-ppm", float),
    ("group_mz_tolerance_da", "--group-mz-da", float),
    ("group_rt_tolerance", "--group-rt-tolerance", float),
    ("group_min_scans", "--group-min-scans", int),
    ("group_min_correlation", "--group-min-correlation", float),
    ("group_min_overlap", "--group-min-overlap", float),
    ("isf_min_relative_intensity", "--isf-min-relative-intensity", float),
    ("isf_max_peaks", "--isf-max-peaks", int),
)
_CHOICES = {
    "quant_method": ("peak_height", "peak_area"),
    "noise_mode": ("adaptive", "none"),
    "noise_estimator": ("low_tail", "clipped"),
    "noise_sigma": (3.0,),
    "profile_mode": ("centroid", "reject"),
    "best_ms2": ("fragment_intensity", "precursor_intensity", "apex"),
}


def _config_options(parser):
    parser.add_argument(
        "--config", type=Path, help="JSON processing configuration; explicit flags override it"
    )
    parser.add_argument("--adducts", help="Ion IDs/ranges or exact names, e.g. 1-15,22,23")
    parser.add_argument(
        "--adduct-file",
        type=Path,
        help="JSON list of additional ion definitions with unique IDs/names",
    )


def _resolve_config(args):
    from .adducts import Adduct, parse_selection

    values = json.loads(args.config.read_text()) if args.config else {}
    if not isinstance(values, dict):
        raise ValueError("Configuration must be a JSON object")
    if args.adduct_file:
        values["custom_adducts"] = json.loads(args.adduct_file.read_text())
    for name in Config.__dataclass_fields__:
        if hasattr(args, name):
            values[name] = getattr(args, name)
    locked = set(values.get("tuning_locked", []))
    locked.update(name for name in TUNABLE_FIELDS if hasattr(args, name))
    if hasattr(args, "smoothing_sigma"):
        locked.add("smoothing_seconds")
    if hasattr(args, "peak_distance"):
        locked.add("peak_distance_seconds")
    values["tuning_locked"] = sorted(locked)
    if args.adducts:
        custom = tuple(Adduct(**row) for row in values.get("custom_adducts", []))
        values["adduct_ids"] = list(parse_selection(args.adducts, custom))
    return Config.from_dict(values)


def _parser():
    parser = argparse.ArgumentParser(
        description="Stream mzML/mzXML to feature Parquet/TSV tables with calibration, MS2 references and ion grouping."
    )
    parser.add_argument("--version", action="version", version=__version__)
    sub = parser.add_subparsers(dest="command", required=True)
    extract = sub.add_parser("extract", help="Extract files, directory trees, or a path manifest")
    extract.add_argument("inputs", nargs="*", help="Files or recursively scanned directories")
    extract.add_argument("-o", "--output-dir")
    extract.add_argument("--manifest", help="UTF-8 paths, one per line, relative to the manifest")
    extract.add_argument(
        "--retry-errors", type=Path, help="Retry inputs from a previous errors.tsv report"
    )
    extract.add_argument(
        "--report-dir",
        type=Path,
        help="Parent for unique run reports (default: OUTPUT_DIR/reports)",
    )
    extract.add_argument("--workers", type=int, default=1, help="File workers (default: 1)")
    extract.add_argument("--executor", choices=("process", "thread"), default="process")
    extract.add_argument(
        "--format",
        choices=("parquet", "tsv"),
        default="parquet",
        help="Feature table format (default: parquet, Zstd compressed)",
    )
    mode = extract.add_mutually_exclusive_group()
    mode.add_argument(
        "--resume", action="store_true", help="Verify and skip completed jobs; repair stale results"
    )
    mode.add_argument("--overwrite", action="store_true")
    _config_options(extract)
    extract.add_argument(
        "--print-config", action="store_true", help="Print resolved JSON settings and exit"
    )
    defaults = Config()
    for name, flag, type_ in _SETTING_FLAGS:
        help_ = name.replace("_", " ") + f" (default: {getattr(defaults, name)})"
        if name in TUNABLE_FIELDS:
            help_ += "; calibrated unless overridden; shown value is the fallback"
        if name in {"ms2_rt_tolerance", "group_rt_tolerance"}:
            help_ += "; seconds"
        if name == "ms2_max_isolation_width":
            help_ += "; full width in Da, 0 disables the optional MS2 association cutoff"
        if name == "ms2_gapfill_rt_window":
            help_ += "; seconds on each side of the first unmatched MS2 in a seeded feature"
        if name == "noise_fraction":
            help_ += "; weakest positive fraction in (0, 1], low_tail only"
        if name == "noise_sigma":
            help_ += "; fixed at 3; tune the background subset with --noise-fraction"
        if name == "scan_cache_mb":
            help_ += "; per-worker decoded scan budget in MiB; 0 disables replay"
        if name == "noise_min_neighbors":
            help_ += "; minimum selected background count, excluding the candidate"
        if name in {"noise_clip_sigma", "noise_iterations"}:
            help_ += "; clipped estimator only"
        extract.add_argument(
            flag,
            dest=name,
            type=type_,
            default=argparse.SUPPRESS,
            choices=_CHOICES.get(name),
            help=help_,
        )
    for name, flag, help_ in (
        (
            "auto_tune",
            "--auto-tune",
            "Estimate mass/timing settings independently for each polarity",
        ),
        ("segment", "--segment", "Split chromatographic peaks"),
        (
            "assume_centroid",
            "--assume-centroid",
            "Allow missing mode metadata; never relabel profile scans",
        ),
        ("associate_ms2", "--ms2", "Associate MS2 scans and select a representative scan"),
        (
            "ms2_seeded_features",
            "--ms2-seeded-features",
            "Create unmatched MS2 precursor features and gap-fill their MS1 measurements",
        ),
        (
            "group_ions",
            "--group-ions",
            "Annotate candidate adduct and source-fragment relationships",
        ),
        (
            "group_isf",
            "--isf",
            "Include water-loss and MS2-supported fragment candidates during grouping",
        ),
    ):
        extract.add_argument(
            flag,
            dest=name,
            action=argparse.BooleanOptionalAction,
            default=argparse.SUPPRESS,
            help=help_ + f" (default: {getattr(defaults, name)})",
        )
    sub.add_parser("warmup", help="Compile/load numerical kernels before a large run")
    sub.add_parser("defaults", help="Print all default processing parameters as JSON")
    adducts = sub.add_parser("adducts", help="List numbered ion forms and their default selection")
    _config_options(adducts)
    adducts.add_argument(
        "--json", action="store_true", help="Print machine-readable ion definitions"
    )
    return parser


def main(argv=None):
    parser = _parser()
    args = parser.parse_args(argv)
    if args.command == "warmup":
        from .features import warmup

        warmup()
        print("MassCubeX kernels ready", file=sys.stderr)
        return 0
    if args.command == "defaults":
        print(json.dumps(Config().to_dict(), indent=2))
        return 0
    try:
        config = _resolve_config(args)
    except (OSError, ValueError, TypeError) as error:
        parser.error(str(error))
    if args.command == "adducts":
        from .adducts import catalog_rows

        rows = catalog_rows(config.custom_adducts, config.adduct_ids)
        if args.json:
            print(json.dumps(rows, indent=2))
        else:
            print("id\tname\tpolarity\tcharge\tmolecules\tionic_mass_shift\tkind\tenabled")
            for row in rows:
                print(
                    f"{row['id']}\t{row['name']}\t{row['polarity']}\t{row['charge']}\t{row['molecules']}\t{row['mass_shift']:.9f}\t{row['kind']}\t{str(row['enabled']).lower()}"
                )
        return 0
    if args.print_config:
        print(json.dumps(config.to_dict(), indent=2))
        return 0
    if not args.inputs and not args.manifest and not args.retry_errors:
        parser.error("Provide inputs, --manifest, or --retry-errors")
    if not args.output_dir:
        parser.error("Provide --output-dir")
    if args.workers < 1:
        parser.error("--workers must be positive")
    report_root = args.report_dir or Path(args.output_dir) / "reports"
    options = {
        "inputs": args.inputs,
        "manifest": args.manifest,
        "retry_errors": str(args.retry_errors) if args.retry_errors else None,
        "output_dir": str(Path(args.output_dir).absolute()),
        "workers": args.workers,
        "executor": args.executor,
        "output_format": args.format,
        "resume": args.resume,
        "overwrite": args.overwrite,
    }
    try:
        with BatchReport(report_root, config, options) as report:
            print(f"Batch report: {report.summary_path}", file=sys.stderr)
            with closing(
                run_batch(
                    iter_inputs(args.inputs, args.manifest, args.retry_errors),
                    args.output_dir,
                    config,
                    workers=args.workers,
                    executor=args.executor,
                    output_format=args.format,
                    resume=args.resume,
                    overwrite=args.overwrite,
                )
            ) as results:
                for result in results:
                    report.record(result)
                    print(json.dumps(result), flush=True)
            if report.counts["total"] == 0:
                raise ValueError("No supported input files found")
    except KeyboardInterrupt:
        print(
            "Interrupted; completed outcomes and error reports have been preserved.",
            file=sys.stderr,
        )
        return 130
    except (OSError, ValueError, RuntimeError, csv.Error) as error:
        print(str(error), file=sys.stderr)
        return 1
    total, errors = report.counts["total"], report.counts["failed"]
    print(f"{total} outcomes; {errors} errors; summary: {report.summary_path}", file=sys.stderr)
    return 1 if errors else 0
