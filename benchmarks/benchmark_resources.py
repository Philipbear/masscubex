"""Measure file throughput and the resident memory of explicitly owned workers.

Workers are warmed before the steady-state timer starts. A separate startup
timer includes imports, JIT/cache loading and one whole-file warmup per worker.
No system process enumeration, cache purging or changes to scientific settings
are used. psutil is needed only for this benchmark, not for extraction.
"""

import argparse
from collections import Counter
import hashlib
import json
import multiprocessing
import os
from pathlib import Path
import platform
import queue
import resource
import sys
import time
import traceback


def digest(path):
    checksum = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            checksum.update(block)
    return checksum.hexdigest()


def scientific(value):
    if isinstance(value, dict):
        return {
            k: scientific(v)
            for k, v in value.items()
            if k not in {"elapsed_seconds", "timings_seconds", "implementation_sha256"}
        }
    if isinstance(value, list):
        return [scientific(v) for v in value]
    return value


def remove_result(output):
    output.unlink(missing_ok=True)
    Path(str(output) + ".json").unlink(missing_ok=True)


def relieve_allocators(process):
    """Measure best-effort release of this process's unused allocations only."""
    import ctypes
    import gc
    import pyarrow as pa

    start = time.perf_counter()
    result = {"before_bytes": process.memory_info().rss}
    result["gc_collected"] = gc.collect()
    result["after_gc_bytes"] = process.memory_info().rss
    pool = pa.default_memory_pool()
    result["arrow_live_bytes"] = pool.bytes_allocated()
    pool.release_unused()
    result["after_arrow_bytes"] = process.memory_info().rss
    libc = ctypes.CDLL(None)
    if platform.system() == "Darwin":
        function = libc.malloc_zone_pressure_relief
        function.argtypes = [ctypes.c_void_p, ctypes.c_size_t]
        function.restype = ctypes.c_size_t
        result["native_released_bytes"] = function(None, 0)
    elif platform.system() == "Linux" and hasattr(libc, "malloc_trim"):
        function = libc.malloc_trim
        function.argtypes = [ctypes.c_size_t]
        function.restype = ctypes.c_int
        result["malloc_trim_result"] = function(0)
    result["after_native_bytes"] = process.memory_info().rss
    result["seconds"] = time.perf_counter() - start
    return result


def worker(
    incoming, outgoing, package_path, directory, config_dict, first_input,
    release_every=0, worker_id=0,
):
    """Keep scientific objects inside a worker; send compact result records."""
    if package_path:
        sys.path.insert(0, package_path)
    import psutil
    from masscubex import Config, process_file
    from masscubex.features import warmup
    from masscubex.pipeline import implementation_id

    process = psutil.Process()
    config = Config.from_dict(config_dict) if config_dict else Config()
    output = Path(directory) / f"warmup-{os.getpid()}-{worker_id}.parquet"
    try:
        warmup()
        process_file(first_input, output, config=config, overwrite=True)
        remove_result(output)
        outgoing.put({
            "kind": "ready", "pid": os.getpid(),
            "rss_bytes": process.memory_info().rss,
            "implementation_sha256": implementation_id()[0],
            "dependencies": implementation_id()[1], "config": config.to_dict(),
        })
        for task_number, task in enumerate(iter(incoming.get, None), 1):
            index, source = task
            output = Path(directory) / (Path(source).name + ".parquet")
            try:
                result = process_file(source, output, config=config, overwrite=True)
                metadata = json.loads(Path(str(output) + ".json").read_text())
                encoded = json.dumps(scientific(metadata), sort_keys=True).encode()
                row = {
                    "kind": "file", "index": index, "pid": os.getpid(),
                    "name": Path(source).name, "seconds": result.elapsed_seconds,
                    "features": result.features, "stages": metadata["timings_seconds"],
                    "input_sha256": digest(Path(source)),
                    "output_sha256": metadata["output_sha256"],
                    "scientific_sha256": hashlib.sha256(encoded).hexdigest(),
                    "output_bytes": metadata["output_bytes"],
                    "rss_after_bytes": process.memory_info().rss,
                    "peak_rss_bytes": resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
                    * (1 if platform.system() == "Darwin" else 1024),
                }
                remove_result(output)
                if release_every and task_number % release_every == 0:
                    row["allocator_relief"] = relieve_allocators(process)
                    row["rss_after_bytes"] = process.memory_info().rss
                outgoing.put(row)
            except Exception:
                outgoing.put({"kind": "error", "name": Path(source).name,
                              "error": traceback.format_exc()})
    except Exception:
        outgoing.put({"kind": "fatal", "error": traceback.format_exc()})
    finally:
        remove_result(output)


def thread_workers(count, arguments):
    """Host threads in one owned process so the parent can enforce its RAM budget."""
    from concurrent.futures import ThreadPoolExecutor

    with ThreadPoolExecutor(max_workers=count) as pool:
        futures = [pool.submit(worker, *arguments, worker_id=i) for i in range(count)]
        for future in futures:
            future.result()


def sample_memory(processes, monitor):
    rss = {}
    for process in processes:
        if process.pid is None:
            continue
        try:
            rss[str(process.pid)] = monitor.Process(process.pid).memory_info().rss
        except monitor.NoSuchProcess:
            pass
    return rss


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--directory", type=Path, default=Path("validation/MSV000086989"))
    parser.add_argument("--manifests", nargs="+", default=["download_manifest.json"])
    parser.add_argument("--package-path", type=Path)
    parser.add_argument("--output-directory", type=Path, required=True)
    parser.add_argument("--report", type=Path, required=True)
    parser.add_argument("--reference", type=Path)
    parser.add_argument("--config", type=Path)
    parser.add_argument("--workers", type=int, default=1)
    parser.add_argument("--executor", choices=["process", "thread"], default="process")
    parser.add_argument("--limit", type=int)
    parser.add_argument("--memory-budget-mib", type=float, default=3072)
    parser.add_argument("--sample-seconds", type=float, default=0.2)
    parser.add_argument("--release-every", type=int, default=0, help="Experimental allocator relief every N files (0: off)")
    args = parser.parse_args()
    if args.workers < 1 or args.memory_budget_mib <= 0 or args.sample_seconds <= 0:
        parser.error("Workers, memory budget and sampling interval must be positive")
    if args.release_every < 0:
        parser.error("--release-every cannot be negative")
    if args.release_every and args.executor == "thread":
        parser.error("Experimental allocator relief is measured in isolated processes only")
    import psutil

    rows = []
    for manifest in args.manifests:
        rows.extend(json.loads((args.directory / manifest).read_text()))
    rows = rows[:args.limit]
    if not rows or len({r["name"] for r in rows}) != len(rows):
        parser.error("Input manifests must identify distinct files")
    reference = {}
    if args.reference:
        reference = {r["name"]: r for r in json.loads(args.reference.read_text())["files"]}
    config = json.loads(args.config.read_text()) if args.config else None
    args.output_directory.mkdir(parents=True, exist_ok=True)
    args.report.parent.mkdir(parents=True, exist_ok=True)
    tasks = [(i, str((args.directory / r["name"]).absolute())) for i, r in enumerate(rows)]
    report = {
        "python": platform.python_version(), "platform": platform.platform(),
        "workers": args.workers, "executor": args.executor, "input_files": len(rows),
        "input_bytes": sum(Path(p).stat().st_size for _, p in tasks),
        "memory_budget_mib": args.memory_budget_mib, "files": [], "samples": [],
        "allocator_relief_every": args.release_every,
        "notes": [
            "Warm filesystem/JIT; startup and one warmup per worker reported separately.",
            "RSS is summed over these workers only, including shared pages in each process.",
            "Output hashing, provenance checks and deletion add small measured bookkeeping costs.",
            "Raw files, source code and persistent reports are never deleted.",
        ],
    }
    ctx = multiprocessing.get_context("spawn")
    incoming, outgoing = ctx.Queue(maxsize=2 * args.workers), ctx.Queue()
    worker_arguments = (
        incoming, outgoing, str(args.package_path.absolute()) if args.package_path else None,
        str(args.output_directory.absolute()), config, tasks[0][1], args.release_every,
    )
    processes = (
        [ctx.Process(target=worker, args=worker_arguments) for _ in range(args.workers)]
        if args.executor == "process"
        else [ctx.Process(target=thread_workers, args=(args.workers, worker_arguments))]
    )
    start = time.perf_counter()
    last_sample, steady_start = -1.0, None
    ready, sent, completed = [], 0, 0
    peak_rss = peak_total = peak_controller = 0
    controller = psutil.Process()
    report["system_available_before_bytes"] = psutil.virtual_memory().available
    try:
        for process in processes:
            process.start()
        while completed < len(tasks):
            now = time.perf_counter()
            if now - last_sample >= args.sample_seconds:
                rss = sample_memory(processes, psutil)
                total = sum(rss.values())
                coordinator_rss = controller.memory_info().rss
                peak_rss = max(peak_rss, total)
                peak_controller = max(peak_controller, coordinator_rss)
                peak_total = max(peak_total, total + coordinator_rss)
                report["samples"].append({"seconds": now - start, "rss_bytes": total,
                                          "controller_rss_bytes": coordinator_rss})
                if total > args.memory_budget_mib * 1024**2:
                    raise MemoryError(f"Worker RSS exceeded {args.memory_budget_mib:g} MiB budget")
                last_sample = now
            try:
                message = outgoing.get(timeout=args.sample_seconds)
            except queue.Empty:
                failed = [p.exitcode for p in processes if p.exitcode is not None]
                if failed:
                    raise RuntimeError(f"Worker exited before finishing: {failed}")
                continue
            if message["kind"] in {"fatal", "error"}:
                raise RuntimeError(message["error"])
            if message["kind"] == "ready":
                ready.append(message)
                if len(ready) < args.workers:
                    continue
                steady_start = time.perf_counter()
                report["startup_seconds"] = steady_start - start
                report["worker_ready"] = ready
                report.update({k: ready[0][k] for k in ["implementation_sha256", "dependencies", "config"]})
            else:
                if reference:
                    old = reference[message["name"]]
                    for key in ["input_sha256", "output_sha256", "scientific_sha256"]:
                        if message[key] != old[key]:
                            raise AssertionError(f"{key} changed: {message['name']}")
                    message["exact_reference_match"] = True
                report["files"].append(message)
                completed += 1
                print(json.dumps({"done": completed, "total": len(tasks),
                                  "name": message["name"], "seconds": message["seconds"],
                                  "rss_mib": round(peak_rss / 1024**2, 1)}), flush=True)
            if steady_start is not None:
                while sent < len(tasks) and sent - completed < 2 * args.workers:
                    incoming.put(tasks[sent])
                    sent += 1
        report["wall_seconds"] = time.perf_counter() - steady_start
        report["files_per_hour"] = 3600 * len(tasks) / report["wall_seconds"]
        report["status"] = "ok"
    except Exception:
        report["status"] = "failed"
        report["error"] = traceback.format_exc()
        raise
    finally:
        report["peak_worker_rss_bytes"] = peak_rss
        report["peak_controller_rss_bytes"] = peak_controller
        report["peak_total_rss_bytes"] = peak_total
        report["system_available_after_bytes"] = psutil.virtual_memory().available
        if report.get("status") == "ok":
            for _ in range(args.workers):
                incoming.put(None)
        else:
            for process in processes:
                if process.is_alive():
                    process.terminate()
        for process in processes:
            if process.pid is None:
                continue
            process.join(timeout=10)
            if process.is_alive():
                process.kill()
                process.join(timeout=5)
        # A budget cancellation can interrupt finally blocks inside a worker.
        # Clean only these known benchmark targets after every worker has exited.
        targets = {Path(source).name + ".parquet" for _, source in tasks}
        targets.update(
            f"warmup-{p.pid}-{i}.parquet" for p in processes for i in range(args.workers)
        )
        for path in args.output_directory.iterdir():
            if path.name in targets or any(
                path.name in {name + ".json", name + ".lock"}
                or path.name.startswith("." + name + ".") and path.name.endswith(".partial")
                for name in targets
            ):
                path.unlink(missing_ok=True)
        stages = Counter()
        for row in report["files"]:
            stages.update(row["stages"])
        report["stage_cpu_wall_sum_seconds"] = dict(stages)
        report["files"].sort(key=lambda r: r["index"])
        args.report.write_text(json.dumps(report, indent=2) + "\n")
        incoming.close()
        outgoing.close()


if __name__ == "__main__":
    main()
