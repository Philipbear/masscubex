"""Exercise the real CLI scheduler, including pool recycling, with owned-PID RSS.

Repeated cycles reuse the same immutable files to test long-running memory use.
Output checksums and scientific provenance are checked after every outcome, then
the generated table is removed. Pool observation is confined to this benchmark;
it does not change scheduling, worker initialization or extraction functions.
"""

import argparse
from concurrent.futures import ProcessPoolExecutor
import hashlib
import itertools
import json
import multiprocessing
import os
from pathlib import Path
import platform
import queue
import signal
import sys
import time
import traceback

from benchmark_resources import digest, remove_result, scientific


_observer = None


class ObservedPool(ProcessPoolExecutor):
    """Expose only this executor's known child PIDs to the resource monitor."""

    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        _observer.put({"kind": "pool_created"})

    def submit(self, fn, *args, **kwargs):
        future = super().submit(fn, *args, **kwargs)
        _observer.put({"kind": "pids", "pids": list(self._processes or {})})
        return future

    def shutdown(self, **kwargs):
        pids = list(self._processes or {})
        try:
            super().shutdown(**kwargs)
        finally:
            _observer.put({"kind": "retired", "pids": pids})


def cli_host(outgoing, package_path, sources, directory, workers, cycles, recycle_files):
    global _observer
    if hasattr(os, "setsid"):
        os.setsid()  # Give this benchmark's descendants their own cancellation group.
    if package_path:
        sys.path.insert(0, package_path)
    from masscubex import Config
    from masscubex import cli

    _observer = outgoing
    cli.ProcessPoolExecutor = ObservedPool
    runtime = Path(cli.__file__).with_name("worker_runtime.py")
    outgoing.put({"kind": "ready", "process_group": hasattr(os, "setsid"),
                  "cli_sha256": digest(Path(cli.__file__)), "config": Config().to_dict(),
                  "worker_runtime_sha256": digest(runtime) if runtime.exists() else None})
    try:
        paths = itertools.chain.from_iterable(sources for _ in range(cycles))
        results = cli.run_batch(
            paths, directory, Config(), workers=workers, overwrite=True,
            worker_recycle_files=recycle_files,
        )
        for result in results:
            if result["status"] != "ok":
                raise AssertionError(result)
            source, output = Path(result["input"]), Path(result["output"])
            metadata = json.loads(Path(str(output) + ".json").read_text())
            assert digest(output) == metadata["output_sha256"]
            encoded = json.dumps(scientific(metadata), sort_keys=True).encode()
            outgoing.put({
                "kind": "file", "name": source.name, "seconds": result["elapsed_seconds"],
                "features": result["features"], "input_sha256": digest(source),
                "output_sha256": metadata["output_sha256"],
                "scientific_sha256": hashlib.sha256(encoded).hexdigest(),
                "implementation_sha256": metadata["implementation_sha256"],
                "stages": metadata["timings_seconds"],
            })
            remove_result(output)
        outgoing.put({"kind": "complete"})
    except Exception:
        outgoing.put({"kind": "error", "error": traceback.format_exc()})


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--directory", type=Path, default=Path("validation/MSV000086989"))
    parser.add_argument("--manifests", nargs="+", default=["download_manifest.json"])
    parser.add_argument("--package-path", type=Path)
    parser.add_argument("--output-directory", type=Path, required=True)
    parser.add_argument("--report", type=Path, required=True)
    parser.add_argument("--reference", type=Path)
    parser.add_argument("--workers", type=int, default=2)
    parser.add_argument("--cycles", type=int, default=3)
    parser.add_argument("--limit", type=int)
    parser.add_argument("--worker-recycle-files", type=int, default=8)
    parser.add_argument("--memory-budget-mib", type=float, default=3072)
    args = parser.parse_args()
    if args.workers < 2 or args.cycles < 1 or args.worker_recycle_files < 0 or args.memory_budget_mib <= 0:
        parser.error("Use at least two process workers, positive cycles and nonnegative recycling")
    import psutil

    rows = []
    for manifest in args.manifests:
        rows.extend(json.loads((args.directory / manifest).read_text()))
    rows = rows[:args.limit]
    sources = [str((args.directory / row["name"]).absolute()) for row in rows]
    if not sources or len(set(sources)) != len(sources):
        parser.error("Manifests must contain distinct input paths")
    reference = {}
    if args.reference:
        reference = {r["name"]: r for r in json.loads(args.reference.read_text())["files"]}
    args.output_directory.mkdir(parents=True, exist_ok=True)
    args.report.parent.mkdir(parents=True, exist_ok=True)
    report = {
        "python": platform.python_version(), "platform": platform.platform(),
        "harness_sha256": digest(Path(__file__)),
        "workers": args.workers, "cycles": args.cycles, "unique_inputs": len(sources),
        "total_jobs": len(sources) * args.cycles,
        "worker_recycle_files": args.worker_recycle_files, "files": [], "samples": [],
        "pool_count": 0, "memory_budget_mib": args.memory_budget_mib,
        "notes": [
            "Real run_batch scheduler; pool startup, retirement and verification bookkeeping included.",
            "Cycles reuse the same cohort; they are not additional independent input files.",
            "RSS covers only the owned coordinator and observed child workers; shared pages are counted per process.",
            "The monitor reads this executor's private process map; it never enumerates system processes.",
        ],
    }
    ctx = multiprocessing.get_context("spawn")
    outgoing = ctx.Queue()
    host = ctx.Process(target=cli_host, args=(
        outgoing, str(args.package_path.absolute()) if args.package_path else None,
        sources, str(args.output_directory.absolute()), args.workers, args.cycles,
        args.worker_recycle_files,
    ))
    started = time.perf_counter()
    host.start()
    owned = {host.pid: psutil.Process(host.pid)}
    seen_workers = set()
    peak = 0
    process_group, run_started = False, None
    last_sample = -1.0
    try:
        while True:
            now = time.perf_counter()
            if now - last_sample >= 0.2:
                total, host_rss = 0, 0
                for pid, process in list(owned.items()):
                    try:
                        rss = process.memory_info().rss
                        total += rss
                        if pid == host.pid:
                            host_rss = rss
                    except psutil.NoSuchProcess:
                        owned.pop(pid)
                peak = max(peak, total)
                report["samples"].append({"seconds": now - started, "rss_bytes": total,
                                          "coordinator_rss_bytes": host_rss})
                if total > args.memory_budget_mib * 1024**2:
                    raise MemoryError("CLI processes exceeded the benchmark RAM budget")
                last_sample = now
            try:
                message = outgoing.get(timeout=0.2)
            except queue.Empty:
                if host.exitcode is not None:
                    raise RuntimeError(f"CLI host exited without completion: {host.exitcode}")
                continue
            kind = message["kind"]
            if kind == "ready":
                process_group = message["process_group"]
                report["cli_sha256"] = message["cli_sha256"]
                report["worker_runtime_sha256"] = message.get("worker_runtime_sha256")
                report["config"] = message["config"]
                run_started = time.perf_counter()
            elif kind == "pool_created":
                report["pool_count"] += 1
            elif kind == "pids":
                for pid in message["pids"]:
                    if pid not in owned:
                        seen_workers.add(pid)
                        try:
                            owned[pid] = psutil.Process(pid)
                        except psutil.NoSuchProcess:
                            pass  # A tiny smoke job can finish before observation.
            elif kind == "retired":
                for pid in message["pids"]:
                    owned.pop(pid, None)
            elif kind == "file":
                if "implementation_sha256" not in report:
                    report["implementation_sha256"] = message["implementation_sha256"]
                assert message["implementation_sha256"] == report["implementation_sha256"]
                old = reference.get(message["name"])
                if old is None:
                    reference[message["name"]] = message
                else:
                    for key in ["input_sha256", "output_sha256", "scientific_sha256"]:
                        if message[key] != old[key]:
                            raise AssertionError(f"{key} changed: {message['name']}")
                message["exact_reference_match"] = True
                report["files"].append(message)
                if len(report["files"]) % 10 == 0:
                    print(json.dumps({"done": len(report["files"]), "total": report["total_jobs"],
                                      "pools": report["pool_count"], "rss_mib": peak / 1024**2}), flush=True)
            elif kind == "error":
                raise RuntimeError(message["error"])
            elif kind == "complete":
                assert len(report["files"]) == report["total_jobs"]
                report["wall_seconds"] = time.perf_counter() - run_started
                report["files_per_hour"] = 3600 * report["total_jobs"] / report["wall_seconds"]
                report["status"] = "ok"
                break
    except Exception:
        report["status"] = "failed"
        report["error"] = traceback.format_exc()
        raise
    finally:
        if report.get("status") != "ok":
            if process_group:
                try:
                    os.killpg(host.pid, signal.SIGTERM)
                except ProcessLookupError:
                    pass
            else:
                for process in owned.values():
                    try:
                        process.terminate()
                    except psutil.NoSuchProcess:
                        pass
        host.join(timeout=10)
        if host.is_alive():
            host.kill()
            host.join(timeout=5)
        report["peak_pipeline_rss_bytes"] = peak
        report["worker_processes_created"] = len(seen_workers)
        # Every output is known in advance; preserve unrelated files if a user
        # reuses a directory. Partial files/locks can survive forced cancellation.
        from masscubex.cli import output_for

        for source in sources:
            output = output_for(source, args.output_directory)
            remove_result(output)
            Path(str(output) + ".lock").unlink(missing_ok=True)
            if output.parent.exists():
                for partial in output.parent.glob("." + output.name + ".*.partial"):
                    partial.unlink()
                try:
                    output.parent.rmdir()
                except OSError:
                    pass
        args.report.write_text(json.dumps(report, indent=2) + "\n")
        outgoing.close()


if __name__ == "__main__":
    main()
