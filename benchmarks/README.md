# Benchmarks and validation

Start with the current reports:

- [Speed versus pinned MassCube 1.2.22](results/UPSTREAM_SPEED.md): identical adaptive and fixed-floor workloads, scientific parity, serialization boundaries.
- [100-file extraction speed and RAM](results/OVERNIGHT_OPTIMIZATION.md): sustained worker recycling, stage costs, capacity estimates and raw measurements.
- [Algorithm parity](results/UPSTREAM_RESTORATION.md) and [upstream comparison](../docs/UPSTREAM_COMPARISON.md).
- [Code audit](../docs/CODE_AUDIT.md): implementation scope and validation history.

## Reproduce measurements

Install development dependencies with `python -m pip install -e '.[dev]'`.
Download inputs using `download_massive.py`; the current portable input manifest
is [here](results/overnight-input-manifest-100.json). Raw downloads and generated
tables belong under ignored `validation/`, not in Git.

| Script | Purpose |
| --- | --- |
| `benchmark_exact_speed.py` | Repeated single-process extraction; byte-identical tables and scientific-sidecar comparison against a frozen implementation. |
| `benchmark_upstream_speed.py` | Five-file comparison against pinned original functions, with equivalent preprocessing and grouping. Requires Pyteomics and the pinned upstream Git objects. |
| `benchmark_resources.py` | Warm process/thread throughput and owned-process RSS. |
| `benchmark_cli_resources.py` | Actual scheduler startup/recycling, per-file parity and a benchmark-only RSS stop budget. |
| `benchmark_coordinator_memory.py` | Coordinator imports/startup RSS. |
| `audit_warmup.py` | Production kernel signatures before/after extraction. |
| `profile_extraction.py` | Inclusive extraction substeps; profile separately from timed benchmarks. |

Each script has `--help`. Reports retain configurations, input/output hashes and
implementation fingerprints. [The source archive](baselines/README.md) preserves
the pre-optimization reference without generated caches. Scientific-sidecar
hashes include local input paths/mtimes; regenerate the reference on another
checkout. Input and table checksums are portable.

## Historical investigations

Older reports in `results/` describe the implementation recorded in each report.
They are retained as evidence and may have different scientific defaults:

- [Earlier validation](results/VALIDATION.md)
- [Feature counts and output precision](results/OUTPUT_AND_COUNTS.md)
- [Blank/QC EIC review](results/EIC_REVIEW.md)
- [Peak-shape investigation](results/PEAK_SHAPES.md)
- [Pre-restoration performance](results/PERFORMANCE_CURRENT.md)

Use the current reports above for the present workflow. Generated tables,
binary profiles, temporary environments and source snapshots are disposable;
retain compact measurements, independent reference fixtures and input manifests.
