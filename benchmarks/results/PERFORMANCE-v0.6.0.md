# Extraction performance: MassCubeX 0.6.0

This compares the previous MassCubeX 0.5.1 extraction workflow at commit
`3bf31e4ac4609eb6f80435969bfec935f1c46bd5` with 0.6.0. It is not a timing
comparison against upstream MassCube's different scientific workflow.

## Results

| Workflow | 50-file median | Speed versus 0.5.1 | Feature tables |
| --- | ---: | ---: | ---: |
| 0.5.1, TSV, one worker | 253.23 s | 1.00× | 431.85 MB |
| 0.6.0, TSV, one worker | 163.75 s | 1.55× | 431.85 MB |
| 0.6.0, Parquet, one worker | 145.60 s | 1.74× | 85.72 MB |
| 0.6.0, Parquet, four process workers | 42.19 s | 6.00× including parallelism | 85.72 MB |

The extraction optimizations alone reduce elapsed time by **35.3%**. With
Parquet as the default, elapsed time falls **42.5%** and feature-table size falls
**80.2%** (5.04× smaller). MB means 1,000,000 bytes. JSON sidecars add about
2.4 MB per 50-file set; including them, output size falls from 434.3 to 88.2 MB.

The three complete-pass times were:

- Previous TSV: 252.48, 253.24, 253.23 seconds.
- Optimized TSV: 163.46, 163.75, 164.50 seconds.
- Optimized Parquet: 145.71, 145.47, 145.60 seconds.
- Optimized Parquet, four workers: 42.25, 41.88, 42.19 seconds.

Four workers give 3.45× the throughput of one optimized worker on this machine,
including process startup. All 50 parallel Parquet tables are byte-identical to
their single-worker counterparts. This scaling result is separate from the
1.74× improvement attributable to the code and default output format.

Raw reports: [previous TSV](extraction-v0.5.1-tsv.json),
[optimized TSV](extraction-v0.6.0-tsv.json),
[optimized Parquet](extraction-v0.6.0-parquet.json),
[four workers](extraction-v0.6.0-parquet-4-workers.json).

Single-process peak RSS was 472 MiB with optimized TSV and 589 MiB with Parquet.
Parquet uses less disk space but adds Arrow/library and conversion memory in this
run. Choose worker counts with per-file memory needs in mind.

## Where time goes

The production timers from the final 50-file, one-worker Parquet pass show:

| Stage | Total for 50 files | Share of elapsed time |
| --- | ---: | ---: |
| MS1 extraction, including input decoding, noise and peak summarization | 81.08 s | 55.7% |
| Automatic calibration, including its input pass | 21.13 s | 14.5% |
| Ion grouping | 16.64 s | 11.4% |
| MS2-seeded raw-MS1 gap filling | 12.74 s | 8.8% |
| Parquet table writing, including fsync | 9.47 s | 6.5% |
| MS2 association | 3.78 s | 2.6% |
| Checksums, metadata, locks and remaining overhead | 0.76 s | 0.5% |

The corresponding optimized TSV pass spent 27.13 seconds writing tables.
Serialization falls by about 65% with Parquet; peak extraction now dominates
the remaining runtime. Stage numbers describe one pass, while the end-to-end
comparison above uses the median of three.

Separate instrumentation on the same biological sample, blank and QC (three
repeats each, both versions writing TSV) isolates the main hotspots:

| Operation | Previous mean/file | Optimized mean/file | Speedup |
| --- | ---: | ---: | ---: |
| Ion grouping | 1.896 s | 0.345 s | 5.49× |
| Chromatographic peak summarization (inside MS1 extraction) | 1.454 s | 1.080 s | 1.35× |
| MS1 noise estimation (inside MS1 extraction) | 0.192 s | 0.192 s | unchanged |
| Extraction-pass decoding (inside MS1 extraction) | 0.236 s | 0.233 s | essentially unchanged |

Ion grouping was the largest individual hotspot; **chromatographic peak
summarization is now the main remaining hotspot**. The noise algorithm was not
the dominant cost and was left unchanged. These substeps overlap the parent
MS1-extraction stage and must not be added to that stage. See the
[previous profile](profile-v0.5.1-tsv.json) and
[optimized profile](profile-v0.6.0-tsv.json); a separate
[Parquet profile](profile-v0.6.0-parquet.json) is also included.

## Correctness

The [complete output comparison](extraction-equivalence-v0.6.0.json) verifies all
50 files: **1,540,612 features**, **125,266 associated MS2 scans**, **5,143 seeded
features** and **205,868 relationship entries**. Every optimized TSV is
byte-identical to 0.5.1. Every decoded Parquet scalar, scan list and relationship
matches its TSV counterpart, in the same row order. Scientific sidecar content,
including calibration decisions and unresolved-MS2 bookkeeping, is unchanged.
This preserves the results of the [prior independent raw-data audit](MSV000086989-ms2-v0.5.1.json).

Single-worker reports were generated at commit `6bdbc54`. The subsequent
Windows-only fix opens the final flush handle as writable on Windows, where
`fsync` requires it; the measured macOS path and numerical algorithms are
unchanged. Reports retain the implementation fingerprints actually measured.

## Measurement conditions

- Apple M2 Max, 32 GiB RAM, macOS 26.6.2 ARM64, Python 3.10.12.
- The same 50 public [MSV000086989 files](https://massive.ucsd.edu/ProteoSAFe/dataset.jsp?task=2f73277b4e034948acebfdf1edab17ed):
  40 biological samples, six blanks and four QCs; 604,070,716 input bytes.
- Default scientific settings throughout: automatic calibration, adaptive MS1
  noise, segmentation, MS2 assignment/seeding/gap filling and ion grouping all
  enabled. No thresholds, adduct rules or precision were relaxed for speed.
- Three complete passes per benchmark, with one worker unless stated otherwise.
  Use the median elapsed time for the entire 50-file pass, not the sum of median
  individual-file times. Inputs and output directories are on the local SSD.
- Warm filesystem cache and JIT. `warmup()` plus one complete file run are
  excluded; measured runs overwrite outputs and include parsing, computation,
  serialization, fsync, checksums and completion-sidecar publication. Downloads
  and correctness audits run separately. No simultaneous tests or benchmarks.
- BLAS/OpenMP/Numba thread counts are limited to one per worker. Multiple-worker
  measurements use spawned processes and include their startup for every pass.
- Dependencies and the exact implementation fingerprint are recorded in each
  benchmark JSON. The local run used NumPy 1.26.0, Numba 0.60.0, SciPy 1.13.1,
  lxml 5.2.2, PyArrow 17.0.0 and pynumpress 0.0.9.

## What changed

1. **Compiled ion-candidate rejection.** Vectorized mass searches feed a Numba
   iterator that rejects incompatible apex times, empty traces and self-pairs
   before entering Python. It streams compatible pairs instead of allocating a
   large pair matrix. Per-feature charges and timing data are computed once.
   MS2 fragment searches use the same filter and vectorized mass bounds.
2. **Less peak-summarization overhead.** Gaussian kernels are cached, padding is
   allocated directly, and median spacing uses a compiled helper. A plateau-aware
   local-maximum check skips SciPy prominence/distance work only when there cannot
   be two peaks. The multiple-peak path and quantitative arithmetic are preserved.
3. **Parquet output.** Zstd level 3, float64 measurements, nullable int64 scan IDs,
   native lists and typed relationship structs replace repeated text/JSON
   serialization. Row groups contain at most 8,192 features; categorical strings
   use dictionary encoding. `--format tsv` retains text output.
4. **Production timings.** Every new sidecar includes `timings_seconds`, so future
   datasets can be profiled without enabling cProfile or editing the source.

## Reproduction

Run from the repository root in an installed environment. Raw data, archived
source and generated tables stay under the ignored `validation/` directory.

```bash
python benchmarks/download_massive.py --directory validation/MSV000086989
mkdir -p validation/performance/baseline-v0.5.1
git archive 3bf31e4ac4609eb6f80435969bfec935f1c46bd5 | tar -x -C validation/performance/baseline-v0.5.1
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 NUMBA_NUM_THREADS=1

PYTHONPATH=validation/performance/baseline-v0.5.1/src python benchmarks/benchmark_extraction.py \
  --output-directory validation/performance/baseline-tsv \
  --report benchmarks/results/extraction-v0.5.1-tsv.json
PYTHONPATH=src python benchmarks/benchmark_extraction.py \
  --output-directory validation/performance/optimized-tsv \
  --report benchmarks/results/extraction-v0.6.0-tsv.json
PYTHONPATH=src python benchmarks/benchmark_extraction.py --format parquet \
  --output-directory validation/performance/optimized-parquet \
  --report benchmarks/results/extraction-v0.6.0-parquet.json
PYTHONPATH=src python benchmarks/benchmark_extraction.py --format parquet --workers 4 \
  --output-directory validation/performance/parquet-4-workers \
  --report benchmarks/results/extraction-v0.6.0-parquet-4-workers.json
PYTHONPATH=src python benchmarks/compare_extraction_outputs.py
```

`profile_extraction.py` separately times one biological sample, one blank and one
QC, each three times. Select the archived/current source with `PYTHONPATH` as
above and pass `--output-directory`, `--report` and `--format`. Its extraction
substeps are included in the parent extraction time; do not add them twice.
Its shared trace-tracking timer includes calibration and extraction and is
reported separately, so it must not be attributed entirely to extraction.
Lightweight timers still add overhead, so use the unprofiled 50-file measurements
for end-to-end speed comparisons.

## Scope

These are positive-mode, centroid mzXML files from one study. Genuine profile
mzML, polarity-switching acquisitions, substantially larger files, network
storage and cold-start latency are not represented by these timings. Synthetic
tests exercise those acquisition modes, but do not establish their throughput.
No million-file run was performed. Peak RSS, where reported, includes imports,
JIT/library caches and all processing stages in the single benchmark process;
it is not just the raw-spectrum buffer or the sum across parallel workers.
