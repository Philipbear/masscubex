# Extraction performance: MassCubeX 0.8.0

**MS1 extraction time falls 12.0%; median complete single-worker processing falls 2.8%.**
The same 50 public files, scientific defaults and Parquet output were measured
three times per configuration. The previous version was freshly re-run from
`ca63c8c`; the optimized processing code is `6a78440`.
The earlier requested 10× end-to-end code speedup remains unachieved.

## Complete command and two million similar files

Apple M2 Max, 32 GiB RAM, Python 3.10.12, macOS 26.6.2 ARM64. This cohort contains
50 positive-mode **centroid mzXML** files: 604.07 MB total, **12.08 MB/file** on
average. All steps remain enabled: calibration, adaptive noise, peak detection,
MS2 association/seeding/raw gap filling, ion grouping and compressed output.
The cohort is [MSV000086989](https://massive.ucsd.edu/ProteoSAFe/dataset.jsp?task=2f73277b4e034948acebfdf1edab17ed):
40 biological samples, six blanks and four QCs.

The full **CLI command**, including startup, bounded scheduling and batch
reports, was run three times at each worker count. These measurements supply
the capacity projections:

| Process workers on this machine | Median / 50 files | Observed three-run range | Projection / 2 million |
| --- | ---: | ---: | ---: |
| 1 | 56.67 s | 56.59–57.22 s | 26.24 days |
| 4 | 17.80 s | 17.54–18.32 s | 8.24 days |
| 8 | 13.71 s | 12.51–13.98 s | 6.35 days |
| 12 | 19.81 s | 13.03–20.44 s | 9.17 days |

Eight workers gave the best CLI median. Twelve workers were more variable and
slower at the median; no linear scaling with local worker count is assumed.
The separate four-worker corruption/resume audit took 29.16 seconds for its
initial 50-valid-plus-five-invalid batch. It is a correctness run rather than
the matched throughput benchmark, but illustrates that operational times vary.

For eight workers, the central estimate is **6.35 days**;
an explicit allowance for 20% unavailable capacity increases it to
**7.93 days**. This allowance is
a planning assumption, not a confidence interval or a bound on slow files.
The formula is `median_batch_seconds / 50 × 2,000,000`. Startup is included in
these short batches and would be amortized over a longer-lived worker pool.
These are projections from 50 files, not a measured two-million-file run.

For an independent cluster of equivalent machines, a **linear, unmeasured**
projection from the eight-worker CLI result is:

| Total workers | Ideal processing estimate | With 20% unavailable capacity |
| --- | ---: | ---: |
| 32 (4 machines × 8) | 38.1 h | 47.6 h |
| 64 (8 machines × 8) | 19.0 h | 23.8 h |
| 128 (16 machines × 8) | 9.5 h | 11.9 h |

Cluster estimates require each machine to sustain the measured CPU and storage
throughput. A central disk or network bottleneck can erase projected scaling.
At this cohort's size/compression, two million files imply approximately
**24.16 TB input, 3.43 TB Parquet tables and 98.1 GB JSON sidecars** (decimal units).
Transferring the input alone at an assumed sustained 100 MB/s takes about
67.1 hours; at 1 GB/s it takes 6.7 hours. These are bandwidth calculations,
not measured transfer rates. Downloads, staging, retries and alignment are
excluded. Hot filesystem-cache measurements are optimistic for cold remote
storage. Larger/profile files and other acquisition patterns need separate
throughput measurements; file count alone does not predict their runtime.

[Machine-readable capacity calculations](capacity-v0.8.0-2million.json).

## Matched before/after processing benchmark

The processing API comparison uses the same 50-file manifest and order, three
passes each, with all scientific stages and output publication included.
The one-worker median falls **60.60 → 58.92 s**, a **2.8% time reduction**
(1.029× throughput). These measurements support the before/after claim;
the separately run CLI measurements above support the command capacity table.
Their differences are not an isolated estimate of command overhead.

Current API process-pool medians for 1/4/8/12 workers were respectively
**58.92 / 16.91 / 12.29 / 11.90 s** per 50 files. The CLI's greater high-worker
variability is retained in the report; the best isolated API timing is not used
to project cluster capacity. No unrelated benchmark, test or audit was run
concurrently with either measurement series.

## Where processing time goes

Exclusive production timers from the last single-worker pass. The per-file
column is the 50-file total divided by 50; shares use that same pass's elapsed
time, **not** a sum of overlapping profiler substeps.

| Stage | 0.7.0 / 50 files | 0.8.0 / 50 files | Current / file | Current share |
| --- | ---: | ---: | ---: | ---: |
| MS1 extraction | 22.03 s | 19.39 s | 0.388 s | 32.9% |
| Calibration, including XML/binary decoding | 14.54 s | 14.67 s | 0.293 s | 24.9% |
| Parquet writing and fsync | 9.54 s | 9.29 s | 0.186 s | 15.8% |
| Ion grouping | 8.20 s | 8.44 s | 0.169 s | 14.3% |
| Raw-MS1 seeded gap filling | 3.67 s | 3.66 s | 0.073 s | 6.2% |
| MS2 association | 2.55 s | 2.78 s | 0.056 s | 4.7% |
| Checksums, completion records and orchestration | — | 0.71 s | 0.014 s | 1.2% |

MS1 extraction remains the largest stage. Calibration/decoding and output
remain substantial; optimizing MS1 alone cannot produce a 10× complete-workflow
speedup. A separate sample/blank/QC profile (three repetitions each) measures
noise estimation at **0.141 → 0.106 s/file**, a
**1.32×** kernel-stage speedup.
Peak-summarization timer boundaries changed: 0.8 includes native trace draining,
which 0.7 charged to tracking. Those isolated timer values are not directly
comparable. Shared reader/tracker timers span several stages and must not be
added to the stage table.

## Implementation and equivalence

- Noise moments are evaluated lazily only for queried low-tail subtrees.
  Counts remain current. Prefix reuse guards the first unselected rank too,
  preserving the original merge order and float64 results.
- Completed traces are drained, segmented and quantified in native scan batches.
  Calibration retains its raw-trace interface. Named settings and candidate
  records make options, measurements and rare SciPy fallback cases explicit.
- Gaussian reflection is materialized once per trace instead of repeated inside
  every convolution pair. Pool copying uses direct scalar loops. Scientific
  defaults and float64 precision are unchanged; fast-math remains disabled.

The [full table comparison](extraction-equivalence-v0.8.0.json) finds **zero
numerical or structural differences** in 1,540,612 features, 125,266 associated
MS2 scans, 5,143 seeded features and 205,868 ion-relationship entries. All
scientific sidecar values match. All 50 Parquet files are also byte-identical
to 0.7.0 and to every tested API/CLI worker configuration
([byte comparison](parallel-equivalence-v0.8.0.json)). The [noise audit](noise-equivalence-v0.8.0.json)
also finds exact means, standard deviations and thresholds at **34,905,337 raw
MS1 points across 55,198 scans**, including points removed before detection.
The 50 Da, lowest 10%, minimum five, leave-one-out μ + 3σ rule is unchanged.

Whole-process peak RSS, including warmup/JIT work, was
**673.3 → 722.7 MiB** in these runs.
These are high-water marks, not per-file or arbitrary-input memory guarantees.
The decoded scan-cache budget remains 64 MiB per worker; traces, features,
fragment evidence and library/compiler buffers consume additional memory.

See [validation](VALIDATION.md), [methods](../../docs/METHODS.md) and
[code structure](../../docs/DEVELOPMENT.md). The previous release's broader
optimization history remains in the [0.7.0 report](PERFORMANCE-v0.7.0.md).

## Repetition details and reproduction

Baseline one-worker passes: 60.00, 60.60, 61.24 s.
Current **API** passes:

- 1 worker(s): 58.75, 58.92, 58.94 s.
- 4 worker(s): 17.58, 16.91, 16.87 s.
- 8 worker(s): 10.73, 12.29, 12.44 s.
- 12 worker(s): 11.54, 11.90, 12.46 s.

Current **CLI** passes:

- 1 worker(s): 56.59, 56.67, 57.22 s.
- 4 worker(s): 17.54, 18.32, 17.80 s.
- 8 worker(s): 12.51, 13.71, 13.98 s.
- 12 worker(s): 13.03, 19.81, 20.44 s.

Runs were serial, with no simultaneous benchmark/profiler/test process.
Every pass overwrites outputs and decodes the input files again. The bounded
scan cache replays scans only between stages of the same file, never between
benchmark repetitions. JIT warmup and a complete warmup file are excluded;
checksums, Parquet publication and fsync are included. The 50-file benchmark
uses the processing API and process pool; CLI reporting is checked separately.

```bash
# Archive the exact previous source without changing the working checkout.
mkdir -p validation/baseline-0.7
git archive ca63c8c | tar -x -C validation/baseline-0.7

# Use the same 50-file manifest downloaded for MSV000086989.
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 \
PYTHONPATH=validation/baseline-0.7/src .venv/bin/python benchmarks/benchmark_extraction.py \
  --format parquet --repeats 3 --output-directory validation/baseline-parquet \
  --report validation/baseline-timing.json

OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 \
.venv/bin/python benchmarks/benchmark_extraction.py --format parquet --repeats 3 \
  --workers 8 --output-directory validation/current-parquet \
  --report validation/current-timing.json
```

Use workers 1, 4, 8 and 12 in separate runs. Clear generated Numba caches after
changing imported numerical helpers; see the development guide. The baseline
noise audit imports a temporary source copy to avoid contaminating its normal
Numba cache with an aliased module name.

Raw reports: [baseline](extraction-v0.7.0-recheck-parquet.json),
[one worker](extraction-v0.8.0-parquet.json),
[four](extraction-v0.8.0-parquet-4-workers.json),
[eight](extraction-v0.8.0-parquet-8-workers.json),
[twelve](extraction-v0.8.0-parquet-12-workers.json),
[baseline profile](profile-v0.7.0-recheck-parquet.json),
[current profile](profile-v0.8.0-parquet.json).

Full-command reports: [one worker](extraction-cli-v0.8.0-1-workers.json),
[four](extraction-cli-v0.8.0-4-workers.json),
[eight](extraction-cli-v0.8.0-8-workers.json),
[twelve](extraction-cli-v0.8.0-12-workers.json).
Reproduce with the same environment variables as the API benchmark:

```bash
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 \
.venv/bin/python benchmarks/benchmark_cli.py --workers 8 --repeats 3 \
  --output-directory validation/cli-benchmark --report validation/cli-timing.json
```
