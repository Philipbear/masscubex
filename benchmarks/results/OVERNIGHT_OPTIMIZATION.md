# Extraction speed and RAM — overnight 2026-11-04/05

Completed locally. Recommended for the measured workload: **four processes,
recycle setting 8**, at **1,921 files/hour** and **2.853 GiB peak summed RSS**.
Three processes with the same setting use **2.191 GiB** at **1,494 files/hour**.
All 100 unique default tables and all 1,800 successful sustained jobs match the
frozen pre-round output. Scientific defaults are unchanged.

## Workload and exactness

Fifty additional public files were downloaded from MSV000086989, giving **100
unique inputs**: 54 biological samples, 42 blanks and four QCs. They total
1,319,964,953 bytes (1.32 GB), all positive centroid QExactive mzXML. Default
extraction produces **10,039,342 features**, with 442,656,376 bytes of Parquet
output. This cohort does not benchmark profile detection or all repository types.

Every table is **byte-identical** to the frozen pre-round implementation at one,
two, three and four processes and two threads. Scientific JSON sidecars also
match after excluding only elapsed/stage timings and the changed implementation
fingerprint. No scientific defaults, peak decisions, scan assignments, grouping
rules, output schema or precision changed in this round.

## Single-file improvements

The continuous serial 100-file run fell from **817.51 s to 644.31 s**: 21.2% less
wall time, or **26.9% more files/hour**. Peak total RSS fell from **1.021 GiB to
0.953 GiB** (6.6%). The same original and additional manifests are used in the
same order, with identical settings and bookkeeping. The five-file repeated
check separately fell from 36.40 s to 29.02 s (sum of two-repeat file medians).

- Batch native raw-median gating/peak selection, preserving float32 median
  arithmetic and the existing SciPy fallback for ambiguous distance ties.
- Reuse known apex state in Gaussian/noise diagnostics and specialize the
  exact NumPy area/correlation expressions. Reductions, BLAS order and integer
  quantification stay the same; no fastmath or approximate scores.
- Build bounded sparse grouping candidates in one packed pass, use int32
  indices within their existing size limit, and filter each anchor snapshot
  once. Keep greedy consumption and overwrite order unchanged.
- Reuse EIC time joins and store already-float32 processed spectra without
  redundant float64 storage. Numerical query outputs retain their previous dtype.
- Prime all production JIT signatures/layouts and initialize provenance/writer
  imports before the first large feature graph. A fresh-cache audit found no
  extra numerical signatures compiled by the checked real centroid file.

The process coordinator no longer imports a second numerical runtime. Its
measured startup RSS fell from **161.5 MiB to 28.0 MiB**. The writer initializer
uses a discarded one-row memory buffer and publishes no extra file.

## Fresh-worker throughput matrix

Apple M2 Max, 32 GiB RAM, macOS 26.6.2, Python 3.10.12, NumPy 1.26.0, Numba
0.60.0, SciPy 1.13.1, PyArrow 17.0.0. The normal defaults include a 64 MiB raw
scan cache, calibration, adaptive MS1 noise, MS2 seeding and ion grouping.
These timings exclude per-worker full-path warmup. They include extraction,
Parquet writing, checksum/provenance comparisons and deletion of checked outputs.
Benchmarks/tests were not run concurrently with each other.

| Execution | 100 files | Files/hour | Peak summed RSS |
| --- | ---: | ---: | ---: |
| 1 process | 10.74 min | 559 | 0.95 GiB |
| 2 processes | 5.56 min | 1,079 | 1.53 GiB |
| 3 processes | 3.79 min | 1,585 | 2.20 GiB |
| 4 processes | 2.96 min | 2,030 | 2.62 GiB |
| 2 threads | 12.93 min | 464 | 1.16 GiB |

Processes outperform threads for this workload; two threads are slower than one
process. These are short, fresh-worker peaks. They underestimate long-running
allocator retention, so four processes are not recommended solely from this table.

## Sustained CLI runs

The best measured throughput under the chosen 3 GiB budget uses **four process
workers** and `--worker-recycle-files 8`. **400 jobs** (four cycles of the 100-file
cohort) completed in **749.73 s (12.50 min)** at **1,921 files/hour**, with
**2.853 GiB peak summed pipeline RSS**. All 400 outputs match the frozen baseline.
The run created 13 successive pools, confirming that restarts release retained
allocator memory rather than just delaying the peak to the next pool.

Use this configuration for this measured workload on the Mac:

```bash
masscubex extract --manifest files.txt -o features --workers 4 --resume
# Explicit equivalent of the new recycling default:
masscubex extract --manifest files.txt -o features --workers 4 --resume \
  --worker-recycle-files 8
```

The pool drains after `8 * workers` extraction attempts (32 with four workers),
then restarts only if work remains. The setting bounds the pool's total attempts,
not the exact count received by each individual process. Resume hits do not
count, no file is interrupted, and failed files keep their own error record.
`--worker-recycle-files 0` disables retirement. This policy applies to process
pools with more than one worker; one-worker/direct and thread execution do not
recycle. Worker count still defaults to one.

A separate **1,000-job** run used three processes with setting 32: **2,281.34 s**,
**1,578 files/hour**, **2.535 GiB peak RSS**, 11 pools and identical outputs.
Four workers with setting 8 provide **21.7% more throughput** than that
longer-interval three-worker run, using **326 MiB more peak RSS**.

| Processes | Recycle setting | Verified jobs | Wall time, including restarts | Files/hour | Peak pipeline RSS | 2 million similar files |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| 3 | 8 | 400 | 16.06 min | 1,494 | 2.191 GiB | 55.8 days |
| 3 | 32 | 1,000 | 38.02 min | 1,578 | 2.535 GiB | 52.8 days |
| 4 | 8 | 400 | 12.50 min | 1,921 | 2.853 GiB | 43.4 days |

The three-process/setting-8 run completed 17 pool lifetimes. It uses 13.6% less
peak RSS than setting 32 with 5.3% lower throughput, which supports the shorter
recycling default. Four processes fit the selected budget and give the highest
measured sustained throughput; three provide a lower-memory option. All **1,800
jobs** in these successful sustained runs match the reference.

Four workers exceeded the **3 GiB benchmark guard** in sustained trials with
recycling settings 128, 32 and 16, despite the lower short-batch peak. The
successful setting-8 run is why the default was reduced. Guard stops are
controlled cancellation of this benchmark's owned processes, not scientific
output mismatches. They do not establish that four workers are unsafe for the
whole Mac; they establish that those intervals miss the chosen footprint target.

Repeated cycles reuse the same raw cohort; they are not additional independent
inputs. The initializer's numerical/writer code is unchanged in the final CLI;
only its recycling default/help changed after measuring explicit settings.

At the measured sustained four-process rate, **two million similar files need
about 43.4 days** on this Mac running continuously, excluding download time.
The tested three-process/setting-32 run projects 52.8 days. These are cohort-based
capacity estimates, not measurements of two million files.

## Remaining stage costs

Serial wall-time contributions across all 100 files:

| Stage | Seconds | Share of measured wall time |
| --- | ---: | ---: |
| MS1 feature extraction | 254.07 | 39.4% |
| Ion grouping | 228.16 | 35.4% |
| Parquet writing | 68.18 | 10.6% |
| Calibration | 48.89 | 7.6% |
| MS2 association | 20.32 | 3.2% |
| MS2-seeded gap filling | 10.82 | 1.7% |
| Other orchestration/checking | 13.87 | 2.2% |

Extraction and grouping still dominate. These are stage wall times, not CPU-user
time; their sum should not be added to inclusive profiler timings.

## Original MassCube control

The fresh five-file [upstream comparison](UPSTREAM_SPEED.md) is **16.56×** faster
in aggregate on identical adaptive-noise workloads (12.75–23.74× per file), with
zero mismatched core rows/MS2 references/grouping labels. On the fixed Orbitrap
floor control it is **5.03×**, also with zero mismatched rows. Upstream lacks the
adaptive estimator, calibration and seeding extensions, so those comparisons
share preprocessing and explicitly disable unsupported extensions. Selected
single-file grouping is enabled in both; this is not stock upstream's default.

The original timer excludes serialization; MassCubeX includes its complete
Parquet transaction. See the linked comparison for the pinned commit, adapter,
parameters and timing boundaries. Do not use the matched speedup as a guarantee
for every instrument, preset or file.

## Reproducibility and scope

Numerical implementation SHA-256:
`fca3d76939fc32766261a865252ba33710f8f90f5abbe89ff5790748268d8a8e`.
The [before-round source archive](../baselines/README.md) retains the exact
baseline without generated Python/Numba caches. Its fingerprint is
`da4d382e1609a392fa09765a31c96c851e1de21ccab6eb8defa7281d81167a68`.
CLI/runtime/harness hashes are recorded separately in the sustained reports.

RSS sums count shared library pages separately in each process. They are
conservative process-footprint measurements, not unique private memory or total
Mac memory use. Only owned/known worker PIDs are inspected; unrelated apps and
system caches are untouched. The 3 GiB guard belongs to the benchmark; normal
CLI runs control workers/recycling and do not enforce a total-RAM cap. Per-file
feature/relationship memory still scales with feature count.

Downloads, first-ever compilation and independent test/audit runs are excluded
from throughput projections. Warm filesystem/cache results on these relatively
small centroid files should not be extrapolated to larger/profile files or
network storage without measurement. Sustained tests include repeated JIT-cache
loading at pool restarts. PyArrow reports restricted sysctl capability probes
in this sandbox; all comparisons use the same environment.

## Validation and cleanup

**265 targeted checks pass**, including full-precision numerical comparisons,
extreme-float32/NaN median behavior, sparse-query limits/ties, worker recycling,
resume behavior, failed-file continuation and atomic output failure. Ruff and
`git diff --check` pass. The refreshed example TSV remains byte-identical
(SHA-256 `15a83773ef2b2bdbb7e0dcc119de4ca503589e61592c5e6caf641eaa864e91f9`)
with current implementation provenance.

Older MS2/seeding/polarity fixture failures were reproduced in the earlier
baseline/candidate audit. This round also reproduced the raw-reader duplicate-
merging assertion against the frozen baseline; reader parsing bodies are
unchanged. See [the code audit](../../docs/CODE_AUDIT.md) and
[validation record](overnight-validation-2026-11-05.json) for the exact scope.
No scientific-policy changes were made to alter those assertions. This is not a
clean full release-suite claim.

Before the round, 32,262,353 bytes of obsolete wheel/smoke/build results were
removed. Fresh-cache validation additionally removed 157 stale Numba cache files
(9,866,135 bytes). Raw inputs, compact reports and the example table are retained.
Checked repeated feature tables are deleted as each job finishes. Final cleanup
removed another **151,423,389 bytes** of task-owned source snapshots, temporary
outputs/logs and binary profiles: **193,551,877 bytes removed across this round**.
Text profiles, compact measurements, raw inputs, the example and the 83 KB
baseline source archive remain. All benchmark pools have shut down; no background
extraction job is left running.

Changes are local. No commit, push, PR, release or full platform/wheel matrix.

[100-input checksum/URL manifest](overnight-input-manifest-100.json) ·
[Before 100-file report](overnight-baseline-continuous100-2026-11-05.json) ·
[Final serial report](overnight-final-100-w1-2026-11-05.json) ·
[Five-file exact comparisons](overnight-final-five-2026-11-05.json) ·
[Final profile](overnight-final-five-2026-11-05.profile.txt) ·
[JIT audit](overnight-warmup-audit-final-2026-11-05.json) ·
[Machine-readable summary](overnight-summary-2026-11-05.json) ·
[Validation record](overnight-validation-2026-11-05.json) ·
[Run inventory](overnight-run-inventory-2026-11-05.json) ·
[Three-process sustained run](overnight-runtime-3-recycle8-2026-11-05.json) ·
[Four-process sustained run](overnight-runtime-4-recycle8-2026-11-05.json) ·
[Cleanup inventory](overnight-cleanup-2026-11-04.json)


Reproduce after downloading the two manifests' immutable inputs:

```bash
python benchmarks/benchmark_resources.py --directory validation/MSV000086989 \
  --manifests download_manifest.json additional_manifest.json --workers 1 \
  --reference benchmarks/results/overnight-reference-100-2026-11-04.json \
  --output-directory /tmp/mcx-resources --report /tmp/mcx-resources.json
python benchmarks/benchmark_cli_resources.py --directory validation/MSV000086989 \
  --manifests download_manifest.json additional_manifest.json --workers 4 \
  --cycles 4 --worker-recycle-files 8 --memory-budget-mib 3072 \
  --reference benchmarks/results/overnight-reference-100-2026-11-04.json \
  --output-directory /tmp/mcx-cli --report /tmp/mcx-cli.json
```

Scientific-sidecar hashes include local input paths/mtimes. On another checkout,
first generate its local baseline reference from the source archive, then supply
that report to `--reference`; input and table checksums remain portable. Install
the `dev` dependencies for the RSS monitor. The harness default was updated to 8
after the explicit-setting trials; recorded harness/CLI hashes describe the code
at measurement time, while the runtime and numerical fingerprints are unchanged.

Publication follow-up: the complete local source and installed-wheel suites each pass 637 tests.
Older expectations were reconciled against the pinned upstream contract without changing
the numerical implementation. See the [publication record](publication-validation-2026-11-05.json)
and [code audit](../../docs/CODE_AUDIT.md). Statements above about deferred publication
describe the original measurement round.
