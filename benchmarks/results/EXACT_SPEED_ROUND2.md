# Exact-output optimization, round 2 — 2026-11-04

This is the round-2 record. The [later five-file upstream comparison](UPSTREAM_SPEED.md)
contains the current timings and direct original-MassCube benchmark.

This round reduces processing time by **17.2%**, or **1.21× additional speedup**,
against the already optimized implementation from round 1. Both feature tables
are byte-for-byte identical, including all grouping and MS2 information.

The same two positive-mode centroided mzXML files were measured serially with
default settings and compressed Parquet. Values are medians of **three warm
runs per file and version**; JIT/full-path warmup and profiling are excluded.

| File | Features, unchanged | Before this round (s) | After (s) | Speedup |
| --- | ---: | ---: | ---: | ---: |
| NIST_POS_Samp_14-03.mzXML | 93,731 | 15.31 | 12.71 | 1.20× |
| QC_Mix_01.mzXML | 65,652 | 10.34 | 8.52 | 1.21× |
| Both files | 159,383 | 25.64 | 21.23 | 1.21× |

Compared with the original pre-round-1 measurements, the cumulative speedup is
**3.53×**. Both rounds retain exactly the same Parquet hashes. These are warm
single-worker results on this Mac, not a throughput guarantee for other files,
profile input or multiple workers.

## Stage totals across both files

| Stage | Before (s) | After (s) |
| --- | ---: | ---: |
| Calibration | 1.06 | 1.01 |
| MS1 extraction, including splitting/finalization | 12.72 | 10.94 |
| Ion grouping | 9.14 | 7.55 |
| MS2 association and seeded gap filling | 0.70 | 0.56 |
| Table writing | 1.67 | 1.24 |

Extraction remains the largest stage. Per-stage medians do not necessarily sum
to the total median, and small transaction/sorting costs are outside these
stages. Calibration and MS2 algorithms were untouched; changes in their small
timings should not be interpreted as targeted speedups.

## What changed

- Short **trimmed spans under five samples**, including any interior zeros,
  directly apply upstream's existing length guards and measurement arithmetic.
  They cannot split and have zero Gaussian/noise scores; singletons also have
  zero area. This removes temporary upstream objects and unused shape arrays,
  while preserving float32 casts, first-apex ties, NumPy area arithmetic and
  array-based RT conversion. Longer spans retain the original-function path.
- Grouping combines binary searches, mass checks and eligibility filtering in
  one compiled query. NumPy resolves the subtraction/comparison dtypes; explicit
  rounding preserves its float32/float64 behavior and strict boundaries. There
  is no fastmath. Unsupported or nonfinite query cases retain the NumPy path.
  Per-polarity adduct rules are cached without changing their order.
- Parquet batching resolves each Arrow column name once and uses bounded
  iterator slices. Rounding, schema, row-group size and compression are unchanged.

Only `features.py`, `grouping.py` and `output.py` changed in the scientific package
during this round. Configuration and pinned MassCube 1.2.22 source bodies remain
unchanged.

## Verification and timing variability

All six measured optimized runs match the complete baseline Parquet bytes and
scientific JSON metadata. Only timing fields and the source fingerprint are
excluded from the JSON comparison. All 4,654 MS2 scans remain assigned. The
example TSV was also regenerated and checked byte-for-byte.

**151 focused tests passed**, covering reference summaries, numeric boundaries,
mass-promotion combinations, grouping, segmentation and output formatting.
No full release suite, additional real-data cohort, push or release was performed.

The first optimized QC timing was 11.48 seconds and the second 8.35 seconds, so a
third repetition was added for **both** implementations; the third optimized QC
run was 8.52 seconds. All measurements, including the slow run, are retained in
the JSON report. Temporary benchmark outputs and binary profiles are removed
after retaining these results.

[Complete measurements, hashes and settings](exact-speed-round2-2026-11-04.json) ·
[Before profile](exact-speed-round2-before-2026-11-04.txt) ·
[After profile](exact-speed-round2-after-2026-11-04.txt) ·
[Round 1](EXACT_SPEED.md) · [Benchmark command](../benchmark_exact_speed.py)
