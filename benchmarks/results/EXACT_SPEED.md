# Exact-output optimization — 2026-11-04

This is the round-1 record. See [round 2](EXACT_SPEED_ROUND2.md) for the latest
timings; the measurements below are retained unchanged.

The local extraction workflow is **2.93× faster across two files**, with
byte-identical Parquet tables. Scientific algorithms and default parameters are
unchanged from the frozen source snapshot taken before this optimization.

Single-worker medians of two measured runs per file, using default settings and
compressed Parquet. JIT/full-path warmup and profiling are excluded from timings.

| Input | Features, before = after | Before (s) | After (s) | Speedup |
| --- | ---: | ---: | ---: | ---: |
| NIST_POS_Samp_14-03.mzXML | 93,731 | 46.56 | 15.29 | 3.04× |
| QC_Mix_01.mzXML | 65,652 | 28.34 | 10.24 | 2.77× |
| Both files | 159,383 | 74.90 | 25.53 | 2.93× |

## Where the time went

Stage totals below sum the per-file medians; small transaction/sorting overhead
is outside these stages. Small differences in untouched stages are timing variation.

| Stage, both files | Before (s) | After (s) |
| --- | ---: | ---: |
| Calibration | 1.02 | 1.06 |
| MS1 extraction, including splitting/finalization | 20.27 | 12.92 |
| Ion grouping | 51.23 | 9.12 |
| MS2 association and seeded gap filling | 0.65 | 0.58 |
| Table writing | 1.59 | 1.72 |

Ion grouping improved **5.62×** and extraction **1.57×**. Extraction is now the
largest stage. The original-function finalization/splitting remains a substantial
cost; this optimization does not replace those scientific calculations.

## Changes that preserve the algorithms

- Build the sorted mass-search index once in float64, an exact embedding of the
  float32 masses, instead of repeatedly converting the entire index for each
  query. Actual mass membership still uses the original float32 expression.
- Compile candidate eligibility/RT filtering with Numba, without fastmath.
  Preserve strict boundaries, priority ordering and nearest-RT ties. Assignment
  generations preserve the per-anchor eligibility snapshot without copying a
  whole mask each time. Empty searches return immediately.
- For EICs already sampled on the anchor's ordered time grid, skip redundant
  time intersections; retain the identical NumPy Pearson calculation.
- Reuse finalized measurements for unchanged traces, stop repeating an unchanged
  splitting result, and skip unconsumed top-average/asymmetry calculations.
  Avoid deep-copying full parents immediately before the original subset/finalize
  operations rebuild child buffers. Parent/child buffer independence is tested.
- Use scalar bisections for adjacent acquired MS1 time lookup, preserving the
  minute/second roundtrip, nearest-time comparison and ties.

The pinned MassCube 1.2.22 scientific source bodies remain unchanged. Tolerances,
adaptive noise, splitting rules, quantification, ion catalog and MS2 selection
are unchanged. No feature filters or lower-precision calculations were introduced.

## Equality checks and scope

Every measured optimized run has the same complete Parquet SHA256 as its
baseline: all 25 columns, row order, measurements, grouping relationships and
MS2 references. Scientific JSON sidecars also match; only source fingerprints
and elapsed/stage timings are excluded from that comparison. All 4,654 MS2 scans
remain assigned. The example TSV was also regenerated and compared byte-for-byte.

120 distinct focused tests passed, including original-function parity,
full-precision summaries, buffer independence, RT boundaries, assignment
snapshots and correlation equivalence. This was not a full release-suite run.

Only these **two positive-mode centroided mzXML files** from MSV000086989 were
processed. The timings are warm-cache results on this Mac, not a guarantee for
larger files, profile data, other instruments or multiple workers. No additional
real-data cohort was run. Temporary benchmark tables and binary profiles are removed after
saving the report; the example TSV remains.

[Measurements, hashes, settings and environment](exact-speed-2026-11-04.json) ·
[Before profile](exact-speed-before-2026-11-04.txt) ·
[After profile](exact-speed-after-2026-11-04.txt) ·
[Benchmark command](../benchmark_exact_speed.py)

To compare another saved baseline, run the command with exactly two input paths,
`--package-path /path/to/baseline-src`, `--output-directory /tmp/mcx-before` and
`--report /tmp/mcx-before.json`. Then run the current implementation with a new
output/report location and `--compare /tmp/mcx-before`. Both runs use two measured
repeats by default. `--profile-first` profiles the separate warmup only.
