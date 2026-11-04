# Regional pilot sampling — MassCubeX 0.14.0

Calibration now selects the strongest **128 centroids per m/z region per MS1
scan**, using the fixed default regions **[0,200), [200,500), [500,∞)**. Gaussian
similarity ≥0.90, at least 20 accepted anchors, and mean + 3 sample SD are unchanged.

## Coverage on the same 50 raw files

The comparison uses the same 6 blanks, 4 QCs and 40 biological samples from
MSV000086989, all positive centroid mzXML. The baseline is the unchanged v0.13
implementation at commit `1d8487b`, with the same three regions and quality rules.

| Evidence | Global strongest 128 | Strongest 128 per region |
| --- | ---: | ---: |
| Files with independent MS1 fit: 0–200 | 0/50 | 8/50 |
| Files with independent MS1 fit: 200–500 | 49/50 | 50/50 |
| Files with independent MS1 fit: 500+ | 0/50 | 8/50 |
| Independent MS1 file-region fits | 49/150 | 66/150 |
| Median accepted anchors by region | 9 / 45 / 5 | 15 / 72.5 / 11.5 |
| Files with a supported global MS1 fit | 49/50 | 50/50 |
| Files with independent precursor fit, by region | 0 / 32 / 0 | 0 / 47 / 0 |
| Files with a supported global precursor fit | 45/50 | 49/50 |
| Reference detected features within locally fitted regions | 428,886 (57.7%) | 482,532 (64.9%) |

The increase is **17 locally fitted regions**, or 34.7% relative to the old 49.
Coverage remains limited in the end regions: 42/50 files in each still use the
file-wide fallback. More sampled centroids cannot guarantee enough isolated,
well-shaped chromatographic peaks. These counts measure available calibration
evidence, not absolute mass accuracy or chemical identity.
Feature coverage uses the same **743,611 archived v0.13 detected feature masses**
for both versions, so differences in output feature counts do not change the
denominator. The end regions remain under-supported despite the larger sample.

## Paired processing time

Three warm serial passes over all 50 files, on the existing M2 Max / Python
3.10.12 environment, produced these results:

| Time for all 50 files | Global strongest 128 | Strongest 128 per region |
| --- | ---: | ---: |
| Full pass 1 | 49.62 s | 53.25 s |
| Full pass 2 | 48.20 s | 53.43 s |
| Full pass 3 | 47.82 s | 55.41 s |
| **Median full workflow** | **48.20 s** | **53.43 s** |
| Calibration stage, median | 14.98 s | 19.87 s |
| Feature extraction stage, median | 15.74 s | 15.88 s |
| Ion grouping stage, median | 5.28 s | 5.19 s |
| MS2 association stage, median | 1.46 s | 1.47 s |
| MS2 gap filling stage, median | 3.25 s | 3.56 s |
| Table writing stage, median | 6.92 s | 7.05 s |

Full processing is **10.9% slower**, adding 5.23 seconds per 50 files, or
approximately **0.105 seconds per file**. Calibration is 32.6% slower and is now
the largest stage; its 4.89-second increase accounts for most of the additional
time. Stage medians are computed independently and need not sum to the full-pass
median, which also includes pipeline overhead. These are warm serial processing
times, not download times or a concurrent-worker scaling test.

The new median is **1.069 seconds per file**. A direct extrapolation to two
million files of this size and type is **24.7 days with one worker**, before
downloads and other repository overhead. This is not a prediction for larger,
longer, profile-mode or otherwise heterogeneous public files. Parallel throughput
also depends on memory, storage and available CPU resources.

## Sampling and memory details

- `--tuning-peak-limit 128` now limits each region, allowing at most 384 centroids
  per scan with the default three regions. Sparse regions keep all available
  points and do not lend spare capacity to other regions. Positive, negative
  and unknown streams remain separate.
- Partial intensity selection avoids a full intensity sort. Cutoff ties retain
  the earlier mass-order indices, and selected indices remain in mass order.
  Exact boundary masses enter the higher region. Profile data are centroided
  before this selection; the adaptive noise filter is not applied to the pilot.
- Competing-ion checks still inspect the entire centroided scan, including
  unselected ions and neighbors across region boundaries.
- Qualified complete candidates have separate per-region pools, each bounded
  at `10 × tuning_max_peaks` (1,000 by default). This prevents a crowded region
  from displacing another region's candidates. The maximum total pool remains
  3,000 per polarity. None of these 50 files reaches a regional pool limit;
  the largest retained pool total is 153, so this safeguard does not explain
  the observed coverage changes in this cohort.
- Fixed regions need no new reading pass. The optional equal-width layout
  (`--tuning-mz-regions`) first discovers the complete acquired range separately
  by polarity, then uses the same fixed boundaries for sampling and mass fitting.
  This uses bounded scan replay when available, or another streaming read after
  cache overflow. No single early scan determines the boundaries for the run.
- Sidecars record selected-centroid and qualified/retained-candidate counts per
  region. The feature-table schema and output precision are unchanged.

## Quantification and implementation checks

**708 source tests and 708 isolated installed-wheel tests pass**, plus lint.
Tests compare partial selection with an independent full-sort oracle and cover
ties, empty/sparse regions, exact boundaries, unselected competing ions across
boundaries, crowded candidate pools, expanding acquisition ranges, polarity
switching, bounded/disabled cache replay and retained MS2 fragments.

The independent raw-data audit preserves **125,266/125,266 MS2 scans** exactly
once: 120,953 match detected features and 4,313 have seeded features. It checks
175,287 raw seeded-quantification sample points. The output has **748,062 rows**:
743,749 detected and 4,313 seeded, versus 747,913 rows before (+149, or 0.020%).
More calibration anchors also change learned time settings, so changes in
feature counts are not attributable solely to mass tolerances.

| Group | Previously detected | Newly detected | Previously seeded | Newly seeded |
| --- | ---: | ---: | ---: | ---: |
| Blank, 6 files | 74,972 | 74,949 | 180 | 193 |
| QC, 4 files | 56,994 | 57,660 | 141 | 141 |
| Biological, 40 files | 611,645 | 611,140 | 3,981 | 3,979 |

The small net count change masks changes in individual files: Blank_00 gains
310 total rows and QC_Mix_01 gains 326, while NIST_POS_Samp_17-02 loses 392.
The noise filter and peak-shape guard are unchanged; these results do not prove
that individual feature identities are unchanged or chemically validated.

The paired benchmark uses three warm serial passes per implementation, alternating
version order across passes, with the same files, dependencies, options and host.
JIT warmup, startup, downloads and independent audits are excluded. Baseline
outputs must reproduce the archived v0.13 Parquet checksums; repeated outputs
must be identical within each implementation, and new Parquet rows must agree
with the independently audited TSV rows. All those checks passed, including
**748,062 typed Parquet/TSV row comparisons**. See
`benchmarks/benchmark_regional_sampling.py` for the complete procedure.

Evidence: [paired coverage and timing measurements](regional-sampling-v0.14.0.json),
[independent raw-data audit](MSV000086989-ms2-v0.14.0.json),
[release fingerprint and tests](release-v0.14.0.json), and
[complete calibration algorithm](../../docs/AUTO_TUNING.md).

The source and wheel, full output tables, test XML and individual pass reports
remain under `validation/regional-sampling-v0.14/` rather than Git. Reproduce by
copying the v0.13 `src/masscubex` package into its `baseline/src/` directory,
then running `benchmarks/validate_ms2_coverage.py --workers 4` followed by
`benchmarks/benchmark_regional_sampling.py` serially. The benchmark also uses the
archived v0.13 tables in `validation/calibration-coverage-v0.13/three/` and the
baseline fingerprint in `baseline/identity.json`.
