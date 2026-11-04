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

The increase is **17 locally fitted regions**, or 34.7% relative to the old 49.
Coverage remains limited in the end regions: 42/50 files in each still use the
file-wide fallback. More sampled centroids cannot guarantee enough isolated,
well-shaped chromatographic peaks. These counts measure available calibration
evidence, not absolute mass accuracy or chemical identity.

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

The paired benchmark uses three warm serial passes per implementation, alternating
version order across passes, with the same files, dependencies, options and host.
JIT warmup, startup, downloads and independent audits are excluded. Baseline
outputs must reproduce the archived v0.13 Parquet checksums; repeated outputs
must be identical within each implementation, and new Parquet rows must agree
with the independently audited TSV rows. See
`benchmarks/benchmark_regional_sampling.py` for the complete procedure.

The source and wheel, full output tables, test XML and individual pass reports
remain under `validation/regional-sampling-v0.14/` rather than Git. Reproduce by
copying the v0.13 `src/masscubex` package into its `baseline/src/` directory,
then running `benchmarks/validate_ms2_coverage.py --workers 4` followed by
`benchmarks/benchmark_regional_sampling.py` serially. The benchmark also uses the
archived v0.13 tables in `validation/calibration-coverage-v0.13/three/` and the
baseline fingerprint in `baseline/identity.json`.
