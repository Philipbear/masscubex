# Development calibration policy 5

Local comparison on 2026-11-04. This is a development change, not a release.
No PR, push, wheel build, full test suite or CI workflow was run for this initial comparison.
Subsequent 0.15.0 release validation added a √2 MS1 fallback transfer margin
after detecting trace fragmentation in a polarity-switching fixture; independent
regional fits are unchanged. See [current validation](VALIDATION.md) and
[current performance](PERFORMANCE_CURRENT.md) for the final published version.
Earlier policy-4 release/benchmark records remain historical snapshots.

Implementation SHA256: `3ac1728c1140deec4b60ff2af30e83829422920825152d5394dea4c695d4e289`.

## Method changes

- Pilot radius: ±max(0.025 Da, 50 ppm), previously 0.05 Da / 100 ppm.
- At least five **observed MS1 scans**, previously seven. Tracker zeros do not
  count. At least three central mass measurements remain required, previously five.
- Gaussian similarity ≥0.80, previously 0.90. Median selected-anchor similarity
  across this cohort remains approximately 0.91–0.92 in each region.
- Remove separate relative-contrast, multiple-maximum and total-variation vetoes.
  Keep plateau rejection, competing-mass checks, central mass stability and a
  measured rise/fall around the apex. Flat traces still fail.
- Central mass measurements are restricted to the main apex's half-height
  interval. Secondary peaks elsewhere cannot supply central observations.
- A five-scan peak can use raw 10%-height crossings when smoothing raises its
  endpoints above that level. Missing-scan zeros cannot manufacture the observed
  flanks required below half height. Completely clipped peaks remain excluded.
- Screen retained candidate pools **before** the final intensity/stratification
  quota. Mass errors above `median + max(4.5 × 1.4826 × MAD, 1 ppm)` are removed.
  This is a one-sided robust screen, distinct from the unchanged mean + 3 sample
  SD tolerance estimate. FWHM retains the existing log-MAD screen. Regional pools
  are screened independently. The precursor-error screen uses the same mass rule.
- For an unsupported region, transfer the largest fitted **ppm** component and
  use a 0.0005 Da floor. Do not copy a high-mass absolute Da estimate to low mass.
  Supported local fits retain their own Da/ppm components. If file-wide evidence
  is insufficient, configured fallbacks remain. Explicit overrides are preserved.
- The minimum independent evidence remains **20 anchors**. The per-scan sampling
  cap remains 128 per region; the final representative cap remains 100 per region.

Detailed rules: [AUTO_TUNING.md](../../docs/AUTO_TUNING.md).

## Same 50-file cohort

The cohort comprises 40 biological files, six blanks and four QCs from
MSV000086989. These are positive-mode centroid mzXML files.

| Region | Independent fits: before → after | Median anchors: before → after | Fewest anchors after |
|---|---:|---:|---:|
| 0–200 | 8/50 → 50/50 | 15 → 100 | 45 |
| 200–500 | 50/50 → 50/50 | 72.5 → 100 | 100 |
| 500+ | 8/50 → 50/50 | 11.5 → 100 | 27 |

Overall: 66/150 → 150/150 independent regional MS1 fits. Median retained
quality-qualified candidates before the robust screen/final quota are
133.5 / 555.5 / 143.0; the reported final medians reach the existing 100-anchor cap.

All low and middle regions now resolve to ±5 ppm over the acquired m/z range
(approximately 100+). High regions retain their independent estimates: median
Da component 0.002367, range 0.001555–0.005975; ppm median 5, range 5–6.391.
The actual radius remains `max(Da, m/z × ppm × 1e-6)`.

In `NIST_POS_Samp_14-03`, the low region changes from 16 anchors and a fallback
of 0.0042545 Da / 5.6609 ppm to 100 anchors and a local fit of 0.0005 Da / 5 ppm.
At m/z 150, the radius therefore changes from ±0.0042545 Da (28.36 ppm) to
±0.0007500 Da (5 ppm). Its high-region Da component falls to 0.0024352 Da.

## Focused validation

90 calibration/regional-sampling tests pass, including five-real-scan support,
missing-flank rejection, competing masses, outlier selection order, independent
NumPy/SciPy measurement agreement, polarity separation, regional boundaries,
protected overrides and prevention of cross-region Da transfer. Changed Python
files pass Ruff and the diff whitespace check.

Three extraction smoke checks also complete:

| File | Total feature rows: before → after | Assigned MS2 scans |
|---|---:|---:|
| Blank_00 | 11,814 → 11,787 | 2,268 |
| QC_Mix_01 | 12,623 → 12,598 | 2,058 |
| NIST_POS_Samp_14-03 | 16,766 → 16,486 | 2,596 |

All 6,922 MS2 scans are reported assigned, with zero unassigned. These smoke
checks do not replace the independent raw-scan audit used for release validation.
Timing settings also change when calibration anchors change; row-count differences
cannot be attributed solely to mass tolerances.

The calibration-only pass took 20.55 seconds across 50 files, excluding warmup.
This is one diagnostic pass, not a paired speed benchmark or full extraction time.
Greater coverage does not establish absolute mass accuracy or chemical purity;
this comparison contains no new independent holdout validation and no real
profile-mode or negative-mode cohort.

Reproduce with `.venv/bin/python benchmarks/audit_calibration_policy.py`.
Full per-file diagnostics and the three Parquet outputs are saved locally under
`validation/calibration-policy5/`; raw inputs and policy-4 sidecars are required.
