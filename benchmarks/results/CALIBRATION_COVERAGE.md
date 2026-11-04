# Calibration coverage and fixed mass regions — MassCubeX 0.13.0

The default is now **[0,200), [200,500), [500,∞)**. On the same 50 files, this
layout puts more detected features in a region with sufficient local calibration
evidence than the proposed four-region layout. It does not solve the shortage
of accepted anchors at either end of the mass range. Gaussian similarity ≥0.90,
the 20-peak evidence minimum, and mean + 3 sample SD remain unchanged.

## Why only 57/250 of the old regions had local fits

The cohort contains 6 blanks, 4 QCs and 40 biological samples from MSV000086989,
all positive centroid mzXML. The old five equal-width regions were approximately
100–380, 380–660, 660–940, 940–1220 and 1220–1500. Local MS1 support was
49, 8, 0, 0 and 0 files, respectively.

There are three contributors:

1. **High-mass signal is sparse.** Centroids at m/z ≥600 account for 7.07% of raw
   centroid observations and 4.38% of summed MS1 intensity. These are observation
   and intensity shares, not counts of distinct chemical compounds.
2. **The pilot selects only the strongest 128 centroids per scan before tracing.**
   Only 3.29% of those selected observations are at m/z ≥600. Weaker peaks can
   drop out of this selection in successive scans; later representative
   balancing cannot recover a trace that was never assembled.
3. **Shape and support screening remove additional candidates.** At the default
   Gaussian threshold, the median accepted counts for the four fixed regions
   are 9, 39, 11 and 1. The dominant high-mass rejection is fewer than seven
   retained observations, followed by shape/support failures. A short pilot
   trace is not proof that the underlying signal is noise: sampling can break it.

Consequently, this is also a low-mass calibration problem. The 0–200 region has
many observations, but fewer than 20 accepted anchors in every file with the
current sampler and quality rules. In the original five-region layout, even
removing the Gaussian gate leaves fewer than ten accepted anchors in each of
the upper three regions of every file. Lowering that gate alone cannot solve it.

## Three versus four regions, with identical quality rules

Both layouts were run through extraction on all 50 files. Feature coverage uses
the **same 743,635 archived v0.12 detected feature masses** as its denominator,
so changes in output feature count do not bias the comparison.

| Measure | 0–200 / 200–500 / 500+ | 0–200 / 200–400 / 400–600 / 600+ |
| --- | ---: | ---: |
| Files with a local MS1 fit, by region | 0 / 49 / 0 | 0 / 49 / 2 / 0 |
| Median accepted anchors, by region | 9 / 45 / 5 | 9 / 39 / 11 / 1 |
| Locally fitted file-region pairs | 49 / 150 | 51 / 200 |
| Reference features lying in a supported region | 428,886 (57.7%) | 361,008 (48.5%) |
| Held-out peak-error summaries within the effective tolerance | 3,127 / 3,136 | 3,127 / 3,136 |
| Output rows | 747,913 | 747,905 |
| Associated MS2 scans | 125,266 | 125,266 |
| Warm serial extraction time, all 50 files | 45.70 s | 46.16 s |

Three regions are the preferred starting default **for this cohort at the
current quality settings**: splitting the middle mass range creates an extra region that
rarely has its own fit, without improving the held-out check. Broader-region
coverage alone does not establish better mass accuracy. Blank_00 has only 19
accepted global anchors and retains configured mass settings throughout.
With three regions, local precursor fits are supported in 0 / 32 / 0 files;
global precursor fits are supported in 45 files.

The held-out check orders Gaussian-qualified candidates by RT/mass, uses five
interleaved folds, and fits on the other four folds. Each of 3,136 candidates
is held out once. It compares that peak's Q95 absolute mass-error summary with
the effective MS1 radius at its mass. It does **not** measure raw-observation
coverage, weak-ion errors, absolute mass accuracy or compound identity. The
5 ppm floor dominates many estimated windows, limiting this comparison's power.

Timing includes extraction stages and output writing, after warmup, with the
layouts interleaved file by file on the same M2 Max. The independent audits and
holdout calculations are excluded. This single paired pass does not establish
a speed difference between the layouts.

## Would 0.8 or ten peaks help?

Replaying the same candidates gives the following numbers of supported
file-region pairs. The global evidence minimum stays at 20 in every case;
the ten-peak column is an exploratory regional availability count, not a
production setting or a validated ten-peak mass fit.

| Gaussian minimum | Regional minimum | Three regions | Four regions | Old five equal-width regions |
| --- | ---: | ---: | ---: | ---: |
| 0.90 | 20 | 49/150 | 51/200 | 57/250 |
| 0.90 | 10 | 81/150 | 104/200 | 87/250 |
| 0.80 | 20 | 60/150 | 76/200 | 85/250 |
| 0.80 | 10 | 118/150 | 144/200 | 97/250 |

More supported regions do not by themselves show that a relaxed model is
accurate. A ten-peak fit also has a less stable SD estimate. These experiments
do not justify silently lowering the default quality threshold or evidence
minimum. Gaussian 0.8 remains available explicitly through the existing CLI.

An additional pilot-sampling experiment kept Gaussian ≥0.90 and 20 anchors:

| Sampling per scan | Maximum centroids | Supported files in four fixed regions |
| --- | ---: | ---: |
| Strongest globally, current default | 128 | 0 / 49 / 2 / 0 |
| Strongest 32 in each region | 128 | 0 / 45 / 14 / 0 |
| Strongest globally | 512 | 16 / 50 / 34 / 0 |
| Strongest 128 in each region | 512 | 8 / 50 / 37 / 1 |

Sampling more centroids recovers additional anchors, particularly at low and
intermediate mass, while the highest region remains poorly supported. With
regional sampling at 128 per bin, the median accepted count in 600+ rises from
1 to 5.5, still below 20. These are exploratory pilot runs; **production retains
the global 128-centroid sampler**. A future sampling change needs its own
accuracy and throughput validation rather than treating fit counts as the goal.

## Implementation and validation

- New `--tuning-mz-breaks 200 500` default; use `200 400 600` to request four
  regions. Exact boundary masses enter the higher region. An unbounded final
  upper limit is JSON `null`. `--tuning-mz-regions 5` retains the previous
  equal-width acquired-range option.
- Sparse regions retain the documented file-wide tolerance envelope. MS1,
  precursor matching and gap filling use the same region lookup; gap filling
  never gets a wider physical MS1 window. Noise filtering, segmentation and
  feature-table schema are unchanged.
- **689 source tests and 689 installed-wheel tests pass**, plus lint. New cases
  cover fixed/custom boundaries, an unbounded tail, metadata independence,
  JSON round trips, invalid layouts and legacy CLI behavior.
  All [12 CI configurations](https://github.com/Philipbear/masscubex/actions/runs/41779617877)
  also pass source and installed-wheel tests: Linux, macOS and Windows with
  Python 3.10–3.13, on implementation commit `697b0f4`.
- The independent Pyteomics audit accounts for **125,266/125,266 MS2 scans**:
  120,964 match detected features and 4,302 have seeded features. It also checks
  157,324 raw seeded-quantification sample points. There are 4,143 measured seeds,
  122 flagged noise-height imputations, 35 with no local peak and 2 with no local
  MS1 ion. Imputed areas remain missing.
- All **747,913 typed Parquet and TSV rows agree** for the default layout:
  743,611 detected + 4,302 seeded. This is 23 fewer total rows than v0.12;
  a small count difference does not imply that individual rows are identical.

Evidence: [layout comparison and sampling summary](calibration-layouts-v0.13.0.json),
[independent raw scan/quantification audit](MSV000086989-ms2-v0.13.0.json),
[release fingerprint and tests](release-v0.13.0.json), and
[complete calibration methods](../../docs/AUTO_TUNING.md).

Reproduce from the existing 50-file manifest with
`benchmarks/audit_calibration_coverage.py --workers 4`, then
`benchmarks/validate_ms2_coverage.py --workers 4`, followed by
`benchmarks/compare_calibration_layouts.py` serially. The last command also needs
the archived v0.12 Parquet tables under `validation/regional-tuning-v0.12/parquet/`.
Pilot candidates, full tables, installed wheel and test XML are retained locally
under `validation/calibration-coverage-v0.13/` and
`validation/MSV000086989/v0.13.0-ms2/`. The original diagnostic pilot was recorded
before the layout change; its implementation hash is retained separately in the
comparison report, which verifies its projected support against production fits.
