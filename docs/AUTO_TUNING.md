# Per-file calibration and polarity switching

Version 0.13 uses calibration policy 3: Gaussian-qualified representatives,
three fixed m/z regions (0–200 / 200–500 / 500+), and mean + 3 sample SD mass-error margins.
Output defaults to Parquet; the calculations apply equally to TSV.

Automatic calibration and ion grouping are **on by default**. Calibration is a
new MassCubeX implementation. It estimates repeatability and chromatography from
the file; it does not calibrate absolute mass accuracy against known compounds.
The MS1 noise default is **fixed 100 Da bins / lowest 10% / minimum 5 / μ + 3σ**.
Neither the noise settings nor the ion-grouping mass tolerance is auto-tuned.

## Broad search, then qualified representatives

An additional streaming pass precedes feature detection. The pilot mass radius is
`max(0.05 Da, m/z × 100 ppm)` on each side, configurable with
`--tuning-search-da` and `--tuning-search-ppm`. Positive, negative and unknown
MS1 streams have independent trackers. A scan of the other polarity never
increments a trace's missed-scan count.

The pilot operates on reader-normalized MS1 without the adaptive filter. Profile
MS1 is centroided using zero baseline/threshold, the configured minimum profile
support (3 samples) and spacing-gap factor (4). Profile rejection still applies
if explicitly requested. For each scan, at most the strongest 128 centroids are
tracked, with deterministic mass-order ties. A candidate is marked as competing
if another centroid in its broad mass interval has at least 20% of its intensity;
this check uses the complete scan, including centroids outside the strongest 128.

Pilot traces need at least 7 observations and allow one missing same-polarity
scan. They use the configured trace time-gap limit. Trackers flush after 512
same-polarity scans or 120 seconds; a peak cut at a block boundary cannot be
accepted without both flanks. Limits are configurable. This bounds pilot trace
memory; it deliberately allows fallback when peaks are too broad or sparsely
sampled to characterize safely.

Each candidate is screened as follows:

1. Require finite values and a positive maximum. Smooth intensity with a
   Gaussian σ of 0.6 observations for this pilot quality check only.
2. Set the baseline to the raw intensity 10th percentile. Require the smoothed
   apex above baseline by more than half the raw maximum.
3. Require one internal maximum with prominence at least 30% of that amplitude.
   Reject plateaus where the raw maximum repeats at least `max(3, n // 3)` times.
4. Require at least 5 raw observations above baseline plus half the amplitude,
   with no competing-mass flag among them.
5. Interpolate both half-height and 10%-height crossings on each side of the
   smoothed peak. Reject missing flanks and total smoothed variation above
   four times the amplitude.
6. Use the median mass of the central observations as the representative mass.
   Its mass error is the 95th percentile of absolute deviations from that
   median. Reject error above 40% of the pilot mass radius.

7. Require symmetric Gaussian similarity **≥0.90**, configurable with
   `--tuning-min-gaussian-similarity`. Use the original intensities after
   subtracting the baseline above and clipping negative excess to zero; do not
   use the smoothed curve. The Gaussian is centered on the raw excess apex;
   its sigma is the bracketed half-height width divided by 2.355. Similarity is
   the Pearson correlation, clamped to [0,1], over the complete pilot trace.
   This is the same symmetric template construction as feature-table Gaussian
   similarity, but the table score uses uncorrected intensities. Calibration is
   deliberately more selective than detection: a low score excludes an anchor,
   not the corresponding detected feature. Shape cannot establish ion purity;
   the competing-mass and central-mass-stability checks still apply.

Keep a bounded pool of the strongest `10 × tuning_max_peaks × region_count`
candidates per polarity: 3,000 by default. At the end of the pass, select up to
100 representatives globally for timing and the file-wide fallback, and up to
100 independently in each m/z region. Each selection uses m/z and RT quartiles
(up to 16 cells), taking the strongest in each cell in rounds. The pilot still
tracks only the top 128 centroids per scan, so sparse high-mass coverage is
reported rather than assumed. `--tuning-max-peaks` controls each selection, not
an overall cap of 100 across all regions. The global top-128 selection can
under-sample weaker mass regions before this balanced representative selection.
The [50-file coverage audit](../benchmarks/results/CALIBRATION_COVERAGE.md)
separates this sampling effect from sparse signal and shape rejection.

Before estimation, screen FWHM and per-peak mass error in ppm independently in
each selection using log-space MAD: retain values within
`max(4.5 × 1.4826 × MAD, log(4))` of the median log value. For mass error,
values below 0.5 ppm share a floor so quantized exports do not reject every
nonzero deviation. Regional screening is independent of global screening:
a genuinely less precise m/z region is not discarded simply because other
regions have lower errors. Require **at least 20 surviving peaks** for each
local MS1 fit. The minimum can be increased, never reduced below 20.

## Fixed mass regions and optional equal-width regions

The default regions are **[0,200), [200,500), [500,∞)**. Exact boundary masses
200 and 500 belong to the higher region. These bounds are shared across files
and polarities; each polarity learns its own tolerances. A target mass selects
the region, but its matching interval may cross a region boundary.

Use `--tuning-mz-breaks 200 400 600` for four regions or supply other strictly
increasing positive boundaries. The last region is unbounded, represented as
`"upper": null` in JSON; infinity is used only inside numerical lookup arrays.
JSON `tuning_mz_breaks` defaults to `[200, 500]`; Python uses the immutable tuple
`(200.0, 500.0)`. An empty JSON list or Python tuple creates one unbounded region.
`tuning_mz_regions` defaults to `null` in JSON (`None` in Python), and an explicit
integer selects the older equal-width layout instead. A CLI boundary list
overrides a region count loaded from JSON; the two CLI options are mutually
exclusive.

`--tuning-mz-regions 5` preserves the five equal-width acquired-range layout.
For that mode, the acquired span is determined as follows. Range provenance is
also recorded with fixed boundaries, but does not move those boundaries.

For each polarity, collect the union of MS1 acquisition windows from mzML
`scanWindow` lower/upper limits (MS:1000501 / MS:1000500), or mzXML
`startMz` / `endMz`. These are acquisition limits; mzXML `lowMz` / `highMz`
are observed extrema. See the [PSI controlled vocabulary](https://github.com/HUPO-PSI/psi-ms-CV)
and [mzXML schema documentation](https://tools.proteomecenter.org/formats/mzXML/mzXML_xmlspy_docs.html).

Use the union of all observed reader-normalized MS1 array extents when explicit
windows are absent, before selecting strong centroids, centroiding profile
arrays or applying the adaptive filter. Do not infer range from representative
masses. Invalid optional window metadata is counted and ignored; if observed
signals extend outside valid metadata, extend the union to include them and
label that source in provenance. Empty scans can still supply valid metadata.
Different scan windows are combined within each polarity; this is not a
separate model for each acquisition regime or RT interval.

Split this span into the requested number of equal-width regions. Regions are
[lower, upper), except the final finite upper bound is included. Outside this
finite range, use the file-wide fallback. With no usable range, this optional
mode retains a global model. `--tuning-mz-regions 1` uses one finite region.
The fixed-boundary default does not require acquisition-range metadata.

## Separate MS1 and MS/MS mass estimates

Each representative contributes one value: its **95th percentile absolute
mass deviation from its own median**, measured using the raw observations in
its central half-height portion. For the collection of these per-peak errors:

```text
candidate Da  = mean(error_Da)  + k × sample_SD(error_Da)
candidate ppm = mean(error_ppm) + k × sample_SD(error_ppm)
MS1 Da        = max(0.0005 Da, candidate Da)
MS1 ppm       = max(5 ppm, candidate ppm)
```

The default **k = 3** is configurable with `--tuning-mass-sigma`; `5` is a wider
alternative. SD uses `ddof=1`, after the robust screen above. These are nonnegative
per-peak Q95 errors, not pooled signed scan residuals. This rule is a safety
margin, **not a 99.7% coverage guarantee or an estimate of absolute mass accuracy**.
Equal peak weights prevent a strong or densely sampled peak from dominating.
Both components are capped at the corresponding broad pilot-search component.
The final half-window is `max(Da, m/z × ppm × 1e−6)`.

MS/MS precursor tolerance has independent evidence requirements: single
selected-ion precursor, known matching polarity, RT inside an accepted peak's
FWHM, and exactly one candidate within the broad mass search. Match against the
union of global and regional representatives, including across region edges;
selecting a region cannot hide a competitor. Isolation-target-only masses do
not calibrate precursor accuracy. For each distinct representative, take Q95
of its absolute precursor-to-MS1-median mass differences. Screen those per-peak
ppm errors with the same log-MAD rule and 0.5 ppm floor, independently per region.
Require **20 distinct peaks with MS2 evidence**, not 20 spectra from one peak.
Then compute mean + k sample SD in Da and ppm, with each component at least the
corresponding local MS1 estimate, and cap at the pilot components.

A region with insufficient evidence does **not** lower its shape threshold or
minimum count. Its fallback is the component-wise maximum of the successful
file-wide estimate and all supported regional estimates. If the file-wide fit
is unsupported, retain configured fallback values as well. MS1 and precursor
fallbacks are independent. An unsupported MS2 region retains at least the local
MS1 estimate unless explicitly overridden. With fewer than 20 global anchors,
unknown polarity or no usable cadence, retain configured mass settings instead
of fitting regions. Default unsupported MS1 is ±max(0.002 Da, 10 ppm); default
unsupported MS2 is ±max(0.01 Da, 20 ppm).

The lookup applies to MS1 trace tracking, MS2 association, preceding-MS1
precursor-intensity lookup, seeded-feature reuse, and raw gap filling. A trace
looks up its running intensity-weighted m/z; MS2 and gap-fill queries look up
precursor m/z. Gap filling inherits the local MS1 pair, with its physical radius
always capped at the local MS1 radius; narrower overrides are honored. Ion
adduct/fragment grouping still uses its explicit 0.002 Da / 10 ppm defaults.

Global scalar settings in each effective config are the **fallback envelope**;
`mass_regions` contains the actual per-region settings. The sidecar records
layout, boundaries, range source, invalid-metadata counts, shape rejections, selected representative
indices and scores, local counts, means, sample SDs, estimated/effective settings
and fallback status. Explicit scalar overrides apply in every automatically
learned region. A nonempty `mass_regions` supplied explicitly through JSON/Python
is preserved and takes precedence over scalars inside its bounds;
use an empty list for automatic regional estimation. Disabling auto-tuning
keeps configured settings, including any explicitly supplied regions.

## Time settings in seconds

Collect actual same-polarity MS1 intervals in a deterministic reservoir of at
most 4,096 values, using separate fixed seeds per polarity. Screen intervals
with the log-MAD rule above. Let `d` be their median, `c` their 95th percentile,
`f10`/`f50` the accepted FWHM 10th percentile/median, and `w90` the accepted
10%-height width 90th percentile.

| Effective setting | Estimate |
| --- | --- |
| Gaussian smoothing σ | `max(0.5 × d, 0.08 × f10)` seconds |
| Minimum separation of chromatographic maxima | `max(d, 0.5 × f10)` seconds |
| Maximum elapsed trace gap | `1.25 × c × (max_gap_scans + 1)` |
| Maximum interval integrated during gap fill | Same as elapsed trace gap |
| Maximum age of a preceding MS1 used for precursor intensity | `2 × c` |
| MS2-to-feature boundary extension | `max(c, 0.15 × f50)` |
| Seed search half-window | `min(tuning_block_seconds, max(3 × c, w90))` |
| Fixed-window seed-reuse tolerance | `max(d, 0.5 × f10)`; primarily relevant to fallback mode |
| Ion-grouping maximum apex-time difference | `max(d, 0.25 × f50)` |

For feature segmentation, interpolate onto a uniform RT grid at approximately
the median interval, capped at four times the original observation count.
Convert seconds to grid samples for Gaussian smoothing and peak separation;
map valleys back to original observations. **Height and trapezoidal area use
original measurements and actual times**, never the interpolated/smoothed curve.
The same-polarity scan-gap count remains 2 by default and is not learned.

When there are fewer than 20 representatives, only convert the legacy smoothing
σ (0.6 scans) and separation (3 scans) to seconds using the observed median
cadence. Peak-width-dependent settings and mass tolerances retain their
configured fallbacks. With no usable cadence, retain the scan-unit settings.
Unknown polarity never supplies learned mass/width settings.

## Refined raw-MS1 gap filling

For a polarity with supported calibration, each unmatched usable MS2 initially
creates its own query. Search the learned time half-window using one EIC
anchored at the exact precursor m/z. Its radius is the smaller of the configured
gap-fill radius and the normal MS1 trace radius, each calculated as
`max(Da, m/z × ppm × 1e−6)`. There is no broad pilot-mass search, recentering,
or padding outside that interval. Keep only the raw slices needed by active queries; for profile data,
centroid again with zero baseline/threshold. Thus a below-noise ion can still
be measured. Do not apply the detection floor or minimum detected-scan count.

Take the maximum centroid intensity inside the anchored interval in each scan
(zero when absent). This is a forced-window measurement, not deconvolution of
unresolved or coisolated ions.

Require positive evidence within `max(ms2_rt_tolerance, 2 × local median MS1
interval)` of the original seed time. No local ion gives
`gapfill_status=no_local_ms1_signal`; no ion anywhere gives `no_ms1_signal`.
Version 0.10 removes the former broad-search / 3× candidate-ranking rule:
overlapping candidate windows could count the same raw centroid twice. The
[blank/QC audit](../benchmarks/results/EIC_REVIEW.md) documents this defect.

For the selected narrow EIC, interpolate and smooth in seconds as above. Find
maxima at the configured prominence and separation, then valleys between them.
Select the basin containing the seed time; ties use apex distance then height.
Within it, trim toward the apex until the smoothed signal reaches the basin's
10th-percentile baseline plus 5% of its amplitude. Map these bounds outward to
original MS1 times. Flat traces retain the basin. If the seed time is outside
the resulting bounds plus `ms2_rt_tolerance`, report `no_local_peak` with blank
quantification. These are forced measurements, not noise-qualified detections.

Integrate original, unfiltered EIC values within those bounds, excluding
intervals above `gapfill_max_gap_seconds`. Do not subtract baseline, extrapolate
to window edges, or interpolate missing intensity for quantification. Store the
unchanged precursor mass, original seed time, actual mass radius and integration
bounds in the TSV. Repeated queries merge only when they recover the same
known polarity, mass within the narrow radius, identical original MS1 apex time
(within 1 microsecond), and each seed time lies within the other query's bounds
plus the MS2 time tolerance. The first retained measurement wins. Different
recovered apices remain separate; no chain of broad time windows can merge them.

If a query exceeds 64 centroids in any slice, 100,000 retained mass values or
10,000 MS1 observations, stop collecting it, keep its MS2 references, and label
`sampling_limit` with blank quantification. No truncated-sample estimate is
published. A query is refined as soon as its time window finishes, and its raw
slices are released immediately. Only the narrow intensity/time trace needed
by enabled ion grouping survives. `ms2_gapfill_peak_retained_values` reports the
maximum simultaneously retained centroid pairs. These limits bound each query,
not the number of concurrent seeds.

With auto-tuning off or insufficient representatives, use the fixed-window
fallback: ±5 seconds, raw EIC radius ±max(0.002 Da, 10 ppm), no local mass or
boundary refinement. Profile fallback takes the maximum raw profile sample.
Seed reuse has its own ±5-second anchor tolerance, independent of the search
window; changing one does not silently change the other. All these fallback
values are configurable. Older 0.4.0 used the broader MS2 mass tolerance for
this measurement. If no ion is found, use the median positive regional noise
threshold from up to the three closest compatible MS1 scans in the query
window. Default imputed height is 0.1× that reference; area remains missing,
scan count is zero, and `ms1_quant_status=imputed_noise` labels it explicitly.
The reference is in `apex_noise_threshold`, although there is no measured apex.
`gapfill_status` retains the original absence reason. `--gapfill-noise-fraction 0`
disables imputation. Missing scans, ambiguous polarity, unsupported noise,
peak-bound mismatch and collection limits remain unquantified. See the output
status reference; imputed height is not an observed signal or evidence for grouping.

## Polarity and reproducibility

- Extraction, pilot fitting, timing, gap filling and ion relations are separated
  by recorded polarity. Negative-mode intervals are measured between negative
  MS1 scans, including time spent on positive acquisition.
- An unknown-polarity MS2 can inherit polarity from an explicit reference to
  the latest retained MS1 of a known polarity; record that inference in JSON.
  It is never inferred from ion intensity or which feature happens to match.
- In mixed-polarity files, unresolved unknown MS2 gets its own seeded row with
  blank MS1 quantification. If all MS1 polarity metadata is absent but MS2
  explicitly contains both polarities, neither may borrow the pooled MS1.
  A single known polarity can still use entirely unknown-polarity MS1 exports.
- Preceding-MS1 precursor-intensity fallbacks are also rechecked after all
  polarities are known; ambiguous sources cannot drive best-MS2 ranking.
- Conflicting positive and negative flags on an mzML spectrum are a file error.
  Missing/incorrect vendor metadata cannot be reconstructed with certainty.

`parameter_tuning` in every JSON sidecar records the policy version, elapsed
time, search settings, rejected-candidate counts, chosen representatives,
per-polarity support/status, estimated/applied parameters and full effective
configs. Original requested settings are retained separately. Resume fingerprints
include this implementation and invalidate earlier outputs after changes.

Use `--no-auto-tune` for fixed settings. Explicit CLI overrides are protected,
even when equal to a built-in default. Non-default JSON/Python values are also
protected. To lock a default-valued JSON/Python setting while keeping the rest
automatic, include its name in `tuning_locked`. Exported default JSON values
remain fallback values unless changed or listed as locked. A legacy scan-unit
override protects its seconds counterpart. For example:

```bash
masscubex extract study -o results --workers 8
masscubex extract study -o results --mz-ppm 10 --ms2-mz-ppm 20
masscubex extract study -o fixed --no-auto-tune
masscubex extract study -o ungrouped --no-group-ions
```

## Limits of the evidence

Strong isolated peaks are useful calibration anchors, but weak-ion mass jitter,
saturation, unusual chromatography, heterogeneous acquisition regimes and export
rounding may differ. The Gaussian cutoff, error margins and floors above are heuristics,
not universal instrument specifications. Calibration does not prove DDA, resolve
DIA mixtures or certify ion identity. Narrow windows can miss signal; broad
windows can mix ions. Inspect the saved evidence and use overrides when needed.

Known-truth tests cover alternating polarities with different cadence, mass
jitter and peak width, both centroid/profile mzML/mzXML, missing polarity,
parent-reference inference, below-noise seeds, nearby interferents and separate
chromatographic peaks. Real-data validation uses 50 positive-mode centroid
mzXML files; it is **not genuine polarity-switching LCMS validation**. See
[the versioned validation report](../benchmarks/results/VALIDATION.md).
