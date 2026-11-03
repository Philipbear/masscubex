# Per-file calibration and polarity switching (0.5.1)

Automatic calibration and ion grouping are **on by default**. Calibration is a
new MassCubeX implementation. It estimates repeatability and chromatography from
the file; it does not calibrate absolute mass accuracy against known compounds.
The **50 Da / lowest 10% / minimum 5 / μ + 3σ MS1 noise rule is unchanged**.
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

At most the strongest `10 × tuning_max_peaks` candidates (1,000 by default) are
retained per polarity. Divide that pool by m/z and RT quartiles into up to 16
cells; take the strongest from each cell in rounds, up to 100 representatives.
This favors intensity while reducing domination by a single m/z or RT region.
The pool itself is intensity-selected, so complete coverage is not guaranteed.

Before estimation, screen FWHM and mass error in ppm using log-space MAD:
retain values within `max(4.5 × 1.4826 × MAD, log(4))` of the median log value.
For mass error, values below 0.5 ppm share a floor to avoid treating all
nonzero deviations as outliers in quantized exports. **At least 20 accepted
representative peaks per known polarity** are required after screening. The
minimum can be increased, but not reduced below 20.

## Separate MS1 and MS/MS mass estimates

Each representative gets one vote, regardless of intensity or number of scans.
If `Q95` is the 95th percentile across accepted representatives:

```text
MS1 Da floor = max(0.0005 Da, 2 × Q95(per-peak mass error in Da))
MS1 ppm     = max(5 ppm,      2 × Q95(per-peak mass error in ppm))
```

Each component is capped at its pilot-search component. Trace matching and raw
gap-fill extraction use `max(Da floor, m/z × ppm × 1e−6)` as their radius.
The gap-fill parameters are separate settings, so either may be overridden.

MS/MS precursor tolerance requires its own evidence: a single selected-ion
precursor, known matching polarity and RT within a representative's FWHM, with
exactly one candidate in the broad mass search. Isolation-target-only masses
do not calibrate precursor accuracy. For each distinct representative, take
Q95 of its absolute MS2-precursor-to-MS1-median mass differences. Screen these
errors in ppm using the same log-MAD rule and 0.5 ppm floor. Require **at least
20 distinct MS1 representatives with accepted MS2 evidence**, not 20 repeated
spectra from one ion. Then:

```text
MS2 Da floor = max(learned MS1 Da floor, 1.5 × Q95(precursor error in Da))
MS2 ppm      = max(learned MS1 ppm,      1.5 × Q95(precursor error in ppm))
```

Again cap at the pilot components. If MS1 estimation succeeds but precursor
evidence is insufficient, keep the configured MS2 fallback, normally
±max(0.01 Da, 20 ppm). Do not substitute the MS1 estimate without MS2 evidence.

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
creates its own query. Search the learned time half-window, with a broad mass
radius equal to the larger of the pilot-search and gap-fill radii. Retain two
additional narrow radii of context beyond each search edge, so a center near
the edge still has its complete narrow EIC. Candidate centers must remain
inside the original search. Keep only the raw slices needed by active queries; for profile data,
centroid again with zero baseline/threshold. Thus a below-noise ion can still
be measured. Do not apply the detection floor or minimum detected-scan count.

Sort the collected masses and partition into clusters whose full diameter is
at most twice the narrow gap-fill radius at the first mass; this prevents
unlimited chaining. Use each cluster's median mass as a candidate center and
reconstruct its exact narrow EIC from the retained slices, taking the maximum
centroid intensity in each scan (zero when absent).

Require positive evidence within `max(ms2_rt_tolerance, 2 × local median MS1
interval)` of the original seed time. Prefer centers inside the independently
resolved MS2 mass tolerance. If none qualify, consider the broad-search
candidates. Rank by local maximum intensity, then distance to precursor mass.
If the best is less than three times the second-best eligible candidate,
report `ambiguous_mass` with blank quantification. This factor is an explicit
heuristic; MS2 triggering alone does not prove which nearby ion was isolated.

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
refined mass, original seed precursor/time, actual mass radius and integration
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
this measurement. Zero means compatible scans existed with no signal; blank
means unavailable, ambiguous, or unsupported. See the output status reference.

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
rounding may differ. The multipliers/floors above are conservative heuristics,
not universal instrument specifications. Calibration does not prove DDA, resolve
DIA mixtures or certify ion identity. Narrow windows can miss signal; broad
windows can mix ions. Inspect the saved evidence and use overrides when needed.

Known-truth tests cover alternating polarities with different cadence, mass
jitter and peak width, both centroid/profile mzML/mzXML, missing polarity,
parent-reference inference, below-noise seeds, nearby interferents and separate
chromatographic peaks. Real-data validation uses 50 positive-mode centroid
mzXML files; it is **not genuine polarity-switching LCMS validation**. See
[the versioned validation report](../benchmarks/results/VALIDATION.md).
