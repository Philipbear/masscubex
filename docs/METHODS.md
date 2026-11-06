# MassCubeX methods

MassCubeX retains the core tracing, splitting, quantification, MS2 assignment and
ion-grouping algorithms from MassCube 1.2.22 at commit
`356748ede691862b9a87f67812db5e728a703b9d`. Its default additions are adaptive
MS1 noise, automatic mass/timing parameter estimation, profile MS1 centroiding
and MS2-seeded gap filling. These additions can change the detected/seeded rows;
default output is not claimed to be identical to upstream. The selected adduct
catalog and repository-oriented input/output are also MassCubeX adaptations.
See [the comparison table](UPSTREAM_COMPARISON.md).

## Input and scan preprocessing

Before scientific processing, [acquisition screening](ACQUISITION.md) identifies
confirmed DIA using metadata and repeated window evidence. Such files use
MS1-only extraction by default: MS1 isotope/adduct grouping is retained, while
DDA-style MS/MS association, MS2 seeding and fragment-based grouping are disabled.
Unknown modes remain labelled unknown and use the configured ordinary processing.
No DIA deconvolution is performed; explicit skip/error policies remain available.

The streaming readers support mzML/mzXML, gzip, float32/float64 binary arrays,
zlib and optional MS-Numpress. They validate XML, lengths, encodings, finite
values, scan times and polarity metadata; no network vocabulary lookup is used.
Explicit optical spectrum types are skipped before MS-level or array validation;
their wavelength arrays never enter mass-feature processing.
Explicit zero-length mzML arrays may have empty binary payloads even when zlib
compression is declared. Nonempty compressed payloads still require successful
decompression and matching array lengths; corrupted data is not treated as empty.
Profile MS1 is centroided by default; `--profile-mode reject` disables this extension.

Scientific scan processing uses float32 m/z/intensity arrays, matching upstream.
For MS1, estimate adaptive thresholds on those arrays, retain intensities strictly
above threshold and masses strictly between 0 and 100,000 Da, then apply original
`centroid_signals`. This sorts by m/z, joins consecutive neighbors separated by
less than 0.002 Da, sums each group's intensities and calculates its
intensity-weighted mass in float32. This scan-level merging is distinct from
chromatographic feature mass, which is the apex mass. The displayed threshold
for a merged signal is the maximum threshold of its contributing centroids;
that diagnostic does not introduce another detection gate.

Retention times are converted to seconds for I/O and upstream minutes for
scientific finalization. The upstream acquisition range is 0–10,000 minutes.
Polarity streams are processed independently: a scan of another polarity does
not advance the missing-scan counter. Unknown polarity is not silently assigned
to a known stream. Ion-mobility-resolved MS1 inputs are not collapsed implicitly.

## Adaptive MS1 noise

Each scan is divided into fixed 100 Da bins anchored at zero: [0,100), [100,200),
etc. Boundary masses enter the higher bin. For N positive intensities in a bin,
select the lowest `max(5, ceil(0.10*N))` values, limited by availability. With at
least five positives, calculate their mean μ and population standard deviation
σ (`ddof=0`), and use `max(min_intensity, μ + 3σ)`. All points in that bin share
the threshold; there is no leave-one-out exclusion and no borrowing from adjacent bins.
With fewer than five positives, use only the configured absolute floor (zero by default).

The compiled estimator partitions once per bin. Shifted moments avoid cancellation
on large baselines. It does not smooth, subtract or rescale retained intensities.
A low-tail standard deviation is not a calibrated false-positive probability.
No chromatographic shape gate is added to compensate for a permissive noise threshold.

`--noise none --min-intensity VALUE` provides a fixed-floor control. Historical
sliding windows and clipped estimators remain explicit parameter alternatives;
they are not default. `noise_sigma` remains fixed at 3.

`--noise-estimator consistency` selects the revised population/evidence model
described in [the estimator design and experiments](NOISE_ESTIMATOR_DESIGN.md).
Its approved μ + 3σ decisions are unchanged by the optimized implementation:
compiled searches preserve strict mass boundaries and ties; stable intensity
grouping preserves fixed-diameter bands; background reductions retain NumPy's
original float64 accumulation order. A histogram kernel reproduces explicit
NumPy bin edges exactly. Filtered training surveys are reused and released after
their last temporal neighbor. Persistent centroids and decisions use temporary
disk-backed arrays, while MassCube tracing, splitting, quantification and ion
grouping remain unchanged. This mode is opt-in and does not replace the default
low-tail estimator. Raw-gap-fill extraction and imputation retain their existing
raw-spectrum rules, including below-threshold searches.

## Trace building

The optimized array-pool tracker implements original `detect_features` rules:

1. Start a trace for each first-scan signal in mass order.
2. For each later scan, process active traces in stable descending intensity
   order from the previous scan. Missed observations have intensity zero.
3. Query the closest mass to the last observed mass, within the applicable
   MS1 tolerance (absolute difference strictly less than the radius).
   An exact midpoint between neighboring masses is unassigned. If the closest
   signal was already claimed, do not try a second mass.
4. Append a matching signal and reset the gap counter; otherwise append zero
   at the previous mass. Finish after the eleventh consecutive miss (ten allowed).
   There is no default elapsed-time gap cutoff.
5. Start new traces for unclaimed signals. Expired-trace removal order and stable
   ties follow upstream. Retain detected features with at least five nonzero MS1
   observations by default. `--min-scans 1` retains short traces.

Positive/negative streams use the same rules separately. Automatic calibration
is on by default and supplies each stream's regional MS1 tolerance:
`max(Da floor, last observed m/z × ppm × 1e−6)`. It changes the search radius,
not the last-observed-mass target, nearest-mass decision or missing-scan rule.
With insufficient file-wide evidence, retain the configured fallback, originally
±0.01 Da with no ppm term. `--no-auto-tune` uses configured values directly.
Calibration never adds an elapsed-time cutoff to the ten-missing-scan rule;
`--max-gap-seconds` remains an explicit optional limit.

## Splitting and measurements

Upstream rules finalize and split each trace. Trim leading and trailing zeros.
Splitting requires at least five nonzero observations. The configurable
`min_scans` threshold also applies to detected features after splitting.
A trace-noise score greater than 3 prevents splitting. The score
zeros observations below 5% of the trace maximum, zero-pads the ends, and computes
half the total absolute first difference divided by the maximum, minus one,
with the original arithmetic; spans shorter than five samples have score zero.

For trimmed spans shorter than five samples (including interior zeros), a direct
equivalent avoids temporary upstream feature objects: no split is possible,
Gaussian/noise scores are zero, and a singleton has zero area. Longer spans are
represented by intervals until final measurements are needed. This avoids
deep-copying parent objects and calculating discarded parent summaries.

The interval adapter retains raw valleys, shared valley samples, trimming,
first-apex ties and both passes. Measurement specializations reuse the known
apex and omit duplicate copies/reductions, retaining the original noise/Gaussian
arithmetic and exact area expression below. Gaussian smoothing uses the
installed SciPy C kernel with cached coefficients; a compiled peak selector falls
back to the original SciPy path for ambiguous distance-priority ties. A public
SciPy fallback covers changes to its private kernel API. No fastmath or extra
shape/intensity filters are used.

The raw-median gate uses native order statistics. An even float32 median adds
the two middle float32 observations before division, retaining NumPy's rounding
and overflow behavior. This changes execution only; it does not change which
maxima qualify. Full-precision comparisons cover both contiguous and strided
float32/float64 traces.

For eligible traces, zero-pad and Gaussian-smooth (sigma 0.6 scans). Find peaks
with prominence 3% of the smoothed maximum and distance three scans. Retain only
maxima whose raw intensity exceeds the raw median. Split at raw minima between
adjacent selected maxima; the valley sample is shared by the two intervals.
Apply two passes, re-finalizing each subset as upstream.

Feature m/z and RT come from the first maximum. Height is the original integer
conversion of that float32 maximum. Area is the original integer conversion of
`trapz(float32 intensities, RT in minutes) * 60`. Gaussian similarity uses
upstream's apex-centered Gaussian with FWHM-derived sigma and Pearson similarity.
It is a diagnostic, not a detection/splitting gate. Original source bodies and
attribution are in `upstream.py`; parameters are bound without mutating module globals.

Final rows are ordered by m/z, preserving upstream's stable parent-trace ordering.
Native IDs and numeric scan references are unchanged. Outputs round m/z to five
decimals, times in seconds to two, intensity values to one and correlations to four.
These formatting rules are separate from upstream's measurement arithmetic.

## MS/MS handling

Use upstream fixed MS2 intensity filtering and 0.002 Da adjacent-centroid merging.
The default floor is 10,000 for an `orbitrap` header keyword, otherwise 500,
matching upstream's first-300-lines instrument selection; `--ms2-min-intensity`
overrides it. The relative threshold defaults to zero. MS2 filtering is not adaptive.
The summed fragment score is computed on the processed float32 spectrum using
upstream's sum. Ties select the first acquisition, without an apex-distance tie-break.

Associate using the applicable precursor mass tolerance and the feature's
adjacent-MS1 scan edges. These association edges extend one acquired MS1 scan outside each
trimmed trace endpoint, and are represented in float32 minutes as upstream.
Among overlapping candidates choose the largest peak height, not the nearest
mass or apex time. Auto-tuning supplies precursor mass tolerances when supported
and a boundary-time allowance from measured scan cadence/peak widths. Without
calibration, the defaults remain strict ±0.1 Da and no extra RT allowance, as
upstream. There is no isolation-width cutoff. The first precursor entry is used,
as upstream. Unmatched usable precursors are passed to the enabled seeding
extension; scans still ineligible for assignment remain in JSON diagnostics.

Every assigned feature exports all scan numbers/native IDs/source spectrum
indices plus the best scan number. Sparse sidecar fallbacks preserve the best
native ID/index when scan numbers are absent or duplicated. Original scan
indices remain zero-based positions in the complete source spectrum sequence.

## Ion grouping

Follow upstream strongest-first grouping with the selected ion catalog. Each
new positive/negative anchor starts as [M+H]+/[M-H]-. Enumerate isotopes, allowed
ions from the anchor neutral mass, then its best-MS2 products. Grouped features
do not start competing anchors. Within an anchor, upstream's candidate-mask
snapshot permits later assignments to replace earlier ones.

Use mass error <0.01 Da, apex RT difference <3 seconds, and partner EICs within
±0.005 Da over the anchor scan interval. Check only the nearest-RT non-isotope
candidate. Pearson correlation on shared positive observations must exceed 0.9;
fewer than three returns 1, as upstream. Constant longer traces fail naturally.
Isotope centers are `mz + 1 + (i-1)*1.003355` for i=1..3. Child-ion isotopes
additionally require height below 1.5 times their parent. No extra overlap or
minimum-five-scan rule applies. The isotope/source-fragment flags follow the
single-file upstream assignment semantics. See [the catalog and labels](ION_GROUPING.md).

The optimized grouping index caches only immutable mass/RT/polarity candidates.
An anchor's initial assignment mask is applied to all its cached queries in one
call, preserving the greedy order, current-anchor snapshot and later replacement
behavior. Dynamic child-isotope checks still use current assignments. Sparse
storage and target tables have explicit size limits; broad or unsupported cases
fall back to the scalar query path. The common NumPy 1.x path batches queries;
NumPy 2.x keeps its general scalar-promotion path. These are execution changes,
not different grouping criteria.

Sparse IDs and offsets use signed 32-bit storage under the existing index size
bound. Packed storage grows within that bound in one discovery pass, rather
than repeating every binary search. Relationship lists are allocated only for
features with evidence. Processed grouping spectra retain their original
float32 values in a temporary memory map; EIC results remain float64. An anchor's
scan-time join is reused for its partner queries. Pearson specialization retains
NumPy's float64 row means, matrix dot product, covariance scaling, sequential
standard-deviation divisions and clipping; it is not an approximate formula.

## Extensions enabled by default

These stages extend the pinned upstream algorithms:

- Automatic calibration (`--no-auto-tune` to disable): region-stratified pilot peaks, Gaussian similarity ≥0.80,
  at least 20 representatives, robust error screening and mean + 3 SD margins;
  regions 0–200, 200–500 and 500+. Mass/timing estimates and fallbacks are recorded.
  [Full calibration details](AUTO_TUNING.md).
- Profile centroiding (`--profile-mode reject` to disable): split profile mass signals at local valleys,
  require at least three samples, reject unsupported zero plateaus and acquisition
  spacing gaps, estimate centroids using numerical first moments and preserve
  apex intensities. This is additional MS1 mass detection, followed by the default
  upstream preprocessing/tracing. It does not imply profile-MS2 deconvolution.
- MS2 seeding (`--no-ms2-seeded-features` to disable): cluster unmatched precursors into explicitly labelled
  rows and measure raw MS1 EICs. Cap the gap-fill radius at the normal MS1 mass
  tolerance; a narrower configured gap-fill radius is honored. With supported calibration, use local chromatographic bounds;
  otherwise use the configured ±5-second window. If no ion is found, impute height
  as 0.1× regional noise, leave area missing and label the value. Seed provenance
  remains internal/JSON; redundant seed m/z and RT columns are not exported.

Other explicit parameter overrides can change settings. JSON `algorithm_policy`
records whether these extensions were enabled; `config` records all overrides.

## Repository execution and verification

The file transaction owns validation, extraction, optional extensions, grouping,
atomic Parquet/TSV publication and a JSON completion record. `--resume` checks
input identity, configuration, implementation/dependency fingerprints and output
checksums. File workers process independent files; file errors are recorded in
persistent summary/error manifests and do not stop other inputs. Global setup or
storage failures can still abort a batch. Temporary grouping storage closes on
success or error. A budgeted decoded-scan cache supports optional repeat passes.

[Parity checks](UPSTREAM_COMPARISON.md) explicitly disable these three extensions
and compare complete measurements and MS2 references with original functions on
fixed-floor and adaptive-noise inputs.
Counts alone are insufficient. Historical throughput and count reports retain
their original implementation hashes and do not describe these new defaults.
Alignment, library annotation and DIA deconvolution remain outside this command.
