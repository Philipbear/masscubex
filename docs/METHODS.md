# MassCubeX: implementation and methods

This document describes the implemented MassCubeX 0.13.0 processing path,
including the low-tail MS1 estimator, MS2-seeded features and persistent batch reports. It is a description of
the code and its defaults, not a claim of a newly published scientific method.
For a quick start, see the [README](../README.md). Commands, Python usage and
the output column reference are in [USAGE.md](USAGE.md).
Automatic calibration, its constants and polarity-switching rules are detailed in
[AUTO_TUNING.md](AUTO_TUNING.md). For every numbered ion and custom-ion definitions, see [ION_GROUPING.md](ION_GROUPING.md).

## 1. Relationship to MassCube

The starting point was a Git clone of
[Philipbear/masscube at `356748ede691862b9a87f67812db5e728a703b9d`](https://github.com/Philipbear/masscube/tree/356748ede691862b9a87f67812db5e728a703b9d).
The source history, Apache 2.0 license and upstream attribution are retained.
The runtime has been rewritten as `masscubex`; it is not API-compatible with
the old `masscube` package. The comparisons below refer to that exact commit.

| Component | Relationship to the cloned MassCube implementation |
| --- | --- |
| Input scope | Replaced the vendor/cache-oriented readers with streaming mzML/mzXML readers, including gzip. The cloned commit's reader explicitly rejected mzXML; this fork implements it directly. |
| MS1 noise | Replaced fixed absolute cutoffs with regional m/z background. Version 0.10 defaults to fixed 100 Da bins and the lowest-intensity 10% of positive samples, expanded to at least 5, followed by mean + 3 standard deviations of that subset (§3). Historical sliding windows, clipped estimates and an absolute floor remain selectable. |
| Profile mass detection | Added explicit profile peak detection and numerical first-moment centroiding (§4). This differs from upstream `centroid_signals`, which merged adjacent masses separated by less than 0.002 Da and summed their intensities. |
| Trace building | Retained the concept of matching successive MS1 peaks into ion traces; replaced Python per-centroid objects/list updates with compiled array pools, changed matching to a running intensity-weighted mass, and added a time-gap limit (§5). |
| Chromatographic splitting | Reimplemented Gaussian smoothing, prominence and area boundaries. Relative prominence defaults to 0.1 rather than upstream's 0.03. Version 0.11 restores upstream's noise score as a trigger for a new per-peak shape guard, adds Gaussian-similarity reporting, and permits supported multiple peaks in noisy traces (§6). The above-median peak requirement is not retained. |
| Quantification | Retained peak height as the default and trapezoidal area in intensity × seconds. Added explicit selection and kept both measurements without integer truncation; detected-feature m/z is intensity-weighted rather than the upstream apex mass (§6). MS2 seeds record the original precursor and may refine their mass from raw MS1 (§7). |
| MS/MS | Retained summed MS2 intensity as the default best-spectrum criterion. Added all associated native identifiers, best-scan identifiers and fragment/precursor intensities, optional ranking methods, MS2-seeded features with forced MS1 quantification, and per-scan unresolved reasons. Isolation-width exclusion is off by default (§7). |
| Automatic settings | Added broad-search pilot fitting with at least 20 qualified representatives per polarity, independent MS1/MS2 tolerance estimates, physical-time smoothing and local gap-fill boundaries ([details](AUTO_TUNING.md)). |
| Ion relationships | Replaced greedy forced group membership with default-on, evidence-labelled pairwise hypotheses; corrected signed-charge mass arithmetic and undersupported correlation; added user-selected forms and water losses (§8). |
| Repository execution | Added lazy discovery, bounded independent file workers, compressed Parquet by default (optional TSV), atomic table/JSON publication, checksum-verified resume, persistent error reports, and targeted retries (§9–10). |
| Removed scope | No vendor SDKs, raw cache formats, spectral-library annotation, isotope assignment, normalization, statistics, plotting, retention-time correction, or integrated dataset alignment. Alignment is a separate future process. |

In the audited [upstream parameters](https://github.com/Philipbear/masscube/blob/356748ede691862b9a87f67812db5e728a703b9d/src/masscube/params.py),
the actual MS1 absolute default is 1,000 and the MS2 default is 500. Comments
recommend MS1 thresholds of 30,000 for Orbitrap and 1,000 for QTOF. These should
not be confused with universal instrument noise levels. MassCubeX does not
choose an MS1 threshold from the instrument name. The historical comparisons
against 500 and 30,000 are controlled comparisons, not a claim that both were
upstream MS1 defaults.

Implementation map:

| Source | Responsibility |
| --- | --- |
| [readers.py](../src/masscubex/readers.py) | XML, binary arrays, metadata, scan identifiers |
| [scan_cache.py](../src/masscubex/scan_cache.py) | Budgeted decoded-scan replay and streaming fallback |
| [chromatograms.py](../src/masscubex/chromatograms.py) | Compiled peak segmentation and measurements |
| [peak_shapes.py](../src/masscubex/peak_shapes.py) | Gaussian similarity, trace roughness and evidence for splits in noisy traces |
| [calibration.py](../src/masscubex/calibration.py) | Compiled quality screening; reproducible accepted-candidate measurements |
| [noise.py](../src/masscubex/noise.py) | Adaptive MS1 background and threshold |
| [centroid.py](../src/masscubex/centroid.py) | Profile MS1 mass detection |
| [features.py](../src/masscubex/features.py) | Trace pools, segmentation, quantification |
| [ms2.py](../src/masscubex/ms2.py) | Precursor-intensity fallback, association, representative selection |
| [tuning.py](../src/masscubex/tuning.py) | Per-polarity pilot selection, robust parameter estimates, fallbacks |
| [gapfill_refinement.py](../src/masscubex/gapfill_refinement.py) | Narrow raw EICs, local boundaries and refined-seed deduplication |
| [gapfill.py](../src/masscubex/gapfill.py) | New MS2-seeded features and streaming raw-MS1 EIC measurements |
| [adducts.py](../src/masscubex/adducts.py), [grouping.py](../src/masscubex/grouping.py) | Ion chemistry and relationship evidence |
| [pipeline.py](../src/masscubex/pipeline.py) | Per-file orchestration and atomic outputs |
| [cli.py](../src/masscubex/cli.py), [reporting.py](../src/masscubex/reporting.py) | Discovery, scheduling, continuation, reports and retries |
| [config.py](../src/masscubex/config.py) | Validated, immutable processing configuration |

## 2. Scan reading and normalization

The default sequence starts with a separate streaming calibration pass
([full method](AUTO_TUNING.md)), then: read a scan → normalize arrays/metadata → estimate
MS1 background → retain centroids or detect profile masses → extend ion traces
→ split and quantify completed traces → associate MS2 → create unmatched
MS2 seeds and, when needed, stream MS1 again for targeted gap filling → propose
ion relationships → publish one Parquet/TSV table and one JSON sidecar.

### Formats and streaming

File extensions are case-insensitive: `.mzML`, `.mzXML`, and either with `.gz`.
The reader uses `lxml.iterparse`. Calibration, extraction and raw-MS1 gap filling
remain separate logical passes. Version 0.7.0 replays decoded scans when their
conservatively charged array/metadata storage fits `scan_cache_mb` (64 MiB per
worker by default). On overflow, the cache is discarded and later passes stream
the input again; `--scan-cache-mb 0` always streams. A cache is usable only after
the first pass reaches valid XML EOF. MS2 metadata is copied per pass so polarity
inference, precursor intensity estimation and fragment pruning cannot contaminate
subsequent passes. Raw MS1 arrays stay unfiltered for gap filling. Cached arrays
are released before ion grouping. Sidecar `scan_cache` records the budget, peak
charged bytes and physical-read/replay counts. XML elements and preceding siblings
are released after use; no byte-offset index is constructed. External entity resolution, DTD loading and network access
are disabled. There is no runtime controlled-vocabulary download.

mzML supports plain and indexed wrappers, referenceable parameter groups,
little-endian float32/float64 arrays, uncompressed/zlib data, and optional
MS-Numpress linear, positive-integer and short-logged-float encodings, including
their zlib combinations. Numpress requires the optional `pynumpress` dependency.
Unrelated binary array types are ignored. Required arrays, binary types,
compression terms and declared per-spectrum lengths are checked.

mzXML supports interleaved mass/intensity pairs, 32/64-bit precision, network
(big) or little endian, and uncompressed/zlib data. Per-scan `peaksCount` is
checked. A stack handles nested scans. MS1 is emitted when its `<peaks>` ends,
before nested MS2, preserving acquisition order and the preceding-MS1 fallback.
An mzXML processing-level centroid flag is used only when a scan has no flag.

Only MS1 and optionally MS2 are processed. Higher levels and chromatograms are
not feature inputs. Their spectrum positions still count toward the source
spectrum index. Unsupported encodings fail that file; there is no silent
fallback to a different numerical interpretation.

### Numerical and metadata conventions

- Arrays become contiguous float64. Masses must be positive and finite;
  intensities must be finite. Array dimensions and lengths must agree.
- Centroid intensities must be nonnegative; zero-intensity entries are removed.
  Profile zeros are retained as boundaries and negative samples are clipped
  to zero. Duplicate m/z coordinates are summed after a stable mass sort when
  needed; this applies to both modes. There is no additional near-mass merging
  of already centroided peaks.
- mzML time is converted from declared seconds, minutes, hours or milliseconds
  to seconds. mzXML day/hour/minute/second durations are converted to seconds.
  Missing/unsupported MS1 units, negative/nonfinite times, and nonincreasing
  MS1 times within a polarity are errors.
- Positive, negative and unknown polarities have separate trace trackers.
  Interleaving polarities does not count as a missing scan of the other polarity.
- Explicit profile MS1 is centroided automatically, or rejected with
  `profile_mode="reject"`. Missing mode is an error unless
  `assume_centroid=True`; that setting never overrides explicit profile mode.
- Recognized ion-mobility-resolved mzML MS1 coordinates are rejected. The
  current implementation expects one ordinary mass spectrum per MS1 time and
  polarity, without an additional mobility dimension.

`input_points` counts normalized MS1 points before adaptive filtering: after
centroid zeros/duplicate coordinates have been handled, but including profile
zero samples. `retained_centroids` counts accepted centroid entries or detected
profile masses. Their ratio on profile data is not a noise-removal percentage.

## 3. New MS1 algorithm: low-intensity background across m/z

### Default since 0.10: fixed 100 Da bins

Each MS1 scan is processed independently. Bin edges are anchored at zero:
`[0, 100)`, `[100, 200)`, `[200, 300)`, etc. A point exactly at 200 belongs
to `[200, 300)`. `--noise-window-da` controls bin width; `--noise-layout bins`
is the default. Bins never pool scans, retention times or polarities.

For each bin with `N` positive intensities, select its lowest
`k = max(5, ceil(0.10 × N))` values. Compute the mean and population standard
deviation (`ddof=0`) of **that same subset**, then use
`threshold = max(min_intensity, μ + 3σ)`. All points in the bin share one
threshold, including points in the selected subset; there is no per-point
leave-one-out exclusion. Zero padding does not contribute. Fewer than five
positive samples gives an unsupported estimate and uses only `min_intensity`
(zero by default). Empty neighboring bins are not borrowed.

Numba partitions each bin once to select the low tail; it does not sort a
separate neighborhood for every mass. Moments are calculated after shifting
by a selected value to preserve small variation on a large baseline. Overflow
raises a file error. Centroid intensities must be strictly above the threshold;
profile samples use the complete grid and the same thresholds for centroiding.
Intensity itself is not smoothed or baseline-subtracted.

This change simplifies the regional noise estimate; it does **not** guarantee
stronger filtering or correct over-splitting of persistent chromatographic
background. A standard deviation of a truncated low tail is not a calibrated
false-positive probability. The [raw EIC audit](../benchmarks/results/EIC_REVIEW.md)
shows both plausible weak peaks and excessive splitting. Peak tracking and
chromatographic segmentation have not been changed in 0.10.

### Historical sliding estimator, still selectable

Use `--noise-layout sliding --noise-window-da 50` to reproduce the v0.9 noise
neighborhood (±25 Da). A width of 100 gives the optional ±50 Da sliding window.
The following subsections document this historical sliding implementation.

This estimator operates separately on every MS1 scan. It estimates local
background along the **physical m/z axis**, not along retention time. It does
not smooth the reported chromatographic intensities or subtract an estimated
baseline from accepted centroid heights.

Let the sorted scan be `(m_i, I_i)`, `i = 0…N−1`, with nonnegative intensity.
The default `noise_estimator="low_tail"` uses the **lowest-intensity** fraction
of the neighboring positive points. "Lowest 10%" refers to rank by intensity,
not the last 10% of m/z coordinates or a percentage of cumulative intensity.
There is no smoothing of the whole spectrum or averaging of all local peaks.

### Sliding reference subset, moments and threshold

Let `W` be the full window width (v0.9: 50 Da, ±25 Da), `p` the low-intensity fraction
(0.10), `B` the minimum selected background count (5), `s` the detection
multiplier (fixed at 3), and `F` the absolute intensity floor (0). For each candidate:

```text
J_i = {j : |m_j − m_i| ≤ W/2, I_j > 0, j ≠ i}
C_i = number of points in J_i

if C_i < B:
    mu_i = sd_i = 0
    threshold_i = F
else:
    k_i = max(B, ceil(p * C_i))
    S_i = the k_i lowest-intensity entries of J_i
    mu_i = sum(I_j for j in S_i) / k_i
    sd_i = sqrt(sum((I_j − mu_i)² for j in S_i) / k_i)
    threshold_i = max(F, mu_i + s * sd_i)
```

**Both the mean and standard deviation come from the same selected subset.**
The detection multiplier is fixed at 3 in both estimators. `Config`, JSON and
the CLI reject any other `noise_sigma`; the field is retained so explicit
`noise_sigma=3` configurations remain valid. Adjust `noise_fraction` to compare
background subsets, such as 0.05 versus 0.10, while retaining the same 3σ rule.
The standard deviation is a population value (`ddof=0`). The candidate never
influences its own estimate. Window endpoints are inclusive, and windows
truncate at spectrum boundaries without expanding across missing m/z regions.
Each reference point has equal weight. Ties at the subset boundary have the
same intensity, so which tied point is selected does not change the moments.

With 40 eligible neighbors, `ceil(0.10 * 40)=4`, so the minimum expands the
subset to the weakest **5 of 40 (12.5%)**. With 51 neighbors, the subset has
6 points. With only 4 neighbors, the estimate is unsupported: it uses `F`
rather than pretending the minimum was met. Thus the effective fraction can
exceed 10%, especially in sparse centroid exports. With the zero default floor,
unsupported windows retain positive centroids. Profile zero padding is never
a reference point; it remains in the mass-detection grid as a boundary.

For centroid MS1, retain an entry only when `I_i > threshold_i`; original mass
and intensity are unchanged. For profile MS1, estimate from positive profile
samples, evaluate the background/threshold at **every** grid point including
zeros, and pass the complete grid to §4. The background for profile scans
therefore describes low positive **samples**, not already detected mass peaks.
`noise_mode="none"` skips either estimator and returns zero background/std
plus the configured floor everywhere.

### Efficient rank selection and stable variance

The Numba kernel sorts positive intensities once per scan to assign stable
intensity ranks. A segment tree over those ranks stores the active count,
mean and sum of squared deviations from the mean (`M2`) for each node. Two
pointers advance an inclusive physical m/z window; adding/removing a sample
changes its leaf and updates ancestor counts. The candidate's leaf is temporarily
removed while querying the lowest `k_i` active ranks, then restored.

An order-statistic query descends the tree using subtree counts. Complete
low-rank subtrees are combined until exactly `k_i` references are covered.
It never sorts or copies each overlapping window. Centered moments of groups
`a` and `b` are merged as:

```text
n = n_a + n_b
d = mean_b − mean_a
mean = mean_a + d * n_b / n
M2 = M2_a + M2_b + d² * n_a * n_b / n
```

Empty nodes contribute nothing; population `sd = sqrt(max(0, M2 / k_i))`.
This avoids subtraction of two nearly equal large raw second moments, and
rebuilding nodes avoids accumulating subtractive error as strong peaks leave
the window. Before building the tree, an intensity range whose squared span
times the positive-point count would overflow float64 raises a file error.

Version 0.7.0 first queries the lowest-k prefix with the candidate present.
If the candidate rank is outside that prefix, excluding it cannot change the
selected subset, so no temporary tree updates are needed. Otherwise it is
removed, the prefix is queried again, and it is restored. The selection count
always excludes the candidate. Centered moments and μ + 3σ are unchanged.

Version 0.8.0 marks ancestor moments dirty instead of recalculating moments
that no low-tail query uses. A query refreshes its requested subtrees in
postorder; clean subtrees are reused and empty subtrees need no descent.
The centered merge order is unchanged. A cached prefix result is reusable only
when its size, selected ranks **and first unselected rank** have not changed.
The last condition preserves the query's subtree partition and float64 merge
order even when brighter points leave the window. Candidate exclusion is
temporary and never changes the cached result. These are computation shortcuts,
not changes to the selected background or threshold rule.

The one-time sort plus window updates and prefix queries take **O(N log N)**
amortized time and O(N) temporary storage: each dirty node is refreshed at most
once between updates. Tree storage is three numeric arrays and one boolean
dirty array with twice the next power of two of the positive-point count,
plus rank/order/output arrays.
Kernels are compiled and cached with `nogil=True`; there is no Python loop
per m/z point, full sliding-window matrix, or nested Numba thread pool.

### Previous optional estimator: iteratively clipped moments

`noise_estimator="clipped"` with `noise_layout="sliding"` preserves the 0.2.0 algorithm for comparison.
`noise_clip_sigma` and `noise_iterations` apply **only** to this estimator;
`noise_fraction` applies **only** to `low_tail`. The legacy estimator starts
from all points, including profile zeros, before clipping strong candidates.

Let `W` be the configured full window width, `K` the clipping-iteration count
(default 2), `c` the clipping multiplier (2), `s` the fixed detection multiplier (3),
`B` the configured minimum background-neighbor count, and `F` the absolute floor (0).
The current shared defaults are `W=100` and `B=5`. For the historical 0.2.0
behavior, select `clipped` with `noise_layout="sliding"`, `noise_window_da=25` and
`noise_min_neighbors=8`; the diagnostic script pins those historical settings.

#### Legacy window, moments and clipping

Start with background membership `v_i = true` for every point. At each point,
the background reference set is:

```text
J_i(v) = {j : |m_j − m_i| ≤ W/2, v_j is true, j ≠ i}
n_i    = number of points in J_i
mu_i   = max(0, sum(I_j for j in J_i) / n_i)
sd_i   = sqrt(max(0, sum(I_j² for j in J_i) / n_i − mu_i²))
```

Moments are supported only when `n_i ≥ B`; otherwise both are zero. The
variance is a **population** variance (divisor `n_i`, not `n_i−1`). Window
endpoints are inclusive; windows simply truncate at the spectrum edges.
Every retained sample has equal weight, so uneven profile sampling can affect
the estimate even though the window is measured in Da.

For each of the `K` iterations, calculate these moments using the current
membership mask and update all memberships simultaneously:

```text
v_i = v_i AND (window_i unsupported OR I_i ≤ mu_i + c * sd_i)
```

This is a scan-wide membership mask determined by each point's own local
window. It is not a separate iterative sorting/clipping operation for every
overlapping window. Exclusions from the background are monotonic: a discarded
background candidate does not re-enter on a later iteration. A strong peak
excluded from the background is still eligible to become a detected signal.

Recalculate moments once more after clipping and set:

```text
threshold_i = max(F, mu_i + s * sd_i)   if n_i ≥ B
threshold_i = F                        otherwise
```

For centroid MS1, retain an entry only when `I_i > threshold_i` (strictly
greater). Its original m/z and intensity remain unchanged. For profile MS1,
pass the background `mu_i` and threshold arrays to the detector in §4 without
first deleting samples. `noise_mode="none"` sets background/std to zero and
uses the configured floor everywhere.

#### Legacy rolling implementation

Two pointers move the window bounds monotonically through sorted m/z.
Running count, sum and sum of squares are updated when points enter or leave.
Each point is removed from its own reference set, even if it remains in the
background mask. Compensated additions/subtractions retain small contributions
when very large intensities enter or leave a window. Candidate subtraction is
performed before applying the stored compensation. Negative variance caused
by rounding is clamped to zero.

Each moments pass is O(N); total estimator work is O((K+1)N), with O(N)
temporary storage. Kernels are Numba-compiled, cached and release the GIL.
There is no per-mass Python window loop and no independent Numba thread pool.

### Interpretation and limits

Neither estimator chooses an instrument-specific absolute cutoff. A point below threshold
is not automatically known to be electronic noise; chemical background and
real weak analytes can also be removed. Conversely, intense persistent
electronic signals can survive.

Sparse windows intentionally use `F` alone. With the default zero floor they
retain positive centroids because there is insufficient local background
evidence to estimate noise reliably. The lowest-intensity subset is a **truncated distribution**;
its standard deviation is narrower than that of the full background. Applying
3 times this value does not estimate an upper three-sigma boundary for all
noise. For example, in the fixed-seed synthetic uniform-background test
(`U[8,12]`), this rule retains more than 80% of noise-only centroids. It also
produces many background traces on the old synthetic chromatographic fixture;
the previous "three known peaks, no extras" result applies to `clipped` only.
This more permissive behavior is part of the requested low-tail rule, not
evidence that the extra retained peaks are real analytes.

Dense unresolved chemical signal, broad profile peaks, changing sampling
density and unusual noise distributions can bias either estimate. Positive
profile tails can dominate the low subset. The multiplier is empirical:
`s=3` is not a calibrated significance test, signal-to-noise ratio or false-discovery rate.
No blank subtraction or cross-file normalization is performed.

The [same-scan comparison](../benchmarks/figures/README.md) records successive default settings and the historical estimator,
selected-count diagnostics, before/after plots and all per-centroid values.
Version 0.3.0 introduced the low-tail default (25 Da, 5%, minimum 8); 0.3.1
changed it to 50 Da, 1%, minimum 5. Version 0.3.2 uses 50 Da, 10%, minimum 5
and enforces a detection multiplier of 3. These changes update extraction
behavior and provenance; old output
sidecars fail the configuration/implementation match and are not resumed as
current results. Those noise-setting revisions did not change the feature-table columns.

## 4. New profile MS1 mass detector and centroiding

This detector receives the complete normalized profile grid, background and
threshold arrays from §2–3. Already centroided scans bypass this step.

### Locate maxima and integration intervals

1. Walk through internal samples. A peak must rise strictly above its left
   neighbor and eventually fall strictly below its right neighbor. Equal-height
   plateaus are treated as one maximum; the apex is the integer midpoint,
   rounded down. A maximum extending to either spectrum edge is excluded.
2. Find the lowest raw-intensity sample between adjacent maxima/plateaus;
   the first minimum wins a tie. These valleys, together with the spectrum
   ends, bound the interval available to each peak.
3. Reject an apex at or below its local threshold. Define its reference grid
   spacing as the smaller of the gaps immediately outside the plateau's left
   and right edges. The maximum allowed sampling gap is this spacing times
   `profile_gap_factor` (default 4).
4. Expand from the apex toward its bounding valleys. Stop before a neighbor
   at or below **that neighbor's** threshold, or before an m/z gap above the
   allowed gap. No interpolation through a threshold crossing or missing
   sampling region is performed.
5. Require at least `profile_min_points` (default 3) contiguous accepted samples
   in the resulting interval. This rejects isolated one-bin spikes even if
   their intensity is large.

### Calculate the mass centroid

Within the accepted interval let `y_i = max(0, I_i − mu_i)` and choose the raw
apex mass as origin `m_0`. For each adjacent pair use `x_i = m_i − m_0` and
`dx = x_(i+1) − x_i`. Integrate the piecewise-linear signal exactly:

```text
A_i = dx * (y_i + y_(i+1)) / 2
Q_i = dx * [x_i*(2*y_i + y_(i+1)) + x_(i+1)*(y_i + 2*y_(i+1))] / 6

centroid_mz = m_0 + sum(Q_i) / sum(A_i)
```

Require positive total area. Relative mass coordinates reduce cancellation
when the absolute m/z is large. Actual `dx` values account for nonuniform
sampling; this is not an intensity-weighted average that assumes equal-width
bins. Output intensity is the **original raw profile apex height**, and the
output noise threshold is the threshold at that apex. The integrated mass
profile area is used only to locate m/z; it is not substituted for intensity.

Maxima, valleys and accepted intervals are traversed in compiled loops with
linear work in the profile size and O(N) temporary storage. Valleys prevent
merging separately resolved maxima; sampling-gap checks prevent centroiding
across absent parts of a profile. The method does not deconvolve unresolved
shoulders, fit instrument-specific line shapes, or reconstruct truncated
tails. Saturation, undersampling and threshold-truncated asymmetric peaks can
affect mass accuracy. Profile-apex heights are not guaranteed to equal a
vendor's centroid intensities.

The real 50-file validation corpus contains centroid exports. Profile accuracy
has been checked with synthetic Gaussian peaks on nonuniform grids, resolved
neighbors, plateaus, discontinuities, spikes, and complete profile-to-feature
fixtures; it has not been validated across real vendor profile datasets.

## 5. Accelerated ion-trace construction

Each polarity has a separate tracker. On each MS1 scan:

1. Visit existing traces in descending **most recently observed positive
   intensity**, with stable ties. A missed scan does not replace that priority
   value with zero.
2. Close a trace if the elapsed time since its last matched positive point is
   greater than `max_gap_seconds` (15 s).
3. Search sorted centroid masses around the trace's running intensity-weighted
   m/z `m`, within `max(mz_tolerance_da, m * mz_tolerance_ppm * 1e−6)`.
   Defaults are 0.002 Da and 10 ppm. Choose the nearest **unused** centroid;
   equal mass-error ties choose the lower mass/earlier sorted entry. If the
   nearest centroid is already taken, another available candidate in the
   window can still be used. Each centroid joins at most one trace.
4. On a match of mass `x`, intensity `I`, and previous total weight `S`, update
   `m ← m + (x − m)*I/(S+I)` and `S ← S+I`. Save the original point, time and
   threshold; reset the consecutive-gap counter.
5. Otherwise increment the gap counter. Up to `max_gap_scans` (2) missing MS1
   scans are represented by zero intensity at their actual acquisition times.
   Close the trace when the counter exceeds that limit. The scan that exceeds
   the limit is not appended. The time-gap rule still applies independently.
6. Start new traces from every unclaimed centroid. Flush remaining traces at
   end of file. Discard a trace with fewer than `min_scans` observed positive
   points (5 by default).

Active-trace state and linked histories live in contiguous numeric arrays.
Integer free lists recycle both trace slots and history nodes immediately on
trace closure. Capacity grows geometrically when needed; it does not grow
merely because another completed feature was encountered. Numba handles
matching, node insertion and draining without a Python object per centroid.

Version 0.8.0 uses a batched extraction interface: copy and recycle a completed
trace, segment and quantify it, then retain only its feature summaries and
optional two-column grouping trace. The raw chromatogram stays inside compiled
code; there is no Python/Numba round trip per completed trace. Calibration
retains its raw-trace interface. Pool reuse cannot modify already emitted
traces. Rare ambiguous peak priorities carry their original trace to the SciPy
fallback in the same output position.

For F active traces and N current centroids, the scan requires a priority sort
(approximately F log F), binary-search starts (F log N), candidate-window
inspection, and linear insertion/draining work. It avoids an all-traces by
all-centroids distance matrix. Pools keep their high-water capacity until the
file ends. A persistent ion can still retain a long history: memory is not
constant in acquisition length.

Compared with upstream, the default gap allowance is 2 rather than 10 scans,
the target is a weighted mass rather than the last mass, and both ppm and
absolute limits are available. Therefore this is not guaranteed to reproduce
all upstream feature boundaries, even with adaptive noise disabled.

## 6. Chromatographic segmentation and quantification

Remove leading/trailing zero observations from a completed trace. With
`segment=True` and positive `smoothing_seconds` or `peak_distance_seconds`,
interpolate onto a uniform RT grid near the median interval, capped at four
times the original observation count. Convert seconds to grid samples for
Gaussian smoothing and minimum peak distance; map valleys back to original
observations. Calibration sets physical-time parameters by default. Otherwise
use legacy `smoothing_sigma=0.6` and `peak_distance=3` samples.

Pad the intensity vector with one zero at each end and apply the order-zero
Gaussian filter with reflected boundaries, equivalent to SciPy's
`gaussian_filter1d`. Find maxima with prominence
`peak_prominence × maximum_smoothed_intensity` (default 0.1). Resampling and
smoothing affect segmentation only, never original-sample quantification.

Version 0.7.0 fuses resampling, reflected Gaussian smoothing (radius
`int(4*sigma + 0.5)`), plateau midpoints, peak selection and quantification in
Numba. Minimum distance is applied **before** prominence, matching SciPy's
selection order; prominence searches extend until a strictly higher sample or
the signal boundary. Exact or near-equal distance priorities (within eight
float64 epsilons relative to the larger magnitude) use the retained SciPy path
so its unstable-sort tie policy is preserved. The named summary record carries
original trace bounds, avoiding implicit positional fields in the Python layer.
No fast-math or float32 conversion is used. Strided measurement loops avoid
allocating several temporary arrays per peak. Compiled exp/sum/dot operation
order can change final float64 bits: numerical parity is tested with tight
tolerances, while feature boundaries and MS/MS identities must match exactly.
See the [performance and equivalence report](../benchmarks/results/PERFORMANCE.md).
The SciPy compatibility source and BSD notice are credited in `NOTICE`.

Version 0.8.0 materializes reflected Gaussian boundary samples once per trace
instead of repeating reflection arithmetic for every convolution pair. The
kernel, outside-in summation order, padding and measured observations stay
unchanged. Named `PeakSettings` and `PeakCandidate` records carry the batch
options, summaries and fallback status across the native/Python boundary.

Between consecutive detected maxima, cut at the first minimum of the smoothed
signal. Add the trace's first and last points as outer boundaries. Adjacent
segments share their valley point but never share an integration interval.
Each segment must contain at least `min_scans` nonzero samples. If fewer than
two maxima exist, retain the unsplit trace if it passes the length criterion.
`--no-segment` skips peak splitting entirely.

### Shape guard introduced in 0.11

Original MassCube computes Gaussian similarity for final features, but does
**not** use that score to gate feature splitting. Its splitter returns an
unsplit trace when the whole-trace noise score exceeds 3; a Gaussian-similarity
threshold of 0.7 is used in its later alignment workflow. MassCubeX versions
through 0.10 used Gaussian **smoothing**, without either that noise-score guard
or final Gaussian-similarity reporting. This allowed long, fluctuating background
traces to yield many intervals. See the [before/after review](../benchmarks/results/PEAK_SHAPES.md).

With `peak_shape_guard=True` (default), calculate upstream's noise score before
segmentation. For height `H=max(I)`, set `z_i=I_i/H` when `I_i >= 0.05H`, otherwise
zero. The score is

```text
N = max(0, (sum(abs(z_(i+1) - z_i)) + z_0 + z_last) / 2 - 1)
```

Traces with fewer than five points or nonpositive height have score zero.
`N <= peak_noise_score_limit` (default 3) keeps the existing splitting behavior.
For `N > 3`, apply the following additional evidence checks. A large score
is a trigger for scrutiny, not a reason to delete the trace or automatically
merge every peak: five separated ideal peaks can also have `N > 3`.

1. Linearly interpolate only the tracker's internal missing-scan zeros for
   **split detection**. These runs are already bounded by `max_gap_scans` and
   `max_gap_seconds`. This prevents tolerated missing observations from
   manufacturing deep valleys. Do not join independently completed traces.
   Keep original intensities, zeros, timestamps and MS1 thresholds for output
   measurements and ion-grouping traces.
2. Resample/smooth/find candidate maxima as described above. Use their original
   adjacent valleys as fixed candidate basins. Rejected neighbors do not enlarge
   a basin during scoring. The outer basins end at the actual trace boundaries;
   zero padding supplies no baseline evidence.
3. Estimate a robust temporal fluctuation scale from triples of consecutive
   positive, original observations. For each interior sample, use

   ```text
   w_i = (t_i - t_(i-1)) / (t_(i+1) - t_(i-1))
   e_i = [I_i - (1-w_i) I_(i-1) - w_i I_(i+1)] / sqrt(1 + (1-w_i)^2 + w_i^2)
   sigma_time = median(abs(e_i - median(e))) / 0.6744897501960817
   ```

   Linear trends cancel, unequal scan spacing is accounted for, and triples
   containing missing observations are excluded. Fewer than five residuals
   leave the estimate unsupported (`sigma_time=0`), so only the other shape
   checks apply. The normalization estimates per-observation sigma under
   independent additive Gaussian noise; it is a robust heuristic for real LC-MS
   fluctuations, not a calibrated false-positive probability.
4. Set the candidate baseline to the higher of its two smoothed valley/boundary
   intensities. Require positive smoothed apex excess at least
   `peak_min_snr × sigma_time` (default 5). This chromatographic check is separate
   from the unchanged **mass-axis MS1 threshold μ + 3σ**.
5. Inside that basin, take positive original observations and subtract the
   baseline, clipping negative excesses to zero. Require at least five observed
   samples and at least three above **half the maximum raw excess height**, with
   support on both sides of the smoothed candidate apex. Interpolated points do
   not count. This rejects isolated spikes even when smoothing spreads them.
6. Require a Gaussian-template Pearson correlation of at least
   `peak_min_gaussian_similarity` (default 0.7). Accept either the symmetric
   apex/FWHM template below or a two-sided template whose left and right sigma
   equal their respective half-height distances divided by `sqrt(2 log(2))`.
   The latter accommodates peak fronting and tailing. These guard scores use
   baseline-corrected positive observations; they are not the exported score.
7. Split only between the surviving candidates, using the first smoothed valley
   between each pair. Fewer than two survivors retains the entire original
   trace as one feature, subject to the existing minimum observation count.

`--no-peak-shape-guard` reproduces the previous segmentation. The limit, S/N
multiplier and similarity threshold are explicit CLI/JSON settings, not
automatically tuned. The compiled and independent NumPy/SciPy paths implement
the same rules, including the equal-priority fallback.

### Shape diagnostics and interpretation

`gaussian_similarity` reports the symmetric MassCube-style shape score on each
final detected feature's **original** trace, including missing-scan zeros and
without baseline subtraction. Set the Gaussian center to the first raw apex;
find the first sample below half height on each side, using the trace boundary
if no crossing exists. Set sigma to that full width divided by 2.355 (minimum
1e-6 seconds), generate the Gaussian template, and calculate Pearson correlation
clamped to [0, 1]. Normalize both profiles to unit apex before arithmetic to
avoid scale-dependent overflow. Fewer than five points, zero height or zero
variance gives zero. Output rounds the score to two decimals.

`segmentation_status` is `split` or `unsplit`, prefixed with `guarded_` when the
whole-trace noise score triggered the extra checks. `disabled` means segmentation
was switched off. MS2-seeded forced measurements are not segmented by this
algorithm and have `not_applicable` with a null Gaussian score. A low Gaussian
score does not remove a row. A retained background-like trace retains its
original height, integrated area and all assigned MS2 references; its area
should not be interpreted as a clean chromatographic compound peak.

These diagnostics do not identify contamination or determine biological truth.
Strongly structured background can resemble elution peaks. Conversely, weak,
poorly sampled or irregular true peaks in a noisy trace may remain combined;
the guard does not promise to resolve every coeluting isomer. Acquisition gaps
and failed mass continuity still produce separate traces. Neither raw-data
baseline subtraction nor a cross-gap trace-merging algorithm is introduced.

All quantitative measurements use original, unsmoothed intensity observations,
including zeros for tolerated missing scans. For a segment `(m_i, I_i, t_i)`:

```text
feature_mz   = sum(m_i * I_i) / sum(I_i)
peak_height  = max(I_i)
rt_apex      = t_i at the first maximum original intensity
peak_area    = sum((t_(i+1) − t_i) * (I_i + I_(i+1)) / 2)
scan_count   = count(I_i > 0)
```

Times and integration boundaries are seconds; area is intensity × seconds.
Neither endpoint is extrapolated to an invented zero outside the feature.
For profile input, each `I_i` is the raw mass-peak apex height from §4.
Reported area and height receive no chromatographic baseline subtraction.
Filtering still affects which observations and boundaries survive, so this
does not imply unbiased recovery of signal below the detection threshold.

`quant_method` records the preferred measurement (`peak_height` by default,
or `peak_area`) once in the JSON sidecar. Both measurement columns remain;
schema 6 omits the redundant `intensity` and `quant_method` columns.
`apex_noise_threshold` is the retained MS1 threshold at the chromatographic
apex. The derived `apex_threshold_ratio` is also omitted; dividing height by
threshold is not a statistical S/N estimate. Measurements are rounded only
at serialization: mass coordinates/radii to 5 decimal places, all RTs in
seconds and intensities/areas/thresholds to 2. Rounding can change interval
membership if downstream code reconstructs a query from table coordinates;
use source scan identifiers and full-precision settings for exact provenance.
Non-finite feature values fail the file instead of publishing invalid quantification.
Calculations remain float64. Feature IDs follow trace-completion and
segment order, starting at 1 separately in each file, rather than mass order.

## 7. MS/MS provenance, association and representative selection

### Read metadata and scores

For mzML, read the selected ion's m/z, intensity, charge and parent spectrum
reference, plus the isolation window's lower/upper offsets. If selected-ion
m/z is absent, an isolation-target m/z is usable only when both offsets are
known. More than one precursor/selected ion is treated as multiplexed.
For mzXML the analogous fields come from `<precursorMz>` and its attributes,
including `precursorScanNum`, `precursorIntensity`, charge and window width.

Save the exact native ID. A numeric native scan number is extracted only from
an entirely numeric identifier or a standalone `scan=N` field. Other native
formats retain their ID with a missing numeric scan number. Independently,
the zero-based spectrum index counts every source spectrum, including ignored
MS levels. An index is never silently substituted for a native number.

By default, decode each MS2 array and save its float64 intensity sum. The sum
includes all exported nonnegative peaks, including residual precursor peaks.
Sorting, removing zeros and merging identical masses preserve the sum; explicit
profile negatives are clipped to zero. Profile MS2 is scored by the sum of
samples, not a profile mass integral, so different sampling densities can
affect those scores. The vendor's TIC may describe pre-export data and is not
substituted for the calculated sum.

Without ion-fragment analysis, arrays are then discarded. With that analysis,
only explicitly centroided MS2 retains product masses: first retain peaks
with intensity/base-peak intensity ≥ `isf_min_relative_intensity` (0.05), then
keep at most `isf_max_peaks` (100), ranked by intensity with stable mass-order
ties. **The representative score was calculated before this pruning.** MS1
adaptive noise estimation is never applied to MS2.

Precursor intensity is also recorded when positive. If absent or zero, inspect
the latest MS1 of compatible polarity; an unknown polarity can fall back only
when one previous-polarity stream exists. Require the supplied parent reference
to match, when present, a scan age from 0 through `precursor_max_age_seconds`, and a
mass match within the MS2 tolerance. Use the nearest centroid's original
intensity, or the maximum original profile sample in that interval. The source
is labelled `metadata`, `preceding_ms1` or `unavailable`. Older MS1 scans are
not retained just to resolve an earlier parent reference.
After the complete file is read, mixed-polarity runs invalidate preceding-MS1
intensity fallbacks whose precursor or source MS1 polarity was ambiguous.
Explicit precursor intensities from the source metadata are retained.

### Assign MS2 after chromatographic splitting

Build a sorted feature-mass index. Skip and count multiplexed precursors,
missing precursor/time metadata, and known isolation widths greater than
an explicitly enabled `ms2_max_isolation_width`. The default is **0 (disabled)**;
versions through 0.3.2 defaulted to 4 Da. Width means the full isolation span:
the sum of lower and upper offsets in mzML, or `windowWideness` in mzXML.
It is unrelated to the MS1 noise window or precursor-to-feature mass tolerance.
Neither a small width nor a missing width proves DDA acquisition. A saved
configuration with a positive limit retains it unless explicitly overridden.

For each remaining MS2 at precursor mass `p` and time `t`, find features within
`max(ms2_mz_tolerance_da, p * ms2_mz_tolerance_ppm * 1e−6)` (configured fallback 0.01 Da / 20 ppm)
whose integration interval, extended by `ms2_rt_tolerance` (fallback 2 s) at both ends,
contains `t`. Known opposite polarities are incompatible. An unknown-polarity MS2 in a mixed-polarity file remains initially unassigned,
even if only one feature happens to match. An explicit reference to the latest
retained known-polarity MS1 may resolve its polarity, with the inference recorded.
Learned mass and time settings are resolved separately per polarity.

Assign to at most one remaining feature, in this priority order:

1. Inside the actual, unextended feature interval.
2. Smaller absolute precursor mass error.
3. Smaller distance from the feature apex time.
4. Higher feature peak height.
5. Lower feature ID.

No hit initially means unassigned, not a file failure. The next step creates
MS2-seeded features by default. The initial outcome for each scan distinguishes
missing mass/time, multiplexed precursors, an optional width limit, no mass
candidate, no compatible polarity, no compatible time, and ambiguous polarity.
These rules do not deconvolve coisolated precursors or DIA windows.

### New MS2-seeded features and forced MS1 gap filling

**Upstream audit.** In the cloned commit, [`allocate_ms2_to_features`](https://github.com/Philipbear/masscube/blob/356748ede691862b9a87f67812db5e728a703b9d/src/masscube/raw_data_utils/core.py#L227)
matches existing features using a 0.1 Da precursor tolerance and feature time
boundaries, selects the highest peak if several match, and skips unmatched
scans. It does not create features for them. Upstream's
[`gap_filling`](https://github.com/Philipbear/masscube/blob/356748ede691862b9a87f67812db5e728a703b9d/src/masscube/alignment.py#L274)
instead fills missing measurements after cross-file alignment, from the
unfiltered cache using a local EIC maximum and trapezoidal area. MassCubeX
adopts that forced-measurement concept for **new single-file MS2 seeds**;
dataset alignment is still absent.

With `ms2_seeded_features=True`, process initially unmatched MS2 in source
spectrum order. Keep all existing detected-feature assignments unchanged.
A seed requires one precursor, a usable precursor m/z and retention time,
and no exclusion by an explicitly configured isolation-width limit.

A polarity with at least 20 qualified representatives uses refined time bounds:
extract one narrow EIC anchored at the exact precursor mass, find local
chromatographic bounds, then merge queries only when they recover the same
mass/polarity and original MS1 apex. Version 0.10 removes the broad mass search
and its 3× candidate-ranking rule (`ambiguous_mass` in older tables).
The complete algorithm, quality rules, safety limits and constants are in
[AUTO_TUNING.md](AUTO_TUNING.md#refined-raw-ms1-gap-filling).

With automatic calibration off or insufficient representatives, use a fixed
window `[max(0, t − w), t + w]`, where `w=5` seconds by default. Reuse a seed
only within its separate `ms2_seed_rt_tolerance` (5 seconds) of the original
anchor and within the resolved MS2 mass tolerance. Ties use mass difference,
anchor-time difference and feature order. The window never expands through
repeated MS2 scans. Both modes cap the mass radius at the normal MS1 radius:
`min(max(gapfill_Da, p × gapfill_ppm × 1e−6), max(MS1_Da, p × MS1_ppm × 1e−6))`.
The two configured radii default to 0.002 Da / 10 ppm, and calibration normally
sets both to the same learned values. A wider explicit gap-fill override is
capped; a narrower one is honored. The precursor center never moves.

Both modes stream raw MS1 once for all seeds using per-polarity active-query
schedules. No adaptive detection threshold or minimum feature scan count is
applied. Reader normalization still sorts, promotes m/z to float64, sums exact
duplicates and clamps negative profile samples to zero. Refined profile queries
centroid with zero baseline/threshold; fixed-window profile queries take the
maximum raw sample. Neither sums intensity across a profile peak.

```text
J_s = maximum unfiltered intensity in the narrow mass interval, or zero
height = max(J_s) within the reported bounds
area = sum((t_(s+1) − t_s) × (J_s + J_(s+1)) / 2) within those bounds
```

Use original scan times in seconds. Exclude intervals longer than
`gapfill_max_gap_seconds`; do not extrapolate to window edges or subtract a
baseline. A single observation has zero area. Refinement smoothing determines
bounds only. Overflow fails the file. These forced measurements can recover
below-noise signal, but can also measure background or interference; an MS2
trigger is not proof of MS1 identity. Separate seeds need not represent
separate compounds or disjoint integrations.

Positive and negative MS1 are never pooled. Unknown MS2 can inherit a unique
file MS1 polarity, or an explicit known-polarity parent reference. Unresolved
polarity in mixed files gives unavailable quantification. Entirely unknown MS1
may support one known MS2 polarity, but never both positive and negative MS2
streams. Known polarity without corresponding MS1 scans stays unavailable.

Rows have `feature_origin=ms2_seeded` or `ms1_detected`. Seeds retain original
`seed_precursor_mz` and `seed_rt_seconds`; seeded `mz` stays at the precursor
anchor while refined time bounds may change. `gapfill_mz_radius_da` records the measurement
radius. `gapfill_status` distinguishes `fixed_window`, `peak_bounds`,
`single_observation` and unsuccessful refinement. `scan_count` counts positive
MS1 observations. Measured seeds have no detection threshold. When no MS1 ion
is present (including no local signal near the seed time), estimate regional
noise independently in up to the three closest compatible MS1 scans in the
query window. Use the median of their positive thresholds; each threshold is
the same bin's μ + 3σ, including any explicit minimum floor. Unsupported zero
thresholds supply no imputation evidence. Set height to 0.1× that reference,
area to missing, `scan_count=0`, and `ms1_quant_status=imputed_noise`.
`apex_noise_threshold` holds the reference used; there is no measured apex for
this row. `gapfill_status` retains the absence reason. No trace is manufactured
for grouping. `--gapfill-noise-fraction 0` disables imputation. Unknown polarity,
missing scans, unsupported noise, local-peak mismatch and sampling limits do
not receive fabricated heights. An actual below-threshold MS1 ion is measured
normally, even if weaker than this possible placeholder.

| `ms1_quant_status` | Meaning and quantitative output |
| --- | --- |
| `detected` | Normal MS1 detection |
| `gap_filled` | Raw MS1 measurement; height/area populated |
| `imputed_noise` | No MS1 signal; height is 0.1× regional threshold, area missing |
| `no_ms1_signal` | Compatible scans but no signal and no usable/enabled imputation; blank |
| `no_ms1_scans` | No compatible MS1 scans; quantitative fields blank |
| `ambiguous_polarity` | No defensible MS1 polarity; blank |
| `ambiguous_mass` | Historical v0.9 candidate ranking only; no longer emitted |
| `no_local_ms1_signal` | Broad query has ions but none near the seed time; blank |
| `no_local_peak` | Recovered peak bounds do not support the seed time; blank |
| `sampling_limit` | Query exceeds the documented bounded collection limits; blank |

Every assigned scan belongs to exactly one detected or seeded row, and all
native IDs, numeric scan numbers when available, and source spectrum indices
are exported. The same representative-selection rules apply to seeded rows.
Use `(source file, zero-based spectrum index)` for unique reverse lookup;
native numeric scan numbers may repeat across controllers or be absent.

The sidecar stores `ms2_matched_detected`, `ms2_initially_unassigned`, initial
reason counts, `ms2_seeded_features`, `ms2_seeded_scans`, gap-fill status counts,
and final assigned/unassigned counts. `unassigned_ms2` records every unresolved
scan's native ID, source index, numeric number, available precursor/time/polarity
metadata and reason. Those include missing precursor m/z/time, multiprecursor
spectra and a user-enabled width limit; no fictional single precursor is
created. An MS1-free file with usable MS2 seeds still produces rows, with blank
MS1 quantification and `no_ms1_scans`. A file with neither MS1 scans nor usable
MS2 seeds is an error. Corrupted inputs follow the ordinary file-error policy.

Under default settings, every parsed MS2 with usable single-precursor m/z/time
maps to one feature. This is a bookkeeping guarantee, not a guarantee that a
coisolated or DIA spectrum has a unique biological precursor. `--no-ms2-seeded-features`
turns off row creation and records unmatched scans instead; `--no-ms2` opts out
of MS2 parsing/association entirely. File-level DIA detection is not implemented.

### Choose the representative and write references

Within each feature, save every associated scan number, native ID and spectrum
index in acquisition order. Their arrays correspond position by position;
missing numeric scan numbers are JSON `null`, not removed entries.

| `best_ms2` | Primary ranking | Tie rules |
| --- | --- | --- |
| `fragment_intensity` (default) | Largest complete exported MS2 intensity sum | Nearest feature apex, then earliest spectrum index |
| `precursor_intensity` | Largest metadata/fallback precursor intensity | Nearest feature apex, then earliest spectrum index |
| `apex` | Smallest absolute time distance to feature apex | Higher available precursor intensity, then earliest spectrum index |

A missing intensity ranks below a known zero. When an intensity-based method
has no available score in any associated scan, selection becomes
`apex_fallback`. Default pipeline MS2 summation normally supplies a score even
for an empty MS2 array (zero). The sidecar records `config.best_ms2` and a sparse
list `best_ms2_apex_fallback_feature_ids` identifying these exceptional rows.
Best-scan columns retain the identifiers, full fragment sum when decoded, and
precursor intensity/source. Schema 6 omits the repeated method and derived
score columns. Selection still uses full-precision values before rounding.

Choosing precursor/apex ranking while fragment grouping is off avoids MS2
binary decoding. `--no-ms2` disables both association and those metadata
records. Unknown/profile MS2 may still be associated and scored, but does not
supply centroid product-mass evidence. Complete MS2 spectra remain in the
source file rather than being duplicated into every feature table.

## 8. Adducts, water losses and in-source-fragment candidates

Grouping is on by default. With `group_ions=True`, keep feature intensity/time
traces until relationships have been evaluated. The enabled catalog is ions
**1–15 plus 22–23**; 16–21 (multiple charges and trimers) remain available but
disabled. The full forms, atomic-mass shifts and custom schema are documented
in [ION_GROUPING.md](ION_GROUPING.md) and emitted by `masscubex adducts --json`.

For molecule count `n`, signed charge `q`, neutral mass `M`, and total ionic
shift `s`, use `m/z = (n*M+s)/abs(q)`, with inverse
`M = (abs(q)*m/z−s)/n`. Built-in shifts sum the added/subtracted neutral atomic
masses and subtract `q * electron_mass`. The sign of `q` specifies polarity;
it never makes the m/z divisor negative. Custom shifts are supplied directly
as total ionic mass shifts, not per-charge offsets.

For each pair of enabled ordinary ion forms in a known polarity, derive
neutral masses from each potential source feature and find matching predicted
partner masses. Require a positive neutral mass and candidate m/z. Use the
larger of 10 ppm and 0.002 Da by default. Known associated precursor-charge
magnitudes restrict allowed forms; if several magnitudes are recorded, a form
must match at least one. Unknown charge permits hypotheses, not confirmation.

Adduct mass bounds are found by a linear sweep of sorted query bounds and sorted
feature masses. Nonmonotone bounds from custom ion forms/extreme tolerances fall
back to independent binary searches. Small MS2 fragment queries retain binary
searches. The compiled candidate iterator rejects self-pairs,
empty traces and incompatible apex times in a compiled loop before creating
Python correlation/evidence objects. It reuses per-feature apex/count/charge
data and vectorizes MS2 product-mass searches. Candidate ordering, inclusive
tolerance boundaries, correlation calculations and evidence ranking remain
unchanged.

Every relationship also requires:

- Apex time difference ≤ `group_rt_tolerance` (3 s).
- At least `group_min_scans` (5) paired **positive** observations at exactly
  equal source MS1 times. Missing samples are not interpolated.
- Shared count / larger feature's positive scan count ≥ `group_min_overlap`
  (0.5). This guards against tiny subsets of otherwise different traces.
- Pearson correlation ≥ `group_min_correlation` (0.9), calculated with online
  mean/variance/covariance updates. Fewer than three paired points or constant
  traces return zero correlation, never a perfect match. Correlation is
  clamped to [−1, 1] for rounding. Pairwise evidence is cached during the file.

Ordinary adduct links record both hypothesized forms, average inferred neutral
mass, shared count, correlation and mass error. That error is
`1e6 * (observed_partner_mz−predicted_partner_mz)/predicted_partner_mz` in the
enumerated direction; the symmetric link retains the same evidence values.
No tallest feature is automatically declared to be protonated/deprotonated.

Water-loss forms 22/23 predict `[M+H-H2O]+` and `[M+H-2H2O]+` from a putative
`[M+H]+`, corresponding to approximately 18.010565 and 36.021129 Da. Mass and
chromatographic evidence can propose these as `neutral_loss` links without
MS2 support. Each explicitly reports `ms2_supported`; a positive value records
the actual supporting MS2 scan and product. Successive water-loss chains can
produce more than one valid hypothesis.

Generic `in_source_fragment` links require a matching centroid product from
any associated parent MS2, not only the best MS2. The current interpretation
uses singly charged or unknown-charge parent spectra and products below the
parent m/z by more than the mass tolerance. It searches child features using
the same mass/RT/correlation requirements. For each parent–child pair, prefer
the highest product/base-peak relative intensity, then smaller product mass
error, then earlier evidence spectrum index. A pair already labelled a water
loss does not get a redundant generic fragment link. Ordinary adduct
alternatives remain visible. `--no-isf` disables both types of fragment/loss
analysis while retaining ordinary adduct candidates.

Each link appears on both feature rows and uses local feature IDs. Fragment
links also identify parent/fragment roles. Sidecar relationship counts count
each hypothesis once, not both row appearances. Traces are released after
grouping. Every link is labelled `status="candidate"`: correlation and shared
MS2 products do not prove chemical identity or where fragmentation happened.
Rows are never collapsed, removed or assigned a confirmed compound label.

## 9. File outputs, identity and restarting

The CLI names outputs using the first 128 bits (32 hex characters) of the
SHA-256 of the resolved absolute source path, and shards them by the first two
hex characters. Up to 100 UTF-8 bytes of the source basename are included for
readability. Repeated basenames in different paths therefore do not overwrite
one another; moving a source changes its output ID. The Python single-file
API also accepts an exact `.parquet` or `.tsv` path; the extension selects the
format. Without a path, Parquet is the default. The CLI selects formats with
`--format parquet|tsv`.

Parquet uses an explicit, stable Arrow schema: float64 quantitative values,
nullable int64 identifiers, native MS2 lists and a list of typed structs for
ion evidence. Absent relationship-kind fields are null, not fabricated evidence.
Unknown evidence fields are rejected instead of silently discarded. Zstd level 3
compression is lossless **after output rounding**. At most 8,192 features are
buffered per row group; flat categorical strings use dictionary encoding.
TSV uses fixed decimal places and compact JSON arrays. The 29-column schema
omits five redundant fields listed in [USAGE.md](USAGE.md#feature-table-schema).
Both writers round binary64 values with Python's half-even rule and never
mutate the feature objects. Nested relationship masses have 5 decimals,
ppm errors and relative intensities have 2, and correlation has 6. Parquet
retains float64 storage because float32 cannot preserve five decimal places
throughout a typical mass range. The `output_precision` metadata describes
rounding; identifiers and lists remain exact. `iter_feature_rows` streams typed rows
from either format and normalizes absent relationship fields for comparison.

An exclusive `.parquet.lock` or `.tsv.lock` records PID, hostname and creation time. Write data and
metadata into temporary `.partial` files in the output directory; flush and
fsync them before atomic replacement. Publish the table first and its JSON
sidecar last as the completion marker. Ordinary exceptions remove temporary
files and release the lock. Parsing errors leave pre-existing outputs intact.
A failure between the two final renames can leave mismatched files, which
cannot pass resume validation. A hard-killed process can leave locks/partials;
inspect its recorded owner and remove a lock only after that worker has stopped.

The per-file sidecar uses schema version 7 and stores all processing settings, source
path/size/nanosecond mtime, counts and memory diagnostics, output format/size/SHA-256 and
an implementation/dependency fingerprint. The fingerprint hashes the numerical
and per-file processing source modules and versions of NumPy, Numba, SciPy,
lxml, PyArrow and optional pynumpress. CLI/report-only edits do not change that numerical
fingerprint. Resume also checks the configured software/schema version.
`output_bytes` and `output_sha256` apply to both formats; TSV keeps its legacy
`tsv_bytes` and `tsv_sha256` aliases. Old sidecars do not pass the new resume
validation. A different format uses a separate output path.

`timings_seconds` measures calibration, feature extraction, MS2 association,
MS2 gap filling, ion grouping and table writing. Feature extraction includes
decoding, adaptive noise, centroiding, trace tracking and peak summarization;
calibration and gap filling include their own source passes. Generator timing
excludes the consumer's work, including for the lazy MS1-only writer. Checksums,
locks, metadata publication and other setup are outside these stage timers.

`--resume` skips only when source identity, configuration, implementation,
version, output format, size and checksum all match and the completion metadata is
valid. Otherwise it recomputes. `--overwrite` always recomputes; without either,
existing outputs become a reported file error. Source identity is rechecked
after reading, but the entire input is not content-hashed. Same-size content
changes with restored timestamps are not detectable by this identity scheme;
keep inputs immutable. CSV/TSV quoting and JSON array encoding preserve IDs
and relationships without flattening or silently dropping missing values.

A valid MS1 run with zero retained features writes an empty Parquet table with
the complete schema, or a header-only TSV. Files
without MS1 may instead produce usable MS2-seeded rows with blank quantitative
values; without either MS1 or usable MS2 seeds they are errors. Scalar results
and sidecar statistics provide m/z,
RT, polarity and quantification for a future independent alignment process;
they do not implement alignment, RT correction or cross-file feature identity.

## 10. Millions of files: scheduling, failures and reports

### Bounded scheduling

Directory traversal uses a stack of open `scandir` iterators, without collecting
all filenames from a directory. Nested symlinked directories are not followed.
Manifests and error reports are read lazily. Duplicate paths are not globally
deduplicated, because a million-entry set would grow with the job; avoid
duplicate inputs or use separate batches with checksum-verified resume.

One worker processes files sequentially. With multiple workers, processes are
the default; threads are optional. At most `2 * workers` futures are submitted.
New inputs are consumed as outcomes become available. File state is never
shared between workers. The parent validates/warmups kernels, and spawned
process workers load/initialize their own caches. Numerical kernels use
`cache=True, nogil=True` without a nested Numba thread pool. Set a writable
`NUMBA_CACHE_DIR` on clusters, and control BLAS/OpenMP threads when running
many file workers. The benchmark found processes faster for this XML/Python
workload, despite GIL-free numerical kernels.

Default memory includes up to 64 MiB of conservatively charged decoded-scan
cache per worker (configurable/disableable), the current raw spectrum, active trace histories,
completed feature summaries, compact MS2 metadata and their associations.
Calibration adds bounded pilot traces/candidate pools and MS2 metadata
during its logical pass. Its one-trace quality screening runs in a separate
compiled module; all selection thresholds and robust per-polarity policies are
unchanged. Accepted candidates are re-measured with the previous NumPy/SciPy
arithmetic so calibrated mass/time bounds remain bit-reproducible; the numerous
rejected traces incur only compiled screening. Gap filling packs each bounded raw query into contiguous arrays once
so repeated candidate EIC queries avoid Python calls per sample. Refined gap filling retains bounded complete slices
per active query; total memory still depends on simultaneous queries.
Grouping additionally retains feature intensity/time histories and bounded
product-mass lists. Completed rows can be streamed directly when both grouping
and MS2 association are off. Report memory consists of counters, configuration
and the current outcome; full success/error lists are never accumulated.

### Corrupted inputs do not terminate the batch

The file-job boundary catches ordinary exceptions, including malformed or
truncated XML, invalid base64/zlib/gzip, invalid arrays/metadata, unsupported
formats, missing/unreadable inputs, existing-output conflicts and file-output
errors. Record the affected absolute input path, exception class and message,
then proceed to other inputs. Do not silently salvage a partially parsed file
as a complete feature table. Directory/manifest discovery errors are also
recorded while other accessible inputs continue. A discovery error names the
inaccessible directory or manifest; the contents that could not be enumerated
are unknown and are not presented as an exhaustive list of failed files.

An abruptly lost process worker can break its pool. Pending tasks then receive
`stage="worker"` errors if their outcomes cannot be recovered, and the scheduler
creates a new pool for undispatched inputs. An affected input is **not thereby
proven corrupt**; other tasks can be collateral failures. They are not blindly
resubmitted to the failed pool. One-worker mode runs directly in the caller;
neither it nor thread execution can isolate a native crash of the whole Python
interpreter. Multiple process workers provide child-process crash containment.
There is no per-file timeout in this version.

Global failures are distinct: invalid CLI settings, an unusable compilation
environment, an invalid retry-report format, inability to write the run report,
or user interruption can prevent completion. Where reporting remains writable,
the summary records `aborted` or `interrupted`, not successful completion.

### Automatically saved reports

Every `extract` run that reaches processing creates a unique UTC timestamp/UUID
run directory (`--print-config` and argument validation do not create reports):

```text
OUTPUT_DIR/reports/<run-id>/
    summary.json
    errors.tsv
    files.jsonl
```

`--report-dir PATH` changes the parent directory; it still creates a unique
child for each invocation. Earlier reports are never overwritten, including
when `--resume` or `--retry-errors` is used. The path is printed to stderr.

`files.jsonl` saves each success, verified skip or error as it is observed in
completion order; stdout still carries the same per-input JSON records.
`errors.tsv` contains only failures, with columns `input`, `input_kind`, `stage`,
`error_type`, and `error_message`. TSV quoting preserves embedded tabs/newlines;
use a TSV parser rather than splitting text lines. Input kinds distinguish
files from discovery failures on directories/manifests.

`summary.json` records run status, UTC start/update/end times, elapsed time,
counts of total/succeeded/skipped/failed outcomes, newly written feature count,
exception/stage counts, fatal error if any, command options and processing
configuration. Totals include discovery-failure records, not an invented count
of unseen files. Report paths are relative to the summary directory.

Logs are flushed after every outcome. They and an atomically replaced summary
are fsynced initially, at checkpoints when at least 30 seconds have elapsed
and another outcome arrives, and at normal/interrupted shutdown. There is no
unbounded in-memory error list. During one long-running file the summary may
remain at its previous checkpoint. After a hard kill it may still say `running`,
and flushed log records can be newer than the summary; a partially written
final log record is not a completed outcome. Per-file completion sidecars
remain the authority for `--resume`.

Final CLI exit codes:

| Code | Meaning |
| ---: | --- |
| 0 | Every discovered input finished or was verified/skipped |
| 1 | All possible work finished with one or more file/discovery errors, or a global run failure occurred; distinguish these using summary status/fatal error |
| 2 | Invalid command-line/configuration arguments before extraction |
| 130 | Keyboard interruption |

A nonzero final code after corrupt files does **not** mean the batch stopped
at the first bad file. Inspect `status="completed_with_errors"` and its counts.

### Retry and Python usage

After repairing/redownloading inputs, retry just the recorded sources:

```bash
masscubex extract --retry-errors features/reports/RUN_ID/errors.tsv \
  -o features --config params.json --workers 8 --resume
```

Replace `RUN_ID` with the run folder printed by the CLI. Reuse the original
processing configuration; it is recorded in the summary
but is not silently imported by `--retry-errors`. File rows are retried, directory
rows are rediscovered, and manifest rows are reopened. A second report records
the retry outcomes; the original error list remains intact. For killed workers,
resolve confirmed stale locks before retrying. To recheck the whole dataset,
rerun the original input command with `--resume`.

`process_file(...)` intentionally raises errors to its Python caller. The
batch API catches per-file failures and yields outcome dictionaries. To use
the same persistence from a Python scheduler:

```python
from masscubex import Config
from masscubex.cli import iter_inputs, run_batch
from masscubex.reporting import BatchReport

def main():
    config = Config()
    with BatchReport("features/reports", config, {"workers": 4}) as report:
        for result in run_batch(iter_inputs(["/data/study"]), "features", config, workers=4):
            report.record(result)

if __name__ == "__main__":  # required for spawned process workers
    main()
```

## 11. Complete processing-parameter reference

`masscubex defaults` prints the authoritative JSON defaults. Explicit flags
override `--config` values; unknown JSON keys and invalid types/ranges are
rejected. Python configuration objects are immutable. All processing settings,
including exact custom-ion definitions, are preserved in per-file provenance.

Automatic calibration is on by default; values below are configured fallback
settings unless explicitly protected. Per-polarity effective settings, support,
selected representatives and fallback reasons appear in `parameter_tuning`.
CLI overrides protect a value even when it equals its default. Non-default
JSON/Python values are also protected; `tuning_locked` protects default-valued
settings. Full derivations and override semantics: [AUTO_TUNING.md](AUTO_TUNING.md).

| JSON/Python key | CLI flag | Default / interpretation |
| --- | --- | --- |
| `auto_tune` | `--auto-tune` / `--no-auto-tune` | `true`; per-file, per-polarity calibration |
| `tuning_min_peaks` | `--tuning-min-peaks` | 20; minimum accepted representatives, cannot be lowered |
| `tuning_max_peaks` | `--tuning-max-peaks` | 100 per global/regional selection; at least the minimum |
| `tuning_mz_breaks` | `--tuning-mz-breaks` | `[200, 500]`; fixed regions [0,200), [200,500), [500,∞); empty JSON/Python sequence gives one unbounded region |
| `tuning_mz_regions` | `--tuning-mz-regions` | `null`; an explicit count (1–100) selects equal-width acquired-m/z regions instead of fixed boundaries |
| `tuning_min_gaussian_similarity` | `--tuning-min-gaussian-similarity` | 0.90; raw baseline-corrected symmetric Gaussian correlation for calibration anchors |
| `tuning_mass_sigma` | `--tuning-mass-sigma` | 3; mean + k sample SD of screened per-peak Q95 absolute errors |
| `mass_regions` | JSON/Python only | Empty by default; learned local MS1/MS2/gap-fill settings in effective configs; a final `upper: null` means unbounded |
| `tuning_search_ppm`, `tuning_search_da` | `--tuning-search-ppm`, `--tuning-search-da` | 100 ppm / 0.05 Da broad pilot radius components |
| `tuning_peak_limit` | `--tuning-peak-limit` | 128 strongest MS1 centroids per scan; minimum 32 |
| `tuning_block_scans`, `tuning_block_seconds` | `--tuning-block-scans`, `--tuning-block-seconds` | 512 scans / 120 s per pilot block |
| `tuning_locked` | JSON/Python only | Empty list; tunable names protected from automatic updates |
| `smoothing_seconds` | `--smoothing-seconds` | 0 means use scan units; calibrated when enabled |
| `peak_distance_seconds` | `--peak-distance-seconds` | 0 means use scan units; calibrated when enabled |
| `ms2_seed_rt_tolerance` | `--ms2-seed-rt-tolerance` | 5 s fallback; seed reuse independently of the search half-window |
| `gapfill_mz_tolerance_ppm`, `gapfill_mz_tolerance_da` | `--gapfill-mz-ppm`, `--gapfill-mz-da` | 10 ppm / 0.002 Da fallback; raw-MS1 measurement radius |
| `gapfill_noise_fraction` | `--gapfill-noise-fraction` | 0.1; imputed height / regional threshold; 0 disables |
| `precursor_max_age_seconds` | `--precursor-max-age-seconds` | 15 s fallback; preceding-MS1 precursor-intensity age limit |
| `gapfill_max_gap_seconds` | `--gapfill-max-gap-seconds` | 15 s fallback; maximum integrated gap-fill interval |
| `mz_tolerance_ppm` | `--mz-ppm` | 10; trace-matching ppm |
| `mz_tolerance_da` | `--mz-da` | 0.002 Da; trace-matching floor |
| `noise_mode` | `--noise` | `adaptive`; alternative `none` |
| `noise_estimator` | `--noise-estimator` | `low_tail`; alternative `clipped` (the 0.2.0 estimator) |
| `noise_fraction` | `--noise-fraction` | 0.10; positive fraction ≤ 1, used by `low_tail` |
| `noise_layout` | `--noise-layout` | `bins`; optional `sliding` |
| `noise_window_da` | `--noise-window-da` | 100 Da bin width; full width in sliding mode |
| `noise_sigma` | `--noise-sigma` | Fixed at 3; other values rejected by API, JSON and CLI |
| `noise_clip_sigma` | `--noise-clip-sigma` | 2; background-clipping multiplier, `clipped` only |
| `noise_iterations` | `--noise-iterations` | 2; integer ≥ 0, `clipped` only |
| `noise_min_neighbors` | `--noise-min-neighbors` | 5; integer ≥ 2; minimum selected low-tail count, or minimum supported count after legacy clipping, excluding the candidate |
| `min_intensity` | `--min-intensity` | 0; nonnegative absolute MS1 floor |
| `min_scans` | `--min-scans` | 5; integer ≥ 2, per detected trace/final segment; seeded forced measurements bypass it |
| `max_gap_scans` | `--max-gap-scans` | 2; integer ≥ 0 |
| `max_gap_seconds` | `--max-gap-seconds` | 15 s fallback; trace gaps only, learned with supported calibration |
| `segment` | `--segment` / `--no-segment` | `true` |
| `peak_shape_guard` | `--peak-shape-guard` / `--no-peak-shape-guard` | `true`; additional evidence for splits in noisy traces |
| `peak_noise_score_limit` | `--peak-noise-score-limit` | 3; nonnegative whole-trace score above which the guard applies |
| `peak_min_snr` | `--peak-min-snr` | 5; nonnegative smoothed peak excess / temporal noise requirement within the guard |
| `peak_min_gaussian_similarity` | `--peak-min-gaussian-similarity` | 0.7; [0,1], symmetric or two-sided local-template correlation within the guard |
| `smoothing_sigma` | `--smoothing-sigma` | 0.6 MS1 samples |
| `peak_prominence` | `--peak-prominence` | 0.1; positive fraction ≤ 1 |
| `peak_distance` | `--peak-distance` | 3 samples; integer ≥ 1 |
| `assume_centroid` | `--assume-centroid` / `--no-assume-centroid` | `false`; affects only absent MS1 mode metadata |
| `profile_mode` | `--profile-mode` | `centroid`; alternative `reject` |
| `profile_min_points` | `--profile-min-points` | 3; integer ≥ 3 |
| `profile_gap_factor` | `--profile-gap-factor` | 4; local spacing multiplier |
| `quant_method` | `--quant-method` | `peak_height`; alternative `peak_area` |
| `associate_ms2` | `--ms2` / `--no-ms2` | `true` |
| `ms2_mz_tolerance_ppm` | `--ms2-mz-ppm` | 20 |
| `ms2_mz_tolerance_da` | `--ms2-mz-da` | 0.01 Da |
| `ms2_rt_tolerance` | `--ms2-rt-tolerance` | 2 s; nonnegative boundary extension |
| `ms2_max_isolation_width` | `--ms2-max-isolation-width` | 0 (disabled); a positive value limits full isolation width in Da |
| `ms2_seeded_features` | `--ms2-seeded-features` / `--no-ms2-seeded-features` | `true`; create unmatched precursor features and measure raw MS1 |
| `ms2_gapfill_rt_window` | `--ms2-gapfill-rt-window` | 5 s on each side of a new seed; positive half-window |
| `best_ms2` | `--best-ms2` | `fragment_intensity`; alternatives `precursor_intensity`, `apex` |
| `group_ions` | `--group-ions` / `--no-group-ions` | `true` |
| `group_isf` | `--isf` / `--no-isf` | `true`; active only during grouping |
| `adduct_ids` | `--adducts` | 1–15, 22, 23; CLI accepts IDs, ranges, exact names |
| `custom_adducts` | `--adduct-file` or config JSON | Empty list; unique validated ion definitions |
| `group_mz_tolerance_ppm` | `--group-mz-ppm` | 10 |
| `group_mz_tolerance_da` | `--group-mz-da` | 0.002 Da |
| `group_rt_tolerance` | `--group-rt-tolerance` | 3 s; nonnegative maximum apex difference |
| `group_min_scans` | `--group-min-scans` | 5; integer ≥ 3 |
| `group_min_correlation` | `--group-min-correlation` | 0.9; positive fraction ≤ 1 |
| `group_min_overlap` | `--group-min-overlap` | 0.5; positive fraction ≤ 1 |
| `isf_min_relative_intensity` | `--isf-min-relative-intensity` | 0.05; fraction from 0 to 1 |
| `isf_max_peaks` | `--isf-max-peaks` | 100; integer ≥ 1 |

All numeric settings must be finite. Unless explicitly allowed to be zero above,
scale/window/time settings are positive. Each pair of Da/ppm tolerances is
nonnegative and must have at least one positive member. Neutral-loss selection
requires its compatible parent ion. `--no-ms2` still permits adduct and
water-loss hypotheses; it cannot supply MS2-supported fragment evidence.
Execution-only options such as input paths, worker count, executor, report
location, retry/resume and overwrite belong to the run summary rather than
the per-file processing configuration.

## 12. Verification and remaining scope

The [validation report](../benchmarks/results/VALIDATION.md) records the
50-file MassIVE corpus, per-file checksums, timings, MS2 reference checks and
three optional grouping examples. Numerical tests cover m/z-dependent noise,
scale changes, sparse spectra, strong-peak interference, profile mass truth,
trace competition, missing scans, splitting/area conservation, polarity,
native IDs, ranking, signed-charge masses and insufficient grouping evidence.
Batch tests deliberately interleave corrupted/missing inputs with valid files,
verify continued execution in all worker modes, and retry only the failures.

The core-engine speed comparison excludes XML, adaptive filtering and optional
grouping; the 50-file timings are host/workload-specific and not million-file
throughput guarantees. Feature counts, filtering fractions and candidate links
are not validated metabolite counts or measured false-discovery rates.
Real profile validation, instrument-specific mass-accuracy studies, dilution
linearity, DIA deconvolution, mobility-aware features and dataset alignment
remain outside the demonstrated scope.
