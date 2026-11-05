# Optional intensity-consistency MS1 noise estimator

Enable with `--noise-estimator consistency`. The default remains the regional
low-tail estimator described in [methods](METHODS.md). Both use μ + 3σ.
This MS1 adaptation is distinct from the inherited MassCube feature algorithms.

## Rationale

The [Spectral Denoising paper](https://pmc.ncbi.nlm.nih.gov/articles/PMC12087458/)
identifies electronic structure through unusually frequent, near-identical
intensities across different m/z values. Its MS/MS formula-based chemical
filter is not applied here. Repeated intensity across unrelated masses can be
background; persistence at the same mass across RT can instead support an ion.

## Implemented algorithm

1. Use unrounded positive centroid intensities. Sort them stably and form fixed
   intensity bands within 0.1% of the band's maximum, without chaining bands.
   A candidate needs at least `max(8, ceil(log2(scan_points)))` members,
   a 50 Da span and support in at least three distinct 10 Da mass cells.
2. Identify supported background populations in fixed 100 Da bins. Separate a
   repeated minimum-intensity atom from a broader background distribution.
   A low continuous density mode is an alternative when repeated bands do not
   identify the full background. Density fits use a smoothed log-intensity
   histogram; the final mean and standard deviation use linear intensities.
   Unsupported populations abstain rather than invent a threshold.
3. Pool background evidence from up to two same-polarity survey scans on each
   side, within 15 seconds. Exclude protected ions from background fitting.
4. Search up to ten same-polarity survey scans on each side, within 15 seconds,
   for physical evidence. Two other scans with a unique mass match protect a
   centroid. Use the normal regional MS1 tolerance, without widening it.
   Original precursor metadata also protects matching MS1 ions near MS/MS events.
5. Reject only identified background candidates at or below their regional
   μ + 3σ ceiling. A repeated intensity above that ceiling remains. Distinct weak
   ions and points with MS1 continuity or MS/MS support remain eligible.

Background consistency alone cannot distinguish every true ion from noise.
Temporal protection can retain persistent contaminants. A score based on
`(intensity - μ) / σ` describes local across-m/z variation, not chromatographic
EIC signal-to-noise. Thresholds are heuristic and do not establish ion identity.

## Profile MS1 handling

For this estimator, dense profile segments receive a conditional detection
stage before the common centroid filter. Segments need at least 64 samples and
95% positive values. Estimate a rolling-median baseline over 31 samples.
Estimate electronic variation from the larger of robust second-difference
noise and baseline-residual noise. Detect apices above baseline + 3σ and
retain shoulders above baseline + σ, requiring three contiguous samples.
Centroid masses are first moments; exported heights preserve the original apex.

Sparse, zero-padded or unsupported profiles do not receive an inferred amplitude
cutoff. They use the physical mass-detection guards, followed by the common
centroid filter. Profile-MS2 deconvolution is not implemented.

## Execution

Numba kernels accelerate intensity grouping, recurrence, precursor lookup and
histogram counting. NumPy/SciPy reductions retain their arithmetic order;
fast-math is disabled. A per-file temporary disk-backed cache holds centroids,
thresholds and keep masks, and is removed on success or error. Feature and
relationship memory can still grow with the number of extracted features.

The source is in `consistency_noise.py`, `consistency_kernels.py`,
`consistency_cache.py` and `profile_noise.py`. Automated tests use generated
signals, strict mass boundaries, independent calculations, temporary-file
cleanup and thread-isolated output layouts. Raw datasets, figures and trial
results are not included in the distribution.
