# Compact output and original MassCube comparison — 0.9.0

Tested on the same 50 positive-mode, centroid mzXML files from MSV000086989:
40 biological samples, six blanks and four pooled QC files. Raw data and full
tables remain under `validation/`; compact reports are linked below.

## Output changes

Schema 6 has 27 columns instead of 32. Removed `intensity`, `quant_method`,
`best_ms2_selection_method`, `best_ms2_score` and `apex_threshold_ratio`.
Both height and area remain. Quantification preference and best-MS2 ranking
are file-level parameters; exceptional nearest-apex fallbacks are recorded
by feature ID in the JSON sidecar. All native MS2 IDs, scan numbers, spectrum
indices, best-scan references and ion evidence remain.

Both TSV and Parquet round mass coordinates/radii to five decimal places;
RTs in seconds and intensities/areas/thresholds to two. Nested evidence follows
the same mass/intensity policy, with two decimals for ppm error and six for
correlation. Calculations remain float64. Parquet retains float64 storage for
the rounded values; float32 would not preserve five decimals across the mass
range. Rounding does not alter scan association, integration or peak detection.

| All 50 files | Before | After | Reduction |
| --- | ---: | ---: | ---: |
| TSV | 431.85 MB | 251.27 MB | 41.8% |
| Parquet, Zstd | 85.71 MB | 54.96 MB | 35.9% |

Sizes use decimal MB and exclude JSON sidecars. The example
`NIST_POS_Samp_01-01.mzXML` retains 24,980 rows: its TSV shrinks from 7,350,371
to 4,400,441 bytes; Parquet shrinks from 1,459,791 to 965,394 bytes.

## Feature counts: the concern is real

| 50-file total | Original MassCube | MassCubeX default | MassCubeX with a fixed 50,000 MS1 floor |
| --- | ---: | ---: | ---: |
| All output rows | 1,446,350 | 1,540,612 | 551,449 |
| Detected MS1 features with ≥5 nonzero observations | 391,188 | 1,535,469 | 544,514 |
| Added MS2-seeded rows | 0 | 5,143 | 6,935 |

The total row count increases by only 6.5%, but original MassCube writes many
short traces that MassCubeX excludes. Among detected features with at least
five nonzero MS1 observations, MassCubeX reports **3.93× as many**. This common
five-observation diagnostic is not the upstream single-file filter: that path
keeps short traces. Its later alignment count rule is `total_scans > 5 OR MS2
present`, followed by additional alignment filters.

The increase is concentrated in weak signals and in blanks/QCs:

| Group | Files | Original rows | MassCubeX rows | MassCubeX rows with height <50,000 |
| --- | ---: | ---: | ---: | ---: |
| Biological | 40 | 1,361,921 | 1,074,484 | 469,275 |
| Blank | 6 | 44,035 | 281,431 | 224,978 |
| QC | 4 | 40,394 | 184,697 | 140,757 |

These counts do not label individual features as false positives. Blanks can
contain real contaminants, and a fixed intensity floor can discard real weak
ions. Nevertheless, the concentration in weak blank/QC features warrants
review before repository-scale production.

## What explains the difference?

1. **Noise filtering is the largest contributor.** Changing only MassCubeX's
   MS1 noise policy to a fixed 50,000 floor reduces detected features from
   1,535,469 to 544,514 (64.5%). Learned tolerances and timing estimates were
   verified unchanged. The default low-tail filter retains 27,960,819 MS1
   centroids; the fixed-floor variant retains 9,575,398. Original preprocessing
   retains 9,575,050 after its adjacent-centroid merging. The lowest 10% is a
   truncated intensity distribution; its mean + 3σ is not a three-sigma cutoff
   for the entire background population.
2. **Tracing and splitting both contribute after matching the floor.** Before
   segmentation, original MassCube has 284,427 traces with ≥5 observations;
   MassCubeX has 327,158 (+15.0%). After segmentation, these become 391,188 and
   544,514 (+39.2%). Original splitting skips noisy traces and requires maxima
   above the trace median. MassCubeX uses different prominence/boundary logic,
   learned timing and no equivalent noise-score/median gate. This experiment
   separates the stages; it does not isolate every individual matching rule.
   In blanks, fixed-floor MassCubeX grows from 13,773 unsplit traces to 46,462
   detected peaks, versus original's 10,349 to 14,283.
3. **MS2 seeding is a small contribution.** It adds 5,143 rows, 0.33% of the
   default table. It cannot explain the main count increase.

Recommended next work: validate a background estimate that does not mistake
the truncated low tail for the full noise population, and require stronger
chromatographic evidence before splitting weak/background traces. Use blanks,
QC reproducibility and representative raw EICs to assess sensitivity and false
splitting. A universal 50,000 cutoff is only a diagnostic control here, not a
new recommendation for all instruments. Version 0.9.0 changes output only;
the user-selected 50 Da / lowest 10% / μ + 3σ defaults remain in effect.

## Meaning of `ambiguous_mass`

During refined MS2-seeded gap filling, multiple nearby MS1 mass clusters can
fit the precursor search. Candidates must have local time evidence; candidates
within the learned precursor tolerance are preferred. If the strongest
eligible candidate has less than three times the runner-up's local maximum
intensity, the algorithm cannot confidently choose one. It preserves the
MS2 references but leaves MS1 height and area missing. This is an ambiguous
assignment, not an invalid m/z or mass-precision error. There are 720 such
seeded rows in these 50 files. See [the full refinement rule](../../docs/AUTO_TUNING.md).

## Baseline and reproducibility

Original source is MassCube 1.2.22 at
`356748ede691862b9a87f67812db5e728a703b9d`, the pinned original checkout.
Every Python source file in the audit copy was compared byte-for-byte with
that Git commit; hashes are included in the report. Its Orbitrap defaults
use an MS1 floor of **50,000** and MS2 floor of 10,000, 0.002 Da adjacent-centroid
merging, 0.01 Da tracing tolerance, ten missing scans and two segmentation passes.

This pinned version explicitly rejects mzXML. The audit uses Pyteomics to
decode the original mzXML and translates only metadata into the upstream
mzML scan-extractor schema. Arrays and RTs in minutes pass unchanged to that
extractor; preprocessing, detection, segmentation, finalization and MS2
assignment are original code. This is an audit through an input adapter, not
a claim that the upstream command accepts mzXML. It covers this cohort only.

Reproduction scripts:

- [Original count audit](../audit_upstream_features.py): `--upstream`, `--input`, `--output`.
- [Count comparison and fixed-floor/splitting controls](../compare_upstream_counts.py): `--input`, `--upstream`, `--upstream-source`, `--candidate`, `--floor-output`, `--unsegmented-output`, `--report`.
- [Compact-table parity](../validate_compact_output.py): `--reference`, `--candidate`, `--tsv`, `--reference-tsv`, `--report`.
- [Independent raw-data audit](../validate_ms2_coverage.py): `--directory`, `--workers`, `--output`. Audits retain unrounded writer inputs privately so RT/mass rounding cannot change reconstructed interval membership.

The original-only `IsoSpecPy` dependency was installed into an isolated
validation directory; it is not a MassCubeX runtime dependency.

## Validation evidence

- [Compact-output report](compact-output-v0.9.0.json): all 1,540,612 rows equal
  their v0.8 full-precision baseline after only the specified rounding and
  column removals. TSV and Parquet agree; scientific sidecar fields and
  205,868 ion-relationship entries are preserved.
- [MS2/raw-MS1 audit](MSV000086989-ms2-v0.9.0.json): all 125,266 MS2 scans are
  assigned exactly once. All 4,252 numeric seeded measurements were checked
  against independently decoded raw data using unrounded coordinates;
  unavailable measurements remain blank. Their ambiguity/status rules are
  covered by synthetic tests, not independently established biological truth.
- [Original per-file results](upstream-orbitrap-counts-v0.9.0.json) and
  [full comparison](upstream-comparison-v0.9.0.json) include file checksums and
  counts for every input, group and diagnostic control.
- 551 source tests and 551 installed-wheel tests pass, including precision, nulls, large integer IDs,
  boundary-sensitive gap filling and atomic/resumable output.

A [single warm serial Parquet pass](compact-parquet-run-v0.9.0.json) took 64.60 seconds for all 50 files on the
Apple M2 Max; table writing accounted for 12.53 seconds. This is a validation
timing, not a repeated cross-version throughput benchmark. It excludes input
downloads and independent reference checks.
