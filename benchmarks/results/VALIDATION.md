# MassCubeX validation — 2026-11-03

## Version 0.5.1: publication audit

The expanded suite passes **302 tests**. It adds 44 deterministic randomized
checks against independent sorted-window and isolated-ion references; exercises
both file formats, centroid/profile mode, 32/64-bit arrays, compressed/raw arrays
and gzip/plain files; and checks tuned process/thread equivalence, resume,
batched polarity switching and late input failure. The numerical subset also
passes **59 tests with Numba compilation disabled**.
The same **302 tests also pass against the installed wheel**, with the source
path override disabled. See the [release summary](release-v0.5.1.json).

The audit found and fixed four issues:

- TSV floats now use 17 significant digits. Exported mass/time coordinates
  round-trip exactly, preserving inclusive interval membership during raw-data
  reconstruction.
- Non-finite feature values caused by numerical overflow fail the affected
  file instead of publishing a successful table containing `inf`/`nan`.
- Preceding-MS1 precursor-intensity fallbacks are rechecked after all file
  polarities are known. Ambiguous MS1 polarity cannot drive representative-MS2
  ranking in a mixed-polarity file; source metadata intensities are retained.
- Completed refined gap-fill queries release raw slices immediately. Regression
  tests verify early release and actual activation of the dense-query limit.
  The real corpus's peak retained raw slices range from 229 to 1,263 centroid
  pairs per file; optional narrow grouping traces are separate from that count.

The final [independent 50-file audit](MSV000086989-ms2-v0.5.1.json) is tied to the
exact processing implementation fingerprint. It again accounts for all
**125,266 MS2 scans**, with zero unresolved scans, and independently reconstructs
all 4,252 numeric seeded measurements. The 1,540,612 feature rows and seeded
status counts match 0.5.0; TSV checksums change because float precision changed.

Separate [process](batch-errors-process-v0.5.1.json) and
[thread](batch-errors-thread-v0.5.1.json) CLI runs interleaved the 50 real files
with five corrupt/missing inputs. Both completed all 50 good files, reported
exactly the five failures, produced byte-identical good TSVs, cleaned temporary
outputs and skipped all 50 completed tables on resume. Corrupt inputs remained
reported during resume. Thread validation limited BLAS/OpenMP thread pools to
one thread to avoid nested parallelism.

The [CI workflow](../../.github/workflows/tests.yml) tests both the source tree
and installed wheel on Linux, macOS and Windows with Python 3.10–3.13. The raw
corpus is excluded from Git and distributions. This remains one positive-mode
centroid-mzXML cohort; genuine profile/polarity-switching repository data and
biological precursor truth are outside this audit. No finite test suite proves
correctness for every possible acquisition or corrupted file.

## Version 0.5.0: per-polarity calibration and local gap filling

Automatic calibration and ion relationships are enabled by default. The pilot
starts at ±max(0.05 Da, 100 ppm), screens isolated, well-sampled strong peaks,
and needs at least 20 accepted representatives per polarity. MS2 tolerance
needs separate support from at least 20 distinct representatives with MS2.
Physical-time smoothing, peak separation, time tolerances and local raw-MS1
gap-fill bounds are new. The MS1 noise rule is unchanged: 50 Da, lowest 10%,
minimum 5 neighbors, μ + 3σ. See [the complete method](../../docs/AUTO_TUNING.md).

**218 tests passed; Ruff passed.** New tests cover representative quality,
outliers and quantized masses, insufficient evidence, deterministic selection,
protected overrides, time-based segmentation, and interleaved polarity with
different scan cadence, mass jitter, widths and intensities. Known-truth
centroid/profile mzML/mzXML tests verify MS2 association and original-sample
height/area per polarity. Missing polarity, explicit parent inference and
conflicting metadata are tested. Raw below-noise seeds, nearby interfering
ions, distinct chromatographic peaks and repeated-query merging are covered.

The real-data audit exposed a boundary-slice defect: a recovered center near
the broad mass-search edge could lack part of its narrow EIC. Queries now retain
extra mass context while centers remain restricted to the original search.
A regression test reconstructs the full EIC across that boundary.

The final [50-file report](MSV000086989-ms2-v0.5.0.json) verifies original source
checksums and all MS2 identifiers/representative scores through independent
Pyteomics decoding. All 50 files supported both MS1 and MS2 calibration, with
76–100 accepted MS1 representatives. The resulting tables contain 1,540,612 rows.

| Outcome | Count |
| --- | ---: |
| Input MS2 / assigned MS2 | **125,266 / 125,266** |
| MS2 matched to independently detected features | 120,123 |
| MS2 assigned to seeded features | 5,143 |
| Unresolved MS2 | **0** |
| Seeded rows with numeric raw-MS1 measurements | 4,252 |
| Seeded rows with ambiguous mass | 720 |
| Seeded rows without local MS1 signal near the seed | 104 |
| Seeded rows outside recovered local peak support | 67 |

Every numeric seeded height, area, positive scan count and status was checked
against independently reconstructed raw EICs, using the exported mass radius,
bounds and effective per-polarity interval limit. The other 891 seeded rows keep
all MS2 references with explicit status and blank quantification. The union of
assigned indices and unresolved records equals all input MS2 for each file.
Default grouping reported 70,042 adduct, 13,075 water-loss and 19,817 generic
in-source-fragment candidate pairs; these are hypotheses, not validated identities.

In this cohort, learned ppm floors were 5 for both MS1/MS2; Da floors ranged
0.0005615–0.0011902. Learned gap-fill search half-windows ranged 11.67–27.29 s,
while integration uses narrower recovered local bounds. These values are
observations from this dataset, not proposed universal constants. Strong-peak
repeatability does not establish absolute accuracy or weak-ion sensitivity.

The four-worker audit took 91.68 s, including independent reference checks.
Reported per-file processing time had median 4.96 s (range 4.21–10.01 s), with
median calibration time 0.46 s. These concurrent-run timings include ion
grouping and cannot be compared directly with older configurations or treated
as a universal million-file throughput estimate.

The real cohort remains positive-mode centroid mzXML. **Polarity-switching and
profile validation here are synthetic**, with known expected signals; genuine
switching LCMS files still need external validation. Complete MS2 bookkeeping
does not establish biological precursor identity, DIA purity or deconvolution.
Version 0.4.0 results below are retained as historical results with different
parameters, grouping defaults and gap-fill behavior.

## Version 0.4.0: complete usable-MS2 coverage and forced MS1 measurements

The isolation-width cutoff defaults to **0 (disabled)**. Unmatched usable
single-precursor MS2 now seed features, measured from raw MS1 within ±5 seconds.
Existing detected-feature assignments are preserved. The new origin and
quantification-status columns separate detections, forced measurements, measured
absence and unavailable MS1 data. Missing precursor metadata and multiplexed
precursors remain explicit per-scan unresolved records in the JSON sidecar.
MS2-only files with usable precursors produce rows with blank MS1 quantification.

**193 automated tests passed.** New coverage includes wide isolation windows,
centroid/profile mzML and mzXML, gzip, below-threshold and short MS1 traces,
repeated MS2, duplicate/missing native scan numbers, polarity switching,
missing metadata, zero versus unavailable signal, MS2-only input, float32 mass
boundaries, duplicate centroid normalization, exact areas, long acquisition
gaps, configuration overrides, resume, second-pass failures, and identical
seeded TSVs from process/thread workers.

The [50-file audit](MSV000086989-ms2-v0.4.0.json), produced by
[validate_ms2_coverage.py](../validate_ms2_coverage.py), used the same downloaded
MSV000086989 sample and verified source SHA-256 values:

| Outcome | MS2 scans |
| --- | ---: |
| Total input MS2 | 125,266 |
| Matched independently detected features | 122,032 (97.4183%) |
| Recovered through new seeded features | 3,234 (2.5817%) |
| Final assigned scans | **125,266 (100%)** |
| Final unresolved scans | **0** |

Initially, 1,187 scans had no feature in the mass tolerance and 2,047 had no
matching feature interval. All 3,234 added features had measurable raw MS1
signal. The final tables contain 1,714,388 rows, including these 3,234 additions.

Pyteomics independently decoded all 125,266 MS2 scans. The audit checked every
assigned native ID, numeric number, source index, precursor mass/time/polarity,
and representative-score selection, with no scan assigned to multiple rows.
It independently reconstructed all 3,234 seeded EICs, using float64 m/z and
summing exact duplicate masses before testing interval membership, and checked
height, area, positive scan count and status against the exported table.
The union of assigned source indices and unresolved records equaled all input
MS2 indices for every file.

This corpus contains centroid mzXML with recorded isolation widths no greater
than 4 Da. It therefore validates the coverage/gap-fill path, while synthetic
mzML/mzXML tests demonstrate that the disabled cutoff retains wider scans.
Neither the 100% mapping result nor raw-MS1 gap filling proves biological
precursor identity or DIA purity. Seed windows may overlap and do not necessarily
capture an entire chromatographic peak. Real profile-DIA validation remains
outside this audit. The report's wall time includes extraction **and independent
reference checks**; it is not an extraction-throughput benchmark.

## Version 0.3.2: lowest 10%, detection multiplier fixed at 3

Defaults are 50 Da, lowest 10%, minimum 5 positive neighbors. Both API and CLI
reject detection multipliers other than 3, including values loaded from JSON.
An explicit `--noise-sigma 3` can override an old incompatible JSON value.
The low-tail rank-tree algorithm and strict retention comparison are unchanged.

**166 automated tests passed.** Coverage includes the new fraction, exact
minimum and ceiling rules, independent sorted-window numerical agreement,
fixed-multiplier validation, explicit CLI overrides, profile mass detection,
MS2 references, corruption handling and process/thread execution.

The [three-scan comparison](noise-fractions-v0.3.2.json) compares 5% and 10%
at the same 50 Da / minimum 5 / 3σ settings, on a sample, blank and QC scan.
Within m/z 100–250, retained counts are respectively 120→120, 275→254 and
355→332. [Figures and all centroid values](../figures/README.md) are available.
These are centroid-filter sensitivity checks, not identified analytes or
labeled noise-rejection accuracy. The archived 0.3.1 comparisons with 4σ/5σ
were exploratory and are not current defaults or supported runtime options.
The full 50-file timing results below remain historical measurements of their
stated configurations; no throughput claim is transferred to this version.

## Historical version 0.3.1: updated defaults

Default MS1 noise settings changed to a **50 Da window**, **lowest 1%**, and
**minimum 5 positive neighbors**, with the same subset mean + 3σ rule.
The same-scan [figure and per-centroid data](../figures/README.md) were regenerated:
m/z 100–250 retains 120/148 centroids and 96.9% of summed intensity, versus
83/148 and 91.2% under 0.3.0. All displayed windows use the minimum of five.
This is a sensitivity illustration without ground-truth noise labels.

The updated automated suite has **162 passing tests**, including the inclusive
±25 Da boundary, exact five-neighbor support, four-neighbor fallback, minimum
expansion, ceiling at 501 neighbors, CLI defaults/overrides, profile centroiding,
MS2 bookkeeping and corrupt-file continuation. The full 50-file timings below
were measured under the stated older settings; they are not 0.3.1 throughput claims.

## Historical version 0.3.0: lowest-intensity background subset

The same 50 complete files were rerun with the new default estimator:
`noise_estimator=low_tail`, lowest 5% of positive neighbors, expanded to at
least 8, followed by the mean plus 3 population standard deviations **of that
subset**. All other extraction defaults are unchanged. Exact formulas and
implementation are in [Methods §3](../../docs/METHODS.md#3-new-ms1-algorithm-moving-low-intensity-background-across-mz).
The [same-scan figure](../figures/README.md) includes the effective subset sizes.

| Measurement | Previous clipped default (0.2.0) | New low-tail default (0.3.0) |
| --- | ---: | ---: |
| Retained MS1 centroids / 34,905,337 | 7,028,048 (20.1%) | 27,040,414 (77.5%) |
| Chromatographic features | 308,303 | 1,630,998 |
| Associated MS2 / 125,266 | 100,232 | 121,831 |
| Features with MS2 | 76,634 | 92,221 |
| Sequential, 50 files | 35.43 s | 100.15 s |
| Four processes, 50 files | 10.78 s | 28.85 s |
| Four threads, 50 files | 55.00 s | 147.00 s |
| Resume: 50 verified and skipped | 0.076 s | 0.201 s |

These are separate warm runs on the same Apple Silicon host, not a randomized
speed comparison. Inputs were local; initial parent imports/JIT and downloads
are excluded, worker startup is included. No tests or other numerical workloads
ran alongside the new extraction timings. Median file time was 1.890 s; parent
peak RSS was 577.2 MB, including libraries, validation and reference checks,
excluding child processes. This is not per-file RSS or a million-file guarantee.

The new rule is substantially more permissive. More retained centroids, a more
expensive rank-selection kernel, more trace/segmentation work and larger TSVs
contribute to higher runtime. The low-tail standard deviation is **not** a
full-background standard deviation. More features/assignments are not evidence
of better accuracy; many extra features can be background. The synthetic
uniform-noise fixture explicitly demonstrates retention of more than 80% of
noise points. The previous three-analyte/no-extra-feature fixture remains a
regression test of the optional `clipped` method, not a claim for `low_tail`.

All 50 sequential/process/thread TSVs were byte-identical. All outputs passed
resume checksum checks. Independent Pyteomics decoding verified every one of
125,266 MS2 scans and all 121,831 assignments, including native IDs, source
indices, summed exported intensity and best-scan tie rules. The first five MS1
scans in each file (250 total) had identical decoded arrays and times. As before,
all public exports were centroided; real profile data remain unvalidated.

**162 automated tests passed** after the estimator change. Added coverage
compares the compiled rank tree against independently sorted neighborhoods,
including fractions from 0.01 to 1, physical window boundaries, ties, zero
padding, candidate exclusion, ceiling/minimum rules, unsupported windows,
intensity scaling and high-baseline variance precision. Existing profile,
MS2, grouping, resume, corruption and worker-recovery checks also pass.

Representative feature counts, without ground-truth labels:

| File | Low-tail adaptive | Previous clipped adaptive | Fixed 500 | Fixed 30,000 |
| --- | ---: | ---: | ---: | ---: |
| NIST_POS_Samp_01-01.mzXML | 26,144 | 4,714 | 33,057 | 15,427 |
| Blank_00.mzXML | 37,689 | 6,359 | 45,721 | 13,561 |
| QC_Mix_01.mzXML | 40,781 | 7,953 | 50,401 | 13,621 |

The separate [noise-kernel benchmark](noise-v0.3.0.json), after JIT warmup,
measured median low-tail costs of 0.235 ms, 3.062 ms and 40.089 ms for 1,000,
10,000 and 100,000 synthetic points at a 25 Da window (five repetitions).
Widening the 100,000-point window to 500 Da took 42.538 ms: the implementation
does not sort each enlarged neighborhood. The legacy clipped kernel took
2.296 ms for 100,000 points / 25 Da and is cheaper, but selects a different
subset. These timings include validation, sorting and temporary allocation;
they exclude XML, centroiding, chromatographic work and I/O. Reproduce with
`python benchmarks/benchmark_noise.py`; these inputs do not provide noise truth.

Optional grouping was rerun on the sample, blank and QC. All original feature
measurements and MS2 references stayed unchanged; relationships were symmetric
and every claimed fragment was verified in independently decoded spectra.

| File | Extraction + grouping | Adduct candidates | Water-loss candidates | Other fragment candidates |
| --- | ---: | ---: | ---: | ---: |
| NIST_POS_Samp_01-01.mzXML | 2.826 s | 1,478 | 290 | 409 |
| Blank_00.mzXML | 3.958 s | 1,317 | 163 | 323 |
| QC_Mix_01.mzXML | 4.244 s | 954 | 113 | 243 |

The [grouping report](ion-grouping-v0.3.0.json) records these unconfirmed candidates.
A fresh [mixed-batch check](batch-errors-v0.3.0.json) then interleaved the 50
valid files with the same five corrupt/missing inputs. All 50 valid TSVs matched
the new reference hashes, all five failures were logged, and the final status
was `completed_with_errors` with exit code 1. This CLI run took 29.84 s with
four processes, including startup and reports. No partial files or locks remained.

The 0.3.0 extraction outputs are under `validation/MSV000086989/v0.3.0/`.
The [full report](MSV000086989-v0.3.0.json) records input/output hashes,
configuration, implementation fingerprint, dependency versions and per-file
counts/timings. Run `python benchmarks/validate_repository.py --reference-check`
to regenerate a report for the installed source version.

## Historical 0.2.0 batch continuation and persistent error reports

After adding automatic reports, a fresh CLI run interleaved the same **50 valid
MassIVE files** with **five intentionally invalid inputs**: truncated mzML,
truncated mzXML, invalid base64 mzXML, truncated gzip mzML, and a missing file.
All 50 valid files completed and their TSVs matched the reference checksums.
All five bad inputs appeared in the persistent `errors.tsv` and `files.jsonl`;
`summary.json` reported `completed_with_errors`, 50 successes and 5 failures.
No partial TSVs or lock files remained. Exit code 1 was returned after all
55 inputs were handled.

This mixed batch took **12.24 s with four processes**, including CLI startup
and report writes on the same host with cached kernels. It is a different
timing boundary from the extraction-only measurements below. The
[machine-readable check](batch-errors-v0.2.0.json) records the outcomes and
location of the full local reports. Reproduce with:

```bash
python benchmarks/validate_batch_errors.py
```

Automated tests additionally exercise corrupted/nonfinite arrays and
continuation in sequential, process and thread modes, retry-only processing
after repairing failures, preservation of previous reports, checkpointing on
interruption, inaccessible-directory/manifest errors, and recovery after an
actual child process exits abruptly. **134 tests passed** across the suite at that stage. Reports use bounded-memory streaming rather than accumulating failures.

## Version 0.2.0: MS/MS references, ranking and ion relationships

The same 50 complete files were rerun with the 0.2.0 defaults: peak-height
quantification, MS2 association enabled, best MS2 selected by the sum of all
exported MS2 intensities, and optional ion grouping disabled. The run contains
**125,266 MS2 scans**, of which **100,232** were assigned to **76,634 features**.
All scans in this sample had narrow, single-precursor metadata; the remaining
25,034 MS2 did not match a retained feature under the mass/time/polarity rules.
MS1 results remain **308,303 features**, with the same filtering counts as 0.1.0.

| Execution | Files | Wall time |
| --- | ---: | ---: |
| Sequential | 50 | 35.43 s |
| Four processes | 50 | 10.78 s |
| Four threads | 50 | 55.00 s |
| Resume, all outputs verified and skipped | 50 | 0.076 s |

Median file time was **0.691 s**. All sequential/process/thread TSVs were
byte-identical. These are single warm measurements on the host described
below; downloads and initial parent JIT/imports are excluded, worker startup
is included. The test suite did not run alongside these timings. Parent peak
RSS was 376.5 MB, including validation libraries and reference checks, excluding
child processes. It is not a per-file memory figure or a controlled memory
comparison with the historical run.

Independent Pyteomics 4.6.3 decoding checked **all 125,266 MS2 scans**: native
identifiers, zero-based spectrum indices, RT, precursor mass/intensity/reference,
and the sum of the complete exported intensity array. Every reported scan
reference and representative choice in all 50 TSVs was then checked against
those spectra. All 100,232 assignments had consistent mass, time and polarity;
each source MS2 was assigned at most once. Representative scores matched the
largest full-spectrum sum within each feature, with the documented tie rules.
An additional 250 MS1 spectra had equal arrays and times in both readers.
This validates bookkeeping and numerical consistency, not precursor purity
or the biological correctness of feature assignments.

Optional grouping was exercised separately with the user-selected ions 1–15,
22 and 23. Each proposed link was symmetric, and every claimed MS2 product
mass/intensity was verified in independently decoded source spectra. Grouping
preserved all MS1 columns, MS2 references and representative selections.

| File | Extraction + grouping | Adduct candidates | Water-loss candidates | Other MS2-supported fragment candidates |
| --- | ---: | ---: | ---: | ---: |
| NIST_POS_Samp_01-01.mzXML | 0.770 s | 409 | 49 | 220 |
| Blank_00.mzXML | 0.847 s | 323 | 31 | 130 |
| QC_Mix_01.mzXML | 0.914 s | 278 | 13 | 93 |

The 346 distinct evidence MS2 spectra support 487 candidate links (including
some water losses). Candidates in blanks are expected to include background
chemistry. These counts are not confirmed compounds, true-positive rates or
measured source-fragment identities. Unit tests cover exact negative/positive
ion masses, water-loss chains, insufficient/constant/uncorrelated traces,
unknown native scan numbers, missing precursor intensity, configuration and
quantification selection. **125 automated tests passed.**

The 0.2.0 outputs are under `validation/MSV000086989/v0.2.0/`: `features/`,
`process_features/`, `thread_features/`, and `grouped/` for the three optional
grouping examples. Reports:
[full run](MSV000086989-v0.2.0.json) and
[optional grouping](ion-grouping-v0.2.0.json).

```bash
python benchmarks/validate_repository.py --reference-check
python benchmarks/validate_ion_grouping.py
```

## Historical version 0.1.0: MS1-only baseline

The following results were recorded before MS2 association and ion
relationships were added. They are retained to distinguish the earlier
workload from the current defaults.

The requested [MassIVE dataset MSV000086989](https://massive.ucsd.edu/ProteoSAFe/dataset.jsp?task=2f73277b4e034948acebfdf1edab17ed) was tested using
**50 complete public mzXML files**: 40 NIST sample injections, six blanks and
four QC injections. Selection is deterministic and spans the available filenames;
it does not represent 50 independent biological specimens.

All 50 files processed successfully. The sample contains **604.1 MB**,
**55,198 MS1 scans**, and **34,905,337 input centroids**.
Outputs contain **308,303 chromatographic features**. Adaptive filtering
retained 7,028,048 centroids and removed 79.87%
of input points. This fraction is a filtering rate, not a measured noise-removal
accuracy or false discovery rate. Feature counts are not metabolite counts.

## Timing and determinism

Measured on `macOS-26.6.2-arm64-arm-64bit`, Python 3.10.12. Input files were already
local. Parent imports and initial JIT compilation are excluded; process-worker
startup is included. This is one warm run per configuration, not a hardware-
independent throughput guarantee. No test suite ran alongside the reported run.

| Execution | Files | Wall time |
| --- | ---: | ---: |
| Sequential | 50 | 23.72 s |
| Four processes | 50 | 8.12 s |
| Four threads | 50 | 39.05 s |
| Resume, all outputs verified and skipped | 50 | 0.034 s |

Median single-file extraction time: **0.460 s**.
All process/thread TSVs match their sequential counterparts byte for byte.
Processes are the default because threads were slower on this workload, despite
GIL-free numerical kernels. XML handling and Python feature finalization remain
part of the workload.

The parent validation process reached 489.7 MB peak RSS,
including Python, Numba, scientific libraries, thread execution and reference
checks. This does not include child-process RSS and is not a per-file memory
measurement. Per-file live trace-node high-water marks ranged from
16,866 to 47,181;
completed traces are recycled.

## Correctness checks

- The first five MS1 spectra in every file were compared with Pyteomics 4.6.3:
  **250 spectra**, equal decoded m/z and intensity arrays, equal RT after unit conversion.
- All 50 outputs passed SHA-256 verification and were skipped on resume.
- Unit/integration coverage includes mzML and mzXML, both binary precisions,
  zlib, gzip, six MS-Numpress variants, nested scans, time units, malformed data,
  mass matching, gaps, polarity switching, integration, segmentation, atomic
  writes, restart validation and process/thread execution.
- A synthetic chromatographic test embeds three analytes in three different
  m/z-dependent noise bands. The historical clipped filter recovers all three without
  extra features in that fixture; a fixed 30,000 cutoff recovers none.
- Profile tests cover Gaussian masses on nonuniform grids (under 0.1 ppm
  centroid error in those fixtures), resolved neighbors, flat tops, sampling
  discontinuities, one-bin spikes, and profile-to-feature extraction.

**All selected MassIVE exports were already centroided.** They do not validate
profile centroiding on real instrument profiles. Undersampled, saturated or
unresolved mass peaks need further instrument-specific validation. Neither the
synthetic results nor the public-file feature counts establish quantitative
biological accuracy; use known standards, blanks and dilution series for that.

## Fixed-cutoff comparisons

All other settings were held constant. The following are representative file
feature counts, not labeled true/false positives:

| File | Adaptive | Fixed 500 | Fixed 30,000 |
| --- | ---: | ---: | ---: |
| NIST_POS_Samp_01-01.mzXML | 4,714 | 33,057 | 15,427 |
| Blank_00.mzXML | 6,359 | 45,721 | 13,561 |
| QC_Mix_01.mzXML | 7,953 | 50,401 | 13,621 |

## Reproduce

```bash
python -m pip install -e '.[dev,numpress]' pyteomics==4.6.3
python benchmarks/download_massive.py --directory validation/MSV000086989
# The validation script now runs the current version and writes versioned results.
python benchmarks/validate_repository.py --reference-check
python benchmarks/benchmark_engine.py --scans 500 --peaks 500 --repeats 3
```

The downloader uses anonymous **encrypted FTP**, which the service now requires.
The converted mzML directory denied anonymous access during this session, so
original public mzXML exports were used. Downloaded data are excluded from Git.

Per-file outputs are in `validation/MSV000086989/features/`. The complete public
manifest, input checksums, parameters, per-file metrics and output checksums are
in [MSV000086989.json](MSV000086989.json). That report also records the source and
dependency fingerprint used by resume validation.

The separate [engine benchmark](engine.json) measured **16.1× median speedup**
over cloned upstream commit `356748ede691862b9a87f67812db5e728a703b9d` on identical
500-scan × 500-centroid inputs. Both engines returned 500 features. This isolates
trace detection and summary; it excludes XML, adaptive filtering and segmentation,
and does not assert feature equivalence under all upstream processing settings.
