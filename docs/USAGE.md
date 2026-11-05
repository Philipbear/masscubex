# Usage and output

Core algorithms follow pinned MassCube 1.2.22, with adaptive MS1 noise and
auto-tuning, profile centroiding and MS2-seeded gap filling enabled by default.
See [METHODS](METHODS.md) and [the upstream comparison](UPSTREAM_COMPARISON.md).

The installer selects NumPy **>=1.24,<2.4**. The pinned upstream finalization
calls `np.trapz`, which [NumPy 2.4 removed](https://numpy.org/doc/2.4/release/2.4.0-notes.html#removal-of-deprecated-functions-and-arguments).
The dependency bound preserves those original functions and their arithmetic.

## Optional Numpress decoder

Ordinary float32/float64 mzML and mzXML arrays need no additional decoder.
MS-Numpress-encoded mzML needs `pynumpress`; install it on Python 3.10–3.12 with
`python -m pip install -e '.[numpress]'` from this checkout.

On Python 3.13, `pynumpress` 0.0.9's isolated build pins NumPy 2.0, which predates
Python 3.13 wheels. Install MassCubeX first, then build the decoder against the
compatible NumPy already installed in your virtual environment:

```bash
python -m pip install -e .
python -m pip install "setuptools>=77" wheel Cython
python -m pip install --no-build-isolation "pynumpress>=0.0.9,<0.1"
```

A C/C++ compiler is needed when no decoder wheel is available (for example,
Xcode Command Line Tools on macOS or Microsoft C++ Build Tools on Windows).
The explicit build tools above are needed because pip does not install build
dependencies when [build isolation is disabled](https://pip.pypa.io/en/stable/reference/build-system/#disabling-build-isolation).
CI exercises this route on Python 3.13 on all three platforms. Missing Numpress
support is reported as a file error when an input uses that encoding; other
files in the batch continue.

## Extract features

```bash
masscubex extract sample.mzML -o features
masscubex extract /data/study -o features --workers 3 --resume
masscubex extract --manifest files.txt -o features --workers 3 --resume
masscubex extract sample.mzML -o features --format tsv
masscubex defaults > params.json
masscubex extract sample.mzML -o features --config params.json
masscubex extract --help
```

Parquet with Zstd compression is the default; TSV is optional. Input extensions
are case-insensitive and include gzip. Manifest paths are relative to the manifest.
Output names contain a source-path hash. Each successful input receives a table and JSON sidecar.
Independent workers use processes by default; `--executor thread` is available.
For large batches, run `masscubex warmup` once first. The CLI defaults to one
worker. Each worker has its own scan-cache budget and feature graph.
Final feature/relationship memory grows with the number of features in a
file, so measure representative large files before increasing concurrency.

With multiple process workers, pools drain and restart after `8 * workers`
extraction attempts to release accumulated allocator memory. Resume hits do not
count. Use `--worker-recycle-files N` to adjust this interval, or `0` to disable.
Files are never interrupted for recycling; failures still receive individual
error records. The process coordinator keeps numerical libraries in workers.

## Scientific defaults and extensions

| Setting | Default |
| --- | --- |
| MS1 noise | 100 Da bins, lowest 10%, minimum five background points, μ + 3σ |
| Trace tolerance | Auto-tuned by polarity and m/z region; radius `max(Da, m/z × ppm × 1e−6)`; unsupported file-wide calibration falls back to strict ±0.01 Da |
| Trace gaps | Ten missing observations, no elapsed-time limit |
| Minimum output observations | One: short traces are retained as upstream |
| Near-mass merging | 0.002 Da after filtering |
| Splitting | Upstream noise veto >3, sigma 0.6 scans, prominence 3%, distance three, two passes |
| Quantification | Apex m/z; upstream integer height and area |
| MS2 association | Auto-tuned precursor tolerance/RT allowance when supported; fallback strict ±0.1 Da and adjacent acquired MS1 edges; highest-height candidate |
| Best MS2 | Summed retained fragment intensity; first acquisition wins ties |
| MS2 noise | Fixed instrument default: 10,000 for Orbitrap header keyword, otherwise 500 |
| Ion grouping | On, with the selected catalog |
| Confirmed DIA | Skip before extraction; unknown modes are reported and processed |
| Auto-tuning / profile centroiding / MS2 seeding | On / centroid profile / on |

Select the revised signal-preserving estimator explicitly:

```bash
masscubex extract sample.mzML -o features --noise-estimator consistency
masscubex extract /data/study -o features --noise-estimator consistency --workers 3 --resume
```

It uses μ + 3σ in fixed bins (`--noise-window-da`, default 100), with independent
MS1/MS2 protection. The existing calibration pilot runs first; the estimator
then uses its learned, polarity-specific mass radii. Profile MS1 uses the
conditional profile-aware mass detector. The existing `low_tail` default and
`clipped` alternative remain available. `--noise-fraction` controls `low_tail`,
not the consistency population model; `--noise-layout sliding` is unsupported
with consistency. Single-spectrum `estimate_noise` cannot evaluate this mode;
use the extraction pipeline, which supplies neighboring scans.

Centroid arrays and point decisions use temporary memory maps beside the output;
these are removed on success or a handled exception. Temporary storage scales
with the number of detected/raw centroids, while short background-training
copies are released after their last consumer. Final MS2/grouping still requires
the file's feature summaries, so this is not a hard bound on total worker RSS.
The existing raw-MS1 gap-fill search and its raw-scan regional imputation rule
are retained. [Complete estimator method](NOISE_ESTIMATOR_DESIGN.md).

Disable individual additions when needed:

```bash
masscubex extract sample.mzML -o features --no-auto-tune
masscubex extract sample.mzML -o features --profile-mode reject
masscubex extract sample.mzML -o features --no-ms2-seeded-features
# Upstream parameter baseline, retaining adaptive MS1 noise and the selected catalog
masscubex extract sample.mzML -o features --no-auto-tune --profile-mode reject --no-ms2-seeded-features
```

Unmatched usable MS2 precursors create labelled seeded features, with raw-MS1
measurements or labelled noise-derived heights. Scans lacking usable metadata
remain in JSON diagnostics. With seeding disabled, all unmatched scans are
reported there instead of adding rows. [Full extension details](AUTO_TUNING.md).

Use `--mz-da`, `--mz-ppm`, `--max-gap-scans`, `--min-scans`,
`--centroid-mz-da`, `--ms2-min-intensity`, `--ms2-relative-intensity`,
`--ms2-mz-da`, `--ms2-rt-tolerance`, `--smoothing-sigma`, `--peak-distance`,
`--peak-prominence` and `--segmentation-passes` for explicit overrides.
All parameters are available through JSON and `--print-config`. Explicit mass
overrides are protected from automatic calibration. The removed seconds-based
smoothing/separation options cannot silently change upstream splitting.

## Ion forms

```bash
masscubex adducts
masscubex adducts --json
masscubex extract sample.mzML -o features --adducts 1-15,22,23
masscubex extract sample.mzML -o features --no-group-ions
```

Defaults include IDs 1–15 and the two requested water losses, 22–23. Multiply
charged forms and trimers (16–21) remain disabled. Each strongest ungrouped ion
starts as [M+H]+ or [M-H]-. One assigned form replaces competing pair hypotheses;
monomer H/Na pairs do not get duplicate dimer assignments. Actual dimer ions remain
allowed. `--no-isf` disables best-MS2 fragment matching; remove 22–23 from the
catalog to exclude the water-loss forms. [Catalog/custom ions](ION_GROUPING.md).

## Output columns

Schema 13 has 15 default columns. Parquet with Zstd compression is the default;
use `--format tsv` for text output. m/z has five decimals, times are seconds with two,
intensities have one and correlations have four. Original MS1 scan arrays are
float32; height and area follow upstream integer quantification before formatting.
Nested MS2 references and relationships are native Parquet lists and JSON in TSV.

The default measurement is `peak_height`. These switches apply to both formats:

```bash
# Add area alongside height:
masscubex extract sample.mzML -o features --include-peak-area
# Select area alone:
masscubex extract sample.mzML -o features --quant-method peak_area
# Include detailed ion-group evidence:
masscubex extract sample.mzML -o features --include-ion-relations
# Include area and detailed relationships (17 columns total):
masscubex extract sample.mzML -o features --include-peak-area --include-ion-relations
# Include the two diagnostic columns:
masscubex extract sample.mzML -o features --include-apex-noise-threshold --include-gaussian-similarity
```

The same controls are `include_peak_area`, `include_ion_relations`,
`include_apex_noise_threshold`, `include_gaussian_similarity`, and `quant_method`
in JSON/Python configuration. All four include options default to `false`;
enabling all four gives 19 columns with height as the primary measurement.
Omitting detailed relationships preserves grouping and the compact
adduct/isotope/source-fragment labels. The JSON sidecar's `output_columns`
records the exact per-file layout, also embedded in Parquet metadata. Changing
these options invalidates a previous `--resume` result.

| Column | Meaning |
| --- | --- |
| `feature_id` | Sequential identifier within this file |
| `mz` | Upstream apex mass for detected features; precursor m/z for seeded features |
| `rt_apex_seconds` | Time of its highest original intensity |
| `rt_start_seconds`, `rt_end_seconds` | Integration boundaries |
| `peak_height` | Default measurement; upstream integer detected height; seeds use raw-MS1 height or a labelled noise-derived estimate when status is `imputed_noise`; blank if unavailable; omitted when area alone is selected |
| `peak_area` (optional) | With `--include-peak-area` or `--quant-method peak_area`: upstream integer trapezoidal area, intensity × seconds; seeds use refined bounds or the fallback window |
| `apex_noise_threshold` (optional) | With `--include-apex-noise-threshold`: local MS1 threshold at the detected apex; regional reference for an imputed seed |
| `gaussian_similarity` (optional) | With `--include-gaussian-similarity`: symmetric apex/FWHM Gaussian Pearson score on the original detected trace, 0–1, four decimals; null for MS2-seeded forced measurements |
| `polarity` | `positive`, `negative`, or `unknown` |
| `feature_origin` | `ms1_detected` or `ms2_seeded` |
| `ms1_quant_status` | `detected`, `gap_filled`, `imputed_noise`, `no_ms1_signal`, `no_ms1_scans`, `ambiguous_polarity`, `no_local_ms1_signal`, `no_local_peak`, or `sampling_limit` |
| `ms2_scan_numbers` | List of native scan numbers; `null` where the format has no numeric scan number |
| `best_ms2_scan_number` | Native scan number of the selected best spectrum; blank when unavailable |
| `adduct` | One assigned ion form; ungrouped positive/negative anchors start as `[M+H]+`/`[M-H]-`; blank if grouping is off or polarity unknown |
| `group_id` | Greedy ion-group ID within this file, including singleton groups |
| `isotope_state` | Upstream `M+0`/`M+1`/`M+2`/`M+3` state; blank if grouping is off or polarity unknown |
| `is_in_source_fragment` | Assigned source-fragment flag from the selected best parent MS2; null when checks are disabled or polarity unknown |
| `ion_relations` (optional) | With `--include-ion-relations`: assigned isotope/adduct/fragment links with partner feature IDs and available evidence |


Previously removed redundant measurement/selection and seed precursor/RT
columns remain absent. Their settings or exceptional identifiers are in JSON.
Both height and area remain available internally; these output switches do not
change their calculation or the grouping algorithm. MS2 counts are the length
of `ms2_scan_numbers`; isotope membership is
given by `isotope_state`. These derived columns and the internal segmentation
status are omitted from both formats, as is `scan_count`.

Source scan numbers may be absent or duplicated. For those exceptions,
`ms2_identifier_fallbacks` in the JSON sidecar records the `feature_id`, zero-based
`ms2_list_index` within that row's scan list, `scan_number`, original `native_id`,
and zero-based source `spectrum_index`. Duplicates are checked across all MS2
scans in the file, including unassigned scans. Unique numeric scan numbers need
no extra entry. `best_ms2_identifier_fallbacks` preserves the selected spectrum's
identity when its number is absent or repeated. Full relationship evidence is
available through the optional `ion_relations` column.

## Python

```python
from masscubex import Config, process_file, iter_feature_rows
result = process_file("sample.mzML", "sample.parquet", config=Config())
for row in iter_feature_rows(result.output):
    print(row["feature_id"], row["mz"], row["ms2_scan_numbers"])
```

An explicit extension selects output format. For repeated passes, the decoded
scan cache is bounded by `--scan-cache-mb` (64 MiB per worker); zero disables it.
Processed MS1 data used for grouping is file-backed and closed on every exit.

## Failures and resuming

A corrupt or missing input does not stop other files. Batch reports are saved to
`features/reports/<run-id>/summary.json`, `errors.tsv`, `skipped.tsv` and `files.jsonl`. They
record successful/skipped/error files and reasons. Retry repaired inputs with
`--retry-errors PATH`. A batch with file failures returns nonzero after processing
remaining inputs. Global setup/storage failures can still stop the batch.

Confirmed DIA is skipped by default and listed with detection evidence in
`skipped.tsv`; this does not make the exit code nonzero. `summary.json` separates
DIA exclusions from resume hits. Unknown modes remain explicitly unknown and
are processed. `--dia-policy error` treats DIA as a file error;
`--dia-policy process` explicitly permits it. See [acquisition screening](ACQUISITION.md)
for evidence rules, formats and limitations.

Temporary outputs are published atomically. Resume checks input identity,
configuration, implementation/dependency fingerprints and output checksums.
An exclusive lock protects each output; remove a stale lock only after confirming
its worker stopped. Keep raw files immutable during processing. The new upstream
policy changes the implementation fingerprint, so older results are regenerated
under `--resume` rather than silently reused.
