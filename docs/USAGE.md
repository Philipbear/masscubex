# Usage and output

Core algorithms follow pinned MassCube 1.2.22, with adaptive MS1 noise and
auto-tuning, profile centroiding and MS2-seeded gap filling enabled by default.
See [METHODS](METHODS.md) and [the upstream comparison](UPSTREAM_COMPARISON.md).

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
worker; the examples use three processes as a starting point for this 32 GiB
Mac. Four processes gave the highest tested sustained throughput below 3 GiB on
the current cohort; see the [speed/RAM measurements](../benchmarks/results/OVERNIGHT_OPTIMIZATION.md).
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
| Auto-tuning / profile centroiding / MS2 seeding | On / centroid profile / on |

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

Schema 10 retains 25 columns. m/z has five decimals, times are seconds with two,
intensities have one and correlations have four. Original MS1 scan arrays are
float32; height and area follow upstream integer quantification before formatting.
Nested MS2 references and relationships are native Parquet lists and JSON in TSV.

| Column | Meaning |
| --- | --- |
| `feature_id` | Sequential identifier within this file |
| `mz` | Upstream apex mass for detected features; precursor m/z for seeded features |
| `rt_apex_seconds` | Time of its highest original intensity |
| `rt_start_seconds`, `rt_end_seconds` | Integration boundaries |
| `peak_height` | Upstream integer detected height; seeds use raw-MS1 height or a labelled noise-derived estimate when status is `imputed_noise`; blank if unavailable |
| `peak_area` | Upstream integer trapezoidal area, intensity × seconds; seeds use refined bounds or the fallback window |
| `scan_count` | Number of observed nonzero MS1 samples |
| `apex_noise_threshold` | Local MS1 threshold at the detected apex; regional reference for an imputed seed |
| `gaussian_similarity` | Symmetric apex/FWHM Gaussian Pearson score on the original detected trace, 0–1, four decimals; null for MS2-seeded forced measurements |
| `segmentation_status` | `split`, `unsplit`, `guarded_split`, `guarded_unsplit`, `disabled`, or `not_applicable` for seeds; `guarded_` means the upstream noise cutoff prevented further splitting |
| `polarity` | `positive`, `negative`, or `unknown` |
| `feature_origin` | `ms1_detected` or `ms2_seeded` |
| `ms1_quant_status` | `detected`, `gap_filled`, `imputed_noise`, `no_ms1_signal`, `no_ms1_scans`, `ambiguous_polarity`, `no_local_ms1_signal`, `no_local_peak`, or `sampling_limit` |
| `ms2_count` | Number of associated MS2 scans |
| `ms2_scan_numbers` | List of native scan numbers; `null` where the format has no numeric scan number |
| `ms2_native_ids` | List of original spectrum identifiers |
| `ms2_spectrum_indices` | List of zero-based positions among all source spectra |
| `best_ms2_scan_number` | Native scan number of the selected best spectrum; blank when unavailable |
| `adduct` | One assigned ion form; ungrouped positive/negative anchors start as `[M+H]+`/`[M-H]-`; blank if grouping is off or polarity unknown |
| `group_id` | Greedy ion-group ID within this file, including singleton groups |
| `is_isotope`, `isotope_state` | Upstream isotope flag and `M+0`/`M+1`/`M+2`/`M+3` state; blank if grouping is off or polarity unknown |
| `is_in_source_fragment` | Assigned source-fragment flag from the selected best parent MS2; null when checks are disabled or polarity unknown |
| `ion_relations` | Assigned isotope/adduct/fragment links with partner feature IDs and available evidence |


Previously removed redundant measurement/selection and seed precursor/RT
columns remain absent. Their settings or exceptional identifiers are in JSON.
`quant_method` remains a file-level preference; both height and area are always
exported. Source scan numbers may be absent or duplicated; native IDs and source
spectrum indices allow reliable lookup. Best-scan exceptions are recorded in
`best_ms2_identifier_fallbacks` in the sidecar.

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
`features/reports/<run-id>/summary.json`, `errors.tsv` and `files.jsonl`. They
record successful/skipped/error files and reasons. Retry repaired inputs with
`--retry-errors PATH`. A batch with file failures returns nonzero after processing
remaining inputs. Global setup/storage failures can still stop the batch.

Temporary outputs are published atomically. Resume checks input identity,
configuration, implementation/dependency fingerprints and output checksums.
An exclusive lock protects each output; remove a stale lock only after confirming
its worker stopped. Keep raw files immutable during processing. The new upstream
policy changes the implementation fingerprint, so older results are regenerated
under `--resume` rather than silently reused.
