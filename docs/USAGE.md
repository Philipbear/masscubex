# MassCubeX usage and output reference

For installation and a short introduction, see the [README](../README.md).
For algorithm details, see [METHODS.md](METHODS.md).

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
# One file, or all supported files in a directory tree
masscubex extract sample.mzML -o features
masscubex extract /data/MSV000086989 -o features --workers 8 --resume

# Millions of paths without a giant shell argument list
masscubex extract --manifest files.txt -o features --workers 8 --resume > jobs.jsonl

# Threads are also supported; benchmark against processes on your storage
masscubex extract /data/study -o features --workers 8 --executor thread --resume

# Record area as the preferred measurement; both height and area are written
masscubex extract sample.mzML -o features --quant-method peak_area

# Compressed Parquet is the default; TSV remains available
masscubex extract sample.mzML -o features --format tsv

# Adduct / in-source-fragment candidates run by default; opt out when needed
masscubex extract sample.mzML -o ungrouped --no-group-ions
```

Manifest paths are UTF-8, one per line, relative to the manifest's directory
unless absolute. Empty lines and lines starting with `#` are ignored. Directory
input is recursive. Symlinked directories are not followed. Avoid duplicate
inputs in the same job.

Each input produces a Parquet table and a `.parquet.json` sidecar by default.
CLI outputs are sharded as
`features/<hash-prefix>/<input-name>.<path-hash>.features.parquet`. With
`--format tsv`, the endings become `.tsv` and `.tsv.json`. The hash uses the
resolved absolute input path, preventing collisions between repeated basenames
or mzML/mzXML exports. Moving inputs changes their output IDs. A JSON-lines job
report on stdout maps every input to its output; stderr carries the summary.
Corrupted or missing files do not stop other inputs. Every run automatically
saves a summary and an error-file list, as described below. A batch with file
failures returns a nonzero exit code **after processing the remaining files**.

## Memory and scan replay

The default decoded-scan cache budget is **64 MiB per worker**. Calibration,
extraction and raw-MS1 gap filling reuse it within a file. If a file exceeds the
budget, extraction discards the cache and streams subsequent passes from disk.
`--scan-cache-mb 0` disables it; a larger value can help large files when RAM
permits. Feature histories, MS2 metadata and Arrow buffers need additional memory.
The cache does not reuse features across files or runs, and does not replace
checksum-verified `--resume` behavior. Sidecars report read/replay counts and
peak charged cache bytes under `scan_cache`.

## Corrupted files and batch reports

No extra flag is needed to continue after malformed XML, broken binary arrays,
truncated gzip data, missing files or invalid file metadata. The failed file is
reported, its partial result is discarded, and other files continue. The CLI
automatically creates a unique report directory for each invocation:

```text
features/reports/<run-id>/summary.json   # counts, run status, settings and timing
features/reports/<run-id>/errors.tsv     # failed paths, exception types and messages
features/reports/<run-id>/files.jsonl    # every success, verified skip and error
```

The report location is printed to stderr. `--report-dir /path/to/reports`
changes the parent folder. Logs are streamed during processing; the summary
is periodically checkpointed and finalized at the end. Previous reports are
preserved, including when resuming or retrying. Reports do not require holding
millions of filenames or errors in memory.

After fixing or redownloading failed inputs, retry only that report's sources:

```bash
masscubex extract --retry-errors features/reports/RUN_ID/errors.tsv \
  -o features --config params.json --workers 8 --resume
```

Replace `RUN_ID` with the run folder printed by the CLI. Use the same processing
settings as the original run; `--retry-errors` reads
paths, not settings. Include `--format tsv` if that was the original format.
A new report records the retry. Omit `--config params.json`
if the original run used defaults. Use a TSV parser to read error messages or
paths containing tabs/newlines.

`summary.json` distinguishes `completed`, `completed_with_errors`, `interrupted`
and `aborted`. Exit code 1 can mean a completed batch with file errors; it does
not imply that later files were skipped. Directory/manifest discovery failures
are recorded too, with `stage=discovery`; their unseen contents are not counted
as an exhaustive failed-file list. Unexpected worker losses use `stage=worker`
and a new process pool handles undispatched files. An affected file is not
automatically known to be corrupt. Global setup/report-storage failures and
user interruptions can still abort the run; see the
[failure and restart details](METHODS.md#10-millions-of-files-scheduling-failures-and-reports).

## Python API

For an exact output path or your own scheduler:

```python
from masscubex import Config, process_file

result = process_file(
    "sample.mzXML",
    "features/sample.parquet",
    config=Config(noise_layout="bins", noise_window_da=100, noise_sigma=3, min_scans=5),
    resume=True,
)
print(result.features, result.elapsed_seconds)
```

Without an output path, `process_file("sample.mzXML")` writes
`sample.mzXML.features.parquet`. An explicit `.tsv` path selects TSV; alternatively
use `output_format="tsv"`. A path/format mismatch is rejected. The CLI `--format`
option is an execution setting, separate from the scientific JSON configuration.

Read selected Parquet columns without parsing the rest:

```python
import pyarrow.parquet as pq

table = pq.read_table("features/sample.parquet", columns=["mz", "peak_height"])
```

For typed row-by-row access to either format, including normalized relationship
dictionaries, use `from masscubex import iter_feature_rows` and iterate over
`iter_feature_rows("features/sample.parquet")`.

Settings are immutable and file state is independent, so the API can be called
from separate threads or processes. The default CLI uses one worker; multiple
workers use processes unless `--executor thread` is selected.

## Parameters and ion selection

```bash
masscubex extract --help
masscubex defaults > params.json
masscubex extract --config params.json --print-config
masscubex extract /data/study -o features --config params.json --workers 8 --resume

masscubex adducts                       # numbered catalog and enabled state
masscubex adducts --json                # exact ionic masses and definitions
masscubex extract sample.mzML -o grouped --group-ions --adducts 1-15,22,23
```

Every processing parameter is available in JSON; explicit CLI flags override
that file. The resolved configuration is saved beside each table. Custom ion
definitions can be supplied through `--adduct-file ions.json`; see the
[ion-grouping reference](ION_GROUPING.md).

The selected defaults are **1–15**: H, Na, K, NH4, deprotonation, Cl, Br,
formate, acetate, acetonitrile, methanol, and the four listed dimers.
**16–21 are disabled** (multiple charges and trimers). **22 and 23** add
`[M+H-H2O]+` and `[M+H-2H2O]+` as neutral-loss candidates, separately from
ordinary adducts. Grouping is on by default; `--no-group-ions` disables it.

Grouping requires compatible polarity, mass, apex time and correlated
chromatographic traces: defaults are 10 ppm / 0.002 Da, a calibrated apex-time
tolerance (3 seconds when calibration is unsupported), at least
five shared scans, Pearson correlation ≥ 0.9, and ≥ 50% overlap relative to
the longer feature. Water losses can be proposed from this evidence alone;
`ms2_supported` distinguishes those also observed in an associated parent
MS/MS spectrum. Generic in-source-fragment candidates require MS/MS support.
These are competing **candidate relationships**, not confirmed identities or
forced compound groups. `--no-isf` disables both neutral-loss and fragment
analysis while keeping ordinary adduct hypotheses.

## MS/MS scans and representative selection

MS2 association is enabled by default. Each accepted single-precursor MS2 scan
is assigned to at most one feature using precursor m/z, polarity and the
feature's time interval. Calibration estimates these settings independently
per polarity; fallback values are 20 ppm / 0.01 Da and a 2-second boundary
extension. Tolerance uses the larger of the ppm and Da limits.
Overlapping candidates are ranked by containment in the actual interval,
mass error, apex-time distance, peak height, then feature ID.

**The best scan defaults to the largest sum of exported MS2 peak intensities**,
preserving upstream MassCube's ranking criterion. It sums the actual arrays,
including any residual precursor peak, before optional fragment-evidence
pruning. Vendor TIC metadata can describe different pre-export samples, so it
does not replace this sum. Ties prefer proximity to the feature apex, then
earlier acquisition. MS1 adaptive noise filtering is never applied to MS2.

```bash
# Explicit default; alternatives remain available
masscubex extract sample.mzML -o features --best-ms2 fragment_intensity
masscubex extract sample.mzML -o features --best-ms2 precursor_intensity
masscubex extract sample.mzML -o features --best-ms2 apex

# MS1-only throughput comparison
masscubex extract sample.mzML -o ms1_only --no-ms2
```

Precursor ranking uses recorded precursor intensity, then a matching preceding
MS1 measurement if metadata is absent or zero. The fallback checks polarity,
parent reference when supplied, m/z and scan age. It uses the nearest centroid
or the profile apex within the tolerance. If all associated scans lack that
intensity, selection falls back to apex proximity. The sidecar records
`config.best_ms2` and the exceptional `best_ms2_apex_fallback_feature_ids`.
The table retains best-scan identifiers, fragment/precursor intensities and the
precursor-intensity source.

The isolation-width cutoff is **disabled by default** (`ms2_max_isolation_width=0`).
A positive `--ms2-max-isolation-width` opts into excluding wider MS2 scans; this
is a full precursor isolation width, not a matching tolerance or DIA classifier.
Old configuration files containing `4` retain that explicit value; override
with `--ms2-max-isolation-width 0` or update the JSON.

By default, an unmatched single-precursor MS2 with usable m/z and time creates
an **MS2-seeded feature**. A narrow raw EIC anchored at its precursor m/z
supplies height/area. With supported calibration, local peak boundaries
define integration. Repeated scans merge only after recovering the same MS1
apex. Insufficient evidence uses a fixed ±5-second window and separate
±max(0.002 Da, 10 ppm) fallback mass radius. Both modes cap the radius at the
normal MS1 tolerance; a wider gap-fill override cannot broaden it. They bypass the MS1 detection threshold
and minimum scan count. Existing detected-feature assignments are preserved.

Use `--ms2-gapfill-rt-window SECONDS` to change that half-window, or
`--no-ms2-seeded-features` to keep only independently detected MS1 features.
`feature_origin` and `ms1_quant_status` identify the added rows. If no MS1 ion
is found, height is imputed as 0.1× regional noise, with status `imputed_noise`,
missing area and zero observed scans. `apex_noise_threshold` records the noise
reference; it is the median positive threshold from up to the three closest
compatible MS1 scans in the query window. Use `--gapfill-noise-fraction 0` to
leave height missing instead. Missing scans, unsupported noise, unresolved
polarity, local-peak mismatch and exceeded limits remain blank. Unusable precursor metadata,
multiprecursor spectra, and an explicitly enabled width cutoff are recorded
per scan in the JSON sidecar's `unassigned_ms2`, with native IDs, source indices
and reasons. The sidecar counts partition all parsed MS2 into assigned and
unassigned scans. `--no-ms2` explicitly turns this bookkeeping off.

For reverse lookup, use `(input file, ms2_spectrum_index)` as the unique key;
scan numbers can repeat or be absent. Search `ms2_spectrum_indices` in the
feature table, then `unassigned_ms2` if no row matches. A valid usable
single-precursor MS2 maps to exactly one detected or seeded feature under the
default settings. This is a deterministic association, not proof of a unique
biological precursor. DIA classification and fragment deconvolution are not implemented. Profile/unknown-mode MS2 can
still be associated and scored; profile scores sum nonnegative exported samples.
Only explicitly centroided MS2 supplies fragment-mass evidence for grouping.

## Automatic settings and polarity switching

Calibration and ion grouping run by default. The pilot starts at
±max(0.05 Da, 100 ppm), then requires at least 20 qualified strong peaks per
polarity; precursor mass matching needs evidence from at least 20 distinct
MS1 representatives with MS2. Unsupported estimates retain configured fallback
values and report the reason. Noise remains fixed at μ + 3σ.

Use `--no-auto-tune` for static settings, or explicit flags to protect individual
settings, for example `--mz-ppm 10 --ms2-mz-ppm 20`. Non-default JSON/Python
values are protected too; include default-valued names in `tuning_locked` to
protect them. `masscubex defaults` prints fallback settings;
`parameter_tuning.effective_configs` in each result's JSON shows actual values.

Positive and negative MS1 streams use their own cadence, learned tolerances,
traces and gap filling. Missing MS2 polarity can be resolved by an explicit
parent reference; unresolved mixed-polarity cases keep MS2 references with blank
quantification. See [all algorithms and controls](AUTO_TUNING.md).

## Feature-table schema

Detected rows are emitted in trace-completion order, followed by MS2-seeded
rows in seed-creation order. IDs are local to each file. Retention times are
always **seconds**. Parquet stores float64 values, nullable int64 scan IDs, native
lists and relationship structs, using Zstd level 3 compression. Row groups contain
at most 8,192 features. Output schema 7 has **29 columns**. Both formats round
m/z, seed m/z and mass radii to **5 decimal places**, all retention times in
seconds to **2**, and every intensity/area/threshold to **2**. TSV writes trailing
zeros; Parquet stores the same rounded numbers without text formatting. Nested
ion evidence uses 5 decimals for masses, 2 for ppm errors and relative
intensities, and 6 for correlation. The sidecar and Parquet metadata specify
this policy in `output_precision`. Non-finite feature output fails the file.

Computations stay float64; rounding does not affect extraction, grouping,
scan associations or representative selection. Keep scan IDs when retrieving
raw data: rounded coordinates may fall on a different side of a matching or
integration boundary. Do not treat the displayed mass radius as the exact
processing tolerance; full settings remain in the sidecar.
Missing scalars are null in Parquet and blank cells in TSV.

| Column | Meaning |
| --- | --- |
| `feature_id` | Sequential identifier within this file |
| `mz` | Intensity-weighted mean for detected features; original precursor m/z for seeded features |
| `rt_apex_seconds` | Time of its highest original intensity |
| `rt_start_seconds`, `rt_end_seconds` | Integration boundaries |
| `peak_height` | Detected/raw MS1 height, or a labelled noise-derived estimate when status is `imputed_noise`; blank if unavailable |
| `peak_area` | Trapezoidal integral, intensity × seconds; seeded rows integrate refined local bounds or the fallback window |
| `scan_count` | Number of observed nonzero MS1 samples |
| `apex_noise_threshold` | Local MS1 threshold at the detected apex; regional reference for an imputed seed |
| `gaussian_similarity` | Symmetric apex/FWHM Gaussian Pearson score on the original detected trace, 0–1, two decimals; null for MS2-seeded forced measurements |
| `segmentation_status` | `split`, `unsplit`, `guarded_split`, `guarded_unsplit`, `disabled`, or `not_applicable` for seeds; `guarded_` means the noisy trace received extra shape checks |
| `polarity` | `positive`, `negative`, or `unknown` |
| `feature_origin` | `ms1_detected` or `ms2_seeded` |
| `ms1_quant_status` | `detected`, `gap_filled`, `imputed_noise`, `no_ms1_signal`, `no_ms1_scans`, `ambiguous_polarity`, `no_local_ms1_signal`, `no_local_peak`, or `sampling_limit` |
| `seed_precursor_mz`, `seed_rt_seconds` | Original first MS2 anchor for seeded rows; blank for detected rows |
| `gapfill_status` | `fixed_window`, `peak_bounds`, `single_observation`, or a refinement failure status; blank for detected rows |
| `gapfill_mz_radius_da` | Actual narrow EIC radius when measured; blank if unavailable |
| `ms2_count` | Number of associated MS2 scans |
| `ms2_scan_numbers` | List of native scan numbers; `null` where the format has no numeric scan number |
| `ms2_native_ids` | List of original spectrum identifiers |
| `ms2_spectrum_indices` | List of zero-based positions among all source spectra |
| `best_ms2_scan_number`, `best_ms2_native_id`, `best_ms2_spectrum_index` | Identifiers of the selected scan; empty if unavailable |
| `best_ms2_fragment_intensity` | Sum of exported MS2 intensities, when decoded |
| `best_ms2_precursor_intensity`, `best_ms2_precursor_intensity_source` | Precursor intensity and `metadata`, `preceding_ms1`, or `unavailable` |
| `ion_relations` | List of optional adduct, neutral-loss and in-source-fragment candidates |

`intensity`, `quant_method`, `best_ms2_selection_method`, `best_ms2_score` and
`apex_threshold_ratio` are omitted from schema 6. Height and area are both
available; `--quant-method` records the preferred measurement once in
`config.quant_method` (default `peak_height`) without changing the table's
measurements. Best-spectrum ranking is in `config.best_ms2`. A height/threshold
ratio can be calculated when the threshold is nonzero; it is not a statistical S/N.

`ambiguous_mass` is a historical v0.9 status, no longer emitted. Its broad
candidate search compared the two highest local MS1 heights; a ratio below
3 caused missing quantification. Overlapping candidate windows could count
the same centroid twice. Version 0.10 uses one precursor-anchored EIC within
the normal MS1 tolerance and removes that rule. See the [audit](../benchmarks/results/EIC_REVIEW.md).

Schema 7 adds the two shape diagnostics. By default, traces with MassCube's
noise score above 3 need supported peak shapes before splitting. Otherwise
the whole trace is retained, with its original measurements and MS2 references.
A low Gaussian similarity or a guarded status is not a contamination diagnosis.
Use `--no-peak-shape-guard` for previous splitting behavior; thresholds are
`--peak-noise-score-limit 3`, `--peak-min-snr 5`, and
`--peak-min-gaussian-similarity 0.7`. The S/N check is across retention time;
MS1 mass-axis noise remains μ + 3σ. See the [method](METHODS.md#shape-guard-introduced-in-011)
and [real-file comparison](../benchmarks/results/PEAK_SHAPES.md).

MS2 arrays are ordered by acquisition and correspond position by position.
A native number is never invented from a spectrum index. Parquet lists are
already typed; relationship structs have null fields for evidence absent from
that relationship kind. Read TSV with a TSV/CSV parser, then JSON-decode its
list fields. `[]` means no associated scans
or no proposed relationships. A best native ID may exist even when its native
scan number is unavailable.

A threshold ratio is not a statistical S/N estimate. Missing scans within a
trace contribute zero intensity at their actual acquisition times. Integration
trims trailing missing scans. Segments share a valley point but never count a
trapezoid interval twice. `--no-segment` retains entire ion traces. Five nonzero
scans are required per independently detected MS1 feature by default;
MS2-seeded rows bypass this detection requirement. Their apex is the first
maximum raw MS1 measurement, falling back to the seed time if unavailable.
Their noise-threshold fields are blank. `gapfill_status` and the preserved
seed coordinates distinguish refinement from fixed-window fallback.

The JSON sidecar records requested configuration, per-polarity effective
configuration, calibration evidence/fallbacks, version, implementation/dependency
fingerprint, input identity, scan/point counts, MS2 assignments, unresolved
IDs/reasons, polarity inferences, seeded-feature and relationship counts,
profile-scan count, memory-pool diagnostics, elapsed time and output SHA-256.
`output_format`, `output_bytes` and `output_sha256` identify either format;
TSV also retains the compatibility keys `tsv_bytes` and `tsv_sha256`.
`timings_seconds` records calibration, feature extraction (including reading),
MS2 association, raw-MS1 gap filling, ion grouping and table writing. Metadata
publication, checksums and other overhead are outside these stage timers.
A valid MS1 run with no detected or seeded features writes an empty table with
the complete schema (or a header-only TSV).
MS2-only files with usable precursors produce seeded rows with blank MS1
quantification. Inputs with neither MS1 nor usable MS2 seeds are errors.
These tables provide m/z, RT, polarity and quantitative
values for a future independent dataset-level alignment process; alignment and
retention-time correction are not implemented here.

## Input integrity and restarting

- mzML supports indexed wrappers, referenceable parameter groups, little-endian
  float32/float64, uncompressed/zlib arrays, and optional MS-Numpress (including
  zlib combinations). mzXML supports nested scans, 32/64-bit pairs, big/little
  endian and zlib. Other encodings fail explicitly.
- Metadata must identify profile or centroid mode. `--assume-centroid` permits
  **missing** mode metadata only; it never overrides explicit profile metadata.
- MS1 times must increase strictly within each polarity. Ion-mobility-resolved
  frames must be collapsed first. MS2 references and scores are exported;
  complete fragment spectra remain in the source files.
- Arrays, time units and counts are validated. Unsorted masses are sorted and
  duplicate masses summed. Invalid arrays abort that file rather than creating
  a successful but incomplete table.
- Outputs use temporary files and atomic replacements. The JSON sidecar is
  published last as the completion marker. `--resume` checks version, implementation/dependencies, settings,
  input path/size/mtime, and the output checksum; stale/corrupt results are
  recomputed. Inputs are not content-hashed, so keep input files immutable.
- Existing outputs require `--resume` or `--overwrite`. Exclusive `.parquet.lock` / `.tsv.lock`
  files prevent simultaneous writes. After a killed worker, inspect the lock's
  PID/host and remove it only once that worker has stopped. Interrupted
  `.partial` files can then be removed. Rename/lock semantics require a
  filesystem with atomic local-style operations.

