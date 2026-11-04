# Code structure and numerical development

`pipeline.process_file` owns the file transaction: validate/resume, calibrate,
extract, associate/gap-fill MS2, group ions, write a temporary table, then publish
it with a checksum and completion marker. `cli.run_batch` owns scheduling and
per-file failure reporting. Neither contains per-mass numerical loops.

| Module | Responsibility |
| --- | --- |
| `readers.py`, `array_validation.py` | XML/binary decoding and fast normalized-array validation; unusual arrays use the full normalization/error path. |
| `formats.py` | Lightweight filename recognition, shared by discovery and numerical readers. |
| `scan_cache.py` | Budgeted decoded-scan replay with streaming fallback and independent MS2 metadata per pass. |
| `upstream.py`, `upstream_adapter.py` | Unchanged pinned scientific function bodies, preprocessing, reference/fallback adapters. |
| `summary.py` | Equivalent interval splitting, first-apex/count kernels, exact SciPy smoothing and original NumPy measurement arithmetic. |
| `measurements.py` | Exact specializations of upstream Gaussian/noise/area arithmetic and NumPy's two-vector Pearson calculation. |
| `features.py` | Recyclable trace pools, batched native feature extraction and the `Feature` record. Calibration uses the raw-trace interface. `LAST_MZ`, `HEAD_NODE`, etc. name pool columns; `available_nodes` and `active_count` name state counters. |
| `chromatograms.py` | SciPy-equivalent calibration primitives and peak selection; ambiguous distance ties use the original SciPy path. |
| `peak_shapes.py` | Compiled Gaussian diagnostics for calibration and numerical checks. |
| `reference_chromatograms.py` | Independent Gaussian diagnostic for calibration checks. |
| `calibration.py` | Compiled quality screening of one representative trace; accepted candidates retain NumPy/SciPy measurement arithmetic to preserve calibrated boundaries. Named `RepresentativeValues` and explicit rejection reasons. |
| `tuning.py` | Gaussian-qualified representatives, regional per-polarity estimates, range provenance, evidence minima and override policy. |
| `pilot_sampling.py` | Deterministic strongest-centroid selection per mass region and complete acquired-range discovery for optional equal-width layouts. |
| `mass_tolerance.py`, `mass_tolerance_kernels.py` | Lightweight immutable region definitions, Python lookup and lazily loaded compiled lookup. An unbounded upper limit is `None` in records/JSON and infinity in numerical arrays. |
| `noise.py`, `noise_regions.py`, `centroid.py` | Sliding/binned MS1 noise estimation and profile mass detection. |
| `ms2.py`, `gapfill.py`, `gapfill_refinement.py` | MS2 bookkeeping, unmatched precursor seeds, query scheduling and raw-MS1 measurements. |
| `grouping.py`, `grouping_spectra.py` | Upstream greedy assignments and temporary file-backed partner EIC queries. |
| `grouping_queries.py` | Bounded sparse mass/RT candidate discovery; assignment state remains in the greedy consumer. |
| `output.py` | Separate TSV and columnar Parquet writers, shared output schema and typed readers. |
| `worker_runtime.py` | Worker initialization: compile the production layouts and prime provenance/Parquet schema allocation before a large feature graph exists. |

Python handles policies and labelled records; Numba handles repeated numerical
work. Scientific scan arrays and measurement arithmetic follow upstream float32 behavior,
including integer height/area finalization. Kernels avoid fast-math. Format only
at the output boundary according to `output.DECIMAL_PLACES`; table values are
not full-precision interval coordinates. Keep exact pre-export coordinates for
independent raw-data audits. Helpers
called only from a compiled kernel inherit that caller's GIL-free execution;
public heavy kernels use `nogil=True`. Avoid nested numerical thread pools
when distributing files across workers.

Preserve acquisition order, strict upstream mass windows and inclusive MS2 time bounds, first-maximum and
first-valley tie rules, polarity separation and scan provenance when changing
kernels. Summation and transcendental implementations can differ in the final
float64 bits. Differential tests compare numerical measurements tightly and
require identical boundaries, counts and identifiers. Do not replace those
checks with a feature-count-only comparison.

During development, keep changes local and run focused checks appropriate to the
change. Create PRs or publish changes only when requested; run the full release
checks when the user approves publication. For release validation, run
`pytest -q` and `ruff check src tests benchmarks`. Tests include independent
SciPy comparisons, exhaustive small plateau/distance cases, corrupted inputs,
profile scans, polarity switching, cache overflow/disable paths, atomic output,
and Parquet/TSV round trips. `NUMBA_DISABLE_JIT=1` is useful for debugging selected
small numerical fixtures. A full source and installed-wheel matrix runs on
Linux, macOS and Windows with Python 3.10–3.13.

For speed work, benchmark serially on immutable inputs with identical scientific
settings. Exclude JIT warmup, retain full-pass timings and implementation hashes,
and profile separately. `profile_extraction.py` reports inclusive substeps;
shared scan replay and tracker time overlap several stages. Never add them to
parent stages or present a kernel speedup as whole-workflow throughput.

After changing a helper imported by another compiled module, clear generated
Numba `.nbc`/`.nbi` caches under `src/masscubex/__pycache__` before validation.
Numba's disk cache does not reliably invalidate callers after a change in a
different source file. Never ship generated caches in the wheel. The clean
installed-wheel and CI runs verify compilation without local cache assumptions.

`compare_noise_versions.py` checks every raw MS1 mean, standard deviation and
threshold against an archived noise implementation, including discarded points.
`compare_extraction_versions.py` compares complete tables and scientific sidecars,
not just feature counts. Time these audits separately from throughput runs.

Pinned upstream function bodies in `tests/upstream_masscube.py` provide the
independent scientific reference. `tests/upstream_function_hashes.json` verifies
that runtime copied bodies remain unchanged. `verify_upstream_defaults.py`
explicitly disables calibration, profile centroiding and MS2 seeding to compare
full measurements/MS2 references on real files. These extensions are enabled
in normal default runs. Preserve source attribution.

Performance adapters reuse finalized measurements and stop retrying an unchanged
split result. `ExtractionFeature` avoids deep-copying full parent buffers only
inside the pinned splitter's immediate subset/finalize sequence; it is not a
general-purpose shallow copy of a mutable upstream feature. Tests compare every
exported measurement at full precision and verify parent/child buffer independence.
The sorted grouping lookup index embeds float32 masses exactly in float64 to
avoid NumPy recasting the entire index for every query. The compiled mass query
uses NumPy's resolved subtraction/comparison dtypes, including explicit float32
rounding where required; unsupported/nonfinite cases retain the NumPy path.
Assignment generations retain the same per-anchor eligibility snapshot and
replacement rules. The query and eligibility loops use no fastmath; their strict
mass/RT decisions and priority order match the NumPy/Python reference.

Trimmed trace spans under five samples use a direct equivalent of upstream's
length guards and final measurements. Interior missing-scan zeros count toward
span length. Both diagnostic scores are zero and no split is possible; a
singleton additionally has zero area. Float32 conversion, first-apex selection,
NumPy area arithmetic and array-based minute/second conversion are preserved.
Longer spans use equivalent interval execution and measurement specializations.
Parquet batching resolves column
names once and uses bounded iterator slices; schema, rounding and row groups
are unchanged.

`benchmark_exact_speed.py` processes the explicitly specified files with unchanged
defaults, a warmup and measured repeats. Select an archived source tree with
`--package-path`; supply its output directory with `--compare` for the optimized
run. It requires byte-identical Parquet output and equal scientific sidecars,
excluding only elapsed/stage timings and the implementation fingerprint.
`--profile-first` profiles only the separate warmup, never a measured repeat.

`benchmark_resources.py` measures process/thread throughput and explicitly owned
worker RSS under a configured budget. It checks table and scientific-sidecar
hashes against a frozen reference and deletes generated outputs after checking.
`benchmark_cli_resources.py` exercises the actual scheduler over repeated cycles,
including startup, checksum verification and pool recycling. It observes only
the executor's known process IDs. Repeated cycles are not additional independent
raw files. RSS sums include shared pages in each process and are conservative
estimates rather than unique private-memory measurements.

The process coordinator keeps numerical imports in workers. Initialization is
validated before file submission. CLI process pools drain and restart after
`8 * workers` extraction attempts by default; resume hits are excluded. One
deferred path avoids unnecessary final-pool startup without accumulating nested
iterators over long runs. This policy applies to process pools with more than
one worker; `--worker-recycle-files 0` disables it.

`audit_warmup.py` compares loaded kernel signatures before/after a real-file
run. Warmup covers production float32 column views, reversed finished IDs,
calibration's explicit arguments and read-only compact grouping maps, preventing
avoidable compilation in the first normal centroid file.
