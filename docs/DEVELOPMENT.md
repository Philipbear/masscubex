# Code structure and numerical development

`pipeline.process_file` owns the file transaction: validate/resume, calibrate,
extract, associate/gap-fill MS2, group ions, write a temporary table, then publish
it with a checksum and completion marker. `cli.run_batch` owns scheduling and
per-file failure reporting. Neither contains per-mass numerical loops.

| Module | Responsibility |
| --- | --- |
| `readers.py`, `array_validation.py` | XML/binary decoding and fast normalized-array validation; unusual arrays use the full normalization/error path. |
| `scan_cache.py` | Budgeted decoded-scan replay with streaming fallback and independent MS2 metadata per pass. |
| `features.py` | Recyclable trace pools, batched native feature extraction and the `Feature` record. Calibration uses the raw-trace interface. `MEAN_MZ`, `HEAD_NODE`, etc. name pool columns; `available_nodes` and `active_count` name state counters. |
| `chromatograms.py` | Compiled segmentation and original-sample measurements. Named `PeakSettings`, `PeakSummary` and `PeakCandidate` records carry options, measurements and fallback status. |
| `peak_shapes.py` | MassCube-style Gaussian/roughness scores and additional shape evidence for splitting noisy traces. Missing-scan interpolation is used only for segmentation. |
| `reference_chromatograms.py` | Readable SciPy compatibility path for tied peak priorities and a differential-test oracle. |
| `calibration.py` | Compiled quality screening of one representative trace; accepted candidates retain NumPy/SciPy measurement arithmetic to preserve calibrated boundaries. Named `RepresentativeValues` and explicit rejection reasons. |
| `tuning.py` | Gaussian-qualified representatives, regional per-polarity estimates, range provenance, evidence minima and override policy. |
| `mass_tolerance.py` | Immutable region records and consistent Python/compiled target-mass window lookup. |
| `noise.py`, `noise_regions.py`, `centroid.py` | Sliding/binned MS1 noise estimation and profile mass detection. |
| `ms2.py`, `gapfill.py`, `gapfill_refinement.py` | MS2 bookkeeping, unmatched precursor seeds, query scheduling and raw-MS1 measurements. |
| `grouping.py` | Mass-window candidates and evidence-qualified ion relationships. |
| `output.py` | Separate TSV and columnar Parquet writers, shared output schema and typed readers. |

Python handles policies and labelled records; Numba handles repeated numerical
work. All calculations remain float64 and kernels avoid fast-math. Round only
at the output boundary according to `output.DECIMAL_PLACES`; table values are
not full-precision interval coordinates. Keep exact pre-export coordinates for
independent raw-data audits. Helpers
called only from a compiled kernel inherit that caller's GIL-free execution;
public heavy kernels use `nogil=True`. Avoid nested numerical thread pools
when distributing files across workers.

Preserve acquisition order, inclusive mass/time boundaries, first-maximum and
first-valley tie rules, polarity separation and scan provenance when changing
kernels. Summation and transcendental implementations can differ in the final
float64 bits. Differential tests compare numerical measurements tightly and
require identical boundaries, counts and identifiers. Do not replace those
checks with a feature-count-only comparison.

Run `pytest -q` and `ruff check src tests benchmarks`. Tests include independent
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
