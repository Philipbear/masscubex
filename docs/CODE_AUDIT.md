# Extraction performance audit

## Publication checks — 2026-11-05

The complete local suites pass: **637 source tests and 637 installed-wheel
tests**, with no skips or expected failures. Ruff and `git diff --check` pass.
The wheel was imported from its installation directory with the source-path
override disabled. All 35 packaged Python modules match the source files; the
wheel includes the Apache/BSD notices and no generated Python/Numba caches.
The source archive now also includes the independent reference hash fixture.

The first full run found 93 failures. Every failed node also failed against
the frozen pre-optimization source. The fixtures still expected the earlier
weighted-average mass, floating-point area, five-observation minimum,
seconds-based splitting or former MS2 selection/filtering defaults. These tests
were migrated to the current contract: pinned upstream functions supply exact
quantification expectations, toy-data thresholds are explicit, and switching
fixtures check detected and seeded features separately. No tests were skipped
or marked xfail, and no scientific implementation was changed to make them pass.

Two public sample/blank files were processed through both source and wheel;
their Parquet bytes and scientific metadata match the saved reference.
The numerical implementation fingerprint remains
`fca3d76939fc32766261a865252ba33710f8f90f5abbe89ff5790748268d8a8e`.
The [publication record](../benchmarks/results/publication-validation-2026-11-05.json)
retains counts, hashes, the original failed-node audit and smoke results.
Pushing runs the existing Linux/macOS/Windows, Python 3.10–3.13 CI matrix.
These local checks do not themselves establish that remote matrix's outcome.

The first remote run exposed NumPy 2.4's removal of `np.trapz` on newer Python
jobs. The dependency is now constrained to `numpy>=1.24,<2.4`, preserving the
verbatim upstream functions and arithmetic. This packaging correction leaves
the numerical implementation fingerprint and the Mac benchmark environment
unchanged.

## Overnight continuation — 2026-11-05

The continued inspection covers 35 modules, including the new lightweight
format/tolerance imports, measurement specializations and worker initializer. Scientific defaults
are unchanged. Fifty additional MSV000086989 mzXML files were downloaded;
validation now covers 100 unique files (54 biological samples, 42 blanks and
four QCs), totaling 1.32 GB and 10,039,342 output features. This remains one
positive centroid-mzXML cohort, not a representative sample of all repositories.

Every final-code table and scientific sidecar matches the frozen pre-round
baseline at one, two, three and four process workers. Comparisons exclude only
elapsed/stage timings and the implementation fingerprint. The benchmark records
input/table/provenance hashes; generated tables are removed after checking.

The added execution changes are:

- Reuse first-apex state in Gaussian/noise measurements and retain the exact
  one-dimensional NumPy area expression and Pearson arithmetic.
- Apply the raw-median peak gate in native code, with NumPy-equivalent float32
  arithmetic, and preserve the existing SciPy tie fallback.
- Discover bounded sparse candidates once, using compact IDs/offsets, then
  apply each anchor's eligibility snapshot in one call. Reuse partner EIC time
  joins; retain processed float32 spectra without redundant float64 storage.
- Load and validate numerical libraries in workers. Prime the actual production
  array layouts and explicit calibration arguments before file processing.
- Prime provenance and Parquet schema/writer imports before allocating the
  first large feature graph. The discarded one-row buffer is not an output.
- Drain/recycle long-running process pools, excluding resume hits and avoiding
  empty final pools or growing iterator chains.

The fresh-cache focused run passed 262 checks. The final focused run passed
**265 checks**, including the added recycling/resume and extreme-float32
regressions; Ruff and the whitespace check pass. Full release validation remains
deferred. One additional
older reader assertion expects duplicate m/z entries to be merged by the raw
reader. It fails identically in the frozen baseline and current version; all
reader parsing function/class bodies are unchanged. Duplicate handling remains
in upstream preprocessing. This is recorded alongside the earlier test debt
below, rather than changing scientific behavior during optimization.

Measured performance, worker/RAM tradeoffs and sustained-run results are in
[the overnight report](../benchmarks/results/OVERNIGHT_OPTIMIZATION.md).

## Previous five-file round — 2026-11-04

The audit covered all 31 Python modules, following input discovery, validation,
decoding/cache replay, preprocessing, calibration, tracing, splitting,
measurements, MS2 assignment/seeding, grouping and transactional output. It used
the pinned MassCube 1.2.22 reference at commit
`356748ede691862b9a87f67812db5e728a703b9d` and five previously downloaded sample,
QC and blank mzXML files. No additional real-data cohort was processed.

| Area inspected | Finding/action |
| --- | --- |
| CLI and reports | Input discovery and worker submission remain lazy/bounded; file failures continue to be reported independently. No scheduling/policy changes. |
| Readers and array validation | Streaming XML, supported binary encodings, finite-array validation and polarity metadata remain unchanged. |
| Raw scan cache | Budgeted replay and overflow-to-streaming behavior remain unchanged; no unbounded raw-file cache was introduced. |
| Adaptive noise and profile detection | No threshold, bin, centroid or intensity changes. |
| Calibration and mass windows | All selection, fit, fallback and override policies remain unchanged. Effective settings are compared in full-default sidecars. |
| Trace tracker | Last-observed mass, strict nearest/tie rules, stable intensity priority and ten missing scans remain unchanged. |
| Splitting and measurements | Replace temporary parent/child objects with intervals and defer measurements until final leaves. Cache count/apex/trim state. Preserve original score/area arithmetic, SciPy smoothing, raw valleys and both passes. |
| MS2 association | Cache invariant float32 minute bounds/heights instead of reconstructing arrays for each candidate; preserve scalar comparison semantics and selection order. Empty-record runs avoid unnecessary indexing. |
| Seeded gap filling | No measurement, mass-radius, missing-height, deduplication or status changes. |
| Ion grouping | Cache bounded static mass/RT candidates and consume them against the current assignment generation. Preserve strongest-first order, eligibility snapshots, overwrites, isotope rules and exact NumPy correlations. |
| Output/provenance | Schema, rounding, row groups, compression, checksum/resume and atomic publication remain unchanged. New numerical modules enter the implementation fingerprint. |

The smoothing adapter calls SciPy's actual C correlation kernel, so it preserves
the installed library's accumulation behavior. It uses public SciPy when private
helpers are unavailable. The compiled peak selector defers ambiguous distance
priorities to SciPy. Large catalogs, broad query windows, exceptional numeric
ranges and unsupported scalar-promotion paths retain fallbacks. No approximate
Gaussian/correlation scores or relaxed output comparisons were introduced.

Sparse indices and target tables are bounded; final features and relationships
still require memory proportional to a file's feature count. This is streaming
input, not a claim of constant memory for a globally grouped feature table.

## Historical validation and test debt

223 focused checks pass, including full-precision reference measurements,
bit-identical smoothing, plateau/tie cases, sparse-query boundaries, grouping,
calibration, output formats, failed-file continuation and atomic output failure.
The five full-default Parquet tables and scientific sidecars are compared against
the frozen pre-optimization version. A separate benchmark compares original
MassCube measurements/MS2 references/group labels on identical workloads.

A wider check of three older MS2/seeding/polarity fixture files found **52 failed
and 32 passed tests in both the frozen baseline and optimized version**, with
exactly the same failed test names. Several fixtures still assume former rules
(apex tie breaking, dropping short traces, rejecting multiple precursors or
floating rather than upstream integer areas). Other polarity/seed expectations
need a separate scientific-policy review. These failures were not introduced by
this speed change, were not hidden by changing production algorithms, and must
be resolved before release validation. Their names are retained in the benchmark
JSON audit record. This audit is not a clean full-release-suite result.

Those earlier rounds stayed local and deferred publication checks. The
publication section above records the subsequent full-suite reconciliation.
