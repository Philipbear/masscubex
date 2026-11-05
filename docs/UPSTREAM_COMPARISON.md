# Upstream algorithms and MassCubeX defaults

MassCubeX follows [huaxuyu/masscube](https://github.com/huaxuyu/masscube), pinned to
`356748ede691862b9a87f67812db5e728a703b9d` (MassCube 1.2.22).
Core feature algorithms are retained from this exact version. **Adaptive MS1
noise, automatic tolerance tuning, profile MS1 centroiding and MS2-seeded gap
filling are enabled by default**. The selected ion catalog, repository I/O,
polarity separation and compact output are also intentional adaptations.

## Inherited behavior and configured differences

| Step | Current default |
| --- | --- |
| Scan representation | float32 m/z and intensities, as upstream |
| Acquisition acceptance | MassCubeX screens metadata and repeated windows, skipping recognized DIA before extraction; unknown modes remain labelled unknown |
| Scan cleanup | Strict mass/intensity bounds; adjacent-centroid merging at 0.002 Da after filtering |
| Trace matching | Last observed mass, not a running weighted mean; exact midpoint ties do not match. Regional auto-tuned radius by default; upstream 0.01 Da is the uncalibrated fallback |
| Conflicting matches | Test the nearest mass only; if already claimed, insert a missing observation rather than choosing another mass |
| Trace priority | Initial scan in mass order; subsequently stable descending latest-scan intensity, including zero for missed observations |
| Missing scans | Ten allowed; finish on the eleventh miss; no default elapsed-time cutoff |
| Output filtering | Keep short traces, including one observation; no five-observation minimum |
| Segmentation | Original noise >3 veto, sigma 0.6 scans, 3% prominence, distance 3, raw-median apex gate, raw valleys, two passes |
| Feature mass | Apex m/z |
| Quantification | Original integer peak height and trapezoidal area (intensity × seconds), including upstream float32 arithmetic and minute-to-second conversion |
| Gaussian similarity | Original upstream calculation; diagnostic, not an extra feature filter |
| MS2 preprocessing | Original fixed absolute/relative filter and adjacent-centroid merging |
| MS2 association | Adjacent-MS1-scan time edges, strongest feature wins overlaps; calibration supplies mass radii and an RT margin. Without calibration, upstream strict 0.1 Da and zero RT margin |
| Best MS2 | Largest sum of processed fragment intensities; first acquisition wins a tie |
| Unmatched MS2 | MassCubeX creates labelled seeded features for usable unmatched precursors and gap-fills raw MS1; unresolved scans remain in JSON |
| Grouping | Upstream strongest-first isotope/adduct/best-MS2 fragment algorithm, with the selected catalog |

The scientific function bodies for splitting, finalization, scores and centroid
merging are retained unmodified in `src/masscubex/upstream.py`, with source
attribution. The optimized tracker is compared against original `detect_features`.
The interval adapter and compiled peak selector preserve those rules and are
checked against the original functions, including raw valleys and tie behavior.
Feature mass remains the apex mass. Thread/process workers do not mutate shared
upstream globals.

## Enabled extensions

`--auto-tune`, `--profile-mode centroid` and `--ms2-seeded-features` are all on
by default. Calibration changes mass and MS2 timing parameters, without changing
the ten-missing-MS1-scan rule or introducing an elapsed trace-gap cutoff. The
other stages add profile mass detection and unmatched-precursor gap filling.
The actual per-polarity and regional parameters, fallback reasons and enabled
extensions are recorded in JSON provenance. Disable these stages with
`--no-auto-tune --profile-mode reject --no-ms2-seeded-features` for the core
upstream comparison configuration; adaptive noise and the selected catalog remain.
See [methods](METHODS.md) and [calibration](AUTO_TUNING.md).

## Parameters, formats and labels

The default MS2 intensity floor follows upstream header selection: 10,000 when
`orbitrap` occurs in the first 300 XML header lines, otherwise 500. Explicit
`--ms2-min-intensity` overrides this. Adaptive MS1 filtering never uses that
instrument-specific floor. Fixed-floor comparison runs explicitly set the same
MS1/MS2 thresholds in both programs.

The selected catalog remains IDs 1–15 and 22–23. Each strongest ungrouped
positive anchor starts as `[M+H]+`; the negative equivalent is `[M-H]-`.
H/Na monomer pairs do not produce alternative dimer-pair hypotheses. Actual
dimer-mass features remain eligible. Grouping assignments are tentative,
not library identifications.

mzML/mzXML streaming, per-polarity processing, compressed Parquet/TSV output,
scan provenance, failure reports and resumable transactions are repository
adaptations. Outputs have 15 default columns, with peak area, detailed ion
relationships and noise/shape diagnostics available as explicit output options,
and use the requested rounding. The internal
MS1 scan representation is now float32, following upstream; earlier claims that
all computation remained float64 no longer describe this version.

## Verification

Synthetic tests compare tracing, segmentation, measurements, MS2 ownership and
ion grouping with the pinned upstream functions. Function-body hashes guard
the copied scientific code. Other tests cover profile scans, polarity switching,
corrupt inputs, noise-model evidence, conditional output columns and restarts.
See [development checks](DEVELOPMENT.md). Test datasets and generated reports
are not bundled with the repository.
