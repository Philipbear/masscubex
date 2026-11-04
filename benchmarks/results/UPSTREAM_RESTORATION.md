# Upstream restoration — 2026-11-04

This report validates the restored MassCube scientific core, including short traces, last-mass matching, ten missing scans, apex mass, original quantification and MS2 assignment. The runs below disable auto-tuning, profile MS1 centroiding and MS2-seeded gap filling; adaptive MS1 noise and the selected catalog are controlled differences. These three extensions were subsequently enabled by default at the user's request. The figures below remain a core-algorithm comparison, not counts for the current default pipeline. See the [current README](../../README.md).

## Complete-row comparisons

| File | Original, fixed MS1 50,000 | MassCubeX, fixed MS1 50,000 | Both with adaptive MS1 filter | Adaptive rows with ≥5 observations |
| --- | ---: | ---: | ---: | ---: |
| NIST_POS_Samp_01-01.mzXML | 28,492 | 28,492 | 70,776 | 15,754 |
| NIST_POS_Samp_14-03.mzXML | 20,865 | 20,865 | 75,100 | 17,928 |
| Blank_00.mzXML | 12,142 | 12,142 | 55,614 | 10,286 |
| Blank_Post_14-03.mzXML | 5,519 | 5,519 | 114,323 | 10,622 |
| QC_Mix_01.mzXML | 9,225 | 9,225 | 54,872 | 10,675 |
| QC_Mix_04.mzXML | 10,446 | 10,446 | 84,163 | 14,175 |
| **Total** | **86,689** | **86,689** | **454,848** | **79,440** |

All 12 file/mode pairs have zero mismatched rows for exported mass, apex/start/end times, height, area, observation count, Gaussian similarity, all assigned MS2 source indices and best-scan choice. The fixed runs use MS2 floor 10,000; adaptive runs use 500, corresponding to the supplied configuration. Full grouping labels were additionally compared for NIST_POS_Samp_14-03 in both modes, with zero differences.

The adaptive column is an explicit control: the original feature algorithms receive the same adaptive MS1 thresholded scans. Stock MassCube does not implement that adaptive filter. This verifies the intended distinction rather than implying that stock defaults produce 454,848 rows.

375,408 of those adaptive rows have fewer than five nonzero observations and are retained as upstream does. Counts do not establish chemical validity. The former local example had 19,847 rows under altered tracing/quantification/calibration rules; this comparison's example has 75,100 rows, including 57,172 short traces, before enabling automatic calibration and MS2 seeding.

## Methods and limitations

- Reference: huaxuyu/masscube commit `356748ede691862b9a87f67812db5e728a703b9d`; function bodies are retained with attribution and checked by source hashes.
- Input adapter: Pyteomics decodes the public mzXML files because the pinned original reader rejects that extension. Original numeric arrays/time values feed unchanged upstream detection, segmentation, finalization, association and representative selection functions.
- Checks are at exported precision (m/z 5, seconds 2, intensity 1, correlation 4 decimals); randomized synthetic tracker tests also compare full measurements/trace intensities. Output formatting and per-polarity repository support are intentional adapters.
- Grouping comparisons use the selected adduct catalog. The default positive stream assumes [M+H]+ for each ungrouped anchor; no pairwise alternative dimer hypotheses are generated.
- The scientific source functions are currently called directly for splitting/finalization. The earlier throughput reports no longer apply. The recorded elapsed times are single diagnostic runs, not repeated production benchmarks.
- No production parameters were silently retuned to match counts. No push, release or full release suite was performed.

[Detailed JSON report](upstream-restoration-2026-11-04.json) · [Reproduction script](../verify_upstream_defaults.py) · [Current algorithms](../../docs/METHODS.md)
