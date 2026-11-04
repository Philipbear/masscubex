# Feature counts before full upstream restoration — 2026-11-04

This report describes the earlier local implementation. See [the restored upstream defaults and new comparison](UPSTREAM_RESTORATION.md) for the current pipeline.

Fresh runs of six files from MSV000086989: two biological samples, two blanks and two QCs. The local MassCubeX implementation includes the restored upstream splitting and greedy grouping rules. No extraction settings or algorithms were changed for this audit.

| File | MassCube total | MassCube ≥5 scans | MassCubeX detected MS1 | MassCubeX MS2-seeded | MassCubeX total |
| --- | ---: | ---: | ---: | ---: | ---: |
| NIST_POS_Samp_01-01.mzXML | 28,492 | 7,306 | 17,020 | 147 | 17,167 |
| NIST_POS_Samp_14-03.mzXML | 20,865 | 7,109 | 19,790 | 57 | 19,847 |
| Blank_00.mzXML | 12,142 | 3,887 | 13,439 | 63 | 13,502 |
| Blank_Post_14-03.mzXML | 5,519 | 1,932 | 14,188 | 26 | 14,214 |
| QC_Mix_01.mzXML | 9,225 | 2,975 | 14,381 | 46 | 14,427 |
| QC_Mix_04.mzXML | 10,446 | 2,736 | 17,059 | 38 | 17,097 |
| **Total** | **86,689** | **25,945** | **95,877** | **377** | **96,254** |

“≥5 scans” means at least five nonzero MS1 observations. Original MassCube retains shorter traces in its single-file table; MassCubeX requires five observations for detected rows. The common ≥5-observation count is a comparison diagnostic, not a claim about upstream’s output filter or downstream alignment filters. MS2-seeded rows are reported separately.

Using that common observation criterion, MassCubeX detects **3.70×** as many features overall: **2.55×** in samples, **4.75×** in blanks and **5.51×** in QCs. Raw total rows increase only 11.0%, because upstream includes 60,744 shorter traces. The 377 seeded rows are 0.39% of the MassCubeX total and do not explain the excess.

## Intensity-floor diagnostic

The instrument metadata identifies a Thermo Q Exactive. Original MassCube’s Orbitrap default MS1 threshold is 50,000. In an additional diagnostic run, only MassCubeX’s adaptive noise filter was replaced by a fixed 50,000 floor; learned calibration and every other effective processing setting were verified unchanged.

| Group | Original MassCube ≥5 scans | Default MassCubeX detected | MassCubeX detected with fixed 50,000 floor |
| --- | ---: | ---: | ---: |
| biological | 14,415 | 36,810 | 15,024 |
| blank | 5,819 | 27,627 | 6,876 |
| QC | 5,711 | 31,440 | 7,001 |
| **Total** | **25,945** | **95,877** | **28,901** |

The fixed-floor control reduces detected MassCubeX features by **69.9%**, leaving **11.4% more** than original MassCube at the common observation criterion. This strongly implicates the different noise policy in the remaining count difference. It does not establish that the extra features are false positives, or justify adopting one fixed cutoff across repositories. Tracing, calibration and centroid preprocessing still differ between implementations; their residual effects were not isolated here. Production defaults remain unchanged.

## Reproducibility and scope

- Original source: MassCube 1.2.22, commit `356748ede691862b9a87f67812db5e728a703b9d`. The Python sources were copied byte-for-byte from the pinned Git object. Fresh counts exactly reproduce the earlier pinned audit for these six inputs.
- This upstream version rejects mzXML directly. The existing metadata-only adapter uses Pyteomics to pass original arrays and retention times in minutes to upstream’s mzML scan extractor. Preprocessing, tracing, two-pass segmentation, finalization and MS2 allocation use unchanged upstream functions. Grouping changes labels, not the row set, so upstream counting ends after finalization.
- MassCubeX: local schema 10; implementation SHA256 `122cca9e40bf12ce996875ae4c406da2368be91a78881ce50fa00d2177d45cfa`. Grouping remains enabled, adaptive MS1 noise uses fixed 100 Da bins and lowest 10% with μ + 3σ, and splitting retains the restored upstream rules.
- All six fresh original runs and twelve MassCubeX runs completed. No release suite or parameter optimization was performed. Temporary sources, feature tables and isolated dependency files were removed after retaining the compact report.
- These are counts for one positive-mode centroid Orbitrap dataset; they do not measure identification accuracy or biological validity.

[Full counts, input/source hashes, parameters and diagnostics](feature-counts-current-2026-11-04.json). Reproduction: [compare_current_feature_counts.py](../compare_current_feature_counts.py), using the upstream-only dependencies on `PYTHONPATH` as described in [audit_upstream_features.py](../audit_upstream_features.py).

```bash
PYTHONPATH=/path/to/isolated-upstream-dependencies python benchmarks/compare_current_feature_counts.py \
  --report benchmarks/results/feature-counts-current-2026-11-04.json
```
