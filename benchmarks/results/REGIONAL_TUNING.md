# Regional mass calibration — MassCubeX 0.12.0

The default now selects calibration anchors with symmetric Gaussian similarity
≥0.90, divides each polarity's acquired m/z span into five equal-width regions,
and uses mean + 3 sample SD after robust outlier screening. Both MS1 and precursor
fits require 20 distinct qualifying peaks per fitted region. The detailed
[algorithm and all fallbacks](../../docs/AUTO_TUNING.md) specify the error
summaries, acquisition-range metadata, floors/caps and override rules.

## What is supported by the 50-file sample

All 50 raw files are from MSV000086989: 6 blanks, 4 QCs and 40 biological samples,
positive centroid mzXML. Shape screening leaves enough global MS1 evidence in
49 files. Blank_00 has 19 accepted representatives and retains configured mass
settings. Global precursor evidence is sufficient in 45 files.

| Equal-width region, low to high m/z | Files with local MS1 fit | Files with local precursor fit |
| --- | ---: | ---: |
| 1 | 49 | 34 |
| 2 | 8 | 0 |
| 3 | 0 | 0 |
| 4 | 0 | 0 |
| 5 | 0 | 0 |

Thus 57 of 250 MS1 regions have independent fits; 193 use an explicit file-wide
fallback. Five bins do not manufacture evidence in sparse high-mass regions.
These files omit explicit acquisition limits, so bounds come from all observed
MS1 arrays before peak selection. Gaussian shape and isolation checks do not
prove a single chemical identity. Negative-mode and polarity-switching behavior
are tested synthetically, not established by this single repository dataset.

For QC_Mix_01, the observed span is approximately 100–1499.71 m/z:

| Approximate m/z interval | Accepted MS1 / MS2 peaks | MS1 Da component, 3 SD | MS1 Da component, 5 SD | Source |
| --- | ---: | ---: | ---: | --- |
| 100–379.94 | 62 / 24 | 0.000500 | 0.000673 | Regional |
| 379.94–659.88 | 8 / 3 | 0.000542 | 0.000811 | File-wide fallback |
| 659.88–939.83 | 0 / 0 | 0.000542 | 0.000811 | File-wide fallback |
| 939.83–1219.77 | 0 / 0 | 0.000542 | 0.000811 | File-wide fallback |
| 1219.77–1499.71 | 0 / 0 | 0.000542 | 0.000811 | File-wide fallback |

The ppm component is 5 in each of these cases. The actual half-window is the
larger of the Da component and m/z × ppm × 1e−6; the high-mass rows therefore
have identical physical windows at 3 SD and 5 SD. For example, at m/z 500,
both use ±0.0025 Da. Exact region edges are in the JSON report.

## Feature counts and 3 SD versus 5 SD

| Group | Files | v0.11 detected | 3 SD detected | 5 SD detected |
| --- | ---: | ---: | ---: | ---: |
| blank | 6 | 75,146 | 74,972 | 74,971 |
| QC | 4 | 57,625 | 57,002 | 56,983 |
| biological | 40 | 610,654 | 611,661 | 611,659 |

At the 3 SD default, there are **743,635 detected + 4,301 seeded = 747,936 rows**.
The prior version had 747,745 total rows: a 0.026% increase. The Gaussian gate
also changes which peaks set chromatography timing parameters; these count
changes are not attributable solely to mass windows. The shape guard and
adaptive MS1-noise algorithm themselves are unchanged.

Using 5 SD gives 747,912 total rows, only 24 fewer than 3 SD in this dataset.
A wider mass window can join observations into existing traces or assign an
MS2 to an existing feature, so feature counts need not increase monotonically.
The 5 ppm lower bound dominates most of these data. This comparison does not
show that 3 SD is universally sufficient for weak ions or other instruments.
It supports retaining 3 SD as the chosen starting default, with 5 SD available
through `--tuning-mass-sigma 5`. The margin operates on per-peak Q95 absolute
errors; it does not imply a normal-distribution coverage guarantee.

## Validation and processing time

- **676 source tests and 676 isolated installed-wheel tests pass**, plus lint.
  Tests cover Gaussian rejection, heterogeneous regional mass scatter, sample
  SD arithmetic, sparse/unknown-polarity fallbacks, metadata and observed range
  sources, boundaries, matching across edges, precursor-intensity lookup,
  distinct MS2 evidence, manual overrides and strict gap-fill caps.
- Independent Pyteomics decoding verifies all **125,266 MS2 scans exactly once**:
  120,965 attach to detected features and 4,301 to seeded features. No MS2 is lost.
- Raw seeded quantification passes: 4,142 measured seeds, 122 flagged noise-height
  imputations with missing area, 35 with no local peak and 2 with no local ion.
- All **747,936 typed Parquet and TSV rows agree**. Precision/schema stay unchanged.

A warm serial pass over all 50 files took **47.77 s**, excluding
independent audits/downloads. Three separate warm serial passes over the ten
blank/QC files have a median of **9.86 s**.

| Stage, ten blank/QC files | Median seconds across three passes |
| --- | ---: |
| calibration | 3.304 |
| feature_extraction | 4.077 |
| ion_grouping | 0.798 |
| ms2_association | 0.247 |
| ms2_gapfill | 0.316 |
| table_write | 1.129 |

Measurements use the existing M2 Max environment with warmed JIT/cache and no
other benchmark/test process running. Historical v0.11 times were 47.70 s for
50 files and 9.85 s for ten blank/QC files; those are archived runs, not a paired
speed test. No speed improvement or regression is inferred from that difference.

Evidence: [regional comparison](regional-tuning-v0.12.0.json),
[independent scan/quantification audit](MSV000086989-ms2-v0.12.0.json), and
[release fingerprint, wheel checksum and tests](release-v0.12.0.json).
Reproduce with `benchmarks/validate_ms2_coverage.py --workers 4`, followed by
`benchmarks/validate_regional_tuning.py` serially. Raw data, complete sidecars,
installed wheel and test XML remain in `validation/regional-tuning-v0.12/`
and `validation/MSV000086989/v0.12.0-ms2/` rather than Git.
