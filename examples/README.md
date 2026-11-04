# Example feature table

[NIST_POS_Samp_14-03.features.tsv](NIST_POS_Samp_14-03.features.tsv) is the complete
current default output for one public MSV000086989 sample: **93,731 rows**,
comprising **93,722 detected MS1 features and 9 MS2-seeded features**.
Short traces are retained under upstream rules: 74,478 rows have fewer than five
nonzero observations. **All 2,596 MS2 scans are assigned; none are unassigned**.
Auto-tuning, profile centroiding and MS2 seeding are enabled; this particular
input contains centroided scans, so it does not exercise profile mass detection.

The table has 25 columns. m/z has five decimals, RT is seconds with two decimals,
intensities have one decimal and correlations have four. Detected mass comes
from the apex; height and area follow upstream integer quantification.
`adduct`, group/isotope labels and the source-fragment flag expose assignments.
The default [M+H]+ label is an assumption, not a confirmed chemical identification.
Redundant seed precursor/time and best-MS2 measurement columns remain absent.

The scientific core follows MassCube 1.2.22, commit
`356748ede691862b9a87f67812db5e728a703b9d`, with enabled MassCubeX additions.
The explicit upstream-baseline run (auto-tuning/seeding disabled) has 75,100 rows
and matches the original-function reference given the same adaptive MS1 filter.
That parity check does not claim that the calibrated default output matches
stock MassCube. [Core comparison](../benchmarks/results/UPSTREAM_RESTORATION.md).

For this default run, all three positive-polarity MS1 regions have independent
fits, each using 100 screened representatives:

| m/z region | Absolute term (Da) | Relative term (ppm) | Example actual ±radius |
| --- | ---: | ---: | --- |
| 0–200 | 0.000500 | 5 | 0.000750 Da at m/z 150 |
| 200–500 | 0.000500 | 5 | 0.001750 Da at m/z 350 |
| 500+ | 0.002435 | 5 | 0.004000 Da at m/z 800 |

The radius is the maximum of the absolute and relative terms. These are the
values used by trace matching, not the pilot search window. There is no learned
elapsed-time trace cutoff; ten consecutive missing same-polarity MS1 scans are
allowed. Full fitted values, evidence and timings are in the adjacent JSON.

Reproduce the defaults, or disable the added stages explicitly:

```bash
masscubex extract NIST_POS_Samp_14-03.mzXML -o features --format tsv
masscubex extract sample.mzML -o features \
  --no-auto-tune --profile-mode reject --no-ms2-seeded-features
```

Parquet remains the default output. CLI names include an input-path hash; this
example was named explicitly with the Python `process_file(source, output)` API.
Source: [MSV000086989 dataset](https://massive.ucsd.edu/ProteoSAFe/dataset.jsp?task=2f73277b4e034948acebfdf1edab17ed).
Input SHA256: `458edf67e81118f19c9beb4335b80bf6927fa5a8871ba2ae9475e3321537a128`.
Output SHA256: `15a83773ef2b2bdbb7e0dcc119de4ca503589e61592c5e6caf641eaa864e91f9`.
Implementation SHA256: `fca3d76939fc32766261a865252ba33710f8f90f5abbe89ff5790748268d8a8e`.
