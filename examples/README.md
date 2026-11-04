# Example feature table

[NIST_POS_Samp_14-03.features.tsv](NIST_POS_Samp_14-03.features.tsv) is the complete
MassCubeX 0.15.0 output for one public MSV000086989 sample: **16,486 features**, including
**68 MS/MS-seeded features**, with all **2,596 MS/MS scans** assigned.

The table has 29 columns. m/z is rounded to five decimals; retention times are
seconds with two decimals; height, area and other intensities have two decimals.
The adjacent JSON sidecar records parameters, mass calibration, provenance and
MS/MS bookkeeping. The table itself contains all associated MS/MS scan references
and the selected best scan for each feature. Parquet remains the default format;
this example explicitly requests TSV.

Generate a TSV through the CLI:

```bash
masscubex extract NIST_POS_Samp_14-03.mzXML -o features --format tsv
```

CLI output names include an input-path hash. This example was named explicitly
through the Python API's `process_file(source, output)` argument.

Source: [MSV000086989 dataset](https://massive.ucsd.edu/ProteoSAFe/dataset.jsp?task=2f73277b4e034948acebfdf1edab17ed).
Input SHA256: `458edf67e81118f19c9beb4335b80bf6927fa5a8871ba2ae9475e3321537a128`.
Output SHA256: `3117149159a87a8d8092b4ee82709e445626e82104d6ccf009033407be5e26ce`.
Implementation SHA256: `a0912eedb3bec2680077841a7f5ba9b85c8477e4b562d43b7b1e718e46885aa7`.
