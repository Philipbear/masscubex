# MassCubeX

Fast extraction of one feature table per mzML or mzXML file, built on
[MassCube 1.2.22](https://github.com/huaxuyu/masscube/tree/356748ede691862b9a87f67812db5e728a703b9d).
MassCubeX supports profile MS1 data, polarity switching, adaptive MS1 noise
filtering, automatic mass-tolerance tuning, MS/MS scan association and ion grouping.
Dataset alignment is a separate downstream task.

## Install

Requires Python 3.10 or newer.

```bash
git clone https://github.com/Philipbear/masscubex.git
cd masscubex
python -m venv .venv
source .venv/bin/activate
python -m pip install -e .
masscubex warmup
```

On Windows, activate with `.venv\Scripts\activate`. MS-Numpress input additionally
requires `python -m pip install -e '.[numpress]'`; see the
[installation guide](docs/USAGE.md#optional-numpress-decoder) for Python 3.13.

## Command-line usage

Extract one file, a directory tree, or a manifest containing one input path per line:

```bash
masscubex extract sample.mzML -o features
masscubex extract /data/study -o features --workers 4 --resume
masscubex extract --manifest files.txt -o features --workers 4 --resume
```

Directories are searched recursively. Manifest paths are relative to the manifest
file. `.mzML`, `.mzXML`, and their gzip forms are accepted. Each worker processes
an independent file; reduce `--workers` for large files or limited RAM.
`--resume` verifies completed outputs and reprocesses stale or incomplete results.

Parquet with Zstd compression is the **default**. To request TSV:

```bash
masscubex extract sample.mzXML -o features --format tsv
```

### Choose output columns

The default table has 15 columns, including **peak height**, RT boundaries,
all associated MS/MS scan numbers, the best MS/MS scan, feature origin and compact
adduct/isotope/in-source-fragment labels. `scan_count` is not exported.

| Option | Effect |
| --- | --- |
| `--include-peak-area` | Add peak area alongside height |
| `--quant-method peak_area` | Export area instead of height |
| `--include-ion-relations` | Add detailed ion-group relationships |
| `--include-apex-noise-threshold` | Add the local apex noise threshold |
| `--include-gaussian-similarity` | Add the Gaussian peak-shape score |

```bash
masscubex extract sample.mzML -o features --include-peak-area --include-ion-relations
```

RT is in seconds, rounded to two decimals. m/z has five decimals, intensities
have one, and correlation/shape scores have four. Tables have source-path hashes
in their filenames to avoid collisions, plus a JSON sidecar containing parameters,
calibration, acquisition information and exceptional MS/MS identifier lookups.

### Noise, grouping and saved parameters

The default MS1 noise threshold is μ + 3σ from the lowest 10% of positive
intensities in fixed 100 Da bins, with a minimum of five background points.
The optional consistency estimator also uses μ + 3σ, identifies repeated
intensity populations and protects ions supported by neighboring scans or MS/MS:

```bash
masscubex extract sample.mzML -o features --noise-estimator consistency
masscubex adducts
masscubex extract sample.mzML -o features --adducts 1-15,22,23
masscubex defaults > params.json
masscubex extract /data/study -o features --config params.json --workers 4 --resume
masscubex extract --help
```

Explicit flags override JSON settings. Auto-tuning, profile centroiding,
MS/MS-seeded gap filling and ion grouping are enabled by default.

## Failed files and DIA

Corrupt or missing files are recorded and processing continues. Confidently
recognized DIA acquisitions, including AIF/MSᴱ and SWATH, are skipped by default.
Unknown acquisition modes are processed with their uncertainty recorded.

Each run writes `features/reports/<run-id>/` containing `summary.json`,
`errors.tsv`, `skipped.tsv`, and `files.jsonl`. A batch with file errors finishes
with a nonzero exit code; global setup or storage failures can abort a run.
After repairing inputs, retry the error list:

```bash
masscubex extract --retry-errors features/reports/RUN_ID/errors.tsv -o features --resume
```

## Speed

On five centroid mzXML files on an Apple M2 Max, matched extraction was
**16.6× faster with adaptive MS1 noise** and **5.0× with fixed intensity floors**
than MassCube 1.2.22. Calibration, profile centroiding and MS/MS seeding were
disabled for these comparisons; speed depends on input and settings.

## Methods and license

Tracing, peak splitting, quantification and ion grouping follow the pinned
MassCube version. Adaptive MS1 noise, tolerance calibration, profile centroiding,
MS/MS-seeded gap filling and repository execution are MassCubeX additions.
These extensions can change feature counts relative to stock MassCube.

[CLI reference](docs/USAGE.md) · [Methods](docs/METHODS.md) ·
[Calibration](docs/AUTO_TUNING.md) · [Ion grouping](docs/ION_GROUPING.md) ·
[DIA detection](docs/ACQUISITION.md) · [Development](docs/DEVELOPMENT.md)

Apache License 2.0, with retained BSD notices for NumPy/SciPy compatibility
code. See [LICENSE](LICENSE), [NOTICE](NOTICE), and [third-party licenses](licenses).
