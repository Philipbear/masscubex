# MassCubeX

**One mzML or mzXML file → one compressed feature table, with MS/MS scan references.**

MassCubeX adapts [MassCube](https://github.com/huaxuyu/masscube) for large public
LC-MS repositories. It provides accelerated extraction, adaptive MS1 noise
filtering, automatic mass-tolerance calibration, profile MS1 centroiding and
parallel file processing. Dataset alignment is a separate downstream task.

## Install

Python 3.10 or newer:

```bash
git clone https://github.com/Philipbear/masscubex.git
cd masscubex
python -m venv .venv
source .venv/bin/activate
python -m pip install -e .
masscubex warmup
```

MS-Numpress-encoded mzML also needs `python -m pip install -e '.[numpress]'`.
See [installation details](docs/USAGE.md#optional-numpress-decoder), including
the Python 3.13 build instructions. Processing runs offline.

## Extract features

```bash
masscubex extract sample.mzML -o features
masscubex extract /data/study -o features --workers 4 --resume
masscubex extract --manifest files.txt -o features --workers 4 --resume
masscubex extract sample.mzML -o features --format tsv
masscubex defaults > params.json
masscubex extract --help
```

Directories are searched recursively. Manifest entries are relative to the
manifest directory; `.mzML`, `.mzXML` and their gzip forms are supported.
File workers use processes. Pools recycle by default to release retained memory;
three workers are a lower-memory alternative to the four-worker example.

## Speed compared with MassCube

Five public sample/blank/QC files, one process, Apple M2 Max:

| Matched configuration | MassCube 1.2.22 | MassCubeX | Speedup |
| --- | ---: | ---: | ---: |
| Same adaptive MS1 filter | 396.99 s | 23.97 s | **16.56×** |
| Same fixed Orbitrap intensity floors | 34.66 s | 6.89 s | **5.03×** |

Both engines use identical preprocessing and the selected ion catalog/grouping.
Calibration, profile centroiding and MS2 seeding are disabled for this comparison
because upstream lacks those extensions. MassCubeX timings include Parquet writing;
original timings exclude serialization. Scientific results match in both controls.
See [parameters, parity checks and raw measurements](benchmarks/results/UPSTREAM_SPEED.md).

With normal MassCubeX defaults, sustained four-process extraction reached
**1,921 files/hour at 2.85 GiB summed peak RSS** on this Mac. All 100 unique
feature tables and 1,800 successful repeated jobs matched the pre-optimization
version. These measurements use positive centroid mzXML; other file sizes and
acquisition modes need their own capacity measurements.
[Speed, stage timings and RAM](benchmarks/results/OVERNIGHT_OPTIMIZATION.md).

## What is inherited, and what is new?

The scientific reference is **MassCube 1.2.22**, commit
[`356748ede691862b9a87f67812db5e728a703b9d`](https://github.com/huaxuyu/masscube/tree/356748ede691862b9a87f67812db5e728a703b9d).
Trace priority/matching, peak splitting, Gaussian diagnostics, quantification and
isotope/adduct/in-source-fragment grouping follow this pinned version.

Enabled additions are adaptive MS1 noise thresholds, regional mass-tolerance
calibration, profile MS1 centroiding and labelled MS2-seeded features with raw-MS1
gap filling. The selected adduct catalog is configurable. These additions can
change feature counts relative to stock MassCube.
[Algorithm comparison](docs/UPSTREAM_COMPARISON.md) · [Complete methods](docs/METHODS.md).

## Outputs and failed files

Each successful input produces a **Zstd-compressed Parquet table** and JSON
provenance sidecar. TSV is optional. Tables contain height and area, RT in seconds,
MS/MS scan numbers and the best scan, origin labels and ion-group assignments.
Output precision is five decimals for m/z, two for RT, one for intensity and four
for correlations. Output paths include a source-path hash; `--resume` verifies
completed results before skipping them.

Corrupted or missing files are recorded while the remaining inputs continue:

```text
features/reports/<run-id>/summary.json
features/reports/<run-id>/errors.tsv
features/reports/<run-id>/files.jsonl
```

Use `--retry-errors path/to/errors.tsv` after repairing failed inputs. A batch
with file errors finishes with a nonzero exit code. Global setup/storage errors
can abort a run. [Commands, columns and restart behavior](docs/USAGE.md).

## Documentation and development

- [Calibration and gap filling](docs/AUTO_TUNING.md)
- [Adducts and in-source fragments](docs/ION_GROUPING.md)
- [Example TSV](examples/NIST_POS_Samp_14-03.features.tsv)
- [Benchmarks and validation](benchmarks/README.md)
- [Code structure and development checks](docs/DEVELOPMENT.md)

DIA classification/deconvolution, mobility-resolved extraction, library
annotation and dataset alignment are not implemented. The package and command
are named `masscubex`; the Python API differs from `masscube`.

Apache License 2.0, with retained BSD notices for SciPy/NumPy compatibility
portions. See [LICENSE](LICENSE), [NOTICE](NOTICE) and [third-party licenses](licenses).
Original MassCube attribution and Git history are retained.
