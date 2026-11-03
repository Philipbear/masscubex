<!-- Modified for MassCubeX; replaces the upstream MassCube README. -->
# MassCubeX

**One mzML or mzXML file → one feature TSV, with MS/MS scan references.**
MassCubeX rewrites [MassCube](https://github.com/Philipbear/masscube) for public
LC-MS repositories: streaming input, Numba-accelerated extraction, profile MS1
centroiding, adaptive noise filtering, and independent file workers.
Dataset alignment remains a separate downstream task.

## Install

Python 3.10 or newer, with a Python/NumPy combination supported by Numba:

```bash
git clone https://github.com/Philipbear/masscubex.git
cd masscubex
python -m venv .venv
source .venv/bin/activate
python -m pip install -e .
masscubex warmup
```

For MS-Numpress-encoded mzML, install `python -m pip install -e '.[numpress]'`.
Processing runs offline. The package, import and command are named `masscubex`;
this fork is not API-compatible with `masscube`.

## Extract features

```bash
# One file
masscubex extract sample.mzML -o features

# A directory tree, or a manifest with one path per line
masscubex extract /data/study -o features --workers 8 --resume
masscubex extract --manifest files.txt -o features --workers 8 --resume

# Important options
masscubex extract sample.mzML -o features --quant-method peak_area
masscubex extract sample.mzML -o features --no-auto-tune
masscubex extract sample.mzML -o features --no-group-ions
masscubex defaults > params.json
masscubex extract /data/study -o features --config params.json --workers 8
masscubex extract --help
```

Both `.mzML` and `.mzXML`, including `.gz`, are supported. Manifest paths are
relative to the manifest directory unless absolute. Workers use processes by
default; `--executor thread` is also available.

## Defaults and outputs

- **Quantification:** peak height; both height and area are always included.
- **Automatic calibration:** broad initial mass search, then at least 20 qualified
  peaks per polarity to learn mass and time settings; explicit fallbacks when
  evidence is insufficient. Polarity-switching files use separate streams.
- **MS1 noise:** 50 Da window, lowest 10% of positive neighbors, minimum 5,
  threshold μ + 3σ. Profile MS1 is centroided automatically.
- **MS/MS:** all associated scan IDs and the best scan, selected by summed
  fragment intensity. The isolation-width cutoff is disabled.
- **Unmatched MS/MS:** create a marked precursor feature and gap-fill raw MS1
  height and area using local peak boundaries when calibration is supported.
  The fallback search is ±5 seconds. Ambiguous measurements remain labelled. Unusable
  precursor metadata is reported with scan IDs and reasons in the JSON sidecar.
- **Ion grouping:** on by default; selected adducts and water-loss candidates are
  listed by `masscubex adducts`.

Each input produces a feature TSV and a JSON sidecar with settings, provenance,
counts and unresolved MS/MS scans. Output names include an input-path hash to
avoid basename collisions. `--resume` verifies completed outputs before skipping.
See the [usage and output reference](docs/USAGE.md) for columns, scan lookup,
Python usage and configuration options.

Corrupted or missing input files **do not stop the remaining files**. Each run saves:

```text
features/reports/<run-id>/summary.json
features/reports/<run-id>/errors.tsv
features/reports/<run-id>/files.jsonl
```

Retry repaired inputs using `--retry-errors path/to/errors.tsv` with the original
processing settings. A batch with file failures returns a nonzero exit code after
processing the other files. Global setup/storage failures can still abort a run.

## Documentation

- [Implemented methods and changes from MassCube](docs/METHODS.md)
- [Automatic calibration, gap filling and polarity switching](docs/AUTO_TUNING.md)
- [Commands, output columns, MS/MS lookup and restart behavior](docs/USAGE.md)
- [Adduct and in-source-fragment candidates](docs/ION_GROUPING.md)
- [Tests, real-data validation and benchmarks](benchmarks/results/VALIDATION.md)

DIA classification/deconvolution, mobility-resolved feature extraction, library
annotation and dataset alignment are not implemented.

## Development and license

```bash
python -m pip install -e '.[dev,numpress]'
pytest -q
ruff check src tests benchmarks
python -m build
```

Apache License 2.0; see [LICENSE](LICENSE) and [NOTICE](NOTICE). Original MassCube
attribution and Git history are retained.
