# Development

Install the development dependencies and run the checks:

```bash
python -m pip install -e '.[dev]'
ruff check src tests
python -m pytest -q
python -m build
```

Tests generate small mzML/mzXML fixtures in temporary directories. Downloaded
raw datasets, feature tables and benchmark reports are not required or bundled.
The CI matrix tests source and installed-wheel behavior on Linux, macOS and
Windows with Python 3.10–3.13.

## Code structure

- `pipeline.py` owns calibration, extraction, grouping and atomic table/JSON publication.
- `cli.py` schedules independent files and records failures and DIA exclusions.
- `readers.py`, `formats.py` and `acquisition.py` decode and screen input files.
- `upstream.py` contains the pinned MassCube scientific reference bodies;
  `upstream_adapter.py` connects them to streaming input.
- `features.py`, `summary.py`, `measurements.py` and `chromatograms.py` implement
  equivalent trace, splitting and measurement execution.
- `noise.py`, `noise_regions.py` and `centroid.py` implement default noise and mass detection.
- `consistency_noise.py`, `consistency_kernels.py`, `consistency_cache.py` and
  `profile_noise.py` implement the optional intensity-consistency estimator.
- `tuning.py`, `calibration.py`, `pilot_sampling.py` and `mass_tolerance*.py`
  select qualified representatives and estimate regional tolerances.
- `ms2.py`, `gapfill.py` and `gapfill_refinement.py` preserve scan ownership and
  measure unmatched precursors.
- `grouping*.py` implement inherited isotope, adduct and source-fragment assignments.
- `output.py` selects per-file columns and writes rounded TSV/Parquet tables.

## Numerical and publication rules

Keep the pinned upstream bodies in `src/masscubex/upstream.py` and
`tests/upstream_masscube.py` unchanged, including their attribution. The saved
function hashes and independent reference comparisons guard those algorithms.
Preserve strict windows, acquisition and tie order, polarity isolation,
float32 scan representation and original measurement arithmetic. Round only
when serializing. Numba kernels must not enable fast-math.

After changing a helper used by compiled modules, remove the generated Numba
`.nbc`/`.nbi` caches under `src/masscubex/__pycache__` before verification.
Never package caches. For speed measurements, use identical scientific settings,
exclude JIT warmup, time full passes and compare values and scan ownership.
A kernel improvement is not a complete-workflow speedup.

Run focused checks during development. Run the full source suite, lint and
package checks before an explicitly requested publication. CI verifies the
remaining platform matrix after a push.
