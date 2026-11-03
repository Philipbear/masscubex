# Extraction performance: MassCubeX 0.7.0

**The measured single-worker speedup is 2.46×. The requested further 10× end-to-end target was not reached.**
The comparison uses the same 50 public files, scientific defaults, Parquet format,
and machine. The 0.6.0 baseline was re-run from commit `6d70109`; production code
for 0.7.0 is commit `7564ce9`. Timing excludes correctness audits and JIT warmup.

## Complete extraction, including output

| Version / workers | Median for 50 files | Seconds/file | Tables |
| --- | ---: | ---: | ---: |
| 0.6.0, one worker, freshly re-run | 147.64 s | 2.953 | 85.72 MB |
| 0.7.0, one worker | 60.11 s | 1.202 | 85.71 MB |
| 0.7.0, four process workers | 17.73 s | 0.355 | 85.71 MB |

One-worker runtime falls 59.3%.
Four workers provide 3.39× the throughput of one current worker,
including worker startup. This parallel scaling is separate from the code speedup.
The previous four-worker result was 42.19 s; see the [0.6.0 report](PERFORMANCE-v0.6.0.md).
Do not combine parallel scaling with the single-worker result to claim a further
10× algorithmic speedup. MB means 1,000,000 bytes.

Three complete passes: baseline **147.64, 147.63, 148.12 s**; current one-worker
**59.64, 60.11, 60.53 s**; current four-worker **17.52, 17.73, 18.29 s**.
Fresh runs overwrite output and decode each input again; results are never
reused across files or benchmark repetitions. The new scan cache only replays
raw decoded data between stages of the *same* extraction call.

Peak process RSS was 615.3 MiB for the baseline and
594.8 MiB for the current one-worker run. These are whole-process
high-water marks, not memory guarantees for arbitrary inputs. The scan cache
budget defaults to 64 MiB per worker; trace histories, feature records, fragment
evidence and library buffers require additional memory.

Raw reports: [baseline](extraction-v0.6.0-recheck-parquet.json),
[current](extraction-v0.7.0-parquet.json), [four workers](extraction-v0.7.0-parquet-4-workers.json).

## Remaining bottlenecks

Production timers from the last complete single-worker pass:

| Stage | 0.6.0 | 0.7.0 | Current elapsed share |
| --- | ---: | ---: | ---: |
| MS1 extraction | 82.39 s | 21.26 s | 35.1% |
| Calibration, including the physical XML/binary read | 21.41 s | 14.49 s | 23.9% |
| Parquet writing and fsync | 9.52 s | 9.74 s | 16.1% |
| Ion grouping | 17.02 s | 8.45 s | 14.0% |
| Raw-MS1 seeded gap filling | 13.60 s | 3.66 s | 6.0% |
| MS2 association | 3.42 s | 2.19 s | 3.6% |

The remaining time is checksums, completion metadata, locks and orchestration.
MS1 extraction remains the largest aggregate stage; calibration/decoding and
output are substantial too. Reaching 10× against the baseline would require
about 14.76 s for all 50 files on one worker. The current calibration and
output stages alone exceed that budget. The present work does not establish
that target for arbitrary public repositories.

Separate profiling on one biological sample, one blank and one QC, with three
repetitions each, gives:

| Extraction substep | 0.6.0 mean/file | 0.7.0 mean/file | Speedup |
| --- | ---: | ---: | ---: |
| Chromatographic segmentation and peak measurement | 1.103 s | 0.172 s | 6.43× |
| Adaptive MS1 noise estimation | 0.193 s | 0.140 s | 1.37× |

These are substeps of MS1 extraction, **not additional time**. The new shared
scan-reading/replay timer spans calibration, extraction and gap filling; it
cannot be added to their totals. Parquet writing is effectively unchanged at
this measurement scale. The columnar writer refactor improves separation and
avoids per-row dictionaries; it is not credited with a measured speedup.
[Baseline profile](profile-v0.6.0-recheck-parquet.json), [current profile](profile-v0.7.0-parquet.json).

## Implementation and readability

- Fused Numba chromatogram resampling, reflected Gaussian smoothing, plateau,
  distance/prominence selection and original-sample measurements. Exact/near
  distance-priority ties retain the SciPy compatibility path. No fast-math.
- Compiled calibration screening rejects poor candidates cheaply. Accepted
  candidates retain the original NumPy/SciPy measurement arithmetic, preserving
  all calibrated mass/time settings on the full cohort.
- Bounded decoded-scan replay avoids repeated XML/base64/zlib work. Overflow
  discards the cache and falls back to streaming. Calibration drops fragment
  arrays immediately instead of retaining them outside the cache budget.
- Low-tail noise avoids temporary tree removal/restoration when a candidate is
  outside the selected background prefix. The 50 Da, lowest 10%, minimum five,
  leave-one-out μ + 3σ rule is unchanged.
- Linear sweeps find sorted adduct mass windows; unusual nonmonotone bounds use
  binary searches. Gap filling packs each bounded raw query once before its
  compiled EIC searches. Normal spectrum validation uses one compiled pass.
- Empty MS2 provenance is not repeatedly rebuilt. Parquet columns are formed
  directly, with finite-value and evidence-schema validation preserved.
- Numerical modules have explicit responsibilities, named result records and
  documented tracker-state columns. TSV/Parquet writers are separate functions.
  See [code structure](../../docs/DEVELOPMENT.md) and [methods](../../docs/METHODS.md).

## Scientific equivalence

The [full comparison](extraction-equivalence-v0.7.0.json) checks **1,540,612 rows**,
including **5,143 MS2-seeded features**, **125,266 associated MS2 scans**, and
**205,868 ion-relationship entries**, plus every scientific sidecar field.
Feature order/identities, heights, boundaries, counts, MS2 lists/best scans,
relation partners/types, correlation evidence and calibration values match
exactly. No scientific stage or default was disabled for the benchmark.

The outputs are **not byte-identical across versions**: different float64
summation/exp evaluation can change final bits. The largest observed relative
area difference is 3.68e-15; the largest absolute
m/z difference is 4.55e-13 Da. Derived relation mass errors differ by
at most 8.22e-10 ppm. Noise-threshold differences are at most
3.62e-16 relative, without a change in retained centroid counts.
The audit permits 3e-14 relative/near-zero absolute error only in listed
measurement fields, and 2e-7 ppm absolute error in derived mass errors;
identities, boundaries, heights and scan lists require exact equality.

## Conditions and reproduction

Apple M2 Max, 32 GiB RAM, macOS 26.6.2 ARM64, Python 3.10.12. The corpus is
[MSV000086989](https://massive.ucsd.edu/ProteoSAFe/dataset.jsp?task=2f73277b4e034948acebfdf1edab17ed):
40 biological samples, six blanks and four QCs; 604,070,716 bytes. All are
positive-mode centroid mzXML. Profile data, mzML and polarity switching have
synthetic regression coverage; these timings do not establish their throughput.
Local SSD, warm filesystem/JIT caches; worker startup included in parallel runs.
No simultaneous CPU-heavy tests, benchmarks or audits. Every report records the
scientific configuration, dependency versions and implementation fingerprint.

```bash
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 NUMBA_NUM_THREADS=1
mkdir -p validation/optimization-v0.7/baseline
git archive 6d70109 | tar -x -C validation/optimization-v0.7/baseline
PYTHONPATH=validation/optimization-v0.7/baseline/src python benchmarks/benchmark_extraction.py \
  --format parquet --repeats 3 --output-directory validation/optimization-v0.7/baseline-parquet \
  --report benchmarks/results/extraction-v0.6.0-recheck-parquet.json
python benchmarks/benchmark_extraction.py --format parquet --repeats 3 \
  --output-directory validation/optimization-v0.7/final-parquet \
  --report benchmarks/results/extraction-v0.7.0-parquet.json
python benchmarks/compare_extraction_versions.py \
  --reference validation/optimization-v0.7/baseline-parquet \
  --candidate validation/optimization-v0.7/final-parquet \
  --report benchmarks/results/extraction-equivalence-v0.7.0.json
```

Use `--workers 4` with a separate output/report path for parallel measurement.
`profile_extraction.py` takes the same format and directory/report arguments.
See [validation](VALIDATION.md) for corruption, raw-data and cross-platform checks.
