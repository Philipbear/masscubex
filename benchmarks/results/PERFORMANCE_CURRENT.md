# Historical MassCubeX 0.15.0 processing speed

These measurements predate restoration of the pinned upstream tracing,
short-trace retention, splitting and quantification rules. They do **not**
describe the current workflow, despite the unchanged development version number.
See the [current 100-file speed/RAM study](OVERNIGHT_OPTIMIZATION.md) and
[five-file upstream benchmark](UPSTREAM_SPEED.md). The original
measurements and projections below are retained as historical results only.

Measured 2026-11-04 on an Apple M2 Max (32 GB), macOS 26.6.2, Python 3.10.12.
The same 50 public MSV000086989 positive centroid-mzXML files total 604.1 MB.
That implementation's defaults included adaptive MS1 noise, policy-5 calibration, MS2 association
and seeded gap filling, ion grouping, and compressed Parquet output.

| Workers | 50 files, median wall time | Effective seconds/file | Files/hour | 2 million similar files |
|---|---:|---:|---:|---:|
| 1 | 55.39 s | 1.108 | 3,249 | 25.65 days |
| 4 | 17.15 s | 0.343 | 10,493 | 7.94 days |

Three complete serial passes: 55.39 s, 55.55 s, 55.28 s.
Three four-worker passes: 17.15 s, 17.26 s, 17.15 s.
Four workers provide 3.23× measured batch throughput.

## Where time goes

Median sums of per-file stage times across the three serial passes:

| Stage | Seconds / 50 files | Seconds / file |
|---|---:|---:|
| Calibration | 24.31 | 0.486 |
| MS1 feature extraction | 14.69 | 0.294 |
| MS2 association | 1.40 | 0.028 |
| MS2-seeded gap filling | 2.98 | 0.060 |
| Ion grouping | 4.35 | 0.087 |
| Table writing | 6.49 | 0.130 |

Calibration is the largest stage: approximately 44% of serial wall time. Stage medians
are calculated independently and omit small orchestration overhead; they need
not sum exactly to the full workflow median.

## Scope and reproducibility

- Warm JIT and filesystem cache. Downloads, first-run compilation, independent
  raw-data audits and test execution are excluded. No other benchmark or tests
  ran concurrently. Routine documentation work is not a dedicated-machine guarantee.
- Four-worker timings include fresh worker-pool startup for each 50-file batch.
  Long-running workers may amortize that overhead differently.
- Two-million-file figures extrapolate the measured cohort at continuous operation.
  Profile data, larger files, different chromatography, slower storage and network
  transfer can change throughput substantially. These are capacity estimates, not guarantees.
- All repeated outputs are identical within a worker count, and one/four-worker
  tables are identical. The fingerprint and full configuration are saved in the JSON reports.
- Output tables total 30.63 MB, plus 10.31 MB of JSON sidecars.
- Serial peak process RSS was 600 MiB, including warmup and loaded libraries.

[Serial measurements](extraction-v0.15.0-serial.json) · [Four-worker measurements](extraction-v0.15.0-four-workers.json) · [Validation](VALIDATION.md)

```bash
python benchmarks/benchmark_extraction.py --format parquet --workers 1 --repeats 3 \
  --output-directory validation/speed/serial --report validation/speed/serial.json
python benchmarks/benchmark_extraction.py --format parquet --workers 4 --repeats 3 \
  --output-directory validation/speed/four --report validation/speed/four.json
```

Raw inputs and the download manifest remain in `validation/MSV000086989/` locally.
Repeated generated tables were removed after validation; rerun these commands to regenerate them.
