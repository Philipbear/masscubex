# Upstream extraction comparison — 2026-11-05

The latest five-file comparison is **16.56× faster in aggregate** than pinned
MassCube 1.2.22 on the same adaptive-noise workload. Individual files are
12.75–23.74× faster, with zero mismatched scientific rows. The fixed Orbitrap
floor control is **5.03× faster**. These are different workloads; neither number
is a claim about stock MassCube with all its original defaults.

## Identical adaptive workload

Both engines receive the same adaptive MS1 filtering, strict 0.01 Da trace
radius, ten missing scans, two splitting passes, MS2 floor 500 and selected ion
catalog with grouping enabled. Auto-tuning, profile centroiding and MS2 seeding
are disabled only for this comparison because upstream lacks those extensions.
Normal MassCubeX defaults keep them enabled.

| File | Features in each engine | Original processing (s) | MassCubeX including Parquet (s) | Speedup |
| --- | ---: | ---: | ---: | ---: |
| NIST_POS_Samp_14-03.mzXML | 75,100 | 71.38 | 5.26 | 13.58× |
| NIST_POS_Samp_01-01.mzXML | 70,776 | 65.09 | 5.10 | 12.75× |
| QC_Mix_01.mzXML | 54,872 | 54.87 | 3.63 | 15.10× |
| Blank_00.mzXML | 55,614 | 58.66 | 3.79 | 15.50× |
| Blank_Post_14-03.mzXML | 114,323 | 146.99 | 6.19 | 23.74× |
| **Total** | **370,685** | **396.99** | **23.97** | **16.56×** |

The original timer includes decoding/preprocessing, trace detection, splitting,
finalization, MS2 association/selection and grouping, excluding serialization.
MassCubeX includes its complete Parquet transaction and provenance. Both are
single warm runs per file; warmup and independent parity checks are outside the
timers. This measurement design favors the original engine on serialization.

The reference is huaxuyu/masscube 1.2.22, commit
`356748ede691862b9a87f67812db5e728a703b9d`. Function-body hashes are verified
against the pinned Git objects before running. A Pyteomics adapter supplies
mzXML arrays and minute RTs because that upstream reader accepts mzML only.
The scientific function bodies are unchanged. Both engines share the adaptive
MS1 preprocessing; upstream does not provide that estimator. Single-file ion
grouping is explicitly enabled here and is normally disabled upstream.

## Fixed-floor control

The same five files, fixed MS1 floor 50,000 and MS2 floor 10,000 produce 76,243
features in each engine, with zero mismatched rows. Aggregate speedup is 5.03×.
The adaptive workload above has 370,685 features; its greater density changes
the relative cost of Python feature and grouping loops.

For both modes, parity checks compare core measurements, all associated MS2
references, best-scan choices and grouping assignments. Matching feature counts
alone are insufficient. These runs are positive centroid mzXML from
MSV000086989 on the same M2 Max Mac; they do not benchmark profile mass detection,
other instruments, network I/O or every public repository.

## Normal default workflow

Default MassCubeX uses calibration, profile support, MS2 seeding, grouping and
25-column Parquet output. The five files produce 448,891 features in total.
Their final median runtimes sum to **29.02 s**, versus **36.40 s** for the frozen
pre-overnight version. Each median uses two measured warm repeats; all five
Parquet tables and scientific sidecars match that baseline exactly, apart from
timings and implementation fingerprints. The calibrated defaults are not the
output-equivalent upstream configuration above.

The broader [100-file speed/RAM study](OVERNIGHT_OPTIMIZATION.md) covers the
normal default pipeline, sustained worker recycling and the recommended local
concurrency. See [methods](../../docs/METHODS.md) for the exact implementation
specializations and [the code audit](../../docs/CODE_AUDIT.md) for focused checks
and pre-existing test debt. This is development validation, with no PR, push,
release or full platform/wheel matrix.

[Adaptive comparison JSON](overnight-upstream-matched-2026-11-05.json) ·
[Fixed-floor JSON](overnight-upstream-orbitrap-2026-11-05.json) ·
[Default five-file medians](overnight-final-five-2026-11-05.json) ·
[Pre-round five-file baseline](overnight-before-five-2026-11-04.json)

The previous 2026-11-04 measurement (10.81× adaptive, 3.75× fixed floor) is
retained in [its original JSON](upstream-speed-five-files-2026-11-04.json).
Do not combine timings from the two runs.

Reproduce the latest comparison with `benchmarks/benchmark_upstream_speed.py`
and exactly five inputs. Use `--mode matched` for adaptive preprocessing or
`--mode orbitrap` for the fixed-floor control. The script requires Pyteomics and
the pinned upstream Git objects. `benchmark_exact_speed.py` accepts
`--package-path` and `--compare` for frozen-source table/sidecar checks. The
before-round source archive and hashes are in [baselines](../baselines/README.md).

Publication follow-up: the complete local source and installed-wheel suites each pass 637 tests.
Older expectations were reconciled against the pinned upstream contract without changing
the numerical implementation. See the [publication record](publication-validation-2026-11-05.json)
and [code audit](../../docs/CODE_AUDIT.md). Statements above about deferred publication
describe the original measurement round.
