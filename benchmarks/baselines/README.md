# Before-round source snapshot

`masscubex-before-overnight-2026-11-04.tar.gz` retains the exact extraction source
used before the overnight optimization, without generated Python/Numba caches.
It is a benchmark reference, not a release artifact. License notices are included.

The extraction implementation fingerprint on the recorded dependencies is
`da4d382e1609a392fa09765a31c96c851e1de21ccab6eb8defa7281d81167a68`.
The archive SHA-256 is
`f879310dabf467ff21c61b51967c248f5c0dc5706f5df220e8289cfa3df67451`.

Extract into an empty directory and pass that directory to the benchmark's
`--package-path` option. First generate a baseline report on the local inputs;
then run the current implementation with `--reference` pointing to that report.
Raw checksums are in `results/overnight-input-manifest-100.json` relative to
`benchmarks/`. Scientific-sidecar hashes include local input paths/mtimes and
must be regenerated for another checkout. The source snapshot and recorded
dependency versions let the two implementations be compared on the same inputs.
