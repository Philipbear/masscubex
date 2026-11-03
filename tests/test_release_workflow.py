"""Whole-file equivalence, isolation and restart checks for the tuned workflow."""

import json
from pathlib import Path

import pytest

from masscubex import Config, process_file
from masscubex.cli import run_batch
from conftest import mzml_spectrum, write_mzml
from test_ms2 import precursor_mzml
from test_ms2_seeding import outputs
from test_polarity_switching import switching_file


@pytest.mark.parametrize("executor", ["process", "thread"])
def test_supported_tuning_is_identical_in_parallel_and_resumable(tmp_path, executor):
    source = switching_file(tmp_path, profile=True, weak_seeds=True)
    config = Config(
        noise_mode="none", min_intensity=100, tuning_block_scans=2048, tuning_block_seconds=300
    )
    reference = tmp_path / "reference.tsv"
    _, meta = outputs(source, reference, config)
    paths = [source]
    for i in range(3):
        copy = tmp_path / f"copy-{i}.mzML.gz"
        copy.write_bytes(source.read_bytes())
        paths.append(copy)
    results = list(run_batch(paths, tmp_path / "out", config, workers=2, executor=executor))
    assert len(results) == len(paths) and all(r["status"] == "ok" for r in results)
    for row in results:
        assert Path(row["output"]).read_bytes() == reference.read_bytes()
        saved = json.loads(Path(row["output"] + ".json").read_text())
        assert (
            saved["parameter_tuning"]["effective_configs"]
            == meta["parameter_tuning"]["effective_configs"]
        )
    resumed = list(
        run_batch(paths, tmp_path / "out", config, workers=2, executor=executor, resume=True)
    )
    assert all(r["status"] == "skipped" for r in resumed)


@pytest.mark.parametrize("auto", [False, True])
def test_batched_polarity_switching_does_not_bridge_acquisition_gaps(tmp_path, auto):
    events = []
    for times, polarity, height in [
        ([0, 1, 2, 10, 11, 12], "positive", 10),
        ([5, 6, 7, 15, 16, 17], "negative", 100),
    ]:
        for t in times:
            events.append((t, mzml_spectrum(t, t, [100], [height], polarity=polarity)))
        for t in [times[1], times[4]]:
            events.append(
                (
                    t + 0.1,
                    mzml_spectrum(
                        100 + t,
                        t + 0.1,
                        [50],
                        [1],
                        level=2,
                        polarity=polarity,
                        precursor=precursor_mzml(ref=f"scan={t}"),
                    ),
                )
            )
    path = write_mzml(tmp_path / "blocks.mzML", [s for _, s in sorted(events)])
    rows, meta = outputs(
        path,
        tmp_path / "out.tsv",
        Config(auto_tune=auto, noise_mode="none", min_scans=3, max_gap_seconds=4),
    )
    assert len(rows) == 4 and meta["ms2_assigned"] == 4
    assert all(r["feature_origin"] == "ms1_detected" and r["ms2_count"] == "1" for r in rows)
    for row in rows:
        expected = 10 if row["polarity"] == "positive" else 100
        assert float(row["peak_height"]) == expected
        assert float(row["peak_area"]) == 2 * expected


def test_late_parse_error_preserves_previous_outputs_and_removes_partials(tmp_path):
    source = switching_file(tmp_path)
    config = Config(noise_mode="none", tuning_block_scans=2048, tuning_block_seconds=300)
    out = tmp_path / "out.tsv"
    process_file(source, out, config=config)
    metadata = Path(str(out) + ".json")
    old = out.read_bytes(), metadata.read_bytes()
    # Valid content followed by an invalid gzip trailer: fail after processing content.
    source.write_bytes(source.read_bytes()[:-8])
    with pytest.raises((EOFError, OSError)):
        process_file(source, out, config=config, overwrite=True)
    assert (out.read_bytes(), metadata.read_bytes()) == old
    assert not list(tmp_path.glob("*.lock")) and not list(tmp_path.glob("*.partial"))


def test_no_auto_tune_keeps_explicit_fixed_parameters(tmp_path):
    source = switching_file(tmp_path)
    rows, meta = outputs(source, tmp_path / "out.tsv", Config(auto_tune=False, noise_mode="none"))
    assert not meta["parameter_tuning"]["enabled"]
    assert "effective_configs" not in meta["parameter_tuning"]
    assert rows and meta["ms2_assigned"] == 48
