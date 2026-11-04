from dataclasses import replace
from pathlib import Path
import csv
import json
import subprocess
import sys

import pytest

from masscubex import Config, process_file
from masscubex.cli import iter_inputs, output_for, run_batch
from conftest import mzml_spectrum, mzxml_scan, write_mzml, write_mzxml


def test_tsv_and_provenance(sample_file, tmp_path):
    output = tmp_path / "features.tsv"
    result = process_file(sample_file, output)
    assert result.features == 1
    assert result.ms1_scans == 9
    with output.open() as stream:
        (row,) = csv.DictReader(stream, delimiter="\t")
    assert row["mz"] == "100.12346"
    assert float(row["peak_area"]) == 780
    metadata = json.loads(Path(str(output) + ".json").read_text())
    assert metadata["config"]["noise_mode"] == "adaptive"
    assert metadata["retained_centroids"] == 9
    assert len(metadata["tsv_sha256"]) == 64
    assert not list(tmp_path.glob("*.lock"))


def test_resume_validates_config_input_and_table(sample_file, tmp_path):
    output = tmp_path / "features.tsv"
    process_file(sample_file, output)
    before = output.stat().st_mtime_ns
    assert process_file(sample_file, output, resume=True).status == "skipped"
    assert output.stat().st_mtime_ns == before
    output.write_text("truncated")
    assert process_file(sample_file, output, resume=True).status == "ok"
    config = replace(Config(), min_scans=20)
    assert process_file(sample_file, output, resume=True, config=config).features == 0
    metadata = Path(str(output) + ".json")
    metadata.write_text("broken json")
    assert process_file(sample_file, output, resume=True).status == "ok"
    provenance = json.loads(metadata.read_text())
    provenance["schema_version"] = 9  # Previous export shape must be regenerated.
    metadata.write_text(json.dumps(provenance))
    assert process_file(sample_file, output, resume=True).status == "ok"
    provenance = json.loads(metadata.read_text())
    provenance["implementation_sha256"] = "old implementation"
    metadata.write_text(json.dumps(provenance))
    assert process_file(sample_file, output, resume=True).status == "ok"
    sample_file.write_text(sample_file.read_text() + "\n")
    assert process_file(sample_file, output, resume=True).status == "ok"


def test_existing_output_and_locks_protected(sample_file, tmp_path):
    output = tmp_path / "features.tsv"
    output.write_text("user data")
    with pytest.raises(FileExistsError):
        process_file(sample_file, output)
    assert output.read_text() == "user data"
    lock = Path(str(output) + ".lock")
    lock.write_text("a different worker")
    with pytest.raises(FileExistsError, match="locked"):
        process_file(sample_file, output, overwrite=True)
    assert lock.read_text() == "a different worker"


def test_binned_noise_source_change_invalidates_resume(sample_file, tmp_path, monkeypatch):
    from masscubex.pipeline import implementation_id

    output = tmp_path / "noise.parquet"
    process_file(sample_file, output)
    read_bytes = Path.read_bytes

    def changed_source(path):
        original = read_bytes(path)
        return (
            original + b"\n# edited bin estimator\n"
            if path.name == "noise_regions.py"
            else original
        )

    try:
        with monkeypatch.context() as patch:
            patch.setattr(Path, "read_bytes", changed_source)
            implementation_id.cache_clear()
            assert process_file(sample_file, output, resume=True).status == "ok"
    finally:
        implementation_id.cache_clear()


def test_failed_parse_leaves_no_outputs_and_preserves_old_result(sample_file, tmp_path):
    output = tmp_path / "features.tsv"
    process_file(sample_file, output)
    before = output.read_bytes()
    sample_file.write_text(sample_file.read_text()[:-30])
    with pytest.raises(ValueError):
        process_file(sample_file, output, overwrite=True)
    assert output.read_bytes() == before
    assert not list(tmp_path.glob("*.lock"))
    assert not list(tmp_path.glob("*.partial"))


def test_empty_scans_header_only_but_ms2_only_errors(tmp_path):
    path = write_mzml(tmp_path / "empty.mzML", [mzml_spectrum(1, 0, [], [])])
    result = process_file(path, output_format="tsv")
    assert result.features == 0
    assert len(Path(result.output).read_text().splitlines()) == 1
    path = write_mzml(tmp_path / "ms2.mzML", [mzml_spectrum(1, 0, [], [], level=2)])
    with pytest.raises(ValueError, match="No MS1"):
        process_file(path)
    assert not Path(str(path) + ".features.tsv").exists()


def test_mzml_mzxml_outputs_match(sample_file, tmp_path):
    path = write_mzxml(
        tmp_path / "other.mzXML",
        [
            mzxml_scan(i, i * 2, [100.123456789], [v])
            for i, v in enumerate([10, 20, 40, 80, 100, 80, 40, 20, 10])
        ],
    )
    a = process_file(sample_file)
    b = process_file(path)
    assert Path(a.output).read_bytes() == Path(b.output).read_bytes()


@pytest.mark.parametrize("executor", ["thread", "process"])
def test_parallel_cli_manifest_and_resume(sample_file, tmp_path, executor):
    other = tmp_path / "copy.mzML"
    other.write_bytes(sample_file.read_bytes())
    manifest = tmp_path / "paths.txt"
    manifest.write_text("sample.mzML\ncopy.mzML\n")
    output = tmp_path / "out"
    cmd = [
        sys.executable,
        "-m",
        "masscubex",
        "extract",
        "--manifest",
        str(manifest),
        "-o",
        str(output),
        "--workers",
        "2",
        "--executor",
        executor,
    ]
    run = subprocess.run(cmd, capture_output=True, text=True)
    assert run.returncode == 0, run.stderr
    results = [json.loads(line) for line in run.stdout.splitlines()]
    assert len(results) == 2
    assert all(r["status"] == "ok" for r in results)
    assert Path(results[0]["output"]).read_bytes() == Path(results[1]["output"]).read_bytes()
    run = subprocess.run(cmd + ["--resume"], capture_output=True, text=True)
    assert run.returncode == 0, run.stderr
    assert all(json.loads(line)["status"] == "skipped" for line in run.stdout.splitlines())


def test_batch_input_is_bounded_and_name_collisions_avoided(sample_file, tmp_path):
    consumed = []

    def paths():
        for i in range(20):
            path = tmp_path / f"{i}.mzML"
            path.write_bytes(sample_file.read_bytes())
            consumed.append(i)
            yield path

    results = run_batch(paths(), tmp_path / "out", Config(), workers=2, executor="thread")
    assert next(results)["status"] == "ok"
    assert len(consumed) == 4
    results.close()
    assert output_for(tmp_path / "a/run.mzML", "out") != output_for(tmp_path / "b/run.mzML", "out")


def test_cli_continues_after_file_error(sample_file, tmp_path):
    run = subprocess.run(
        [
            sys.executable,
            "-m",
            "masscubex",
            "extract",
            str(tmp_path / "missing.mzML"),
            str(sample_file),
            "-o",
            str(tmp_path / "out"),
        ],
        capture_output=True,
        text=True,
    )
    assert run.returncode == 1
    assert [json.loads(line)["status"] for line in run.stdout.splitlines()] == ["error", "ok"]


def test_directory_discovery_and_relative_manifests(sample_file, tmp_path):
    nested = tmp_path / "nested"
    nested.mkdir()
    other = nested / "another.MZXML.GZ"
    other.write_bytes(b"fixture")
    (nested / "ignored.tsv").write_text("not an input")
    assert set(iter_inputs([tmp_path])) == {sample_file, other}
    manifest = tmp_path / "manifest.txt"
    manifest.write_text("# input files\n\nsample.mzML\nnested/another.MZXML.GZ\n")
    assert set(iter_inputs([], manifest)) == {sample_file, other}
