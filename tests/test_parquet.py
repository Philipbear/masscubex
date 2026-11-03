from dataclasses import replace
import json
from pathlib import Path
import struct

import pyarrow as pa
import pyarrow.parquet as pq
import pytest

from masscubex import Config, iter_feature_rows, process_file
from masscubex.cli import main
from masscubex.features import Feature
from masscubex.output import ROW_GROUP_SIZE, parquet_schema, write_features


def example():
    return Feature(
        mz=100.12345678912345,
        rt_apex_seconds=12.000000000000003,
        rt_start_seconds=10.0,
        rt_end_seconds=14.0,
        peak_height=None,
        peak_area=None,
        scan_count=0,
        apex_noise_threshold=None,
        apex_threshold_ratio=None,
        polarity="positive",
        feature_origin="ms2_seeded",
        ms1_quant_status="no_ms1_signal",
        ms2_count=2,
        ms2_scan_numbers=(None, 9223372036854775807),
        ms2_native_ids=("spectrum=α", "scan=9223372036854775807"),
        ms2_spectrum_indices=(0, 1),
        best_ms2_scan_number=9223372036854775807,
        ion_relations=(
            {
                "relation": "in_source_fragment",
                "status": "candidate",
                "partner_feature_id": 2,
                "parent_feature_id": 1,
                "fragment_feature_id": 2,
                "mass_error_ppm": -1.23456789012345,
                "correlation": 0.99,
                "shared_scans": 5,
                "role": "parent",
                "parent_ion": None,
                "fragment_ion": None,
                "ms2_supported": True,
                "evidence_ms2_scan_number": None,
                "evidence_ms2_native_id": "spectrum=α",
                "evidence_ms2_spectrum_index": 0,
                "fragment_mz": 80.0,
                "fragment_relative_intensity": 0.5,
            },
        ),
    )


def test_parquet_preserves_exact_values_nested_lists_nulls_and_unicode(tmp_path):
    features = [example(), replace(example(), ion_relations=(), ms2_scan_numbers=())]
    a, b = tmp_path / "table.tsv", tmp_path / "table.parquet"
    assert write_features(a, features, "tsv") == write_features(b, features, "parquet") == 2
    assert list(iter_feature_rows(a)) == list(iter_feature_rows(b))
    table = pq.read_table(b)
    assert table.schema.field("mz").type == pa.float64()
    assert pa.types.is_list(table.schema.field("ms2_scan_numbers").type)
    assert pa.types.is_struct(table.schema.field("ion_relations").type.value_type)
    assert table["mz"][0].as_py() == features[0].mz
    assert table["peak_height"][0].as_py() is None
    assert table["ms2_scan_numbers"][0].as_py() == [None, 9223372036854775807]
    assert pq.read_metadata(b).row_group(0).column(0).compression == "ZSTD"


def test_empty_and_multiple_row_groups_have_identical_explicit_schema(tmp_path):
    empty, full = tmp_path / "empty.parquet", tmp_path / "full.parquet"
    write_features(empty, [], "parquet")
    write_features(full, (example() for _ in range(ROW_GROUP_SIZE + 1)), "parquet")
    assert pq.read_schema(empty) == pq.read_schema(full) == parquet_schema()
    assert list(iter_feature_rows(empty)) == []
    assert pq.read_metadata(full).num_row_groups == 2
    ids = [r["feature_id"] for r in iter_feature_rows(full)]
    assert ids == list(range(1, ROW_GROUP_SIZE + 2))


def test_default_parquet_explicit_tsv_and_format_mismatch(sample_file, tmp_path):
    result = process_file(sample_file)
    assert result.output.endswith(".features.parquet")
    metadata = json.loads(Path(result.output + ".json").read_text())
    assert metadata["output_format"] == "parquet"
    assert metadata["output_bytes"] == Path(result.output).stat().st_size
    assert len(metadata["output_sha256"]) == 64
    tsv = tmp_path / "explicit.tsv"
    process_file(sample_file, tsv)
    assert list(iter_feature_rows(result.output)) == list(iter_feature_rows(tsv))
    with pytest.raises(ValueError, match="does not match"):
        process_file(sample_file, tmp_path / "wrong.tsv", output_format="parquet")
    with pytest.raises(ValueError, match="Output format"):
        process_file(sample_file, output_format="unknown")


def test_parquet_resume_checks_content_format_and_configuration(sample_file, tmp_path):
    out = tmp_path / "features.parquet"
    process_file(sample_file, out)
    assert process_file(sample_file, out, resume=True).status == "skipped"
    out.write_bytes(b"PAR1truncated")
    assert process_file(sample_file, out, resume=True).status == "ok"
    metadata_path = Path(str(out) + ".json")
    data = json.loads(metadata_path.read_text())
    data["output_format"] = "tsv"
    metadata_path.write_text(json.dumps(data))
    assert process_file(sample_file, out, resume=True).status == "ok"
    assert process_file(sample_file, out, resume=True, config=Config(min_scans=20)).features == 0
    assert pq.read_table(out).num_rows == 0


def test_failed_parquet_write_keeps_prior_table_and_completion_marker(
    sample_file, tmp_path, monkeypatch
):
    from masscubex import pipeline

    out = tmp_path / "features.parquet"
    process_file(sample_file, out)
    sidecar = Path(str(out) + ".json")
    before = out.read_bytes(), sidecar.read_bytes()

    def fail(path, features, output_format):
        Path(path).write_bytes(b"PAR1incomplete")
        raise OSError("injected output failure")

    monkeypatch.setattr(pipeline, "write_features", fail)
    with pytest.raises(OSError, match="injected"):
        process_file(sample_file, out, overwrite=True)
    assert (out.read_bytes(), sidecar.read_bytes()) == before
    assert not list(tmp_path.glob("*.lock")) and not list(tmp_path.glob("*.partial"))


@pytest.mark.parametrize("output_format", ["parquet", "tsv"])
def test_cli_format_and_nonfinite_validation(sample_file, tmp_path, capsys, output_format):
    assert (
        main(["extract", str(sample_file), "-o", str(tmp_path / "out"), "--format", output_format])
        == 0
    )
    result = json.loads(capsys.readouterr().out)
    assert result["output"].endswith("." + output_format)
    assert len(list(iter_feature_rows(result["output"]))) == 1
    with pytest.raises(ValueError, match="Non-finite"):
        write_features(
            tmp_path / ("bad." + output_format),
            [replace(example(), peak_height=float("inf"))],
            output_format,
        )


def test_no_ms2_no_grouping_streams_to_parquet(sample_file, tmp_path):
    config = Config(associate_ms2=False, group_ions=False)
    out = tmp_path / "stream.parquet"
    process_file(sample_file, out, config=config)
    assert len(list(iter_feature_rows(out))) == 1
    data = json.loads(Path(str(out) + ".json").read_text())
    assert data["timings_seconds"]["feature_extraction"] > 0
    assert data["timings_seconds"]["table_write"] >= 0


def test_unknown_evidence_cannot_be_silently_dropped(tmp_path):
    feature = replace(example(), ion_relations=({"relation": "adduct", "future_field": 1},))
    with pytest.raises(ValueError, match="Unknown ion evidence"):
        write_features(tmp_path / "unknown.parquet", [feature], "parquet")


def test_parquet_preserves_extreme_finite_float_bits(tmp_path):
    values = [-0.0, 5e-324, 1.7976931348623157e308]
    output = tmp_path / "precision.parquet"
    write_features(output, [replace(example(), peak_area=v) for v in values], "parquet")
    actual = pq.read_table(output, columns=["peak_area"])["peak_area"].to_pylist()
    assert [struct.pack("d", v) for v in actual] == [struct.pack("d", v) for v in values]
