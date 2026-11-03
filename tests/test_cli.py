import csv
import json
from pathlib import Path

import pytest

from masscubex import Config
from masscubex.cli import main


def test_default_quantification_and_numbered_catalog(capsys):
    assert main(["defaults"]) == 0
    config = json.loads(capsys.readouterr().out)
    assert config["quant_method"] == "peak_height"
    assert config["best_ms2"] == "fragment_intensity"
    assert config["associate_ms2"] is True
    assert config["ms2_max_isolation_width"] == 0
    assert config["ms2_seeded_features"] is True
    assert config["ms2_gapfill_rt_window"] == 5
    assert config["noise_estimator"] == "low_tail"
    assert config["noise_fraction"] == 0.10
    assert config["noise_sigma"] == 3
    assert config["noise_window_da"] == 50
    assert config["noise_min_neighbors"] == 5
    assert Config.from_dict(config) == Config()
    assert main(["adducts", "--json"]) == 0
    rows = json.loads(capsys.readouterr().out)
    assert [r["id"] for r in rows if r["enabled"]] == list(range(1, 16)) + [22, 23]
    assert main(["adducts", "--adducts", "1,5"]) == 0
    lines = capsys.readouterr().out.splitlines()
    rows = list(csv.DictReader(lines, delimiter="\t"))
    assert [r["name"] for r in rows if r["enabled"] == "true"] == ["[M+H]+", "[M-H]-"]


def test_ms2_seeding_and_width_controls_override_saved_settings(tmp_path, capsys):
    path = tmp_path / "old.json"
    path.write_text(json.dumps({"ms2_max_isolation_width": 4}))
    assert (
        main(
            [
                "extract",
                "--config",
                str(path),
                "--ms2-max-isolation-width",
                "0",
                "--no-ms2-seeded-features",
                "--ms2-gapfill-rt-window",
                "8",
                "--print-config",
            ]
        )
        == 0
    )
    config = json.loads(capsys.readouterr().out)
    assert config["ms2_max_isolation_width"] == 0
    assert config["ms2_seeded_features"] is False
    assert config["ms2_gapfill_rt_window"] == 8


def test_json_configuration_explicit_flags_override_and_custom_ion(tmp_path, capsys):
    path = tmp_path / "settings.json"
    path.write_text(
        json.dumps(
            {"quant_method": "peak_area", "noise_sigma": 5, "segment": False, "adduct_ids": [1, 5]}
        )
    )
    custom = tmp_path / "ions.json"
    custom.write_text(
        json.dumps(
            [{"id": 100, "name": "[M+Li]+", "charge": 1, "molecules": 1, "mass_shift": 7.015454857}]
        )
    )
    assert (
        main(
            [
                "extract",
                "--config",
                str(path),
                "--adduct-file",
                str(custom),
                "--adducts",
                "1,100",
                "--noise-sigma",
                "3",
                "--noise-estimator",
                "low_tail",
                "--noise-fraction",
                "0.05",
                "--noise-min-neighbors",
                "12",
                "--group-ions",
                "--no-isf",
                "--no-ms2",
                "--print-config",
            ]
        )
        == 0
    )
    values = json.loads(capsys.readouterr().out)
    assert values["quant_method"] == "peak_area"
    assert values["noise_sigma"] == 3 and values["segment"] is False
    assert values["noise_estimator"] == "low_tail"
    assert values["noise_fraction"] == 0.05 and values["noise_min_neighbors"] == 12
    assert values["adduct_ids"] == [1, 100]
    assert values["custom_adducts"][0]["name"] == "[M+Li]+"
    assert values["group_ions"] is True and values["group_isf"] is False
    assert values["associate_ms2"] is False
    assert Config.from_dict(values).to_dict() == values


@pytest.mark.parametrize(
    "flags",
    [
        ["--quant-method", "invalid"],
        ["--ms2-mz-da", "0", "--ms2-mz-ppm", "0"],
        ["--group-min-scans", "2"],
        ["--group-min-correlation", "1.1"],
        ["--adducts", "22"],
        ["--noise-estimator", "unknown"],
        ["--noise-fraction", "0"],
        ["--noise-fraction", "1.01"],
        ["--noise-sigma", "4"],
    ],
)
def test_cli_rejects_invalid_settings_before_processing(flags, capsys):
    with pytest.raises(SystemExit) as error:
        main(["extract", "--print-config", *flags])
    assert error.value.code == 2
    assert "error:" in capsys.readouterr().err


def test_cli_area_selection_reaches_output(sample_file, tmp_path, capsys):
    assert (
        main(
            [
                "extract",
                str(sample_file),
                "-o",
                str(tmp_path / "out"),
                "--quant-method",
                "peak_area",
                "--best-ms2",
                "apex",
                "--no-segment",
            ]
        )
        == 0
    )
    result = json.loads(capsys.readouterr().out)
    with Path(result["output"]).open() as stream:
        (row,) = csv.DictReader(stream, delimiter="\t")
    assert row["quant_method"] == "peak_area"
    assert row["intensity"] == row["peak_area"] == "780"
    assert row["peak_height"] == "100"
    metadata = json.loads(Path(result["output"] + ".json").read_text())
    assert metadata["config"]["best_ms2"] == "apex"
    assert metadata["config"]["segment"] is False


def test_cli_rejects_old_non_three_sigma_json_configuration(tmp_path, capsys):
    path = tmp_path / "old_settings.json"
    path.write_text(json.dumps({"noise_sigma": 5}))
    with pytest.raises(SystemExit) as error:
        main(["extract", "--config", str(path), "--print-config"])
    assert error.value.code == 2
    assert "noise_sigma is fixed at 3" in capsys.readouterr().err
