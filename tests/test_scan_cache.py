"""Replay must preserve raw evidence, corruption detection and memory bounds."""

from dataclasses import replace

import numpy as np
import pytest

from masscubex import Config, process_file
from masscubex.output import iter_feature_rows
from masscubex.readers import MS2Scan
from masscubex.scan_cache import ScanSource
from conftest import mzml_spectrum, write_mzml, mzxml_scan, write_mzxml
from test_ms2 import precursor_mzml


def test_ms2_metadata_changes_do_not_leak_between_passes(tmp_path):
    path = write_mzml(
        tmp_path / "scan.mzML",
        [
            mzml_spectrum(1, 0, [100], [40]),
            mzml_spectrum(2, 0.1, [50, 60], [5, 10], level=2, precursor=precursor_mzml()),
        ],
    )
    source = ScanSource(path, Config())
    original = list(source.scans())[-1]
    original.polarity = "negative"
    original.precursor_intensity = 123
    original.mz = np.array([90.0])
    replay = list(source.scans())[-1]
    assert isinstance(replay, MS2Scan) and replay.polarity == "positive"
    assert replay.precursor_intensity != 123
    np.testing.assert_array_equal(replay.mz, [50, 60])
    assert len(list(source.ms1())) == 1
    assert source.stats["reader_passes"] == 1 and source.stats["replay_passes"] == 2


@pytest.mark.parametrize("budget", [0, 0.001, 64])
def test_cache_budget_and_partial_iteration(sample_file, budget):
    source = ScanSource(sample_file, Config(scan_cache_mb=budget))
    partial = source.scans()
    next(partial)
    partial.close()
    assert source._scans == [] and not source._complete
    scans = list(source.scans())
    assert len(scans) == 9 and not source._complete
    assert source.stats["peak_bytes"] <= source.budget


@pytest.mark.parametrize("suffix", ["mzML", "mzXML.gz"])
@pytest.mark.parametrize("profile", [False, True])
@pytest.mark.parametrize("auto_tune", [False, True])
def test_cached_overflow_and_streaming_tables_are_identical(tmp_path, suffix, profile, auto_tune):
    make, write = (mzml_spectrum, write_mzml) if suffix == "mzML" else (mzxml_scan, write_mzxml)
    spectra = []
    for index in range(30):
        mz = np.linspace(99.98, 100.02, 21) if profile else np.array([100.0])
        y = (20 + 100 * np.exp(-(((index - 14) / 5) ** 2))) * np.exp(-(((mz - 100) / 0.008) ** 2))
        mode = (
            {"mode": "profile" if profile else "centroid"}
            if suffix == "mzML"
            else {"centroid": "0" if profile else "1"}
        )
        spectra.append(make(index, float(index), mz, y, **mode))
    path = write(tmp_path / ("input." + suffix), spectra)
    base = Config(auto_tune=auto_tune)
    tables = []
    for budget in (0, 0.001, 64):
        output = tmp_path / f"out-{budget}.parquet"
        process_file(path, output, config=replace(base, scan_cache_mb=budget))
        tables.append(list(iter_feature_rows(output)))
    assert tables[0] == tables[1] == tables[2]
    assert tables[0]


def test_late_xml_corruption_cannot_complete_cache(tmp_path):
    path = write_mzml(tmp_path / "broken.mzML", [mzml_spectrum(1, 0, [100], [20])])
    path.write_text(path.read_text().replace("</mzML>", ""))
    source = ScanSource(path, Config())
    with pytest.raises(Exception):
        list(source.scans())
    assert not source._complete and not source._scans
