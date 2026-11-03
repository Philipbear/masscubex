import numpy as np
import pytest

from masscubex.config import Config
from masscubex.features import TraceTracker, summarize_trace
from masscubex.pipeline import extract_features
from masscubex.readers import Scan


def scan(rt, intensities, mz=None, polarity="positive"):
    intensities = np.array(intensities, dtype=float)
    mz = np.array(mz if mz is not None else [100.0] * len(intensities))
    keep = intensities > 0
    return Scan(rt, mz[keep], intensities[keep], polarity)


def test_apex_mass_and_area_use_original_samples_and_seconds():
    config = Config(min_scans=3, noise_mode="none", segment=False)
    scans = [scan(t, [v], [100 + d]) for t, v, d in [(0, 10, 0), (2, 20, 0.001), (5, 10, 0)]]
    (feature,) = extract_features(scans, config)
    assert feature.mz == pytest.approx(100.0005)
    assert feature.rt_apex_seconds == 2
    assert feature.peak_height == 20
    assert feature.peak_area == 75
    assert feature.scan_count == 3


def test_missing_scans_contribute_zero_area_and_preserve_empty_scan_gaps():
    config = Config(min_scans=2, max_gap_scans=1, noise_mode="none", segment=False)
    (feature,) = extract_features([scan(0, [10]), scan(1, []), scan(2, [10])], config)
    assert feature.scan_count == 2
    assert feature.peak_area == 10
    separate = [scan(i, [10] if i in [0, 1, 4, 5] else []) for i in range(6)]
    assert len(list(extract_features(separate, config))) == 2


def test_long_time_gap_and_polarity_switching_do_not_join_features():
    config = Config(min_scans=2, noise_mode="none", segment=False, max_gap_seconds=3)
    scans = [scan(0, [10]), scan(1, [10]), scan(10, [10]), scan(11, [10])]
    assert len(list(extract_features(scans, config))) == 2
    scans = [scan(i, [10], polarity=p) for i in range(3) for p in ["positive", "negative"]]
    features = list(extract_features(scans, config))
    assert {f.polarity for f in features} == {"positive", "negative"}
    assert all(f.scan_count == 3 for f in features)


def test_equal_rt_is_an_explicit_error():
    with pytest.raises(ValueError, match="retention times"):
        list(extract_features([scan(1, [10]), scan(1, [10])], Config()))


def test_competing_peaks_use_second_available_match():
    config = Config(min_scans=2, noise_mode="none", segment=False, mz_tolerance_da=0.01)
    scans = [scan(0, [100, 90], [100, 100.004]), scan(1, [200, 150], [100.001, 100.008])]
    features = list(extract_features(scans, config))
    assert len(features) == 2
    assert sorted(f.scan_count for f in features) == [2, 2]


def test_ppm_tolerance_and_da_floor():
    config = Config(
        min_scans=2, noise_mode="none", segment=False, mz_tolerance_ppm=10, mz_tolerance_da=0.001
    )
    scans = [scan(0, [10, 20], [100, 1000]), scan(1, [10, 20], [100.002, 1000.008])]
    (feature,) = extract_features(scans, config)
    assert feature.mz == pytest.approx(1000.004)


def test_segment_two_peaks_conserves_area():
    config = Config(min_scans=3, noise_mode="none", peak_prominence=0.1)
    y = np.array([1, 5, 15, 30, 15, 5, 2, 5, 15, 40, 15, 5, 1.0])
    trace = np.column_stack((np.full(len(y), 100.0), y, np.arange(len(y)), np.zeros(len(y))))
    features = list(summarize_trace(trace, config))
    assert len(features) == 2
    assert [f.rt_apex_seconds for f in features] == [3, 9]
    assert sum(f.peak_area for f in features) == np.sum((y[1:] + y[:-1]) / 2)


def test_array_pool_growth_and_recycling():
    config = Config(min_scans=2, max_gap_scans=0, noise_mode="none")
    tracker = TraceTracker(config, capacity=2, node_capacity=2)
    mz = np.linspace(100, 200, 50)
    intensity = np.ones(50)
    threshold = np.zeros(50)
    finished = []
    for cycle in range(20):
        for t in range(3):
            finished.extend(tracker.add(mz, intensity, threshold, cycle * 4 + t))
        finished.extend(tracker.add(np.empty(0), np.empty(0), np.empty(0), cycle * 4 + 3))
    finished.extend(tracker.flush())
    assert len(finished) == 20 * 50
    assert all(len(trace) == 3 for trace in finished)
    assert tracker.available_tracks == len(tracker.track_values)
    assert tracker.available_nodes == len(tracker.nodes)
    assert len(set(tracker.free_tracks[: tracker.available_tracks])) == tracker.available_tracks
    assert len(set(tracker.free_nodes[: tracker.available_nodes])) == tracker.available_nodes
    assert len(tracker.nodes) < 400  # reuse slots instead of retaining previous cycles
