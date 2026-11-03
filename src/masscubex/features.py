"""Compiled online ion traces backed by recyclable array pools.

Only active traces are retained. Finished traces are segmented and emitted;
there is no Python object or Python append call per centroid.
"""

from dataclasses import dataclass, field

import numpy as np
from numba import njit

from .config import Config


# Float state per live trace. The compiled kernels use these named columns
# instead of Python objects for each of the millions of observations.
MEAN_MZ, TOTAL_INTENSITY, LAST_INTENSITY, LAST_OBSERVED_TIME = range(4)
# Integer state: linked-list endpoints and observation/gap counters.
HEAD_NODE, TAIL_NODE, NODE_COUNT, OBSERVATION_COUNT, MISSING_SCANS = range(5)


@njit(cache=True, nogil=True)
def _append_observation(
    slot, mz, intensity, rt, threshold, track_layout, nodes, links, free_nodes, available
):
    """Append one observation to a trace, consuming one free pool node."""
    available -= 1
    node = free_nodes[available]
    nodes[node, 0] = mz
    nodes[node, 1] = intensity
    nodes[node, 2] = rt
    nodes[node, 3] = threshold
    links[node] = -1
    if track_layout[slot, NODE_COUNT] == 0:
        track_layout[slot, HEAD_NODE] = node
    else:
        links[track_layout[slot, TAIL_NODE]] = node
    track_layout[slot, TAIL_NODE] = node
    track_layout[slot, NODE_COUNT] += 1
    return available


@njit(cache=True, nogil=True)
def _advance_tracks(
    mz,
    intensity,
    threshold,
    rt,
    ppm,
    da,
    max_gaps,
    max_seconds,
    active,
    active_count,
    track_values,
    track_layout,
    nodes,
    links,
    free_tracks,
    available_tracks,
    free_nodes,
    available_nodes,
):
    previous_active = active[:active_count].copy()
    # Stronger traces claim peaks first; stable ties preserve acquisition order.
    priority = np.argsort(-track_values[previous_active, LAST_INTENSITY], kind="mergesort")
    used = np.zeros(mz.size, dtype=np.bool_)
    finished_tracks = np.empty(active_count, dtype=np.int64)
    finished_count = 0
    active_count = 0
    for k in priority:
        slot = previous_active[k]
        if rt - track_values[slot, LAST_OBSERVED_TIME] > max_seconds:
            finished_tracks[finished_count] = slot
            finished_count += 1
            continue
        target = track_values[slot, MEAN_MZ]
        tolerance = max(da, target * ppm * 1e-6)
        start = np.searchsorted(mz, target - tolerance)
        best = -1
        error = np.inf
        for j in range(start, mz.size):
            if mz[j] > target + tolerance:
                break
            distance = abs(mz[j] - target)
            if not used[j] and distance < error:
                best, error = j, distance
        if best >= 0:
            used[best] = True
            weight = track_values[slot, TOTAL_INTENSITY] + intensity[best]
            track_values[slot, MEAN_MZ] += (mz[best] - target) * (intensity[best] / weight)
            track_values[slot, TOTAL_INTENSITY] = weight
            track_values[slot, LAST_INTENSITY] = intensity[best]
            track_values[slot, LAST_OBSERVED_TIME] = rt
            track_layout[slot, OBSERVATION_COUNT] += 1
            track_layout[slot, MISSING_SCANS] = 0
            available_nodes = _append_observation(
                slot,
                mz[best],
                intensity[best],
                rt,
                threshold[best],
                track_layout,
                nodes,
                links,
                free_nodes,
                available_nodes,
            )
        else:
            track_layout[slot, MISSING_SCANS] += 1
            if track_layout[slot, MISSING_SCANS] > max_gaps:
                finished_tracks[finished_count] = slot
                finished_count += 1
                continue
            available_nodes = _append_observation(
                slot, target, 0.0, rt, 0.0, track_layout, nodes, links, free_nodes, available_nodes
            )
        active[active_count] = slot
        active_count += 1
    for j in range(mz.size):
        if used[j]:
            continue
        available_tracks -= 1
        slot = free_tracks[available_tracks]
        track_layout[slot, :] = 0
        track_layout[slot, OBSERVATION_COUNT] = 1
        track_values[slot, MEAN_MZ] = mz[j]
        track_values[slot, TOTAL_INTENSITY] = intensity[j]
        track_values[slot, LAST_INTENSITY] = intensity[j]
        track_values[slot, LAST_OBSERVED_TIME] = rt
        available_nodes = _append_observation(
            slot,
            mz[j],
            intensity[j],
            rt,
            threshold[j],
            track_layout,
            nodes,
            links,
            free_nodes,
            available_nodes,
        )
        active[active_count] = slot
        active_count += 1
    return active_count, available_tracks, available_nodes, finished_tracks[:finished_count]


@njit(cache=True, nogil=True)
def _drain_tracks(
    finished_tracks,
    minimum,
    track_layout,
    nodes,
    links,
    free_tracks,
    available_tracks,
    free_nodes,
    available_nodes,
):
    """Materialize qualified completed traces and recycle their pool storage."""
    traces = []
    for slot in finished_tracks:
        keep = track_layout[slot, OBSERVATION_COUNT] >= minimum
        trace = np.empty((track_layout[slot, NODE_COUNT] if keep else 0, 4))
        node = track_layout[slot, HEAD_NODE]
        for j in range(track_layout[slot, NODE_COUNT]):
            if keep:
                trace[j, :] = nodes[node, :]
            free_nodes[available_nodes] = node
            available_nodes += 1
            node = links[node]
        free_tracks[available_tracks] = slot
        available_tracks += 1
        if keep:
            traces.append(trace)
    return traces, available_tracks, available_nodes


class TraceTracker:
    """State for one file and polarity; never shared between workers."""

    def __init__(self, config: Config, capacity=1024, node_capacity=8192):
        self.config = config
        self.track_values = np.zeros((capacity, 4))
        self.track_layout = np.zeros((capacity, 5), dtype=np.int64)
        self.nodes = np.empty((node_capacity, 4))
        self.links = np.empty(node_capacity, dtype=np.int64)
        self.free_tracks = np.arange(capacity, dtype=np.int64)
        self.free_nodes = np.arange(node_capacity, dtype=np.int64)
        self.active = np.empty(capacity, dtype=np.int64)
        self.available_tracks, self.available_nodes, self.active_count = capacity, node_capacity, 0
        self.last_rt = -np.inf
        self.peak_active_traces = 0
        self.peak_live_nodes = 0

    def _ensure(self, peaks):
        if self.available_tracks < peaks:
            old = len(self.track_values)
            size = max(old * 2, old + peaks - self.available_tracks)
            self.track_values = np.vstack((self.track_values, np.zeros((size - old, 4))))
            self.track_layout = np.vstack(
                (self.track_layout, np.zeros((size - old, 5), dtype=np.int64))
            )
            free = np.empty(size, dtype=np.int64)
            free[: self.available_tracks] = self.free_tracks[: self.available_tracks]
            free[self.available_tracks : self.available_tracks + size - old] = np.arange(old, size)
            self.free_tracks = free
            self.available_tracks += size - old
            active = np.empty(size, dtype=np.int64)
            active[: self.active_count] = self.active[: self.active_count]
            self.active = active
        needed = self.active_count + peaks
        if self.available_nodes < needed:
            old = len(self.nodes)
            size = max(old * 2, old + needed - self.available_nodes)
            self.nodes = np.vstack((self.nodes, np.empty((size - old, 4))))
            self.links = np.concatenate((self.links, np.empty(size - old, dtype=np.int64)))
            free = np.empty(size, dtype=np.int64)
            free[: self.available_nodes] = self.free_nodes[: self.available_nodes]
            free[self.available_nodes : self.available_nodes + size - old] = np.arange(old, size)
            self.free_nodes = free
            self.available_nodes += size - old

    def _finish(self, finished_tracks):
        traces, self.available_tracks, self.available_nodes = _drain_tracks(
            finished_tracks,
            self.config.min_scans,
            self.track_layout,
            self.nodes,
            self.links,
            self.free_tracks,
            self.available_tracks,
            self.free_nodes,
            self.available_nodes,
        )
        return traces

    def add(self, mz, intensity, threshold, rt):
        if rt <= self.last_rt:
            raise ValueError(
                "MS1 retention times must increase strictly within each polarity; "
                "sort scans or collapse ion-mobility frames before extraction"
            )
        self.last_rt = rt
        self._ensure(len(mz))
        self.active_count, self.available_tracks, self.available_nodes, finished_tracks = (
            _advance_tracks(
                mz,
                intensity,
                threshold,
                float(rt),
                float(self.config.mz_tolerance_ppm),
                float(self.config.mz_tolerance_da),
                self.config.max_gap_scans,
                float(self.config.max_gap_seconds),
                self.active,
                self.active_count,
                self.track_values,
                self.track_layout,
                self.nodes,
                self.links,
                self.free_tracks,
                self.available_tracks,
                self.free_nodes,
                self.available_nodes,
            )
        )
        self.peak_active_traces = max(self.peak_active_traces, self.active_count)
        self.peak_live_nodes = max(self.peak_live_nodes, len(self.nodes) - self.available_nodes)
        return self._finish(finished_tracks)

    def flush(self):
        finished_tracks = self.active[: self.active_count].copy()
        self.active_count = 0
        return self._finish(finished_tracks)


@dataclass(slots=True)
class Feature:
    mz: float
    rt_apex_seconds: float
    rt_start_seconds: float
    rt_end_seconds: float
    peak_height: float | None
    peak_area: float | None
    scan_count: int
    apex_noise_threshold: float | None
    apex_threshold_ratio: float | None
    polarity: str
    intensity: float | None = 0.0
    quant_method: str = "peak_height"
    ms2_count: int = 0
    ms2_scan_numbers: tuple = ()
    ms2_native_ids: tuple = ()
    ms2_spectrum_indices: tuple = ()
    best_ms2_scan_number: int | None = None
    best_ms2_native_id: str | None = None
    best_ms2_spectrum_index: int | None = None
    best_ms2_selection_method: str | None = None
    best_ms2_score: float | None = None
    best_ms2_fragment_intensity: float | None = None
    best_ms2_precursor_intensity: float | None = None
    best_ms2_precursor_intensity_source: str | None = None
    ion_relations: tuple = ()
    feature_origin: str = "ms1_detected"
    ms1_quant_status: str = "detected"
    seed_precursor_mz: float | None = None
    seed_rt_seconds: float | None = None
    gapfill_status: str | None = None
    gapfill_mz_radius_da: float | None = None
    trace: np.ndarray | None = field(default=None, repr=False, metadata={"tsv": False})
    _ms2_records: tuple = field(default=(), repr=False, metadata={"tsv": False})


def summarize_trace(trace, config: Config, polarity="unknown"):
    """Segment one trace in compiled code; measure the original observations."""
    from .chromatograms import summarize_chromatogram

    summaries, tied = summarize_chromatogram(
        trace,
        config.min_scans,
        config.segment,
        float(config.smoothing_sigma),
        float(config.smoothing_seconds),
        config.peak_distance,
        float(config.peak_distance_seconds),
        float(config.peak_prominence),
    )
    if tied:
        from .reference_chromatograms import summarize_reference_trace

        yield from summarize_reference_trace(trace, config, polarity)
        return
    for peak in summaries:
        yield Feature(
            mz=peak.mz,
            rt_apex_seconds=peak.apex_time,
            rt_start_seconds=peak.start_time,
            rt_end_seconds=peak.end_time,
            peak_height=peak.height,
            peak_area=peak.area,
            scan_count=peak.observations,
            apex_noise_threshold=peak.noise_threshold,
            apex_threshold_ratio=peak.height / peak.noise_threshold
            if peak.noise_threshold > 0
            else None,
            polarity=polarity,
            intensity=peak.height if config.quant_method == "peak_height" else peak.area,
            quant_method=config.quant_method,
            trace=trace[peak.first : peak.last + 1, 1:3].copy() if config.group_ions else None,
        )


def warmup():
    """Compile/load the small fixed kernel set once per long-lived worker."""
    from .noise import estimate_noise
    from .centroid import centroid_profile
    from .grouping import trace_correlation, _candidate_pairs, sorted_mass_windows
    from .gapfill import eic_intensities
    from .calibration import measure_representative
    from .gapfill_refinement import _eic

    config = Config(min_scans=2)
    mz = np.array([100.0, 200.0])
    intensity = np.array([10.0, 20.0])
    threshold = estimate_noise(mz, intensity, config)[2]
    tracker = TraceTracker(config)
    tracker.add(mz, intensity, threshold, 0.0)
    tracker.add(mz, intensity, threshold, 1.0)
    tracker.flush()
    centroid_profile(mz, intensity, np.zeros(2), threshold, 3, 4.0)
    trace = np.column_stack((intensity, np.arange(2, dtype=float)))
    trace_correlation(trace, trace)
    eic_intensities(mz, intensity, mz, 20.0, 0.01)
    list(summarize_trace(np.column_stack((mz, intensity, mz, threshold)), config))
    measure_representative(np.zeros((2, 4)), 0.05, 100.0)
    ids = np.arange(2, dtype=np.int64)
    sorted_mass_windows(mz, mz, mz)
    _eic(mz, intensity, np.array([0, 2], dtype=np.int64), 100.0, 0.01)
    list(_candidate_pairs(ids, ids, ids, ids + 1, mz, intensity, ids + 1, 5.0))
