"""Compiled online ion traces backed by recyclable array pools.

Only active traces are retained. Finished traces are segmented and emitted;
there is no Python object or Python append call per centroid.
"""

from dataclasses import dataclass, field
import math

import numpy as np
from numba import njit

from .config import Config
from .mass_tolerance import ms1_radius, ms1_region_array
from .measurements import integrated_area
from .summary import _span_state, summarize_intervals


# Float state per live trace. The compiled kernels use these named columns
# instead of Python objects for each of the millions of observations.
LAST_MZ, TOTAL_INTENSITY, LAST_INTENSITY, LAST_OBSERVED_TIME = range(4)
# Integer state: linked-list endpoints and observation/gap counters.
HEAD_NODE, TAIL_NODE, NODE_COUNT, OBSERVATION_COUNT, MISSING_SCANS = range(5)


@njit(cache=True, inline="always")
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
    regions,
    max_gaps,
    max_seconds,
    first_scan,
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
    used = np.zeros(mz.size, dtype=np.bool_)
    finished_tracks = np.empty(active_count, dtype=np.int64)
    finished_count = 0
    active_count = 0
    # The initial scan remains in m/z order. Later scans use the previous
    # scan's stable intensity order, including zero intensity for missing ions.
    for slot in previous_active:
        if max_seconds > 0 and rt - track_values[slot, LAST_OBSERVED_TIME] > max_seconds:
            finished_tracks[finished_count] = slot
            finished_count += 1
            continue
        target = np.float32(track_values[slot, LAST_MZ])
        tolerance = ms1_radius(target, da, ppm, regions)
        best = -1
        if mz.size:
            index = np.searchsorted(mz, target)
            if index == 0:
                if mz[0] - target < tolerance:
                    best = 0
            elif index == mz.size:
                if target - mz[-1] < tolerance:
                    best = mz.size - 1
            else:
                before, after = target - mz[index - 1], mz[index] - target
                if after < before and after < tolerance:
                    best = index
                elif before < after and before < tolerance:
                    best = index - 1
                # Upstream leaves an exact midpoint unassigned.
        # Do not search for a second choice if another trace claimed the nearest.
        if best >= 0 and not used[best]:
            used[best] = True
            track_values[slot, LAST_MZ] = mz[best]
            track_values[slot, TOTAL_INTENSITY] += intensity[best]
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
            track_values[slot, LAST_INTENSITY] = 0
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
        track_values[slot, LAST_MZ] = mz[j]
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
    if not first_scan:
        ids = active[:active_count].copy()
        order = np.argsort(-track_values[ids, LAST_INTENSITY], kind="mergesort")
        active[:active_count] = ids[order]
    # Original removes expired traces from the active list in reverse order.
    return active_count, available_tracks, available_nodes, finished_tracks[:finished_count][::-1]


@njit(cache=True)
def _take_trace(
    slot,
    minimum,
    track_layout,
    nodes,
    links,
    free_tracks,
    available_tracks,
    free_nodes,
    available_nodes,
):
    """Copy one qualified trace, then release its nodes and track slot."""
    keep = track_layout[slot, OBSERVATION_COUNT] >= minimum
    trace = np.empty((track_layout[slot, NODE_COUNT] if keep else 0, 4))
    node = track_layout[slot, HEAD_NODE]
    for j in range(track_layout[slot, NODE_COUNT]):
        if keep:
            for column in range(4):
                trace[j, column] = nodes[node, column]
        free_nodes[available_nodes] = node
        available_nodes += 1
        node = links[node]
    free_tracks[available_tracks] = slot
    return trace, available_tracks + 1, available_nodes


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
        trace, available_tracks, available_nodes = _take_trace(
            slot,
            minimum,
            track_layout,
            nodes,
            links,
            free_tracks,
            available_tracks,
            free_nodes,
            available_nodes,
        )
        if len(trace):
            traces.append(trace)
    return traces, available_tracks, available_nodes


class TraceTracker:
    """State for one file and polarity; never shared between workers."""

    def __init__(self, config: Config, capacity=1024, node_capacity=8192):
        self.config = config
        self.mass_regions = ms1_region_array(config)
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

    def _advance(self, mz, intensity, threshold, rt):
        if rt <= self.last_rt:
            raise ValueError(
                "MS1 retention times must increase strictly within each polarity; "
                "sort scans or collapse ion-mobility frames before extraction"
            )
        first_scan = not np.isfinite(self.last_rt)
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
                self.mass_regions,
                self.config.max_gap_scans,
                float(self.config.max_gap_seconds),
                first_scan,
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
        return finished_tracks

    def add(self, mz, intensity, threshold, rt):
        """Return raw completed traces, used by the calibration pilot."""
        return self._finish(self._advance(mz, intensity, threshold, rt))

    def flush(self):
        finished_tracks = self.active[: self.active_count].copy()
        self.active_count = 0
        return self._finish(finished_tracks)

    def _finish_features(self, finished_tracks, polarity):
        for trace in self._finish(finished_tracks):
            yield from summarize_trace(trace, self.config, polarity)

    def add_features(self, mz, intensity, threshold, rt, polarity="unknown"):
        """Return features from the completed traces in this scan."""
        return self._finish_features(self._advance(mz, intensity, threshold, rt), polarity)

    def flush_features(self, polarity="unknown"):
        """Finish the remaining traces at EOF without a per-trace Python loop."""
        finished_tracks = self.active[: self.active_count].copy()
        self.active_count = 0
        return self._finish_features(finished_tracks, polarity)


@dataclass(slots=True)
class Feature:
    mz: float
    rt_apex_seconds: float
    rt_start_seconds: float
    rt_end_seconds: float
    peak_height: float | None
    peak_area: float | None
    scan_count: int = field(metadata={"tsv": False})
    apex_noise_threshold: float | None
    apex_threshold_ratio: float | None = field(metadata={"tsv": False})
    polarity: str
    # Retain the Python API's convenience values without repeating them in
    # every output row. Quantification settings live in the file sidecar.
    intensity: float | None = field(default=0.0, metadata={"tsv": False})
    quant_method: str = field(default="peak_height", metadata={"tsv": False})
    ms2_count: int = field(default=0, metadata={"tsv": False})
    ms2_scan_numbers: tuple = ()
    ms2_native_ids: tuple = field(default=(), metadata={"tsv": False})
    ms2_spectrum_indices: tuple = field(default=(), metadata={"tsv": False})
    best_ms2_scan_number: int | None = None
    best_ms2_native_id: str | None = field(default=None, metadata={"tsv": False})
    best_ms2_spectrum_index: int | None = field(default=None, metadata={"tsv": False})
    best_ms2_selection_method: str | None = field(default=None, metadata={"tsv": False})
    best_ms2_score: float | None = field(default=None, metadata={"tsv": False})
    best_ms2_fragment_intensity: float | None = field(default=None, metadata={"tsv": False})
    best_ms2_precursor_intensity: float | None = field(default=None, metadata={"tsv": False})
    best_ms2_precursor_intensity_source: str | None = field(default=None, metadata={"tsv": False})
    adduct: str | None = None
    group_id: int | None = None
    is_isotope: bool | None = field(default=None, metadata={"tsv": False})
    isotope_state: str | None = None
    is_in_source_fragment: bool | None = None
    ion_relations: tuple = ()
    feature_origin: str = "ms1_detected"
    ms1_quant_status: str = "detected"
    seed_precursor_mz: float | None = field(default=None, metadata={"tsv": False})
    seed_rt_seconds: float | None = field(default=None, metadata={"tsv": False})
    gapfill_status: str | None = field(default=None, metadata={"tsv": False})
    gapfill_mz_radius_da: float | None = field(default=None, metadata={"tsv": False})
    gaussian_similarity: float | None = None
    segmentation_status: str = field(default="not_applicable", metadata={"tsv": False})
    trace: np.ndarray | None = field(default=None, repr=False, metadata={"tsv": False})
    _ms2_records: tuple = field(default=(), repr=False, metadata={"tsv": False})
    _parent_mz: float = field(default=0.0, repr=False, metadata={"tsv": False})
    _peak_edges: tuple | None = field(default=None, repr=False, metadata={"tsv": False})


def summarize_trace(trace, config: Config, polarity="unknown"):
    """Use pinned MassCube splitting/finalization, including short traces.

    The adapter converts internal seconds to upstream minutes and attaches
    provenance. Scientific function bodies live unchanged in upstream.py.
    """
    trace = np.asarray(trace)
    if not len(trace):
        return
    first, stop, observed, _ = _span_state(trace[:, 1], 0, len(trace))
    if observed < config.min_scans:
        return
    if stop - first < 5:
        feature = _short_trace_feature(trace, first, stop, config, polarity)
        if feature is not None:
            if feature.scan_count >= config.min_scans:
                yield feature
            return
    yield from summarize_intervals(trace, config, polarity, Feature, segmentation_status)
    return


def _summarize_trace_reference(trace, config, polarity="unknown"):
    """Previous adapter retained as an independent optimization check."""
    from .upstream_adapter import ExtractionFeature, segmenter

    original = ExtractionFeature()
    original.signals = trace[:, :2]
    original.rt_seq = (trace[:, 2] / 60.0).tolist()
    original.scan_idx_seq = list(range(len(trace)))
    # Most traces are short or remain unsplit. Calculate their exported
    # measurements once, and reuse them when segmentation returns the same
    # unchanged object. Top-average/asymmetry are neither used nor exported.
    original.finalize(ta=False, a_score=False)
    parent_mass = float(original.mz)
    segments = [original]
    if config.segment:
        split = segmenter(config.smoothing_sigma, config.peak_prominence, config.peak_distance)
        terminal = set()
        for _ in range(config.segmentation_passes):
            next_segments = []
            for peak in segments:
                if peak in terminal:
                    next_segments.append(peak)
                    continue
                parts = split(
                    peak,
                    noise_tol=config.peak_noise_score_limit if config.peak_shape_guard else np.inf,
                )
                # The pinned splitter returns [peak] without mutating peak
                # when no split is possible. Repeating it cannot change the
                # result; newly created children still receive all passes.
                if len(parts) == 1 and parts[0] is peak:
                    terminal.add(peak)
                next_segments.extend(parts)
            segments = next_segments
            if all(peak in terminal for peak in segments):
                break
    for peak in segments:
        if peak is not original:
            peak.finalize(ta=False, a_score=False)
        if peak.length < config.min_scans:
            continue
        indices = np.asarray(peak.scan_idx_seq)
        apex_index = peak.scan_idx
        threshold = float(trace[apex_index, 3])
        height, area = float(peak.peak_height), float(peak.peak_area)
        yield Feature(
            mz=float(peak.mz),
            rt_apex_seconds=float(peak.rt * 60.0),
            rt_start_seconds=float(peak.rt_seq[0] * 60.0),
            rt_end_seconds=float(peak.rt_seq[-1] * 60.0),
            peak_height=height,
            peak_area=area,
            scan_count=int(peak.length),
            apex_noise_threshold=threshold,
            apex_threshold_ratio=height / threshold if threshold > 0 else None,
            polarity=polarity,
            intensity=height if config.quant_method == "peak_height" else area,
            quant_method=config.quant_method,
            gaussian_similarity=float(peak.gaussian_similarity),
            segmentation_status=segmentation_status(
                config.segment,
                peak.is_segmented,
                config.peak_shape_guard and peak.noise_score > config.peak_noise_score_limit,
            ),
            trace=np.column_stack((peak.signals[:, 1], trace[indices, 2]))
            if config.group_ions
            else None,
            _parent_mz=parent_mass,
        )


def _short_trace_feature(trace, first, stop, config, polarity):
    """Specialize upstream's length guards, keeping its measurement arithmetic.

    This is about the complete trimmed span, including interior missing-scan
    zeros, not merely the count of positive observations. With fewer than five
    samples both diagnostic scores are exactly zero and splitting always returns
    the same trace. No temporary upstream Feature/peak-shape array is needed.
    """
    if stop - first == 1:
        mass = float(np.float32(trace[first, 0]))
        intensity = np.float32(trace[first, 1])
        if not math.isfinite(mass) or not math.isfinite(float(intensity)) or intensity == 0:
            return None
        apex, count = 0, int(intensity > 0)
        height, area = float(int(intensity)), 0.0
        points = np.array([[intensity, trace[first, 2]]]) if config.group_ions else None
    else:
        signals = np.array(trace[first:stop, :2], dtype=np.float32)
        # Unusual direct callers can supply values that underflow or become
        # invalid on conversion. Leave those cases to the general upstream path,
        # including its existing error behavior.
        if not np.isfinite(signals).all() or not np.all(signals[[0, -1], 1] != 0):
            return None
        apex = int(np.argmax(signals[:, 1]))
        mass = float(signals[apex, 0])
        height = float(int(signals[apex, 1]))
        count = int(np.sum(signals[:, 1] > 0))
        points = (
            np.column_stack((signals[:, 1], trace[first:stop, 2])) if config.group_ions else None
        )
    # NumPy float32 callers can get different precision on scalar division.
    # Retain upstream's array division before
    # converting to Python floats, even for a singleton.
    minutes = (trace[first:stop, 2] / 60.0).tolist()
    # Keep NumPy's float32 y addition and float64 time arithmetic, including
    # upstream's Python integer conversion; do not reassociate this expression.
    if stop - first > 1:
        area = integrated_area(minutes, signals[:, 1])
    threshold = float(trace[first + apex, 3])
    return Feature(
        mz=mass,
        rt_apex_seconds=float(minutes[apex] * 60.0),
        rt_start_seconds=float(minutes[0] * 60.0),
        rt_end_seconds=float(minutes[-1] * 60.0),
        peak_height=height,
        peak_area=area,
        scan_count=count,
        apex_noise_threshold=threshold,
        apex_threshold_ratio=height / threshold if threshold > 0 else None,
        polarity=polarity,
        intensity=height if config.quant_method == "peak_height" else area,
        quant_method=config.quant_method,
        gaussian_similarity=0.0,
        segmentation_status="unsplit" if config.segment else "disabled",
        trace=points,
        _parent_mz=mass,
    )


def segmentation_status(enabled, split, guarded):
    if not enabled:
        return "disabled"
    result = "split" if split else "unsplit"
    return "guarded_" + result if guarded else result


def warmup():
    """Compile/load the small fixed kernel set once per long-lived worker."""
    from .noise import estimate_noise
    from .centroid import centroid_profile
    from .grouping import warmup_matching, trace_correlation
    from .grouping_spectra import extract_eic
    from .gapfill import eic_intensities
    from .noise_regions import regional_thresholds
    from .calibration import measure_representative
    from .gapfill_refinement import _eic
    from .summary import _span_state, maxima
    from .tuning import competing_peaks
    from .peak_shapes import gaussian_similarity as pilot_gaussian_similarity

    config = Config(min_scans=2)
    mz = np.array([100.0, 200.0])
    intensity = np.array([10.0, 20.0])
    threshold = estimate_noise(mz, intensity, config)[2]
    tracker = TraceTracker(config)
    tracker.add(mz, intensity, threshold, 0.0)
    tracker.add(mz, intensity, threshold, 1.0)
    tracker.flush()
    tracker.add_features(mz, intensity, threshold, 2.0)
    tracker.add_features(mz, intensity, threshold, 3.0)
    list(tracker.flush_features())
    centroid_profile(mz, intensity, np.zeros(2), threshold, 3, 4.0)
    centroid_profile(
        mz.astype(np.float32), intensity.astype(np.float32), np.zeros(2), threshold, 3, 4.0
    )
    # Production centroids are float32 column views; calibration uses float64
    # arrays. Expire two traces to also prime the reversed finished-ID view.
    signals = np.column_stack((mz, intensity)).astype(np.float32)
    production = TraceTracker(Config(max_gap_scans=0))
    production.add(signals[:, 0], signals[:, 1], threshold, 0.0)
    production.add(signals[:0, 0], signals[:0, 1], threshold[:0], 1.0)
    trace = np.column_stack((intensity, np.arange(2, dtype=float)))
    trace_correlation(trace, trace)
    pilot_gaussian_similarity(trace[:, 1], intensity)
    warmup_matching()
    shape = np.column_stack((np.arange(7), [1, 2, 4, 8, 4, 2, 1])).astype(np.float32)
    _span_state(shape[:, 1], 0, len(shape))
    _span_state(np.ones(1), 0, 1)
    maxima(shape[:, 1], config)
    points = np.column_stack((mz, intensity))
    points.flags.writeable = False  # Grouping reads a read-only memory map.
    extract_eic(points, np.array([0]), np.array([2]), np.array([0]), 100.0, 0.005)
    compact_points = points.astype(np.float32)
    compact_points.flags.writeable = False
    extract_eic(
        compact_points, np.array([0]), np.array([2]), np.array([0]), np.float32(100.0), 0.005
    )
    eic_intensities(mz, intensity, mz, np.full(len(mz), 0.01))
    regional_thresholds(mz, intensity, mz, config)
    list(summarize_trace(np.column_stack((mz, intensity, mz, threshold)), config))
    measure_representative(np.zeros((2, 4)), 0.05, 100.0, 0.8)
    competing_peaks(mz, intensity, np.arange(2, dtype=np.int64), 50.0, 0.025)
    _eic(mz, intensity, np.array([0, 2], dtype=np.int64), 100.0, 0.01)
