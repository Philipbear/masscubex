"""Compiled online ion traces backed by recyclable array pools.

Only active traces are retained. Finished traces are segmented and emitted;
there is no Python object or Python append call per centroid.
"""

from dataclasses import dataclass, field
from functools import lru_cache

import numpy as np
from numba import njit
from scipy.ndimage import correlate1d
from scipy.signal import find_peaks

from .config import Config


@njit(cache=True, nogil=True)
def _append(slot, mz, intensity, rt, threshold, ti, nodes, links, free_nodes, available):
    available -= 1
    node = free_nodes[available]
    nodes[node, 0] = mz
    nodes[node, 1] = intensity
    nodes[node, 2] = rt
    nodes[node, 3] = threshold
    links[node] = -1
    if ti[slot, 2] == 0:
        ti[slot, 0] = node
    else:
        links[ti[slot, 1]] = node
    ti[slot, 1] = node
    ti[slot, 2] += 1
    return available


@njit(cache=True, nogil=True)
def _advance(
    mz,
    intensity,
    threshold,
    rt,
    ppm,
    da,
    max_gaps,
    max_seconds,
    active,
    nactive,
    tf,
    ti,
    nodes,
    links,
    free_tracks,
    ft,
    free_nodes,
    fn,
):
    old = active[:nactive].copy()
    # Stronger traces claim peaks first; stable ties preserve acquisition order.
    priority = np.argsort(-tf[old, 2], kind="mergesort")
    used = np.zeros(mz.size, dtype=np.bool_)
    done = np.empty(nactive, dtype=np.int64)
    ndone = 0
    nactive = 0
    for k in priority:
        slot = old[k]
        if rt - tf[slot, 3] > max_seconds:
            done[ndone] = slot
            ndone += 1
            continue
        target = tf[slot, 0]
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
            weight = tf[slot, 1] + intensity[best]
            tf[slot, 0] += (mz[best] - target) * (intensity[best] / weight)
            tf[slot, 1] = weight
            tf[slot, 2] = intensity[best]
            tf[slot, 3] = rt
            ti[slot, 3] += 1
            ti[slot, 4] = 0
            fn = _append(
                slot,
                mz[best],
                intensity[best],
                rt,
                threshold[best],
                ti,
                nodes,
                links,
                free_nodes,
                fn,
            )
        else:
            ti[slot, 4] += 1
            if ti[slot, 4] > max_gaps:
                done[ndone] = slot
                ndone += 1
                continue
            fn = _append(slot, target, 0.0, rt, 0.0, ti, nodes, links, free_nodes, fn)
        active[nactive] = slot
        nactive += 1
    for j in range(mz.size):
        if used[j]:
            continue
        ft -= 1
        slot = free_tracks[ft]
        ti[slot, :] = 0
        ti[slot, 3] = 1
        tf[slot, 0] = mz[j]
        tf[slot, 1] = intensity[j]
        tf[slot, 2] = intensity[j]
        tf[slot, 3] = rt
        fn = _append(slot, mz[j], intensity[j], rt, threshold[j], ti, nodes, links, free_nodes, fn)
        active[nactive] = slot
        nactive += 1
    return nactive, ft, fn, done[:ndone]


@njit(cache=True, nogil=True)
def _drain(done, minimum, ti, nodes, links, free_tracks, ft, free_nodes, fn):
    traces = []
    for slot in done:
        keep = ti[slot, 3] >= minimum
        trace = np.empty((ti[slot, 2] if keep else 0, 4))
        node = ti[slot, 0]
        for j in range(ti[slot, 2]):
            if keep:
                trace[j, :] = nodes[node, :]
            free_nodes[fn] = node
            fn += 1
            node = links[node]
        free_tracks[ft] = slot
        ft += 1
        if keep:
            traces.append(trace)
    return traces, ft, fn


class TraceTracker:
    """State for one file and polarity; never shared between workers."""

    def __init__(self, config: Config, capacity=1024, node_capacity=8192):
        self.config = config
        self.tf = np.zeros((capacity, 4))
        self.ti = np.zeros((capacity, 5), dtype=np.int64)
        self.nodes = np.empty((node_capacity, 4))
        self.links = np.empty(node_capacity, dtype=np.int64)
        self.free_tracks = np.arange(capacity, dtype=np.int64)
        self.free_nodes = np.arange(node_capacity, dtype=np.int64)
        self.active = np.empty(capacity, dtype=np.int64)
        self.ft, self.fn, self.nactive = capacity, node_capacity, 0
        self.last_rt = -np.inf
        self.peak_active_traces = 0
        self.peak_live_nodes = 0

    def _ensure(self, peaks):
        if self.ft < peaks:
            old = len(self.tf)
            size = max(old * 2, old + peaks - self.ft)
            self.tf = np.vstack((self.tf, np.zeros((size - old, 4))))
            self.ti = np.vstack((self.ti, np.zeros((size - old, 5), dtype=np.int64)))
            free = np.empty(size, dtype=np.int64)
            free[: self.ft] = self.free_tracks[: self.ft]
            free[self.ft : self.ft + size - old] = np.arange(old, size)
            self.free_tracks = free
            self.ft += size - old
            active = np.empty(size, dtype=np.int64)
            active[: self.nactive] = self.active[: self.nactive]
            self.active = active
        needed = self.nactive + peaks
        if self.fn < needed:
            old = len(self.nodes)
            size = max(old * 2, old + needed - self.fn)
            self.nodes = np.vstack((self.nodes, np.empty((size - old, 4))))
            self.links = np.concatenate((self.links, np.empty(size - old, dtype=np.int64)))
            free = np.empty(size, dtype=np.int64)
            free[: self.fn] = self.free_nodes[: self.fn]
            free[self.fn : self.fn + size - old] = np.arange(old, size)
            self.free_nodes = free
            self.fn += size - old

    def _finish(self, done):
        traces, self.ft, self.fn = _drain(
            done,
            self.config.min_scans,
            self.ti,
            self.nodes,
            self.links,
            self.free_tracks,
            self.ft,
            self.free_nodes,
            self.fn,
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
        self.nactive, self.ft, self.fn, done = _advance(
            mz,
            intensity,
            threshold,
            float(rt),
            float(self.config.mz_tolerance_ppm),
            float(self.config.mz_tolerance_da),
            self.config.max_gap_scans,
            float(self.config.max_gap_seconds),
            self.active,
            self.nactive,
            self.tf,
            self.ti,
            self.nodes,
            self.links,
            self.free_tracks,
            self.ft,
            self.free_nodes,
            self.fn,
        )
        self.peak_active_traces = max(self.peak_active_traces, self.nactive)
        self.peak_live_nodes = max(self.peak_live_nodes, len(self.nodes) - self.fn)
        return self._finish(done)

    def flush(self):
        done = self.active[: self.nactive].copy()
        self.nactive = 0
        return self._finish(done)


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


@lru_cache(maxsize=512)
def _gaussian_weights(sigma):
    radius = int(4.0 * sigma + 0.5)
    positions = np.arange(-radius, radius + 1)
    weights = np.exp((-0.5 / (sigma * sigma)) * positions**2)
    weights /= weights.sum()
    weights.flags.writeable = False
    return weights


def _smooth_padded(signal, sigma):
    """Same zero pad and order-0 reflect Gaussian as scipy.gaussian_filter1d.

    Cache only the small kernel, never a signal or a peak result. A bounded
    cache avoids rebuilding identical kernels for every chromatographic trace.
    """
    padded = np.empty(len(signal) + 2, dtype=np.float64)
    padded[0] = padded[-1] = 0.0
    padded[1:-1] = signal
    return correlate1d(padded, _gaussian_weights(float(sigma)), mode="reflect")


@njit(cache=True, nogil=True)
def _has_multiple_maxima(values):
    """Count local maxima with SciPy's plateau definition, stopping after two.

    With zero/one maximum, prominence and distance cannot produce a split, so
    the expensive peak-property calculation is unnecessary. Multiple maxima
    always use the original SciPy prominence/distance algorithm.
    """
    peaks = 0
    i = 1
    while i < len(values) - 1:
        if values[i] > values[i - 1]:
            end = i
            while end < len(values) - 1 and values[end] == values[end + 1]:
                end += 1
            if end < len(values) - 1 and values[end] > values[end + 1]:
                peaks += 1
                if peaks > 1:
                    return True
            i = end
        i += 1
    return False


@njit(cache=True, nogil=True)
def _median_spacing(rt):
    # Trace columns are strided; Numba's np.diff requires a contiguous input.
    differences = np.empty(len(rt) - 1, dtype=np.float64)
    for i in range(len(differences)):
        differences[i] = rt[i + 1] - rt[i]
    return np.median(differences)


def summarize_trace(trace, config: Config, polarity="unknown"):
    """Split at chromatographic valleys; quantify original, un-smoothed signals.

    Adjacent segments share a valley point, but no trapezoid interval is counted
    twice. Area is intensity * seconds, and missing MS1 scans contribute zeros.
    """
    valid = np.flatnonzero(trace[:, 1] > 0)
    if len(valid) < config.min_scans:
        return
    trace = trace[valid[0] : valid[-1] + 1]
    intensity = trace[:, 1]
    cuts = [0]
    if config.segment:
        rt = trace[:, 2]
        physical = config.smoothing_seconds > 0 or config.peak_distance_seconds > 0
        if physical and len(rt) > 1:
            step = float(_median_spacing(rt))
            size = min(len(rt) * 4, max(2, int(np.ceil((rt[-1] - rt[0]) / step - 1e-9)) + 1))
            grid = np.linspace(rt[0], rt[-1], size)
            step = float(grid[1] - grid[0])
            signal = np.interp(grid, rt, intensity)
            sigma = (
                config.smoothing_seconds / step
                if config.smoothing_seconds > 0
                else config.smoothing_sigma
            )
            distance = (
                max(1, int(np.ceil(config.peak_distance_seconds / step - 1e-9)))
                if config.peak_distance_seconds > 0
                else config.peak_distance
            )
        else:
            grid, signal, sigma, distance = (
                rt,
                intensity,
                config.smoothing_sigma,
                config.peak_distance,
            )
        smooth = _smooth_padded(signal, sigma)
        peaks = np.empty(0, dtype=np.int64)
        if _has_multiple_maxima(smooth):
            peaks, _ = find_peaks(
                smooth, prominence=smooth.max() * config.peak_prominence, distance=distance
            )
        peaks = peaks - 1
        peaks = peaks[(peaks >= 0) & (peaks < len(signal))]
        for left, right in zip(peaks[:-1], peaks[1:]):
            valley = left + int(np.argmin(smooth[left + 1 : right + 2]))
            if physical:
                index = int(np.searchsorted(rt, grid[valley]))
                index = min(index, len(rt) - 1)
                if index and abs(rt[index - 1] - grid[valley]) <= abs(rt[index] - grid[valley]):
                    index -= 1
                valley = index
            if cuts[-1] < valley < len(intensity) - 1:
                cuts.append(valley)
    cuts.append(len(trace) - 1)
    for left, right in zip(cuts[:-1], cuts[1:]):
        part = trace[left : right + 1]
        y = part[:, 1]
        count = int(np.count_nonzero(y))
        if count < config.min_scans:
            continue
        apex = int(np.argmax(y))
        rt = part[:, 2]
        height = float(y[apex])
        threshold = float(part[apex, 3])
        area = float(np.dot(np.diff(rt), (y[:-1] + y[1:]) * 0.5))
        yield Feature(
            mz=float(np.dot(part[:, 0], y / y.sum())),
            rt_apex_seconds=float(rt[apex]),
            rt_start_seconds=float(rt[0]),
            rt_end_seconds=float(rt[-1]),
            peak_height=height,
            peak_area=area,
            scan_count=count,
            apex_noise_threshold=threshold,
            apex_threshold_ratio=height / threshold if threshold > 0 else None,
            polarity=polarity,
            intensity=height if config.quant_method == "peak_height" else area,
            quant_method=config.quant_method,
            trace=part[:, 1:3].copy() if config.group_ions else None,
        )


def warmup():
    """Compile/load the small fixed kernel set once per long-lived worker."""
    from .noise import estimate_noise
    from .centroid import centroid_profile
    from .grouping import trace_correlation, _candidate_pairs
    from .gapfill import eic_intensities

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
    _has_multiple_maxima(intensity)
    _median_spacing(mz)
    ids = np.arange(2, dtype=np.int64)
    list(_candidate_pairs(ids, ids, ids, ids + 1, mz, intensity, ids + 1, 5.0))
