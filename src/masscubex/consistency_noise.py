"""Signal-preserving MS1 noise model with exact reference decision semantics.

Repeated intensity bands identify electronic structure. A local density mode
estimates a broader background without selecting a fixed fraction of points.
Mass recurrence can rescue points that overlap that background. Neither this
model nor its initial parameters are a validated replacement for the defaults.
"""

from collections import Counter
from dataclasses import dataclass
from functools import lru_cache
import math

import numpy as np
from scipy.ndimage import gaussian_filter1d
from scipy.signal import find_peaks

from .mass_tolerance import mass_radius, ms1_region_array
from .consistency_kernels import (
    add_recurrences,
    closest_precursor,
    group_intensities,
    has_distinct,
    histogram_counts,
    resolve_radii,
)


@dataclass(frozen=True)
class Parameters:
    bin_width_da: float = 100.0
    relative_band_width: float = 0.001
    minimum_repeats: int = 8
    repeat_span_da: float = 50.0
    minimum_mass_cells: int = 3
    mass_cell_da: float = 10.0
    pool_half_scans: int = 2
    recurrence_half_scans: int = 10
    maximum_pool_gap_seconds: float = 15.0
    minimum_component_points: int = 20
    minimum_distinct_intensities: int = 5
    log_histogram_width: float = 0.10
    histogram_smoothing: float = 1.0
    maximum_log_sd: float = 0.60
    minimum_recurrences: int = 2  # Other raw MS1 scans, not a feature-count gate.
    require_unique_mass_match: bool = True
    minimum_global_shape_points: int = 5


@dataclass(frozen=True)
class Component:
    lower: float
    upper: float
    mean: float
    sd: float
    points: int
    status: str

    @property
    def threshold(self):
        return self.mean + 3 * self.sd


@dataclass
class Decision:
    candidate_noise: np.ndarray
    coherent: np.ndarray
    threshold: np.ndarray
    components: dict
    ms2_supported: np.ndarray

    def keep(self, rescue=True):
        return (
            ~self.candidate_noise | self.coherent | self.ms2_supported
            if rescue
            else ~self.candidate_noise
        )


def intensity_groups(mz, intensity, parameters, *, published_rule=False):
    """Fixed-diameter groups prevent chaining through the intensity range.

    The published electronic decision uses >=4 points. Our candidate requires
    more support as scan size grows, plus dispersion across independent m/z
    cells. No intensities are rounded to manufacture ties.
    """
    order = np.argsort(intensity, kind="stable")
    minimum = (
        4
        if published_rule
        else max(parameters.minimum_repeats, math.ceil(math.log2(max(2, len(intensity)))))
    )
    rejected, _, lows, highs, sizes = group_intensities(
        mz,
        intensity,
        order,
        np.floor(mz / parameters.mass_cell_da),
        parameters.relative_band_width,
        minimum,
        parameters.repeat_span_da,
        parameters.minimum_mass_cells,
        published_rule,
    )
    groups = [
        (float(lows[i]), float(highs[i]), int(sizes[i])) for i in range(len(lows) - 1, -1, -1)
    ]
    return rejected, groups


def _supported_spread(mz, parameters):
    return len(mz) >= parameters.minimum_component_points and has_distinct(
        np.floor(mz / parameters.mass_cell_da), parameters.minimum_mass_cells
    )


def population_component(
    mz, intensity, parameters, *, status, lower=None, upper=None, minimum=None
):
    """Linear moments of evidence-selected noise; no percentile selection."""
    minimum = parameters.minimum_component_points if minimum is None else minimum
    if (
        len(intensity) < minimum
        or not has_distinct(intensity, min(parameters.minimum_distinct_intensities, minimum))
        or not has_distinct(np.floor(mz / parameters.mass_cell_da), parameters.minimum_mass_cells)
    ):
        return None
    values = intensity.astype(np.float64)
    origin = values[0]
    shifted = values - origin
    center = shifted.mean()
    return Component(
        float(values.min()) if lower is None else lower,
        float(values.max()) if upper is None else upper,
        origin + float(center),
        float(np.sqrt(np.mean((shifted - center) ** 2))),
        len(values),
        status,
    )


def fit_component(mz, intensity, floor_atom, parameters):
    """Estimate a compact, low intensity-density mode, separating floor atoms.

    Discover the mode in log intensity; use its right half-height width to
    estimate the component extent, including a left-censored mode. Moments
    are calculated in linear units from the resulting background population.
    Sparse/broad/non-separated cases explicitly have no fitted component.
    """
    body = intensity[~floor_atom]
    if len(body) < parameters.minimum_component_points or not has_distinct(
        body, parameters.minimum_distinct_intensities
    ):
        return None, "insufficient_body"
    log_values = np.log(body.astype(np.float64))
    width = parameters.log_histogram_width
    first = math.floor(float(log_values.min()) / width) * width - width
    last = math.ceil(float(log_values.max()) / width) * width + 2 * width
    edges = np.arange(first, last + width / 2, width)
    counts = histogram_counts(log_values, edges)
    density = gaussian_filter1d(
        np.pad(counts.astype(float), 2), parameters.histogram_smoothing, mode="constant"
    )
    centers = (edges[:-1] + edges[1:]) / 2
    centers = np.r_[
        centers[0] - 2 * width,
        centers[0] - width,
        centers,
        centers[-1] + width,
        centers[-1] + 2 * width,
    ]
    peaks, _ = find_peaks(density, prominence=max(1.0, float(density.max()) * 0.10))
    if not len(peaks):
        return None, "no_supported_mode"
    # The first supported mode supplies a low-background candidate. Evidence
    # tests below can reject it; a background component is never forced.
    peak = int(peaks[0])
    half_height = density[peak] / 2
    right = peak + 1
    while right < len(density) and density[right] > half_height:
        if right > peak + 1 and density[right] > density[right - 1]:
            return None, "unseparated_mode"
        right += 1
    if right == len(density):
        return None, "unbounded_mode"
    height_drop = density[right - 1] - density[right]
    fraction = (density[right - 1] - half_height) / height_drop if height_drop else 0.0
    crossing = centers[right - 1] + fraction * width
    log_sd = (crossing - centers[peak]) / math.sqrt(2 * math.log(2))
    if not 0 < log_sd <= parameters.maximum_log_sd:
        return None, "diffuse_mode"
    lower = math.exp(float(centers[peak]) - 3 * log_sd)
    upper = math.exp(float(centers[peak]) + 3 * log_sd)
    background = intensity <= upper
    selected_mz, selected = mz[background], intensity[background].astype(np.float64)
    if not _supported_spread(selected_mz, parameters) or not has_distinct(
        selected, parameters.minimum_distinct_intensities
    ):
        return None, "insufficient_spatial_support"
    # Shifted moments remain stable when the export has a large baseline.
    origin = selected[0]
    shifted = selected - origin
    mean = origin + float(shifted.mean())
    sd = float(np.sqrt(np.mean((shifted - shifted.mean()) ** 2)))
    return Component(lower, upper, mean, sd, len(selected), "density_mode"), "density_mode"


def _neighbors(frames, half_window, maximum_gap):
    by_polarity = {}
    for index, frame in enumerate(frames):
        by_polarity.setdefault(frame.polarity, []).append(index)
    result = [[] for _ in frames]
    for indices in by_polarity.values():
        for position, index in enumerate(indices):
            left, right = (
                max(0, position - half_window),
                min(len(indices), position + half_window + 1),
            )
            result[index] = [
                j
                for j in indices[left:right]
                if abs(frames[j].rt - frames[index].rt) <= maximum_gap
            ]
    return result


@lru_cache(maxsize=128)
def _mass_regions(config):
    return ms1_region_array(config)


def mass_recurrence(frame, others, config, minimum, require_unique=True):
    """Match a fixed current mass in other same-polarity raw MS1 scans.

    This protects uncertain points; it does not filter all short features.
    Using a fixed target prevents a sequence of nearby masses from walking.
    """
    count = np.zeros(len(frame.mz), dtype=np.int16)
    radii = resolve_radii(
        frame.mz,
        float(config.mz_tolerance_da),
        float(config.mz_tolerance_ppm),
        _mass_regions(config),
    )
    for other in others:
        add_recurrences(frame.mz, other.mz, radii, count, minimum, require_unique)
    return count >= minimum


def _zeros(frames, dtype, storage, name):
    return (
        storage.zeros(name, dtype)
        if storage is not None
        else [np.zeros(len(f.mz), dtype=dtype) for f in frames]
    )


def precursor_protection(frames, records, configs, parameters=None, *, storage=None):
    """Use raw acquisition evidence, never baseline feature assignments.

    Protect only the closest centroid inside the ordinary learned MS1 radius,
    in the referenced/preceding MS1 survey and bounded neighboring surveys.
    No wider mass window or molecular identity is inferred.
    """
    parameters = parameters or Parameters()
    masks = _zeros(frames, bool, storage, "ms2")
    references = {f.native_id: i for i, f in enumerate(frames) if f.native_id}
    streams = {}
    positions = {}
    for i, frame in enumerate(frames):
        stream = streams.setdefault(frame.polarity, [])
        positions[i] = len(stream)
        stream.append(i)
    times = {p: np.array([frames[j].rt for j in indices]) for p, indices in streams.items()}
    for record in records:
        mass = record.precursor_mz
        if mass is None or not np.isfinite(mass) or mass <= 0 or record.precursor_count != 1:
            continue
        anchor = references.get(record.precursor_ref)
        polarity = frames[anchor].polarity if anchor is not None else record.polarity
        if polarity == "unknown" and len(streams) == 1:
            polarity = next(iter(streams))
        if polarity not in streams or polarity not in configs:
            continue
        if anchor is None:
            if record.rt is None:
                continue
            position = int(np.searchsorted(times[polarity], record.rt, side="right")) - 1
            if position < 0:
                continue
            anchor = streams[polarity][position]
        elif record.polarity not in {"unknown", polarity}:
            continue
        position = positions[anchor]
        indices = streams[polarity]
        radius = mass_radius(mass, configs[polarity])
        for j in indices[
            max(0, position - parameters.recurrence_half_scans) : position
            + parameters.recurrence_half_scans
            + 1
        ]:
            if abs(frames[j].rt - frames[anchor].rt) > parameters.maximum_pool_gap_seconds:
                continue
            if (
                record.rt is not None
                and abs(record.rt - frames[j].rt) > parameters.maximum_pool_gap_seconds
            ):
                continue
            closest = closest_precursor(frames[j].mz, mass, radius)
            if closest >= 0:
                masks[j][closest] = True
    return masks


@dataclass
class _TrainingScan:
    mz: np.ndarray
    intensity: np.ndarray
    floor: np.ndarray
    repeated: np.ndarray
    bins: dict


def _training_scan(frame, coherent, ms2, floor, repeated, parameters):
    keep = ~(coherent | ms2)
    mz = frame.mz[keep]
    keys, starts = np.unique(np.floor(mz / parameters.bin_width_da), return_index=True)
    ends = np.r_[starts[1:], len(mz)]
    return _TrainingScan(
        mz,
        frame.intensity[keep],
        floor[keep],
        repeated[keep],
        {key: (start, end) for key, start, end in zip(keys, starts, ends)},
    )


def _pool_training(scans, key=None):
    """Concatenate in the reference order: survey first, then increasing m/z."""
    slices = [slice(*s.bins.get(key, (0, 0))) if key is not None else slice(None) for s in scans]
    return tuple(
        np.concatenate([getattr(s, name)[part] for s, part in zip(scans, slices)])
        for name in ("mz", "intensity", "floor", "repeated")
    )


def evaluate(
    frames, configs, parameters=None, *, ms2_evidence=None, storage=None, include_published=True
):
    """Return point decisions without modifying scans or their intensities."""
    parameters = parameters or Parameters()
    if ms2_evidence is None:
        ms2_evidence = _zeros(frames, bool, storage, "ms2")
    if len(ms2_evidence) != len(frames) or any(
        mask.shape != f.mz.shape for f, mask in zip(frames, ms2_evidence)
    ):
        raise ValueError("MS2 evidence must match the raw centroid frames")
    for frame in frames:
        if not frame.centroided:
            raise ValueError("The prototype requires mass-detected centroid ions")
        if (
            frame.mz.ndim != 1
            or frame.intensity.shape != frame.mz.shape
            or np.any(np.diff(frame.mz) < 0)
            or np.any(frame.mz <= 0)
            or np.any(frame.intensity <= 0)
            or not np.isfinite(frame.mz).all()
            or not np.isfinite(frame.intensity).all()
        ):
            raise ValueError("Invalid centroid arrays")
    neighbors = _neighbors(frames, parameters.pool_half_scans, parameters.maximum_pool_gap_seconds)
    recurrence_neighbors = _neighbors(
        frames, parameters.recurrence_half_scans, parameters.maximum_pool_gap_seconds
    )
    coherent_masks = _zeros(frames, bool, storage, "coherent")
    for index, frame in enumerate(frames):
        coherent_masks[index][:] = mass_recurrence(
            frame,
            [frames[j] for j in recurrence_neighbors[index] if j != index],
            configs[frame.polarity],
            parameters.minimum_recurrences,
            parameters.require_unique_mass_match,
        )
    repeat_masks = _zeros(frames, bool, storage, "repeated")
    floor_masks = _zeros(frames, bool, storage, "floor")
    published_masks = _zeros(frames, bool, storage, "paper") if include_published else None
    for index, frame in enumerate(frames):
        order = np.argsort(frame.intensity, kind="stable")
        minimum = max(
            parameters.minimum_repeats, math.ceil(math.log2(max(2, len(frame.intensity))))
        )
        repeated, paper, lows, highs, _ = group_intensities(
            frame.mz,
            frame.intensity,
            order,
            np.floor(frame.mz / parameters.mass_cell_da),
            parameters.relative_band_width,
            minimum,
            parameters.repeat_span_da,
            parameters.minimum_mass_cells,
        )
        if len(lows):
            floor_masks[index][:] = (frame.intensity >= lows[-1]) & (frame.intensity <= highs[-1])
        repeat_masks[index][:] = repeated
        if include_published:
            published_masks[index][:] = paper

    # A training survey is filtered once and released after its last neighbor
    # uses it. Working copies stay bounded by the temporal pool, not file size.
    uses = np.zeros(len(frames), dtype=np.int64)
    for indices in neighbors:
        for j in indices:
            uses[j] += 1
    training = {}
    uniform_mass_dtype = len({f.mz.dtype.str for f in frames}) <= 1
    all_candidates = _zeros(frames, bool, storage, "candidates")
    all_thresholds = _zeros(frames, np.float64, storage, "thresholds")
    decisions, counts = [], Counter()
    for index, frame in enumerate(frames):
        pooled = neighbors[index]
        for j in pooled:
            if j not in training:
                training[j] = _training_scan(
                    frames[j],
                    coherent_masks[j],
                    ms2_evidence[j],
                    floor_masks[j],
                    repeat_masks[j],
                    parameters,
                )
        surveys = [training[j] for j in pooled]
        # Mixed-precision callers need the reference's promotion *before*
        # binning. Native extraction scans uniformly use upstream float32.
        global_pool = None if uniform_mass_dtype else _pool_training(surveys)
        global_shape = None
        global_fitted = False
        current_keys, starts = np.unique(
            np.floor(frame.mz / parameters.bin_width_da), return_index=True
        )
        ends = np.r_[starts[1:], len(frame.mz)]
        candidate = all_candidates[index]
        thresholds = all_thresholds[index]
        components = {}
        for key, start, end in zip(current_keys, starts, ends):
            if uniform_mass_dtype:
                mz, values, floor, repeated = _pool_training(surveys, key)
            else:
                region = np.floor(global_pool[0] / parameters.bin_width_da) == key
                mz, values, floor, repeated = (a[region] for a in global_pool)
            component = population_component(
                mz[repeated], values[repeated], parameters, status="repeated_population"
            )
            if component is not None:
                status = component.status
            else:
                component, status = fit_component(mz, values, floor, parameters)
            if component is None:
                if not global_fitted:
                    if global_pool is None:
                        global_pool = _pool_training(surveys)
                    global_shape, _ = fit_component(*global_pool[:3], parameters)
                    global_fitted = True
                if global_shape is not None:
                    compatible = values <= global_shape.upper
                    component = population_component(
                        mz[compatible],
                        values[compatible],
                        parameters,
                        status="global_shape_local_moments",
                        lower=global_shape.lower,
                        upper=global_shape.upper,
                        minimum=parameters.minimum_global_shape_points,
                    )
                    if component is not None:
                        status = component.status
            counts[status] += 1
            values_current = frame.intensity[start:end]
            repeated_current = repeat_masks[index][start:end]
            if component is None:
                if repeated_current.any():
                    atoms = values_current[repeated_current].astype(float)
                    thresholds[start:end] = float(atoms.mean() + 3 * atoms.std())
                candidate[start:end] = repeated_current & (values_current <= thresholds[start:end])
                components[int(key)] = {"status": status}
                continue
            components[int(key)] = {
                "status": status,
                "lower": component.lower,
                "upper": component.upper,
                "mean": component.mean,
                "sd": component.sd,
                "threshold": component.threshold,
                "points": component.points,
            }
            thresholds[start:end] = component.threshold
            in_body = (values_current >= component.lower) & (values_current <= component.upper)
            candidate[start:end] = (repeated_current | in_body) & (
                values_current <= component.threshold
            )
        decisions.append(
            Decision(candidate, coherent_masks[index], thresholds, components, ms2_evidence[index])
        )
        for j in pooled:
            uses[j] -= 1
            if uses[j] == 0:
                del training[j]
    return decisions, published_masks, dict(counts)
