"""MS2 association without loading an entire file's MS1 or fragment spectra."""

from collections import Counter
from dataclasses import dataclass

import numpy as np

from .readers import MS2Scan
from .tuning import config_for
from .mass_tolerance import mass_radius
from .upstream_adapter import preprocess_ms2


class MS2Collector:
    def __init__(self, config, configs=None):
        self.config = config
        self.configs = configs
        self.records = []

    def ms1_stream(self, scans):
        """Collect compact MS2 records while yielding MS1 scans to the extractor."""
        previous = {}
        for scan in scans:
            if not isinstance(scan, MS2Scan):
                previous[scan.polarity] = scan
                yield scan
                continue
            if scan.polarity == "unknown":
                scan.polarity_source = "unknown"
                if scan.precursor_ref:
                    parents = [s for s in previous.values() if s.native_id == scan.precursor_ref]
                    if len(parents) == 1 and parents[0].polarity != "unknown":
                        scan.polarity = parents[0].polarity
                        scan.polarity_source = "precursor_reference"
            effective = config_for(scan.polarity, self.config, self.configs)
            preprocess_ms2(scan, effective)
            if scan.precursor_intensity is not None and scan.precursor_intensity > 0:
                scan.intensity_source = "metadata"
            else:
                scan.precursor_intensity = None
                parent = previous.get(scan.polarity)
                if parent is None and len(previous) == 1:
                    candidate = next(iter(previous.values()))
                    if "unknown" in (candidate.polarity, scan.polarity):
                        parent = candidate
                if (
                    parent is not None
                    and scan.rt is not None
                    and scan.precursor_mz is not None
                    and 0 <= scan.rt - parent.rt <= effective.precursor_max_age_seconds
                    and (not scan.precursor_ref or scan.precursor_ref == parent.native_id)
                ):
                    tolerance = mass_radius(scan.precursor_mz, effective, "ms2")
                    left = np.searchsorted(parent.mz, scan.precursor_mz - tolerance)
                    right = np.searchsorted(parent.mz, scan.precursor_mz + tolerance, side="right")
                    if right > left:
                        if parent.centroided:
                            closest = left + np.argmin(
                                abs(parent.mz[left:right] - scan.precursor_mz)
                            )
                            value = float(parent.intensity[closest])
                        else:
                            value = float(np.max(parent.intensity[left:right]))
                        if value > 0:
                            scan.precursor_intensity = value
                            scan.intensity_source = "preceding_ms1"
                            scan._intensity_polarity = parent.polarity
            if scan.mz is not None and len(scan.mz):
                if scan.fragment_intensity is None:
                    scan.fragment_intensity = float(np.sum(scan.intensity, dtype=np.float64))
                keep = np.flatnonzero(
                    scan.intensity >= scan.intensity.max() * self.config.isf_min_relative_intensity
                )
                if self.config.isf_max_peaks and len(keep) > self.config.isf_max_peaks:
                    order = np.argsort(-scan.intensity[keep], kind="stable")
                    keep = np.sort(keep[order[: self.config.isf_max_peaks]])
                scan.mz, scan.intensity = scan.mz[keep], scan.intensity[keep]
            self.records.append(scan)


@dataclass(frozen=True, slots=True)
class MS2Assignment:
    feature_index: int | None
    reason: str


def ineligible_reason(scan, config):
    """Return only metadata/explicit-policy exclusions, never a guessed DIA label."""
    if scan.precursor_count < 1 or scan.precursor_mz is None:
        return "missing_precursor_mz"
    if scan.rt is None:
        return "missing_retention_time"
    if (
        config.ms2_max_isolation_width > 0
        and scan.isolation_width is not None
        and scan.isolation_width > config.ms2_max_isolation_width
    ):
        return "isolation_width_limit"
    return None


def associate_ms2(features, records, config, stats, configs=None):
    """Assign each single-precursor MS2 to at most one chromatographic feature.

    Isolation width is unrestricted by default. Return an outcome for every
    record so unmatched scans can seed features or be reported explicitly.
    """
    if not records:
        apply_ms2_assignments(features, records, [], config, stats)
        return []
    masses = np.array([feature.mz for feature in features])
    order = np.argsort(masses, kind="stable")
    sorted_masses = masses[order]
    # These float32 minute bounds and heights are invariant across all MS2
    # queries. Keep scalar comparisons below, preserving NumPy's installed
    # scalar-promotion rules, but do not rebuild an array for every candidate.
    edges = np.asarray(
        [
            f._peak_edges or (f.rt_start_seconds / 60.0, f.rt_end_seconds / 60.0)
            for f in features
        ], dtype=np.float32,
    ).reshape(-1, 2)
    heights = np.asarray([f.peak_height or 0.0 for f in features], dtype=np.float32)
    assignments = []
    polarities = set(stats.get("ms1_polarities", [f.polarity for f in features]))
    mixed = len(polarities) > 1 or (
        "unknown" in polarities and len({s.polarity for s in records} - {"unknown"}) > 1
    )
    for scan in records:
        # Earlier scans cannot reveal whether another polarity appears later.
        # Recheck fallback evidence now that the complete file is known.
        if (
            mixed
            and scan.intensity_source == "preceding_ms1"
            and (scan.polarity == "unknown" or scan._intensity_polarity != scan.polarity)
        ):
            scan.precursor_intensity = None
            scan.intensity_source = "unavailable"
        effective = config_for(scan.polarity, config, configs)
        reason = ineligible_reason(scan, effective)
        if reason is not None:
            assignments.append(MS2Assignment(None, reason))
            continue
        if scan.polarity == "unknown" and mixed:
            assignments.append(MS2Assignment(None, "ambiguous_polarity"))
            continue
        tolerance = mass_radius(scan.precursor_mz, effective, "ms2")
        left = np.searchsorted(sorted_masses, scan.precursor_mz - tolerance)
        right = np.searchsorted(sorted_masses, scan.precursor_mz + tolerance, side="right")
        mass_candidates = order[left:right]
        polarity_candidates = [
            int(i)
            for i in mass_candidates
            if (
                features[i].polarity == scan.polarity
                or (not mixed and (scan.polarity == "unknown" or features[i].polarity == "unknown"))
            )
            and abs(features[i].mz - scan.precursor_mz) < tolerance
        ]

        candidates = [
            i
            for i in polarity_candidates
            if edges[i, 0] - effective.ms2_rt_tolerance / 60.0
            <= scan.rt / 60.0
            <= edges[i, 1] + effective.ms2_rt_tolerance / 60.0
        ]
        if (
            scan.polarity == "unknown"
            and len({features[i].polarity for i in candidates} - {"unknown"}) > 1
        ):
            assignments.append(MS2Assignment(None, "ambiguous_polarity"))
            continue
        if not candidates:
            reason = (
                "no_feature_mz"
                if not len(mass_candidates)
                else "no_feature_polarity"
                if not polarity_candidates
                else "no_feature_rt"
            )
            assignments.append(MS2Assignment(None, reason))
            continue

        best = candidates[int(np.argmax(heights[candidates]))]
        assignments.append(MS2Assignment(best, "matched"))
    apply_ms2_assignments(features, records, assignments, config, stats)
    return assignments


def apply_ms2_assignments(features, records, assignments, config, stats):
    """Publish scan references and best spectra from the single authoritative mapping."""
    if len(assignments) != len(records):
        raise ValueError("Every MS2 record must have exactly one assignment outcome")
    assigned = [[] for _ in features]
    unresolved = []
    reasons = Counter()
    for i, (scan, assignment) in enumerate(zip(records, assignments)):
        if assignment.feature_index is not None:
            assigned[assignment.feature_index].append(i)
        else:
            reasons[assignment.reason] += 1
            unresolved.append(
                {
                    "spectrum_index": scan.spectrum_index,
                    "native_id": scan.native_id,
                    "scan_number": scan.scan_number,
                    "rt_seconds": scan.rt,
                    "precursor_mz": scan.precursor_mz,
                    "polarity": scan.polarity,
                    "precursor_count": scan.precursor_count,
                    "isolation_width": scan.isolation_width,
                    "reason": assignment.reason,
                    "polarity_source": scan.polarity_source,
                }
            )
    stats.update(
        ms2_scans=len(records),
        ms2_assigned=len(records) - len(unresolved),
        ms2_unassigned=len(unresolved),
        ms2_skipped_wide_isolation=reasons["isolation_width_limit"],
        ms2_skipped_multiprecursor=reasons["multiple_precursors"],
        ms2_missing_metadata=reasons["missing_precursor_mz"] + reasons["missing_retention_time"],
        ms2_ambiguous_polarity=reasons["ambiguous_polarity"],
        ms2_unassigned_reasons=dict(sorted(reasons.items())),
        unassigned_ms2=unresolved,
        ms2_polarity_inferences=[
            {
                "spectrum_index": s.spectrum_index,
                "native_id": s.native_id,
                "polarity": s.polarity,
                "precursor_ref": s.precursor_ref,
            }
            for s in records
            if s.polarity_source == "precursor_reference"
        ],
    )
    # Configuration is recorded once per file. Record only exceptional rows
    # whose missing intensity required the documented nearest-apex fallback.
    fallback_ids = []
    stats["best_ms2_apex_fallback_feature_ids"] = fallback_ids
    # The table keeps only the best native scan number. Preserve a resolvable
    # reference for the exceptional files where that number is absent/ambiguous.
    identifier_fallbacks = []
    stats["best_ms2_identifier_fallbacks"] = identifier_fallbacks
    for feature_id, (feature, indices) in enumerate(zip(features, assigned), 1):
        if not indices and feature.ms2_count == 0:
            # Fresh unassigned features already have empty provenance. Avoid
            # rebuilding three empty tuples and resetting eight fields twice.
            continue
        indices.sort(key=lambda i: records[i].spectrum_index)
        feature._ms2_records = tuple(indices)
        feature.ms2_count = len(indices)
        feature.ms2_scan_numbers = tuple(records[i].scan_number for i in indices)
        feature.ms2_native_ids = tuple(records[i].native_id for i in indices)
        feature.ms2_spectrum_indices = tuple(records[i].spectrum_index for i in indices)
        feature.best_ms2_scan_number = None
        feature.best_ms2_native_id = None
        feature.best_ms2_spectrum_index = None
        feature.best_ms2_selection_method = None
        feature.best_ms2_score = None
        feature.best_ms2_fragment_intensity = None
        feature.best_ms2_precursor_intensity = None
        feature.best_ms2_precursor_intensity_source = None
        if not indices:
            continue

        def value(scan):
            return (
                scan.fragment_intensity
                if config.best_ms2 == "fragment_intensity"
                else scan.precursor_intensity
            )

        def ranking(i):
            scan = records[i]
            intensity = value(scan)
            intensity = intensity if intensity is not None else -1.0
            distance = abs(scan.rt - feature.rt_apex_seconds)
            return (
                (-distance, intensity, -scan.spectrum_index)
                if config.best_ms2 == "apex"
                else (intensity, -scan.spectrum_index)
            )

        best = records[max(indices, key=ranking)]
        feature.best_ms2_scan_number = best.scan_number
        feature.best_ms2_native_id = best.native_id
        feature.best_ms2_spectrum_index = best.spectrum_index
        if best.scan_number is None or feature.ms2_scan_numbers.count(best.scan_number) > 1:
            identifier_fallbacks.append(
                {
                    "feature_id": feature_id,
                    "native_id": best.native_id,
                    "spectrum_index": best.spectrum_index,
                }
            )
        method = config.best_ms2
        if method != "apex" and value(best) is None:
            method = "apex_fallback"
            fallback_ids.append(feature_id)
        feature.best_ms2_selection_method = method
        feature.best_ms2_score = (
            abs(best.rt - feature.rt_apex_seconds)
            if method in {"apex", "apex_fallback"}
            else value(best)
        )
        feature.best_ms2_fragment_intensity = best.fragment_intensity
        feature.best_ms2_precursor_intensity = best.precursor_intensity
        feature.best_ms2_precursor_intensity_source = best.intensity_source
