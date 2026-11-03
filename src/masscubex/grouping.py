"""Optional, evidence-labelled ion relationships; no forced compound identities."""

from itertools import combinations

import numpy as np
from numba import njit

from .adducts import selected_adducts
from .tuning import config_for


@njit(cache=True, nogil=True)
def sorted_mass_windows(masses, lower, upper):
    """Linear sweep for ordered adduct queries; inclusive upper mass bounds.

    Extreme user tolerances or custom ion forms can produce unordered bounds.
    Those retain independent binary searches rather than assuming monotonicity.
    """
    if np.any(lower[1:] < lower[:-1]) or np.any(upper[1:] < upper[:-1]):
        return np.searchsorted(masses, lower), np.searchsorted(masses, upper, side="right")
    left = np.empty(len(lower), dtype=np.int64)
    right = np.empty(len(upper), dtype=np.int64)
    start = stop = 0
    for index in range(len(lower)):
        while start < len(masses) and masses[start] < lower[index]:
            start += 1
        while stop < len(masses) and masses[stop] <= upper[index]:
            stop += 1
        left[index], right[index] = start, stop
    return left, right


@njit(cache=True, nogil=True)
def trace_correlation(a, b):
    """Pair shared nonzero scans by actual time; constant/short traces give no evidence."""
    i = j = n = 0
    mean_a = mean_b = variance_a = variance_b = covariance = 0.0
    while i < len(a) and j < len(b):
        if a[i, 1] < b[j, 1]:
            i += 1
        elif a[i, 1] > b[j, 1]:
            j += 1
        else:
            x, y = a[i, 0], b[j, 0]
            if x > 0 and y > 0:
                n += 1
                dx, dy = x - mean_a, y - mean_b
                mean_a += dx / n
                mean_b += dy / n
                variance_a += dx * (x - mean_a)
                variance_b += dy * (y - mean_b)
                covariance += dx * (y - mean_b)
            i += 1
            j += 1
    if n < 3 or variance_a <= 0 or variance_b <= 0:
        return n, 0.0
    return n, min(1.0, max(-1.0, covariance / np.sqrt(variance_a * variance_b)))


@njit(cache=True, nogil=True)
def _candidate_pairs(ids, source_ids, left, right, expected, apices, counts, rt_tolerance):
    """Reject incompatible mass-window hits before crossing back into Python.

    Preserve original source/target order, including inclusive RT boundaries.
    Dense chromatograms commonly have thousands of equal-mass peaks at other
    times; only coeluting, nonempty traces need correlation or evidence objects.
    """
    # Yield matches incrementally: dense files must not allocate a full pair list.
    for position in range(len(source_ids)):
        i = source_ids[position]
        if expected[position] <= 0 or counts[i] == 0:
            continue
        for target in range(left[position], right[position]):
            j = ids[target]
            if i != j and counts[j] > 0 and abs(apices[i] - apices[j]) <= rt_tolerance:
                yield position, j


def group_ions(features, records, config, stats, configs=None):
    """Annotate candidate adducts, water losses and MS2-supported source fragments.

    Every hypothesis requires matched mass, polarity, apex RT, shared scans,
    overlap and Pearson correlation. Ambiguous ion forms remain separate
    hypotheses. MS2 support is evidence of plausibility, not proof of identity.
    """
    stats.update(adduct_relations=0, neutral_loss_relations=0, isf_relations=0)
    if not features:
        return
    if any(f.trace is None for f in features):
        raise ValueError("Ion grouping requires traces retained with group_ions=True")
    rules = selected_adducts(config.adduct_ids, config.custom_adducts)
    by_name = {rule.name: rule for rule in rules}
    relations = [[] for _ in features]
    masses = np.array([f.mz for f in features])
    apices = np.array([f.rt_apex_seconds for f in features])
    counts = np.array([f.scan_count for f in features], dtype=np.int64)
    rt_tolerances = np.array(
        [config_for(f.polarity, config, configs).group_rt_tolerance for f in features]
    )
    charges = [
        {abs(records[r].precursor_charge) for r in f._ms2_records if records[r].precursor_charge}
        for f in features
    ]
    cache = {}

    def evidence(i, j):
        if (
            i == j
            or counts[i] == 0
            or counts[j] == 0
            or abs(apices[i] - apices[j]) > rt_tolerances[i]
        ):
            return None
        key = min(i, j), max(i, j)
        if key not in cache:
            count, correlation = trace_correlation(features[i].trace, features[j].trace)
            overlap = count / max(features[i].scan_count, features[j].scan_count)
            cache[key] = (
                (int(count), float(correlation))
                if (
                    count >= config.group_min_scans
                    and overlap >= config.group_min_overlap
                    and correlation >= config.group_min_correlation
                )
                else None
            )
        return cache[key]

    def charge_ok(index, rule):
        return not charges[index] or abs(rule.charge) in charges[index]

    for polarity in ("positive", "negative"):
        ids = np.array([i for i, f in enumerate(features) if f.polarity == polarity], dtype=int)
        ids = ids[np.argsort(masses[ids], kind="stable")]
        mz = masses[ids]
        if not len(ids):
            continue
        rt_tolerance = config_for(polarity, config, configs).group_rt_tolerance

        def candidates(expected):
            tolerances = np.maximum(
                config.group_mz_tolerance_da, expected * config.group_mz_tolerance_ppm * 1e-6
            )
            left, right = sorted_mass_windows(mz, expected - tolerances, expected + tolerances)
            for position, target in _candidate_pairs(
                ids, ids, left, right, expected, apices, counts, rt_tolerance
            ):
                i = int(ids[position])
                j = int(target)
                proof = evidence(i, j)
                if proof is not None:
                    yield i, j, float(expected[position]), proof

        ordinary = [r for r in rules if r.polarity == polarity and r.kind == "adduct"]
        for a, b in combinations(ordinary, 2):
            neutral = a.neutral_mass(mz)
            expected = b.mz(neutral)
            expected[neutral <= 0] = -1
            for i, j, target, (shared, correlation) in candidates(expected):
                if not charge_ok(i, a) or not charge_ok(j, b):
                    continue
                proof = {
                    "relation": "adduct",
                    "status": "candidate",
                    "neutral_mass": float(
                        (a.neutral_mass(masses[i]) + b.neutral_mass(masses[j])) / 2
                    ),
                    "mass_error_ppm": float((masses[j] - target) / target * 1e6),
                    "correlation": correlation,
                    "shared_scans": shared,
                }
                relations[i].append(
                    {**proof, "partner_feature_id": j + 1, "ion": a.name, "partner_ion": b.name}
                )
                relations[j].append(
                    {**proof, "partner_feature_id": i + 1, "ion": b.name, "partner_ion": a.name}
                )
                stats["adduct_relations"] += 1

        if not config.group_isf:
            continue
        # Use any associated, singly charged/unknown-charge centroid MS2 as
        # evidence. Best-MS2 ranking remains independent of this search.
        supported = {}
        for i in ids:
            i = int(i)
            for record_index in features[i]._ms2_records:
                scan = records[record_index]
                if (
                    scan.mz is None
                    or not len(scan.mz)
                    or (scan.precursor_charge and abs(scan.precursor_charge) != 1)
                ):
                    continue
                base = float(scan.intensity.max())
                if base <= 0:
                    continue
                relative = scan.intensity / base
                tolerances = np.maximum(
                    config.group_mz_tolerance_da, scan.mz * config.group_mz_tolerance_ppm * 1e-6
                )
                selected = (relative >= config.isf_min_relative_intensity) & (
                    scan.mz < masses[i] - tolerances
                )
                products, relative = scan.mz[selected], relative[selected]
                tolerances = tolerances[selected]
                left = np.searchsorted(mz, products - tolerances)
                right = np.searchsorted(mz, products + tolerances, side="right")
                sources = np.full(len(products), i, dtype=np.int64)
                for k, j in _candidate_pairs(
                    ids, sources, left, right, products, apices, counts, rt_tolerance
                ):
                    product, fraction = float(products[k]), float(relative[k])
                    j = int(j)
                    if evidence(i, j) is None:
                        continue
                    rank = (-fraction, abs(masses[j] - product), scan.spectrum_index)
                    if (i, j) not in supported or rank < supported[i, j][0]:
                        supported[i, j] = (
                            rank,
                            {
                                "evidence_ms2_scan_number": scan.scan_number,
                                "evidence_ms2_native_id": scan.native_id,
                                "evidence_ms2_spectrum_index": scan.spectrum_index,
                                "fragment_mz": float(product),
                                "fragment_relative_intensity": fraction,
                            },
                        )

        def add_fragment(i, j, relation, target, source_ion=None, fragment_ion=None):
            shared, correlation = evidence(i, j)
            proof = {
                "relation": relation,
                "status": "candidate",
                "parent_feature_id": i + 1,
                "fragment_feature_id": j + 1,
                "correlation": correlation,
                "shared_scans": shared,
                "mass_error_ppm": float((masses[j] - target) / target * 1e6),
                "parent_ion": source_ion,
                "fragment_ion": fragment_ion,
                "ms2_supported": (i, j) in supported,
            }
            if (i, j) in supported:
                proof.update(supported[i, j][1])
            relations[i].append({**proof, "role": "parent", "partner_feature_id": j + 1})
            relations[j].append({**proof, "role": "fragment", "partner_feature_id": i + 1})

        loss_pairs = set()
        for loss in (r for r in rules if r.polarity == polarity and r.kind == "neutral_loss"):
            parent = by_name[loss.parent]
            expected = loss.mz(parent.neutral_mass(mz))
            for i, j, target, _ in candidates(expected):
                if not charge_ok(i, parent) or not charge_ok(j, loss):
                    continue
                add_fragment(i, j, "neutral_loss", target, parent.name, loss.name)
                loss_pairs.add((i, j))
                stats["neutral_loss_relations"] += 1
        for (i, j), (_, proof) in sorted(supported.items()):
            if (i, j) not in loss_pairs:
                add_fragment(i, j, "in_source_fragment", proof["fragment_mz"])
                stats["isf_relations"] += 1

    for feature, links in zip(features, relations):
        feature.ion_relations = tuple(links)
        feature.trace = None  # release optional trace histories once annotated
