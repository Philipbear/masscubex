"""MassCube's strongest-first isotope, adduct and source-fragment grouping.

Adapted from huaxuyu/masscube 356748ede691862b9a87f67812db5e728a703b9d.
Feature IDs stay in extraction order; only the processing order is sorted.
The configurable ion catalog replaces upstream's hard-coded dictionaries.
"""

from math import isfinite
from collections import defaultdict

import numpy as np
from numba import njit

from .adducts import catalog, selected_adducts
from .grouping_queries import anchor_snapshot, snapshot_queries, static_query_index, target_tables
from .measurements import pearson_correlation


FLOAT32 = np.dtype(np.float32)
FLOAT64 = np.dtype(np.float64)


@njit(cache=True, nogil=True)
def _matching_candidates(
    sorted_masses,
    sorted_ids,
    masses,
    lower,
    upper,
    target,
    tolerance,
    subtract_float32,
    compare_float32,
    assigned_epoch,
    epoch,
    snapshot,
    polarities,
    observed,
    apices,
    anchor,
    rt_radius,
    priority,
):
    """Search and filter without temporary mass arrays, using NumPy's dtypes.

    The caller resolves scalar promotion for its installed NumPy version.
    Explicit float32 rounding reproduces subtraction before comparison;
    float64 values are retained when NumPy promotes either operation.
    """
    first = np.searchsorted(sorted_masses, lower)
    stop = np.searchsorted(sorted_masses, upper, side="right")
    if stop <= first:
        return np.empty(0, dtype=np.int64)
    selected = np.empty(stop - first, dtype=np.int64)
    count = 0
    limit = float(np.float32(tolerance)) if compare_float32 else tolerance
    for i in range(first, stop):
        candidate = sorted_ids[i]
        if subtract_float32:
            distance = float(abs(np.float32(masses[candidate] - np.float32(target))))
        else:
            distance = abs(float(masses[candidate]) - target)
        if not distance < limit:
            continue
        eligible = assigned_epoch[candidate] < 0 or (
            snapshot and assigned_epoch[candidate] == epoch
        )
        if (
            eligible
            and polarities[candidate] == polarities[anchor]
            and observed[candidate]
            and abs(apices[candidate] - apices[anchor]) < rt_radius
        ):
            selected[count] = candidate
            count += 1
    selected = selected[:count]
    if count < 2:
        return selected
    return selected[np.argsort(priority[selected])]


@njit(cache=True, nogil=True)
def _eligible_candidates(
    ids, assigned_epoch, epoch, snapshot, polarities, observed, apices, anchor, rt_radius, priority
):
    """Filter the mass-matched slice with upstream eligibility and RT rules.

    Mass arithmetic stays in NumPy. This loop only replaces Python object
    lookups; RT arithmetic is float64 without fastmath. Priority values are
    unique, so their ordering is identical to the stable Python sort.
    """
    selected = np.empty(len(ids), dtype=np.int64)
    count = 0
    for candidate in ids:
        eligible = assigned_epoch[candidate] < 0 or (
            snapshot and assigned_epoch[candidate] == epoch
        )
        if (
            eligible
            and polarities[candidate] == polarities[anchor]
            and observed[candidate]
            and abs(apices[candidate] - apices[anchor]) < rt_radius
        ):
            selected[count] = candidate
            count += 1
    selected = selected[:count]
    if count < 2:
        return selected
    return selected[np.argsort(priority[selected])]


def warmup_matching():
    """Compile both query paths once before a long-lived file worker starts."""
    masses = np.array([100.0], dtype=np.float32)
    ids = np.array([0], dtype=np.int64)
    epochs = np.array([-1], dtype=np.int64)
    polarities = np.array([1], dtype=np.int8)
    observed = np.array([True])
    apices = np.array([0.0])
    _matching_candidates(
        masses.astype(np.float64),
        ids,
        masses,
        99.99,
        100.01,
        100.0,
        0.01,
        True,
        True,
        epochs,
        1,
        True,
        polarities,
        observed,
        apices,
        0,
        0.05,
        ids,
    )
    _eligible_candidates(ids, epochs, 1, True, polarities, observed, apices, 0, 0.05, ids)
    if int(np.__version__.split(".")[0]) == 1:
        static_query_index(
            np.array([100.0]), ids, np.array([100.0], dtype=np.float32),
            np.array([[101.0]]), polarities, observed, apices, 1, 0.05, ids, 0.01, 0.0,
        )
        anchor_snapshot(
            np.array([0], dtype=np.int32), np.array([0, 1], dtype=np.int32), epochs, 0, 1
        )


def trace_correlation(a, b):
    """Upstream shared-positive-scan Pearson rule, paired by acquisition time.

    Upstream returns 1 for fewer than three shared observations. Constant
    traces with >=3 observations return NaN and cannot pass the grouping gate.
    """
    _, ia, ib = np.intersect1d(a[:, 1], b[:, 1], assume_unique=True, return_indices=True)
    x, y = a[ia, 0], b[ib, 0]
    return _positive_correlation(x, y)


def _positive_correlation(x, y):
    """Original Pearson calculation on already aligned acquisition samples."""
    keep = (x > 0) & (y > 0)
    n = int(np.count_nonzero(keep))
    if n < 3:
        return n, 1.0
    with np.errstate(invalid="ignore", divide="ignore"):
        return n, float(pearson_correlation(x[keep], y[keep]))


def isotope_targets(mass):
    """Upstream's M+1..M+3 centers, including its first +1 Da offset."""
    return [(f"M+{i}", mass + 1 + (i - 1) * 1.003355) for i in range(1, 4)]


def group_ions(features, records, config, stats, configs=None, *, spectra=None):
    """Greedily group features, starting from the most intense ungrouped ion.

    Enumerate ordinary/water-loss forms from the anchor's assumed [M+H]+ or
    [M-H]- neutral mass. Only its best MS2 supplies fragment masses. Preserve
    upstream's candidate-mask snapshot within an anchor; a later assignment
    replaces that row's earlier link rather than exporting alternatives.
    """
    stats.update(
        adduct_relations=0,
        neutral_loss_relations=0,
        isf_relations=0,
        isotope_relations=0,
        ion_groups=0,
    )
    if not features:
        return
    if any(f.trace is None for f in features):
        raise ValueError("Ion grouping requires retained feature traces")
    if spectra is None:
        raise ValueError("Ion grouping requires processed MS1 spectra for EIC extraction")
    rules = selected_adducts(config.adduct_ids, config.custom_adducts)
    defaults = {"positive": catalog()[1], "negative": catalog()[5]}
    rules_by_polarity = {
        polarity: [r for r in rules if r.polarity == polarity and r.name != base.name]
        for polarity, base in defaults.items()
    }
    masses = np.array([f.mz for f in features], dtype=np.float32)
    apices = np.array([f.rt_apex_seconds / 60.0 for f in features])
    heights = np.array([f.peak_height or 0.0 for f in features])
    polarity_codes = np.array(
        [{"positive": 1, "negative": -1}.get(f.polarity, 0) for f in features], dtype=np.int8
    )
    observed = np.array([f.scan_count > 0 for f in features])
    order = np.argsort(-heights, kind="stable")
    priority = np.empty(len(order), dtype=np.int64)
    priority[order] = np.arange(len(order))
    sorted_ids = np.argsort(masses, kind="stable")
    # searchsorted promotes a float32 array to each float64 query's dtype.
    # Convert the index once, not the whole feature array for every candidate.
    # This exact float32 -> float64 embedding changes no mass decisions below.
    sorted_masses = masses[sorted_ids].astype(np.float64)
    search_mass = sorted_masses.searchsorted
    # A generation records membership in the anchor's original candidate mask.
    # Children assigned during this anchor remain eligible for later rules;
    # previous anchors and the current anchor itself do not. This is the same
    # snapshot as ~assigned.copy(), without copying every feature per anchor.
    assigned_epoch = np.full(len(features), -1, dtype=np.int64)
    epoch = 0
    # Most ions have no relationships. Allocate lists only when evidence exists.
    links = defaultdict(list)
    parent_of = {}
    index_padding = 2 * np.finfo(np.float32).eps
    rt_radius = config.group_rt_tolerance / 60.0
    tables = target_tables(
        masses, rules_by_polarity, defaults, config, {f.polarity for f in features}
    )
    indices = {}
    for polarity, table in tables.items():
        if (
            np.isfinite(table).all() and np.max(np.abs(table), initial=0) < 1e20
            and config.group_mz_tolerance_ppm < 1e12
        ):
            packed, offsets = static_query_index(
                sorted_masses, sorted_ids, masses, table, polarity_codes, observed,
                apices, 1 if polarity == "positive" else -1, rt_radius, priority,
                config.group_mz_tolerance_da, config.group_mz_tolerance_ppm,
            )
            if len(offsets):
                indices[polarity] = packed, offsets, table.shape[1]
    # NumPy 1.x has value-dependent scalar promotion, while NumPy 2.x treats
    # NumPy scalar dtypes as strong. Ask NumPy rather than assuming either rule.
    constant_comparison = (
        np.result_type(FLOAT32, config.group_mz_tolerance_da)
        if config.group_mz_tolerance_ppm == 0
        else None
    )
    ordinary_scalar_queries = (
        int(np.__version__.split(".")[0]) == 1
        and type(config.group_mz_tolerance_da) is float
        and constant_comparison == FLOAT32
    )

    def matches(target, anchor, snapshot=False):
        tolerance = max(config.group_mz_tolerance_da, target * config.group_mz_tolerance_ppm * 1e-6)
        # Expand only the index preselection; the decision uses upstream's
        # float32 array subtraction and strict tolerance below.
        padding = index_padding * abs(target)
        lower, upper = target - tolerance - padding, target + tolerance + padding
        if ordinary_scalar_queries and type(target) in (float, np.float64) and abs(target) < 1e20:
            # NumPy 1.x safely promotes these bounded scalar values to float32.
            # The default fixed-Da query needs no repeated dtype resolution.
            subtract_dtype = compare_dtype = FLOAT32
        else:
            subtract_dtype = np.result_type(FLOAT32, target)
            compare_dtype = (
                constant_comparison
                if constant_comparison is not None and subtract_dtype == FLOAT32
                else np.result_type(subtract_dtype, tolerance)
            )
        if (
            subtract_dtype in (FLOAT32, FLOAT64)
            and compare_dtype in (FLOAT32, FLOAT64)
            and isfinite(lower)
            and isfinite(upper)
        ):
            return _matching_candidates(
                sorted_masses,
                sorted_ids,
                masses,
                float(lower),
                float(upper),
                float(target),
                float(tolerance),
                subtract_dtype == FLOAT32,
                compare_dtype == FLOAT32,
                assigned_epoch,
                epoch,
                snapshot,
                polarity_codes,
                observed,
                apices,
                anchor,
                rt_radius,
                priority,
            ).tolist()
        # Retain NumPy's behavior for unusual externally supplied scalar dtypes.
        first = search_mass(lower)
        stop = search_mass(upper, side="right")
        if first == stop:
            return []
        ids = sorted_ids[first:stop]
        ids = ids[np.abs(masses[ids] - target) < tolerance]
        if not len(ids):
            return []
        return _eligible_candidates(
            ids,
            assigned_epoch,
            epoch,
            snapshot,
            polarity_codes,
            observed,
            apices,
            anchor,
            rt_radius,
            priority,
        ).tolist()

    def own(child, parent):
        previous = parent_of.get(child)
        if previous is not None:
            links[previous] = [
                link for link in links[previous] if link["partner_feature_id"] != child + 1
            ]
            links[child] = [
                link for link in links[child] if link["partner_feature_id"] != previous + 1
            ]
        parent_of[child] = parent
        assigned_epoch[child] = epoch
        features[child].group_id = features[parent].group_id

    def add_isotope(parent, child, state, target):
        own(child, parent)
        f, p = features[child], features[parent]
        f.adduct, f.is_isotope, f.isotope_state = p.adduct, True, state
        proof = dict(
            relation="isotope",
            status="assigned",
            isotope_state=state,
            mass_error_ppm=float((masses[child] - target) / target * 1e6),
        )
        links[parent].append({**proof, "partner_feature_id": child + 1})
        links[child].append({**proof, "partner_feature_id": parent + 1})
        stats["isotope_relations"] += 1

    for f in features:
        f.adduct = defaults[f.polarity].name if f.polarity in defaults else None
        f.group_id = None
        f.isotope_state = "M+0" if f.polarity in defaults else None
        f.is_isotope = False if f.polarity in defaults else None
        f.is_in_source_fragment = False if config.group_isf and f.polarity in defaults else None

    for anchor in order:
        anchor = int(anchor)
        f = features[anchor]
        if assigned_epoch[anchor] >= 0 or f.polarity not in defaults:
            continue
        assigned_epoch[anchor] = epoch
        epoch += 1
        stats["ion_groups"] += 1
        f.group_id, f.isotope_state = stats["ion_groups"], "M+0"
        if not f.scan_count:
            continue
        base = defaults[f.polarity]
        table = tables.get(f.polarity)
        index = indices.get(f.polarity)
        packed = offsets = None
        if index is not None:
            packed, offsets = anchor_snapshot(*index[:2], assigned_epoch, anchor, index[2])
        neutral = base.neutral_mass(masses[anchor]) if table is None else None
        ions = [
            (i, rule.name, float(rule.mz(neutral)) if table is None else float(table[anchor, i]), None)
            for i, rule in enumerate(rules_by_polarity[f.polarity], 3)
            if index is None or offsets[i] != offsets[i + 1]
        ]
        best = None if not f._ms2_records else next(
            (
                records[i]
                for i in f._ms2_records
                if records[i].spectrum_index == f.best_ms2_spectrum_index
            ),
            None,
        )
        if config.group_isf and best is not None and best.mz is not None:
            fragment_query = 3 + len(rules_by_polarity[f.polarity])
            ions.extend(
                (fragment_query + k, None, float(mass), k) for k, mass in enumerate(best.mz)
            )
        initial_isotopes = isotope_targets(masses[anchor])
        if table is not None and index is None:
            targets = np.array([target for _, target in initial_isotopes] + [ion[2] for ion in ions])
            limits = np.maximum(
                config.group_mz_tolerance_da, targets * config.group_mz_tolerance_ppm * 1e-6
            )
            padding = index_padding * np.abs(targets)
            lower, upper = targets - limits - padding, targets + limits + padding
            if (
                np.isfinite(lower).all() and np.isfinite(upper).all()
                and np.max(np.abs(targets), initial=0) < 1e20
                and np.max(limits, initial=0) < 1e20
            ):
                packed, offsets = snapshot_queries(
                    sorted_masses, sorted_ids, masses, targets, limits, lower, upper,
                    assigned_epoch, polarity_codes, observed, apices, anchor, rt_radius, priority,
                )
                if not len(offsets):
                    offsets = None
        for i, (state, target) in enumerate(initial_isotopes):
            candidates = (
                matches(target, anchor, snapshot=True) if offsets is None
                else packed[offsets[i]:offsets[i + 1]]
            )
            for child in candidates:
                add_isotope(anchor, int(child), state, target)
        times = f.trace[:, 1]
        ordered_times = None
        eic_query = None
        for query, ion, target, product in ions:
            if target <= 0:
                continue
            candidates = (
                matches(target, anchor, snapshot=True)
                if offsets is None or (index is not None and query >= index[2])
                else packed[offsets[query]:offsets[query + 1]]
            )
            if not len(candidates):
                continue
            # Test only the nearest RT candidate. Ties follow descending height.
            child = int(min(candidates, key=lambda j: abs(apices[j] - apices[anchor])))
            if eic_query is None:
                # Reuse the scan-time join for all partners of this anchor.
                prepare = getattr(spectra, "prepare_query", None)
                if prepare is not None:
                    eic_query = prepare(times, f.polarity, config.group_eic_mz_tolerance_da)
                else:
                    def eic_query(mass):
                        return spectra.eic(
                            mass, times, f.polarity, config.group_eic_mz_tolerance_da
                        )
            y = eic_query(masses[child])
            # The EIC is sampled on the anchor's own time grid. For normal,
            # strictly increasing traces the intersection is the identity.
            # Keep the general path for externally supplied unordered traces.
            if ordered_times is None:
                ordered_times = np.all(times[:-1] < times[1:])
            if ordered_times:
                count, correlation = _positive_correlation(f.trace[:, 0], y)
            else:
                count, correlation = trace_correlation(f.trace, np.column_stack((y, times)))
            if not correlation > config.group_min_correlation:
                continue
            own(child, anchor)
            g = features[child]
            g.adduct = ion if ion is not None else f.adduct
            error = float((masses[child] - target) / target * 1e6)
            proof = dict(
                status="assigned",
                correlation=float(correlation),
                shared_scans=int(count),
                mass_error_ppm=error,
            )
            if ion is not None:
                if neutral is None:
                    neutral = base.neutral_mass(masses[anchor])
                proof.update(relation="adduct", neutral_mass=float(neutral))
                links[anchor].append(
                    {**proof, "partner_feature_id": child + 1, "ion": f.adduct, "partner_ion": ion}
                )
                links[child].append(
                    {**proof, "partner_feature_id": anchor + 1, "ion": ion, "partner_ion": f.adduct}
                )
                stats["adduct_relations"] += 1
            else:
                g.is_in_source_fragment = True
                proof.update(
                    relation="in_source_fragment",
                    parent_feature_id=anchor + 1,
                    fragment_feature_id=child + 1,
                    parent_ion=f.adduct,
                    fragment_ion=g.adduct,
                    ms2_supported=True,
                    evidence_ms2_scan_number=best.scan_number,
                    evidence_ms2_native_id=best.native_id,
                    evidence_ms2_spectrum_index=best.spectrum_index,
                    fragment_mz=target,
                    fragment_relative_intensity=float(
                        best.intensity[product] / best.intensity.max()
                    ),
                )
                links[anchor].append({**proof, "role": "parent", "partner_feature_id": child + 1})
                links[child].append({**proof, "role": "fragment", "partner_feature_id": anchor + 1})
                stats["isf_relations"] += 1
            for state, isotope_mass in isotope_targets(masses[child]):
                for isotope in matches(isotope_mass, anchor):
                    if heights[isotope] < heights[child] * config.isotope_rel_int_limit:
                        add_isotope(child, isotope, state, isotope_mass)

    for i, feature in enumerate(features):
        feature.ion_relations = tuple(links.get(i, ()))
        feature.trace = None
    # Count the final relationships, not intermediate overwritten assignments.
    stats.update(adduct_relations=0, isf_relations=0, isotope_relations=0)
    counters = {
        "adduct": "adduct_relations",
        "in_source_fragment": "isf_relations",
        "isotope": "isotope_relations",
    }
    for index, relations in links.items():
        parent = index + 1
        for link in relations:
            if parent < link["partner_feature_id"]:
                stats[counters[link["relation"]]] += 1
