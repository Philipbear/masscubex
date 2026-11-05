# Ion grouping inherited from MassCube

Grouping follows the strongest-first single-file algorithm in
[huaxuyu/masscube](https://github.com/huaxuyu/masscube), using the selected catalog
below. It is enabled by default and changes labels, not intensities or MS2
ownership. See [the complete algorithm and adaptations](UPSTREAM_COMPARISON.md).

Each unassigned anchor starts as `[M+H]+` or `[M-H]-`. Group isotopes, matching
ion forms, and source fragments from the selected best parent MS2. The nearest
RT candidate is checked against the partner EIC over the anchor scan times.
Upstream defaults are mass error <0.01 Da, RT difference <3 seconds, partner
EIC ±0.005 Da and shared-positive-scan Pearson correlation >0.9. Fewer than
three shared observations returns 1, as upstream; constant longer traces fail.
No five-observation or overlap-fraction requirement is added.

The output has one `adduct` string per feature, `group_id`, `isotope_state`
and `is_in_source_fragment`. The default anchor label is an
assumption, not an identification. Related features have one final assignment;
`--include-ion-relations` exports the assigned links and available evidence in
the optional `ion_relations` column. A monomer
H/Na pair does not generate a second dimer H/Na hypothesis. Actual dimer-mass
ions remain allowed. Isotopes use the upstream M+1..M+3 offsets and the 1.5×
relative-height limit for isotopes of assigned child ions.

The named water-loss forms below are handled as ion forms, as upstream. They
are not independently marked as source fragments; a later best-MS2 match can
set that flag under the upstream assignment order. Other fragments inherit the
parent's ion form. Isotope flags do not automatically inherit the source-fragment
flag. Polarity streams are grouped separately; unknown polarity leaves labels
blank. `--no-group-ions` disables grouping; `--no-isf` disables MS2 product matches.

The best spectrum after upstream fixed-intensity filtering and adjacent-centroid
merging supplies fragment masses by default.
`--isf-min-relative-intensity` defaults to 0 and `--isf-max-peaks` to 0 (unlimited);
explicit pruning remains available. MS2 uses upstream adjacent-centroid merging; a separate profile-MS2
mass-detection algorithm is not implemented. Grouping is not DIA deconvolution.

Source-fragment assignment does not use a catalog of allowed neutral losses.
It searches for an MS1 feature matching an observed product-ion m/z from the
anchor's selected best MS2: default mass error <0.01 Da, apex RT difference
<3 seconds, same polarity, and partner EIC correlation >0.9. The EIC is extracted
within ±0.005 Da on the anchor's scan times, and only the nearest-RT mass
candidate is checked. When requested, the output link records the parent/fragment feature IDs,
supporting MS2 identifiers, product-ion mass, mass error and correlation. These
are putative assignments; coeluting unrelated ions or a chimeric MS2 can also
satisfy the rule. As in upstream, fewer than three shared positive scans returns
correlation 1, so that edge case does not establish a well-sampled peak-shape match.

## Selected ion catalog

IDs are stable. `masscubex adducts` shows exact ionic mass shifts; `--json`
also makes the catalog available to scripts. Enabled means considered during
grouping, which runs by default in 0.5.0. Use `--no-group-ions` to disable it.

| ID | Form | Enabled | Interpretation |
| ---: | --- | :---: | --- |
| 1 | `[M+H]+` | yes | Protonated |
| 2 | `[M+Na]+` | yes | Sodium |
| 3 | `[M+K]+` | yes | Potassium |
| 4 | `[M+NH4]+` | yes | Ammonium |
| 5 | `[M-H]-` | yes | Deprotonated |
| 6 | `[M+Cl]-` | yes | Chloride |
| 7 | `[M+Br]-` | yes | Bromide |
| 8 | `[M+HCOO]-` | yes | Formate |
| 9 | `[M+CH3COO]-` | yes | Acetate |
| 10 | `[M+H+CH3CN]+` | yes | Acetonitrile |
| 11 | `[M+H+CH3OH]+` | yes | Methanol |
| 12 | `[2M+H]+` | yes | Protonated dimer |
| 13 | `[2M+Na]+` | yes | Sodium dimer |
| 14 | `[2M+NH4]+` | yes | Ammonium dimer |
| 15 | `[2M-H]-` | yes | Deprotonated dimer |
| 16 | `[M+2H]2+` | no | Doubly protonated |
| 17 | `[M+3H]3+` | no | Triply protonated |
| 18 | `[M-2H]2-` | no | Doubly deprotonated |
| 19 | `[M-3H]3-` | no | Triply deprotonated |
| 20 | `[3M+H]+` | no | Protonated trimer |
| 21 | `[3M-H]-` | no | Deprotonated trimer |
| 22 | `[M+H-H2O]+` | yes | Neutral loss from `[M+H]+` |
| 23 | `[M+H-2H2O]+` | yes | Two water losses from `[M+H]+` |

## Custom ions and reproducibility

For example, save this JSON list to `ions.json`:

```json
[
  {
    "id": 100,
    "name": "[M+Li]+",
    "charge": 1,
    "molecules": 1,
    "mass_shift": 7.015454857
  }
]
```

```bash
masscubex adducts --adduct-file ions.json --adducts 1-15,22,23,100
masscubex extract sample.mzML -o grouped --group-ions \
  --adduct-file ions.json --adducts 1-15,22,23,100
```

IDs and names must be unique, including against disabled built-ins. `mass_shift`
is the **total ionic shift**, not a per-charge m/z offset. Optional `kind` is
`adduct` or `neutral_loss`; a loss also needs `parent`, naming an enabled ordinary
ion of the same charge and molecule count with a higher mass shift. No formula
parser silently infers custom masses. The exact selected definitions and all
settings are saved in each table's JSON sidecar and participate in resume validation.
