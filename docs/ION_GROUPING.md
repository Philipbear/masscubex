# Ion relationships and the upstream audit

MassCubeX retains **peak height** as the default quantitative value and
**summed MS2 peak intensity** as the default representative-spectrum criterion.
Both can be changed independently. Adduct and source-fragment relationships
run by default (`--no-group-ions` disables them) and do not change feature intensities or remove rows.

## What upstream MassCube did

The audit used the cloned commit
[`356748ede691862b9a87f67812db5e728a703b9d`](https://github.com/Philipbear/masscube/tree/356748ede691862b9a87f67812db5e728a703b9d).
These observations describe that commit, not necessarily later upstream versions.

| Area | Audited behavior | MassCubeX behavior |
| --- | --- | --- |
| Quantification | `Params.quant_method = "peak_height"` | Same default; `--quant-method peak_area` selects area. Both columns always survive. |
| Best MS2 | `find_best_ms2` maximizes the sum of its MS2 signal intensities | Same ranking criterion by default; records all associated identifiers and the best scan's score/method. Ties use apex proximity then acquisition order. |
| Adduct candidates | Hardcoded positive/negative dictionaries; descending feature height; an unassigned anchor is assumed to be protonated/deprotonated | User-selected catalog, optional custom ions, and pairwise hypotheses without forcing an anchor identity |
| Negative/charged masses | Enumeration uses signed charge in neutral-mass recovery and in the m/z denominator | Signed charge describes polarity; absolute charge sets the m/z divisor. Ionic mass shifts include electron correction. |
| Chromatographic evidence | Fewer than three shared positive samples return correlation 1.0 | At least five shared samples by default, minimum overlap and correlation; constant or undersupported traces provide no evidence |
| Source fragments | Candidate masses come from the representative parent MS2 spectrum, then mass/RT/correlation checks | Evidence can come from any associated centroid MS2; the supporting scan is recorded independently of the representative scan |
| Group output | Greedy feature-group membership and assigned ion labels | Symmetric candidate links, with the alternatives and evidence preserved |

Source locations: [parameters](https://github.com/Philipbear/masscube/blob/356748ede691862b9a87f67812db5e728a703b9d/src/masscube/params.py),
[MS2 selection](https://github.com/Philipbear/masscube/blob/356748ede691862b9a87f67812db5e728a703b9d/src/masscube/raw_data_utils/core.py),
[grouping and correlation](https://github.com/Philipbear/masscube/blob/356748ede691862b9a87f67812db5e728a703b9d/src/masscube/feature_grouping.py),
and [adduct definitions](https://github.com/Philipbear/masscube/blob/356748ede691862b9a87f67812db5e728a703b9d/src/masscube/utils_functions.py).
The audited grouping tolerances are 0.01 Da, 0.05 minutes and correlation > 0.9.
The rewrite uses seconds consistently. Upstream isotope heuristics have not
been carried over; isotope classification is not implemented in this fork.

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

## Mass and chromatographic evidence

For neutral mass `M`, molecule count `n`, signed charge `q` and ionic mass
shift `s`, `m/z = (n*M + s) / abs(q)` and
`M = (abs(q)*m/z - s) / n`. For the built-in ions, the shift is the sum of
added/subtracted neutral atomic masses minus `q * electron_mass`.
For example, the proton shift is approximately +1.007276452 Da and the
deprotonation shift is its negative. A negative charge never produces a
negative m/z divisor.

For every enabled pair of ordinary ion forms with the same polarity, a
binary search finds features near the predicted mass. A relationship also
requires apex agreement and shared positive intensity observations at the
same actual MS1 acquisition times. A stable Pearson calculation evaluates
those shared observations. Overlap is `shared_count / max(nonzero_count_A,
nonzero_count_B)`. Opposite or unknown polarities are not grouped.
Known associated precursor-charge magnitudes constrain possible ion forms.
Missing charge metadata does not establish a charge state.

Default controls, all available as CLI flags and JSON settings:

| CLI flag | Default | Meaning |
| --- | ---: | --- |
| `--group-mz-ppm`, `--group-mz-da` | 10, 0.002 | Use the larger tolerance |
| `--group-rt-tolerance` | 3 | Maximum apex difference, seconds |
| `--group-min-scans` | 5 | Minimum shared nonzero MS1 observations |
| `--group-min-correlation` | 0.9 | Minimum Pearson correlation |
| `--group-min-overlap` | 0.5 | Shared observations relative to the longer feature |
| `--isf-min-relative-intensity` | 0.05 | Minimum product intensity / MS2 base peak |
| `--isf-max-peaks` | 100 | Maximum retained product masses per MS2 for evidence |

The last two settings limit fragment-evidence memory and work. They **do not
affect the summed-intensity score used to select the best MS2**: the full
exported spectrum is summed first. Scoring does not apply the MS1 noise filter.
Default scoring includes residual precursor signal if present in the MS2 array.

## Water losses and other source fragments

IDs 22 and 23 predict losses of approximately 18.010565 and 36.021129 Da from
a putative `[M+H]+` feature. Their compatible parent ion must be selected too.
They are labelled `neutral_loss`, not ordinary `adduct` relationships.
Mass and chromatographic evidence suffice to propose a candidate;
`ms2_supported: false` explicitly marks the absence of matching MS2 evidence.
When supported, the relationship records the evidence scan number, native ID,
spectrum index, product m/z and relative intensity.

Other `in_source_fragment` candidates require a matching product m/z from
any associated parent MS2. The current interpretation is limited to parents
with singly charged or unknown precursor charge, explicitly centroided MS2,
and product m/z below parent m/z. Profile MS2 contributes scan references and
ranking scores but not centroid product-mass evidence. Higher-charge fragment
interpretation and DIA deconvolution are outside this implementation.

A three-feature water-loss chain can support several hypotheses, including a
successive 18 Da loss as well as a 36 Da loss. Those alternatives remain visible.
Coelution, correlation and even an MS2 match do not prove a source fragment:
coeluting compounds can share fragment masses. Every link retains
`status: "candidate"`; no row is automatically merged, relabelled as a confirmed
compound, or discarded. False-positive rates are not established by these rules.

Each relationship appears in both participating rows of `ion_relations`, with
`partner_feature_id` pointing to the other row. Fragment links also identify
the parent and fragment feature IDs and each row's role. Sidecar relationship
counts count each proposed link once rather than counting its two appearances.

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
