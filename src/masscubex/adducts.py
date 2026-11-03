"""Configurable ion forms, with ionic mass shifts and signed charge metadata.

m/z = (molecules * neutral_mass + mass_shift) / abs(charge).
Mass shifts include the electron correction; negative charge never negates m/z.
"""

from dataclasses import asdict, dataclass
import math
import re

ELECTRON_MASS = 0.000548579909065
ATOMIC_MASS = {
    "H": 1.00782503223,
    "C": 12.0,
    "N": 14.00307400443,
    "O": 15.99491461957,
    "Na": 22.9897692820,
    "K": 38.9637064864,
    "Cl": 34.968852682,
    "Br": 78.9183376,
}


@dataclass(frozen=True)
class Adduct:
    id: int
    name: str
    charge: int
    molecules: int
    mass_shift: float
    kind: str = "adduct"
    parent: str | None = None

    def __post_init__(self):
        if type(self.id) is not int or self.id < 1:
            raise ValueError("Adduct id must be a positive integer")
        if not isinstance(self.name, str) or not self.name or any(c in self.name for c in "\t\r\n"):
            raise ValueError("Adduct name must be a nonempty single-line string")
        if type(self.charge) is not int or self.charge == 0:
            raise ValueError("Adduct charge must be a nonzero signed integer")
        if type(self.molecules) is not int or self.molecules < 1:
            raise ValueError("Adduct molecules must be a positive integer")
        if (
            isinstance(self.mass_shift, bool)
            or not isinstance(self.mass_shift, (int, float))
            or not math.isfinite(self.mass_shift)
        ):
            raise ValueError("Adduct mass_shift must be a finite ionic mass shift")
        if self.kind not in {"adduct", "neutral_loss"}:
            raise ValueError("Adduct kind must be adduct or neutral_loss")
        if self.kind == "neutral_loss" and not self.parent:
            raise ValueError("Neutral-loss forms need a parent ion name")

    @property
    def polarity(self):
        return "positive" if self.charge > 0 else "negative"

    def neutral_mass(self, mz):
        return (abs(self.charge) * mz - self.mass_shift) / self.molecules

    def mz(self, neutral_mass):
        return (self.molecules * neutral_mass + self.mass_shift) / abs(self.charge)


def _ion(id, name, charge, molecules=1, kind="adduct", parent=None, **atoms):
    shift = sum(ATOMIC_MASS[element] * count for element, count in atoms.items())
    return Adduct(id, name, charge, molecules, shift - charge * ELECTRON_MASS, kind, parent)


ADDUCTS = (
    _ion(1, "[M+H]+", 1, H=1),
    _ion(2, "[M+Na]+", 1, Na=1),
    _ion(3, "[M+K]+", 1, K=1),
    _ion(4, "[M+NH4]+", 1, N=1, H=4),
    _ion(5, "[M-H]-", -1, H=-1),
    _ion(6, "[M+Cl]-", -1, Cl=1),
    _ion(7, "[M+Br]-", -1, Br=1),
    _ion(8, "[M+HCOO]-", -1, H=1, C=1, O=2),
    _ion(9, "[M+CH3COO]-", -1, H=3, C=2, O=2),
    _ion(10, "[M+H+CH3CN]+", 1, H=4, C=2, N=1),
    _ion(11, "[M+H+CH3OH]+", 1, H=5, C=1, O=1),
    _ion(12, "[2M+H]+", 1, molecules=2, H=1),
    _ion(13, "[2M+Na]+", 1, molecules=2, Na=1),
    _ion(14, "[2M+NH4]+", 1, molecules=2, N=1, H=4),
    _ion(15, "[2M-H]-", -1, molecules=2, H=-1),
    _ion(16, "[M+2H]2+", 2, H=2),
    _ion(17, "[M+3H]3+", 3, H=3),
    _ion(18, "[M-2H]2-", -2, H=-2),
    _ion(19, "[M-3H]3-", -3, H=-3),
    _ion(20, "[3M+H]+", 1, molecules=3, H=1),
    _ion(21, "[3M-H]-", -1, molecules=3, H=-1),
    _ion(22, "[M+H-H2O]+", 1, kind="neutral_loss", parent="[M+H]+", H=-1, O=-1),
    _ion(23, "[M+H-2H2O]+", 1, kind="neutral_loss", parent="[M+H]+", H=-3, O=-2),
)

# User-selected default: all 1–15, no 16–21, plus the two water-loss hypotheses.
DEFAULT_ADDUCT_IDS = tuple(range(1, 16)) + (22, 23)


def catalog(custom=()):
    result = {ion.id: ion for ion in ADDUCTS}
    names = {ion.name for ion in ADDUCTS}
    for ion in custom:
        if not isinstance(ion, Adduct):
            raise ValueError("custom_adducts must contain Adduct definitions")
        if ion.id in result or ion.name in names:
            raise ValueError(f"Duplicate adduct ID/name: {ion.id} {ion.name}")
        result[ion.id] = ion
        names.add(ion.name)
    return result


def selected_adducts(ids, custom=()):
    available = catalog(custom)
    if not isinstance(ids, tuple) or not ids or any(type(id) is not int for id in ids):
        raise ValueError("adduct_ids must be a nonempty tuple of integer IDs")
    if len(set(ids)) != len(ids) or any(id not in available for id in ids):
        raise ValueError("Unknown or duplicate adduct IDs; see 'masscubex adducts'")
    result = tuple(available[id] for id in ids)
    by_name = {ion.name: ion for ion in result}
    for ion in result:
        if ion.kind == "neutral_loss":
            parent = by_name.get(ion.parent)
            if (
                parent is None
                or parent.kind != "adduct"
                or parent.charge != ion.charge
                or parent.molecules != ion.molecules
                or parent.mass_shift <= ion.mass_shift
            ):
                raise ValueError(
                    f"{ion.name} requires its compatible parent {ion.parent} in adduct_ids"
                )
    return result


def parse_selection(text, custom=()):
    available = catalog(custom)
    by_name = {ion.name: ion.id for ion in available.values()}
    ids = []
    for token in text.split(","):
        token = token.strip()
        if re.fullmatch(r"\d+-\d+", token):
            lo, hi = map(int, token.split("-"))
            if hi < lo or hi - lo > 10000:
                raise ValueError(f"Invalid adduct range: {token}")
            ids.extend(range(lo, hi + 1))
        elif token.isdecimal():
            ids.append(int(token))
        elif token in by_name:
            ids.append(by_name[token])
        else:
            raise ValueError(f"Unknown adduct selection: {token}")
    result = tuple(sorted(set(ids)))
    selected_adducts(result, custom)
    return result


def catalog_rows(custom=(), enabled=DEFAULT_ADDUCT_IDS):
    return [
        {**asdict(ion), "polarity": ion.polarity, "enabled": id in enabled}
        for id, ion in sorted(catalog(custom).items())
    ]
