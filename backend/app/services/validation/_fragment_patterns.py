"""Vendored solvent/salt fragment patterns.

These are the exact ``molvs.fragment.REMOVE_FRAGMENTS`` definitions (name +
SMARTS), vendored here so the deep-composition validation checks keep their
precise fragment-classification behaviour while dropping the unmaintained
``molvs`` dependency. Each pattern's compiled SMARTS reproduces MolVS's
heavy-atom counts exactly (verified by ``test_fragment_patterns.py``).
"""

from __future__ import annotations

from dataclasses import dataclass

from rdkit import Chem

# (name, SMARTS) — verbatim from molvs.fragment.REMOVE_FRAGMENTS.
_RAW_PATTERNS: list[tuple[str, str]] = [
    ("hydrogen", "[H]"),
    ("fluorine", "[F]"),
    ("chlorine", "[Cl]"),
    ("bromine", "[Br]"),
    ("iodine", "[I]"),
    ("lithium", "[Li]"),
    ("sodium", "[Na]"),
    ("potassium", "[K]"),
    ("calcium", "[Ca]"),
    ("magnesium", "[Mg]"),
    ("aluminium", "[Al]"),
    ("barium", "[Ba]"),
    ("bismuth", "[Bi]"),
    ("silver", "[Ag]"),
    ("strontium", "[Sr]"),
    ("zinc", "[Zn]"),
    ("ammonia/ammonium", "[#7]"),
    ("water/hydroxide", "[#8]"),
    ("methyl amine", "[#6]-[#7]"),
    ("sulfide", "S"),
    ("nitrate", "[#7](=[#8])(-[#8])-[#8]"),
    ("phosphate", "[P](=[#8])(-[#8])(-[#8])-[#8]"),
    ("hexafluorophosphate", "[P](-[#9])(-[#9])(-[#9])(-[#9])(-[#9])-[#9]"),
    ("sulfate", "[S](=[#8])(=[#8])(-[#8])-[#8]"),
    ("methyl sulfonate", "[#6]-[S](=[#8])(=[#8])(-[#8])"),
    ("trifluoromethanesulfonic acid", "[#8]-[S](=[#8])(=[#8])-[#6](-[#9])(-[#9])-[#9]"),
    ("trifluoroacetic acid", "[#9]-[#6](-[#9])(-[#9])-[#6](=[#8])-[#8]"),
    ("1,2-dichloroethane", "[Cl]-[#6]-[#6]-[Cl]"),
    ("1,2-dimethoxyethane", "[#6]-[#8]-[#6]-[#6]-[#8]-[#6]"),
    ("1,4-dioxane", "[#6]-1-[#6]-[#8]-[#6]-[#6]-[#8]-1"),
    ("1-methyl-2-pyrrolidinone", "[#6]-[#7]-1-[#6]-[#6]-[#6]-[#6]-1=[#8]"),
    ("2-butanone", "[#6]-[#6]-[#6](-[#6])=[#8]"),
    ("acetate/acetic acid", "[#8]-[#6](-[#6])=[#8]"),
    ("acetone", "[#6]-[#6](-[#6])=[#8]"),
    ("acetonitrile", "[#6]-[#6]#[N]"),
    ("benzene", "[#6]1[#6][#6][#6][#6][#6]1"),
    ("butanol", "[#8]-[#6]-[#6]-[#6]-[#6]"),
    ("t-butanol", "[#8]-[#6](-[#6])(-[#6])-[#6]"),
    ("chloroform", "[Cl]-[#6](-[Cl])-[Cl]"),
    ("cycloheptane", "[#6]-1-[#6]-[#6]-[#6]-[#6]-[#6]-[#6]-1"),
    ("cyclohexane", "[#6]-1-[#6]-[#6]-[#6]-[#6]-[#6]-1"),
    ("dichloromethane", "[Cl]-[#6]-[Cl]"),
    ("diethyl ether", "[#6]-[#6]-[#8]-[#6]-[#6]"),
    ("diisopropyl ether", "[#6]-[#6](-[#6])-[#8]-[#6](-[#6])-[#6]"),
    ("dimethyl formamide", "[#6]-[#7](-[#6])-[#6]=[#8]"),
    ("dimethyl sulfoxide", "[#6]-[S](-[#6])=[#8]"),
    ("ethanol", "[#8]-[#6]-[#6]"),
    ("ethyl acetate", "[#6]-[#6]-[#8]-[#6](-[#6])=[#8]"),
    ("formic acid", "[#8]-[#6]=[#8]"),
    ("heptane", "[#6]-[#6]-[#6]-[#6]-[#6]-[#6]-[#6]"),
    ("hexane", "[#6]-[#6]-[#6]-[#6]-[#6]-[#6]"),
    ("isopropanol", "[#8]-[#6](-[#6])-[#6]"),
    ("methanol", "[#8]-[#6]"),
    ("N,N-dimethylacetamide", "[#6]-[#7](-[#6])-[#6](-[#6])=[#8]"),
    ("pentane", "[#6]-[#6]-[#6]-[#6]-[#6]"),
    ("propanol", "[#8]-[#6]-[#6]-[#6]"),
    ("pyridine", "[#6]-1=[#6]-[#6]=[#7]-[#6]=[#6]-1"),
    ("t-butyl methyl ether", "[#6]-[#8]-[#6](-[#6])(-[#6])-[#6]"),
    ("tetrahydrofurane", "[#6]-1-[#6]-[#6]-[#8]-[#6]-1"),
    ("toluene", "[#6]-[#6]~1~[#6]~[#6]~[#6]~[#6]~[#6]~1"),
    ("xylene", "[#6]-[#6]~1~[#6](-[#6])~[#6]~[#6]~[#6]~[#6]~1"),
]


@dataclass(frozen=True)
class FragmentPattern:
    """A named solvent/salt fragment pattern (drop-in for molvs FragmentPattern)."""

    name: str
    smarts_str: str
    smarts: Chem.Mol | None  # precompiled SMARTS query mol (None if it fails to parse)


def _build() -> list[FragmentPattern]:
    return [
        FragmentPattern(name=name, smarts_str=smarts, smarts=Chem.MolFromSmarts(smarts))
        for name, smarts in _RAW_PATTERNS
    ]


# Precompiled once at import (matches MolVS behaviour).
REMOVE_FRAGMENTS: list[FragmentPattern] = _build()
