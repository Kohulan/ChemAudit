"""Vendored fragment patterns must match MolVS exactly and classify correctly."""

import importlib.util

import pytest
from rdkit import Chem

from app.services.validation._fragment_patterns import (
    REMOVE_FRAGMENTS,
    FragmentPattern,
)

_HAS_MOLVS = importlib.util.find_spec("molvs") is not None


def test_all_patterns_compile():
    assert len(REMOVE_FRAGMENTS) == 61
    for fp in REMOVE_FRAGMENTS:
        assert isinstance(fp, FragmentPattern)
        assert fp.smarts is not None, f"{fp.name} ({fp.smarts_str}) failed to compile"


@pytest.mark.skipif(not _HAS_MOLVS, reason="molvs not installed (dependency removed)")
def test_vendored_patterns_match_molvs_exactly():
    """While molvs is available, prove the vendored list is byte-for-byte faithful."""
    from molvs.fragment import REMOVE_FRAGMENTS as MOLVS_FRAGMENTS

    vendored = [(fp.name, fp.smarts_str) for fp in REMOVE_FRAGMENTS]
    upstream = [(fp.name, fp.smarts_str) for fp in MOLVS_FRAGMENTS]
    assert vendored == upstream

    # Compiled heavy-atom counts (relied on by the matching logic) must agree.
    for v, u in zip(REMOVE_FRAGMENTS, MOLVS_FRAGMENTS):
        vc = v.smarts.GetNumHeavyAtoms() if v.smarts is not None else None
        uc = u.smarts.GetNumHeavyAtoms() if u.smarts is not None else None
        assert vc == uc, f"{v.name}: heavy-atom count {vc} != molvs {uc}"


def test_patterns_match_expected_fragments():
    """Spot-check that representative solvents/salts match their named patterns."""
    by_name = {fp.name: fp for fp in REMOVE_FRAGMENTS}

    water = Chem.MolFromSmiles("O")
    assert water.HasSubstructMatch(by_name["water/hydroxide"].smarts)

    sodium = Chem.MolFromSmiles("[Na+]")
    assert sodium.HasSubstructMatch(by_name["sodium"].smarts)

    methanol = Chem.MolFromSmiles("CO")
    assert methanol.HasSubstructMatch(by_name["methanol"].smarts)
