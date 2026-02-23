"""
Tests for scaffold analysis and R-group decomposition.

Covers: Murcko scaffold grouping, generic scaffold extraction, Shannon entropy
diversity metrics, frequency distribution capping, R-group decomposition, and
edge cases (empty batch, acyclic molecules, invalid SMARTS).
"""

import math

import pytest

from app.services.analytics.scaffold_analysis import (
    compute_rgroup_decomposition,
    compute_scaffold_analysis,
)


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------


def _make_result(smiles: str, index: int) -> dict:
    """Build a minimal batch result dict for testing."""
    return {"smiles": smiles, "index": index}


# ---------------------------------------------------------------------------
# compute_scaffold_analysis tests
# ---------------------------------------------------------------------------


def test_scaffold_groups_benzene_derivatives():
    """Toluene, aniline, phenol all share benzene scaffold c1ccccc1; group count = 3."""
    results = [
        _make_result("Cc1ccccc1", 0),   # toluene
        _make_result("Nc1ccccc1", 1),   # aniline
        _make_result("Oc1ccccc1", 2),   # phenol
    ]
    r = compute_scaffold_analysis(results)

    assert r["unique_scaffold_count"] == 1, "All three should share one scaffold"
    assert len(r["scaffolds"]) == 1
    group = r["scaffolds"][0]
    assert group["count"] == 3
    assert group["scaffold_smiles"] == "c1ccccc1"
    assert set(group["molecule_indices"]) == {0, 1, 2}


def test_scaffold_different_scaffolds():
    """Benzene and pyridine molecules yield two distinct scaffold groups."""
    results = [
        _make_result("c1ccccc1", 0),   # benzene — scaffold = c1ccccc1
        _make_result("c1ccncc1", 1),   # pyridine — scaffold = c1ccncc1
    ]
    r = compute_scaffold_analysis(results)

    assert r["unique_scaffold_count"] == 2
    assert len(r["scaffolds"]) == 2

    scaffold_smiles = {g["scaffold_smiles"] for g in r["scaffolds"]}
    assert "c1ccccc1" in scaffold_smiles
    assert "c1ccncc1" in scaffold_smiles


def test_generic_scaffold():
    """Benzene and pyridine share the same all-carbon generic scaffold (C1CCCCC1)."""
    results = [
        _make_result("c1ccccc1", 0),   # benzene
        _make_result("c1ccncc1", 1),   # pyridine (different heteroatom)
    ]
    r = compute_scaffold_analysis(results)

    generic_smiles = {g["generic_scaffold_smiles"] for g in r["scaffolds"]}
    # Both 6-membered rings should collapse to the same generic scaffold
    assert len(generic_smiles) == 1, (
        f"Expected 1 unique generic scaffold but got: {generic_smiles}"
    )


def test_acyclic_molecules():
    """Ethanol and butane (acyclic) are grouped under the empty scaffold ''."""
    results = [
        _make_result("CCO", 0),    # ethanol
        _make_result("CCCC", 1),   # butane
    ]
    r = compute_scaffold_analysis(results)

    assert r["unique_scaffold_count"] == 1
    group = r["scaffolds"][0]
    assert group["scaffold_smiles"] == "", "Acyclic molecules must use empty string scaffold"
    assert group["count"] == 2
    assert set(group["molecule_indices"]) == {0, 1}


def test_scaffold_shannon_entropy_single():
    """All molecules share the same scaffold — Shannon entropy must be 0."""
    results = [
        _make_result("Cc1ccccc1", 0),
        _make_result("Nc1ccccc1", 1),
        _make_result("Oc1ccccc1", 2),
        _make_result("Fc1ccccc1", 3),
    ]
    r = compute_scaffold_analysis(results)

    assert r["unique_scaffold_count"] == 1
    assert r["shannon_entropy"] == 0.0, (
        f"Single-scaffold entropy must be 0.0, got {r['shannon_entropy']}"
    )


def test_scaffold_shannon_entropy_diverse():
    """N scaffolds with equal frequency (1 each) yield entropy = log2(N)."""
    # 4 distinct single-ring scaffolds, one molecule each
    results = [
        _make_result("c1ccccc1", 0),   # benzene scaffold
        _make_result("c1ccncc1", 1),   # pyridine scaffold
        _make_result("C1CCCC1", 2),    # cyclopentane scaffold
        _make_result("C1CCCCC1", 3),   # cyclohexane scaffold
    ]
    r = compute_scaffold_analysis(results)

    n = r["unique_scaffold_count"]
    assert n == 4

    expected_entropy = math.log2(4)
    assert abs(r["shannon_entropy"] - expected_entropy) < 1e-5, (
        f"Expected entropy {expected_entropy}, got {r['shannon_entropy']}"
    )


def test_scaffold_frequency_distribution_cap():
    """More than 50 distinct scaffolds produce a distribution capped at 50 + 'Other'."""
    # 52 molecules with distinct ring scaffolds
    distinct_scaffold_smiles = [
        # 5-membered aromatics
        "c1ccoc1", "c1ccsc1", "c1cc[nH]c1", "c1cnco1", "c1cncs1",
        "c1cnc[nH]1", "c1cnn[nH]1",
        # 6-membered aromatics
        "c1ccccc1", "c1ccncc1", "c1ccnnc1",
        "c1cncnc1", "c1ccncn1", "c1cnccn1", "c1ccccn1", "c1cnncn1", "c1ncncn1",
        # Saturated rings of varying size and heteroatoms
        "C1CCC1", "C1CCCC1", "C1CCCCC1", "C1CCCCCC1", "C1CCCCCCC1",
        "C1NCC1", "C1CCNC1", "C1CCOC1", "C1CCSC1",
        "C1CCNCC1", "C1CCOCC1", "C1CCSCC1", "C1CNCCN1",
        "C1CNOCC1", "C1CNSCC1", "C1CNSC1", "C1CNOC1",
        # Bicyclics
        "c1ccc2ccccc2c1", "c1ccc2ncccc2c1", "c1ccc2[nH]ccc2c1",
        "c1ccc2[nH]cnc2c1", "c1ccc2occc2c1", "c1ccc2sccc2c1",
        "c1ccc2scnc2c1", "c1ccc2ocnc2c1", "c1ccc2cncnc2c1",
        "c1cncc2ccccc12", "c1nncc2ccccc12",
        "c1ccc2ccc3ccccc3c2c1", "c1ccc2cc3ccccc3cc2c1",
        "c1ccc2nc3ccccc3cc2c1",
        "C1CCC2CCCCC2C1", "C1Cc2ccccc21", "C1CCc2ccccc21",
        "c1ccc2[nH]ncc2c1",
        "c1ccc2c(c1)[nH]c1ccccc12", "c1ccc2c(c1)oc1ccccc12",
        "c1ccc2c(c1)sc1ccccc12",
    ]

    results = [_make_result(smi, i) for i, smi in enumerate(distinct_scaffold_smiles)]
    r = compute_scaffold_analysis(results)

    unique_count = r["unique_scaffold_count"]
    freq_dist = r["frequency_distribution"]

    assert unique_count > 50, (
        f"Expected >50 unique scaffolds for the cap test, got {unique_count}"
    )
    # Capped at 50 labelled entries + "Other"
    non_other = {k: v for k, v in freq_dist.items() if k != "Other"}
    assert len(non_other) == 50, (
        f"Expected exactly 50 non-'Other' entries, got {len(non_other)}"
    )
    assert "Other" in freq_dist, "Expected 'Other' bucket for scaffolds beyond top 50"
    assert freq_dist["Other"] >= 1


def test_scaffold_empty_batch():
    """Empty result list produces empty scaffolds list and entropy = 0."""
    r = compute_scaffold_analysis([])

    assert r["scaffolds"] == []
    assert r["unique_scaffold_count"] == 0
    assert r["shannon_entropy"] == 0.0
    assert r["frequency_distribution"] == {}


# ---------------------------------------------------------------------------
# compute_rgroup_decomposition tests
# ---------------------------------------------------------------------------


def test_rgroup_decomposition_basic():
    """Substituted benzenes decompose with the aromatic core and one R-group each."""
    results = [
        _make_result("Cc1ccccc1", 0),    # methyl substituent
        _make_result("CCc1ccccc1", 1),   # ethyl substituent
        _make_result("CCCc1ccccc1", 2),  # propyl substituent
    ]
    r = compute_rgroup_decomposition(results, core_smarts="[*]c1ccccc1")

    assert r["core_smarts"] == "[*]c1ccccc1"
    assert r["unmatched_count"] == 0
    assert len(r["decomposition"]) == 3

    # Each decomposition entry has molecule_index, core, rgroups
    for entry in r["decomposition"]:
        assert "molecule_index" in entry
        assert "core" in entry
        assert "rgroups" in entry
        assert len(entry["rgroups"]) >= 1, "Expected at least R1 group per molecule"


def test_rgroup_invalid_smarts():
    """Invalid SMARTS returns an error dict without raising an exception."""
    results = [_make_result("Cc1ccccc1", 0)]
    r = compute_rgroup_decomposition(results, core_smarts="[invalid{{")

    assert "error" in r
    assert r["error"] == "Invalid SMARTS pattern"
    assert r["decomposition"] == []
    # unmatched_count equals total submitted
    assert r["unmatched_count"] == len(results)


def test_rgroup_no_matches():
    """Valid SMARTS that matches nothing yields unmatched_count equal to total molecules."""
    results = [
        _make_result("C1CCCCC1", 0),   # cyclohexane — no aromatic ring
        _make_result("CCCC", 1),        # butane — acyclic
    ]
    # This aromatic core pattern won't match aliphatic molecules
    r = compute_rgroup_decomposition(results, core_smarts="c1ccccn1")

    assert r["unmatched_count"] == len(results)
    assert r["decomposition"] == []
    assert "error" not in r


def test_rgroup_none_core_smarts():
    """None/empty core SMARTS returns error dict with all molecules unmatched."""
    results = [_make_result("Cc1ccccc1", 0), _make_result("Nc1ccccc1", 1)]
    r = compute_rgroup_decomposition(results, core_smarts=None)

    assert "error" in r
    assert r["unmatched_count"] == len(results)
    assert r["decomposition"] == []
