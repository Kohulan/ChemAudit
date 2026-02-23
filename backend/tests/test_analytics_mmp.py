"""
Tests for Matched Molecular Pair (MMP), Activity Cliff, and LLE Analytics

Covers:
- MMP detection via BRICS fragmentation
- Pair deduplication and pair limit cap
- Batch size enforcement
- Activity cliff SALI computation
- Lipophilic Ligand Efficiency (LLE) computation
- Edge cases (empty batch, missing activity data)
"""

import pytest
from rdkit import Chem
from rdkit.Chem import Descriptors

from app.services.analytics.mmp import (
    MAX_MMP_BATCH_SIZE,
    MAX_PAIRS_RETURNED,
    _compute_activity_cliffs,
    _detect_mmp_pairs,
    compute_mmp_analysis,
)


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------


def _make_result(smiles: str, index: int, properties: dict | None = None) -> dict:
    """
    Build a minimal batch validation result dict for testing.

    Args:
        smiles: SMILES string for the molecule.
        index: Original batch index (molecule_index).
        properties: Optional extra properties dict (e.g., for activity columns).

    Returns:
        Result dict with status='success' and the supplied fields.
    """
    return {
        "molecule_index": index,
        "status": "success",
        "smiles": smiles,
        "properties": properties or {},
    }


# SMILES that are genuine BRICS-decomposable (aromatic–aliphatic bond) and share
# [16*]c1ccccc1 as the common core fragment.
_PHENYLACETIC_ACID = "c1ccccc1CC(=O)O"  # phenylacetic acid
_PHENYLACETAMIDE = "c1ccccc1CC(=O)N"   # phenylacetamide
_BENZYL_FLUORIDE = "c1ccccc1CF"        # benzyl fluoride


# ---------------------------------------------------------------------------
# MMP pair detection tests
# ---------------------------------------------------------------------------


def test_mmp_finds_pairs():
    """
    Phenylacetic acid and phenylacetamide share [16*]c1ccccc1 as the BRICS core
    and differ in their side chain. The analysis must detect at least one MMP pair.
    """
    results = [
        _make_result(_PHENYLACETIC_ACID, 0),
        _make_result(_PHENYLACETAMIDE, 1),
    ]
    res = compute_mmp_analysis(results)

    assert "pairs" in res
    assert len(res["pairs"]) >= 1, "Expected at least one MMP pair for acid/amide"

    pair = res["pairs"][0]
    # Both molecule indices must be present
    assert {pair["mol_a_index"], pair["mol_b_index"]} == {0, 1}
    # Core must be a valid SMILES containing the benzene ring
    core_mol = Chem.MolFromSmiles(pair["core_smiles"])
    assert core_mol is not None, "core_smiles must be a valid SMILES"
    assert any(
        atom.GetIsAromatic() for atom in core_mol.GetAtoms()
    ), "Core fragment should contain aromatic atoms (benzene ring)"


def test_mmp_no_pairs():
    """
    Structurally unrelated molecules (ethanol, benzene, aspirin) should produce
    no MMP pairs because they share no common BRICS core fragments.
    """
    results = [
        _make_result("CCO", 0),           # ethanol — no BRICS bonds
        _make_result("c1ccccc1", 1),      # benzene — no BRICS bonds
        _make_result("CC(=O)Oc1ccccc1C(=O)O", 2),  # aspirin — BRICS bonds but no shared core with others
    ]
    res = compute_mmp_analysis(results)

    assert "pairs" in res
    assert len(res["pairs"]) == 0, "Unrelated molecules should produce no MMP pairs"


def test_mmp_refuses_large_batch():
    """
    Batches exceeding MAX_MMP_BATCH_SIZE must be refused with a descriptive error,
    not silently processed or crash.
    """
    # Construct batch of MAX_MMP_BATCH_SIZE + 1 results using any valid SMILES
    large_results = [
        _make_result(_PHENYLACETIC_ACID, i) for i in range(MAX_MMP_BATCH_SIZE + 1)
    ]
    res = compute_mmp_analysis(large_results)

    assert res["status"] == "refused"
    assert str(MAX_MMP_BATCH_SIZE) in res["reason"]
    assert res["molecule_count"] == MAX_MMP_BATCH_SIZE + 1


def test_mmp_tanimoto_computed():
    """
    MMP pairs must include a Tanimoto similarity value in the range (0, 1).
    Phenylacetic acid and phenylacetamide are structurally similar but not identical.
    """
    results = [
        _make_result(_PHENYLACETIC_ACID, 0),
        _make_result(_PHENYLACETAMIDE, 1),
    ]
    res = compute_mmp_analysis(results)

    assert len(res["pairs"]) >= 1
    tanimoto = res["pairs"][0]["tanimoto"]
    assert 0 < tanimoto < 1, f"Tanimoto {tanimoto} should be in (0, 1)"


def test_mmp_pair_deduplication():
    """
    Each (mol_a, mol_b) combination must appear at most once in the output, even
    when multiple molecules could share multiple different core fragments.
    """
    # Include three molecules; indices 0 and 2 differ only by molecule_index (same SMILES)
    results = [
        _make_result(_PHENYLACETIC_ACID, 0),
        _make_result(_PHENYLACETAMIDE, 1),
        _make_result(_BENZYL_FLUORIDE, 2),
    ]
    res = compute_mmp_analysis(results)

    pair_keys = [
        (p["mol_a_index"], p["mol_b_index"]) for p in res["pairs"]
    ]
    # All keys must be unique — no duplicates
    assert len(pair_keys) == len(set(pair_keys)), "Each (mol_a, mol_b) pair must appear at most once"


def test_mmp_pairs_limited(monkeypatch):
    """
    When MMP pair detection returns more than MAX_PAIRS_RETURNED pairs,
    the output must be capped at MAX_PAIRS_RETURNED entries.

    The cap is tested by monkeypatching _detect_mmp_pairs to return 1500
    synthetic pairs (well above the 1000-pair limit).
    """
    # Build minimal valid results so compute_mmp_analysis proceeds past all guards
    results = [
        _make_result(_PHENYLACETIC_ACID, 0),
        _make_result(_PHENYLACETAMIDE, 1),
    ]

    # Synthetic pairs: 1500 unique (i, i+1) combinations with fake metadata
    synthetic_pairs = [
        {
            "mol_a_index": i,
            "mol_b_index": i + 1,
            "core_smiles": "[16*]c1ccccc1",
            "rgroup_a": "[8*]CC(=O)O",
            "rgroup_b": "[8*]CC(N)=O",
            "tanimoto": 0.8 - i * 0.0001,  # monotonically decreasing so sort is stable
        }
        for i in range(1500)
    ]

    import app.services.analytics.mmp as mmp_module

    monkeypatch.setattr(mmp_module, "_detect_mmp_pairs", lambda mols, indices: synthetic_pairs)

    res = compute_mmp_analysis(results)

    assert len(res["pairs"]) <= MAX_PAIRS_RETURNED, (
        f"Pairs must be capped at {MAX_PAIRS_RETURNED}, got {len(res['pairs'])}"
    )
    assert len(res["pairs"]) == MAX_PAIRS_RETURNED


# ---------------------------------------------------------------------------
# Activity cliff tests
# ---------------------------------------------------------------------------


def test_activity_cliffs_with_data():
    """
    When activity data is provided, activity_cliffs must be populated with
    SALI values for each MMP pair where both molecules have activity data.
    """
    results = [
        _make_result(_PHENYLACETIC_ACID, 0, properties={"pIC50": "7.5"}),
        _make_result(_PHENYLACETAMIDE, 1, properties={"pIC50": "5.5"}),
    ]
    res = compute_mmp_analysis(results, activity_column="pIC50")

    assert res["activity_cliffs"] is not None, "activity_cliffs should be a list when data is provided"
    assert len(res["activity_cliffs"]) >= 1, "Expected at least one activity cliff entry"

    cliff = res["activity_cliffs"][0]
    assert cliff["sali"] > 0
    assert cliff["activity_diff"] == pytest.approx(2.0, abs=1e-6)
    assert 0 < cliff["tanimoto"] < 1


def test_activity_cliffs_without_column():
    """
    When no activity_column is supplied, activity_cliffs must be None.
    """
    results = [
        _make_result(_PHENYLACETIC_ACID, 0),
        _make_result(_PHENYLACETAMIDE, 1),
    ]
    res = compute_mmp_analysis(results, activity_column=None)

    assert res["activity_cliffs"] is None
    assert res["lle_values"] is None


def test_activity_cliffs_missing_values():
    """
    Activity cliff computation skips pairs where either molecule has no activity
    data in the specified column.
    """
    results = [
        _make_result(_PHENYLACETIC_ACID, 0, properties={"pIC50": "7.5"}),
        _make_result(_PHENYLACETAMIDE, 1, properties={}),  # missing activity
    ]
    res = compute_mmp_analysis(results, activity_column="pIC50")

    # Pairs are found (both BRICS-compatible) but no cliff can be computed —
    # mol_b_index=1 has no activity value
    assert res["activity_cliffs"] is not None
    assert len(res["activity_cliffs"]) == 0, (
        "No cliffs expected when one molecule in a pair lacks activity data"
    )


# ---------------------------------------------------------------------------
# LLE computation tests
# ---------------------------------------------------------------------------


def test_lle_computation():
    """
    LLE = activity - LogP must be computed correctly.

    Toluene (Cc1ccccc1) has RDKit LogP ≈ 1.995. With activity=7.0, the
    expected LLE ≈ 5.005.
    """
    toluene_smiles = "Cc1ccccc1"
    results = [
        _make_result(toluene_smiles, 0, properties={"pIC50": "7.0"}),
    ]
    res = compute_mmp_analysis(results, activity_column="pIC50")

    assert res["lle_values"] is not None
    assert len(res["lle_values"]) == 1

    lle_entry = res["lle_values"][0]
    assert lle_entry["molecule_index"] == 0
    assert lle_entry["activity"] == pytest.approx(7.0, abs=1e-6)

    # Verify LLE = activity - LogP formula using RDKit's own LogP
    mol = Chem.MolFromSmiles(toluene_smiles)
    expected_logp = Descriptors.MolLogP(mol)
    expected_lle = 7.0 - expected_logp

    assert lle_entry["logp"] == pytest.approx(expected_logp, abs=1e-4)
    assert lle_entry["lle"] == pytest.approx(expected_lle, abs=1e-4)
    # LLE should be close to 5.0 (LogP of toluene ≈ 1.995)
    assert lle_entry["lle"] == pytest.approx(5.0, abs=0.1)


def test_lle_no_activity():
    """
    When no activity_column is provided, lle_values must be None.
    """
    results = [
        _make_result("Cc1ccccc1", 0),
    ]
    res = compute_mmp_analysis(results)

    assert res["lle_values"] is None


# ---------------------------------------------------------------------------
# Edge case tests
# ---------------------------------------------------------------------------


def test_mmp_empty_batch():
    """
    An empty results list must return an empty pairs list without crashing.
    """
    res = compute_mmp_analysis([])

    assert res["pairs"] == []
    assert res["activity_cliffs"] is None
    assert res["lle_values"] is None


def test_mmp_activity_column_no_matching_values():
    """
    If the specified activity_column exists but none of the molecules have a
    parseable value under that key, activity_cliffs and lle_values must be None
    and a descriptive note must be included in the result.
    """
    results = [
        _make_result(_PHENYLACETIC_ACID, 0, properties={"pIC50": "not_a_number"}),
        _make_result(_PHENYLACETAMIDE, 1, properties={}),
    ]
    res = compute_mmp_analysis(results, activity_column="pIC50")

    assert res["activity_cliffs"] is None
    assert res["lle_values"] is None
    assert "note" in res
    assert "pIC50" in res["note"]


def test_mmp_activity_cliff_sali_formula():
    """
    SALI = |delta_activity| / (1 - tanimoto) must be computed correctly.

    Uses the internal _compute_activity_cliffs helper with known Tanimoto
    and activity values for a precise arithmetic check.
    """
    # Construct a synthetic MMP pair with known values
    synthetic_pair = {
        "mol_a_index": 0,
        "mol_b_index": 1,
        "core_smiles": "[16*]c1ccccc1",
        "rgroup_a": "[8*]CF",
        "rgroup_b": "[8*]CCl",
        "tanimoto": 0.5,
    }
    activities = {0: 8.0, 1: 5.0}  # delta = 3.0

    cliffs = _compute_activity_cliffs([synthetic_pair], activities)

    assert len(cliffs) == 1
    cliff = cliffs[0]
    # SALI = |8.0 - 5.0| / (1 - 0.5) = 3.0 / 0.5 = 6.0
    assert cliff["sali"] == pytest.approx(6.0, abs=1e-6)
    assert cliff["activity_diff"] == pytest.approx(3.0, abs=1e-6)
    assert cliff["tanimoto"] == pytest.approx(0.5, abs=1e-6)
