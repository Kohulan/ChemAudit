"""
Tests for salt/counterion inventory scoring module.
"""

import pytest
from rdkit import Chem

from app.services.scoring.salt_inventory import (
    SaltInventoryScorer,
    calculate_salt_inventory,
)


class TestSaltInventory:
    """Tests for salt/counterion inventory."""

    def test_single_component_no_salts(self):
        """Single-component molecule should have no salts."""
        mol = Chem.MolFromSmiles("CCO")
        result = calculate_salt_inventory(mol)

        assert not result.has_salts
        assert result.total_fragments == 1
        assert len(result.fragments) == 0

    def test_sodium_salt(self):
        """Sodium salt should be detected."""
        mol = Chem.MolFromSmiles("CC(=O)[O-].[Na+]")
        result = calculate_salt_inventory(mol)

        assert result.has_salts
        assert result.total_fragments == 2

        # Find sodium fragment
        na_fragments = [f for f in result.fragments if f.category == "counterion"]
        assert len(na_fragments) >= 1
        assert any("sodium" in f.name.lower() for f in na_fragments)

    def test_hcl_salt(self):
        """HCl salt should be detected."""
        mol = Chem.MolFromSmiles("CCN.Cl")
        result = calculate_salt_inventory(mol)

        assert result.has_salts
        assert result.total_fragments == 2

        # Find HCl fragment
        salt_fragments = [f for f in result.fragments if f.category != "drug"]
        assert len(salt_fragments) >= 1

    def test_multi_salt(self):
        """Molecule with 2 counterions should detect both."""
        mol = Chem.MolFromSmiles("CC(=O)[O-].[Na+].[Cl-]")
        result = calculate_salt_inventory(mol)

        assert result.has_salts
        assert result.total_fragments == 3

        non_drug = [f for f in result.fragments if f.category != "drug"]
        assert len(non_drug) >= 2

    def test_provenance_fragments(self):
        """Provenance fragments should be used when provided."""
        mol = Chem.MolFromSmiles("CC(=O)O")
        result = calculate_salt_inventory(
            mol, provenance_fragments=["CC(=O)O", "[Na+]"]
        )

        assert result.has_salts
        assert result.total_fragments == 2

    def test_parent_identification(self):
        """Largest fragment should be identified as parent."""
        mol = Chem.MolFromSmiles("CC(=O)Oc1ccccc1C(=O)[O-].[Na+]")
        result = calculate_salt_inventory(mol)

        # Parent should be the larger fragment
        drug_fragments = [f for f in result.fragments if f.category == "drug"]
        assert len(drug_fragments) == 1
        assert drug_fragments[0].heavy_atom_count > 1

    def test_interpretation_present(self):
        """Interpretation should be non-empty."""
        mol = Chem.MolFromSmiles("CC(=O)[O-].[Na+]")
        result = calculate_salt_inventory(mol)

        assert result.interpretation
        assert len(result.interpretation) > 10

    def test_parent_smiles_populated(self):
        """Parent SMILES should be the largest fragment."""
        mol = Chem.MolFromSmiles("CC(=O)Oc1ccccc1C(=O)[O-].[Na+]")
        result = calculate_salt_inventory(mol)

        assert result.parent_smiles
        assert "Na" not in result.parent_smiles

    def test_fragment_mw_populated(self):
        """Each fragment should have MW > 0."""
        mol = Chem.MolFromSmiles("CC(=O)[O-].[Na+]")
        result = calculate_salt_inventory(mol)

        for frag in result.fragments:
            assert frag.mw > 0

    def test_scorer_instantiation(self):
        """SaltInventoryScorer should instantiate."""
        scorer = SaltInventoryScorer()
        assert scorer is not None
