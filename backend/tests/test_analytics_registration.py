"""
Tests for Registration Hash Analytics Service.

Validates compute_registration_hashes with basic hash computation,
tautomer_hash_v2 flag, RDKit version tracking, stereoisomer differentiation,
collision detection, and batch processing.
"""

from __future__ import annotations

from rdkit import rdBase


class TestRegistrationHash:
    """Tests for compute_registration_hashes."""

    def test_basic_hash(self):
        """Ethanol produces a 40-char hex string hash."""
        from app.services.analytics.registration_hash import compute_registration_hashes

        results = [{"smiles": "CCO", "index": 0, "status": "success"}]
        result = compute_registration_hashes(results)

        assert len(result["per_molecule"]) == 1
        mol_hash = result["per_molecule"][0]["hash"]
        assert len(mol_hash) == 40
        assert all(c in "0123456789abcdef" for c in mol_hash)

    def test_tautomer_hash_v2(self):
        """enable_tautomer_hash_v2=True is used, verified by result dict."""
        from app.services.analytics.registration_hash import compute_registration_hashes

        results = [{"smiles": "CCO", "index": 0, "status": "success"}]
        result = compute_registration_hashes(results)

        assert result["tautomer_hash_v2"] is True

    def test_version_stored(self):
        """Result includes rdkit_version field matching rdBase.rdkitVersion."""
        from app.services.analytics.registration_hash import compute_registration_hashes

        results = [{"smiles": "CCO", "index": 0, "status": "success"}]
        result = compute_registration_hashes(results)

        assert "rdkit_version" in result
        assert result["rdkit_version"] == rdBase.rdkitVersion

    def test_meso_tartaric(self):
        """Meso-tartaric acid and L-tartaric acid produce different registration hashes (D-16)."""
        from app.services.analytics.registration_hash import compute_registration_hashes

        results = [
            {
                "smiles": "[C@@H](O)(C(=O)O)[C@H](O)C(=O)O",  # meso-tartaric
                "index": 0,
                "status": "success",
            },
            {
                "smiles": "[C@@H](O)(C(=O)O)[C@@H](O)C(=O)O",  # L-tartaric
                "index": 1,
                "status": "success",
            },
        ]
        result = compute_registration_hashes(results)

        hashes = [m["hash"] for m in result["per_molecule"]]
        assert len(hashes) == 2
        assert hashes[0] != hashes[1], (
            "Meso and L-tartaric acid should have different registration hashes"
        )

    def test_collisions(self):
        """Batch with tautomers may produce collision groups."""
        from app.services.analytics.registration_hash import compute_registration_hashes

        # Keto and enol tautomers should produce the same hash with tautomer_hash_v2
        results = [
            {"smiles": "CC(=O)CC", "index": 0, "status": "success"},      # 2-butanone (keto)
            {"smiles": "C/C(O)=C/C", "index": 1, "status": "success"},    # 2-buten-2-ol (enol)
            {"smiles": "c1ccccc1", "index": 2, "status": "success"},       # benzene (different)
        ]
        result = compute_registration_hashes(results)

        assert "collision_groups" in result
        assert isinstance(result["collision_groups"], list)
        # collision_groups should contain groups with count > 1
        for group in result["collision_groups"]:
            assert "hash" in group
            assert "molecule_indices" in group
            assert "count" in group
            assert group["count"] > 1

    def test_batch_registration(self):
        """Batch of 5 SMILES returns correct structure with all expected fields."""
        from app.services.analytics.registration_hash import compute_registration_hashes

        results = [
            {"smiles": "CCO", "index": 0, "status": "success"},
            {"smiles": "c1ccccc1", "index": 1, "status": "success"},
            {"smiles": "CC(=O)O", "index": 2, "status": "success"},
            {"smiles": "CC(=O)Oc1ccccc1C(=O)O", "index": 3, "status": "success"},
            {"smiles": "c1ccc(O)cc1", "index": 4, "status": "success"},
        ]
        result = compute_registration_hashes(results)

        # Top-level keys
        assert "per_molecule" in result
        assert "unique_count" in result
        assert "total_count" in result
        assert "collision_groups" in result
        assert "rdkit_version" in result
        assert "tautomer_hash_v2" in result

        assert result["total_count"] == 5
        assert result["unique_count"] == 5  # all different molecules
        assert result["tautomer_hash_v2"] is True
        assert len(result["per_molecule"]) == 5

        # Check per_molecule structure
        for entry in result["per_molecule"]:
            assert "index" in entry
            assert "smiles" in entry
            assert "hash" in entry
            assert "canonical_smiles" in entry
            assert "formula" in entry
            assert len(entry["hash"]) == 40
