"""
Tests for Butina Clustering Analytics Service.

Validates compute_butina_clustering with configurable Tanimoto distance cutoff,
invalid SMILES handling, single molecule edge case, and identical molecules.
"""

from __future__ import annotations

import pytest


def _make_result(smiles: str, index: int) -> dict:
    """Create a minimal result dict for clustering input."""
    return {"smiles": smiles, "index": index, "status": "success"}


class TestButinaClustering:
    """Tests for compute_butina_clustering."""

    def test_butina_basic(self):
        """5 similar molecules cluster together with correct output structure."""
        from app.services.analytics.clustering import compute_butina_clustering

        # 5 similar benzoic acid derivatives
        results = [
            _make_result("c1ccccc1C(=O)O", 0),       # benzoic acid
            _make_result("c1ccc(O)cc1C(=O)O", 1),     # salicylic acid
            _make_result("c1ccc(N)cc1C(=O)O", 2),     # aminobenzoic acid
            _make_result("c1ccc(F)cc1C(=O)O", 3),     # fluorobenzoic acid
            _make_result("c1ccc(Cl)cc1C(=O)O", 4),    # chlorobenzoic acid
        ]
        result = compute_butina_clustering(results, distance_cutoff=0.4)

        assert "clusters" in result
        assert "cluster_count" in result
        assert "singleton_count" in result
        assert "largest_cluster_size" in result
        assert "distance_cutoff" in result
        assert result["distance_cutoff"] == 0.4

        # At cutoff 0.4, these similar molecules should form at least 1 cluster
        assert result["cluster_count"] >= 1
        total_members = sum(c["size"] for c in result["clusters"])
        assert total_members == 5

        # Each cluster has required keys
        for cluster in result["clusters"]:
            assert "cluster_id" in cluster
            assert "member_indices" in cluster
            assert "size" in cluster
            assert "representative_index" in cluster
            assert cluster["size"] == len(cluster["member_indices"])

    def test_butina_cutoff_variation(self):
        """Loose cutoff (0.2) produces fewer clusters than strict cutoff (0.6)."""
        from app.services.analytics.clustering import compute_butina_clustering

        results = [
            _make_result("c1ccccc1C(=O)O", 0),
            _make_result("c1ccc(O)cc1C(=O)O", 1),
            _make_result("c1ccc(N)cc1C(=O)O", 2),
            _make_result("c1ccc(F)cc1C(=O)O", 3),
            _make_result("c1ccc(Cl)cc1C(=O)O", 4),
        ]
        # Loose distance cutoff -> more molecules per cluster -> fewer clusters
        loose = compute_butina_clustering(results, distance_cutoff=0.2)
        # Strict distance cutoff -> fewer molecules per cluster -> more clusters
        strict = compute_butina_clustering(results, distance_cutoff=0.6)

        assert loose["cluster_count"] <= strict["cluster_count"]

    def test_butina_invalid_smiles_skipped(self):
        """Invalid SMILES are excluded from clustering; valid_indices maps back to originals."""
        from app.services.analytics.clustering import compute_butina_clustering

        results = [
            _make_result("c1ccccc1", 0),           # valid
            _make_result("INVALID_SMILES", 1),      # invalid
            _make_result("c1ccc(O)cc1", 2),         # valid
            _make_result("NOT_A_MOLECULE", 3),       # invalid
            _make_result("c1ccc(N)cc1", 4),          # valid
        ]
        result = compute_butina_clustering(results, distance_cutoff=0.4)

        # Only 3 valid molecules should be clustered
        total_members = sum(c["size"] for c in result["clusters"])
        assert total_members == 3

        # All member indices should map to original batch indices (0, 2, 4)
        all_indices = set()
        for cluster in result["clusters"]:
            for idx in cluster["member_indices"]:
                all_indices.add(idx)
        assert all_indices == {0, 2, 4}

    def test_butina_single_molecule(self):
        """Single molecule returns 1 cluster with 1 member, 0 singletons."""
        from app.services.analytics.clustering import compute_butina_clustering

        results = [_make_result("c1ccccc1", 0)]
        result = compute_butina_clustering(results, distance_cutoff=0.35)

        assert result["cluster_count"] == 1
        assert result["singleton_count"] == 0
        assert result["largest_cluster_size"] == 1
        assert result["clusters"][0]["member_indices"] == [0]

    def test_butina_all_identical(self):
        """Identical molecules produce 1 cluster containing all molecules."""
        from app.services.analytics.clustering import compute_butina_clustering

        results = [
            _make_result("c1ccccc1", 0),
            _make_result("c1ccccc1", 1),
            _make_result("c1ccccc1", 2),
            _make_result("c1ccccc1", 3),
        ]
        result = compute_butina_clustering(results, distance_cutoff=0.35)

        assert result["cluster_count"] == 1
        assert result["singleton_count"] == 0
        assert result["largest_cluster_size"] == 4
        assert sorted(result["clusters"][0]["member_indices"]) == [0, 1, 2, 3]
