"""
Tests for unified_screen orchestrator.

Covers:
  - Aspirin yields zero Kazius alerts
  - 4-aminobiphenyl triggers Mutagenicity Toxicophore concern group
  - Deduplication: concern_groups deduplicated while raw alerts preserve all
  - total_raw >= total_deduped invariant
  - NIBR screening runs without error
  - Complexity filter is NOT included in unified_screen results
"""

import pytest
from rdkit import Chem

from app.services.alerts.unified_screen import unified_screen

# ---------------------------------------------------------------------------
# Helper
# ---------------------------------------------------------------------------

def _get_mol(smiles: str) -> Chem.Mol:
    mol = Chem.MolFromSmiles(smiles)
    assert mol is not None, f"Could not parse SMILES: {smiles}"
    return mol


# ---------------------------------------------------------------------------
# Tests
# ---------------------------------------------------------------------------

class TestUnifiedScreenAspirin:
    """Aspirin: a clean drug-like molecule with minimal alerts."""

    @pytest.fixture(scope="class")
    def aspirin_result(self):
        mol = _get_mol("CC(=O)Oc1ccccc1C(=O)O")
        return unified_screen(mol)

    def test_aspirin_no_kazius_alerts(self, aspirin_result):
        """Aspirin should produce zero Kazius mutagenicity alerts."""
        result = aspirin_result
        kazius_group = result["concern_groups"].get("Mutagenicity Toxicophore")
        if kazius_group is not None:
            # If the group exists it must have zero alerts from KAZIUS source
            kazius_hits = [
                a for a in kazius_group["alerts"]
                if a.catalog_source == "KAZIUS"
            ]
            assert len(kazius_hits) == 0, (
                f"Aspirin should have 0 Kazius alerts, got {len(kazius_hits)}"
            )

    def test_returns_required_keys(self, aspirin_result):
        """unified_screen must return all required top-level keys."""
        required_keys = {
            "alerts", "concern_groups", "total_raw", "total_deduped",
            "screened_catalogs", "has_critical", "has_warning",
        }
        assert required_keys.issubset(aspirin_result.keys()), (
            f"Missing keys: {required_keys - aspirin_result.keys()}"
        )

    def test_screened_catalogs_includes_all_sources(self, aspirin_result):
        """All four screening sources must appear in screened_catalogs."""
        catalogs = aspirin_result["screened_catalogs"]
        for expected in ["CUSTOM", "KAZIUS", "NIBR"]:
            assert expected in catalogs, f"{expected} missing from screened_catalogs"


class TestUnifiedScreen4Aminobiphenyl:
    """4-Aminobiphenyl: known Kazius mutagenicity hit (aromatic amine)."""

    @pytest.fixture(scope="class")
    def abp_result(self):
        mol = _get_mol("Nc1ccc(-c2ccccc2)cc1")
        return unified_screen(mol)

    def test_4aminobiphenyl_kazius_match(self, abp_result):
        """4-Aminobiphenyl should appear in Mutagenicity Toxicophore group."""
        concern_groups = abp_result["concern_groups"]
        assert "Mutagenicity Toxicophore" in concern_groups, (
            "Mutagenicity Toxicophore group missing for 4-aminobiphenyl"
        )
        group = concern_groups["Mutagenicity Toxicophore"]
        kazius_hits = [a for a in group["alerts"] if a.catalog_source == "KAZIUS"]
        assert len(kazius_hits) > 0, "Should have at least one KAZIUS hit"
        names = [a.pattern_name for a in kazius_hits]
        assert "aromatic_amine" in names, (
            f"aromatic_amine not in Kazius hits: {names}"
        )

    def test_4aminobiphenyl_aromatic_amine_has_atom_indices(self, abp_result):
        """aromatic_amine alert must report matched atom indices."""
        group = abp_result["concern_groups"].get("Mutagenicity Toxicophore", {})
        alerts = group.get("alerts", [])
        aromatic_amine_hits = [
            a for a in alerts
            if a.pattern_name == "aromatic_amine" and a.catalog_source == "KAZIUS"
        ]
        assert len(aromatic_amine_hits) > 0
        assert len(aromatic_amine_hits[0].matched_atoms) > 0, (
            "aromatic_amine must return atom indices"
        )


class TestUnifiedScreenDeduplication:
    """Test that concern-group deduplication works correctly."""

    def test_concern_group_dedup(self):
        """
        Nitrobenzene matches nitro in multiple catalogs; concern_groups
        should show ONE entry per unique (concern_group, pattern) key,
        while raw alerts should show all catalog-source hits.
        """
        # Nitrobenzene — nitro group flagged by multiple catalogs
        mol = _get_mol("O=[N+]([O-])c1ccccc1")
        result = unified_screen(mol)

        # raw_alerts may contain the same concern (nitro) from multiple sources
        # but concern_groups must deduplicate within each group
        for group_name, group_data in result["concern_groups"].items():
            deduped_alerts = group_data["alerts"]
            # Check deduplication: normalized keys must be unique within group
            norm_keys = [
                a.pattern_name.lower().replace(" ", "_")
                for a in deduped_alerts
            ]
            assert len(norm_keys) == len(set(norm_keys)), (
                f"Duplicate pattern_name found in concern_group '{group_name}': "
                f"{[k for k in norm_keys if norm_keys.count(k) > 1]}"
            )

    def test_raw_alerts_preserves_all(self):
        """
        For a molecule with multi-catalog hits, total_raw must be >= total_deduped.
        This verifies the raw list is never trimmed.
        """
        # Nitrobenzene again — hits multiple catalogs
        mol = _get_mol("O=[N+]([O-])c1ccccc1")
        result = unified_screen(mol)
        assert result["total_raw"] >= result["total_deduped"], (
            f"total_raw ({result['total_raw']}) < total_deduped ({result['total_deduped']})"
        )

    def test_total_deduped_equals_sum_of_group_counts(self):
        """total_deduped must equal the sum of all group counts."""
        mol = _get_mol("Nc1ccc(-c2ccccc2)cc1")
        result = unified_screen(mol)
        computed_total = sum(
            g["count"] for g in result["concern_groups"].values()
        )
        assert result["total_deduped"] == computed_total, (
            f"total_deduped {result['total_deduped']} != sum of group counts {computed_total}"
        )


class TestUnifiedScreenNIBR:
    """NIBR screening runs without errors for any valid molecule."""

    def test_nibr_screening_runs(self):
        """Any valid molecule can be screened without error."""
        mol = _get_mol("CC(=O)Oc1ccccc1C(=O)O")  # aspirin
        result = unified_screen(mol)
        # Must return a list (even if empty)
        nibr_raw = [a for a in result["alerts"] if a.catalog_source == "NIBR"]
        assert isinstance(nibr_raw, list), "NIBR alerts must be a list"

    def test_nibr_results_are_alert_result_instances(self):
        """NIBR hits must be AlertResult instances with correct catalog_source."""
        from app.services.alerts.alert_manager import AlertResult
        mol = _get_mol("c1ccc(NC(=O)c2ccco2)cc1")  # simple amine + furanyl
        result = unified_screen(mol)
        nibr_hits = [a for a in result["alerts"] if a.catalog_source == "NIBR"]
        for hit in nibr_hits:
            assert isinstance(hit, AlertResult), (
                f"NIBR hit is not AlertResult: {type(hit)}"
            )
            assert hit.catalog_source == "NIBR"


class TestUnifiedScreenComplexity:
    """Complexity percentile is NOT part of unified_screen output."""

    def test_complexity_not_in_alerts(self):
        """
        unified_screen must NOT call complexity_filter.
        Verify by checking that result keys do not include complexity data.
        """
        mol = _get_mol("CC(=O)Oc1ccccc1C(=O)O")
        result = unified_screen(mol)
        # Top-level keys must not include complexity-related keys
        assert "properties" not in result, (
            "complexity_filter output should not be in unified_screen result"
        )
        assert "n_outliers" not in result, (
            "n_outliers should not be in unified_screen result"
        )
        assert "within_range" not in result, (
            "within_range should not be in unified_screen result"
        )

    def test_complexity_not_in_concern_groups(self):
        """
        Complexity percentile groups (e.g. 'MW outlier') must not appear
        in unified_screen concern_groups — they belong to safety/assess.
        """
        mol = _get_mol("CC[C@@H]1NC(=O)[C@H]([C@H](O)[C@@H](C)CCC=CC)N(C)C(=O)"
                       "[C@H](CC(C)C)N(C)C(=O)[C@H](CC(C)C)N(C)C(=O)[C@H](C)"
                       "NC(=O)C(C)(C)N(C)C(=O)[C@H](CC(C)C)N(C)C(=O)[C@H](CC(C)C)"
                       "N(C)C(=O)[C@@H](C)NC(=O)[C@H](C(C)C)N(C)C1=O")
        result = unified_screen(mol)
        # None of the standard complexity property names should appear as group names
        for group_name in result["concern_groups"]:
            assert group_name not in {"MW", "LogP", "RotBonds", "Rings", "HeavyAtoms", "BertzCT"}, (
                f"Complexity property '{group_name}' should not be a concern group"
            )
