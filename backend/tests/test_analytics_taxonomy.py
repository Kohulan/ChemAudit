"""
Tests for Chemical Taxonomy Classification Service.

Validates taxonomy classification with known molecules, multi-category matching,
batch classification, and SMARTS rule compilation.
"""

from __future__ import annotations

from rdkit import Chem


class TestTaxonomyRules:
    """Tests for CHEMOTYPE_RULES and get_compiled_rules."""

    def test_all_rules_compile(self):
        """Every SMARTS in CHEMOTYPE_RULES compiles via Chem.MolFromSmarts()."""
        from app.services.analytics.taxonomy_rules import CHEMOTYPE_RULES

        for rule in CHEMOTYPE_RULES:
            mol = Chem.MolFromSmarts(rule["smarts"])
            assert mol is not None, (
                f"SMARTS failed to compile for rule '{rule['name']}': {rule['smarts']}"
            )

    def test_rule_count(self):
        """CHEMOTYPE_RULES has between 45 and 55 entries."""
        from app.services.analytics.taxonomy_rules import CHEMOTYPE_RULES

        assert 45 <= len(CHEMOTYPE_RULES) <= 55, (
            f"Expected 45-55 rules, got {len(CHEMOTYPE_RULES)}"
        )

    def test_get_compiled_rules_returns_tuples(self):
        """get_compiled_rules returns list of (rule_dict, compiled_mol) tuples."""
        from app.services.analytics.taxonomy_rules import get_compiled_rules

        compiled = get_compiled_rules()
        assert len(compiled) > 0

        for rule_dict, compiled_mol in compiled:
            assert isinstance(rule_dict, dict)
            assert "name" in rule_dict
            assert "smarts" in rule_dict
            assert "category" in rule_dict
            assert "description" in rule_dict
            assert compiled_mol is not None

    def test_compiled_rules_singleton(self):
        """get_compiled_rules returns the same object on repeated calls (singleton)."""
        from app.services.analytics.taxonomy_rules import get_compiled_rules

        first = get_compiled_rules()
        second = get_compiled_rules()
        assert first is second


class TestTaxonomyClassification:
    """Tests for classify_molecule and classify_batch."""

    def test_classify_aspirin(self):
        """Aspirin matches at least Esters category."""
        from app.services.analytics.taxonomy import classify_molecule
        from app.services.analytics.taxonomy_rules import get_compiled_rules

        mol = Chem.MolFromSmiles("CC(=O)Oc1ccccc1C(=O)O")  # aspirin
        compiled = get_compiled_rules()
        categories = classify_molecule(mol, compiled)

        category_names = [c["category"] for c in categories]
        # Aspirin has an ester group and a phenol-like structure
        assert any(
            cat in category_names for cat in ["Esters", "Phenols"]
        ), f"Expected Esters or Phenols in {category_names}"

    def test_classify_diazepam(self):
        """Diazepam matches Benzodiazepines category."""
        from app.services.analytics.taxonomy import classify_molecule
        from app.services.analytics.taxonomy_rules import get_compiled_rules

        mol = Chem.MolFromSmiles("O=C1CN=C(c2ccccc2)c2cc(Cl)ccc2N1C")  # diazepam
        compiled = get_compiled_rules()
        categories = classify_molecule(mol, compiled)

        category_names = [c["category"] for c in categories]
        assert "Benzodiazepines" in category_names, (
            f"Expected Benzodiazepines in {category_names}"
        )

    def test_multi_category(self):
        """A molecule matching multiple categories returns all of them."""
        from app.services.analytics.taxonomy import classify_molecule
        from app.services.analytics.taxonomy_rules import get_compiled_rules

        # Sulfamethoxazole has sulfonamide + oxazole + aniline groups
        mol = Chem.MolFromSmiles("CC1=CC(NS(=O)(=O)c2ccc(N)cc2)=NO1")
        compiled = get_compiled_rules()
        categories = classify_molecule(mol, compiled)

        # Should match at least 2 categories
        assert len(categories) >= 2, (
            f"Expected multiple categories, got {len(categories)}: "
            f"{[c['category'] for c in categories]}"
        )

    def test_classify_batch(self):
        """Batch classification returns per-molecule assignments and category_counts."""
        from app.services.analytics.taxonomy import classify_batch

        results = [
            {"smiles": "CC(=O)Oc1ccccc1C(=O)O", "index": 0, "status": "success"},  # aspirin
            {"smiles": "c1ccc2[nH]ccc2c1", "index": 1, "status": "success"},        # indole
            {"smiles": "c1ccncc1", "index": 2, "status": "success"},                  # pyridine
            {"smiles": "INVALID", "index": 3, "status": "success"},                    # invalid
        ]
        result = classify_batch(results)

        assert "per_molecule" in result
        assert "category_counts" in result
        assert "total_molecules" in result
        assert "classified_molecules" in result
        assert "unclassified_molecules" in result

        assert result["total_molecules"] == 4
        # Category counts may sum to more than molecule count (D-11)
        total_categorizations = sum(result["category_counts"].values())
        assert total_categorizations >= result["classified_molecules"]

        # Check per_molecule structure
        for entry in result["per_molecule"]:
            assert "index" in entry
            assert "smiles" in entry
            assert "categories" in entry
            assert isinstance(entry["categories"], list)
