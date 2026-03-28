"""
Chemical Taxonomy Classification Service

Classifies molecules against ~50 curated SMARTS chemotype rules with
multi-category matching (D-11: molecules can match multiple categories).

Uses the compiled rules singleton from taxonomy_rules.py.
"""

from __future__ import annotations

import logging
from typing import Any

from rdkit import Chem

from app.services.analytics.taxonomy_rules import get_compiled_rules

logger = logging.getLogger(__name__)


def classify_molecule(
    mol: Chem.Mol,
    compiled_rules: list[tuple[dict, Chem.Mol]],
) -> list[dict]:
    """
    Classify a single molecule against compiled SMARTS chemotype rules.

    A molecule can match multiple categories (D-11). All matching rules are
    returned, not just the first match.

    Args:
        mol: RDKit molecule object to classify.
        compiled_rules: List of (rule_dict, compiled_smarts_mol) tuples
            from ``get_compiled_rules()``.

    Returns:
        List of dicts for matching rules, each with ``name``, ``category``,
        ``description``, and ``smarts`` keys.
    """
    matches: list[dict] = []
    for rule_dict, pattern in compiled_rules:
        try:
            if mol.HasSubstructMatch(pattern):
                matches.append(
                    {
                        "name": rule_dict["name"],
                        "category": rule_dict["category"],
                        "description": rule_dict["description"],
                        "smarts": rule_dict["smarts"],
                    }
                )
        except Exception:
            logger.debug(
                "Substructure match failed for rule '%s'",
                rule_dict["name"],
            )
    return matches


def classify_batch(results: list[dict[str, Any]]) -> dict:
    """
    Classify a batch of molecules against chemotype taxonomy rules.

    For each result with a valid SMILES, classifies against all compiled
    SMARTS rules. Molecules can match multiple categories, so
    ``category_counts`` values may sum to more than ``classified_molecules``
    (D-11).

    Args:
        results: Batch result dicts, each containing at least ``smiles`` and
            ``index`` keys.

    Returns:
        Dict with keys:
        - ``per_molecule``: list of dicts with ``index``, ``smiles``,
          ``categories``
        - ``category_counts``: dict mapping category name to count
        - ``total_molecules``: total input count
        - ``classified_molecules``: count matching at least 1 rule
        - ``unclassified_molecules``: count matching zero rules
    """
    compiled = get_compiled_rules()

    per_molecule: list[dict] = []
    category_counts: dict[str, int] = {}
    classified = 0
    unclassified = 0

    for result in results:
        smiles = result.get("smiles", "")
        index = result.get("index", 0)

        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            per_molecule.append(
                {
                    "index": index,
                    "smiles": smiles,
                    "categories": [],
                }
            )
            unclassified += 1
            continue

        matches = classify_molecule(mol, compiled)

        if matches:
            classified += 1
            for match in matches:
                cat = match["category"]
                category_counts[cat] = category_counts.get(cat, 0) + 1
        else:
            unclassified += 1

        per_molecule.append(
            {
                "index": index,
                "smiles": smiles,
                "categories": [m["category"] for m in matches],
            }
        )

    return {
        "per_molecule": per_molecule,
        "category_counts": category_counts,
        "total_molecules": len(results),
        "classified_molecules": classified,
        "unclassified_molecules": unclassified,
    }
