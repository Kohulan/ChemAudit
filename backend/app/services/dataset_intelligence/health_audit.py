"""
Dataset Health Audit Service Module (DSET-01).

Computes a composite 0-100 health score for a dataset of molecules based on
5 sub-scores: parsability, stereo completeness, uniqueness, alert prevalence,
and standardization consistency. Includes property distribution histograms
with drug-space reference overlay data.

Reference weights: Fourches et al. 2010, Mansouri et al. 2024
"""

from __future__ import annotations

import logging
import random
from collections import defaultdict
from dataclasses import dataclass, field
from typing import Any, Callable

import numpy as np
from rdkit.Chem import Descriptors, FindMolChiralCenters

from app.services.alerts.unified_screen import unified_screen
from app.services.diagnostics.std_comparison import compare_pipelines

logger = logging.getLogger(__name__)

# ---------------------------------------------------------------------------
# Drug-space reference histogram (typical ChEMBL drug-like distributions)
# ---------------------------------------------------------------------------
DRUG_SPACE_REFERENCE: dict[str, dict] = {
    "mw": {
        "bin_edges": [
            100, 135, 170, 205, 240, 275, 310, 345, 380, 415,
            450, 485, 520, 555, 590, 625, 660, 695, 730, 765, 800,
        ],
        "reference_counts": [
            120, 450, 890, 1350, 1780, 2100, 2350, 2200, 1900, 1550,
            1200, 880, 620, 420, 280, 180, 110, 65, 35, 15,
        ],
    },
    "logp": {
        "bin_edges": [
            -2.0, -1.5, -1.0, -0.5, 0.0, 0.5, 1.0, 1.5, 2.0, 2.5,
            3.0, 3.5, 4.0, 4.5, 5.0, 5.5, 6.0, 6.5, 7.0, 7.5, 8.0,
        ],
        "reference_counts": [
            80, 150, 320, 580, 950, 1400, 1850, 2200, 2450, 2500,
            2300, 1950, 1500, 1100, 750, 480, 280, 150, 70, 30,
        ],
    },
    "tpsa": {
        "bin_edges": [
            0, 10, 20, 30, 40, 50, 60, 70, 80, 90,
            100, 110, 120, 130, 140, 150, 160, 170, 180, 190, 200,
        ],
        "reference_counts": [
            200, 450, 850, 1300, 1700, 2000, 2200, 2100, 1850, 1550,
            1250, 950, 700, 500, 350, 230, 140, 80, 40, 15,
        ],
    },
}

# Default weights (literature-backed: Fourches et al. 2010)
DEFAULT_WEIGHTS: dict[str, float] = {
    "parsability": 0.25,
    "stereo": 0.15,
    "uniqueness": 0.20,
    "alerts": 0.20,
    "std_consistency": 0.20,
}


@dataclass
class DatasetHealthResult:
    """Result of a dataset health audit."""

    overall_score: float = 0.0
    parsability_score: float = 0.0
    stereo_score: float = 0.0
    uniqueness_score: float = 0.0
    alert_score: float = 0.0
    std_consistency_score: float = 0.0
    weights: dict[str, float] = field(default_factory=lambda: dict(DEFAULT_WEIGHTS))

    molecule_count: int = 0
    parse_failures: int = 0
    stereo_undefined_count: int = 0
    duplicate_count: int = 0
    alert_hit_count: int = 0
    std_disagreement_count: int = 0
    std_sample_size: int = 0

    issues: list[dict] = field(default_factory=list)
    property_distributions: dict = field(default_factory=dict)
    std_pipeline_comparison: dict = field(default_factory=dict)
    contradictions: list[dict] = field(default_factory=list)
    dedup_groups: list[dict] = field(default_factory=list)
    numeric_columns: list[str] = field(default_factory=list)


def compute_health_score(
    molecules: list[dict],
    weights: dict[str, float] | None = None,
    progress_callback: Callable[[str, float], None] | None = None,
) -> DatasetHealthResult:
    """Compute a composite health score for a dataset of molecules.

    Each molecule in ``molecules`` is a dict:
        ``{"index": int, "smiles": str, "mol": Chem.Mol | None,
           "inchikey": str | None, "properties": dict}``

    Args:
        molecules: List of molecule dicts with pre-parsed mol objects.
        weights: Optional custom weights for sub-scores. Uses defaults if None.
        progress_callback: Optional callback(stage_name, fraction) for progress.

    Returns:
        DatasetHealthResult with 5 sub-scores, overall composite, distributions,
        issues, dedup groups, and standardization comparison data.
    """
    if weights is None:
        weights = dict(DEFAULT_WEIGHTS)

    total_count = len(molecules)
    result = DatasetHealthResult(weights=weights, molecule_count=total_count)

    if total_count == 0:
        return result

    def _progress(stage: str, frac: float) -> None:
        if progress_callback is not None:
            progress_callback(stage, frac)

    # ----- Stage 1: Parsability -----
    _progress("parsability", 0.0)
    parsed_mols: list[dict] = []
    for mol_dict in molecules:
        if mol_dict["mol"] is not None:
            parsed_mols.append(mol_dict)
        else:
            result.parse_failures += 1
            result.issues.append({
                "row_index": mol_dict["index"],
                "smiles": mol_dict["smiles"],
                "issue_type": "parse_failure",
                "severity": "critical",
                "description": f"SMILES could not be parsed: {mol_dict['smiles']!r}",
            })

    parsed_count = len(parsed_mols)
    result.parsability_score = parsed_count / total_count if total_count > 0 else 0.0
    _progress("parsability", 1.0)

    # ----- Stage 2: Stereo completeness -----
    _progress("stereo", 0.0)
    fully_defined_count = 0
    for mol_dict in parsed_mols:
        mol = mol_dict["mol"]
        try:
            centers = FindMolChiralCenters(mol, includeUnassigned=True)
            has_undefined = any(c[1] == "?" for c in centers)
            if has_undefined:
                result.stereo_undefined_count += 1
                result.issues.append({
                    "row_index": mol_dict["index"],
                    "smiles": mol_dict["smiles"],
                    "issue_type": "undefined_stereo",
                    "severity": "warning",
                    "description": "Molecule has undefined stereocenters",
                })
            else:
                fully_defined_count += 1
        except Exception:
            fully_defined_count += 1  # No stereo info available — assume OK

    result.stereo_score = (
        fully_defined_count / parsed_count if parsed_count > 0 else 1.0
    )
    _progress("stereo", 1.0)

    # ----- Stage 3: Uniqueness (InChIKey deduplication) -----
    _progress("uniqueness", 0.0)
    inchikey_groups: dict[str, list[int]] = defaultdict(list)
    molecules_with_ik = 0
    for mol_dict in parsed_mols:
        ik = mol_dict.get("inchikey")
        if ik is not None:
            inchikey_groups[ik].append(mol_dict["index"])
            molecules_with_ik += 1

    unique_count = len(inchikey_groups)
    result.duplicate_count = molecules_with_ik - unique_count

    # Track dedup groups with >1 entry
    for ik, rows in inchikey_groups.items():
        if len(rows) > 1:
            result.dedup_groups.append({"inchikey": ik, "rows": rows})
            for row_idx in rows[1:]:  # First occurrence is not flagged
                result.issues.append({
                    "row_index": row_idx,
                    "smiles": "",
                    "issue_type": "duplicate",
                    "severity": "info",
                    "description": f"Duplicate InChIKey: {ik}",
                })

    result.uniqueness_score = (
        unique_count / molecules_with_ik if molecules_with_ik > 0 else 1.0
    )
    _progress("uniqueness", 1.0)

    # ----- Stage 4: Alert prevalence -----
    _progress("alerts", 0.0)
    alert_hit_count = 0
    for i, mol_dict in enumerate(parsed_mols):
        mol = mol_dict["mol"]
        try:
            screen_result = unified_screen(mol)
            if screen_result.get("total_deduped", 0) > 0:
                alert_hit_count += 1
                result.issues.append({
                    "row_index": mol_dict["index"],
                    "smiles": mol_dict["smiles"],
                    "issue_type": "alert_hit",
                    "severity": "warning",
                    "description": (
                        f"Triggered {screen_result['total_deduped']} structural alert(s)"
                    ),
                })
        except Exception as exc:
            logger.warning("Alert screening failed for row %d: %s", mol_dict["index"], exc)
        if (i + 1) % 50 == 0:
            _progress("alerts", (i + 1) / len(parsed_mols))

    result.alert_hit_count = alert_hit_count
    result.alert_score = (
        1.0 - (alert_hit_count / parsed_count) if parsed_count > 0 else 1.0
    )
    _progress("alerts", 1.0)

    # ----- Stage 5: Standardization consistency (sampled) -----
    _progress("std_consistency", 0.0)
    sample_size = min(parsed_count, 500)
    sample = random.sample(parsed_mols, sample_size) if sample_size > 0 else []
    result.std_sample_size = sample_size

    agree_count = 0
    disagree_count = 0
    _PIPELINE_KEYS = ["rdkit_molstandardize", "chembl_style", "minimal_sanitize"]
    per_pipeline_agree = {k: 0 for k in _PIPELINE_KEYS}
    per_pipeline_disagree = {k: 0 for k in _PIPELINE_KEYS}

    for i, mol_dict in enumerate(sample):
        try:
            comp = compare_pipelines(mol_dict["smiles"])
            if comp.get("all_agree", False):
                agree_count += 1
                for k in _PIPELINE_KEYS:
                    per_pipeline_agree[k] += 1
            else:
                disagree_count += 1
                result.issues.append({
                    "row_index": mol_dict["index"],
                    "smiles": mol_dict["smiles"],
                    "issue_type": "std_disagreement",
                    "severity": "info",
                    "description": "Standardization pipelines disagree",
                })
                pipelines = comp.get("pipelines", [])
                smiles_vals = [p.get("smiles") for p in pipelines]
                for idx, k in enumerate(_PIPELINE_KEYS):
                    smi = smiles_vals[idx] if idx < len(smiles_vals) else None
                    others = [
                        s for j, s in enumerate(smiles_vals)
                        if j != idx and s is not None
                    ]
                    if smi is not None and others and any(s == smi for s in others):
                        per_pipeline_agree[k] += 1
                    else:
                        per_pipeline_disagree[k] += 1
        except Exception as exc:
            disagree_count += 1
            for k in _PIPELINE_KEYS:
                per_pipeline_disagree[k] += 1
            logger.warning(
                "compare_pipelines failed for row %d: %s", mol_dict["index"], exc
            )
        if (i + 1) % 50 == 0:
            _progress("std_consistency", (i + 1) / sample_size)

    result.std_disagreement_count = disagree_count
    result.std_consistency_score = (
        agree_count / sample_size if sample_size > 0 else 1.0
    )
    result.std_pipeline_comparison = {
        "sample_size": sample_size,
        "agree_count": agree_count,
        "disagree_count": disagree_count,
        **{k: {"agree": per_pipeline_agree[k], "disagree": per_pipeline_disagree[k]}
           for k in _PIPELINE_KEYS},
    }
    _progress("std_consistency", 1.0)

    # ----- Stage 6: Property distributions -----
    _progress("property_distributions", 0.0)
    mw_values: list[float] = []
    logp_values: list[float] = []
    tpsa_values: list[float] = []

    for mol_dict in parsed_mols:
        mol = mol_dict["mol"]
        try:
            mw_values.append(Descriptors.MolWt(mol))
            logp_values.append(Descriptors.MolLogP(mol))
            tpsa_values.append(Descriptors.TPSA(mol))
        except Exception:
            pass

    property_distributions: dict[str, Any] = {}
    dist_configs = [
        ("mw", mw_values, 100.0, 800.0),
        ("logp", logp_values, -2.0, 8.0),
        ("tpsa", tpsa_values, 0.0, 200.0),
    ]
    for name, values, lo, hi in dist_configs:
        if values:
            counts, bin_edges = np.histogram(values, bins=20, range=(lo, hi))
            property_distributions[name] = {
                "bins": bin_edges.tolist(),
                "counts": counts.tolist(),
            }
        else:
            property_distributions[name] = {"bins": [], "counts": []}

    result.property_distributions = property_distributions
    _progress("property_distributions", 1.0)

    # ----- Overall composite score -----
    sub_scores = {
        "parsability": result.parsability_score,
        "stereo": result.stereo_score,
        "uniqueness": result.uniqueness_score,
        "alerts": result.alert_score,
        "std_consistency": result.std_consistency_score,
    }
    result.overall_score = (
        sum(sub_scores[k] * weights[k] for k in sub_scores) * 100
    )

    return result
