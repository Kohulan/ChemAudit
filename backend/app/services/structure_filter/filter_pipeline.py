"""
Structure Filter Pipeline.

Implements a 6-stage funnel pipeline for filtering SMILES lists through
structure-based criteria. Each stage rejects molecules that fail specific criteria
and passes remaining molecules to the next stage.

Stages:
  1. parse     — SMILES parsability and non-empty atom count
  2. valence   — Explicit valence check (handled by RDKit during parse)
  3. alerts    — Structural alert screening (PAINS/Brenk/Kazius/NIBR)
  4. property  — Physicochemical property thresholds (MW, LogP, TPSA, RotBonds, Rings)
  5. sa        — SA Score threshold
  6. dedup     — InChIKey-based deduplication
  7. novelty   — ChEMBL Tanimoto similarity (only if config.enable_novelty)
"""

from __future__ import annotations

import logging
import os
import sys
from dataclasses import dataclass, field
from typing import Optional

from rdkit import Chem, RDConfig
from rdkit.Chem import Descriptors, rdMolDescriptors
from rdkit.Chem import inchi as rdkit_inchi

from app.services.alerts.kazius_rules import screen_kazius
from app.services.alerts.nibr_filters import screen_nibr
from app.services.structure_filter.filter_config import FilterConfig
from app.services.scoring.safety_filters import _scorer

# SA Score via RDKit Contrib (Phase 7 pattern)
sys.path.append(os.path.join(RDConfig.RDContribDir, "SA_Score"))
import sascorer  # type: ignore  # noqa: E402

logger = logging.getLogger(__name__)

# Stage name constants
STAGE_PARSE = "parse"
STAGE_VALENCE = "valence"
STAGE_ALERTS = "alerts"
STAGE_PROPERTY = "property"
STAGE_SA = "sa"
STAGE_DEDUP = "dedup"
STAGE_NOVELTY = "novelty"


@dataclass
class StageResult:
    """Result summary for a single pipeline stage."""

    stage_name: str
    stage_index: int
    input_count: int
    passed_count: int
    rejected_count: int
    enabled: bool = True


@dataclass
class MoleculeResult:
    """Per-molecule result tracking status and rejection reason."""

    smiles: str
    status: str  # "passed" | "rejected" | "duplicate" | "error"
    failed_at: Optional[str] = None  # stage name where rejection occurred
    rejection_reason: Optional[str] = None


@dataclass
class FilterResult:
    """Aggregated result from the full filter pipeline."""

    input_count: int
    stages: list[StageResult] = field(default_factory=list)
    molecules: list[MoleculeResult] = field(default_factory=list)
    output_count: int = 0


def filter_batch(smiles_list: list[str], config: FilterConfig) -> FilterResult:
    """
    Run the 6-stage structure filter pipeline on a list of SMILES.

    Stages: parse → valence → alerts → property → sa → dedup (+ novelty if enabled).
    Each stage tracks per-molecule status and aggregate counts.

    Args:
        smiles_list: Input SMILES strings to filter.
        config: FilterConfig instance (use PRESETS dict for standard configurations).

    Returns:
        FilterResult with per-stage StageResult counts and per-molecule MoleculeResult list.
    """
    input_count = len(smiles_list)

    # Track per-molecule state: mol object and result
    mol_objects: list[Optional[Chem.Mol]] = [None] * input_count
    results: list[MoleculeResult] = [
        MoleculeResult(smiles=smi, status="passed") for smi in smiles_list
    ]

    # Active indices = molecules not yet rejected
    active = list(range(input_count))

    stages: list[StageResult] = []

    # -------------------------------------------------------------------
    # Stage 1: Parse
    # -------------------------------------------------------------------
    stage_input = len(active)
    parse_rejected: set[int] = set()
    for idx in active:
        smi = smiles_list[idx]
        mol = Chem.MolFromSmiles(smi)
        if mol is None:
            results[idx].status = "rejected"
            results[idx].failed_at = STAGE_PARSE
            results[idx].rejection_reason = "Invalid SMILES"
            parse_rejected.add(idx)
        elif mol.GetNumAtoms() == 0:
            results[idx].status = "rejected"
            results[idx].failed_at = STAGE_PARSE
            results[idx].rejection_reason = "Empty molecule"
            parse_rejected.add(idx)
        else:
            mol_objects[idx] = mol

    active = [i for i in active if i not in parse_rejected]
    stages.append(
        StageResult(
            stage_name=STAGE_PARSE,
            stage_index=1,
            input_count=stage_input,
            passed_count=len(active),
            rejected_count=len(parse_rejected),
        )
    )

    # -------------------------------------------------------------------
    # Stage 2: Valence
    # (RDKit sanitizes by default in MolFromSmiles; explicit valence errors
    # cause parse failure, so valence survivors == parse survivors)
    # -------------------------------------------------------------------
    stage_input = len(active)
    stages.append(
        StageResult(
            stage_name=STAGE_VALENCE,
            stage_index=2,
            input_count=stage_input,
            passed_count=stage_input,
            rejected_count=0,
        )
    )

    # -------------------------------------------------------------------
    # Stage 3: Alerts
    # -------------------------------------------------------------------
    stage_input = len(active)
    alerts_rejected: set[int] = set()

    alerts_enabled = config.use_pains or config.use_brenk or config.use_kazius or config.use_nibr

    for idx in active:
        mol = mol_objects[idx]
        assert mol is not None  # guaranteed by parse stage

        if alerts_enabled:
            alert_names: list[str] = []

            if config.use_pains:
                alert_names.extend(_scorer.get_pains_alerts(mol))
            if config.use_brenk:
                alert_names.extend(_scorer.get_brenk_alerts(mol))
            if config.use_kazius:
                kazius_hits = screen_kazius(mol)
                alert_names.extend(hit.pattern_name for hit in kazius_hits)
            if config.use_nibr:
                nibr_hits = screen_nibr(mol)
                alert_names.extend(hit.pattern_name for hit in nibr_hits)

            if alert_names:
                first = alert_names[0]
                n_more = len(alert_names) - 1
                reason = f"Alert match: {first}" + (f" (+{n_more} more)" if n_more > 0 else "")
                results[idx].status = "rejected"
                results[idx].failed_at = STAGE_ALERTS
                results[idx].rejection_reason = reason
                alerts_rejected.add(idx)

    active = [i for i in active if i not in alerts_rejected]
    stages.append(
        StageResult(
            stage_name=STAGE_ALERTS,
            stage_index=3,
            input_count=stage_input,
            passed_count=len(active),
            rejected_count=len(alerts_rejected),
            enabled=alerts_enabled,
        )
    )

    # -------------------------------------------------------------------
    # Stage 4: Property
    # -------------------------------------------------------------------
    stage_input = len(active)
    property_rejected: set[int] = set()

    for idx in active:
        mol = mol_objects[idx]
        assert mol is not None

        mw = Descriptors.MolWt(mol)
        logp = Descriptors.MolLogP(mol)
        tpsa = Descriptors.TPSA(mol)
        rot_bonds = rdMolDescriptors.CalcNumRotatableBonds(mol)

        rejection_reason: Optional[str] = None

        if not (config.min_mw <= mw <= config.max_mw):
            rejection_reason = f"MW {mw:.1f} outside range [{config.min_mw}, {config.max_mw}]"
        elif not (config.min_logp <= logp <= config.max_logp):
            rejection_reason = (
                f"LogP {logp:.1f} outside range [{config.min_logp}, {config.max_logp}]"
            )
        elif tpsa > config.max_tpsa:
            rejection_reason = f"TPSA {tpsa:.1f} outside range [0, {config.max_tpsa}]"
        elif rot_bonds > config.max_rot_bonds:
            rejection_reason = (
                f"RotBonds {rot_bonds} outside range [0, {config.max_rot_bonds}]"
            )
        elif config.max_rings is not None:
            rings = rdMolDescriptors.CalcNumRings(mol)
            if rings > config.max_rings:
                rejection_reason = (
                    f"Rings {rings} outside range [0, {config.max_rings}]"
                )

        if rejection_reason is not None:
            results[idx].status = "rejected"
            results[idx].failed_at = STAGE_PROPERTY
            results[idx].rejection_reason = rejection_reason
            property_rejected.add(idx)

    active = [i for i in active if i not in property_rejected]
    stages.append(
        StageResult(
            stage_name=STAGE_PROPERTY,
            stage_index=4,
            input_count=stage_input,
            passed_count=len(active),
            rejected_count=len(property_rejected),
        )
    )

    # -------------------------------------------------------------------
    # Stage 5: SA Score
    # -------------------------------------------------------------------
    stage_input = len(active)
    sa_rejected: set[int] = set()

    for idx in active:
        mol = mol_objects[idx]
        assert mol is not None

        sa_score = sascorer.calculateScore(mol)
        if sa_score > config.max_sa_score:
            results[idx].status = "rejected"
            results[idx].failed_at = STAGE_SA
            results[idx].rejection_reason = (
                f"SA Score {sa_score:.1f} > {config.max_sa_score}"
            )
            sa_rejected.add(idx)

    active = [i for i in active if i not in sa_rejected]
    stages.append(
        StageResult(
            stage_name=STAGE_SA,
            stage_index=5,
            input_count=stage_input,
            passed_count=len(active),
            rejected_count=len(sa_rejected),
        )
    )

    # -------------------------------------------------------------------
    # Stage 6: Dedup (InChIKey-based)
    # -------------------------------------------------------------------
    stage_input = len(active)
    dedup_rejected: set[int] = set()
    seen_inchikeys: dict[str, int] = {}  # inchikey -> first seen index (0-based)

    for idx in active:
        mol = mol_objects[idx]
        assert mol is not None

        inchikey = rdkit_inchi.MolToInchiKey(mol)
        if inchikey in seen_inchikeys:
            first_idx = seen_inchikeys[inchikey]
            results[idx].status = "duplicate"
            results[idx].failed_at = STAGE_DEDUP
            results[idx].rejection_reason = f"Duplicate of row {first_idx + 1}"
            dedup_rejected.add(idx)
        else:
            seen_inchikeys[inchikey] = idx

    active = [i for i in active if i not in dedup_rejected]
    stages.append(
        StageResult(
            stage_name=STAGE_DEDUP,
            stage_index=6,
            input_count=stage_input,
            passed_count=len(active),
            rejected_count=len(dedup_rejected),
        )
    )

    # -------------------------------------------------------------------
    # Stage 7: Novelty (only if enabled)
    # -------------------------------------------------------------------
    if config.enable_novelty:
        from app.services.structure_filter.novelty import check_novelty  # noqa: PLC0415

        stage_input = len(active)
        novelty_rejected: set[int] = set()

        for idx in active:
            mol = mol_objects[idx]
            assert mol is not None

            is_novel, max_sim = check_novelty(mol, config.novelty_threshold)
            if not is_novel:
                sim_str = f"{max_sim:.2f}" if max_sim is not None else "N/A"
                results[idx].status = "rejected"
                results[idx].failed_at = STAGE_NOVELTY
                results[idx].rejection_reason = (
                    f"Max Tanimoto {sim_str} > {config.novelty_threshold}"
                )
                novelty_rejected.add(idx)

        active = [i for i in active if i not in novelty_rejected]
        stages.append(
            StageResult(
                stage_name=STAGE_NOVELTY,
                stage_index=7,
                input_count=stage_input,
                passed_count=len(active),
                rejected_count=len(novelty_rejected),
            )
        )

    return FilterResult(
        input_count=input_count,
        stages=stages,
        molecules=results,
        output_count=len(active),
    )
