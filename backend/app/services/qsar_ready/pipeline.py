"""
QSAR-Ready Pipeline Service

10-step configurable curation pipeline for QSAR-ready molecule preparation.
Each step records before/after SMILES, status, and optional detail for full provenance.

Steps (in canonical order):
  1. parse       — SMILES → RDKit mol
  2. metals      — Metal disconnection (MetalDisconnector)
  3. desalt      — Largest fragment selection (LargestFragmentChooser)
  4. normalize   — Functional group normalization (Normalizer)
  5. neutralize  — Charge neutralization (Uncharger)
  6. tautomer    — Canonical tautomer (TautomerEnumerator)
  7. stereo      — Stereochemistry removal (RemoveStereochemistry)
  8. isotope     — Isotope label removal (atom.SetIsotope(0))
  9. filter      — Composition filter (HA count, MW, inorganic rejection)
 10. canonical   — Canonical SMILES generation
"""

import logging
from dataclasses import dataclass, field
from typing import List, Optional

from rdkit import Chem
from rdkit.Chem import inchi as rdkit_inchi
from rdkit.Chem.MolStandardize import rdMolStandardize

logger = logging.getLogger(__name__)


# =============================================================================
# Data Classes
# =============================================================================


@dataclass
class QSARStepResult:
    """Result for a single pipeline step."""

    step_name: str  # "parse", "metals", "desalt", etc.
    step_index: int  # 1-10 (canonical position)
    enabled: bool  # was this step enabled in config?
    status: str = "error"  # "applied" | "no_change" | "skipped" | "error"
    before_smiles: Optional[str] = None
    after_smiles: Optional[str] = None
    detail: Optional[str] = None  # rejection reason, error message, or change description


@dataclass
class QSARReadyConfig:
    """Configuration for the QSAR-Ready pipeline. Each boolean enables/disables a step."""

    enable_metals: bool = True
    enable_desalt: bool = True
    enable_normalize: bool = True
    enable_neutralize: bool = True
    enable_tautomer: bool = True
    enable_stereo_strip: bool = False  # QSAR-2D: True; QSAR-3D: False
    enable_isotope_strip: bool = True
    # Composition filter thresholds (Step 9)
    min_heavy_atoms: int = 3
    max_heavy_atoms: int = 100
    max_mw: float = 1500.0
    remove_inorganics: bool = True

    @classmethod
    def qsar_2d(cls) -> "QSARReadyConfig":
        """QSAR-2D preset: all steps on, strip stereo + isotopes."""
        return cls(enable_stereo_strip=True, enable_isotope_strip=True)

    @classmethod
    def qsar_3d(cls) -> "QSARReadyConfig":
        """QSAR-3D preset: all steps on, preserve stereo, strip isotopes."""
        return cls(enable_stereo_strip=False, enable_isotope_strip=True)

    @classmethod
    def minimal(cls) -> "QSARReadyConfig":
        """Minimal preset: parse + desalt + normalize + canonical only."""
        return cls(
            enable_metals=False,
            enable_neutralize=False,
            enable_tautomer=False,
            enable_stereo_strip=False,
            enable_isotope_strip=False,
            remove_inorganics=False,
        )


@dataclass
class QSARReadyResult:
    """Result of running a SMILES through the QSAR-Ready pipeline."""

    original_smiles: str
    original_inchikey: Optional[str] = None
    curated_smiles: Optional[str] = None
    standardized_inchikey: Optional[str] = None
    inchikey_changed: bool = False  # LOCKED: STATE.md D-14
    status: str = "error"  # "ok" | "rejected" | "duplicate" | "error"
    rejection_reason: Optional[str] = None
    steps: List[QSARStepResult] = field(default_factory=list)


# =============================================================================
# Helper: InChIKey computation
# =============================================================================


def _compute_inchikey(mol: Chem.Mol) -> Optional[str]:
    """
    Compute InChIKey from an RDKit mol object.

    Args:
        mol: RDKit molecule (must not be None).

    Returns:
        InChIKey string or None if computation fails.
    """
    try:
        inchi = rdkit_inchi.MolToInchi(mol)
        if inchi:
            return rdkit_inchi.MolToInchiKey(mol)
    except Exception:
        pass
    return None


# =============================================================================
# Pipeline Steps
# =============================================================================


def _step_parse(smiles: str) -> tuple[Optional[Chem.Mol], QSARStepResult]:
    """Step 1: Parse SMILES string into RDKit mol."""
    step = QSARStepResult(
        step_name="parse",
        step_index=1,
        enabled=True,
        status="error",
        before_smiles=smiles,
    )
    if not smiles or not smiles.strip():
        step.detail = "Empty SMILES string"
        return None, step

    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        step.detail = "Failed to parse SMILES"
        return None, step

    step.status = "applied"
    step.after_smiles = Chem.MolToSmiles(mol)
    return mol, step


def _step_metals(mol: Chem.Mol, enabled: bool) -> tuple[Chem.Mol, QSARStepResult]:
    """Step 2: Disconnect metals using MetalDisconnector."""
    step = QSARStepResult(
        step_name="metals",
        step_index=2,
        enabled=enabled,
    )
    if not enabled:
        step.status = "skipped"
        return mol, step

    try:
        before_smiles = Chem.MolToSmiles(mol)
        disconnector = rdMolStandardize.MetalDisconnector()
        new_mol = disconnector.Disconnect(mol)
        if new_mol is None:
            step.status = "no_change"
            step.before_smiles = before_smiles
            step.after_smiles = before_smiles
            return mol, step
        after_smiles = Chem.MolToSmiles(new_mol)
        step.before_smiles = before_smiles
        step.after_smiles = after_smiles
        step.status = "applied" if before_smiles != after_smiles else "no_change"
        return new_mol, step
    except Exception as exc:
        step.status = "error"
        step.detail = str(exc)
        return mol, step


def _step_desalt(mol: Chem.Mol, enabled: bool) -> tuple[Optional[Chem.Mol], QSARStepResult]:
    """Step 3: Select largest fragment using LargestFragmentChooser."""
    step = QSARStepResult(
        step_name="desalt",
        step_index=3,
        enabled=enabled,
    )
    if not enabled:
        step.status = "skipped"
        return mol, step

    try:
        before_smiles = Chem.MolToSmiles(mol)
        chooser = rdMolStandardize.LargestFragmentChooser()
        new_mol = chooser.choose(mol)
        if new_mol is None:
            # Fragment chooser returned None — degrade gracefully (Pitfall 7)
            step.status = "no_change"
            step.before_smiles = before_smiles
            step.after_smiles = before_smiles
            return mol, step
        after_smiles = Chem.MolToSmiles(new_mol)
        step.before_smiles = before_smiles
        step.after_smiles = after_smiles
        step.status = "applied" if before_smiles != after_smiles else "no_change"
        return new_mol, step
    except Exception as exc:
        step.status = "error"
        step.detail = str(exc)
        return mol, step


def _step_normalize(mol: Chem.Mol, enabled: bool) -> tuple[Chem.Mol, QSARStepResult]:
    """Step 4: Normalize functional groups using Normalizer."""
    step = QSARStepResult(
        step_name="normalize",
        step_index=4,
        enabled=enabled,
    )
    if not enabled:
        step.status = "skipped"
        return mol, step

    try:
        before_smiles = Chem.MolToSmiles(mol)
        normalizer = rdMolStandardize.Normalizer()
        new_mol = normalizer.normalize(mol)
        if new_mol is None:
            step.status = "no_change"
            step.before_smiles = before_smiles
            step.after_smiles = before_smiles
            return mol, step
        after_smiles = Chem.MolToSmiles(new_mol)
        step.before_smiles = before_smiles
        step.after_smiles = after_smiles
        step.status = "applied" if before_smiles != after_smiles else "no_change"
        return new_mol, step
    except Exception as exc:
        step.status = "error"
        step.detail = str(exc)
        return mol, step


def _step_neutralize(mol: Chem.Mol, enabled: bool) -> tuple[Chem.Mol, QSARStepResult]:
    """Step 5: Neutralize charges using Uncharger."""
    step = QSARStepResult(
        step_name="neutralize",
        step_index=5,
        enabled=enabled,
    )
    if not enabled:
        step.status = "skipped"
        return mol, step

    try:
        before_smiles = Chem.MolToSmiles(mol)
        uncharger = rdMolStandardize.Uncharger()
        new_mol = uncharger.uncharge(mol)
        if new_mol is None:
            step.status = "no_change"
            step.before_smiles = before_smiles
            step.after_smiles = before_smiles
            return mol, step
        after_smiles = Chem.MolToSmiles(new_mol)
        step.before_smiles = before_smiles
        step.after_smiles = after_smiles
        step.status = "applied" if before_smiles != after_smiles else "no_change"
        return new_mol, step
    except Exception as exc:
        step.status = "error"
        step.detail = str(exc)
        return mol, step


def _step_tautomer(mol: Chem.Mol, enabled: bool) -> tuple[Chem.Mol, QSARStepResult]:
    """Step 6: Canonicalize tautomer using TautomerEnumerator."""
    step = QSARStepResult(
        step_name="tautomer",
        step_index=6,
        enabled=enabled,
    )
    if not enabled:
        step.status = "skipped"
        return mol, step

    try:
        before_smiles = Chem.MolToSmiles(mol)
        enumerator = rdMolStandardize.TautomerEnumerator()
        new_mol = enumerator.Canonicalize(mol)
        if new_mol is None:
            step.status = "no_change"
            step.before_smiles = before_smiles
            step.after_smiles = before_smiles
            return mol, step
        after_smiles = Chem.MolToSmiles(new_mol)
        step.before_smiles = before_smiles
        step.after_smiles = after_smiles
        step.status = "applied" if before_smiles != after_smiles else "no_change"
        if step.status == "applied" and not enabled:
            # Pitfall 3: TautomerEnumerator may remove E/Z stereo — note in detail
            step.detail = (
                "Tautomer canonicalization may alter E/Z stereochemistry on some bonds"
            )
        return new_mol, step
    except Exception as exc:
        step.status = "error"
        step.detail = str(exc)
        return mol, step


def _step_stereo(mol: Chem.Mol, enabled: bool) -> tuple[Chem.Mol, QSARStepResult]:
    """Step 7: Strip stereochemistry using RemoveStereochemistry (in-place)."""
    step = QSARStepResult(
        step_name="stereo",
        step_index=7,
        enabled=enabled,
    )
    if not enabled:
        step.status = "skipped"
        return mol, step

    try:
        # Capture before_smiles BEFORE in-place mutation (Pitfall 2 analog for stereo)
        before_smiles = Chem.MolToSmiles(mol)
        Chem.RemoveStereochemistry(mol)  # In-place mutation
        after_smiles = Chem.MolToSmiles(mol)
        step.before_smiles = before_smiles
        step.after_smiles = after_smiles
        step.status = "applied" if before_smiles != after_smiles else "no_change"
        return mol, step
    except Exception as exc:
        step.status = "error"
        step.detail = str(exc)
        return mol, step


def _step_isotope(mol: Chem.Mol, enabled: bool) -> tuple[Chem.Mol, QSARStepResult]:
    """Step 8: Strip isotope labels using atom.SetIsotope(0) (in-place)."""
    step = QSARStepResult(
        step_name="isotope",
        step_index=8,
        enabled=enabled,
    )
    if not enabled:
        step.status = "skipped"
        return mol, step

    try:
        # Capture before_smiles BEFORE in-place mutation (Pitfall 2)
        before_smiles = Chem.MolToSmiles(mol)
        for atom in mol.GetAtoms():
            atom.SetIsotope(0)
        after_smiles = Chem.MolToSmiles(mol)
        step.before_smiles = before_smiles
        step.after_smiles = after_smiles
        step.status = "applied" if before_smiles != after_smiles else "no_change"
        return mol, step
    except Exception as exc:
        step.status = "error"
        step.detail = str(exc)
        return mol, step


def _step_filter(
    mol: Chem.Mol,
    config: QSARReadyConfig,
) -> tuple[Optional[Chem.Mol], QSARStepResult, Optional[str]]:
    """
    Step 9: Composition filter — HA count, MW, and inorganic rejection.

    Returns:
        Tuple of (mol_or_None, step_result, rejection_reason).
        mol_or_None is None when the molecule is rejected.
    """
    step = QSARStepResult(
        step_name="filter",
        step_index=9,
        enabled=True,  # Filter always runs (thresholds are always active)
    )

    try:
        from rdkit.Chem import Descriptors

        smiles = Chem.MolToSmiles(mol)
        step.before_smiles = smiles

        heavy_atoms = mol.GetNumHeavyAtoms()
        mw = Descriptors.MolWt(mol)

        # Check min heavy atoms
        if heavy_atoms < config.min_heavy_atoms:
            reason = (
                f"Too few heavy atoms: {heavy_atoms} < {config.min_heavy_atoms}"
            )
            step.status = "applied"
            step.detail = reason
            step.after_smiles = None
            return None, step, reason

        # Check max heavy atoms
        if heavy_atoms > config.max_heavy_atoms:
            reason = (
                f"Too many heavy atoms: {heavy_atoms} > {config.max_heavy_atoms}"
            )
            step.status = "applied"
            step.detail = reason
            step.after_smiles = None
            return None, step, reason

        # Check max MW
        if mw > config.max_mw:
            reason = f"Molecular weight too high: {mw:.1f} > {config.max_mw}"
            step.status = "applied"
            step.detail = reason
            step.after_smiles = None
            return None, step, reason

        # Check inorganics (remove_inorganics enabled — molecule has no carbon)
        if config.remove_inorganics:
            has_carbon = any(
                atom.GetAtomicNum() == 6 for atom in mol.GetAtoms()
            )
            if not has_carbon:
                reason = "Inorganic molecule: no carbon atoms"
                step.status = "applied"
                step.detail = reason
                step.after_smiles = None
                return None, step, reason

        step.status = "no_change"
        step.after_smiles = smiles
        return mol, step, None

    except Exception as exc:
        step.status = "error"
        step.detail = str(exc)
        return mol, step, None


def _step_canonical(mol: Chem.Mol) -> QSARStepResult:
    """Step 10: Generate canonical SMILES (always runs)."""
    step = QSARStepResult(
        step_name="canonical",
        step_index=10,
        enabled=True,
        status="applied",
    )
    try:
        canonical = Chem.MolToSmiles(mol)
        step.before_smiles = canonical
        step.after_smiles = canonical
    except Exception as exc:
        step.status = "error"
        step.detail = str(exc)
    return step


# =============================================================================
# Public API
# =============================================================================


def qsar_ready_single(smiles: str, config: QSARReadyConfig) -> QSARReadyResult:
    """
    Run a single SMILES through the 10-step QSAR-Ready pipeline.

    Computes original_inchikey from the parsed mol BEFORE any pipeline steps.
    Computes standardized_inchikey from the final mol AFTER all steps.
    Sets inchikey_changed = (original_inchikey != standardized_inchikey).

    Args:
        smiles: Input SMILES string.
        config: QSARReadyConfig controlling which steps to enable.

    Returns:
        QSARReadyResult with status "ok", "rejected", or "error".
        Status "duplicate" is set by the caller (batch_processor) after deduplication.
    """
    result = QSARReadyResult(original_smiles=smiles)
    steps: List[QSARStepResult] = []

    # Step 1: Parse
    mol, parse_step = _step_parse(smiles)
    steps.append(parse_step)

    if mol is None:
        result.status = "rejected"
        result.rejection_reason = parse_step.detail or "Failed to parse SMILES"
        result.steps = steps
        # Fill in skipped steps for completeness (indices 2-10)
        _add_skipped_steps(steps, start_index=2, config=config)
        return result

    # Compute original InChIKey AFTER parse, BEFORE any pipeline transformation (D-14)
    result.original_inchikey = _compute_inchikey(mol)

    # Step 2: Metal disconnection
    mol, step = _step_metals(mol, config.enable_metals)
    steps.append(step)
    if mol is None:
        result.status = "error"
        result.rejection_reason = "Null molecule after metals step"
        result.steps = steps
        return result

    # Step 3: Desalt
    mol, step = _step_desalt(mol, config.enable_desalt)
    steps.append(step)
    if mol is None:
        result.status = "error"
        result.rejection_reason = "Null molecule after desalt step"
        result.steps = steps
        return result

    # Step 4: Normalize
    mol, step = _step_normalize(mol, config.enable_normalize)
    steps.append(step)
    if mol is None:
        result.status = "error"
        result.rejection_reason = "Null molecule after normalize step"
        result.steps = steps
        return result

    # Step 5: Neutralize
    mol, step = _step_neutralize(mol, config.enable_neutralize)
    steps.append(step)
    if mol is None:
        result.status = "error"
        result.rejection_reason = "Null molecule after neutralize step"
        result.steps = steps
        return result

    # Step 6: Tautomer
    mol, step = _step_tautomer(mol, config.enable_tautomer)
    steps.append(step)
    if mol is None:
        result.status = "error"
        result.rejection_reason = "Null molecule after tautomer step"
        result.steps = steps
        return result

    # Step 7: Stereo strip
    mol, step = _step_stereo(mol, config.enable_stereo_strip)
    steps.append(step)
    if mol is None:
        result.status = "error"
        result.rejection_reason = "Null molecule after stereo step"
        result.steps = steps
        return result

    # Step 8: Isotope strip
    mol, step = _step_isotope(mol, config.enable_isotope_strip)
    steps.append(step)
    if mol is None:
        result.status = "error"
        result.rejection_reason = "Null molecule after isotope step"
        result.steps = steps
        return result

    # Step 9: Composition filter
    mol, step, rejection_reason = _step_filter(mol, config)
    steps.append(step)
    if mol is None:
        result.status = "rejected"
        result.rejection_reason = rejection_reason
        result.steps = steps
        # Step 10 not reachable — add as skipped
        steps.append(
            QSARStepResult(
                step_name="canonical",
                step_index=10,
                enabled=True,
                status="skipped",
                detail="Skipped: molecule rejected in filter step",
            )
        )
        return result

    # Step 10: Canonical SMILES
    canonical_step = _step_canonical(mol)
    steps.append(canonical_step)

    # Compute standardized InChIKey AFTER all steps (Pitfall 1 avoidance)
    result.standardized_inchikey = _compute_inchikey(mol)
    result.inchikey_changed = result.original_inchikey != result.standardized_inchikey
    result.curated_smiles = canonical_step.after_smiles
    result.status = "ok"
    result.steps = steps

    return result


def _add_skipped_steps(
    steps: List[QSARStepResult],
    start_index: int,
    config: QSARReadyConfig,
) -> None:
    """Append skipped-step placeholders for steps that couldn't run."""
    step_defs = [
        (2, "metals", config.enable_metals),
        (3, "desalt", config.enable_desalt),
        (4, "normalize", config.enable_normalize),
        (5, "neutralize", config.enable_neutralize),
        (6, "tautomer", config.enable_tautomer),
        (7, "stereo", config.enable_stereo_strip),
        (8, "isotope", config.enable_isotope_strip),
        (9, "filter", True),
        (10, "canonical", True),
    ]
    for idx, name, enabled in step_defs:
        if idx >= start_index:
            steps.append(
                QSARStepResult(
                    step_name=name,
                    step_index=idx,
                    enabled=enabled,
                    status="skipped",
                    detail="Skipped: prior step failed",
                )
            )


def qsar_ready_batch(
    smiles_list: list[str],
    config: QSARReadyConfig,
) -> list[QSARReadyResult]:
    """
    Run a list of SMILES strings through the QSAR-Ready pipeline.

    Note: Deduplication by InChIKey is handled here; use this function
    for in-process batch calls. For Celery-based async batch processing,
    use the batch_processor module which also handles progress tracking.

    Args:
        smiles_list: List of SMILES strings.
        config: QSARReadyConfig controlling pipeline behavior.

    Returns:
        List of QSARReadyResult, one per input SMILES.
    """
    results: list[QSARReadyResult] = []
    seen_inchikeys: dict[str, int] = {}

    for i, smiles in enumerate(smiles_list):
        result = qsar_ready_single(smiles, config)
        # Deduplication by standardized InChIKey
        if result.status == "ok" and result.standardized_inchikey:
            key = result.standardized_inchikey
            if key in seen_inchikeys:
                result.status = "duplicate"
                result.rejection_reason = f"Duplicate of molecule {seen_inchikeys[key]}"
            else:
                seen_inchikeys[key] = i
        results.append(result)

    return results
