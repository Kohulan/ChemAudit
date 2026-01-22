"""
Scoring API Routes

Endpoints for molecule scoring including ML-readiness, NP-likeness, and scaffold extraction.
"""
from fastapi import APIRouter, HTTPException
from rdkit import Chem
from rdkit.Chem import Descriptors, rdMolDescriptors
import time

from app.schemas.scoring import (
    ScoringRequest,
    ScoringResponse,
    MoleculeInfoSchema,
    MLReadinessResultSchema,
    MLReadinessBreakdownSchema,
    NPLikenessResultSchema,
    ScaffoldResultSchema,
)
from app.services.parser.molecule_parser import parse_molecule, MoleculeFormat
from app.services.scoring import (
    calculate_ml_readiness,
    calculate_np_likeness,
    extract_scaffold,
)


router = APIRouter()


@router.post("/score", response_model=ScoringResponse)
async def score_molecule(request: ScoringRequest):
    """
    Calculate scores for a molecule.

    Supported scoring types:
    - ml_readiness: ML-readiness score (0-100) with breakdown
    - np_likeness: Natural product likeness score (-5 to +5)
    - scaffold: Murcko scaffold extraction

    Args:
        request: ScoringRequest with molecule and options

    Returns:
        ScoringResponse with requested scoring results

    Raises:
        HTTPException: If molecule cannot be parsed
    """
    start_time = time.time()

    # Parse molecule
    format_map = {
        "smiles": MoleculeFormat.SMILES,
        "inchi": MoleculeFormat.INCHI,
        "mol": MoleculeFormat.MOL,
        "auto": None
    }
    input_format = format_map.get(request.format)

    parse_result = parse_molecule(request.molecule, input_format)

    if not parse_result.success or parse_result.mol is None:
        raise HTTPException(
            status_code=400,
            detail={
                "error": "Failed to parse molecule",
                "errors": parse_result.errors,
                "warnings": parse_result.warnings,
                "format_detected": parse_result.format_detected.value
            }
        )

    mol = parse_result.mol

    # Extract basic molecule info
    mol_info = _extract_molecule_info(mol, request.molecule)

    # Initialize response components
    ml_readiness_result = None
    np_likeness_result = None
    scaffold_result = None

    # Calculate requested scores
    if "ml_readiness" in request.include:
        ml_result = calculate_ml_readiness(mol)
        ml_readiness_result = MLReadinessResultSchema(
            score=ml_result.score,
            breakdown=MLReadinessBreakdownSchema(
                descriptors_score=ml_result.breakdown.descriptors_score,
                descriptors_max=ml_result.breakdown.descriptors_max,
                descriptors_successful=ml_result.breakdown.descriptors_successful,
                descriptors_total=ml_result.breakdown.descriptors_total,
                fingerprints_score=ml_result.breakdown.fingerprints_score,
                fingerprints_max=ml_result.breakdown.fingerprints_max,
                fingerprints_successful=ml_result.breakdown.fingerprints_successful,
                fingerprints_failed=ml_result.breakdown.fingerprints_failed,
                size_score=ml_result.breakdown.size_score,
                size_max=ml_result.breakdown.size_max,
                molecular_weight=ml_result.breakdown.molecular_weight,
                num_atoms=ml_result.breakdown.num_atoms,
                size_category=ml_result.breakdown.size_category,
            ),
            interpretation=ml_result.interpretation,
            failed_descriptors=ml_result.failed_descriptors,
        )

    if "np_likeness" in request.include:
        np_result = calculate_np_likeness(mol)
        np_likeness_result = NPLikenessResultSchema(
            score=np_result.score,
            interpretation=np_result.interpretation,
            caveats=np_result.caveats,
            details=np_result.details,
        )

    if "scaffold" in request.include:
        scaffold_res = extract_scaffold(mol)
        scaffold_result = ScaffoldResultSchema(
            scaffold_smiles=scaffold_res.scaffold_smiles,
            generic_scaffold_smiles=scaffold_res.generic_scaffold_smiles,
            has_scaffold=scaffold_res.has_scaffold,
            message=scaffold_res.message,
            details=scaffold_res.details,
        )

    execution_time = int((time.time() - start_time) * 1000)

    return ScoringResponse(
        molecule_info=mol_info,
        ml_readiness=ml_readiness_result,
        np_likeness=np_likeness_result,
        scaffold=scaffold_result,
        execution_time_ms=execution_time,
    )


def _extract_molecule_info(mol: Chem.Mol, input_string: str) -> MoleculeInfoSchema:
    """
    Extract basic molecule information.

    Args:
        mol: RDKit molecule object
        input_string: Original input string

    Returns:
        MoleculeInfoSchema with basic properties
    """
    try:
        canonical = Chem.MolToSmiles(mol)
        formula = rdMolDescriptors.CalcMolFormula(mol)
        mw = Descriptors.MolWt(mol)
    except Exception:
        canonical = None
        formula = None
        mw = None

    return MoleculeInfoSchema(
        input_string=input_string,
        canonical_smiles=canonical,
        molecular_formula=formula,
        molecular_weight=mw,
    )
