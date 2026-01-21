"""
Validation API Routes

Endpoints for single molecule validation.
"""
from fastapi import APIRouter, HTTPException
from rdkit import Chem
from rdkit.Chem import Descriptors, inchi as rdkit_inchi, rdMolDescriptors
import time

from app.schemas.validation import ValidationRequest, ValidationResponse, MoleculeInfo, CheckResultSchema
from app.services.parser.molecule_parser import parse_molecule, MoleculeFormat
from app.services.validation.engine import validation_engine


router = APIRouter()


@router.post("/validate", response_model=ValidationResponse)
async def validate_molecule(request: ValidationRequest):
    """
    Validate a single molecule.

    Args:
        request: Validation request with molecule and options

    Returns:
        ValidationResponse with validation results and score

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

    # Extract molecule info
    mol_info = extract_molecule_info(mol, request.molecule)

    # Run validation
    results, score = validation_engine.validate(mol, request.checks)

    # Convert to schema
    check_results = [
        CheckResultSchema(
            check_name=r.check_name,
            passed=r.passed,
            severity=r.severity,
            message=r.message,
            affected_atoms=r.affected_atoms,
            details=r.details
        )
        for r in results
    ]

    execution_time = int((time.time() - start_time) * 1000)

    return ValidationResponse(
        molecule_info=mol_info,
        overall_score=score,
        issues=[c for c in check_results if not c.passed],
        all_checks=check_results,
        execution_time_ms=execution_time
    )


@router.get("/checks")
async def list_checks():
    """
    List available validation checks.

    Returns:
        Dictionary mapping check categories to check names
    """
    return validation_engine.list_checks()


def extract_molecule_info(mol: Chem.Mol, input_smiles: str) -> MoleculeInfo:
    """
    Extract molecule properties.

    Args:
        mol: RDKit molecule object
        input_smiles: Original input string

    Returns:
        MoleculeInfo with molecular properties
    """
    try:
        canonical = Chem.MolToSmiles(mol)
        mol_inchi = rdkit_inchi.MolToInchi(mol)
        mol_inchikey = rdkit_inchi.MolToInchiKey(mol) if mol_inchi else None
        formula = rdMolDescriptors.CalcMolFormula(mol)
        mw = Descriptors.MolWt(mol)
        num_atoms = mol.GetNumAtoms()
    except Exception:
        canonical = None
        mol_inchi = None
        mol_inchikey = None
        formula = None
        mw = None
        num_atoms = None

    return MoleculeInfo(
        input_smiles=input_smiles,
        canonical_smiles=canonical,
        inchi=mol_inchi,
        inchikey=mol_inchikey,
        molecular_formula=formula,
        molecular_weight=mw,
        num_atoms=num_atoms
    )
