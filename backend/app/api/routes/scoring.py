"""
Scoring API Routes

Endpoints for molecule scoring including:
- ML-readiness
- NP-likeness
- Scaffold extraction
- Drug-likeness (Lipinski, QED, Veber, Ro3, Ghose, Egan, Muegge)
- Safety filters (PAINS, Brenk, NIH, ZINC, ChEMBL)
- ADMET predictions (SAscore, ESOL, Fsp3, CNS MPO, Pfizer/GSK rules)
- Aggregator likelihood prediction
"""

from __future__ import annotations

import time
from typing import TYPE_CHECKING, Optional

if TYPE_CHECKING:
    pass

from fastapi import APIRouter, Depends, HTTPException, Request

from app.api.routes.validation import extract_molecule_info
from app.core.rate_limit import get_rate_limit_key, limiter
from app.core.security import get_api_key
from app.schemas.scoring import (
    ComparisonRequest,
    MLReadinessBreakdownSchema,
    MLReadinessResultSchema,
    NPLikenessResultSchema,
    RadarAxisSchema,
    RadarComparisonSchema,
    RadarProfileSchema,
    ScaffoldResultSchema,
    ScoringRequest,
    ScoringResponse,
)
from app.services.parser.molecule_parser import MoleculeFormat, parse_molecule
from app.services.scoring import (
    calculate_ml_readiness,
    calculate_np_likeness,
    calculate_radar_comparison,
    extract_scaffold,
)
from app.services.scoring.score_builders import (
    _calculate_admet,
    _calculate_aggregator,
    _calculate_bertz_detail,
    _calculate_bioavailability_radar_route,
    _calculate_boiled_egg_route,
    _calculate_consensus,
    _calculate_druglikeness,
    _calculate_fsp3_detail,
    _calculate_lead_likeness,
    _calculate_ligand_efficiency,
    _calculate_logp_breakdown,
    _calculate_np_breakdown_route,
    _calculate_safety_filters,
    _calculate_salt_inventory,
    _calculate_tpsa_breakdown,
    _ml_dimension_to_schema,
)

router = APIRouter()


@router.post("/score", response_model=ScoringResponse)
@limiter.limit("10/minute", key_func=get_rate_limit_key)
async def score_molecule(
    request: Request,
    body: ScoringRequest,
    api_key: Optional[str] = Depends(get_api_key),
):
    """
    Calculate scores for a molecule.

    Supported scoring types:
    - ml_readiness: ML-readiness score (0-100) with breakdown
    - np_likeness: Natural product likeness score (-5 to +5)
    - scaffold: Murcko scaffold extraction
    - druglikeness: Drug-likeness filters (Lipinski, QED, Veber, Ro3, etc.)
    - safety_filters: Safety alerts (PAINS, Brenk, NIH, ZINC, ChEMBL)
    - admet: ADMET predictions (SAscore, solubility, Pfizer/GSK rules, etc.)
    - aggregator: Aggregator likelihood prediction

    Args:
        body: ScoringRequest with molecule and options

    Returns:
        ScoringResponse with requested scoring results

    Raises:
        HTTPException: If molecule cannot be parsed
    """
    start_time = time.time()

    # IUPAC/common name resolution (try OPSIN then PubChem)
    molecule_input = body.molecule
    if body.format in ("auto", "iupac"):
        from app.services.iupac.converter import detect_input_type, name_to_smiles

        detected = detect_input_type(molecule_input) if body.format == "auto" else "iupac"
        if detected == "iupac" or (detected == "ambiguous" and body.format == "auto"):
            converted, _source = name_to_smiles(molecule_input)
            if converted:
                molecule_input = converted

    # Parse molecule
    format_map = {
        "smiles": MoleculeFormat.SMILES,
        "inchi": MoleculeFormat.INCHI,
        "mol": MoleculeFormat.MOL,
        "iupac": MoleculeFormat.SMILES,
        "auto": None,
    }
    input_format = format_map.get(body.format)

    parse_result = parse_molecule(molecule_input, input_format)

    if not parse_result.success or parse_result.mol is None:
        raise HTTPException(
            status_code=400,
            detail={
                "error": "Failed to parse molecule",
                "errors": parse_result.errors,
                "warnings": parse_result.warnings,
                "format_detected": parse_result.format_detected.value,
            },
        )

    mol = parse_result.mol

    # Extract basic molecule info
    mol_info = extract_molecule_info(mol, body.molecule)

    # Initialize response components
    ml_readiness_result = None
    np_likeness_result = None
    scaffold_result = None
    druglikeness_result = None
    safety_filters_result = None
    admet_result = None
    aggregator_result = None
    consensus_result = None
    lead_likeness_result = None
    salt_inventory_result = None
    ligand_efficiency_result = None
    tpsa_breakdown_result = None
    logp_breakdown_result = None
    bertz_detail_result = None
    fsp3_detail_result = None
    np_breakdown_result = None
    bioavailability_radar_result = None
    boiled_egg_result = None

    # Calculate requested scores
    if "ml_readiness" in body.include:
        ml_result = calculate_ml_readiness(mol)
        ml_readiness_result = MLReadinessResultSchema(
            score=ml_result.score,
            label=ml_result.label,
            color=ml_result.color,
            breakdown=MLReadinessBreakdownSchema(
                structural_quality=_ml_dimension_to_schema(ml_result.breakdown.structural_quality),
                property_profile=_ml_dimension_to_schema(ml_result.breakdown.property_profile),
                complexity_feasibility=_ml_dimension_to_schema(
                    ml_result.breakdown.complexity_feasibility,
                ),
                representation_quality=_ml_dimension_to_schema(
                    ml_result.breakdown.representation_quality,
                ),
            ),
            caveats=ml_result.caveats,
            supplementary=ml_result.supplementary,
            interpretation=ml_result.interpretation,
        )

    if "np_likeness" in body.include:
        np_result = calculate_np_likeness(mol)
        np_likeness_result = NPLikenessResultSchema(
            score=np_result.score,
            interpretation=np_result.interpretation,
            caveats=np_result.caveats,
            details=np_result.details,
        )

    if "scaffold" in body.include:
        scaffold_res = extract_scaffold(mol)
        scaffold_result = ScaffoldResultSchema(
            scaffold_smiles=scaffold_res.scaffold_smiles,
            generic_scaffold_smiles=scaffold_res.generic_scaffold_smiles,
            has_scaffold=scaffold_res.has_scaffold,
            message=scaffold_res.message,
            details=scaffold_res.details,
        )

    if "druglikeness" in body.include:
        druglikeness_result = _calculate_druglikeness(mol)

    if "safety_filters" in body.include:
        safety_filters_result = _calculate_safety_filters(mol)

    if "admet" in body.include:
        admet_result = _calculate_admet(mol)

    if "aggregator" in body.include:
        aggregator_result = _calculate_aggregator(mol)

    if "consensus" in body.include:
        consensus_result = _calculate_consensus(mol)

    if "lead_likeness" in body.include:
        lead_likeness_result = _calculate_lead_likeness(mol)

    if "salt_inventory" in body.include:
        salt_inventory_result = _calculate_salt_inventory(mol)

    if "ligand_efficiency" in body.include:
        ligand_efficiency_result = _calculate_ligand_efficiency(mol)

    if "tpsa_breakdown" in body.include:
        tpsa_breakdown_result = _calculate_tpsa_breakdown(mol)

    if "logp_breakdown" in body.include:
        logp_breakdown_result = _calculate_logp_breakdown(mol)

    if "bertz_detail" in body.include:
        bertz_detail_result = _calculate_bertz_detail(mol)

    if "fsp3_detail" in body.include:
        fsp3_detail_result = _calculate_fsp3_detail(mol)

    if "np_breakdown" in body.include:
        np_breakdown_result = _calculate_np_breakdown_route(mol)

    if "bioavailability_radar" in body.include:
        bioavailability_radar_result = _calculate_bioavailability_radar_route(mol)

    if "boiled_egg" in body.include:
        boiled_egg_result = _calculate_boiled_egg_route(mol)

    execution_time = int((time.time() - start_time) * 1000)

    return ScoringResponse(
        molecule_info=mol_info,
        ml_readiness=ml_readiness_result,
        np_likeness=np_likeness_result,
        scaffold=scaffold_result,
        druglikeness=druglikeness_result,
        safety_filters=safety_filters_result,
        admet=admet_result,
        aggregator=aggregator_result,
        consensus=consensus_result,
        lead_likeness=lead_likeness_result,
        salt_inventory=salt_inventory_result,
        ligand_efficiency=ligand_efficiency_result,
        tpsa_breakdown=tpsa_breakdown_result,
        logp_breakdown=logp_breakdown_result,
        bertz_detail=bertz_detail_result,
        fsp3_detail=fsp3_detail_result,
        np_breakdown=np_breakdown_result,
        bioavailability_radar=bioavailability_radar_result,
        boiled_egg=boiled_egg_result,
        execution_time_ms=execution_time,
    )


@router.post("/score/compare", response_model=RadarComparisonSchema)
@limiter.limit("10/minute", key_func=get_rate_limit_key)
async def compare_molecules(
    request: Request,
    body: ComparisonRequest,
    api_key: Optional[str] = Depends(get_api_key),
):
    """
    Compare multiple molecules using bioavailability radar profiles.

    Accepts up to 10 SMILES strings and returns radar profiles for each,
    plus a drug-like reference profile.

    Args:
        body: ComparisonRequest with list of SMILES

    Returns:
        RadarComparisonSchema with profiles and reference
    """
    result = calculate_radar_comparison(body.smiles_list)

    profiles = [
        RadarProfileSchema(
            smiles=p.smiles,
            axes=[
                RadarAxisSchema(
                    name=a.name,
                    actual_value=a.actual_value,
                    normalized=a.normalized,
                    optimal_min=a.optimal_min,
                    optimal_max=a.optimal_max,
                    in_range=a.in_range,
                    property_name=a.property_name,
                    unit=a.unit,
                )
                for a in p.axes
            ],
            is_reference=p.is_reference,
        )
        for p in result.profiles
    ]

    reference = None
    if result.reference:
        reference = RadarProfileSchema(
            smiles=result.reference.smiles,
            axes=[
                RadarAxisSchema(
                    name=a.name,
                    actual_value=a.actual_value,
                    normalized=a.normalized,
                    optimal_min=a.optimal_min,
                    optimal_max=a.optimal_max,
                    in_range=a.in_range,
                    property_name=a.property_name,
                    unit=a.unit,
                )
                for a in result.reference.axes
            ],
            is_reference=result.reference.is_reference,
        )

    return RadarComparisonSchema(
        profiles=profiles,
        reference=reference,
    )
