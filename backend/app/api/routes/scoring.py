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
    from app.services.scoring.safety_filters import FilterResult

from fastapi import APIRouter, Depends, HTTPException, Request
from rdkit import Chem

from app.api.routes.validation import extract_molecule_info
from app.core.rate_limit import get_rate_limit_key, limiter
from app.core.security import get_api_key
from app.schemas.scoring import (
    ADMETResultSchema,
    AggregatorLikelihoodSchema,
    AtomContributionSchema,
    BertzDetailSchema,
    BioavailabilityRadarSchema,
    BioavailabilitySchema,
    BoiledEggSchema,
    CarbonHybridizationSchema,
    ChEMBLAlertsSchema,
    CNSMPOSchema,
    ComparisonRequest,
    ComplexitySchema,
    ConsensusScoreSchema,
    DrugLikenessResultSchema,
    EganSchema,
    EllipseParamsSchema,
    FilterAlertSchema,
    Fsp3DetailSchema,
    FunctionalGroupContributionSchema,
    GhoseSchema,
    GoldenTriangleSchema,
    GSKRuleSchema,
    LeadLikenessSchema,
    LigandEfficiencySchema,
    LipinskiSchema,
    LogPBreakdownSchema,
    MLDimensionItemSchema,
    MLDimensionSchema,
    MLReadinessBreakdownSchema,
    MLReadinessResultSchema,
    MueggeSchema,
    NPBreakdownSchema,
    NPFragmentSchema,
    NPLikenessResultSchema,
    PfizerRuleSchema,
    QEDSchema,
    RadarAxisSchema,
    RadarComparisonSchema,
    RadarProfileSchema,
    RuleOfThreeSchema,
    RuleSetDetailSchema,
    RuleViolationSchema,
    SafetyFilterResultSchema,
    SaltFragmentSchema,
    SaltInventorySchema,
    ScaffoldResultSchema,
    ScoringRequest,
    ScoringResponse,
    SolubilitySchema,
    SyntheticAccessibilitySchema,
    TPSABreakdownSchema,
    VeberSchema,
)
from app.services.parser.molecule_parser import MoleculeFormat, parse_molecule
from app.services.scoring import (
    calculate_admet,
    calculate_aggregator_likelihood,
    calculate_bertz_detail,
    calculate_bioavailability_radar,
    calculate_boiled_egg,
    calculate_consensus,
    calculate_druglikeness,
    calculate_fsp3_detail,
    calculate_lead_likeness,
    calculate_ligand_efficiency,
    calculate_logp_breakdown,
    calculate_ml_readiness,
    calculate_np_breakdown,
    calculate_np_likeness,
    calculate_radar_comparison,
    calculate_safety_filters,
    calculate_salt_inventory,
    calculate_tpsa_breakdown,
    extract_scaffold,
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


def _calculate_druglikeness(mol: Chem.Mol) -> DrugLikenessResultSchema:
    """Calculate drug-likeness scores and convert to schema."""
    result = calculate_druglikeness(mol, include_extended=True)

    # Build Lipinski schema
    lipinski = LipinskiSchema(
        passed=result.lipinski.passed,
        violations=result.lipinski.violations,
        mw=result.lipinski.mw,
        logp=result.lipinski.logp,
        hbd=result.lipinski.hbd,
        hba=result.lipinski.hba,
        details=result.lipinski.details,
    )

    # Build QED schema
    qed = QEDSchema(
        score=result.qed.score,
        properties=result.qed.properties,
        interpretation=result.qed.interpretation,
    )

    # Build Veber schema
    veber = VeberSchema(
        passed=result.veber.passed,
        rotatable_bonds=result.veber.rotatable_bonds,
        tpsa=result.veber.tpsa,
    )

    # Build Ro3 schema
    ro3 = RuleOfThreeSchema(
        passed=result.ro3.passed,
        violations=result.ro3.violations,
        mw=result.ro3.mw,
        logp=result.ro3.logp,
        hbd=result.ro3.hbd,
        hba=result.ro3.hba,
        rotatable_bonds=result.ro3.rotatable_bonds,
        tpsa=result.ro3.tpsa,
    )

    # Build optional extended filter schemas
    ghose = None
    if result.ghose:
        ghose = GhoseSchema(
            passed=result.ghose.passed,
            violations=result.ghose.violations,
            mw=result.ghose.mw,
            logp=result.ghose.logp,
            atom_count=result.ghose.atom_count,
            molar_refractivity=result.ghose.molar_refractivity,
        )

    egan = None
    if result.egan:
        egan = EganSchema(
            passed=result.egan.passed,
            logp=result.egan.logp,
            tpsa=result.egan.tpsa,
        )

    muegge = None
    if result.muegge:
        muegge = MueggeSchema(
            passed=result.muegge.passed,
            violations=result.muegge.violations,
            details=result.muegge.details,
        )

    return DrugLikenessResultSchema(
        lipinski=lipinski,
        qed=qed,
        veber=veber,
        ro3=ro3,
        ghose=ghose,
        egan=egan,
        muegge=muegge,
        interpretation=result.interpretation,
    )


def _filter_to_schema(fr: FilterResult) -> FilterAlertSchema:
    """Convert a FilterResult dataclass to its Pydantic schema."""
    return FilterAlertSchema(
        passed=fr.passed,
        alerts=fr.alerts,
        alert_details=fr.alert_details,
        alert_count=fr.alert_count,
    )


def _calculate_safety_filters(mol: Chem.Mol) -> SafetyFilterResultSchema:
    """Calculate safety filter results and convert to schema."""
    result = calculate_safety_filters(mol, include_extended=True, include_chembl=True)

    pains = _filter_to_schema(result.pains)
    brenk = _filter_to_schema(result.brenk)
    nih = _filter_to_schema(result.nih) if result.nih else None
    zinc = _filter_to_schema(result.zinc) if result.zinc else None

    # ChEMBL alerts
    chembl = None
    if result.chembl:
        chembl = ChEMBLAlertsSchema(
            passed=result.chembl.passed,
            total_alerts=result.chembl.total_alerts,
            bms=_filter_to_schema(result.chembl.bms) if result.chembl.bms else None,
            dundee=_filter_to_schema(result.chembl.dundee) if result.chembl.dundee else None,
            glaxo=_filter_to_schema(result.chembl.glaxo) if result.chembl.glaxo else None,
            inpharmatica=(
                _filter_to_schema(result.chembl.inpharmatica)
                if result.chembl.inpharmatica
                else None
            ),
            lint=_filter_to_schema(result.chembl.lint) if result.chembl.lint else None,
            mlsmr=_filter_to_schema(result.chembl.mlsmr) if result.chembl.mlsmr else None,
            schembl=_filter_to_schema(result.chembl.schembl) if result.chembl.schembl else None,
        )

    return SafetyFilterResultSchema(
        pains=pains,
        brenk=brenk,
        nih=nih,
        zinc=zinc,
        chembl=chembl,
        all_passed=result.all_passed,
        total_alerts=result.total_alerts,
        interpretation=result.interpretation,
    )


def _calculate_admet(mol: Chem.Mol) -> ADMETResultSchema:
    """Calculate ADMET predictions and convert to schema."""
    result = calculate_admet(mol, include_cns_mpo=True, include_rules=True)

    synthetic_accessibility = SyntheticAccessibilitySchema(
        score=result.synthetic_accessibility.score,
        classification=result.synthetic_accessibility.classification,
        interpretation=result.synthetic_accessibility.interpretation,
    )

    solubility = SolubilitySchema(
        log_s=result.solubility.log_s,
        solubility_mg_ml=result.solubility.solubility_mg_ml,
        classification=result.solubility.classification,
        interpretation=result.solubility.interpretation,
    )

    complexity = ComplexitySchema(
        fsp3=result.complexity.fsp3,
        num_stereocenters=result.complexity.num_stereocenters,
        num_rings=result.complexity.num_rings,
        num_aromatic_rings=result.complexity.num_aromatic_rings,
        bertz_ct=result.complexity.bertz_ct,
        classification=result.complexity.classification,
        interpretation=result.complexity.interpretation,
    )

    cns_mpo = None
    if result.cns_mpo:
        cns_mpo = CNSMPOSchema(
            score=result.cns_mpo.score,
            components=result.cns_mpo.components,
            cns_penetrant=result.cns_mpo.cns_penetrant,
            interpretation=result.cns_mpo.interpretation,
        )

    bioavailability = BioavailabilitySchema(
        tpsa=result.bioavailability.tpsa,
        rotatable_bonds=result.bioavailability.rotatable_bonds,
        hbd=result.bioavailability.hbd,
        hba=result.bioavailability.hba,
        mw=result.bioavailability.mw,
        logp=result.bioavailability.logp,
        oral_absorption_likely=result.bioavailability.oral_absorption_likely,
        cns_penetration_likely=result.bioavailability.cns_penetration_likely,
        interpretation=result.bioavailability.interpretation,
    )

    # Pfizer 3/75 Rule
    pfizer_rule = None
    if result.pfizer_rule:
        pfizer_rule = PfizerRuleSchema(
            passed=result.pfizer_rule.passed,
            logp=result.pfizer_rule.logp,
            tpsa=result.pfizer_rule.tpsa,
            interpretation=result.pfizer_rule.interpretation,
        )

    # GSK 4/400 Rule
    gsk_rule = None
    if result.gsk_rule:
        gsk_rule = GSKRuleSchema(
            passed=result.gsk_rule.passed,
            mw=result.gsk_rule.mw,
            logp=result.gsk_rule.logp,
            interpretation=result.gsk_rule.interpretation,
        )

    # Golden Triangle
    golden_triangle = None
    if result.golden_triangle:
        golden_triangle = GoldenTriangleSchema(
            in_golden_triangle=result.golden_triangle.in_golden_triangle,
            mw=result.golden_triangle.mw,
            logd=result.golden_triangle.logd,
            interpretation=result.golden_triangle.interpretation,
        )

    return ADMETResultSchema(
        synthetic_accessibility=synthetic_accessibility,
        solubility=solubility,
        complexity=complexity,
        cns_mpo=cns_mpo,
        bioavailability=bioavailability,
        pfizer_rule=pfizer_rule,
        gsk_rule=gsk_rule,
        golden_triangle=golden_triangle,
        molar_refractivity=result.molar_refractivity,
        interpretation=result.interpretation,
    )


def _calculate_aggregator(mol: Chem.Mol) -> AggregatorLikelihoodSchema:
    """Calculate aggregator likelihood and convert to schema."""
    result = calculate_aggregator_likelihood(mol)

    return AggregatorLikelihoodSchema(
        likelihood=result.likelihood,
        risk_score=result.risk_score,
        logp=result.logp,
        tpsa=result.tpsa,
        mw=result.mw,
        aromatic_rings=result.aromatic_rings,
        risk_factors=result.risk_factors,
        interpretation=result.interpretation,
        confidence=result.confidence,
        evidence=result.evidence,
    )


def _calculate_consensus(mol: Chem.Mol) -> ConsensusScoreSchema:
    """Calculate consensus drug-likeness score and convert to schema."""
    result = calculate_consensus(mol)

    rule_sets = [
        RuleSetDetailSchema(
            name=rs.name,
            passed=rs.passed,
            violations=[
                RuleViolationSchema(
                    property=v.property,
                    value=v.value,
                    threshold=v.threshold,
                    result=v.result,
                )
                for v in rs.violations
            ],
        )
        for rs in result.rule_sets
    ]

    return ConsensusScoreSchema(
        score=result.score,
        total=result.total,
        rule_sets=rule_sets,
        interpretation=result.interpretation,
    )


def _calculate_lead_likeness(mol: Chem.Mol) -> LeadLikenessSchema:
    """Calculate lead-likeness and convert to schema."""
    result = calculate_lead_likeness(mol)

    return LeadLikenessSchema(
        passed=result.passed,
        violations=result.violations,
        properties=result.properties,
        thresholds=result.thresholds,
        violation_details=[
            RuleViolationSchema(
                property=v.property,
                value=v.value,
                threshold=v.threshold,
                result=v.result,
            )
            for v in result.violation_details
        ],
    )


def _calculate_salt_inventory(mol: Chem.Mol) -> SaltInventorySchema:
    """Calculate salt inventory and convert to schema."""
    result = calculate_salt_inventory(mol)

    return SaltInventorySchema(
        has_salts=result.has_salts,
        parent_smiles=result.parent_smiles,
        fragments=[
            SaltFragmentSchema(
                smiles=f.smiles,
                name=f.name,
                category=f.category,
                mw=f.mw,
                heavy_atom_count=f.heavy_atom_count,
            )
            for f in result.fragments
        ],
        total_fragments=result.total_fragments,
        interpretation=result.interpretation,
    )


def _calculate_ligand_efficiency(mol: Chem.Mol) -> LigandEfficiencySchema:
    """Calculate ligand efficiency and convert to schema."""
    result = calculate_ligand_efficiency(mol)

    return LigandEfficiencySchema(
        le=result.le,
        heavy_atom_count=result.heavy_atom_count,
        activity_value=result.activity_value,
        activity_type=result.activity_type,
        proxy_used=result.proxy_used,
        interpretation=result.interpretation,
    )


def _calculate_tpsa_breakdown(mol: Chem.Mol) -> TPSABreakdownSchema:
    """Calculate TPSA breakdown and convert to schema."""
    result = calculate_tpsa_breakdown(mol)

    return TPSABreakdownSchema(
        total=result.total,
        atom_contributions=[
            AtomContributionSchema(
                atom_index=ac.atom_index,
                symbol=ac.symbol,
                contribution=ac.contribution,
            )
            for ac in result.atom_contributions
        ],
        functional_group_summary=[
            FunctionalGroupContributionSchema(
                group_name=fg.group_name,
                contribution=fg.contribution,
                atom_indices=fg.atom_indices,
            )
            for fg in result.functional_group_summary
        ],
    )


def _calculate_logp_breakdown(mol: Chem.Mol) -> LogPBreakdownSchema:
    """Calculate LogP breakdown and convert to schema."""
    result = calculate_logp_breakdown(mol)

    return LogPBreakdownSchema(
        total=result.total,
        atom_contributions=[
            AtomContributionSchema(
                atom_index=ac.atom_index,
                symbol=ac.symbol,
                contribution=ac.contribution,
            )
            for ac in result.atom_contributions
        ],
        functional_group_summary=[
            FunctionalGroupContributionSchema(
                group_name=fg.group_name,
                contribution=fg.contribution,
                atom_indices=fg.atom_indices,
            )
            for fg in result.functional_group_summary
        ],
    )


def _calculate_bertz_detail(mol: Chem.Mol) -> BertzDetailSchema:
    """Calculate Bertz detail and convert to schema."""
    result = calculate_bertz_detail(mol)

    return BertzDetailSchema(
        bertz_ct=result.bertz_ct,
        num_bonds=result.num_bonds,
        num_atoms=result.num_atoms,
        num_rings=result.num_rings,
        num_aromatic_rings=result.num_aromatic_rings,
        ring_complexity=result.ring_complexity,
        interpretation=result.interpretation,
    )


def _calculate_fsp3_detail(mol: Chem.Mol) -> Fsp3DetailSchema:
    """Calculate Fsp3 detail and convert to schema."""
    result = calculate_fsp3_detail(mol)

    return Fsp3DetailSchema(
        fsp3=result.fsp3,
        total_carbons=result.total_carbons,
        sp3_count=result.sp3_count,
        sp2_count=result.sp2_count,
        sp_count=result.sp_count,
        per_carbon=[
            CarbonHybridizationSchema(
                atom_index=ch.atom_index,
                symbol=ch.symbol,
                hybridization=ch.hybridization,
            )
            for ch in result.per_carbon
        ],
        interpretation=result.interpretation,
    )


def _calculate_np_breakdown_route(mol: Chem.Mol) -> NPBreakdownSchema:
    """Calculate NP-likeness fragment breakdown and convert to schema."""
    result = calculate_np_breakdown(mol)

    return NPBreakdownSchema(
        score=result.score,
        confidence=result.confidence,
        fragments=[
            NPFragmentSchema(
                smiles=f.smiles,
                contribution=f.contribution,
                bit_id=f.bit_id,
                radius=f.radius,
                center_atom_idx=f.center_atom_idx,
                classification=f.classification,
            )
            for f in result.fragments
        ],
        total_fragments=result.total_fragments,
        np_fragment_count=result.np_fragment_count,
        synthetic_fragment_count=result.synthetic_fragment_count,
        interpretation=result.interpretation,
    )


def _calculate_bioavailability_radar_route(mol: Chem.Mol) -> BioavailabilityRadarSchema:
    """Calculate bioavailability radar and convert to schema."""
    result = calculate_bioavailability_radar(mol)

    return BioavailabilityRadarSchema(
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
            for a in result.axes
        ],
        overall_in_range_count=result.overall_in_range_count,
        interpretation=result.interpretation,
    )


def _calculate_boiled_egg_route(mol: Chem.Mol) -> BoiledEggSchema:
    """Calculate BOILED-Egg classification and convert to schema."""
    result = calculate_boiled_egg(mol)

    gi_ellipse = None
    if result.gi_ellipse:
        gi_ellipse = EllipseParamsSchema(
            cx=result.gi_ellipse.cx,
            cy=result.gi_ellipse.cy,
            a=result.gi_ellipse.a,
            b=result.gi_ellipse.b,
        )

    bbb_ellipse = None
    if result.bbb_ellipse:
        bbb_ellipse = EllipseParamsSchema(
            cx=result.bbb_ellipse.cx,
            cy=result.bbb_ellipse.cy,
            a=result.bbb_ellipse.a,
            b=result.bbb_ellipse.b,
        )

    return BoiledEggSchema(
        wlogp=result.wlogp,
        tpsa=result.tpsa,
        gi_absorbed=result.gi_absorbed,
        bbb_permeant=result.bbb_permeant,
        region=result.region,
        gi_ellipse=gi_ellipse,
        bbb_ellipse=bbb_ellipse,
        interpretation=result.interpretation,
    )


def _ml_dimension_to_schema(dim) -> MLDimensionSchema:
    """Convert an MLDimension dataclass to its Pydantic schema."""
    return MLDimensionSchema(
        name=dim.name,
        score=dim.score,
        max_score=dim.max_score,
        items=[
            MLDimensionItemSchema(
                name=item.name,
                score=item.score,
                max_score=item.max_score,
                passed=item.passed,
                subtitle=item.subtitle,
                tooltip=item.tooltip,
            )
            for item in dim.items
        ],
        details=dim.details,
    )


