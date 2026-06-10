"""Scoring schema builders.

Pure ``mol -> schema`` construction helpers extracted from the scoring API
routes to keep the route handlers thin. Each builder wraps a service-layer
calculator and maps its result onto the response schema.
"""

from __future__ import annotations

from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from app.services.scoring.safety_filters import FilterResult

from rdkit import Chem

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
    MueggeSchema,
    NPBreakdownSchema,
    NPFragmentSchema,
    PfizerRuleSchema,
    QEDSchema,
    RadarAxisSchema,
    RuleOfThreeSchema,
    RuleSetDetailSchema,
    RuleViolationSchema,
    SafetyFilterResultSchema,
    SaltFragmentSchema,
    SaltInventorySchema,
    SolubilitySchema,
    SyntheticAccessibilitySchema,
    TPSABreakdownSchema,
    VeberSchema,
)
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
    calculate_np_breakdown,
    calculate_safety_filters,
    calculate_salt_inventory,
    calculate_tpsa_breakdown,
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


