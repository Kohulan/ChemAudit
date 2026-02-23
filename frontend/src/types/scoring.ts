/**
 * Scoring Types
 *
 * TypeScript interfaces for ML-readiness, NP-likeness, scaffold scoring,
 * drug-likeness, safety filters, and ADMET predictions.
 */

/**
 * Breakdown of ML-readiness score components.
 */
export interface MLReadinessBreakdown {
  // Standard descriptors (CalcMolDescriptors - 217 descriptors)
  descriptors_score: number;
  descriptors_max: number;
  descriptors_successful: number;
  descriptors_total: number;

  // Additional descriptors (AUTOCORR2D + MQN)
  additional_descriptors_score: number;
  additional_descriptors_max: number;
  autocorr2d_successful: number;
  autocorr2d_total: number;
  mqn_successful: number;
  mqn_total: number;

  // Fingerprints (7 types)
  fingerprints_score: number;
  fingerprints_max: number;
  fingerprints_successful: string[];
  fingerprints_failed: string[];

  // Size constraints
  size_score: number;
  size_max: number;
  molecular_weight: number | null;
  num_atoms: number | null;
  size_category: 'optimal' | 'acceptable' | 'out_of_range' | 'error' | 'unknown';
}

/**
 * ML-readiness scoring result.
 */
export interface MLReadinessResult {
  score: number;
  breakdown: MLReadinessBreakdown;
  interpretation: string;
  failed_descriptors: string[];
}

/**
 * NP-likeness scoring result.
 */
export interface NPLikenessResult {
  score: number;
  interpretation: string;
  caveats: string[];
  details: Record<string, unknown>;
}

/**
 * Scaffold extraction result.
 */
export interface ScaffoldResult {
  scaffold_smiles: string;
  generic_scaffold_smiles: string;
  has_scaffold: boolean;
  message: string;
  details: Record<string, unknown>;
}

// =============================================================================
// Drug-likeness Types
// =============================================================================

/**
 * Lipinski's Rule of Five results.
 */
export interface LipinskiResult {
  passed: boolean;
  violations: number;
  mw: number;
  logp: number;
  hbd: number;
  hba: number;
  details: Record<string, boolean>;
}

/**
 * QED score results.
 */
export interface QEDResult {
  score: number;
  properties: Record<string, number | string>;
  interpretation: string;
}

/**
 * Veber rules results.
 */
export interface VeberResult {
  passed: boolean;
  rotatable_bonds: number;
  tpsa: number;
}

/**
 * Rule of Three (fragment-likeness) results.
 */
export interface RuleOfThreeResult {
  passed: boolean;
  violations: number;
  mw: number;
  logp: number;
  hbd: number;
  hba: number;
  rotatable_bonds: number;
  tpsa: number;
}

/**
 * Ghose filter results.
 */
export interface GhoseResult {
  passed: boolean;
  violations: number;
  mw: number;
  logp: number;
  atom_count: number;
  molar_refractivity: number;
}

/**
 * Egan filter results.
 */
export interface EganResult {
  passed: boolean;
  logp: number;
  tpsa: number;
}

/**
 * Muegge filter results.
 */
export interface MueggeResult {
  passed: boolean;
  violations: number;
  details: Record<string, boolean>;
}

/**
 * Complete drug-likeness scoring results.
 */
export interface DrugLikenessResult {
  lipinski: LipinskiResult;
  qed: QEDResult;
  veber: VeberResult;
  ro3: RuleOfThreeResult;
  ghose: GhoseResult | null;
  egan: EganResult | null;
  muegge: MueggeResult | null;
  interpretation: string;
}

// =============================================================================
// Safety Filters Types
// =============================================================================

/**
 * Result for a single filter category.
 */
export interface FilterAlertResult {
  passed: boolean;
  alerts: string[];
  alert_count: number;
}

/**
 * ChEMBL structural alerts results.
 */
export interface ChEMBLAlertsResult {
  passed: boolean;
  total_alerts: number;
  bms: FilterAlertResult | null;
  dundee: FilterAlertResult | null;
  glaxo: FilterAlertResult | null;
  inpharmatica: FilterAlertResult | null;
  lint: FilterAlertResult | null;
  mlsmr: FilterAlertResult | null;
  schembl: FilterAlertResult | null;
}

/**
 * Complete safety filter results.
 */
export interface SafetyFilterResult {
  pains: FilterAlertResult;
  brenk: FilterAlertResult;
  nih: FilterAlertResult | null;
  zinc: FilterAlertResult | null;
  chembl: ChEMBLAlertsResult | null;
  all_passed: boolean;
  total_alerts: number;
  interpretation: string;
}

// =============================================================================
// ADMET Types
// =============================================================================

/**
 * Synthetic accessibility score result.
 */
export interface SyntheticAccessibilityResult {
  score: number;
  classification: 'easy' | 'moderate' | 'difficult' | 'unknown';
  interpretation: string;
}

/**
 * ESOL solubility prediction result.
 */
export interface SolubilityResult {
  log_s: number;
  solubility_mg_ml: number;
  classification: 'highly_soluble' | 'soluble' | 'moderate' | 'poor' | 'insoluble' | 'unknown';
  interpretation: string;
}

/**
 * Molecular complexity metrics.
 */
export interface ComplexityResult {
  fsp3: number;
  num_stereocenters: number;
  num_rings: number;
  num_aromatic_rings: number;
  bertz_ct: number;
  classification: 'flat' | 'moderate' | '3d' | 'unknown';
  interpretation: string;
}

/**
 * CNS MPO score result.
 */
export interface CNSMPOResult {
  score: number;
  components: Record<string, number>;
  cns_penetrant: boolean;
  interpretation: string;
}

/**
 * Bioavailability indicators.
 */
export interface BioavailabilityResult {
  tpsa: number;
  rotatable_bonds: number;
  hbd: number;
  hba: number;
  mw: number;
  logp: number;
  oral_absorption_likely: boolean;
  cns_penetration_likely: boolean;
  interpretation: string;
}

/**
 * Pfizer 3/75 Rule result.
 */
export interface PfizerRuleResult {
  passed: boolean;
  logp: number;
  tpsa: number;
  interpretation: string;
}

/**
 * GSK 4/400 Rule result.
 */
export interface GSKRuleResult {
  passed: boolean;
  mw: number;
  logp: number;
  interpretation: string;
}

/**
 * Golden Triangle (Abbott) analysis.
 */
export interface GoldenTriangleResult {
  in_golden_triangle: boolean;
  mw: number;
  logd: number;
  interpretation: string;
}

/**
 * Complete ADMET prediction results.
 */
export interface ADMETResult {
  synthetic_accessibility: SyntheticAccessibilityResult;
  solubility: SolubilityResult;
  complexity: ComplexityResult;
  cns_mpo: CNSMPOResult | null;
  bioavailability: BioavailabilityResult;
  pfizer_rule: PfizerRuleResult | null;
  gsk_rule: GSKRuleResult | null;
  golden_triangle: GoldenTriangleResult | null;
  molar_refractivity: number | null;
  interpretation: string;
}

// =============================================================================
// Aggregator Likelihood Types
// =============================================================================

/**
 * Aggregator likelihood prediction result.
 */
export interface AggregatorLikelihoodResult {
  likelihood: 'low' | 'moderate' | 'high';
  risk_score: number;
  logp: number;
  tpsa: number;
  mw: number;
  aromatic_rings: number;
  risk_factors: string[];
  interpretation: string;
  confidence: number;
  evidence: Array<Record<string, unknown>>;
}

// =============================================================================
// Consensus Drug-Likeness Types
// =============================================================================

/**
 * Per-property result within a rule set.
 */
export interface RuleViolation {
  property: string;
  value: number;
  threshold: string;
  result: 'pass' | 'fail';
}

/**
 * Detail for a single rule set in consensus scoring.
 */
export interface RuleSetDetail {
  name: string;
  passed: boolean;
  violations: RuleViolation[];
}

/**
 * Consensus drug-likeness score across 5 rule sets.
 */
export interface ConsensusScore {
  score: number;
  total: number;
  rule_sets: RuleSetDetail[];
  interpretation: string;
}

// =============================================================================
// Lead-Likeness Types
// =============================================================================

/**
 * Lead-likeness assessment result.
 */
export interface LeadLikeness {
  passed: boolean;
  violations: number;
  properties: Record<string, number>;
  thresholds: Record<string, string>;
  violation_details: RuleViolation[];
}

// =============================================================================
// Salt Inventory Types
// =============================================================================

/**
 * A single fragment from a salt-form molecule.
 */
export interface SaltFragment {
  smiles: string;
  name: string;
  category: string;
  mw: number;
  heavy_atom_count: number;
}

/**
 * Salt/counterion inventory result.
 */
export interface SaltInventory {
  has_salts: boolean;
  parent_smiles: string;
  fragments: SaltFragment[];
  total_fragments: number;
  interpretation: string;
}

// =============================================================================
// Ligand Efficiency Types
// =============================================================================

/**
 * Ligand efficiency result.
 */
export interface LigandEfficiency {
  le: number | null;
  heavy_atom_count: number;
  activity_value: number | null;
  activity_type: string | null;
  proxy_used: boolean;
  interpretation: string;
}

// =============================================================================
// Property Breakdown Types
// =============================================================================

/**
 * Per-atom property contribution.
 */
export interface AtomContribution {
  atom_index: number;
  symbol: string;
  contribution: number;
}

/**
 * Per-functional-group property contribution.
 */
export interface FunctionalGroupContribution {
  group_name: string;
  contribution: number;
  atom_indices: number[];
}

/**
 * TPSA per-atom breakdown.
 */
export interface TPSABreakdown {
  total: number;
  atom_contributions: AtomContribution[];
  functional_group_summary: FunctionalGroupContribution[];
}

/**
 * LogP per-atom breakdown.
 */
export interface LogPBreakdown {
  total: number;
  atom_contributions: AtomContribution[];
  functional_group_summary: FunctionalGroupContribution[];
}

/**
 * Bertz complexity detail.
 */
export interface BertzDetail {
  bertz_ct: number;
  num_bonds: number;
  num_atoms: number;
  num_rings: number;
  num_aromatic_rings: number;
  ring_complexity: number;
  interpretation: string;
}

/**
 * Per-carbon hybridization data.
 */
export interface CarbonHybridization {
  atom_index: number;
  symbol: string;
  hybridization: string;
}

/**
 * Fsp3 per-carbon detail.
 */
export interface Fsp3Detail {
  fsp3: number;
  total_carbons: number;
  sp3_count: number;
  sp2_count: number;
  sp_count: number;
  per_carbon: CarbonHybridization[];
  interpretation: string;
}

// =============================================================================
// NP-Likeness Breakdown Types
// =============================================================================

/**
 * A single NP fragment contribution.
 */
export interface NPFragment {
  smiles: string;
  contribution: number;
  bit_id: number;
  radius: number;
  center_atom_idx: number;
  classification: 'np_characteristic' | 'synthetic_characteristic' | 'neutral';
}

/**
 * NP-likeness fragment breakdown.
 */
export interface NPBreakdown {
  score: number;
  confidence: number;
  fragments: NPFragment[];
  total_fragments: number;
  np_fragment_count: number;
  synthetic_fragment_count: number;
  interpretation: string;
}

// =============================================================================
// Bioavailability Radar & BOILED-Egg Types
// =============================================================================

/**
 * A single axis of the bioavailability radar.
 */
export interface RadarAxis {
  name: string;
  actual_value: number;
  normalized: number;
  optimal_min: number;
  optimal_max: number;
  in_range: boolean;
  property_name: string;
  unit: string;
}

/**
 * Bioavailability radar result.
 */
export interface BioavailabilityRadar {
  axes: RadarAxis[];
  overall_in_range_count: number;
  interpretation: string;
}

/**
 * Ellipse parameters for BOILED-Egg model.
 */
export interface EllipseParams {
  cx: number;
  cy: number;
  a: number;
  b: number;
}

/**
 * BOILED-Egg classification result.
 */
export interface BoiledEgg {
  wlogp: number;
  tpsa: number;
  gi_absorbed: boolean;
  bbb_permeant: boolean;
  region: 'yolk' | 'white' | 'grey';
  gi_ellipse: EllipseParams | null;
  bbb_ellipse: EllipseParams | null;
  interpretation: string;
}

/**
 * A single molecule's radar profile.
 */
export interface RadarProfile {
  smiles: string;
  axes: RadarAxis[];
  is_reference: boolean;
}

/**
 * Multi-molecule radar comparison.
 */
export interface RadarComparison {
  profiles: RadarProfile[];
  reference: RadarProfile | null;
}

// =============================================================================
// Core Types
// =============================================================================

/**
 * Basic molecule information.
 */
export interface ScoringMoleculeInfo {
  input_string: string;
  canonical_smiles: string | null;
  molecular_formula: string | null;
  molecular_weight: number | null;
}

/**
 * Available scoring types.
 */
export type ScoringType =
  | 'ml_readiness'
  | 'np_likeness'
  | 'scaffold'
  | 'druglikeness'
  | 'safety_filters'
  | 'admet'
  | 'aggregator'
  | 'consensus'
  | 'lead_likeness'
  | 'salt_inventory'
  | 'ligand_efficiency'
  | 'tpsa_breakdown'
  | 'logp_breakdown'
  | 'bertz_detail'
  | 'fsp3_detail'
  | 'np_breakdown'
  | 'bioavailability_radar'
  | 'boiled_egg';

/**
 * Request for molecule scoring.
 */
export interface ScoringRequest {
  molecule: string;
  format?: 'auto' | 'smiles' | 'inchi' | 'mol';
  include?: ScoringType[];
}

/**
 * Response containing scoring results.
 */
export interface ScoringResponse {
  status: string;
  molecule_info: ScoringMoleculeInfo;
  ml_readiness: MLReadinessResult | null;
  np_likeness: NPLikenessResult | null;
  scaffold: ScaffoldResult | null;
  druglikeness: DrugLikenessResult | null;
  safety_filters: SafetyFilterResult | null;
  admet: ADMETResult | null;
  aggregator: AggregatorLikelihoodResult | null;
  consensus: ConsensusScore | null;
  lead_likeness: LeadLikeness | null;
  salt_inventory: SaltInventory | null;
  ligand_efficiency: LigandEfficiency | null;
  tpsa_breakdown: TPSABreakdown | null;
  logp_breakdown: LogPBreakdown | null;
  bertz_detail: BertzDetail | null;
  fsp3_detail: Fsp3Detail | null;
  np_breakdown: NPBreakdown | null;
  bioavailability_radar: BioavailabilityRadar | null;
  boiled_egg: BoiledEgg | null;
  execution_time_ms: number;
}

/**
 * Scoring error response.
 */
export interface ScoringError {
  error: string;
  details?: {
    errors?: string[];
    warnings?: string[];
    format_detected?: string;
  };
}
