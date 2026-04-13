/**
 * Profiler Types
 *
 * TypeScript interfaces for the Compound Profiler API responses.
 * All shapes match backend/app/schemas/profiler.py.
 */

// =============================================================================
// PFI (Property Forecast Index)
// =============================================================================

/**
 * PFI result with risk classification.
 * PFI = cLogP + number of aromatic rings. Low PFI → better developability.
 */
export interface PFIResult {
  pfi: number;
  clogp: number;
  aromatic_rings: number;
  risk: 'low' | 'moderate' | 'high';
}

// =============================================================================
// #Stars (QikProp-equivalent outlier count)
// =============================================================================

/**
 * A single property detail for #stars evaluation.
 * Shows whether a property falls within the acceptable range.
 */
export interface StarDetail {
  property: string;
  value: number;
  range_low: number;
  range_high: number;
  violated: boolean;
}

/**
 * #Stars result. Each violated property contributes one star.
 * Lower star count indicates a better overall property profile.
 */
export interface StarsResult {
  stars: number;
  details: StarDetail[];
}

// =============================================================================
// Abbott Bioavailability Score
// =============================================================================

/**
 * Abbott Bioavailability Score result.
 * Estimates probability of oral bioavailability using a 4-class model.
 */
export interface AbbottResult {
  abbott_score: number;
  probability_pct: number;
  tpsa: number;
  lipinski_violations: number;
}

// =============================================================================
// Consensus LogP
// =============================================================================

/**
 * Consensus LogP result combining multiple calculation methods.
 * Uses Wildman-Crippen and XLOGP3 to produce a consensus estimate.
 */
export interface ConsensusLogPResult {
  consensus_logp: number;
  wildman_crippen: number;
  xlogp3_approx: number;
  xlogp3_is_approximation: boolean;
}

// =============================================================================
// Skin Permeation
// =============================================================================

/**
 * Skin permeation result using the Potts-Guy model.
 * log Kp predicts dermal penetration rate.
 */
export interface SkinPermeationResult {
  log_kp: number;
  classification: 'low' | 'moderate' | 'high';
}

// =============================================================================
// 3D Shape Descriptors
// =============================================================================

/**
 * 3D shape descriptor result using PMI/NPR/PBF metrics.
 * Requires conformer generation — may fail for complex molecules.
 */
export interface Shape3DResult {
  pmi1: number;
  pmi2: number;
  pmi3: number;
  npr1: number;
  npr2: number;
  pbf: number;
  shape_class: 'rod' | 'disc' | 'sphere';
  '3d_conformer_failed': boolean;
}

// =============================================================================
// SA Score Comparison
// =============================================================================

/**
 * Detail for a single synthetic accessibility scoring method.
 */
export interface SAScoreDetail {
  score: number;
  scale: string;
  classification: 'easy' | 'moderate' | 'difficult';
  available: boolean;
  error?: string;
}

/**
 * SA Comparison result combining SA Score, SCScore, SYBA, and RAscore slot.
 * RAscore is not available (incompatible Python/TF version requirements).
 */
export interface SAComparisonResult {
  sa_score: SAScoreDetail;
  scscore: SAScoreDetail;
  syba: SAScoreDetail;
  rascore: { available: false; note: string };
}

// =============================================================================
// CNS MPO Score
// =============================================================================

/**
 * CNS MPO (Multi-Parameter Optimization) score result.
 * Uses Wager 2010 definition with 6 CNS-relevant properties.
 * Score out of max_score (typically 6).
 */
export interface CNSMPOResult {
  score: number;
  max_score: number;
  components: Record<string, number>;
}

// =============================================================================
// Custom MPO
// =============================================================================

/**
 * A single component in a custom MPO score computation.
 */
export interface CustomMPOComponent {
  property: string;
  raw_value: number;
  desirability: number;
  weight: number;
}

/**
 * Custom MPO score result.
 * Normalized score ranges from 0 to 1.
 */
export interface CustomMPOResult {
  score: number;
  max_score: number;
  normalized: number;
  components: CustomMPOComponent[];
}

// =============================================================================
// Ligand Efficiency
// =============================================================================

/**
 * Extended ligand efficiency metrics.
 * Requires an activity value (IC50/Ki/pIC50/pKd) to compute.
 */
export interface LEResult {
  pIC50: number;
  LE: number;
  LLE: number;
  LELP: number;
  BEI: number;
  SEI: number;
}

// =============================================================================
// MPO Property Configuration
// =============================================================================

/**
 * A single property configuration for custom MPO scoring.
 * Defines the desirability function shape and acceptable range.
 */
export interface MPOProperty {
  property: string;
  low: number;
  high: number;
  weight: number;
  shape: 'sigmoid' | 'ramp' | 'step';
}

// =============================================================================
// Full Profile Response
// =============================================================================

/**
 * Complete profiler response from /api/v1/profiler/full.
 * Contains results for all computed profiling metrics.
 */
export interface ProfileResponse {
  pfi: PFIResult;
  stars: StarsResult;
  abbott: AbbottResult;
  consensus_logp: ConsensusLogPResult;
  skin_permeation: SkinPermeationResult;
  sa_comparison: SAComparisonResult;
  cns_mpo: CNSMPOResult;
  druglikeness?: Record<string, unknown>;
}

// =============================================================================
// Identifier Resolution
// =============================================================================

/**
 * Response from /api/v1/resolve for identifier-to-SMILES conversion.
 */
export interface ResolveResponse {
  smiles: string;
  identifier: string;
  source: string;
}
