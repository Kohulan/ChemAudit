/**
 * Scoring Types
 *
 * TypeScript interfaces for ML-readiness, NP-likeness, and scaffold scoring.
 */

/**
 * Breakdown of ML-readiness score components.
 */
export interface MLReadinessBreakdown {
  descriptors_score: number;
  descriptors_max: number;
  descriptors_successful: number;
  descriptors_total: number;

  fingerprints_score: number;
  fingerprints_max: number;
  fingerprints_successful: string[];
  fingerprints_failed: string[];

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
 * Request for molecule scoring.
 */
export interface ScoringRequest {
  molecule: string;
  format?: 'auto' | 'smiles' | 'inchi' | 'mol';
  include?: ('ml_readiness' | 'np_likeness' | 'scaffold')[];
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
