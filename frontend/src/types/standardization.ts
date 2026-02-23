/**
 * Standardization Types
 *
 * TypeScript interfaces for standardization requests and responses.
 */

// ---------------------------------------------------------------------------
// Provenance types (Phase 02 — STD-01 through STD-06)
// ---------------------------------------------------------------------------

export interface ChargeChange {
  atom_idx: number;
  element: string;
  before_charge: number;
  after_charge: number;
  rule_name: string;
  smarts: string;
}

export interface BondChange {
  bond_idx: number;
  atom1_idx: number;
  atom2_idx: number;
  before_type: string;
  after_type: string;
  rule_name: string;
}

export interface RadicalChange {
  atom_idx: number;
  element: string;
  before_radicals: number;
  after_radicals: number;
}

export interface RingChange {
  ring_atoms: number[];
  ring_size: number;
  before_type: string;
  after_type: string;
}

export interface FragmentRemoval {
  smiles: string;
  name: string | null;
  role: 'salt' | 'solvent' | 'counterion' | 'unknown';
  mw: number;
}

export interface TautomerProvenance {
  input_smiles: string;
  canonical_smiles: string;
  num_tautomers_enumerated: number;
  modified_atoms: number[];
  modified_bonds: number[];
  stereo_stripped: boolean;
  complexity_flag: boolean;
}

export interface StereoCenterDetail {
  atom_idx: number;
  type: string;
  before_config: string;
  after_config: string;
  reason: string;
}

export interface StereoProvenance {
  stereo_stripped: boolean;
  centers_lost: number;
  bonds_lost: number;
  per_center: StereoCenterDetail[];
}

export interface ProvStageRecord {
  stage_name: string;
  input_smiles: string;
  output_smiles: string;
  applied: boolean;
  charge_changes: ChargeChange[];
  bond_changes: BondChange[];
  radical_changes: RadicalChange[];
  ring_changes: RingChange[];
  fragment_removals: FragmentRemoval[];
  dval_cross_refs: string[];
}

export interface StandardizationProvenance {
  stages: ProvStageRecord[];
  tautomer: TautomerProvenance | null;
  stereo_summary: StereoProvenance | null;
}

// ---------------------------------------------------------------------------
// Core standardization types
// ---------------------------------------------------------------------------

export interface StandardizationOptions {
  /**
   * Include tautomer canonicalization.
   * WARNING: May lose E/Z double bond stereochemistry.
   */
  include_tautomer: boolean;

  /**
   * Attempt to preserve stereochemistry during standardization.
   */
  preserve_stereo: boolean;

  /**
   * Include detailed per-stage provenance records in the response.
   */
  include_provenance?: boolean;
}

export interface StandardizationStep {
  step_name: string;
  applied: boolean;
  description: string;
  changes: string;
}

export interface StereoComparison {
  before_count: number;
  after_count: number;
  lost: number;
  gained: number;
  double_bond_stereo_lost: number;
  warning: string | null;
}

export interface StructureComparison {
  original_atom_count: number;
  standardized_atom_count: number;
  original_formula: string | null;
  standardized_formula: string | null;
  original_mw: number | null;
  standardized_mw: number | null;
  mass_change_percent: number;
  is_identical: boolean;
  diff_summary: string[];
}

export interface CheckerIssue {
  penalty_score: number;
  message: string;
}

export interface StandardizationResult {
  original_smiles: string;
  standardized_smiles: string | null;
  success: boolean;
  error_message: string | null;
  steps_applied: StandardizationStep[];
  checker_issues: CheckerIssue[];
  excluded_fragments: string[];
  stereo_comparison: StereoComparison | null;
  structure_comparison: StructureComparison | null;
  mass_change_percent: number;
  provenance?: StandardizationProvenance | null;
}

export interface MoleculeInfo {
  input_smiles: string;
  canonical_smiles: string | null;
  inchi: string | null;
  inchikey: string | null;
  molecular_formula: string | null;
  molecular_weight: number | null;
  num_atoms: number | null;
}

export interface StandardizeRequest {
  molecule: string;
  format?: 'auto' | 'smiles' | 'inchi' | 'mol';
  options?: StandardizationOptions;
}

export interface StandardizeResponse {
  molecule_info: MoleculeInfo;
  result: StandardizationResult;
  execution_time_ms: number;
}

export interface StandardizeError {
  error: string;
  details?: {
    errors?: string[];
    warnings?: string[];
    format_detected?: string;
  };
}

export interface StandardizeOptionsResponse {
  options: {
    include_tautomer: {
      type: string;
      default: boolean;
      description: string;
      warning: string;
    };
    preserve_stereo: {
      type: string;
      default: boolean;
      description: string;
    };
  };
  pipeline_steps: {
    name: string;
    description: string;
    always_run: boolean;
    requires_option?: string;
  }[];
}
