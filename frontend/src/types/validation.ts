export type Severity = 'critical' | 'error' | 'warning' | 'info' | 'pass';

export interface CheckResult {
  check_name: string;
  passed: boolean;
  severity: Severity;
  message: string;
  affected_atoms: number[];
  details: Record<string, unknown>;
}

export interface MoleculeInfo {
  input_smiles: string;
  canonical_smiles: string | null;
  canonical_smiles_source?: string | null;
  inchi: string | null;
  inchikey: string | null;
  molecular_formula: string | null;
  molecular_weight: number | null;
  num_atoms: number | null;
  num_bonds: number | null;
  num_rings: number | null;
  num_aromatic_rings: number | null;
  num_stereocenters: number | null;
  has_stereochemistry: boolean | null;
}

export interface ValidationRequest {
  molecule: string;
  format?: 'auto' | 'smiles' | 'inchi' | 'mol';
  checks?: string[];
  preserve_aromatic?: boolean;
  input_type?: 'auto' | 'smiles' | 'iupac';
}

export interface InputInterpretation {
  detected_as: 'smiles' | 'iupac';
  original_input: string;
  converted_smiles: string | null;
  conversion_source: string | null;
}

export interface ValidationResponse {
  status: string;
  molecule_info: MoleculeInfo;
  overall_score: number;
  issues: CheckResult[];
  all_checks: CheckResult[];
  execution_time_ms: number;
  input_interpretation?: InputInterpretation;
}

export interface ValidationError {
  error: string;
  details?: {
    errors?: string[];
    warnings?: string[];
    format_detected?: string;
  };
}

export interface ChecksResponse {
  [category: string]: string[];
}

// === Deep Validation Detail Types ===

// DVAL-01/02: Stereoisomer enumeration
export interface StereoisomerDetail {
  undefined_count: number;
  total_centers: number;
  atom_indices: number[];
  stereoisomer_smiles: string[];
  enumeration_cap: number;
  cap_exceeded: boolean;
}

// DVAL-03: Tautomer detection
export interface TautomerDetail {
  tautomer_count: number;
  canonical_smiles: string;
  is_canonical_form: boolean;
  tautomer_smiles: string[];
}

// DVAL-04: Aromatic system validation
export interface AromaticSystemDetail {
  unusual_ring_sizes: Array<{ ring_size: number; atom_indices: number[] }>;
  charged_aromatics: Array<{ atom_idx: number; symbol: string; charge: number }>;
}

// DVAL-05: Coordinate dimension
export interface CoordinateDimensionDetail {
  dimension: '2d' | '3d' | 'no_coordinates' | 'degenerate';
  num_conformers: number;
}

// DVAL-06: Mixture detection
export interface FragmentDetail {
  smiles: string;
  molecular_weight: number;
  heavy_atom_count: number;
  classification: 'drug' | 'salt' | 'solvent' | 'unknown';
  pattern_name: string | null;
}

export interface MixtureDetail {
  num_fragments: number;
  fragments: FragmentDetail[];
}

// DVAL-07: Solvent contamination
export interface SolventDetail {
  solvents_found: Array<{ name: string; smiles: string; molecular_weight: number }>;
  is_pure_solvent: boolean;
}

// DVAL-08: Inorganic filter
export interface InorganicDetail {
  has_carbon: boolean;
  is_inorganic: boolean;
  is_organometallic: boolean;
  metal_atoms: Array<{ atom_idx: number; symbol: string; atomic_num: number }>;
  element_counts: Record<string, number>;
}

// DVAL-09: Radical detection
export interface RadicalDetail {
  radical_atoms: Array<{ atom_idx: number; symbol: string; num_radical_electrons: number }>;
  total_radical_electrons: number;
}

// DVAL-10: Isotope label detection
export interface IsotopeDetail {
  labeled_atoms: Array<{ atom_idx: number; symbol: string; isotope_mass: number; common_name: string | null }>;
  total_labeled: number;
}

// DVAL-11: Trivial molecule
export interface TrivialMoleculeDetail {
  heavy_atom_count: number;
  num_bonds: number;
  is_single_atom: boolean;
  threshold: number;
}

// DVAL-12: Hypervalent atoms
export interface HypervalentDetail {
  hypervalent_atoms: Array<{ atom_idx: number; symbol: string; actual_valence: number; allowed_valences: number[] }>;
}

// DVAL-13: Polymer detection
export interface PolymerDetail {
  has_sgroup_markers: boolean;
  sgroup_types: string[];
  molecular_weight: number;
  exceeds_mw_threshold: boolean;
  has_dummy_atoms: boolean;
  dummy_atom_count: number;
}

// DVAL-14: Ring strain
export interface RingStrainDetail {
  strained_rings: Array<{ ring_size: number; atom_indices: number[] }>;
  total_strained_rings: number;
}

// DVAL-15: Macrocycle detection
export interface MacrocycleDetail {
  macrocycles: Array<{ ring_size: number; atom_indices: number[] }>;
  total_macrocycles: number;
  sssr_note: string;
}

// DVAL-16: Charged species
export interface ChargedSpeciesDetail {
  net_charge: number;
  positive_atoms: Array<{ atom_idx: number; symbol: string; charge: number }>;
  negative_atoms: Array<{ atom_idx: number; symbol: string; charge: number }>;
  is_zwitterion: boolean;
  total_charged_atoms: number;
}

// DVAL-17: Explicit hydrogen audit
export interface ExplicitHydrogenDetail {
  atoms_with_explicit_h: Array<{ atom_idx: number; symbol: string; explicit_h_count: number }>;
  total_explicit_h: number;
  has_h_atom_objects: boolean;
  h_atom_object_count: number;
}

// Severity override type
export type SeverityOverride = 'error' | 'warning' | 'info';

export interface DeepValidationConfig {
  severityOverrides: Record<string, SeverityOverride>;
}

// Domain grouping for deep checks
export type DeepValidationDomain = 'stereo_tautomer' | 'chemical_composition' | 'structural_complexity';

export const DEEP_CHECK_DOMAINS: Record<DeepValidationDomain, { label: string; checks: string[] }> = {
  stereo_tautomer: {
    label: 'Stereo & Tautomers',
    checks: ['stereoisomer_enumeration', 'tautomer_detection', 'aromatic_system_validation', 'coordinate_dimension'],
  },
  chemical_composition: {
    label: 'Chemical Composition',
    checks: ['mixture_detection', 'solvent_contamination', 'inorganic_filter', 'radical_detection', 'isotope_label_detection', 'trivial_molecule'],
  },
  structural_complexity: {
    label: 'Structural Complexity',
    checks: ['hypervalent_atoms', 'polymer_detection', 'ring_strain', 'macrocycle_detection', 'charged_species', 'explicit_hydrogen_audit'],
  },
};
