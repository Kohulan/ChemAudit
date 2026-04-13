/**
 * TypeScript interfaces for the Structure Quality Diagnostics API.
 * Covers DIAG-01 through DIAG-05 endpoints under /api/v1/diagnostics/.
 *
 * All shapes match backend Pydantic schemas in backend/app/schemas/diagnostics.py.
 */

// ====================== SMILES Diagnostics (DIAG-01) ======================

export interface FixSuggestion {
  description: string;
  corrected_smiles: string | null;
  confidence: number;
}

export interface SMILESError {
  raw_message: string;
  position: number | null; // zero-indexed character position
  error_type: string;
  message: string;
  suggestions: FixSuggestion[];
}

export interface SMILESWarning {
  atom_index: number | null;
  type: string;
  message: string;
}

export interface SMILESDiagnosticsResponse {
  valid: boolean;
  canonical_smiles: string | null;
  warnings: SMILESWarning[];
  errors: SMILESError[];
}

// ====================== InChI Layer Diff (DIAG-02) ======================

export interface LayerRow {
  layer: string;
  value_a: string | null;
  value_b: string | null;
  match: boolean;
}

export interface InChIDiffResponse {
  identical: boolean;
  layer_rows: LayerRow[];
  layers_a: Record<string, string>;
  layers_b: Record<string, string>;
}

// ====================== Round-Trip Lossiness (DIAG-03) ======================

export interface RoundTripLoss {
  type: string; // "stereo", "charge", "isotope"
  description: string;
  before: number;
  after: number;
}

export interface RoundTripResponse {
  route: string;
  original_smiles: string;
  intermediate: string | null;
  roundtrip_smiles: string | null;
  lossy: boolean;
  losses: RoundTripLoss[];
  error: string | null;
}

// ====================== Cross-Pipeline Comparison (DIAG-04) ======================

export interface PipelineResult {
  name: string;
  smiles: string;
  inchikey: string;
  mw: number;
  formula: string;
  charge: number;
  stereo_count: number;
}

export interface PropertyComparison {
  property: string;
  values: (string | number)[];
  agrees: boolean;
  structural: boolean;
}

export interface CrossPipelineResponse {
  pipelines: PipelineResult[];
  disagreements: number;
  structural_disagreements: number;
  all_agree: boolean;
  property_comparison: PropertyComparison[];
}

// ====================== File Pre-Validation (DIAG-05) ======================

export interface FileIssue {
  block: number | null;
  line: number;
  issue_type: string;
  severity: string;
  description: string;
}

export interface FilePreValidationResponse {
  file_type: string;
  total_blocks: number | null;
  total_rows: number | null;
  encoding: string | null;
  issue_count: number;
  issues: FileIssue[];
  valid: boolean;
}
