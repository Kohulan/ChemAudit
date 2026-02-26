/**
 * Phase 6: Export, API & Workflow TypeScript Interfaces
 *
 * Types for scoring profiles, bookmarks, audit trail, permalinks,
 * IUPAC input detection, and enhanced export formats.
 */

// =============================================================================
// Scoring Profiles
// =============================================================================

export interface ThresholdRange {
  min?: number;
  max?: number;
}

export interface ScoringProfile {
  id: number;
  name: string;
  description: string | null;
  thresholds: Record<string, ThresholdRange>;
  weights: Record<string, number>;
  is_preset: boolean;
  is_active: boolean;
  created_at: string;
  updated_at: string;
}

export interface ScoringProfileCreate {
  name: string;
  description?: string | null;
  thresholds: Record<string, ThresholdRange>;
  weights: Record<string, number>;
}

export interface ScoringProfileUpdate {
  name?: string;
  description?: string | null;
  thresholds?: Record<string, ThresholdRange>;
  weights?: Record<string, number>;
}

export interface ScoringProfileExport {
  name: string;
  description: string | null;
  thresholds: Record<string, ThresholdRange>;
  weights: Record<string, number>;
}

// =============================================================================
// Bookmarks
// =============================================================================

export interface Bookmark {
  id: number;
  smiles: string;
  name: string | null;
  inchikey: string | null;
  tags: string[];
  notes: string | null;
  source: string | null;
  job_id: string | null;
  created_at: string;
}

export interface BookmarkCreate {
  smiles: string;
  name?: string | null;
  tags?: string[];
  notes?: string | null;
  source?: string | null;
  job_id?: string | null;
}

export interface BookmarkUpdate {
  name?: string | null;
  tags?: string[];
  notes?: string | null;
}

// =============================================================================
// Audit Trail / History
// =============================================================================

export interface AuditEntry {
  id: number;
  smiles: string;
  inchikey: string | null;
  outcome: 'pass' | 'warn' | 'fail';
  score: number | null;
  job_id: string | null;
  molecule_count: number | null;
  pass_count: number | null;
  fail_count: number | null;
  source: 'single' | 'batch';
  created_at: string;
}

export interface AuditHistoryResponse {
  entries: AuditEntry[];
  total: number;
  page: number;
  page_size: number;
}

export interface AuditHistoryParams {
  page?: number;
  page_size?: number;
  date_from?: string;
  date_to?: string;
  outcome?: string;
  source?: string;
  smiles_search?: string;
}

export interface AuditHistoryStats {
  total_validations: number;
  outcome_distribution: Record<string, number>;
  source_distribution: Record<string, number>;
}

// =============================================================================
// Permalinks
// =============================================================================

export interface PermalinkResponse {
  short_id: string;
  job_id: string;
  url: string;
  created_at: string;
  expires_at: string | null;
}

export interface PermalinkResolveResponse {
  job_id: string;
  snapshot_data: Record<string, unknown> | null;
  settings: Record<string, unknown> | null;
}

// =============================================================================
// IUPAC Input Detection
// =============================================================================

export interface InputInterpretation {
  detected_as: 'smiles' | 'iupac';
  original_input: string;
  converted_smiles: string | null;
  conversion_source: string | null;
}

// =============================================================================
// Export Formats
// =============================================================================

export type ExportFormat =
  | 'csv'
  | 'excel'
  | 'sdf'
  | 'json'
  | 'pdf'
  | 'fingerprint'
  | 'dedup'
  | 'scaffold'
  | 'property_matrix';

export interface PDFSectionOption {
  id: string;
  label: string;
  description: string;
  defaultChecked: boolean;
}

export const PDF_SECTION_OPTIONS: PDFSectionOption[] = [
  { id: 'validation_summary', label: 'Validation Summary', description: 'Overall validation statistics and pass/fail counts', defaultChecked: true },
  { id: 'score_distribution', label: 'Score Distribution', description: 'Histogram of validation scores across all molecules', defaultChecked: true },
  { id: 'alert_frequency', label: 'Alert Frequency', description: 'Bar chart of structural alert catalog hits', defaultChecked: true },
  { id: 'chemical_space', label: 'Chemical Space', description: 'PCA/t-SNE chemical space projection', defaultChecked: true },
  { id: 'scaffold_treemap', label: 'Scaffold Treemap', description: 'Murcko scaffold frequency distribution', defaultChecked: true },
  { id: 'statistics', label: 'Statistics', description: 'Property distribution statistics and correlations', defaultChecked: true },
  { id: 'correlation_matrix', label: 'Correlation Matrix', description: 'Property correlation heatmap', defaultChecked: true },
  { id: 'mmp_pairs', label: 'MMP Pairs', description: 'Matched molecular pair analysis results', defaultChecked: true },
];

export interface ExportFormatInfo {
  value: ExportFormat;
  label: string;
  description: string;
  extension: string;
  isNew?: boolean;
}

export const EXPORT_FORMATS: ExportFormatInfo[] = [
  { value: 'csv', label: 'CSV', description: 'Comma-separated values with all validation data', extension: 'csv' },
  { value: 'excel', label: 'Excel', description: 'Formatted spreadsheet with conditional coloring and summary', extension: 'xlsx' },
  { value: 'sdf', label: 'SDF', description: 'Structure-data file with properties attached to molecules', extension: 'sdf' },
  { value: 'json', label: 'JSON', description: 'Full result objects with metadata', extension: 'json' },
  { value: 'pdf', label: 'PDF Report', description: 'Professional report with charts, statistics, and molecule images', extension: 'pdf' },
  { value: 'fingerprint', label: 'Fingerprints', description: 'Morgan/MACCS/RDKit fingerprints as CSV, npy, and npz (ZIP)', extension: 'zip', isNew: true },
  { value: 'dedup', label: 'Deduplicated', description: 'Deduplication summary and annotated molecule CSVs (ZIP)', extension: 'zip', isNew: true },
  { value: 'scaffold', label: 'Scaffold-Grouped', description: 'CSV with Murcko scaffold SMILES and scaffold group assignment', extension: 'csv', isNew: true },
  { value: 'property_matrix', label: 'Property Matrix', description: 'All computed properties in flat CSV and multi-sheet Excel (ZIP)', extension: 'zip', isNew: true },
];
