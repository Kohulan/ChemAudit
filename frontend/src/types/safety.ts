/**
 * Safety Types
 *
 * TypeScript interfaces for the Enhanced Structural Alerts & Safety API.
 * All shapes match backend/app/schemas/alerts.py and backend/app/schemas/safety.py.
 *
 * NOTE: These interfaces cover the new Phase 08 endpoints (/alerts/screen,
 * /safety/assess, /safety/summary) and are distinct from the existing
 * alerts.ts types which cover the legacy /alerts endpoint.
 */

// ====================== Shared Types ======================

export interface SafetyMoleculeInfo {
  smiles: string;
  inchi?: string;
  inchi_key?: string;
  molecular_formula?: string;
  molecular_weight?: number;
  heavy_atom_count?: number;
}

// ====================== Alert Screen Types ======================

export interface AlertResult {
  pattern_name: string;
  description: string;
  severity: 'critical' | 'warning' | 'info';
  matched_atoms: number[];
  catalog_source: string;
  smarts?: string;
  reference?: string;
  scope?: string;
  filter_set?: string;
  catalog_description?: string;
  category?: string;
  concern_group?: string;
}

export interface ConcernGroup {
  name: string;
  count: number;
  severity: string;
  alerts: AlertResult[];
}

export interface AlertScreenResponse {
  status: string;
  molecule_info: SafetyMoleculeInfo;
  alerts: AlertResult[];
  concern_groups: Record<string, ConcernGroup>;
  total_raw: number;
  total_deduped: number;
  screened_catalogs: string[];
  has_critical: boolean;
  has_warning: boolean;
  execution_time_ms: number;
}

// ====================== Safety Assess Types ======================

export interface CypSite {
  site_name: string;
  reaction_type: string;
  matched_atoms: number[];
}

export interface CypResult {
  sites: CypSite[];
  n_sites: number;
}

export interface HergResult {
  herg_risk: 'low' | 'moderate' | 'high';
  risk_score: number;
  max_score: number;
  flags: string[];
  descriptors: Record<string, number | boolean>;
}

export interface Bro5Violation {
  property: string;
  value: number;
  threshold: number;
  direction: string;
}

export interface Bro5Result {
  applicable: boolean;
  passed: boolean;
  message?: string;
  violations: Bro5Violation[];
  values: Record<string, number>;
}

export interface ReosViolation {
  property: string;
  value: number;
  range: [number, number];
  exceeded: boolean;
}

export interface ReosResult {
  passed: boolean;
  violations: ReosViolation[];
  n_violations: number;
  descriptors: Record<string, number>;
}

export interface ComplexityProperty {
  value: number;
  p5: number;
  p95: number;
  outlier: boolean;
  direction: string | null;
}

export interface ComplexityResult {
  properties: Record<string, ComplexityProperty>;
  n_outliers: number;
  outlier_properties: string[];
  within_range: boolean;
}

export interface SafetyAssessResponse {
  status: string;
  molecule_info: SafetyMoleculeInfo;
  cyp_softspots: CypResult;
  herg: HergResult;
  bro5: Bro5Result;
  reos: ReosResult;
  complexity: ComplexityResult;
  execution_time_ms: number;
}

// ====================== Safety Summary Types ======================

export interface SafetySummaryResponse {
  status: string;
  total_alerts: number;
  has_critical: boolean;
  cyp_status: 'default' | 'success' | 'warning' | 'error';
  herg_status: 'success' | 'warning' | 'error';
  bro5_status: 'success' | 'error' | 'default';
  reos_status: 'success' | 'warning' | 'error';
  complexity_outliers: number;
}
