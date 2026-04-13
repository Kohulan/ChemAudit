/**
 * Types for external database integrations (PubChem, ChEMBL, COCONUT, Wikidata, SureChEMBL)
 */

// Request types
export interface IntegrationRequest {
  smiles?: string;
  inchikey?: string;
}

// PubChem
export interface PubChemResult {
  found: boolean;
  cid?: number;
  iupac_name?: string;
  molecular_formula?: string;
  molecular_weight?: number;
  canonical_smiles?: string;
  inchi?: string;
  inchikey?: string;
  synonyms?: string[];
  url?: string;
}

// ChEMBL
export interface BioactivityData {
  target_chembl_id: string;
  target_name?: string;
  target_type?: string;
  activity_type: string;
  activity_value?: number;
  activity_unit?: string;
  assay_chembl_id: string;
  document_chembl_id?: string;
}

export interface ChEMBLResult {
  found: boolean;
  chembl_id?: string;
  pref_name?: string;
  molecule_type?: string;
  max_phase?: number;
  molecular_formula?: string;
  molecular_weight?: number;
  canonical_smiles?: string;
  inchi?: string;
  inchikey?: string;
  bioactivities: BioactivityData[];
  bioactivity_count: number;
  url?: string;
}

// COCONUT
export interface COCONUTResult {
  found: boolean;
  coconut_id?: string;
  name?: string;
  smiles?: string;
  inchi?: string;
  inchikey?: string;
  molecular_formula?: string;
  molecular_weight?: number;
  organism?: string;
  organism_type?: string;
  nplikeness?: number;
  url?: string;
}

// Wikidata
export interface WikidataResult {
  found: boolean;
  label?: string;
  smiles?: string;
  inchi?: string;
  inchikey?: string;
  cas?: string;
  molecular_formula?: string;
  molecular_weight?: number;
  url?: string;
}

// SureChEMBL
export interface SureChEMBLResult {
  found: boolean;
  schembl_id?: string;
  url?: string;
  patent_count?: number;
  source?: string;
  smiles?: string;
  inchi?: string;
  inchikey?: string;
  molecular_weight?: number;
}

// Combined lookup result
export interface DatabaseLookupResult {
  pubchem?: PubChemResult;
  chembl?: ChEMBLResult;
  coconut?: COCONUTResult;
  wikidata?: WikidataResult;
  surechembl?: SureChEMBLResult;
}

export interface IntegrationError {
  error: string;
  details?: Record<string, unknown>;
}

// Universal Identifier Resolution
export interface ResolveRequest {
  identifier: string;
  identifier_type?: string;
}

export interface CrossReferences {
  pubchem_cid?: number;
  chembl_id?: string;
  coconut_id?: string;
  drugbank_id?: string;
  chebi_id?: string;
  unii?: string;
  cas?: string;
  wikipedia_url?: string;
  kegg_id?: string;
}

export interface ResolvedCompound {
  resolved: boolean;
  identifier_type_detected: string;
  canonical_smiles?: string;
  inchi?: string;
  inchikey?: string;
  molecular_formula?: string;
  molecular_weight?: number;
  iupac_name?: string;
  resolution_source: string;
  resolution_chain: string[];
  cross_references: CrossReferences;
  confidence: string;
}

// Cross-Database Comparison
export interface CompareRequest {
  smiles?: string;
  inchikey?: string;
}

export interface DatabaseEntry {
  database: string;
  found: boolean;
  canonical_smiles?: string;
  kekulized_smiles?: string;
  inchi?: string;
  inchikey?: string;
  molecular_formula?: string;
  molecular_weight?: number;
  name?: string;
  url?: string;
}

export interface PropertyComparison {
  property_name: string;
  values: Record<string, string | null>;
  status: string;
  detail?: string;
}

export interface ConsistencyResult {
  entries: DatabaseEntry[];
  comparisons: PropertyComparison[];
  overall_verdict: string;
  summary: string;
}
