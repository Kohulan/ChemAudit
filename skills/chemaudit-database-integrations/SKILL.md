---
name: chemaudit-database-integrations
description: >
  Look up molecules across PubChem, ChEMBL, COCONUT, Wikidata, and SureChEMBL
  patent databases, resolve chemical identifiers (CAS, ChEMBL ID, PubChem CID,
  ChEBI, UNII, DrugBank, names), and run cross-database consistency checks. Use
  when user says "look up in PubChem", "ChEMBL bioactivity", "natural product
  search", "COCONUT lookup", "Wikidata lookup", "resolve CAS number", "resolve
  ChEMBL ID", "cross-database check", "find this compound", "patent search",
  "SureChEMBL", or "identifier resolution". Supports SMILES, InChI, InChIKey,
  CAS, ChEMBL ID, PubChem CID, ChEBI, UNII, DrugBank ID, Wikipedia URL, and
  compound names (OPSIN + PubChem fallback).
license: MIT
allowed-tools: "Bash(curl:*) Bash(python:*)"
compatibility: >
  Requires a running ChemAudit instance with outbound network access to PubChem,
  ChEMBL, COCONUT, Wikidata, and SureChEMBL APIs.
metadata:
  author: Kohulan Rajan
  version: 1.0.0
  mcp-server: chemaudit
  category: cheminformatics
  tags: [PubChem, ChEMBL, COCONUT, Wikidata, identifier-resolution, databases, SureChEMBL]
---

# ChemAudit Database Integrations

## Overview

Five external-database lookups plus two cross-database tools:

| Endpoint | Purpose |
|---|---|
| `POST /api/v1/integrations/pubchem/lookup` | PubChem compound info (CID, IUPAC, synonyms, formula, MW) |
| `POST /api/v1/integrations/chembl/bioactivity` | ChEMBL molecule + bioactivity records |
| `POST /api/v1/integrations/coconut/lookup` | COCONUT natural-product database (>400K entries, organism, NP-likeness) |
| `POST /api/v1/integrations/wikidata/lookup` | Wikidata via SPARQL (label, CAS, Wikipedia link) |
| `POST /api/v1/integrations/surechembl/lookup` | SureChEMBL patent presence |
| `POST /api/v1/integrations/compare` | Consistency check across PubChem + ChEMBL + COCONUT |
| `POST /api/v1/resolve` | Universal identifier resolver for 11 identifier types |

All lookup endpoints accept `{smiles: "...", inchikey: "..."}` — provide at least one. Identifier resolution accepts anything.

## Workflow

### 1. PubChem

```bash
curl -sS -X POST http://localhost:8000/api/v1/integrations/pubchem/lookup \
  -H 'Content-Type: application/json' \
  -d '{"smiles": "CC(=O)Oc1ccccc1C(=O)O"}'
```

Response:

```json
{
  "found": true,
  "cid": 2244,
  "iupac_name": "2-acetyloxybenzoic acid",
  "molecular_formula": "C9H8O4",
  "molecular_weight": 180.16,
  "canonical_smiles": "CC(=O)OC1=CC=CC=C1C(=O)O",
  "inchi": "InChI=1S/C9H8O4/c1-6(10)13-8-5-3-2-4-7(8)9(11)12/h2-5H,1H3,(H,11,12)",
  "inchikey": "BSYNRYMUTXBXSQ-UHFFFAOYSA-N",
  "synonyms": ["aspirin", "acetylsalicylic acid", "2-acetoxybenzoic acid", ...],
  "url": "https://pubchem.ncbi.nlm.nih.gov/compound/2244"
}
```

### 2. ChEMBL bioactivity

```bash
curl -sS -X POST http://localhost:8000/api/v1/integrations/chembl/bioactivity \
  -H 'Content-Type: application/json' \
  -d '{"inchikey": "BSYNRYMUTXBXSQ-UHFFFAOYSA-N"}'
```

Response:

```json
{
  "found": true,
  "chembl_id": "CHEMBL25",
  "pref_name": "ASPIRIN",
  "molecule_type": "Small molecule",
  "max_phase": 4,
  "molecular_formula": "C9H8O4",
  "molecular_weight": 180.16,
  "canonical_smiles": "...",
  "inchikey": "BSYNRYMUTXBXSQ-UHFFFAOYSA-N",
  "bioactivities": [
    {
      "target_chembl_id": "CHEMBL240",
      "target_name": "Cyclooxygenase-1",
      "target_type": "SINGLE PROTEIN",
      "activity_type": "IC50",
      "activity_value": 100.0,
      "activity_unit": "nM",
      "assay_chembl_id": "CHEMBL615139",
      "document_chembl_id": "CHEMBL1127040"
    }
  ],
  "bioactivity_count": 1234,
  "url": "https://www.ebi.ac.uk/chembl/compound_report_card/CHEMBL25"
}
```

`max_phase`: 0 (preclinical), 1-3 (clinical trials), 4 (approved).

### 3. COCONUT natural products

```bash
curl -sS -X POST http://localhost:8000/api/v1/integrations/coconut/lookup \
  -H 'Content-Type: application/json' \
  -d '{"smiles": "CN1C=NC2=C1C(=O)N(C(=O)N2C)C"}'
```

Response:

```json
{
  "found": true,
  "coconut_id": "CNP0123456",
  "name": "Caffeine",
  "smiles": "CN1C=NC2=C1C(=O)N(C(=O)N2C)C",
  "inchi": "...",
  "inchikey": "RYYVLZVUVIJVGH-UHFFFAOYSA-N",
  "molecular_formula": "C8H10N4O2",
  "molecular_weight": 194.19,
  "organism": "Coffea arabica",
  "organism_type": "Plantae",
  "nplikeness": 2.4,
  "url": "https://coconut.naturalproducts.net/compounds/CNP0123456"
}
```

### 4. Wikidata

```bash
curl -sS -X POST http://localhost:8000/api/v1/integrations/wikidata/lookup \
  -H 'Content-Type: application/json' \
  -d '{"inchikey": "BSYNRYMUTXBXSQ-UHFFFAOYSA-N"}'
```

Returns label, canonical SMILES, InChIKey, CAS, formula, MW, and a Wikidata entity URL. SPARQL-based — typically slower than the other integrations. For a Wikipedia URL, use `/resolve` and read `cross_references.wikipedia_url`.

### 5. SureChEMBL patent lookup

```bash
curl -sS -X POST http://localhost:8000/api/v1/integrations/surechembl/lookup \
  -H 'Content-Type: application/json' \
  -d '{"inchikey": "BSYNRYMUTXBXSQ-UHFFFAOYSA-N"}'
```

Response:

```json
{
  "found": true,
  "schembl_id": "SCHEMBL25",
  "url": "https://www.surechembl.org/chemical/SCHEMBL25",
  "patent_count": 8432,
  "source": "unichem+surechembl",
  "smiles": "...",
  "inchikey": "BSYNRYMUTXBXSQ-UHFFFAOYSA-N",
  "molecular_weight": 180.16
}
```

Binary present/absent with direct link; patent count populated when SureChEMBL's enrichment API is reachable.

### 6. Cross-database comparison

```bash
curl -sS -X POST http://localhost:8000/api/v1/integrations/compare \
  -H 'Content-Type: application/json' \
  -d '{"inchikey": "BSYNRYMUTXBXSQ-UHFFFAOYSA-N"}'
```

At least one of `smiles` or `inchikey` required. Runs PubChem, ChEMBL, COCONUT in parallel and aligns outputs.

Response:

```json
{
  "entries": [
    {"database": "PubChem", "found": true, "canonical_smiles": "...", "inchikey": "...", "molecular_formula": "C9H8O4", "molecular_weight": 180.16, "name": "...", "url": "..."},
    {"database": "ChEMBL", "found": true, ...},
    {"database": "COCONUT", "found": false, ...}
  ],
  "comparisons": [
    {"property_name": "canonical_smiles", "values": {"PubChem": "...", "ChEMBL": "..."}, "status": "agree", "detail": null},
    {"property_name": "molecular_formula", "values": {"PubChem": "C9H8O4", "ChEMBL": "C9H8O4"}, "status": "agree", "detail": null},
    {"property_name": "molecular_weight", "values": {"PubChem": "180.16", "ChEMBL": "180.159"}, "status": "minor_mismatch", "detail": "Difference within rounding tolerance"}
  ],
  "overall_verdict": "consistent_minor",
  "summary": "All databases agree on structural identity; minor numeric differences within tolerance."
}
```

Verdicts: `consistent`, `consistent_minor` (minor numeric differences), `inconsistent` (structural disagreement), `no_data` (no database found the compound).

### 7. Universal identifier resolver

Accepts **11 identifier types** and resolves to canonical SMILES with cross-references. Use this before any other lookup when you don't know the structural representation.

```bash
# CAS number
curl -sS -X POST http://localhost:8000/api/v1/resolve \
  -H 'Content-Type: application/json' \
  -d '{"identifier": "50-78-2", "identifier_type": "auto"}'

# ChEMBL ID
curl -sS -X POST http://localhost:8000/api/v1/resolve \
  -H 'Content-Type: application/json' \
  -d '{"identifier": "CHEMBL25", "identifier_type": "auto"}'

# PubChem CID
curl -sS -X POST http://localhost:8000/api/v1/resolve \
  -H 'Content-Type: application/json' \
  -d '{"identifier": "2244"}'

# Compound name (OPSIN then PubChem)
curl -sS -X POST http://localhost:8000/api/v1/resolve \
  -H 'Content-Type: application/json' \
  -d '{"identifier": "acetylsalicylic acid"}'

# Explicit type hint
curl -sS -X POST http://localhost:8000/api/v1/resolve \
  -H 'Content-Type: application/json' \
  -d '{"identifier": "DB00945", "identifier_type": "drugbank_id"}'
```

Supported types: `auto`, `smiles`, `inchi`, `inchikey`, `pubchem_cid`, `chembl_id`, `cas`, `drugbank_id`, `chebi_id`, `unii`, `wikipedia`, `name`.

Response:

```json
{
  "resolved": true,
  "identifier_type_detected": "cas",
  "canonical_smiles": "CC(=O)OC1=CC=CC=C1C(=O)O",
  "inchi": "...",
  "inchikey": "BSYNRYMUTXBXSQ-UHFFFAOYSA-N",
  "molecular_formula": "C9H8O4",
  "molecular_weight": 180.16,
  "iupac_name": "2-acetyloxybenzoic acid",
  "resolution_source": "pubchem",
  "resolution_chain": ["cas", "pubchem_cid", "smiles"],
  "cross_references": {
    "pubchem_cid": 2244,
    "chembl_id": "CHEMBL25",
    "coconut_id": null,
    "drugbank_id": "DB00945",
    "chebi_id": "CHEBI:15365",
    "unii": "R16CO5Y76E",
    "cas": "50-78-2",
    "wikipedia_url": "https://en.wikipedia.org/wiki/Aspirin",
    "kegg_id": "C01405"
  },
  "confidence": "high"
}
```

`resolution_chain` shows the UniChem/PubChem hops taken. `confidence` values: `high`, `medium`, `low`, `none`.

## Examples

### Example 1 — "Look up aspirin in all databases"

1. Resolve the name first: `POST /resolve {"identifier": "aspirin"}` → get InChIKey.
2. `POST /integrations/compare` with the InChIKey → one consolidated response covering PubChem, ChEMBL, COCONUT.
3. For patent info: `POST /integrations/surechembl/lookup` separately.

### Example 2 — "Convert this CAS number to SMILES"

Just one call: `POST /resolve {"identifier": "50-78-2"}` → returns canonical SMILES plus UniChem cross-references to all other databases.

### Example 3 — "Has anyone filed a patent on this compound?"

`POST /integrations/surechembl/lookup {"inchikey": "..."}` → `found` tells you yes/no, `patent_count` tells you how many when enrichment API is available.

### Example 4 — "Get all ChEMBL bioactivity records for aspirin and filter to IC50 values in nM"

1. Resolve to ChEMBL ID or InChIKey.
2. `POST /integrations/chembl/bioactivity {"inchikey": "..."}`.
3. `bioactivities[]` has per-record `activity_type`, `activity_value`, `activity_unit`. Filter client-side to `activity_type=="IC50"` and `activity_unit=="nM"`.

### Example 5 — "Are the PubChem and ChEMBL SMILES consistent for this compound?"

`POST /integrations/compare` with the InChIKey → `comparisons[]` has per-property agreement status. `overall_verdict = "inconsistent"` means structures disagree — worth investigating.

## Rate limits

- `/integrations/pubchem/lookup`, `/chembl/bioactivity`, `/coconut/lookup`, `/wikidata/lookup`, `/surechembl/lookup`: 30/min.
- `/integrations/compare`: 10/min.
- `/resolve`: 30/min.

Upstream APIs (PubChem, ChEMBL, etc.) have their own rate limits. Excessive queries may trigger temporary blocks from the upstream source — the ChemAudit service respects `Retry-After` headers where provided.

## Troubleshooting

### `found: false` but the compound definitely exists in PubChem

Try the other identifier. Some compounds are indexed by InChIKey but not easily queried by SMILES (and vice-versa) due to tautomer/isomer normalization on the upstream side. If `/resolve` finds it via `name`, use the returned InChIKey.

### `/resolve`: `resolved: false`, `confidence: "none"`

None of the 11 resolvers matched. Common causes:
- Typo in the identifier.
- Non-standard format (e.g. CAS with no dashes).
- Compound is too new and hasn't been indexed in UniChem yet.

### Slow Wikidata responses

SPARQL endpoint latency is highly variable. Wikidata is typically the slowest integration; consider skipping it for bulk workflows.

### `/integrations/compare`: `verdict: "inconsistent"`

Structural disagreement across databases — usually a charge state, tautomer, or salt-form difference. Read `comparisons[]` to isolate which property differs. Often resolved by standardizing first (`chemaudit-standardization`) and re-querying.

### ChEMBL returns 0 bioactivities for a known drug

Rare but possible when the compound is registered under a different ChEMBL ID (e.g. a salt form vs parent). `/resolve` the name, get the `chembl_id` from `cross_references`, then query directly.

### SureChEMBL `found: false` but you know the compound is patented

UniChem's SureChEMBL cross-reference isn't comprehensive. When found-by-UniChem misses, the SureChEMBL direct API may still succeed — the service attempts both in sequence.

## Further reading

- `chemaudit-molecule-validation` — validate the structure once resolved.
- `chemaudit-standardization` — before cross-database comparison, standardize to reduce salt-form mismatches.
- `chemaudit-dataset-intelligence` — for dataset-level cross-database consistency checks.
