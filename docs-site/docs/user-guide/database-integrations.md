---
sidebar_position: 6
title: Database Integrations
description: Cross-reference molecules against PubChem, ChEMBL, COCONUT, and Wikidata databases
---

# Database Integrations

ChemAudit integrates with major chemical databases to enrich your molecular data with cross-references, bioactivity information, natural product annotations, and structural consistency checks.

## Supported Databases

### PubChem

| Data Available | Rate Limit |
|----------------|------------|
| Properties, synonyms, IUPAC names, CID | 30 req/min |

PubChem is the world's largest collection of freely accessible chemical information. Use it to:

- Find common names and synonyms
- Get IUPAC names
- Link to PubChem compound pages
- Retrieve molecular properties

### ChEMBL

| Data Available | Rate Limit |
|----------------|------------|
| Bioactivity, targets, clinical phase, ChEMBL ID | 30 req/min |

ChEMBL is a database of bioactive drug-like small molecules. Use it to:

- Find bioactivity data (IC50, Ki, etc.)
- Identify protein targets
- Check clinical development status
- Access published assay results

### COCONUT

| Data Available | Rate Limit |
|----------------|------------|
| Natural product data, organism source, COCONUT ID | 30 req/min |

COCONUT (COlleCtion of Open Natural ProdUcTs) contains natural product structures. Use it to:

- Identify natural product origins
- Find organism sources
- Check natural vs. synthetic classification
- Link to COCONUT database entries

### Wikidata

| Data Available | Rate Limit |
|----------------|------------|
| SMILES (isomeric), InChI, InChIKey, CAS, formula, mass | 30 req/min |

Wikidata is a free, open knowledge base linked to Wikipedia. Use it to:

- Get an independent structural representation from community-curated data
- Cross-check identifiers against a non-cheminformatics source
- Access CAS numbers and molecular mass
- Link to Wikidata entity pages

:::info Isomeric SMILES
ChemAudit fetches **isomeric SMILES** (Wikidata property P2017) when available, which preserves stereochemistry (R/S, E/Z). If isomeric SMILES is not available, it falls back to canonical SMILES (P233).
:::

## How to Search

### Web Interface

1. Enter your molecule on the Single Validation page
2. Navigate to the **Database Lookup** tab
3. View cross-references from all four databases:
   - PubChem properties and synonyms
   - ChEMBL bioactivity and targets
   - COCONUT natural product information
   - Wikidata identifiers and mass
4. The **Cross-Database Comparison** panel runs automatically (controlled by the **Auto-Compare** toggle), comparing structural identifiers across all databases that returned results.

Results load automatically after validation completes.

## Cross-Database Comparison

After database lookup, ChemAudit compares how each database represents your molecule. This helps catch registration errors, stereochemistry inconsistencies, or identifier mismatches across public databases.

### What Is Compared

The comparison checks three structural identifiers:

| Identifier | Comparison Method |
|------------|-------------------|
| **SMILES** | RDKit canonicalization, then tautomer-invariant comparison |
| **InChIKey** | Layer-by-layer analysis (connectivity, stereo, version) |
| **InChI** | Layer analysis (formula, connectivity, hydrogen, charge, stereo) |

### Verdict System

| Verdict | Meaning |
|---------|---------|
| **Consistent** | All compared identifiers match across databases |
| **Minor Differences** | Tautomeric or stereochemical variants of the same base structure |
| **Major Discrepancies** | Structural differences found -- different connectivity or genuinely different compounds |
| **No Data** | Fewer than 2 databases returned data, comparison not possible |

### The "Resolved" Entry

The comparison panel includes a **Resolved** entry alongside the external databases. This is ChemAudit's own RDKit-computed canonical representation of your input SMILES, providing a local reference point. It shows the **kekulized SMILES** form (explicit single/double bonds, no aromaticity notation).

### Stereochemistry Handling

ChemAudit always prefers **isomeric SMILES** (preserving `@`/`@@` stereocenters and `E`/`Z` double bonds) from each database:

- **PubChem**: Fetches `IsomericSMILES` over `CanonicalSMILES`/`ConnectivitySMILES`
- **ChEMBL**: `canonical_smiles` field includes stereochemistry by default
- **COCONUT**: `canonical_smiles` field includes stereochemistry by default
- **Wikidata**: Fetches P2017 (isomeric SMILES) over P233 (canonical SMILES)

This ensures that stereochemical information is preserved in all comparisons.

### API - Cross-Database Comparison

```bash
curl -X POST http://localhost:8001/api/v1/integrations/compare \
  -H "Content-Type: application/json" \
  -d '{"smiles": "C[C@H](N)C(=O)O"}'
```

Response:
```json
{
  "entries": [
    {
      "database": "PubChem",
      "found": true,
      "canonical_smiles": "C[C@@H](C(=O)O)N",
      "inchikey": "QNAYBMKLOCPYGJ-REOHCLBHSA-N",
      "inchi": "InChI=1S/C3H7NO2/..."
    },
    { "database": "ChEMBL", "found": true, "..." : "..." },
    { "database": "COCONUT", "found": true, "..." : "..." },
    { "database": "Wikidata", "found": true, "..." : "..." },
    { "database": "Resolved", "found": true, "..." : "..." }
  ],
  "comparisons": [
    {
      "property_name": "canonical_smiles",
      "status": "match",
      "detail": "Raw SMILES differ but resolve to same canonical form after RDKit canonicalization."
    }
  ],
  "overall_verdict": "consistent",
  "summary": "All compared properties match across databases"
}
```

## Identifier Resolution

ChemAudit can resolve a wide range of chemical identifiers to a canonical SMILES structure.

### Supported Identifier Types

| Type | Example |
|------|---------|
| SMILES | `CCO`, `C[C@H](N)C(=O)O` |
| InChI | `InChI=1S/C2H6O/c1-2-3/h3H,2H2,1H3` |
| InChIKey | `LFQSCWFLJHTTHZ-UHFFFAOYSA-N` |
| PubChem CID | `702` |
| ChEMBL ID | `CHEMBL545` |
| CAS Number | `64-17-5` |
| DrugBank ID | `DB00898` |
| ChEBI ID | `CHEBI:16236` |
| UNII | `3K9958V90M` |
| Wikipedia URL | `https://en.wikipedia.org/wiki/Ethanol` |
| Compound name | `aspirin`, `caffeine` |

### API - Resolve Identifier

```bash
curl -X POST http://localhost:8001/api/v1/integrations/resolve \
  -H "Content-Type: application/json" \
  -d '{"identifier": "CHEMBL25"}'
```

Response:
```json
{
  "resolved": true,
  "identifier_type_detected": "chembl_id",
  "canonical_smiles": "CC(=O)Oc1ccccc1C(=O)O",
  "inchi": "InChI=1S/C9H8O4/...",
  "inchikey": "BSYNRYMUTXBXSQ-UHFFFAOYSA-N",
  "iupac_name": "aspirin",
  "resolution_source": "chembl",
  "resolution_chain": ["chembl_id → ChEMBL API → SMILES"],
  "cross_references": {
    "pubchem_cid": 2244,
    "chembl_id": "CHEMBL25"
  },
  "confidence": "high"
}
```

### API - PubChem Lookup

```bash
curl -X POST http://localhost:8001/api/v1/integrations/pubchem/lookup \
  -H "Content-Type: application/json" \
  -d '{"molecule": "CCO", "format": "smiles"}'
```

Response:
```json
{
  "found": true,
  "cid": 702,
  "iupac_name": "ethanol",
  "molecular_formula": "C2H6O",
  "molecular_weight": 46.07,
  "synonyms": ["ethanol", "ethyl alcohol", "alcohol"],
  "url": "https://pubchem.ncbi.nlm.nih.gov/compound/702"
}
```

### API - ChEMBL Bioactivity

```bash
curl -X POST http://localhost:8001/api/v1/integrations/chembl/bioactivity \
  -H "Content-Type: application/json" \
  -d '{"molecule": "CCO", "format": "smiles"}'
```

Response:
```json
{
  "found": true,
  "chembl_id": "CHEMBL545",
  "pref_name": "ETHANOL",
  "max_phase": 4,
  "bioactivity_count": 1250,
  "bioactivities": [
    {
      "target_chembl_id": "CHEMBL240",
      "target_name": "GABA receptor",
      "activity_type": "IC50",
      "activity_value": 100.0,
      "activity_unit": "nM"
    }
  ],
  "url": "https://www.ebi.ac.uk/chembl/compound_report_card/CHEMBL545"
}
```

### API - COCONUT Lookup

```bash
curl -X POST http://localhost:8001/api/v1/integrations/coconut/lookup \
  -H "Content-Type: application/json" \
  -d '{"molecule": "CCO", "format": "smiles"}'
```

Response:
```json
{
  "found": true,
  "coconut_id": "CNP0123456",
  "name": "Compound Name",
  "smiles": "CCO",
  "molecular_weight": 46.07,
  "organism": "Genus species",
  "url": "https://coconut.naturalproducts.net/compounds/CNP0123456"
}
```

## Understanding Results

### PubChem Results

**CID (Compound ID)**: Unique PubChem identifier

**IUPAC Name**: Systematic chemical name

**Synonyms**: Common names, trade names, registry numbers

**URL**: Direct link to PubChem compound page

:::tip Synonym Search
PubChem synonyms help identify whether a molecule is a known drug, reagent, or natural product.
:::

### ChEMBL Results

**ChEMBL ID**: Unique ChEMBL identifier

**Preferred Name**: ChEMBL's canonical name for the compound

**Max Phase**: Clinical development status:
- **0**: Preclinical
- **1**: Phase I clinical trials
- **2**: Phase II clinical trials
- **3**: Phase III clinical trials
- **4**: Approved drug

**Bioactivity Count**: Number of bioactivity records

**Bioactivities**: Sample bioactivity data (up to 10 results):
- Target ID and name
- Activity type (IC50, Ki, EC50, etc.)
- Activity value and unit

:::info Clinical Phase
Max phase 4 indicates an approved drug. Phase 0 means the compound has bioactivity data but is not in clinical development.
:::

### COCONUT Results

**COCONUT ID**: Unique identifier in the COCONUT database

**Name**: Natural product name

**Organism**: Source organism (plant, fungi, bacteria, etc.)

**URL**: Direct link to COCONUT database entry

## Rate Limits

All database integrations are rate-limited to prevent abuse of external services:

| Database | Anonymous Limit | Authenticated Limit |
|----------|----------------|---------------------|
| PubChem | 30 req/min | 30 req/min |
| ChEMBL | 30 req/min | 30 req/min |
| COCONUT | 30 req/min | 30 req/min |
| Wikidata | 30 req/min | 30 req/min |

:::warning Rate Limiting
If you exceed rate limits, you'll receive a 429 error. Wait before retrying, or use batch processing which handles rate limiting automatically.
:::

## Handling Missing Data

Not all molecules are in all databases. The API returns `found: false` when no match is found:

```json
{
  "found": false,
  "message": "No match found in PubChem"
}
```

This is normal for:
- Novel compounds
- Synthetic intermediates
- Proprietary structures
- Very rare natural products

## Use Cases

### Compound Identification

Lookup molecules to confirm identity:

1. Validate structure
2. Cross-reference against databases
3. Compare synonyms and IDs with your records
4. Verify molecular formula and weight

### Bioactivity Screening

Check if molecules have known bioactivity:

1. Search ChEMBL
2. Review targets and activity types
3. Check clinical development status
4. Access published data for context

### Natural Product Analysis

Identify natural product origins:

1. Search COCONUT
2. Check organism source
3. Review natural product classification
4. Link to literature references

## Integration in Workflows

### Batch Processing

Database lookups can slow batch processing. Use them selectively:

```bash
# Process batch without database lookups (faster)
curl -X POST http://localhost:8001/api/v1/batch/upload \
  -F "file=@molecules.sdf"

# Then lookup specific molecules of interest individually
```

### Compound Curation

Enrich your database with cross-references:

1. Validate and standardize structures
2. Lookup each database
3. Store IDs (CID, ChEMBL ID, COCONUT ID)
4. Link to external pages for detailed information

## Best Practices

1. **Cache results**: Database APIs are rate-limited, cache lookups to avoid repeated requests
2. **Handle failures gracefully**: External APIs can be temporarily unavailable
3. **Verify matches**: Always verify the returned structure matches your input
4. **Respect rate limits**: Implement exponential backoff for retries
5. **Use InChIKey for lookups**: More reliable than SMILES for exact matching

## Next Steps

- **[IUPAC Name Conversion](/docs/user-guide/iupac-conversion)** — Enter chemical names directly instead of SMILES
- **[Exporting Results](/docs/user-guide/exporting-results)** — Export data with database cross-references
- **[Batch Processing](/docs/user-guide/batch-processing)** — Process large datasets
- **[API Reference](/docs/api/endpoints)** — Full integrations API documentation
