---
sidebar_position: 15
title: Identifier Resolution
description: Resolve any chemical identifier to canonical SMILES with cross-database linking
---

# Identifier Resolution

ChemAudit's universal identifier resolver converts any chemical identifier into a canonical SMILES structure with cross-database linking. It supports 10+ identifier types with automatic detection.

## Supported Identifier Types

| Type | Example | Detection Pattern |
|------|---------|-------------------|
| **SMILES** | `CCO`, `C[C@H](N)C(=O)O` | Valid SMILES characters |
| **InChI** | `InChI=1S/C2H6O/c1-2-3/h3H,2H2,1H3` | Starts with `InChI=` |
| **InChIKey** | `LFQSCWFLJHTTHZ-UHFFFAOYSA-N` | 27-character pattern with hyphens |
| **PubChem CID** | `702` | Numeric, resolved via PubChem |
| **ChEMBL ID** | `CHEMBL545` | Starts with `CHEMBL` |
| **CAS Number** | `64-17-5` | Digit-digit-digit pattern |
| **DrugBank ID** | `DB00898` | Starts with `DB` |
| **ChEBI ID** | `CHEBI:16236` | Starts with `CHEBI:` |
| **UNII** | `3K9958V90M` | 10-character alphanumeric |
| **Wikipedia URL** | `https://en.wikipedia.org/wiki/Ethanol` | Wikipedia URL pattern |
| **Compound name** | `aspirin`, `caffeine` | Fallback — resolved via PubChem/OPSIN |

## How It Works

1. **Auto-detection**: The identifier type is automatically detected from the input pattern
2. **Resolution chain**: The resolver follows a chain of lookups to convert the identifier to SMILES
3. **Cross-references**: After resolution, UniChem and database APIs are queried to populate cross-references (PubChem CID, ChEMBL ID, DrugBank ID, ChEBI ID, KEGG ID)

:::info Resolution Sources
Different identifier types use different resolution sources:
- **SMILES/InChI/InChIKey**: Direct RDKit conversion or UniChem lookup
- **PubChem CID**: PubChem PUG REST API
- **ChEMBL ID**: ChEMBL API
- **CAS/DrugBank/ChEBI/UNII**: UniChem cross-reference service
- **Wikipedia URL**: Wikidata SPARQL query
- **Compound names**: PubChem name resolution, then OPSIN as fallback
:::

## Using the Web Interface

The Identifier Resolver card appears on the **Database Lookup** tab of Single Validation:

1. Enter any supported identifier in the resolver input field
2. The identifier type is auto-detected and shown
3. Click **Resolve** (or it resolves automatically)
4. View results:
   - Resolved canonical SMILES
   - InChI and InChIKey
   - Molecular formula and weight
   - IUPAC name
   - Resolution source and chain
   - Cross-references with links to external databases

The resolved SMILES can be used directly for validation by clicking **Validate This Structure**.

## API Reference

### Resolve Identifier

```bash
curl -X POST http://localhost:8001/api/v1/integrations/resolve \
  -H "Content-Type: application/json" \
  -d '{
    "identifier": "CHEMBL25",
    "identifier_type": "auto"
  }'
```

| Field | Type | Required | Description |
|-------|------|----------|-------------|
| `identifier` | string | Yes | The identifier to resolve |
| `identifier_type` | string | No | `auto` (default), or specify: `smiles`, `inchi`, `inchikey`, `pubchem_cid`, `chembl_id`, `cas`, `drugbank_id`, `chebi_id`, `unii`, `wikipedia`, `name` |

**Response:**

```json
{
  "resolved": true,
  "identifier_type_detected": "chembl_id",
  "canonical_smiles": "CC(=O)Oc1ccccc1C(=O)O",
  "inchi": "InChI=1S/C9H8O4/c1-6(10)13-8-5-3-2-4-7(8)9(11)12/h2-5H,1H3,(H,11,12)",
  "inchikey": "BSYNRYMUTXBXSQ-UHFFFAOYSA-N",
  "molecular_formula": "C9H8O4",
  "molecular_weight": 180.16,
  "iupac_name": "aspirin",
  "resolution_source": "chembl",
  "resolution_chain": ["chembl_id → ChEMBL API → SMILES"],
  "cross_references": {
    "pubchem_cid": 2244,
    "chembl_id": "CHEMBL25",
    "drugbank_id": "DB00945",
    "chebi_id": "CHEBI:15365",
    "kegg_id": "D00109"
  }
}
```

If the identifier cannot be resolved:

```json
{
  "resolved": false,
  "identifier_type_detected": "unknown",
  "canonical_smiles": null,
  "resolution_source": null,
  "resolution_chain": ["Unable to detect identifier type"]
}
```

```python
import requests

# Resolve a ChEMBL ID
response = requests.post(
    "http://localhost:8001/api/v1/integrations/resolve",
    json={"identifier": "CHEMBL25"}
)
result = response.json()
if result["resolved"]:
    print(f"SMILES: {result['canonical_smiles']}")
    print(f"Source: {result['resolution_source']}")
    print(f"Cross-refs: {result['cross_references']}")

# Resolve a CAS number
response = requests.post(
    "http://localhost:8001/api/v1/integrations/resolve",
    json={"identifier": "50-78-2"}
)

# Resolve a compound name
response = requests.post(
    "http://localhost:8001/api/v1/integrations/resolve",
    json={"identifier": "ibuprofen"}
)
```

## Rate Limits

| Endpoint | Limit |
|----------|-------|
| `POST /integrations/resolve` | 30 req/min |

## Next Steps

- **[Database Integrations](/docs/user-guide/database-integrations)** — Full database lookup and comparison
- **[Single Validation](/docs/user-guide/single-validation)** — Validate resolved structures
- **[IUPAC Conversion](/docs/user-guide/iupac-conversion)** — Name-to-SMILES conversion
