<div align="center">

<img src="assets/logo.png" alt="ChemVault" width="80" />

# API Reference

### Complete REST API Documentation

[![OpenAPI](https://img.shields.io/badge/OpenAPI-3.0-85EA2D?logo=swagger&logoColor=white)](http://localhost:8000/docs)

**Interactive Docs:** http://localhost:8000/docs

</div>

---

## 📋 Table of Contents

- [Overview](#-overview)
- [Authentication](#-authentication)
- [Endpoints](#-endpoints)
  - [Validation](#validation)
  - [Alerts](#alerts)
  - [Scoring](#scoring)
  - [Standardization](#standardization)
  - [Batch Processing](#batch-processing)
  - [Database Integrations](#database-integrations)
- [Error Handling](#-error-handling)
- [Rate Limits](#-rate-limits)

---

## 🌐 Overview

### Base URL

```
http://localhost:8000/api/v1
```

### Content Type

All requests must include:
```
Content-Type: application/json
```

### Response Format

All responses follow this structure:

```json
{
  "data": { },      // Response payload
  "error": null,    // Error details (if any)
  "meta": {         // Metadata
    "request_id": "abc123",
    "processing_time_ms": 45
  }
}
```

---

## 🔐 Authentication

Authentication is **optional** for development. For production, use API keys:

```bash
curl -H "X-API-Key: your-api-key" \
  http://localhost:8000/api/v1/validate
```

---

## 📡 Endpoints

### Validation

#### `POST /validate`

Validate a single molecule.

**Request:**
```json
{
  "molecule": "CCO",
  "format": "auto",
  "checks": ["valence", "aromaticity", "stereo"]
}
```

| Field | Type | Required | Description |
|-------|------|----------|-------------|
| `molecule` | string | ✅ | SMILES, InChI, or MOL block |
| `format` | string | ❌ | `auto`, `smiles`, `inchi`, `mol` (default: `auto`) |
| `checks` | array | ❌ | Specific checks to run (default: all) |

**Response:**
```json
{
  "valid": true,
  "validation_score": 95,
  "checks": [
    {
      "name": "valence",
      "passed": true,
      "severity": "critical",
      "message": "All atoms have valid valence",
      "details": {}
    }
  ],
  "molecule_info": {
    "formula": "C2H6O",
    "molecular_weight": 46.07,
    "heavy_atom_count": 3
  },
  "standardized_smiles": "CCO",
  "inchi": "InChI=1S/C2H6O/c1-2-3/h3H,2H2,1H3",
  "inchi_key": "LFQSCWFLJHTTHZ-UHFFFAOYSA-N"
}
```

---

#### `GET /checks`

List available validation checks.

**Response:**
```json
{
  "checks": [
    {
      "name": "valence",
      "description": "Validates atom valence",
      "severity": "critical",
      "enabled_by_default": true
    }
  ]
}
```

---

### Alerts

#### `POST /alerts`

Screen molecule for structural alerts.

**Request:**
```json
{
  "molecule": "c1ccc2c(c1)nc(n2)Sc3nnnn3C",
  "format": "smiles",
  "catalogs": ["PAINS", "BRENK"]
}
```

| Field | Type | Required | Description |
|-------|------|----------|-------------|
| `molecule` | string | ✅ | Molecule to screen |
| `format` | string | ❌ | Input format (default: `auto`) |
| `catalogs` | array | ❌ | Alert catalogs (default: all) |

**Response:**
```json
{
  "has_alerts": true,
  "total_alerts": 2,
  "alerts": [
    {
      "catalog": "PAINS",
      "rule_id": "pains_a_123",
      "name": "thiazole_amine_A",
      "severity": "high",
      "description": "May cause assay interference",
      "matched_smarts": "[#7]-c1nc2ccccc2s1",
      "matched_atoms": [4, 5, 6, 7, 8]
    }
  ],
  "screened_catalogs": ["PAINS", "BRENK"]
}
```

---

#### `GET /alerts/catalogs`

List available alert catalogs.

**Response:**
```json
{
  "catalogs": [
    {
      "name": "PAINS",
      "description": "Pan-Assay Interference Compounds",
      "version": "1.0",
      "rule_count": 480
    }
  ]
}
```

---

### Scoring

#### `POST /score`

Calculate molecular scores.

**Request:**
```json
{
  "molecule": "CCO",
  "format": "smiles",
  "include": ["ml_readiness", "np_likeness", "scaffold"]
}
```

| Field | Type | Required | Description |
|-------|------|----------|-------------|
| `molecule` | string | ✅ | Molecule to score |
| `include` | array | ❌ | Score types to include |

**Response:**
```json
{
  "ml_readiness": {
    "score": 85,
    "components": {
      "descriptor_calculability": 90,
      "fingerprint_generation": 95,
      "structural_complexity": 70,
      "standardization": 85
    },
    "recommendations": []
  },
  "np_likeness": {
    "score": -1.2,
    "interpretation": "Synthetic-like",
    "percentile": 25
  },
  "scaffold": {
    "murcko_scaffold": "c1ccccc1",
    "generic_scaffold": "*c1ccccc1",
    "scaffold_class": "Benzene"
  }
}
```

---

### Standardization

#### `POST /standardize`

Standardize a molecule.

**Request:**
```json
{
  "molecule": "CC(=O)Oc1ccccc1C(=O)[O-].[Na+]",
  "format": "smiles",
  "options": {
    "normalize": true,
    "reionize": true,
    "remove_salts": true,
    "neutralize": true,
    "remove_stereo": false,
    "canonicalize": true
  }
}
```

**Response:**
```json
{
  "original_smiles": "CC(=O)Oc1ccccc1C(=O)[O-].[Na+]",
  "standardized_smiles": "CC(=O)Oc1ccccc1C(=O)O",
  "standardized_inchi": "InChI=1S/C9H8O4/...",
  "standardized_inchi_key": "BSYNRYMUTXBXSQ-UHFFFAOYSA-N",
  "changes_applied": [
    {
      "step": "remove_salts",
      "description": "Removed fragment: [Na+]"
    },
    {
      "step": "neutralize",
      "description": "Neutralized carboxylate anion"
    }
  ]
}
```

---

### Batch Processing

#### `POST /batch/upload`

Upload file for batch processing.

**Request:** `multipart/form-data`

| Field | Type | Required | Description |
|-------|------|----------|-------------|
| `file` | file | ✅ | SDF or CSV file |
| `smiles_column` | string | ❌ | Column name for CSV (default: `SMILES`) |

**Response:**
```json
{
  "job_id": "550e8400-e29b-41d4-a716-446655440000",
  "status": "pending",
  "total_molecules": 10000,
  "message": "Job submitted. Processing 10000 molecules."
}
```

---

#### `GET /batch/{job_id}`

Get batch job results.

**Query Parameters:**

| Param | Type | Description |
|-------|------|-------------|
| `page` | int | Page number (default: 1) |
| `page_size` | int | Results per page (default: 50, max: 100) |
| `status_filter` | string | Filter by `success` or `error` |
| `min_score` | int | Minimum validation score |
| `max_score` | int | Maximum validation score |

**Response:**
```json
{
  "job_id": "550e8400-e29b-41d4-a716-446655440000",
  "status": "complete",
  "statistics": {
    "total": 10000,
    "successful": 9850,
    "errors": 150,
    "avg_validation_score": 87.5,
    "processing_time_seconds": 125.4
  },
  "results": [
    {
      "index": 0,
      "smiles": "CCO",
      "name": "Ethanol",
      "status": "success",
      "validation": { ... },
      "alerts": { ... },
      "scoring": { ... }
    }
  ],
  "page": 1,
  "page_size": 50,
  "total_pages": 200
}
```

---

#### `GET /batch/{job_id}/status`

Get job status (lightweight).

**Response:**
```json
{
  "job_id": "550e8400-e29b-41d4-a716-446655440000",
  "status": "processing",
  "progress": 45.5,
  "processed": 4550,
  "total": 10000,
  "eta_seconds": 68
}
```

---

#### `DELETE /batch/{job_id}`

Cancel a running batch job.

**Response:**
```json
{
  "job_id": "550e8400-e29b-41d4-a716-446655440000",
  "status": "cancelling",
  "message": "Job cancellation requested"
}
```

---

### Database Integrations

#### `POST /integrations/pubchem/lookup`

Search PubChem database.

**Request:**
```json
{
  "molecule": "CCO",
  "format": "smiles"
}
```

**Response:**
```json
{
  "found": true,
  "cid": 702,
  "name": "Ethanol",
  "synonyms": ["Ethyl alcohol", "Alcohol"],
  "molecular_formula": "C2H6O",
  "molecular_weight": 46.07,
  "iupac_name": "ethanol",
  "canonical_smiles": "CCO",
  "properties": {
    "xlogp": -0.3,
    "hbd": 1,
    "hba": 1,
    "tpsa": 20.2,
    "rotatable_bonds": 0
  }
}
```

---

#### `POST /integrations/chembl/bioactivity`

Search ChEMBL for bioactivity data.

**Response:**
```json
{
  "found": true,
  "chembl_id": "CHEMBL545",
  "pref_name": "ETHANOL",
  "molecule_type": "Small molecule",
  "max_phase": 4,
  "bioactivities_count": 1250,
  "targets": [
    {
      "target_chembl_id": "CHEMBL240",
      "target_name": "GABA receptor",
      "activity_count": 45
    }
  ]
}
```

---

#### `POST /integrations/coconut/lookup`

Search COCONUT natural products database.

**Response:**
```json
{
  "found": true,
  "coconut_id": "CNP0123456",
  "name": "Compound Name",
  "np_likeness_score": 2.1,
  "source_organisms": ["Genus species"],
  "molecular_framework": "Terpenoid"
}
```

---

## ⚠️ Error Handling

### Error Response Format

```json
{
  "error": {
    "code": "VALIDATION_ERROR",
    "message": "Invalid SMILES string",
    "details": {
      "molecule": "Invalid character at position 5"
    }
  }
}
```

### Error Codes

| Code | HTTP Status | Description |
|------|-------------|-------------|
| `VALIDATION_ERROR` | 400 | Invalid input data |
| `PARSE_ERROR` | 400 | Cannot parse molecule |
| `NOT_FOUND` | 404 | Resource not found |
| `RATE_LIMITED` | 429 | Too many requests |
| `INTERNAL_ERROR` | 500 | Server error |

---

## 🚦 Rate Limits

| Endpoint | Limit |
|----------|-------|
| `/validate` | 100 req/min |
| `/alerts` | 100 req/min |
| `/score` | 100 req/min |
| `/batch/upload` | 10 req/min |
| `/integrations/*` | 30 req/min |

**Rate Limit Headers:**
```
X-RateLimit-Limit: 100
X-RateLimit-Remaining: 95
X-RateLimit-Reset: 1706500800
```

---

<div align="center">

**Interactive API Explorer:** http://localhost:8000/docs

</div>
