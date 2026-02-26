---
sidebar_position: 3
title: Endpoints
description: Complete ChemAudit API endpoint reference with request/response examples
---

# API Endpoints

Complete reference for all ChemAudit API endpoints, organized by feature.

## Health and Configuration

### GET /health

Health check endpoint. Returns system status and RDKit version.

**Response:**

```json
{
  "status": "healthy",
  "app_name": "ChemAudit",
  "app_version": "1.0.0",
  "rdkit_version": "2025.09.3"
}
```

### GET /config

Get public application configuration including deployment limits.

**Response:**

```json
{
  "app_name": "ChemAudit",
  "app_version": "1.0.0",
  "deployment_profile": "medium",
  "limits": {
    "max_batch_size": 10000,
    "max_file_size_mb": 500,
    "max_file_size_bytes": 524288000
  }
}
```

## Validation

### POST /validate

Validate a single molecule. Results are cached by InChIKey for 1 hour.

**Request:**

```json
{
  "molecule": "CCO",
  "format": "auto",
  "checks": ["all"],
  "preserve_aromatic": false
}
```

| Field | Type | Required | Description |
|-------|------|----------|-------------|
| `molecule` | string | Yes | SMILES, InChI, or MOL block (max 10,000 chars) |
| `format` | string | No | `auto`, `smiles`, `inchi`, `mol` (default: `auto`) |
| `checks` | array | No | Specific checks to run (default: `["all"]`) |
| `preserve_aromatic` | bool | No | Keep aromatic notation in output SMILES (default: `false`) |

**Response:**

```json
{
  "status": "completed",
  "molecule_info": {
    "input_smiles": "CCO",
    "canonical_smiles": "CCO",
    "inchi": "InChI=1S/C2H6O/c1-2-3/h3H,2H2,1H3",
    "inchikey": "LFQSCWFLJHTTHZ-UHFFFAOYSA-N",
    "molecular_formula": "C2H6O",
    "molecular_weight": 46.07,
    "num_atoms": 3
  },
  "overall_score": 95,
  "issues": [],
  "all_checks": [
    {
      "check_name": "valence_check",
      "passed": true,
      "severity": "critical",
      "message": "All atoms have valid valence"
    }
  ],
  "execution_time_ms": 12,
  "cached": false
}
```

### POST /validate/async

Validate using the Celery high-priority queue. Use when batch jobs are running.

**Query Parameters:**

| Param | Type | Description |
|-------|------|-------------|
| `timeout` | int | Max wait time in seconds (1-60, default: 30) |

**Request:** Same as POST /validate

**Response:** Same as POST /validate with additional `queue` field

### GET /checks

List available validation checks grouped by category.

**Response:**

```json
{
  "basic": ["valence_check", "kekulize_check", "sanitization_check"],
  "stereo": ["undefined_stereo_check", "stereo_consistency_check"],
  "representation": ["smiles_length_check", "inchi_generation_check"]
}
```

## Structural Alerts

### POST /alerts

Screen molecule for structural alerts.

**Request:**

```json
{
  "molecule": "c1ccc2c(c1)nc(n2)Sc3nnnn3C",
  "format": "smiles",
  "catalogs": ["PAINS"]
}
```

**Response:**

```json
{
  "status": "completed",
  "alerts": [
    {
      "pattern_name": "thiazole_amine_A",
      "description": "Potential assay interference",
      "severity": "warning",
      "matched_atoms": [4, 5, 6, 7],
      "catalog_source": "PAINS_A",
      "smarts": "[#7]-c1nc2ccccc2s1"
    }
  ],
  "total_alerts": 1,
  "has_critical": false,
  "has_warning": true
}
```

### POST /alerts/quick-check

Fast check for presence of alerts (no detailed results).

**Request:** Same as POST /alerts

**Response:**

```json
{
  "has_alerts": true,
  "checked_catalogs": ["PAINS"]
}
```

### GET /alerts/catalogs

List available structural alert catalogs.

**Response:**

```json
{
  "catalogs": {
    "PAINS": {
      "name": "PAINS",
      "description": "Pan-Assay Interference Compounds",
      "pattern_count": 480,
      "severity": "warning"
    }
  },
  "default_catalogs": ["PAINS"]
}
```

## Scoring

### POST /score

Calculate comprehensive molecular scores.

**Request:**

```json
{
  "molecule": "CCO",
  "format": "smiles",
  "include": ["ml_readiness", "druglikeness", "admet"]
}
```

**Available score types:**

- `ml_readiness`: ML-readiness score (0-100)
- `np_likeness`: Natural product likeness
- `scaffold`: Murcko scaffold extraction
- `druglikeness`: Drug-likeness filters
- `safety_filters`: Safety alerts summary
- `admet`: ADMET predictions
- `aggregator`: Aggregator likelihood

**Response:** See [User Guide - Scoring](/docs/user-guide/scoring/overview) for detailed response schemas.

## Standardization

### POST /standardize

Standardize a molecule using ChEMBL-compatible pipeline.

**Request:**

```json
{
  "molecule": "CC(=O)Oc1ccccc1C(=O)[O-].[Na+]",
  "format": "smiles",
  "options": {
    "include_tautomer": false,
    "preserve_stereo": true
  }
}
```

**Response:**

```json
{
  "result": {
    "original_smiles": "CC(=O)Oc1ccccc1C(=O)[O-].[Na+]",
    "standardized_smiles": "CC(=O)OC1=CC=CC=C1C(O)=O",
    "success": true,
    "steps_applied": [...],
    "excluded_fragments": ["[Na+]"],
    "structure_comparison": {
      "mass_change_percent": -10.87
    }
  }
}
```

### GET /standardize/options

Get available standardization options with descriptions.

## Batch Processing

### POST /batch/upload

Upload file for batch processing.

**Request:** `multipart/form-data`

| Field | Type | Required | Description |
|-------|------|----------|-------------|
| `file` | file | Yes | SDF, CSV, TSV, or TXT file |
| `smiles_column` | string | No | Column name for SMILES in CSV |
| `name_column` | string | No | Column name for molecule names/IDs |
| `include_extended_safety` | bool | No | Include NIH and ZINC filters |
| `include_chembl_alerts` | bool | No | Include ChEMBL pharma filters |
| `include_standardization` | bool | No | Run standardization pipeline |

**Response:**

```json
{
  "job_id": "550e8400-e29b-41d4-a716-446655440000",
  "status": "pending",
  "total_molecules": 10000,
  "message": "Job submitted. Processing 10000 molecules."
}
```

### POST /batch/detect-columns

Detect columns in delimited text file for SMILES and Name selection.

### GET /batch/\{job_id\}

Get batch job results with pagination and filtering.

**Query Parameters:**

| Param | Type | Description |
|-------|------|-------------|
| `page` | int | Page number (default: 1) |
| `page_size` | int | Results per page (1-100, default: 50) |
| `status_filter` | string | Filter by `success` or `error` |
| `min_score` | int | Minimum validation score (0-100) |
| `max_score` | int | Maximum validation score (0-100) |
| `sort_by` | string | Sort field |
| `sort_dir` | string | Sort direction (`asc`, `desc`) |

### GET /batch/\{job_id\}/status

Lightweight status check for a batch job.

**Response:**

```json
{
  "job_id": "550e8400-...",
  "status": "processing",
  "progress": 45.5,
  "processed": 455,
  "total": 1000,
  "eta_seconds": 68
}
```

### GET /batch/\{job_id\}/stats

Get aggregate statistics for a batch job (without individual results).

### DELETE /batch/\{job_id\}

Cancel a running batch job.

## Export

### GET /batch/\{job_id\}/export

Export batch results to a file.

**Query Parameters:**

| Param | Type | Required | Description |
|-------|------|----------|-------------|
| `format` | string | Yes | `csv`, `excel`, `sdf`, `json`, or `pdf` |
| `score_min` | int | No | Minimum validation score filter |
| `score_max` | int | No | Maximum validation score filter |
| `status` | string | No | Filter by `success`, `error`, or `warning` |
| `indices` | string | No | Comma-separated molecule indices |

**Response:** File download with appropriate `Content-Disposition` header

### POST /batch/\{job_id\}/export

Export with molecule selection via request body (for large selections).

**Request:**

```json
{
  "indices": [0, 1, 5, 23, 42]
}
```

Query parameters same as GET version.

## Database Integrations

### POST /integrations/pubchem/lookup

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
  "iupac_name": "ethanol",
  "synonyms": ["ethanol", "ethyl alcohol"],
  "url": "https://pubchem.ncbi.nlm.nih.gov/compound/702"
}
```

### POST /integrations/chembl/bioactivity

Search ChEMBL for bioactivity data.

### POST /integrations/coconut/lookup

Search COCONUT natural products database.

## API Key Management

All endpoints require admin authentication via `X-Admin-Secret` header.

### POST /api-keys

Create a new API key.

**Request:**

```json
{
  "name": "my-application",
  "description": "Key for production use",
  "expiry_days": 90
}
```

**Response (201):**

```json
{
  "key": "the-full-api-key-shown-only-once",
  "name": "my-application",
  "created_at": "2026-02-01T00:00:00Z",
  "expires_at": "2026-05-02T00:00:00Z"
}
```

### GET /api-keys

List all API keys (metadata only, not the actual key values).

### DELETE /api-keys/\{key_id\}

Revoke an API key. Returns 204 No Content on success.

## Bookmarks

### GET /bookmarks

List bookmarks with pagination and optional filters.

**Query Parameters:**

| Param | Type | Description |
|-------|------|-------------|
| `page` | int | Page number (default: 1) |
| `page_size` | int | Results per page (1-200, default: 50) |
| `tag` | string | Filter by tag |
| `search` | string | SMILES substring search |
| `source` | string | Filter by source |

### POST /bookmarks

Create a new bookmark. Returns 201.

**Request:**

```json
{
  "smiles": "CCO",
  "name": "Ethanol",
  "tags": ["simple", "alcohol"],
  "notes": "Test molecule",
  "source": "single_validation"
}
```

### PUT /bookmarks/\{id\}

Update bookmark tags and notes. Returns 404 if not found.

### DELETE /bookmarks/\{id\}

Delete a single bookmark. Returns 204 on success.

### DELETE /bookmarks/bulk

Bulk delete bookmarks by IDs. Query parameter: `ids` (list of integers). Returns 204.

### POST /bookmarks/batch-submit

Submit selected bookmarks as a new batch job.

**Request:**

```json
{
  "bookmark_ids": [1, 2, 3]
}
```

## History

### GET /history

Paginated validation audit trail with filters.

**Query Parameters:**

| Param | Type | Description |
|-------|------|-------------|
| `page` | int | Page number (default: 1) |
| `page_size` | int | Results per page (1-200, default: 50) |
| `date_from` | ISO8601 | Filter entries after this date |
| `date_to` | ISO8601 | Filter entries before this date |
| `outcome` | string | `pass`, `warn`, or `fail` |
| `source` | string | `single` or `batch` |
| `smiles_search` | string | SMILES substring match |

### GET /history/stats

Get summary statistics: total validations, outcome distribution, source distribution.

## Scoring Profiles

### GET /profiles

List all scoring profiles (8 presets + user-created).

### POST /profiles

Create a custom scoring profile. Returns 201.

### GET /profiles/\{id\}

Get a single profile by ID.

### PUT /profiles/\{id\}

Update a custom profile. Returns 400 if attempting to modify a preset.

### DELETE /profiles/\{id\}

Soft-delete a custom profile. Returns 400 for presets. Returns 204 on success.

### POST /profiles/\{id\}/duplicate

Duplicate any profile (including presets) as a new custom profile. Returns 201.

### GET /profiles/\{id\}/export

Export a profile as JSON for sharing.

### POST /profiles/import

Import a profile from JSON. Returns 201.

## Batch Analytics

### GET /batch/\{job_id\}/analytics

Get analytics status and cached results for a completed batch job.

### POST /batch/\{job_id\}/analytics/\{type\}

Trigger an on-demand analytics computation. Types: `scaffold`, `chemical_space`, `mmp`, `similarity_search`, `rgroup`.

**Optional body parameters:**

| Param | Type | Used By |
|-------|------|---------|
| `method` | string | chemical_space (`pca` or `tsne`) |
| `activity_column` | string | mmp (property key for cliff detection) |
| `query_smiles` | string | similarity_search |
| `query_index` | int | similarity_search |
| `top_k` | int | similarity_search (default: 10) |
| `core_smarts` | string | rgroup (SMARTS pattern) |

## Batch Subset Actions

### POST /batch/\{job_id\}/subset/revalidate

Re-validate a subset as a new batch job. Body: `{"indices": [0, 5, 12]}`.

### POST /batch/\{job_id\}/subset/rescore

Re-score a subset with a different profile. Body: `{"indices": [0, 5], "profile_id": 3}`.

### POST /batch/\{job_id\}/subset/score-inline

Synchronous inline scoring. Body: `{"indices": [0, 5], "profile_id": 4}`.

### POST /batch/\{job_id\}/subset/export

Export selected molecules only. Body: `{"indices": [0, 5], "format": "csv"}`.

## Similarity

### POST /validate/similarity

Calculate ECFP4 Tanimoto similarity between two molecules.

**Request:**

```json
{
  "smiles_a": "CC(=O)Oc1ccccc1C(=O)O",
  "smiles_b": "CC(=O)Nc1ccc(O)cc1"
}
```

## Permalinks

### POST /permalinks

Create a shareable permalink for a batch report. Returns short_id, URL, and 30-day expiry.

### GET /report/\{short_id\}

Resolve a permalink. Returns job data and snapshot. Returns 410 Gone if expired.

## Session

### DELETE /me/data

GDPR erasure — permanently deletes all bookmarks and history for the current session or API key.

**Response:**

```json
{
  "status": "purged",
  "deleted": {"bookmarks": 42, "history": 150}
}
```

## Next Steps

- **[Authentication](/docs/api/authentication)** - API key setup
- **[WebSocket](/docs/api/websocket)** - Real-time batch updates
- **[Error Handling](/docs/api/error-handling)** - Error responses
- **[Rate Limits](/docs/api/rate-limits)** - Rate limit details
- **[Interactive Docs](http://localhost:8001/api/v1/docs)** - Try the API
