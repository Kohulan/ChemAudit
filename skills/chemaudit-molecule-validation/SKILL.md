---
name: chemaudit-molecule-validation
description: >
  Validate and score chemical structures using ChemAudit's 16 deep validation
  checks, 1,500+ structural alerts, and multi-rule drug-likeness scoring. Use
  when user says "validate this molecule", "check this SMILES", "is this
  compound valid", "ML-readiness score", "drug-likeness", "deep validation",
  "quality score", "PAINS check", or asks about stereochemistry, valence,
  tautomers, or structural issues. Supports SMILES, InChI, MOL blocks, IUPAC
  names (OPSIN + PubChem), and database identifiers (ChEMBL ID, CAS, PubChem
  CID, InChIKey) via the /resolve endpoint.
license: MIT
allowed-tools: "Bash(curl:*) Bash(python:*) Bash(chemaudit:*)"
compatibility: >
  Requires a running ChemAudit instance (docker compose up) or use the CLI
  with --local flag for offline mode. Python 3.11+.
metadata:
  author: Kohulan Rajan
  version: 1.0.0
  mcp-server: chemaudit
  category: cheminformatics
  tags: [chemistry, validation, drug-discovery, SMILES, molecular-properties]
---

# ChemAudit Molecule Validation

## Overview

Single-molecule validation with three complementary layers:

1. **Validation** (`/validate`) — 16 deep cheminformatics checks, 0-100 score.
2. **Scoring** (`/score`) — drug-likeness (Lipinski, QED, Veber, Ghose, Egan, Muegge, Ro3), ADMET, ML-readiness, NP-likeness, safety filters.
3. **Alerts** (`/alerts`, `/alerts/screen`) — PAINS, Brenk, NIH, ZINC, Kazius, NIBR, plus 7 ChEMBL pharma filter sets.

Score interpretation: **≥70 pass, 40–69 warn, <40 fail**.

## Access modes

- **API**: HTTP calls to a running ChemAudit instance (default `http://localhost:8000`).
- **CLI**: `chemaudit validate|score|standardize|profile` with `--local` for offline mode.
- **MCP**: `chemaudit` MCP server exposes all tagged endpoints at `/mcp`.

Input formats auto-detected: SMILES, InChI, MOL block, IUPAC name (OPSIN then PubChem fallback). Identifier resolution (ChEMBL ID, CAS, PubChem CID, ChEBI, UNII, DrugBank, etc.) goes through `/api/v1/resolve` first.

## Workflow

### 1. Quick validation (default 16 checks + scoring)

```bash
curl -sS -X POST http://localhost:8000/api/v1/validate \
  -H 'Content-Type: application/json' \
  -d '{"molecule": "CC(=O)Oc1ccccc1C(=O)O", "format": "auto"}'
```

Response (abridged):

```json
{
  "status": "completed",
  "molecule_info": {
    "canonical_smiles": "CC(=O)OC1=CC=CC=C1C(=O)O",
    "inchikey": "BSYNRYMUTXBXSQ-UHFFFAOYSA-N",
    "molecular_formula": "C9H8O4",
    "molecular_weight": 180.16,
    "num_stereocenters": 0,
    "has_stereochemistry": false
  },
  "overall_score": 95,
  "issues": [],
  "all_checks": [...],
  "execution_time_ms": 42,
  "cached": false
}
```

### 2. Targeted deep checks

Pass specific check names instead of `"all"` (the default). Full list in `references/validation-checks.md`:

```bash
curl -sS -X POST http://localhost:8000/api/v1/validate \
  -H 'Content-Type: application/json' \
  -d '{
    "molecule": "CC(C)(C)C1=CC=C(C=C1)[N+](=O)[O-]",
    "checks": ["stereoisomer_enumeration", "tautomer_detection",
               "hypervalent_atoms", "radical_detection", "charged_species"]
  }'
```

Each failing check returns `severity` (info/warning/error), `message`, `affected_atoms`, and a `details` dict.

### 3. Drug-likeness and ADMET scoring

```bash
curl -sS -X POST http://localhost:8000/api/v1/score \
  -H 'Content-Type: application/json' \
  -d '{
    "molecule": "CC(=O)Oc1ccccc1C(=O)O",
    "include": ["druglikeness", "admet", "np_likeness", "ml_readiness",
                "safety_filters", "consensus", "bioavailability_radar",
                "boiled_egg"]
  }'
```

Druglikeness returns Lipinski (0–4 violations), QED (0–1), Veber, Ro3, Ghose, Egan, Muegge. ADMET returns SAscore, ESOL solubility, Fsp3, CNS MPO, Pfizer 3/75, GSK 4/400, Golden Triangle. ML-readiness returns 0–100 with 4-dimension breakdown.

### 4. Scoring profile presets (8 built-in)

Profiles live in the database and are seeded at startup. Retrieve or use by ID:

```bash
curl -sS http://localhost:8000/api/v1/profiles
```

The 8 presets are: **Drug-like (Lipinski)**, **Lead-like**, **Fragment-like (Rule of 3)**, **CNS-penetrant**, **Ghose (Amgen)**, **Veber (GSK)**, **PPI-like**, **NP-like**. Presets are immutable — duplicate before editing.

### 5. Structural alerts

Targeted screen against a chosen set of catalogs:

```bash
curl -sS -X POST http://localhost:8000/api/v1/alerts \
  -H 'Content-Type: application/json' \
  -d '{
    "molecule": "CC(=O)Oc1ccccc1C(=O)O",
    "catalogs": ["PAINS", "BRENK", "NIH", "ZINC"]
  }'
```

Unified one-shot screen across **all** catalogs (PAINS A/B/C, BRENK, NIH, ZINC, Kazius, NIBR, 7 ChEMBL pharma filters, 21 custom SMARTS) with deduplicated concern groups:

```bash
curl -sS -X POST http://localhost:8000/api/v1/alerts/screen \
  -H 'Content-Type: application/json' \
  -d '{"molecule": "CCOC(=O)N=Nc1ccccc1"}'
```

List available catalogs:

```bash
curl -sS http://localhost:8000/api/v1/alerts/catalogs
```

**Critical caveat**: alerts are investigation flags, not rejections. 87 FDA-approved drugs contain PAINS patterns.

### 6. Compound profiling (PFI, #stars, bioavailability, SA comparison)

```bash
curl -sS -X POST http://localhost:8000/api/v1/profiler/full \
  -H 'Content-Type: application/json' \
  -d '{"smiles": "CC(=O)Oc1ccccc1C(=O)O"}'
```

Returns PFI, Abbott #stars, consensus LogP, skin permeation, CNS MPO, SA comparison (SA Score + SCScore + SYBA when installed). 3D shape descriptors are lazy — call `/profiler/shape-3d` separately.

### 7. CLI usage

```bash
chemaudit validate --smiles "CC(=O)Oc1ccccc1C(=O)O"
chemaudit validate --smiles "CCO" --local                # offline
chemaudit validate --file compounds.csv --format json
echo "CCO" | chemaudit validate                           # stdin
chemaudit score --smiles "CCO" --local
chemaudit profile --smiles "CCO"
```

`--local` bypasses the HTTP server and calls service functions directly.

### 8. Identifier resolution (CAS, ChEMBL ID, PubChem CID, names)

Before validating a compound given as a non-structural identifier, resolve it first:

```bash
curl -sS -X POST http://localhost:8000/api/v1/resolve \
  -H 'Content-Type: application/json' \
  -d '{"identifier": "50-78-2", "identifier_type": "auto"}'
```

Supported types: `auto`, `smiles`, `inchi`, `inchikey`, `pubchem_cid`, `chembl_id`, `cas`, `drugbank_id`, `chebi_id`, `unii`, `wikipedia`, `name`. Response includes `canonical_smiles`, `cross_references` (CIDs, ChEMBL, DrugBank, UNII, CAS), and `resolution_chain` provenance.

## Examples

### Example 1 — "Validate this SMILES and tell me if it's ML-ready"

1. `POST /api/v1/validate` with `{"molecule": "<smi>"}` → check `overall_score`, list `issues`.
2. `POST /api/v1/score` with `{"molecule": "<smi>", "include": ["ml_readiness", "druglikeness"]}` → read `ml_readiness.score` (0–100) and `interpretation`.
3. Summarise: parseable, passes Lipinski with X violations, ML-readiness = Y/100 with breakdown across 4 dimensions.

### Example 2 — "Run deep validation and explain all issues"

1. `POST /api/v1/validate` with `{"molecule": "<smi>", "checks": ["all"]}`.
2. Filter `all_checks` where `passed=false`, group by `severity`.
3. For each issue, explain the `check_name` using `references/validation-checks.md`, show `affected_atoms` and the `details` dict.

### Example 3 — "Score against Lipinski, Ghose, and CNS-penetrant profiles"

1. `GET /api/v1/profiles` → find IDs for `Drug-like (Lipinski)`, `Ghose (Amgen)`, `CNS-penetrant`.
2. `POST /api/v1/score` with `{"molecule": "<smi>", "include": ["druglikeness"]}` gets Lipinski and Ghose directly.
3. For CNS-penetrant, use `/profiler/full` → read `cns_mpo.score` and `cns_mpo.cns_penetrant` boolean.

## Rate limits

Anonymous (no API key):
- 10/min: `/validate`, `/checks`, `/score`, `/score/compare`, `/alerts`, `/alerts/catalogs`, `/alerts/quick-check`, `/alerts/screen`, `/standardize`.
- 30/min: `/validate/async`, `/validate/similarity`, `/standardize/options`, `/profiler/full`, `/profiler/efficiency`, `/profiler/mpo`.
- 20/min: `/profiler/sa-comparison`.
- 10/min: `/profiler/shape-3d`.

API key via `X-API-Key` header raises the tier (300/min on most endpoints).

## Troubleshooting

### 400 "Failed to parse molecule"

The detail dict includes `errors[]`, `warnings[]`, and `format_detected`. Common fixes:
- Unbalanced ring closures — use `POST /api/v1/diagnostics/smiles` for position-specific fix suggestions.
- IUPAC name with unusual punctuation — pass `"format": "iupac"` explicitly to force OPSIN/PubChem.
- MOL block — ensure `M  END` terminator is present and on its own line.

### 400 "Invalid characters in molecule string"

The validator blocks `< > & ; | $ \``. Strip these client-side before calling.

### Input auto-classified as IUPAC when it's a valid SMILES

Short valid SMILES like `CO`, `O`, `CCO` are ambiguous. Pass `"input_type": "smiles"` or `"format": "smiles"` to skip name resolution.

### 504 on `/validate/async`

The high-priority Celery worker is unavailable. Either use the sync `/validate` endpoint or check that the `high_priority` queue worker is running (`docker compose logs worker-priority`).

### 429 rate limit

Use `X-API-Key` header for the 300/min tier, or stagger requests client-side. Batch endpoints (`/batch/upload`) have their own limits — see the `chemaudit-batch-validation` skill.

### InChIKey-based caching

`/validate` caches by `(inchikey, checks)`; repeated calls return `cached: true` with the same score. To bypass, add a meaningless suffix to `checks` or flush Redis.

## Further reading

- `references/validation-checks.md` — catalog of all 16 deep validation checks with severity, rationale, and fix recipes.
- `chemaudit-batch-validation` — upload CSV/SDF for bulk validation with analytics.
- `chemaudit-diagnostics` — when you need to know *why* a SMILES fails, not just that it did.
- `chemaudit-standardization` — before validating heterogeneous datasets, standardize first.
