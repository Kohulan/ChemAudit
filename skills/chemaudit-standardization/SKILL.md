---
name: chemaudit-standardization
description: >
  Standardize chemical structures using ChemAudit's ChEMBL-compatible pipeline
  (Checker → Standardizer → GetParent → optional Tautomer) with full per-stage
  provenance, charge/bond/ring-change tracking, and stereochemistry comparison.
  Use when user says "standardize this molecule", "normalize structure", "strip
  salts", "neutralize charges", "canonical SMILES", "ChEMBL standardization",
  "parent molecule", "tautomer canonicalization", or "show me what changed
  after standardization". Shows every transformation applied with before/after
  InChIKey and atom-level diffs.
license: MIT
allowed-tools: "Bash(curl:*) Bash(python:*) Bash(chemaudit:*)"
compatibility: >
  Requires a running ChemAudit instance, or use the CLI with --local flag for
  offline mode (no HTTP server required).
metadata:
  author: Kohulan Rajan
  version: 1.0.0
  mcp-server: chemaudit
  category: cheminformatics
  tags: [standardization, ChEMBL, salt-stripping, normalization, provenance]
---

# ChemAudit Standardization

## Overview

ChEMBL-compatible standardization pipeline with four stages:

1. **Checker** — runs `ChEMBL structure checker`; detects structural problems with penalty scores (penalty>0 means the checker flagged something).
2. **Standardizer** — fixes nitro groups, metal bonds, sulphoxides, etc. via normalization SMARTS rules.
3. **GetParent** — extracts the largest organic fragment (removes salts, counterions, solvents).
4. **Tautomer** (optional, off by default) — canonicalizes tautomers via `TautomerEnumerator`. **Warning**: may destroy E/Z double-bond stereochemistry.

Optional detailed provenance records atom-level charge/bond/ring changes per stage.

## Access modes

- **API**: `POST /api/v1/standardize`.
- **CLI**: `chemaudit standardize --smiles "CCO"` (server) or `--local` (offline).
- **Cross-pipeline diagnostic**: `POST /api/v1/diagnostics/cross-pipeline` runs three pipelines (RDKit MolStandardize, ChEMBL-style, minimal) and compares outputs.

## Workflow

### 1. Default standardization

```bash
curl -sS -X POST http://localhost:8000/api/v1/standardize \
  -H 'Content-Type: application/json' \
  -d '{"molecule": "CC(=O)Oc1ccccc1C(=O)[O-].[Na+]", "format": "auto"}'
```

Response:

```json
{
  "molecule_info": {
    "input_smiles": "CC(=O)Oc1ccccc1C(=O)[O-].[Na+]",
    "canonical_smiles": "...",
    "inchikey": "...",
    "molecular_formula": "C9H7NaO4",
    "molecular_weight": 203.14
  },
  "result": {
    "original_smiles": "CC(=O)Oc1ccccc1C(=O)[O-].[Na+]",
    "standardized_smiles": "CC(=O)Oc1ccccc1C(=O)O",
    "success": true,
    "steps_applied": [
      {"step_name": "checker", "applied": true, "description": "ChEMBL structure checker", "changes": "Found 1 issue"},
      {"step_name": "standardizer", "applied": true, "description": "Normalize functional groups", "changes": "No changes"},
      {"step_name": "get_parent", "applied": true, "description": "Remove salts and counterions", "changes": "Removed 1 fragment: Na+"}
    ],
    "checker_issues": [
      {"penalty_score": 5, "message": "Salt form detected"}
    ],
    "excluded_fragments": ["[Na+]"],
    "stereo_comparison": {
      "before_count": 0, "after_count": 0, "lost": 0, "gained": 0,
      "double_bond_stereo_lost": 0, "warning": null
    },
    "structure_comparison": {
      "original_atom_count": 14, "standardized_atom_count": 13,
      "original_formula": "C9H7NaO4", "standardized_formula": "C9H8O4",
      "original_mw": 203.14, "standardized_mw": 180.16,
      "mass_change_percent": -11.31, "is_identical": false,
      "diff_summary": ["Removed: Na+"]
    },
    "mass_change_percent": -11.31,
    "provenance": null
  },
  "execution_time_ms": 87
}
```

### 2. Enable tautomer canonicalization (opt-in)

```bash
curl -sS -X POST http://localhost:8000/api/v1/standardize \
  -H 'Content-Type: application/json' \
  -d '{
    "molecule": "Oc1nc(N)nc2c1ncn2",
    "options": {
      "include_tautomer": true,
      "preserve_stereo": true
    }
  }'
```

Adds a `tautomer_canonicalization` step. Watch `stereo_comparison.double_bond_stereo_lost` — tautomer enumeration can collapse E/Z bonds.

### 3. Full provenance (atom-level diff)

```bash
curl -sS -X POST http://localhost:8000/api/v1/standardize \
  -H 'Content-Type: application/json' \
  -d '{
    "molecule": "CC(=O)Oc1ccccc1C(=O)[O-].[Na+]",
    "options": {
      "include_provenance": true
    }
  }'
```

The `result.provenance` field fills with per-stage records:

```json
{
  "stages": [
    {
      "stage_name": "standardizer",
      "input_smiles": "...", "output_smiles": "...",
      "applied": true,
      "charge_changes": [
        {"atom_idx": 11, "element": "O", "before_charge": -1, "after_charge": 0, "rule_name": "carboxylate", "smarts": "..."}
      ],
      "bond_changes": [...],
      "radical_changes": [],
      "ring_changes": [],
      "fragment_removals": [],
      "dval_cross_refs": []
    },
    {
      "stage_name": "get_parent",
      "applied": true,
      "fragment_removals": [
        {"smiles": "[Na+]", "name": "Sodium", "role": "counterion", "mw": 22.99}
      ],
      ...
    }
  ],
  "tautomer": null,
  "stereo_summary": {
    "stereo_stripped": false, "centers_lost": 0, "bonds_lost": 0,
    "per_center": [], "dval_cross_refs": []
  }
}
```

Cross-references to Deep Validation check IDs (DVAL-01, DVAL-03) appear in `dval_cross_refs` when the molecule fails those checks — useful for understanding whether standardization addressed a flagged issue.

### 4. Pass prior validation results for cross-referencing

If you already ran deep validation and want the standardization provenance to cross-reference DVAL IDs:

```bash
curl -sS -X POST http://localhost:8000/api/v1/standardize \
  -H 'Content-Type: application/json' \
  -d '{
    "molecule": "...",
    "options": {
      "include_provenance": true,
      "dval_results": {
        "undefined_stereo": {"count": 2},
        "tautomer_detection": {"count": 1}
      }
    }
  }'
```

### 5. Get available options

```bash
curl -sS http://localhost:8000/api/v1/standardize/options
```

Returns the canonical list of pipeline steps, their defaults, and warnings (e.g. tautomer stereo-loss caveat).

### 6. Compare three pipelines side-by-side

```bash
curl -sS -X POST http://localhost:8000/api/v1/diagnostics/cross-pipeline \
  -H 'Content-Type: application/json' \
  -d '{"molecule": "CC(=O)Oc1ccccc1C(=O)[O-].[Na+]", "format": "auto"}'
```

Runs:
- RDKit `MolStandardize` (Cleanup + LargestFragment + Uncharger + TautomerEnumerator).
- ChEMBL-style pipeline (the same as `/standardize`).
- Minimal pipeline (parse + sanitize only).

Compares SMILES, InChIKey, MW, formula, charge, stereo count across all three. Surfaces disagreements — critical when picking a pipeline for a new dataset.

### 7. CLI usage

```bash
chemaudit standardize --smiles "CC(=O)Oc1ccccc1C(=O)[O-].[Na+]"
chemaudit standardize --smiles "..." --local --format json   # offline
chemaudit standardize --file compounds.csv                    # first line
echo "CCO" | chemaudit standardize                            # stdin
```

`--local` uses `app.services.standardization.chembl_pipeline` directly — no HTTP. Outputs a Rich table by default, JSON when piped or `--format json`.

## Stereo comparison

`result.stereo_comparison` summarizes:

| Field | Meaning |
|---|---|
| `before_count` / `after_count` | Defined stereocenters before and after |
| `lost` / `gained` | Net change in defined stereocenters |
| `double_bond_stereo_lost` | E/Z bonds lost to tautomerization |
| `warning` | Populated when `lost > 0` or `double_bond_stereo_lost > 0` |

`preserve_stereo: true` (default) makes the pipeline attempt to retain stereochemistry through each stage. Tautomer canonicalization is the main offender; the warning surfaces when it strips stereo.

## Examples

### Example 1 — "Standardize and show what changed"

1. `POST /standardize` (defaults).
2. Print `result.original_smiles` vs `result.standardized_smiles`.
3. Loop `result.steps_applied[]`; print each step's `description` + `changes`.
4. Print `result.excluded_fragments[]` for the salt/solvent story.
5. Show `structure_comparison.mass_change_percent` — usually negative (lost a counterion).

### Example 2 — "Why do these three SMILES have different InChIKeys?"

1. `POST /diagnostics/cross-pipeline` for each — read each pipeline's InChIKey.
2. Look at `disagreements` and `property_comparison[]` to see which property differs.
3. If `structural_disagreements > 0`, the three standardization approaches produced different molecular graphs — investigate the per-pipeline `error` fields.

### Example 3 — "Strip salts and neutralize without touching tautomers"

Defaults (`include_tautomer: false`) already do this. The pipeline runs Checker → Standardizer (normalizes functional groups) → GetParent (salt stripping) — leaving tautomers alone. The result's canonical SMILES is safe to use as a dedup key without risking E/Z loss.

### Example 4 — "Preserve stereochemistry in a hard case"

Pass `{"options": {"preserve_stereo": true, "include_tautomer": false}}` — default. If you need tautomer canonicalization and stereo, accept that some E/Z information may be lost; the response surfaces exactly how much in `stereo_comparison.double_bond_stereo_lost`.

## Rate limits

- `/standardize`: 10/min.
- `/standardize/options`: 30/min.
- `/diagnostics/cross-pipeline`: 10/min.

API-key tier raises to 300/min.

## Troubleshooting

### `result.success = false`

The pipeline raised an exception; `result.error_message` has details. Common causes:
- Unparseable input (also surfaces as 400 with `format_detected`).
- Tautomer enumeration timed out (100+ tautomers) — the engine flags `complexity_flag=true`.
- Pipeline chose a fragment with unusual chemistry (rare; open an issue with the InChIKey).

### `standardized_smiles = null` but `success = true`

Stripping all fragments left nothing. Means the input was entirely a salt/counterion with no parent (e.g. `[Na+].[Cl-]`). Not an error — genuinely no parent exists.

### `mass_change_percent = 0` but `steps_applied` shows applied stages

The standardizer applied normalizations that don't change mass (e.g. resonance structure normalization, redrawing nitro groups). Check `structure_comparison.diff_summary[]` for the qualitative story.

### `double_bond_stereo_lost > 0` after enabling tautomer

Expected. Tautomer canonicalization breaks and reforms double bonds. If E/Z geometry matters, set `include_tautomer: false` — accept that different tautomeric inputs will produce different canonical SMILES.

### InChIKey changed but structures "look the same"

Likely tautomer canonicalization, or a charge/protonation normalization. Run with `include_provenance: true` and inspect `charge_changes[]` and `bond_changes[]` to see the atom-level story.

### Checker issues with `penalty_score = 0`

Zero-penalty entries are informational. Non-zero penalties mean the ChEMBL checker considers the input substandard.

## Further reading

- `chemaudit-qsar-ready` — when you want the full 10-step ML-curation pipeline, not just the 3-4 ChEMBL steps.
- `chemaudit-diagnostics` — for cross-pipeline comparison and SMILES round-trip checking.
- `chemaudit-molecule-validation` — for deep checks that complement standardization.
