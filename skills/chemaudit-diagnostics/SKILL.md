---
name: chemaudit-diagnostics
description: >
  Diagnose SMILES parse errors with position and fix suggestions, compare InChI
  strings layer-by-layer, check format round-trip lossiness (SMILES→InChI→SMILES
  and SMILES→MOL→SMILES), compare three standardization pipelines side-by-side,
  and pre-validate SDF/CSV files for structural integrity. Use when user says
  "why does this SMILES fail", "diagnose this structure", "InChI diff", "layer
  comparison", "round-trip check", "compare InChI strings", "file pre-validation",
  "SDF integrity", "M END missing", or "fix this SMILES error".
license: MIT
allowed-tools: "Bash(curl:*) Bash(python:*)"
compatibility: >
  Requires a running ChemAudit instance. No external dependencies beyond the API.
metadata:
  author: Kohulan Rajan
  version: 1.0.0
  mcp-server: chemaudit
  category: cheminformatics
  tags: [SMILES, InChI, diagnostics, error-detection, file-validation]
---

# ChemAudit Diagnostics

## Overview

Five targeted diagnostic endpoints. Each answers a different "why didn't this work?" question.

| Endpoint | Answers |
|---|---|
| `POST /api/v1/diagnostics/smiles` | "Why does this SMILES fail to parse?" |
| `POST /api/v1/diagnostics/inchi-diff` | "How do these two InChI strings differ?" |
| `POST /api/v1/diagnostics/roundtrip` | "Does this SMILES survive a round-trip through InChI (or MOL)?" |
| `POST /api/v1/diagnostics/cross-pipeline` | "Do three different standardization pipelines agree?" |
| `POST /api/v1/diagnostics/file-prevalidate` | "Is this SDF/CSV structurally sound before batch upload?" |

## Workflow

### 1. SMILES error diagnostics (DIAG-01)

Uses a dual strategy: RDKit `DetectChemistryProblems` for parseable-but-problematic SMILES, log-capture for unparseable SMILES. Returns error type, character position, and ranked fix suggestions.

```bash
curl -sS -X POST http://localhost:8000/api/v1/diagnostics/smiles \
  -H 'Content-Type: application/json' \
  -d '{"smiles": "C(C)(C)(C)(C)C"}'
```

Response:

```json
{
  "valid": false,
  "canonical_smiles": null,
  "warnings": [],
  "errors": [
    {
      "raw_message": "Explicit valence for atom ... is 5, is greater than permitted",
      "position": 0,
      "error_type": "valence_error",
      "message": "Carbon at position 0 has 5 bonds (max: 4)",
      "suggestions": [
        {
          "description": "Remove one neighbor branch",
          "corrected_smiles": "C(C)(C)(C)C",
          "confidence": 0.9
        }
      ]
    }
  ]
}
```

Error types returned: `unmatched_bracket`, `valence_error`, `ring_closure_mismatch`, `unknown_atom_symbol`, `invalid_charge`, `parse_error`.

For valid-but-suspicious SMILES, `valid: true` with populated `warnings[]` (RDKit chemistry warnings, e.g. kekulization issues).

### 2. InChI layer diff (DIAG-02)

Pure string comparison — no RDKit. Parses each InChI into its constituent layers and produces a per-layer diff table.

```bash
curl -sS -X POST http://localhost:8000/api/v1/diagnostics/inchi-diff \
  -H 'Content-Type: application/json' \
  -d '{
    "inchi_a": "InChI=1S/C9H8O4/c1-6(10)13-8-5-3-2-4-7(8)9(11)12/h2-5H,1H3,(H,11,12)",
    "inchi_b": "InChI=1S/C9H8O4/c1-6(10)13-8-5-3-2-4-7(8)9(11)12/h2-5H,1H3,(H,11,12)/t7-/m1/s1"
  }'
```

Response:

```json
{
  "identical": false,
  "layer_rows": [
    {"layer": "formula", "value_a": "C9H8O4", "value_b": "C9H8O4", "match": true},
    {"layer": "connections", "value_a": "c1-6(...)", "value_b": "c1-6(...)", "match": true},
    {"layer": "hydrogens", "value_a": "h2-5H,1H3,(H,11,12)", "value_b": "h2-5H,1H3,(H,11,12)", "match": true},
    {"layer": "stereo_tetrahedral", "value_a": null, "value_b": "t7-", "match": false},
    {"layer": "stereo_parity", "value_a": null, "value_b": "m1", "match": false},
    {"layer": "stereo_marker", "value_a": null, "value_b": "s1", "match": false}
  ],
  "layers_a": {"formula": "C9H8O4", "connections": "...", "hydrogens": "..."},
  "layers_b": {"formula": "C9H8O4", "connections": "...", "hydrogens": "...", "stereo_tetrahedral": "t7-", ...}
}
```

Great for answering "same compound, stereo-defined vs racemic?" at a glance — rows with `match=false` isolate the exact disagreement.

### 3. Format round-trip (DIAG-03)

Does the molecule survive conversion to an intermediate and back?

```bash
# SMILES → InChI → SMILES (detects stereo/isotope loss)
curl -sS -X POST http://localhost:8000/api/v1/diagnostics/roundtrip \
  -H 'Content-Type: application/json' \
  -d '{"smiles": "C[C@H](N)C(=O)O", "route": "smiles_inchi_smiles"}'

# SMILES → MOL block → SMILES (detects stereo/charge loss)
curl -sS -X POST http://localhost:8000/api/v1/diagnostics/roundtrip \
  -H 'Content-Type: application/json' \
  -d '{"smiles": "...", "route": "smiles_mol_smiles"}'
```

Response:

```json
{
  "route": "smiles_inchi_smiles",
  "original_smiles": "C[C@H](N)C(=O)O",
  "intermediate": "InChI=1S/C3H7NO2/c1-2(4)3(5)6/h2H,4H2,1H3,(H,5,6)/t2-/m0/s1",
  "roundtrip_smiles": "C[C@H](N)C(=O)O",
  "lossy": false,
  "losses": [],
  "error": null
}
```

When lossy:

```json
{
  "lossy": true,
  "losses": [
    {"type": "stereo", "description": "2 stereocenters lost", "before": 2, "after": 0},
    {"type": "isotope", "description": "Isotope label stripped", "before": 1, "after": 0}
  ]
}
```

Loss types: `stereo`, `charge`, `isotope`.

### 4. Cross-pipeline standardization comparison (DIAG-04)

Runs the molecule through **three** pipelines and compares outputs. Useful for picking which pipeline to use on a new dataset, or diagnosing why standardized SMILES disagree between sources.

```bash
curl -sS -X POST http://localhost:8000/api/v1/diagnostics/cross-pipeline \
  -H 'Content-Type: application/json' \
  -d '{"molecule": "CC(=O)Oc1ccccc1C(=O)[O-].[Na+]", "format": "auto"}'
```

Response:

```json
{
  "pipelines": [
    {
      "name": "RDKit MolStandardize",
      "smiles": "CC(=O)Oc1ccccc1C(=O)O",
      "inchikey": "BSYNRYMUTXBXSQ-UHFFFAOYSA-N",
      "mw": 180.157,
      "formula": "C9H8O4",
      "charge": 0,
      "stereo_count": 0,
      "error": null,
      "highlight_atoms": [],
      "highlight_bonds": []
    },
    {"name": "ChEMBL Pipeline", ...},
    {"name": "Minimal", "smiles": "CC(=O)Oc1ccccc1C(=O)[O-].[Na+]", ...}
  ],
  "disagreements": 3,
  "structural_disagreements": 1,
  "all_agree": false,
  "property_comparison": [
    {"property": "smiles", "values": [...], "agrees": false, "structural": true},
    {"property": "inchikey", "values": [...], "agrees": false, "structural": true},
    {"property": "mw", "values": [180.157, 180.157, 203.14], "agrees": false, "structural": false},
    ...
  ]
}
```

`highlight_atoms` / `highlight_bonds` mark atoms outside the cross-pipeline MCS — useful for rendering a side-by-side visual diff in the UI.

### 5. File pre-validation (DIAG-05)

Catches file-level problems before you pay the cost of a batch upload.

```bash
curl -sS -X POST http://localhost:8000/api/v1/diagnostics/file-prevalidate \
  -F "file=@compounds.sdf"
```

Response (SDF):

```json
{
  "file_type": "sdf",
  "total_blocks": 1000,
  "total_rows": null,
  "encoding": null,
  "issue_count": 3,
  "issues": [
    {"block": 42, "line": 1250, "issue_type": "missing_m_end", "severity": "error",
     "description": "Block 42 is missing the 'M  END' terminator"},
    {"block": 87, "line": 2310, "issue_type": "malformed_count_line", "severity": "error",
     "description": "Block 87 has a malformed counts line"},
    {"block": null, "line": null, "issue_type": "suspicious_content", "severity": "warning",
     "description": "Found pattern that could indicate embedded code"}
  ],
  "valid": false
}
```

Response (CSV):

```json
{
  "file_type": "csv",
  "total_blocks": null,
  "total_rows": 500,
  "encoding": "utf-8",
  "issue_count": 1,
  "issues": [
    {"block": null, "line": null, "issue_type": "missing_smiles_column", "severity": "error",
     "description": "No column matches 'SMILES' (case-insensitive)"}
  ],
  "valid": false
}
```

Issue types: `missing_m_end`, `malformed_count_line`, `encoding_fallback`, `encoding_error`, `suspicious_content`, `missing_smiles_column`, `empty_rows`, `empty_file`, `duplicate_columns`.

Severity levels: `error`, `warning`, `info`. `valid=false` only when at least one `error` severity appears — warnings don't block validity.

Max upload size: 50 MB (hard-coded in this endpoint).

## Examples

### Example 1 — "This SMILES fails to parse; what's wrong?"

1. `POST /diagnostics/smiles` with the bad SMILES.
2. Read `errors[0].error_type` and `message` — tells you *what*.
3. Read `errors[0].position` — tells you *where* in the string.
4. Try `errors[0].suggestions[0].corrected_smiles` — first suggestion is ranked highest confidence.

### Example 2 — "I have two InChI strings for the 'same' compound — are they really the same?"

1. `POST /diagnostics/inchi-diff` with both.
2. Check `identical` — if `true`, you're done.
3. If `false`, scan `layer_rows[]` for rows with `match=false`. Formula differences mean different compounds; only stereo/isotope differences mean "same skeleton, different specification".

### Example 3 — "Will my dataset survive saving as an SDF?"

1. For each distinct SMILES, `POST /diagnostics/roundtrip` with `route="smiles_mol_smiles"`.
2. Any entry with `lossy=true` means something will be lost when saving. The `losses[]` array tells you what.
3. Stereo loss is the most common — MOL v2000 handles it poorly for some R/S cases. If it matters, save as SDF v3000 or export as SMILES instead.

### Example 4 — "Pre-flight an SDF before a 10k-molecule batch"

1. `POST /diagnostics/file-prevalidate` with the file.
2. If `valid=true`, go straight to `/batch/upload`.
3. If `valid=false`, fix each `severity=error` issue; warnings can be ignored at your discretion. Common fix: append `M  END` to blocks missing it, or re-save from RDKit.

## Rate limits

- `/diagnostics/smiles`, `/diagnostics/inchi-diff`, `/diagnostics/roundtrip`: 30/min.
- `/diagnostics/cross-pipeline`, `/diagnostics/file-prevalidate`: 10/min.

## Troubleshooting

### 400 "Invalid characters in SMILES string"

The validator blocks `< > & ; | $ \``. These chars are never valid in SMILES — strip client-side.

### `errors[]` is empty but `valid=false`

Rare but possible for cases where the parse error is in the RDKit C++ layer with no accompanying log message. Try `/standardize` to see if the ChEMBL pipeline's error message is more informative.

### Round-trip shows lossy for SMILES that "should" work

Check `losses[].type`:
- `stereo`: the intermediate format doesn't encode your stereo (e.g. atropisomerism).
- `isotope`: InChI's optional isotope layer was dropped in the round trip (shouldn't happen with standard InChI but occasionally does on edge cases).
- `charge`: MOL v2000 loses some charge states.

### Cross-pipeline: all three disagree

Usually means the input is unusual (organometallic, radical, mixture). Treat each pipeline's output as "one opinion of many" and pick the one matching the source database's convention (ChEMBL for ChEMBL-derived data, RDKit for exotic cases).

### File pre-validation: encoding_fallback warning

The file wasn't valid UTF-8; the validator fell back to Latin-1 and succeeded. Usually fine but worth noting — re-save as UTF-8 if you'll reuse the file.

### File pre-validation: suspicious_content

Pattern-matched against script/macro patterns. False positives happen on legitimate SMILES containing `<`, `>`, or similar characters (rare). Audit-logged server-side regardless.

## Further reading

- `chemaudit-molecule-validation` — for check-level issues once a SMILES parses.
- `chemaudit-standardization` — for the ChEMBL-style standardization pipeline.
- `chemaudit-batch-validation` — pair `/file-prevalidate` with `/batch/upload`.
