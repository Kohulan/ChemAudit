# ChemAudit Batch Export Formats

Nine export formats available at `GET /api/v1/batch/{job_id}/export?format=<fmt>` (or POST for large index selections).

## Common filters

All formats accept: `score_min`, `score_max`, `status` (`success`/`error`/`warning`), `indices` (comma-separated for GET, array body for POST).

---

## `csv`

Flat CSV with one row per molecule. Columns include:

- `index`, `name`, `smiles`, `canonical_smiles`, `inchikey`
- `status`, `error`
- Validation: `overall_score`, `issues` (semicolon-joined check names), `num_issues`
- Scoring (when enabled at upload): `qed`, `sa_score`, `mw`, `logp`, `hbd`, `hba`, `tpsa`, `rotatable_bonds`, `num_aromatic_rings`, `fsp3`, `lipinski_passed`
- Alerts (when enabled): `alert_pains`, `alert_brenk`, `alert_nih`, `alert_zinc`, `alert_count`
- Standardization (when `include_standardization=true`): `standardized_smiles`, `std_changes`

Media type: `text/csv`. Extension: `.csv`.

---

## `excel`

Query params:
- `include_images=true|false` — embed 2D depictions of each molecule (slower, larger files).
- `sheet_layout=single|multi` — one sheet with everything, or one sheet per section (Validation / Scoring / Alerts / Standardization / Properties).

Styling: conditional colouring on score columns (green ≥70, amber 40–69, red <40), frozen header row, column-width auto-fit.

Media type: `application/vnd.openxmlformats-officedocument.spreadsheetml.sheet`. Extension: `.xlsx`.

---

## `sdf`

MOL blocks concatenated with `$$$$` terminators. Each molecule has properties attached:

- `<SMILES>`, `<InChIKey>`, `<Validation_Score>`, `<Issues>` (newline-joined)
- `<QED>`, `<SA_Score>`, `<MW>`, `<LogP>`, `<TPSA>` (when scoring enabled)
- `<PAINS_Alerts>`, `<BRENK_Alerts>`, etc. (when alerts enabled)

Query param `include_audit=true` adds full JSON-serialized validation and scoring payloads as single-line properties for downstream parsing.

Media type: `chemical/x-mdl-sdfile`. Extension: `.sdf`.

---

## `json`

Newline-delimited array of full result objects. Structure per molecule:

```json
{
  "index": 0,
  "smiles": "...",
  "name": null,
  "status": "success",
  "validation": {"overall_score": 87, "issues": [...], "all_checks": [...]},
  "alerts": {...},
  "scoring": {"druglikeness": {...}, "admet": {...}, "np_likeness": {...}},
  "standardization": {...},
  "profiling": {...},
  "safety_assessment": {...}
}
```

Includes everything the backend stored — the most faithful representation.

Media type: `application/json`. Extension: `.json`.

---

## `pdf`

Query params:
- `sections=validation_summary,score_distribution,alert_breakdown,top_issues,property_histograms` — pick a subset. Default is all sections.
- `include_audit=true` — add a per-molecule audit appendix.

Contains:
- Executive summary (counts, pass rate, processing time).
- Score distribution histogram.
- Alert breakdown bar chart by catalog.
- Property histograms (MW, LogP, TPSA, QED, SA score).
- Top-N failing-check frequency table.
- Optional: per-molecule thumbnail cards with SMILES, score, issues.

Media type: `application/pdf`. Extension: `.pdf`.

---

## `fingerprint`

ZIP archive containing fingerprint vectors computed from the canonical SMILES:

- `morgan_r2_2048.csv`, `morgan_r2_2048.npy`, `morgan_r2_2048.npz` — Morgan radius=2, 2048 bits.
- `maccs.csv`, `maccs.npy`, `maccs.npz` — 166-bit MACCS keys.
- `rdkit_default.csv`, `rdkit_default.npy`, `rdkit_default.npz` — RDKit path-based fingerprint.
- `index.csv` — mapping row index → `smiles`, `inchikey`, `name`.

`.csv` = human-readable one-row-per-molecule with bit string; `.npy` = dense NumPy; `.npz` = sparse CSR.

Media type: `application/zip`. Extension: `.zip`.

---

## `dedup`

ZIP with deduplication analysis keyed by InChIKey:

- `summary.csv` — one row per unique InChIKey with cluster size and representative SMILES.
- `duplicate_groups/<inchikey>.csv` — per-group CSV listing all rows, their original indices, and any property deltas.
- `README.txt` — summary counts (input, unique, duplicates removed).

Media type: `application/zip`. Extension: `.zip`.

---

## `scaffold`

Single CSV with Murcko scaffold analysis:

- Columns: `index`, `smiles`, `scaffold_smiles` (Murcko), `generic_scaffold_smiles` (carbon-only skeleton), `scaffold_group_id`, `scaffold_group_size`.
- Molecules sharing the same Murcko scaffold receive the same `scaffold_group_id` (1-indexed, ordered by frequency).

Media type: `text/csv`. Extension: `.csv`.

---

## `property_matrix`

ZIP with molecule × property matrices:

- `properties_flat.csv` — molecules as rows, every numeric property as a column.
- `properties_multi.xlsx` — Excel with one sheet per property family (Physicochemical / Druglikeness / ADMET / Safety).
- `schema.json` — column definitions (units, source, computation method).

Media type: `application/zip`. Extension: `.zip`.

---

## Subset export (hand-picked indices)

POST when indices would exceed URL length on GET:

```bash
curl -sS -X POST "http://localhost:8000/api/v1/batch/<job_id>/export?format=excel&include_images=true" \
  -H 'Content-Type: application/json' \
  -d '{"indices": [0, 3, 17, 42, 101, 303, ...]}' -o selected.xlsx
```

Also available as a dedicated endpoint:

```bash
curl -sS -X POST "http://localhost:8000/api/v1/batch/<job_id>/subset/export" \
  -H 'Content-Type: application/json' \
  -d '{"indices": [0, 3, 17], "format": "sdf"}' -o subset.sdf
```
