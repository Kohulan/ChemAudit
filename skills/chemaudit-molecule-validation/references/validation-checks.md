# ChemAudit Deep Validation Checks

The 16 deep validation checks registered at backend startup (see
`EXPECTED_DEEP_VALIDATION_CHECKS` in `backend/app/main.py`). Each check returns
`{check_name, passed, severity, message, affected_atoms, details}` via
`POST /api/v1/validate`.

Pass specific names in the `checks` field, or `["all"]` for every check.

## M1.1 — Stereo & Tautomer

| Check name | DVAL ID | What it flags |
|---|---|---|
| `stereoisomer_enumeration` | DVAL-01/02 | Undefined stereocenters; reports the count of possible stereoisomers |
| `tautomer_detection` | DVAL-03 | Multiple tautomers possible (input may not be the canonical form) |
| `aromatic_system_validation` | DVAL-04 | Incorrect aromaticity assignment or Kekulé inconsistencies |
| `coordinate_dimension` | DVAL-05 | Missing or 2D-only coordinates when 3D are required |

## M1.2 — Chemical Composition

| Check name | DVAL ID | What it flags |
|---|---|---|
| `mixture_detection` | DVAL-06 | Multiple disconnected components (likely salt or mixture) |
| `solvent_contamination` | DVAL-07 | Common solvents (water, DMSO, MeOH, EtOH, etc.) as separate fragments |
| `inorganic_filter` | DVAL-08 | No carbon atoms in the largest fragment |
| `radical_detection` | DVAL-09 | Unpaired electrons (radicals) on atoms |
| `isotope_label_detection` | DVAL-10 | Non-default isotopes (²H, ¹³C, ¹⁵N, etc.) |
| `trivial_molecule` | DVAL-11 | Fewer heavy atoms than a minimum threshold (noise) |

## M1.3 — Structural Complexity

| Check name | DVAL ID | What it flags |
|---|---|---|
| `hypervalent_atoms` | DVAL-12 | Atoms exceeding standard valence (e.g. 5-valent carbon) |
| `polymer_detection` | DVAL-13 | Polymer repeat markers or implausibly large homopolymers |
| `ring_strain` | DVAL-14 | Highly strained small rings (3- and 4-membered) in unusual contexts |
| `macrocycle_detection` | DVAL-15 | Rings of ≥ 12 atoms |
| `charged_species` | DVAL-16 | Net charge ≠ 0; lists per-atom formal charges |
| `explicit_hydrogen_audit` | DVAL-17 | Unexpected explicit H atoms that should be implicit |

## Severity semantics

- `info` — reported for transparency (e.g. isotope labels), not a quality issue on its own.
- `warning` — resolvable curation issue (mixture, unstandardized tautomer, undefined stereo).
- `error` — structural problem that usually means the molecule should not enter an ML dataset as-is (hypervalent atoms, polymer markers, radicals).

Score mapping: each failing check removes points in the scoring engine according
to its severity. Final `overall_score` is clamped to [0, 100]. See
`backend/app/services/validation/engine.py` for the exact weights.

## Common fix recipes

| Failing check | Typical fix |
|---|---|
| `mixture_detection` | Run the `chemaudit-standardization` pipeline; `get_parent` keeps the largest fragment. |
| `solvent_contamination` | Same — standardization strips common solvents automatically. |
| `tautomer_detection` | Enable `include_tautomer` in `/standardize` **only** if losing E/Z stereo is acceptable. |
| `stereoisomer_enumeration` | Assign stereochemistry explicitly, or drop stereo via QSAR-Ready pipeline with `enable_stereo_strip: true`. |
| `charged_species` | `/standardize` with default `Uncharger` neutralizes where chemistry permits. |
| `inorganic_filter` | Exclude from ML dataset or confirm the compound is intentional (metal complexes). |
| `radical_detection` | Usually an input error — re-parse from the original source. |
| `polymer_detection` | Exclude; polymers are out of scope for small-molecule ML models. |
