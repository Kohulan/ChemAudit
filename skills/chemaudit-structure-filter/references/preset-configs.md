# Structure Filter Preset Configurations

Source of truth: `frontend/src/types/structure_filter.ts` (`PRESET_CONFIGS`) and
`backend/app/services/structure_filter/filter_config.py` (`PRESETS`). Values below
match the locked D-12 thresholds and D-15 differentiated weight vectors.

## `drug_like`

Balanced default for general drug-discovery screening. Balanced composite weights.

```json
{
  "min_mw": 200, "max_mw": 500,
  "min_logp": -1, "max_logp": 5,
  "max_tpsa": 140,
  "max_rot_bonds": 10,
  "max_rings": null,
  "max_sa_score": 5.0,
  "use_pains": true, "use_brenk": true, "use_kazius": true, "use_nibr": false,
  "enable_novelty": false, "novelty_threshold": 0.85,
  "weight_validity": 0.3, "weight_qed": 0.3, "weight_alert_free": 0.2, "weight_sa": 0.2
}
```

**Rationale**: Lipinski-compatible MW/LogP, Veber-compatible TPSA/RotBonds, PAINS + Brenk + Kazius on. NIBR off because it's aggressive for general HTS libraries. Balanced weights to reward overall drug-likeness without over-weighting any single component.

## `lead_like`

Tighter than drug-like — used during lead optimization where MW and lipophilicity should drift upward as potency is added. QED weighting emphasized.

```json
{
  "min_mw": 200, "max_mw": 350,
  "min_logp": -1, "max_logp": 3.5,
  "max_tpsa": 140,
  "max_rot_bonds": 7,
  "max_rings": null,
  "max_sa_score": 4.0,
  "use_pains": true, "use_brenk": true, "use_kazius": true, "use_nibr": true,
  "enable_novelty": false, "novelty_threshold": 0.85,
  "weight_validity": 0.2, "weight_qed": 0.4, "weight_alert_free": 0.2, "weight_sa": 0.2
}
```

**Rationale**: Hann & Keserű lead-like space (MW 200–350, LogP ≤ 3.5, RotBonds ≤ 7) gives room for lead→drug MW growth. All alert catalogs including NIBR on because Novartis deck filters catch many problems that survive PAINS/Brenk. Weight QED higher because drug-likeness is the lead-optimization goal.

## `fragment_like`

Astex Rule of 3 for fragment-based drug discovery. Emphasizes clean building blocks (alert-free, easy to synthesize).

```json
{
  "min_mw": 100, "max_mw": 300,
  "min_logp": -1, "max_logp": 3,
  "max_tpsa": 100,
  "max_rot_bonds": 3,
  "max_rings": 3,
  "max_sa_score": 3.0,
  "use_pains": true, "use_brenk": false, "use_kazius": false, "use_nibr": false,
  "enable_novelty": false, "novelty_threshold": 0.85,
  "weight_validity": 0.2, "weight_qed": 0.2, "weight_alert_free": 0.3, "weight_sa": 0.3
}
```

**Rationale**: Tight MW/LogP window; ring cap of 3 keeps fragments small; SA ≤ 3.0 ensures synthesizability. Only PAINS on because Brenk and Kazius patterns are rare at fragment size and false-positive rates climb. Weights emphasize alert-free and SA — the two most important qualities for a fragment-screen library.

## `permissive`

Escape hatch for "just check parseability and deduplicate". Alerts all off.

```json
{
  "min_mw": 100, "max_mw": 800,
  "min_logp": -5, "max_logp": 8,
  "max_tpsa": 200,
  "max_rot_bonds": 15,
  "max_rings": null,
  "max_sa_score": 7.0,
  "use_pains": false, "use_brenk": false, "use_kazius": false, "use_nibr": false,
  "enable_novelty": false, "novelty_threshold": 0.85,
  "weight_validity": 0.4, "weight_qed": 0.3, "weight_alert_free": 0.1, "weight_sa": 0.2
}
```

**Rationale**: Most-permissive property windows; skip the alert stages entirely; SA ≤ 7 excludes only the hardest-to-synthesize molecules. Heavy validity weight because with alerts off, "it parsed and sanitized" is the main signal.

## Config field reference

| Field | Type | Default (drug_like) | Meaning |
|---|---|---|---|
| `min_mw` / `max_mw` | float | 200 / 500 | Molecular weight window (Da) |
| `min_logp` / `max_logp` | float | -1 / 5 | Crippen LogP window |
| `max_tpsa` | float | 140 | Topological polar surface area (Å²) |
| `max_rot_bonds` | int | 10 | Rotatable bonds |
| `max_rings` | int \| null | null | Ring count cap; `null` = no limit |
| `max_sa_score` | float | 5.0 | Ertl SA score (1 easy, 10 hard) |
| `use_pains` | bool | true | Screen PAINS A/B/C catalogs |
| `use_brenk` | bool | true | Screen Brenk unwanted-reactivity catalog |
| `use_kazius` | bool | true | Screen Kazius mutagenicity toxicophores |
| `use_nibr` | bool | false | Screen NIBR Novartis deck filters |
| `enable_novelty` | bool | false | Add ChEMBL Tanimoto novelty stage |
| `novelty_threshold` | float | 0.85 | Max Tanimoto to any ChEMBL compound (reject above) |
| `weight_validity` | float | 0.3 | Composite weight: validity component |
| `weight_qed` | float | 0.3 | Composite weight: QED component |
| `weight_alert_free` | float | 0.2 | Composite weight: alert-free binary |
| `weight_sa` | float | 0.2 | Composite weight: inverted SA score |

Weights do not need to sum to 1; the score is normalized internally. But keeping them normalized makes comparisons across presets clearer.

## Choosing between preset and custom config

**Use a preset when**:
- You want the built-in scientifically-motivated defaults.
- You want differentiated weight vectors without computing them yourself.
- Reproducibility with other REINVENT users matters (preset name is self-documenting).

**Use a custom config when**:
- Your project has specific MW/LogP constraints (e.g. oncology tends to relax LogP).
- You need novelty filtering (no preset enables it).
- You want specific alert catalogs only (e.g. Brenk on, Kazius off).
