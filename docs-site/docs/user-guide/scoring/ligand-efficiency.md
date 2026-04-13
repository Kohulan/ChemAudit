---
sidebar_position: 10
title: Ligand Efficiency
description: Ligand efficiency (LE) and lipophilic ligand efficiency (LLE) metrics for lead optimization
---

# Ligand Efficiency

Ligand efficiency metrics quantify how effectively a molecule uses its size to achieve binding potency. These are key metrics in lead optimization, helping medicinal chemists prioritize compounds that are potent relative to their molecular size and lipophilicity.

## Metrics

### Ligand Efficiency (LE)

**Formula:** LE = pActivity / N_HA

Where:
- **pActivity** = −log₁₀(IC₅₀ or Kd in molar)
- **N_HA** = number of heavy (non-hydrogen) atoms

| LE Value | Interpretation |
|----------|----------------|
| **≥ 0.3** | Efficient — acceptable for drug-like leads |
| **0.2–0.3** | Moderate — may benefit from optimization |
| **< 0.2** | Inefficient — consider fragment-based starting points |

:::info Proxy Activity
When no experimental activity data is available, ChemAudit uses the QED (Quantitative Estimate of Drug-likeness) score as a proxy pIC50 value. The response indicates when a proxy was used via the `proxy_used` field.
:::

### Lipophilic Ligand Efficiency (LLE)

**Formula:** LLE = pActivity − LogP

LLE (also known as LipE) separates potency from lipophilicity. A high LLE indicates that binding potency is driven by specific molecular interactions rather than non-specific hydrophobic binding.

| LLE Value | Interpretation |
|-----------|----------------|
| **≥ 5** | Excellent — potency is not driven by lipophilicity |
| **3–5** | Good — reasonable balance of potency and lipophilicity |
| **< 3** | Potency may be primarily driven by lipophilicity |

## Why It Matters

During lead optimization, molecular weight and lipophilicity tend to increase ("molecular obesity"). Tracking LE and LLE helps:

- **Prioritize efficient leads**: Compounds with high LE use fewer atoms to achieve the same potency
- **Avoid lipophilicity traps**: High LLE compounds are less likely to have off-target effects or poor ADMET properties
- **Guide SAR decisions**: When adding substituents, check whether LE is maintained or declining

## API Usage

Request ligand efficiency as part of the scoring endpoint:

```bash
curl -X POST http://localhost:8001/api/v1/score \
  -H "Content-Type: application/json" \
  -d '{
    "molecule": "c1ccc(cc1)c2cc(nn2c3ccc(cc3)S(=O)(=O)N)C(F)(F)F",
    "include": ["ligand_efficiency"]
  }'
```

**Response:**

```json
{
  "molecule_info": { "..." : "..." },
  "ligand_efficiency": {
    "le": 0.28,
    "heavy_atom_count": 27,
    "activity_value": 7.5,
    "activity_type": "pIC50",
    "proxy_used": true
  },
  "execution_time_ms": 12
}
```

| Field | Type | Description |
|-------|------|-------------|
| `le` | float | Ligand Efficiency value |
| `heavy_atom_count` | int | Number of heavy atoms |
| `activity_value` | float | pIC50 value (real or proxy) |
| `activity_type` | string | Activity type used |
| `proxy_used` | bool | Whether QED was used as a proxy |

## References

- Hopkins, A. L., Groom, C. R. & Alex, A. (2004). Ligand efficiency: a useful metric for lead selection. *Drug Discovery Today*, 9(10), 430–431.
- Leeson, P. D. & Springthorpe, B. (2007). The influence of drug-like concepts on decision-making in medicinal chemistry. *Nature Reviews Drug Discovery*, 6(11), 881–890.

## Next Steps

- **[Bioavailability Radar](/docs/user-guide/scoring/bioavailability-permeation)** — Oral bioavailability prediction
- **[Drug-Likeness](/docs/user-guide/scoring/drug-likeness)** — Lipinski, QED, and consensus scoring
- **[ADMET](/docs/user-guide/scoring/admet)** — Absorption, distribution, metabolism, excretion predictions
