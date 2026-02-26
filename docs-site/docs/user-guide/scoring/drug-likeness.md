---
sidebar_position: 3
title: Drug Likeness
description: Evaluate compliance with Lipinski, QED, Veber, and other drug-likeness rules
---

# Drug-Likeness Scoring

Drug-likeness scores evaluate whether a molecule has properties consistent with known drugs. ChemAudit implements multiple established drug-likeness filters.

## Available Filters

| Filter | Criteria | Pass Rate in DrugBank |
|--------|----------|----------------------|
| **Lipinski (Ro5)** | MW ≤ 500, LogP ≤ 5, HBD ≤ 5, HBA ≤ 10 | ~85% |
| **QED** | Quantitative Estimate (0-1 scale) | Mean ~0.6 |
| **Veber** | Rotatable bonds ≤ 10, TPSA ≤ 140 | ~90% |
| **Rule of Three** | MW ≤ 300, LogP ≤ 3, HBD ≤ 3, HBA ≤ 3 | Fragment-like |
| **Ghose** | MW 160-480, LogP -0.4-5.6, atoms 20-70, MR 40-130 | ~75% |
| **Egan** | LogP ≤ 5.88, TPSA ≤ 131.6 | ~85% |
| **Muegge** | Multiple criteria | ~80% |

## Lipinski's Rule of Five

The most widely used drug-likeness filter:

| Property | Criterion | Rationale |
|----------|-----------|-----------|
| **Molecular Weight** | ≤ 500 Da | Oral bioavailability correlation |
| **LogP** | ≤ 5 | Membrane permeability |
| **H-Bond Donors** | ≤ 5 | Solubility and permeability |
| **H-Bond Acceptors** | ≤ 10 | Solubility and permeability |

**Violations allowed:** Up to 1 violation is often acceptable

:::info Exceptions
Many successful drugs violate Lipinski: antibiotics, antifungals, natural products, and biologics-inspired molecules.
:::

## QED (Quantitative Estimate of Drug-likeness)

Continuous score from 0 to 1 based on 8 molecular properties:

- Molecular weight
- LogP
- H-bond donors
- H-bond acceptors
- Polar surface area
- Rotatable bonds
- Aromatic rings
- Structural alerts

| QED Score | Interpretation |
|-----------|---------------|
| **0.8-1.0** | Highly drug-like |
| **0.6-0.8** | Moderately drug-like |
| **0.4-0.6** | Low drug-likeness |
| **0.0-0.4** | Not drug-like |

## Veber Rules

Focus on oral bioavailability:

| Property | Criterion | Purpose |
|----------|-----------|---------|
| **Rotatable Bonds** | ≤ 10 | Molecular flexibility |
| **TPSA** | ≤ 140 Å² | Membrane permeability |

Developed from analysis of rat oral bioavailability data.

## Rule of Three (Ro3)

Fragment-like criteria for screening libraries:

- MW ≤ 300 Da
- LogP ≤ 3
- H-bond donors ≤ 3
- H-bond acceptors ≤ 3

Used for fragment-based drug discovery.

## Ghose Filter

Defines drug-like chemical space:

- MW: 160-480 Da
- LogP: -0.4 to 5.6
- Atom count: 20-70
- Molar refractivity: 40-130

Based on analysis of the World Drug Index.

## Egan Rules

Simple bioavailability filter:

- LogP ≤ 5.88
- TPSA ≤ 131.6 Å²

Derived from 1,2000 compounds with known oral bioavailability.

## Muegge Filter

Comprehensive 9-parameter filter for druglike chemical space:

| Property | Threshold |
|----------|-----------|
| Molecular Weight | 200–600 Da |
| LogP | −2 to 5 |
| TPSA | ≤ 150 A² |
| Ring Count | ≤ 7 |
| Carbon Count | > 4 |
| Heteroatom Count | > 1 |
| Rotatable Bonds | ≤ 15 |
| H-Bond Donors | ≤ 5 |
| H-Bond Acceptors | ≤ 10 |

**Pass criterion:** All 9 parameters must pass (0 violations).

**Reference:** Muegge et al. (2001). Simple selection criteria for drug-like chemical matter. *Journal of Medicinal Chemistry*, 44(12), 1841–1846.

## Consensus Drug-Likeness (0–5)

A composite score counting how many of the major filter sets the molecule passes:

1. Lipinski (≤ 1 violation)
2. Veber (both parameters pass)
3. Egan (both parameters pass)
4. Ghose (all parameters in range)
5. Muegge (all 9 parameters pass)

| Score | Interpretation |
|-------|----------------|
| **5** | Excellent drug-likeness |
| **4** | Very good |
| **3** | Acceptable |
| **2** | Borderline |
| **0–1** | Poor drug-likeness |

## Lead-Likeness

Identifies molecules in lead-like chemical space, which provides room for optimization during hit-to-lead campaigns.

| Property | Threshold |
|----------|-----------|
| Molecular Weight | 200–350 Da |
| LogP | −1 to 3 |
| Rotatable Bonds | ≤ 7 |

**Pass criterion:** All must pass.

## API Usage

```bash
curl -X POST http://localhost:8001/api/v1/score \
  -H "Content-Type: application/json" \
  -d '{
    "molecule": "CC(=O)Oc1ccccc1C(=O)O",
    "include": ["druglikeness"]
  }'
```

Response:
```json
{
  "druglikeness": {
    "lipinski": {
      "passed": true,
      "violations": 0,
      "mw": 180.16,
      "logp": 1.19,
      "hbd": 1,
      "hba": 4
    },
    "qed": {
      "score": 0.71,
      "interpretation": "Moderately drug-like"
    },
    "veber": {
      "passed": true,
      "rotatable_bonds": 3,
      "tpsa": 63.6
    },
    "ro3": {
      "passed": true,
      "violations": 0
    },
    "ghose": {
      "passed": true,
      "violations": 0
    },
    "egan": {
      "passed": true
    },
    "muegge": {
      "passed": true,
      "violations": 0
    },
    "interpretation": "Passes Lipinski and Veber rules"
  }
}
```

## Interpreting Results

### All Filters Pass

Excellent drug-likeness. Molecule has properties consistent with most known oral drugs.

### Lipinski Passes, Others Fail

Reasonable drug-likeness. Review specific failures to understand limitations.

### Multiple Failures

May still be drug-like if:
- Targeting specific therapeutic areas (antibiotics often fail)
- Natural product (often fail due to complexity)
- Prodrug or special formulation

### Low QED (less than 0.4)

Indicates poor overall drug-likeness, but context matters. Some successful drugs have low QED.

## Use Cases

### Virtual Screening

Pre-filter compound libraries:

```
lipinski_passed = true AND
veber_passed = true AND
qed >= 0.5
```

### Lead Optimization

Track drug-likeness during optimization:

- Monitor QED score changes
- Avoid introducing Lipinski violations
- Maintain favorable TPSA and rotatable bonds

### Library Design

Design screening libraries with good drug-like properties:

- Ghose filter for general drug-like space
- Ro3 for fragment libraries
- Lipinski + Veber for lead-like compounds

## Limitations

**These rules are guidelines, not absolutes:**

- Many drugs violate these rules
- Rules derived from oral drugs (not applicable to all routes)
- Don't account for specific mechanisms or targets
- Don't predict efficacy or safety

:::warning Context Matters
Use drug-likeness filters for prioritization, not as strict cutoffs. Always consider therapeutic area, target, and route of administration.
:::

## Best Practices

1. **Use multiple filters**: No single filter is perfect
2. **Consider therapeutic area**: Some areas (antibiotics, antivirals) often violate rules
3. **Track QED over time**: Monitor during optimization
4. **Don't over-optimize**: Chasing perfect drug-likeness can sacrifice potency
5. **Validate experimentally**: Predicted properties should be confirmed

## References

1. Bickerton, G. R. et al. (2012). Quantifying the chemical beauty of drugs. *Nature Chemistry*, 4(2), 90–98.
2. Congreve, M. et al. (2003). A 'rule of three' for fragment-based lead discovery? *Drug Discovery Today*, 8(19), 876–877.
3. Egan, W. J. et al. (2000). Prediction of drug absorption using multivariate statistics. *Journal of Medicinal Chemistry*, 43(21), 3867–3877.
4. Ghose, A. K. et al. (1999). A knowledge-based approach in designing combinatorial or medicinal chemistry libraries. *Journal of Combinatorial Chemistry*, 1(1), 55–68.
5. Lipinski, C. A. et al. (2001). Experimental and computational approaches to estimate solubility and permeability. *Advanced Drug Delivery Reviews*, 46(1-3), 3–26.
6. Muegge, I. et al. (2001). Simple selection criteria for drug-like chemical matter. *Journal of Medicinal Chemistry*, 44(12), 1841–1846.
7. Veber, D. F. et al. (2002). Molecular properties that influence the oral bioavailability of drug candidates. *Journal of Medicinal Chemistry*, 45(12), 2615–2623.

## Next Steps

- **[Scoring Overview](/docs/user-guide/scoring/overview)** — All scoring systems
- **[ADMET](/docs/user-guide/scoring/admet)** — ADMET predictions
- **[Safety Filters](/docs/user-guide/scoring/safety-filters)** — Structural alerts
