---
sidebar_position: 5
title: ADMET
description: Predict Absorption, Distribution, Metabolism, Excretion, and Toxicity properties
---

# ADMET Predictions

ADMET scoring predicts key pharmacokinetic and physicochemical properties: Absorption, Distribution, Metabolism, Excretion, and Toxicity.

## Available Predictions

| Property | Method | Output | Interpretation |
|----------|--------|--------|----------------|
| **Synthetic Accessibility** | SAscore | 1-10 (1=easy, 10=hard) | Ease of synthesis |
| **Aqueous Solubility** | ESOL (Delaney) | LogS, mg/mL, classification | Water solubility |
| **Complexity** | Fsp3, stereocenters, rings, Bertz CT | Various metrics | Molecular complexity |
| **CNS MPO** | Multi-parameter optimization | Score 0-6 | CNS penetration likelihood |
| **Bioavailability** | TPSA, rotatable bonds, Lipinski | Predictions | Oral absorption, CNS flags |
| **Pfizer 3/75 Rule** | LogP < 3, TPSA > 75 | Pass/fail | Toxicity risk reduction |
| **GSK 4/400 Rule** | MW ≤ 400, LogP ≤ 4 | Pass/fail | Lead-like properties |
| **Golden Triangle** | MW vs LogD plot | In/out | Optimal drug-like space |

## Synthetic Accessibility (SAscore)

Estimates how difficult a molecule is to synthesize:

| Score | Classification | Interpretation |
|-------|---------------|----------------|
| **1-3** | Easy | Simple synthesis, few steps |
| **4-6** | Moderate | Standard synthetic methods |
| **7-9** | Difficult | Complex, many steps |
| **10** | Very Difficult | Extremely challenging |

Based on fragment contributions and complexity penalties.

## Aqueous Solubility (ESOL)

Predicts water solubility using the Delaney ESOL linear regression model.

### ESOL Equation

```
LogS = 0.16 − 0.63 × LogP − 0.0062 × MW + 0.066 × RotBonds − 0.74 × AP
```

Where:
- **LogS** = log₁₀(aqueous solubility in mol/L)
- **LogP** = Wildman-Crippen LogP
- **MW** = molecular weight (Da)
- **RotBonds** = rotatable bond count
- **AP** = aromatic proportion (aromatic atoms / heavy atoms)

**Conversion:** `solubility_mg_mL = 10^LogS × MW / 1000`

### Classification

| LogS | Category | Approx. mg/mL |
|------|----------|---------------|
| **≥ −1** | Highly soluble | > 100 |
| **−1 to −3** | Soluble | 1–100 |
| **−3 to −4** | Moderately soluble | 0.1–1 |
| **−4 to −5** | Poorly soluble | < 0.1 |
| **< −5** | Insoluble | Very low |

Good aqueous solubility (LogS > −4) is favorable for oral drugs.

**Reference:** Delaney (2004). ESOL: Estimating aqueous solubility directly from molecular structure. *JCICS*, 44(3), 1000–1005.

## Molecular Complexity

Multiple complexity metrics:

**Fsp3** (Fraction sp3 carbons):
- Higher Fsp3 = more saturated, better 3D character
- < 0.25: Flat/aromatic
- 0.25–0.42: Moderate 3D
- \> 0.42: Good 3D character (target for drug-likeness)

**Reference (Fsp3):** Lovering et al. (2009). Escape from flatland. *Journal of Medicinal Chemistry*, 52(21), 6752–6756.

**Stereocenters:**
- Count of chiral centers
- More stereocenters = higher synthetic complexity

**Ring systems:**
- Total rings and aromatic rings
- Complex scaffolds have many rings

**Bertz CT** (Complexity index):
- Higher values = more complex
- Accounts for branching and symmetry

## CNS MPO (Multiparameter Optimization)

Pfizer's Central Nervous System Multiparameter Optimization score predicts CNS penetration likelihood on a 0–6 scale.

### Component Scoring

Each component scores 0–1, and the total is the sum of all 6 components:

| Parameter | Score = 1.0 | Linear decrease | Score = 0 |
|-----------|-------------|----------------|-----------|
| **MW** | ≤ 360 Da | 360–500 Da | > 500 Da |
| **LogP** | ≤ 3.0 | 3.0–5.0 | > 5.0 |
| **TPSA** | ≤ 40 A² | 40–90 A² | > 90 A² |
| **HBD** | 0 | 1–3 (−0.25 per donor) | > 3 |
| **LogD** | 0.5 (estimated) | — | — |
| **pKa** | 0.5 (estimated) | — | — |

:::info Placeholder Values
LogD and pKa components use placeholder values (0.5 each) as these require experimental data or more complex models to calculate accurately.
:::

### Interpretation

| CNS MPO Score | Interpretation |
|---------------|---------------|
| **≥ 5** | Excellent CNS penetration |
| **4–5** | Good CNS penetration |
| **3–4** | Moderate |
| **< 3** | Poor CNS penetration |

**Reference:** Wager et al. (2010). Moving beyond rules: the development of a CNS MPO approach. *ACS Chemical Neuroscience*, 1(6), 435–449.

## Bioavailability Predictions

**Oral Absorption:**

```
oral_absorption_likely = lipinski_ok AND veber_ok
```

Where `lipinski_ok` = MW ≤ 500, LogP ≤ 5, HBD ≤ 5, HBA ≤ 10 and `veber_ok` = rotatable bonds ≤ 10, TPSA ≤ 140.

**CNS Penetration:**

```
cns_penetration_likely = TPSA ≤ 90 AND MW ≤ 450 AND HBD ≤ 3 AND 1 ≤ LogP ≤ 4
```

## Pfizer 3/75 Rule

A toxicity risk indicator:

```
at_risk = LogP > 3 AND TPSA < 75
```

Compounds meeting both criteria have statistically higher rates of toxicity and off-target promiscuity.

**Reference:** Hughes et al. (2008). Physiochemical drug properties associated with in vivo toxicological outcomes. *Bioorganic & Medicinal Chemistry Letters*, 18(17), 4872–4875.

## GSK 4/400 Rule

GlaxoSmithKline's compound quality guideline:

```
favorable = MW ≤ 400 AND LogP ≤ 4
```

Provides room for optimization while maintaining favorable ADMET properties.

**Reference:** Gleeson (2008). Generation of a set of simple, interpretable ADMET rules of thumb. *Journal of Medicinal Chemistry*, 51(4), 817–834.

## Golden Triangle

Abbott's optimal property space for balanced permeability and metabolic stability:

```
in_triangle = 200 ≤ MW ≤ 450 AND −0.5 ≤ LogP ≤ 5
```

Molecules within this space tend to have favorable permeability-metabolism balance.

**Reference:** Johnson et al. (2009). Using the Golden Triangle to optimize clearance and oral absorption. *Bioorganic & Medicinal Chemistry Letters*, 19(17), 5560–5564.

## API Usage

```bash
curl -X POST http://localhost:8001/api/v1/score \
  -H "Content-Type: application/json" \
  -d '{
    "molecule": "CC(=O)Oc1ccccc1C(=O)O",
    "include": ["admet"]
  }'
```

Response:
```json
{
  "admet": {
    "synthetic_accessibility": {
      "score": 1.5,
      "classification": "Easy",
      "interpretation": "Very easy to synthesize"
    },
    "solubility": {
      "log_s": -2.1,
      "solubility_mg_ml": 1.43,
      "classification": "Soluble",
      "interpretation": "Good aqueous solubility"
    },
    "complexity": {
      "fsp3": 0.11,
      "num_stereocenters": 0,
      "num_rings": 1,
      "num_aromatic_rings": 1,
      "bertz_ct": 245.2,
      "classification": "Low complexity"
    },
    "cns_mpo": {
      "score": 3.8,
      "cns_penetrant": false
    },
    "bioavailability": {
      "oral_absorption_likely": true,
      "cns_penetration_likely": false
    },
    "pfizer_rule": {
      "passed": true,
      "logp": 1.19,
      "tpsa": 63.6
    },
    "gsk_rule": {
      "passed": true,
      "mw": 180.16,
      "logp": 1.19
    },
    "golden_triangle": {
      "in_golden_triangle": true,
      "mw": 180.16,
      "logd": 1.19
    }
  }
}
```

## Interpretation Guidelines

**Good ADMET Profile:**
- SAscore < 5 (easy to synthesize)
- LogS > -4 (good solubility)
- Fsp3 > 0.4 (sufficient saturation)
- Oral absorption likely
- Passes Pfizer and GSK rules
- In golden triangle

**Poor ADMET Profile:**
- SAscore > 7 (hard to synthesize)
- LogS < -6 (poorly soluble)
- Very complex (many stereocenters, rings)
- CNS penetration when not desired
- Fails multiple rules

## Limitations

All predictions are computational estimates:

- Based on training data (mostly drug-like molecules)
- May not generalize to unusual scaffolds
- Don't replace experimental measurements
- Use for prioritization, not absolute decisions

:::warning Experimental Validation Required
ADMET predictions guide early decisions but must be validated experimentally for lead candidates.
:::

## Use Cases

### Lead Optimization

Track ADMET during optimization:

- Monitor solubility changes
- Avoid increasing synthetic complexity
- Maintain favorable bioavailability
- Stay in golden triangle

### Compound Prioritization

Rank compounds by ADMET profile:

```
SAscore < 5 AND
LogS > -4 AND
oral_absorption_likely = true AND
pfizer_rule_passed = true
```

### Library Design

Design screening libraries with good ADMET:

- Target SAscore < 5
- Ensure solubility > -4
- Maintain Fsp3 > 0.4
- Stay in golden triangle

## Next Steps

- **[Drug-likeness](/docs/user-guide/scoring/drug-likeness)** - Lipinski, QED, Veber rules
- **[Scoring Overview](/docs/user-guide/scoring/overview)** - All scoring systems
- **[API Reference](/docs/api/endpoints)** - Full scoring API
