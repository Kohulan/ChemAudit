---
sidebar_position: 8
title: Aggregator Likelihood
description: Predict colloidal aggregation risk in biochemical assays
---

# Aggregator Likelihood

Aggregator likelihood scoring predicts whether a molecule is likely to form colloidal aggregates in biochemical assays, which can lead to false positive results.

## What Are Aggregators?

Aggregators are small molecules that form colloidal particles in aqueous solution, leading to:

- Non-specific enzyme inhibition
- False positives in biochemical assays
- Artifacts in high-throughput screening
- Misleading structure-activity relationships

Aggregation is concentration-dependent and can be disrupted by detergents.

## How It's Calculated

The aggregator risk score (0–1) is the sum of triggered risk increments, capped at 1.0:

### Risk Indicators

| Indicator | Condition | Risk Increment |
|-----------|-----------|---------------|
| High lipophilicity | LogP > 4.0 | +0.30 |
| Moderate lipophilicity | LogP 3.0–4.0 | +0.15 |
| Low polar surface | TPSA < 40 A² | +0.20 |
| Extended aromatics | ≥ 4 aromatic rings | +0.20 |
| Moderate aromatics | 3 aromatic rings | +0.10 |
| Large molecular weight | MW > 500 Da | +0.10 |
| High conjugation | Aromatic fraction > 0.7 and > 20 heavy atoms | +0.15 |
| Known scaffold match | SMARTS pattern match | +0.20 per match |

### Known Aggregator Scaffolds (SMARTS-based)

The following empirically known aggregator scaffolds are matched:

- Rhodanines
- Quinones
- Catechols
- Curcumin-like structures
- Phenol-sulfonamides
- Flavonoids
- Acridines
- Phenothiazines
- Extended conjugated polyenes
- Azo compounds

**Confidence:** `triggered_indicators / total_indicators`

## Output

The score provides:

| Field | Description |
|-------|-------------|
| `likelihood` | Risk category: Low, Medium, High |
| `risk_score` | Numerical score (0.0-1.0) |
| `risk_factors` | List of contributing factors |
| `interpretation` | Human-readable explanation |

## Risk Categories

| Likelihood | Risk Score | Interpretation |
|------------|------------|----------------|
| **Low** | 0.0-0.3 | Unlikely to aggregate |
| **Medium** | 0.3-0.7 | May aggregate at high concentration |
| **High** | 0.7-1.0 | Likely aggregator |

## API Usage

```bash
curl -X POST http://localhost:8001/api/v1/score \
  -H "Content-Type: application/json" \
  -d '{
    "molecule": "c1ccc2c(c1)ccc3c2ccc4c3cccc4",
    "include": ["aggregator"]
  }'
```

Response:
```json
{
  "aggregator": {
    "likelihood": "High",
    "risk_score": 0.85,
    "logp": 6.2,
    "tpsa": 0.0,
    "mw": 228.29,
    "aromatic_rings": 4,
    "risk_factors": [
      "High LogP (6.2)",
      "Low TPSA (0.0)",
      "Many aromatic rings (4)"
    ],
    "interpretation": "High risk of aggregation - likely false positive in assays"
  }
}
```

## Interpreting Results

### Low Risk

Molecule is unlikely to aggregate:

- Proceed with assay screening
- Standard assay conditions should work
- False positive risk minimal

### Medium Risk

Molecule may aggregate at high concentrations:

- Use appropriate concentration range
- Include detergent in assay (0.01% Triton X-100 or similar)
- Verify hits with dose-response curves
- Check for non-specific effects

### High Risk

Molecule likely aggregates:

- High risk of false positives
- Always include detergent in assay
- Perform counterscreening
- Use light scattering or turbidity checks
- Consider alternative screening methods
- Verify hits with orthogonal assays

## Common Risk Factors

**High LogP (> 5):**

Hydrophobic molecules partition poorly in water and tend to self-associate.

**Low TPSA (< 40 Å²):**

Poor polar surface area means low water solubility, promoting aggregation.

**Multiple Aromatic Rings (> 3):**

Planar aromatic systems stack via pi-pi interactions, forming aggregates.

**Moderate Molecular Weight (200-500 Da):**

The size range where aggregation is most common.

## Use Cases

### Hit Validation

Check aggregation risk before spending resources on hits:

```
aggregator_likelihood != "High"
```

### Assay Design

Design assays with aggregation in mind:

- Include detergent for high-risk compounds
- Use dynamic light scattering to detect aggregates
- Run turbidity checks

### Library Curation

Filter aggregators from screening libraries:

```
aggregator_likelihood = "Low" AND
risk_score < 0.3
```

### Dose-Response Analysis

Interpret dose-response curves:

- Aggregators often show Hill slopes > 1
- Non-specific, steep inhibition
- Detergent disrupts activity

## Experimental Validation

Test for aggregation:

**Detergent sensitivity:**
- Run assay with and without detergent (0.01-0.1% Triton X-100)
- Aggregators lose activity with detergent

**Dynamic light scattering (DLS):**
- Measure particle size in solution
- Aggregators show particles > 100 nm

**Turbidity:**
- Visual inspection or spectrophotometry (320-340 nm)
- Aggregators cause cloudiness

**Concentration dependence:**
- Test multiple concentrations
- Aggregators show non-linear behavior

## False Positives vs True Hits

| Observation | Aggregator | True Hit |
|-------------|------------|----------|
| **Detergent effect** | Activity lost | Activity maintained |
| **Hill slope** | > 1 (steep) | ~1 (normal) |
| **Light scattering** | Positive | Negative |
| **Concentration curve** | Non-linear | Linear/saturable |
| **SAR** | Flat (no clear SAR) | Clear SAR |

## Limitations

**Prediction accuracy:**
- Based on empirical patterns
- Cannot predict all aggregators
- Experimental validation required

**Context dependent:**
- Aggregation depends on assay conditions
- Buffer composition affects aggregation
- Temperature and concentration matter

**Not absolute:**
- Some high-risk molecules don't aggregate
- Some low-risk molecules do aggregate
- Use as a guide, not absolute truth

:::warning Experimental Confirmation
Always validate aggregation experimentally for critical compounds. Computational prediction guides prioritization but doesn't replace experimental testing.
:::

## Best Practices

1. **Check all screening hits**: Especially those with steep dose-response
2. **Use detergent controls**: Include detergent in assay for high-risk compounds
3. **Orthogonal assays**: Confirm hits in different assay formats
4. **Document risk scores**: Track aggregator scores for future reference
5. **Filter libraries**: Remove high-risk aggregators before screening

## Reference

Irwin, J. J. et al. (2015). An aggregation advisor for ligand discovery. *Journal of Medicinal Chemistry*, 58(17), 7076–7087.

## Next Steps

- **[Safety Filters](/docs/user-guide/scoring/safety-filters)** — Structural alert screening
- **[Drug-likeness](/docs/user-guide/scoring/drug-likeness)** — Drug-likeness filters
- **[Scoring Overview](/docs/user-guide/scoring/overview)** — All scoring systems
