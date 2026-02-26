---
sidebar_position: 2
title: ML Readiness
description: Evaluate molecular suitability for machine learning applications
---

# ML-Readiness Scoring

The ML-readiness score evaluates how suitable a molecule is for machine learning workflows by assessing structural quality, physicochemical properties, complexity, and representation quality across four dimensions totaling 100 points.

## Score Dimensions

| Dimension | Max Points | What It Measures |
|-----------|-----------|------------------|
| **Structural Quality** | 20 | Structural soundness for ML pipelines |
| **Property Profile** | 35 | Physicochemical property desirability |
| **Complexity & Feasibility** | 25 | Synthetic tractability and complexity |
| **Representation Quality** | 20 | Numerical representability for ML models |

## 1. Structural Quality (20 points)

Binary pass/fail checks on fundamental structural soundness. Each item either passes (full points) or fails (0 points).

| Item | Points | Pass Condition |
|------|--------|----------------|
| Single component | 5 | Exactly 1 fragment — no mixtures, salts, or solvents |
| Standard organic elements | 5 | No metal atoms present (no organometallics) |
| No radicals | 3 | No atoms with unpaired electrons |
| Reasonable charge | 3 | Net formal charge between −2 and +2 |
| No dummy atoms | 4 | No R-groups or attachment points (atomic number ≠ 0) |
| **Total** | **20** | Sum of passed items |

**Non-scored caveats** (reported as warnings):
- Isotope labels detected (deuterium, ¹³C, tritium)
- Trivial molecules (≤ 3 heavy atoms)

## 2. Property Profile (35 points)

Desirability-scored physicochemical properties measuring how well a molecule fits within typical ML training set distributions. Each property is scored using a trapezoidal desirability function.

| Property | Ideal Range | Max Points |
|----------|------------|------------|
| Molecular Weight | 200–500 Da | 6 |
| LogP (Wildman-Crippen) | 0.5–4.0 | 6 |
| TPSA | 40–120 A² | 5 |
| H-Bond Donors | 0–3 | 4 |
| H-Bond Acceptors | 2–8 | 4 |
| Rotatable Bonds | 1–8 | 5 |
| Aromatic Rings | 1–3 | 5 |
| **Total** | | **35** |

### Desirability Function

Each property is scored using a trapezoidal desirability function:

```
If min ≤ value ≤ max:  d = 1.0  (ideal range → full score)
If value < min:        d = max(0, 1.0 − (min − value) / range)
If value > max:        d = max(0, 1.0 − (value − max) / range)

where range = max − min
points = d × max_points
```

**Example:** A molecule with MW = 600 Da:
- range = 500 − 200 = 300
- d = max(0, 1.0 − (600 − 500) / 300) = 0.667
- points = 0.667 × 6 = 4.0

:::info Property Ranges
These ranges reflect the most common distributions in drug-like compound datasets used for ML training. Molecules outside these ranges aren't necessarily bad — they just fall outside the typical training distribution.
:::

## 3. Complexity & Feasibility (25 points)

Assesses synthetic tractability and structural complexity, which affect practical utility in ML-driven drug discovery campaigns.

| Component | Max Points | Calculation |
|-----------|-----------|-------------|
| QED | 8 | `QED.qed(mol) × 8` |
| SA Score | 8 | Inverse mapping (see below) |
| Fsp3 | 4 | `desirability(Fsp3, 0.2, 0.6) × 4` |
| Stereocenters | 5 | Complexity-adjusted scoring (see below) |
| **Total** | **25** | |

### SA Score → Points

Synthetic Accessibility Score (1–10) is inversely mapped to points:

| SA Score | Points | Interpretation |
|----------|--------|----------------|
| ≤ 3 | 8.0 | Easy to synthesize |
| 3–5 | 8.0 − (SA − 3) × 2.0 | Moderate complexity |
| 5–7 | 4.0 − (SA − 5) × 2.0 | Difficult synthesis |
| > 7 | 0.0 | Very difficult |

### Stereocenter Scoring

| Stereocenters | Base Score | Notes |
|---------------|-----------|-------|
| 0–4 | 5.0 | Manageable complexity |
| 5–8 | 5.0 − (n − 4) × 0.75 | Decreasing linearly |
| > 8 | 0.0 | Too complex |

**Penalty:** If more than 50% of stereocenters are undefined, the base score is halved (×0.5).

## 4. Representation Quality (20 points)

Measures how well the molecule can be numerically represented for ML models — the core requirement for any ML application.

| Component | Max Points | What It Tests |
|-----------|-----------|---------------|
| Descriptor completeness | 5 | Fraction of 451 RDKit descriptors computed successfully |
| Fingerprint generation | 5 | Weighted success across 7 fingerprint types |
| Fingerprint informativeness | 5 | Ideal bit density between 1% and 30% |
| Conformer generation | 5 | 3D coordinate generation via ETKDGv3 |
| **Total** | **20** | |

### Descriptors Tested (451 total)

| Descriptor Set | Count | Method |
|----------------|-------|--------|
| Standard RDKit | 217 | `Descriptors.CalcMolDescriptors()` |
| AUTOCORR2D | 192 | `rdMolDescriptors.CalcAUTOCORR2D()` |
| MQN | 42 | `rdMolDescriptors.MQNs_()` |

Score: `round(5.0 × (successful / 451), 2)`

### Fingerprint Types & Weights

| Fingerprint | Bits | Weight | Description |
|-------------|------|--------|-------------|
| Morgan (radius=2) | 2048 | 8 | Circular fingerprints (ECFP4-like) |
| Morgan Features | 2048 | 8 | Feature-based Morgan |
| MACCS Keys | 167 | 8 | 166-bit MACCS structural keys |
| Atom Pair | 2048 | 4 | Atom pair descriptors |
| Topological Torsion | 2048 | 4 | Topological torsion descriptors |
| RDKit FP | 2048 | 4 | Daylight-like path fingerprints |
| Avalon | 512 | 4 | Avalon toolkit fingerprints |
| **Total weight** | | **40** | |

Score: `round(5.0 × (sum of successful weights / 40), 2)`

### Fingerprint Informativeness

Measures whether fingerprints have useful information content (not too sparse, not too dense):

| Avg Bit Density | Score | Interpretation |
|-----------------|-------|----------------|
| 1–30% | 5.0 | Ideal information content |
| < 1% | Proportional | Too sparse (molecule too simple) |
| 30–45% | Decreasing | Too dense (losing discriminative power) |
| > 45% | 0.0 | Not informative |

### Conformer Generation

| Method | Points |
|--------|--------|
| ETKDGv3 success (seed=42, maxIter=500) | 5 |
| Random coordinate fallback | 3 |
| Complete failure | 0 |

## Overall Score & Tiers

```
Total = Structural Quality + Property Profile + Complexity & Feasibility + Representation Quality
```

| Score | Tier | Interpretation |
|-------|------|----------------|
| **85–100** | Excellent | Suitable for most ML workflows without modification |
| **70–84** | Good | Minor limitations; generally suitable with standard preprocessing |
| **50–69** | Moderate | Usable but may need careful feature selection or preprocessing |
| **30–49** | Limited | Significant challenges; consider alternatives or specialized models |
| **0–29** | Poor | Not recommended for standard ML pipelines |

The UI displays an interpretation banner with:
- Overall score badge
- Tier-specific guidance text
- Per-dimension health tags showing completion percentage

## API Usage

```bash
curl -X POST http://localhost:8001/api/v1/score \
  -H "Content-Type: application/json" \
  -d '{
    "molecule": "CC(=O)Oc1ccccc1C(=O)O",
    "include": ["ml_readiness"]
  }'
```

Response:
```json
{
  "ml_readiness": {
    "score": 88,
    "dimensions": [
      {
        "name": "Structural Quality",
        "score": 20.0,
        "max_score": 20,
        "items": [
          {"name": "Single component", "score": 5.0, "max_score": 5, "passed": true},
          {"name": "Standard organic elements", "score": 5.0, "max_score": 5, "passed": true},
          {"name": "No radicals", "score": 3.0, "max_score": 3, "passed": true},
          {"name": "Reasonable charge", "score": 3.0, "max_score": 3, "passed": true},
          {"name": "No dummy atoms", "score": 4.0, "max_score": 4, "passed": true}
        ]
      },
      {
        "name": "Property Profile",
        "score": 30.5,
        "max_score": 35,
        "items": [
          {"name": "MW", "score": 6.0, "max_score": 6, "passed": true},
          {"name": "LogP", "score": 6.0, "max_score": 6, "passed": true}
        ]
      },
      {
        "name": "Complexity & Feasibility",
        "score": 21.0,
        "max_score": 25
      },
      {
        "name": "Representation Quality",
        "score": 16.5,
        "max_score": 20
      }
    ],
    "interpretation": "Good ML candidate with minor limitations...",
    "caveats": []
  }
}
```

## Use Cases

### Dataset Curation for ML

Filter molecules before creating training datasets:

```
ml_readiness_score >= 80 AND validation_score >= 90
```

This ensures:
- Structurally sound for ML (no mixtures, metals, radicals)
- Properties within typical training distributions
- All required descriptors and fingerprints generate successfully
- 3D conformer can be generated

### Model Applicability Domain

Use ML-readiness to define applicability domain:

- Train models only on molecules with score >= 80
- Flag predictions on molecules with score < 80 as uncertain
- Exclude molecules with score < 60 from predictions

### Dimension-Specific Analysis

Use individual dimension scores to diagnose issues:

| Low Dimension | Likely Issue | Action |
|---------------|-------------|--------|
| Structural Quality | Mixtures, metals, radicals | Clean up structure first |
| Property Profile | Unusual MW, LogP, TPSA | May be outside typical drug-like space |
| Complexity & Feasibility | Hard to synthesize, many stereocenters | Consider if practical for your use case |
| Representation Quality | Descriptor failures, no 3D | May need specialized featurization |

## Limitations

**Does not test:**
- Descriptor quality or relevance to specific models
- Model-specific feature requirements
- Chemical space coverage of your training set
- Experimental measurability

**Assumes:**
- Standard RDKit descriptors are sufficient
- Common fingerprint types are appropriate
- Property ranges derived from typical drug-like datasets

:::tip Custom Requirements
ML-readiness tests standard descriptors and fingerprints. If your model uses custom features (e.g., graph neural network features), you'll need additional validation.
:::

## References

1. Bickerton, G. R. et al. (2012). Quantifying the chemical beauty of drugs. *Nature Chemistry*, 4(2), 90–98.
2. Ertl, P. & Schuffenhauer, A. (2009). Estimation of synthetic accessibility score. *Journal of Cheminformatics*, 1(1), 8.
3. Lovering, F. et al. (2009). Escape from flatland. *Journal of Medicinal Chemistry*, 52(21), 6752–6756.
4. Rogers, D. & Hahn, M. (2010). Extended-connectivity fingerprints. *Journal of Chemical Information and Modeling*, 50(5), 742–754.

## Next Steps

- **[Scoring Overview](/docs/user-guide/scoring/overview)** — All scoring systems
- **[Batch Processing](/docs/user-guide/batch-processing)** — Score large datasets
- **[API Reference](/docs/api/endpoints)** — Full scoring API
