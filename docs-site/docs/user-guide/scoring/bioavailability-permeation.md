---
sidebar_position: 9
title: Bioavailability & Permeation
description: Bioavailability Radar and BOILED-Egg classification for oral absorption and brain penetration
---

# Bioavailability & Permeation

ChemAudit evaluates a molecule's potential for oral bioavailability and membrane permeation using two complementary models: the **Bioavailability Radar** (6-axis physicochemical profile) and the **BOILED-Egg** classification (GI absorption and BBB penetration).

Both models appear in the **Scoring Profiles** tab under the "Bioavailability & Permeation" card after validating a molecule.

## Bioavailability Radar

The Bioavailability Radar is a 6-axis normalized assessment of oral bioavailability. Each axis represents a physicochemical property critical for a drug to be absorbed orally. If all six properties fall within their optimal ranges, the molecule has an excellent oral bioavailability profile.

### The Six Axes

| Axis | Property | Optimal Range | Unit | Why It Matters |
|------|----------|---------------|------|----------------|
| **LIPO** | WLOGP (Lipophilicity) | -0.7 to 5.0 | — | Governs membrane permeability. Too hydrophilic = poor absorption; too lipophilic = poor solubility and toxicity risk. Upper bound from Lipinski's Rule of Five. |
| **SIZE** | Molecular Weight | 150 to 500 | Da | Larger molecules diffuse poorly across membranes. Upper bound of 500 Da from Lipinski's Rule of Five. |
| **POLAR** | TPSA (Topological Polar Surface Area) | 20 to 130 | A² | Reflects surface polarity. Values above ~140 A² correlate with poor intestinal absorption. Based on Veber's rules. |
| **INSOLU** | LogS (ESOL aqueous solubility) | -6 to 0 | log mol/L | A drug must dissolve in GI fluid before absorption. LogS below -6 indicates very poor solubility. |
| **INSATU** | Fraction Csp3 (sp3 carbon fraction) | ≥ 0.25 | ratio | Measures 3D character. Higher Fsp3 improves solubility and clinical success. From Lovering et al.'s "Escape from Flatland." |
| **FLEX** | Rotatable Bonds | ≤ 9 | count | Excessive flexibility imposes entropic penalties on binding and reduces oral bioavailability. Based on Veber's rules. |

### Interpretation

| In-Range Count | Interpretation |
|----------------|----------------|
| 6/6 | Excellent bioavailability profile |
| 4–5/6 | Good profile |
| 2–3/6 | Moderate profile — some properties outside optimal ranges |
| 0–1/6 | Poor profile — significant bioavailability concerns |

### Normalization

Each axis is normalized to a 0–1 scale:
- **Within optimal range**: normalized value = 1.0
- **Outside range**: linearly decreases toward 0 based on distance from the range boundary, reaching 0 at a distance equal to the range width

### Implementation Note

ChemAudit uses **WLOGP** (Wildman-Crippen LogP via RDKit's `Crippen.MolLogP`) for the lipophilicity axis. The original SwissADME paper uses XLOGP3 for the radar. Both are well-validated lipophilicity descriptors with closely correlated values, and the same optimal ranges apply.

## BOILED-Egg Classification

**BOILED-Egg** stands for **B**rain **O**r **I**ntestina**L** **E**stimate**D** permeation. It is a graphical model that simultaneously predicts passive gastrointestinal (GI) absorption and blood-brain barrier (BBB) penetration using only two molecular descriptors:

- **X-axis:** TPSA (Topological Polar Surface Area)
- **Y-axis:** WLOGP (Wildman-Crippen LogP)

### Classification Regions

The model defines three zones on a TPSA vs WLOGP scatter plot:

| Region | Color | Meaning | Approximate Boundaries |
|--------|-------|---------|----------------------|
| **White** | Green | High probability of passive GI absorption | TPSA < ~142 A², WLOGP between ~-2.3 and ~+6.8 |
| **Yolk** | Amber | High probability of BBB permeation (also GI absorbed) | TPSA < ~79 A², WLOGP between ~+0.4 and ~+6.0 |
| **Grey** | Grey | Predicted as neither passively GI-absorbed nor BBB-permeant | Outside both ellipses |

A molecule in the **yolk** is always also in the **white** — BBB-permeant compounds are predicted to be GI-absorbed as well, consistent with biological reality.

### Elliptical Model

Each region is defined by an ellipse in the TPSA-WLOGP plane. A molecule's (TPSA, WLOGP) point is tested against each ellipse using the standard point-in-ellipse equation:

```
((TPSA - cx) / a)² + ((WLOGP - cy) / b)² ≤ 1
```

Where `(cx, cy)` is the ellipse center and `(a, b)` are the semi-axes. The ellipse parameters were optimized using Monte-Carlo optimization evaluated by the Matthews Correlation Coefficient (MCC).

### Accuracy

| Metric | GI Absorption (White) | BBB Permeation (Yolk) |
|--------|----------------------|----------------------|
| Internal Accuracy | 93% | 90% |
| Internal MCC | 0.70 | 0.79 |
| 10-fold CV Accuracy | 92% | 88% |
| 10-fold CV MCC | 0.65 | 0.75 |
| External Validation (46 FDA NCEs) | 83% | — |

### Training Data

- **GI absorption dataset:** 660 molecules (567 well-absorbed, 93 poorly absorbed)
- **BBB permeation dataset:** 260 molecules (156 BBB-permeant, 104 non-permeant)

## References

1. **Bioavailability Radar / SwissADME:** Daina A, Michielin O, Zoete V. SwissADME: a free web tool to evaluate pharmacokinetics, drug-likeness and medicinal chemistry friendliness of small molecules. *Sci. Rep.* 2017, **7**, 42717. [DOI: 10.1038/srep42717](https://doi.org/10.1038/srep42717)

2. **BOILED-Egg:** Daina A, Zoete V. A BOILED-Egg to predict gastrointestinal absorption and brain penetration of small molecules. *ChemMedChem* 2016, **11**(11), 1117–1121. [DOI: 10.1002/cmdc.201600182](https://doi.org/10.1002/cmdc.201600182)

3. **Wildman-Crippen LogP (WLOGP):** Wildman SA, Crippen GM. Prediction of physicochemical parameters by atomic contributions. *J. Chem. Inf. Comput. Sci.* 1999, **39**(5), 868–873. [DOI: 10.1021/ci990307l](https://doi.org/10.1021/ci990307l)

4. **Lipinski's Rule of Five:** Lipinski CA, Lombardo F, Dominy BW, Feeney PJ. Experimental and computational approaches to estimate solubility and permeability in drug discovery and development settings. *Adv. Drug Deliv. Rev.* 2001, **46**(1-3), 3–26. [DOI: 10.1016/S0169-409X(00)00129-0](https://doi.org/10.1016/S0169-409X(00)00129-0)

5. **Veber Rules:** Veber DF, Johnson SR, Cheng HY, Smith BR, Ward KW, Kopple KD. Molecular properties that influence the oral bioavailability of drug candidates. *J. Med. Chem.* 2002, **45**(12), 2615–2623. [DOI: 10.1021/jm020017n](https://doi.org/10.1021/jm020017n)

6. **Escape from Flatland (Fsp3):** Lovering F, Bikker J, Humblet C. Escape from flatland: increasing saturation as an approach to improving clinical success. *J. Med. Chem.* 2009, **52**(21), 6752–6756. [DOI: 10.1021/jm901241e](https://doi.org/10.1021/jm901241e)

7. **ESOL Solubility:** Delaney JS. ESOL: estimating aqueous solubility directly from molecular structure. *J. Chem. Inf. Comput. Sci.* 2004, **44**(3), 1000–1005. [DOI: 10.1021/ci034243x](https://doi.org/10.1021/ci034243x)

## Next Steps

- **[ADMET](/docs/user-guide/scoring/admet)** — SA Score, ESOL, CNS MPO, and pharmaceutical rules
- **[Drug-likeness](/docs/user-guide/scoring/drug-likeness)** — Lipinski, QED, Veber, and consensus scoring
- **[Scoring Overview](/docs/user-guide/scoring/overview)** — All scoring systems at a glance
- **[Scoring Profiles](/docs/user-guide/scoring/profiles)** — Custom thresholds and desirability scoring
