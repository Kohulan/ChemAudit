<div align="center">

<img src="assets/logo.png" alt="ChemAudit" width="80" />

# Scoring Methodology

### Technical Reference for All Calculations, Formulas, and Thresholds

</div>

---

## Table of Contents

- [Overview](#overview)
- [Validation Checks](#validation-checks)
  - [Basic Checks](#basic-checks)
  - [Deep Validation](#deep-validation)
    - [Chemical Composition](#chemical-composition)
    - [Structural Complexity](#structural-complexity)
    - [Stereo & Tautomers](#stereo--tautomers)
- [ML-Readiness Scoring](#ml-readiness-scoring)
  - [Structural Quality](#1-structural-quality-20-points)
  - [Property Profile](#2-property-profile-35-points)
  - [Complexity & Feasibility](#3-complexity--feasibility-25-points)
  - [Representation Quality](#4-representation-quality-20-points)
  - [Overall Score & Tiers](#overall-score--tiers)
- [Drug-Likeness](#drug-likeness)
  - [QED](#qed-quantitative-estimate-of-drug-likeness)
  - [Lipinski Rule of Five](#lipinski-rule-of-five)
  - [Veber Rules](#veber-rules)
  - [Rule of Three](#rule-of-three)
  - [Ghose Filter](#ghose-filter)
  - [Egan Filter](#egan-filter)
  - [Muegge Filter](#muegge-filter)
  - [Consensus Drug-Likeness](#consensus-drug-likeness)
  - [Lead-Likeness](#lead-likeness)
- [ADMET Predictions](#admet-predictions)
  - [Synthetic Accessibility](#synthetic-accessibility-sa-score)
  - [Aqueous Solubility (ESOL)](#aqueous-solubility-esol)
  - [Molecular Complexity](#molecular-complexity)
  - [CNS MPO Score](#cns-mpo-score)
  - [Bioavailability Indicators](#bioavailability-indicators)
  - [Pfizer 3/75 Rule](#pfizer-375-rule)
  - [GSK 4/400 Rule](#gsk-4400-rule)
  - [Golden Triangle](#golden-triangle)
- [NP-Likeness](#np-likeness)
- [Aggregator Likelihood](#aggregator-likelihood)
- [Safety Filters](#safety-filters)
- [Standardization Pipeline](#standardization-pipeline)
- [References](#references)

---

## Overview

ChemAudit evaluates chemical structures through a multi-layered pipeline:

1. **Validation** — Structural integrity checks (basic + deep)
2. **Scoring** — Quantitative assessment across 7 dimensions (ML-Readiness, Drug-Likeness, ADMET, NP-Likeness, Aggregator, Safety, Scaffold)
3. **Standardization** — Structure normalization using a ChEMBL-compatible pipeline
4. **Database Lookup** — Cross-referencing with PubChem, ChEMBL, and COCONUT

All molecular property calculations use RDKit (2025.9.5). Scoring modules are deterministic — the same input always produces the same output.

---

## Validation Checks

### Basic Checks

These checks run on every molecule and assess fundamental structural validity.

| ID | Check | Severity | What It Tests |
|----|-------|----------|---------------|
| BASIC-01 | **Parsability** | Critical | Can RDKit parse the input into a valid molecule object? |
| BASIC-02 | **Sanitization** | Error | Does the molecule pass RDKit sanitization (valence correction, aromaticity perception, kekulization)? |
| BASIC-03 | **Valence** | Critical | Do all atoms have chemically valid bond counts? Uses `Chem.DetectChemistryProblems()`. |
| BASIC-04 | **Aromaticity** | Error | Can aromatic systems be kekulized (assigned explicit single/double bonds)? |
| BASIC-05 | **Connectivity** | Warning | Is the molecule a single connected component? Multiple fragments indicate mixtures or salts. |

**Scoring:** Each check contributes to the overall validation score (0–100). Critical failures result in a score of 0; warnings reduce the score proportionally.

### Deep Validation

Deep validation runs 17 specialized checks organized into three domains.

#### Chemical Composition

Checks that examine what the molecule is made of.

| ID | Check | Default Severity | Logic | Details |
|----|-------|-----------------|-------|---------|
| DVAL-06 | **Mixture Detection** | Warning | Counts disconnected fragments via `Chem.GetMolFrags()`. Each fragment classified as drug, salt, solvent, or unknown using MolVS fragment patterns. | Reports fragment SMILES, MW, heavy atom count, classification |
| DVAL-07 | **Solvent Contamination** | Warning | Matches fragments against 15+ known solvents (water, DMSO, DMF, acetonitrile, methanol, ethanol, acetone, THF, DCM, chloroform, toluene, hexane, etc.) via canonical SMILES match + substructure match + MolVS patterns. | Lists detected solvent names |
| DVAL-08 | **Inorganic Filter** | Warning/Error | Checks for carbon atoms and scans for metals. Error if no carbon (inorganic); Warning if carbon + metal (organometallic). | Reports element counts, is_inorganic, is_organometallic flags |
| DVAL-09 | **Radical Detection** | Warning | Checks `atom.GetNumRadicalElectrons() > 0` for any atom. | Lists radical atoms with electron counts |
| DVAL-10 | **Isotope Labels** | Info | Checks `atom.GetIsotope() > 0` for any atom. Maps to common names (deuterium, carbon-13, tritium, etc.). | Lists isotope-labeled atoms |
| DVAL-11 | **Trivial Molecule** | Error | Flags molecules with ≤ 3 heavy atoms as too small for meaningful analysis. | Reports heavy atom count |

#### Structural Complexity

Checks that assess structural features that may complicate analysis.

| ID | Check | Default Severity | Logic | Details |
|----|-------|-----------------|-------|---------|
| DVAL-12 | **Hypervalent Atoms** | Warning | Compares actual valence against maximum allowed valence per element. Transition metals with unrestricted valence are skipped. | Lists hypervalent atoms |
| DVAL-13 | **Polymer Detection** | Info | Three heuristics: (1) SGroup markers (SRU/COP type), (2) MW > 1500 Da, (3) dummy atoms (Z = 0). | Reports detection method |
| DVAL-14 | **Ring Strain** | Warning | Flags 3- or 4-membered rings, which carry significant angle strain. Uses `ring_info.AtomRings()`. | Lists strained ring sizes and atom indices |
| DVAL-15 | **Macrocycle Detection** | Info | Flags rings with > 12 atoms. Uses SSSR (Smallest Set of Smallest Rings). | Reports ring sizes |
| DVAL-16 | **Charged Species** | Info | Calculates net charge, identifies positive/negative atoms, detects zwitterions (`net_charge == 0` with both positive and negative atoms). | Reports net charge, charged atom positions, zwitterion flag |
| DVAL-17 | **Explicit Hydrogen Audit** | Info | Counts explicit H atoms and checks for mixed explicit/implicit H representation (a common RDKit pitfall). | Reports explicit H count, consistency |

#### Stereo & Tautomers

Checks related to stereochemistry and tautomeric forms.

| ID | Check | Default Severity | Logic | Details |
|----|-------|-----------------|-------|---------|
| DVAL-01 | **Stereoisomer Enumeration** | Warning | Enumerates possible stereoisomers from undefined centers using `EnumerateStereoisomers()`. Cap: 128 (2^7). | Count, SMILES list (if ≤ cap) |
| DVAL-02 | **Undefined Stereocenters** | Warning | Uses `Chem.FindMolChiralCenters(includeUnassigned=True)`. Counts centers marked `?`. | Undefined/total counts, atom indices |
| DVAL-03 | **Tautomer Detection** | Info | Uses `rdMolStandardize.TautomerEnumerator()` to enumerate tautomers and check if input matches the canonical form. | Count, canonical SMILES, is_canonical flag |
| DVAL-04 | **Aromatic System Validation** | Warning | Flags unusual aromatic ring sizes (not 5 or 6) and charged aromatic atoms. | Ring sizes, charged atom positions |
| DVAL-05 | **Coordinate Dimension** | Info | Inspects conformer data: 2D, 3D, no_coordinates, or degenerate. | Dimension classification |

**Severity Override:** Users can customize the severity level of any deep validation check through the severity configuration panel.

---

## ML-Readiness Scoring

Evaluates how suitable a molecule is for machine learning workflows. The score (0–100) is computed across four dimensions.

### 1. Structural Quality (20 points)

Binary pass/fail checks assessing structural soundness for ML pipelines.

| Item | Points | Pass Condition |
|------|--------|----------------|
| Single component | 5 | Exactly 1 fragment (no mixtures/salts) |
| Standard organic elements | 5 | No metal atoms present |
| No radicals | 3 | No unpaired electrons on any atom |
| Reasonable charge | 3 | Net formal charge ≤ ±2 |
| No dummy atoms | 4 | No R-groups or attachment points (Z ≠ 0) |
| **Total** | **20** | Sum of passed items |

**Caveats (non-scored, reported as warnings):**
- Isotope labels detected
- Trivial molecules (≤ 3 heavy atoms)

### 2. Property Profile (35 points)

Desirability-scored physicochemical properties measuring how well a molecule fits within typical ML training set distributions.

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

**Desirability function:**

Each property is scored using a trapezoidal desirability function:

```
If min ≤ value ≤ max:  d = 1.0  (ideal range, full score)
If value < min:        d = max(0, 1.0 − (min − value) / range)
If value > max:        d = max(0, 1.0 − (value − max) / range)

where range = max − min
```

The desirability `d` (0–1) is multiplied by the maximum points for that property:

```
points = round(d × max_points, 2)
```

### 3. Complexity & Feasibility (25 points)

Assesses synthetic tractability and structural complexity, which affect practical utility in ML campaigns.

| Component | Max Points | Calculation |
|-----------|-----------|-------------|
| QED | 8 | `QED.qed(mol) × 8` |
| SA Score | 8 | See mapping below |
| Fsp3 | 4 | `desirability(Fsp3, 0.2, 0.6) × 4` |
| Stereocenters | 5 | See mapping below |
| **Total** | **25** | |

**SA Score → Points mapping:**

```
SA ≤ 3:              8.0 points  (easy to synthesize)
3 < SA ≤ 5:          8.0 − (SA − 3) × 2.0  (linear interpolation)
5 < SA ≤ 7:          4.0 − (SA − 5) × 2.0  (linear interpolation)
SA > 7:              0.0 points  (very difficult)
```

**Stereocenter scoring:**

```
total_centers ≤ 4:      base = 5.0
4 < total_centers ≤ 8:  base = 5.0 − (total_centers − 4) × 0.75
total_centers > 8:      base = 0.0

If undefined / total > 0.5:  base × 0.5  (50% penalty for majority undefined)
```

### 4. Representation Quality (20 points)

Measures how well the molecule can be numerically represented for ML models.

| Component | Max Points | Logic |
|-----------|-----------|-------|
| Descriptor completeness | 5 | Fraction of 451 descriptors computed successfully × 5 |
| Fingerprint generation | 5 | Weighted success across 7 fingerprint types |
| Fingerprint informativeness | 5 | Ideal bit density 1–30% |
| Conformer generation | 5 | ETKDGv3 success = 5 pts; random fallback = 3 pts |
| **Total** | **20** | |

**Descriptors tested (451 total):**
- 217 standard RDKit descriptors (`Descriptors.CalcMolDescriptors()`)
- 192 AUTOCORR2D descriptors (`rdMolDescriptors.CalcAUTOCORR2D()`)
- 42 MQN descriptors (`rdMolDescriptors.MQNs_()`)

**Fingerprint types & weights:**

| Fingerprint | Bits | Weight |
|-------------|------|--------|
| Morgan (radius=2) | 2048 | 8 |
| Morgan Features | 2048 | 8 |
| MACCS Keys | 167 | 8 |
| Atom Pair | 2048 | 4 |
| Topological Torsion | 2048 | 4 |
| RDKit FP | 2048 | 4 |
| Avalon | 512 | 4 |
| **Total weight** | | **40** |

Score = `round(5.0 × (sum of successful weights / 40), 2)`

**Fingerprint informativeness:**

```
avg_density = mean(on_bits / total_bits) across all generated fingerprints

0.01 ≤ density ≤ 0.30:  5.0  (ideal range)
density < 0.01:          5.0 × (density / 0.01)  (too sparse)
0.30 < density ≤ 0.45:  5.0 × (1.0 − (density − 0.30) / 0.15)  (too dense)
density > 0.45:          0.0
```

**Conformer generation:**
- ETKDGv3 (`seed=42`, `maxIterations=500`): 5 points if successful
- Random coordinate fallback: 3 points if ETKDGv3 fails
- Complete failure: 0 points

### Overall Score & Tiers

```
Total Score = Structural Quality + Property Profile + Complexity & Feasibility + Representation Quality
```

| Score | Tier | Interpretation |
|-------|------|----------------|
| 85–100 | **Excellent** | Suitable for most ML workflows without modification |
| 70–84 | **Good** | Minor limitations; generally suitable with standard preprocessing |
| 50–69 | **Moderate** | Usable but may need careful feature selection or preprocessing |
| 30–49 | **Limited** | Significant challenges; consider alternatives or specialized models |
| 0–29 | **Poor** | Not recommended for standard ML pipelines |

---

## Drug-Likeness

Evaluates compliance with established pharmaceutical drug-likeness rules.

### QED (Quantitative Estimate of Drug-likeness)

A composite score (0–1) integrating eight molecular properties using a desirability function approach.

**Components:** MW, ALogP, HBA, HBD, PSA, rotatable bonds, aromatic rings, structural alerts

| QED Score | Interpretation |
|-----------|----------------|
| ≥ 0.67 | Favorable |
| 0.49–0.66 | Moderate |
| < 0.49 | Unfavorable |

**Reference:** Bickerton et al. (2012). Quantifying the chemical beauty of drugs. *Nature Chemistry*, 4(2), 90–98.

### Lipinski Rule of Five

The most widely used drug-likeness filter for oral bioavailability prediction.

| Property | Threshold |
|----------|-----------|
| Molecular Weight | ≤ 500 Da |
| LogP (Wildman-Crippen) | ≤ 5.0 |
| H-Bond Donors | ≤ 5 |
| H-Bond Acceptors | ≤ 10 |

**Pass criterion:** ≤ 1 violation allowed.

**Reference:** Lipinski et al. (2001). Experimental and computational approaches to estimate solubility and permeability in drug discovery and development settings. *Advanced Drug Delivery Reviews*, 46(1-3), 3–26.

### Veber Rules

Predicts oral bioavailability based on molecular flexibility and polar surface area.

| Property | Threshold |
|----------|-----------|
| Rotatable Bonds | ≤ 10 |
| TPSA | ≤ 140 A² |

**Pass criterion:** Both must pass.

**Reference:** Veber et al. (2002). Molecular properties that influence the oral bioavailability of drug candidates. *Journal of Medicinal Chemistry*, 45(12), 2615–2623.

### Rule of Three

Identifies fragment-like molecules suitable for fragment-based drug design.

| Property | Threshold |
|----------|-----------|
| Molecular Weight | < 300 Da |
| LogP | ≤ 3.0 |
| H-Bond Donors | ≤ 3 |
| H-Bond Acceptors | ≤ 3 |
| Rotatable Bonds | ≤ 3 |
| TPSA | ≤ 60 A² |

**Pass criterion:** All must pass.

**Reference:** Congreve et al. (2003). A 'rule of three' for fragment-based lead discovery? *Drug Discovery Today*, 8(19), 876–877.

### Ghose Filter

Drug-likeness criteria based on analysis of known drugs.

| Property | Accepted Range |
|----------|---------------|
| Molecular Weight | 160–480 Da |
| LogP | -0.4 to 5.6 |
| Atom Count | 20–70 |
| Molar Refractivity | 40–130 |

**Pass criterion:** All properties must fall within range.

**Reference:** Ghose et al. (1999). A knowledge-based approach in designing combinatorial or medicinal chemistry libraries for drug discovery. *Journal of Combinatorial Chemistry*, 1(1), 55–68.

### Egan Filter

Predicts intestinal absorption using a simple two-parameter model.

| Property | Threshold |
|----------|-----------|
| LogP | ≤ 5.88 |
| TPSA | ≤ 131.6 A² |

**Pass criterion:** Both must pass.

**Reference:** Egan et al. (2000). Prediction of drug absorption using multivariate statistics. *Journal of Medicinal Chemistry*, 43(21), 3867–3877.

### Muegge Filter

A comprehensive 9-parameter filter for druglike chemical space.

| Property | Threshold |
|----------|-----------|
| Molecular Weight | 200–600 Da |
| LogP | -2 to 5 |
| TPSA | ≤ 150 A² |
| Ring Count | ≤ 7 |
| Carbon Count | > 4 |
| Heteroatom Count | > 1 |
| Rotatable Bonds | ≤ 15 |
| H-Bond Donors | ≤ 5 |
| H-Bond Acceptors | ≤ 10 |

**Pass criterion:** All 9 parameters must pass (0 violations).

**Reference:** Muegge et al. (2001). Simple selection criteria for drug-like chemical matter. *Journal of Medicinal Chemistry*, 44(12), 1841–1846.

### Consensus Drug-Likeness

A 0–5 score counting how many of the following filter sets the molecule passes:

1. Lipinski (≤ 1 violation)
2. Veber (both parameters pass)
3. Egan (both parameters pass)
4. Ghose (all parameters in range)
5. Muegge (all 9 parameters pass)

| Score | Interpretation |
|-------|----------------|
| 5 | Excellent drug-likeness |
| 4 | Very good |
| 3 | Acceptable |
| 2 | Borderline |
| 0–1 | Poor drug-likeness |

### Lead-Likeness

Identifies molecules in lead-like chemical space, which provides room for optimization.

| Property | Threshold |
|----------|-----------|
| Molecular Weight | 200–350 Da |
| LogP | -1 to 3 |
| Rotatable Bonds | ≤ 7 |

**Pass criterion:** All must pass.

---

## ADMET Predictions

Calculated molecular properties predictive of ADMET (Absorption, Distribution, Metabolism, Excretion, Toxicity) behavior. All values are computed from molecular structure, not experimental measurements.

### Synthetic Accessibility (SA Score)

Estimates how difficult a molecule is to synthesize using fragment-based complexity analysis.

**Range:** 1–10 (lower = easier to synthesize)

| SA Score | Classification |
|----------|---------------|
| < 4 | Easy to synthesize |
| 4–6 | Moderate complexity |
| > 6 | Difficult to synthesize |

**Method:** Combines fragment contributions (how common molecular fragments are in known molecules) with a complexity penalty (stereocenters, macrocycles, spiro centers, ring fusions).

**Reference:** Ertl & Schuffenhauer (2009). Estimation of synthetic accessibility score of drug-like molecules based on molecular complexity and fragment contributions. *Journal of Cheminformatics*, 1(1), 8.

### Aqueous Solubility (ESOL)

Predicts aqueous solubility using the ESOL linear regression model.

**Equation:**

```
LogS = 0.16 − 0.63 × LogP − 0.0062 × MW + 0.066 × RotBonds − 0.74 × AP
```

Where:
- `LogS` = log₁₀(aqueous solubility in mol/L)
- `LogP` = Wildman-Crippen LogP
- `MW` = molecular weight (Da)
- `RotBonds` = rotatable bond count
- `AP` = aromatic proportion (aromatic atoms / heavy atoms)

**Conversion to mg/mL:**

```
solubility_mol_L = 10^LogS
solubility_mg_mL = solubility_mol_L × MW / 1000
```

**Classification:**

| LogS | Category | Approx. mg/mL |
|------|----------|---------------|
| ≥ -1 | Highly soluble | > 100 |
| -1 to -3 | Soluble | 1–100 |
| -3 to -4 | Moderately soluble | 0.1–1 |
| -4 to -5 | Poorly soluble | < 0.1 |
| < -5 | Insoluble | Very low |

**Reference:** Delaney (2004). ESOL: Estimating aqueous solubility directly from molecular structure. *Journal of Chemical Information and Computer Sciences*, 44(3), 1000–1005.

### Molecular Complexity

Multiple metrics characterizing 3D character and structural complexity.

| Metric | Calculation | Interpretation |
|--------|-------------|----------------|
| **Fsp3** | Fraction of sp³ carbons | Higher = more 3D character |
| **Stereocenters** | Count of chiral centers | Higher = more complex |
| **Ring Count** | Total rings (SSSR) | Complexity indicator |
| **Aromatic Rings** | Aromatic ring count | Flatness indicator |
| **Bertz Complexity** | Graph-based topological entropy | Overall complexity metric |

**Fsp3 classification:**

| Fsp3 | Character |
|------|-----------|
| < 0.25 | Flat/aromatic |
| 0.25–0.42 | Moderate 3D |
| > 0.42 | Good 3D character |

**Reference (Fsp3):** Lovering et al. (2009). Escape from flatland: increasing saturation as an approach to improving clinical success. *Journal of Medicinal Chemistry*, 52(21), 6752–6756.

### CNS MPO Score

Pfizer's Central Nervous System Multiparameter Optimization score (0–6 scale). Predicts likelihood of CNS penetration.

**Component scoring:**

| Parameter | Score = 1.0 | Linear decrease | Score = 0 |
|-----------|-------------|----------------|-----------|
| MW | ≤ 360 Da | 360–500 | > 500 |
| LogP | ≤ 3.0 | 3.0–5.0 | > 5.0 |
| TPSA | ≤ 40 A² | 40–90 | > 90 |
| HBD | 0 | 1–3 (−0.25 per donor) | > 3 |
| LogD | 0.5 (placeholder) | — | — |
| pKa | 0.5 (placeholder) | — | — |

```
CNS MPO = sum of all 6 component scores
```

| Score | Interpretation |
|-------|----------------|
| ≥ 5 | Excellent CNS penetration |
| 4–5 | Good CNS penetration |
| 3–4 | Moderate |
| < 3 | Poor CNS penetration |

**Reference:** Wager et al. (2010). Moving beyond rules: the development of a central nervous system multiparameter optimization (CNS MPO) approach to enable alignment of druglike properties. *ACS Chemical Neuroscience*, 1(6), 435–449.

### Bioavailability Indicators

Composite assessments combining multiple property thresholds.

**Oral absorption:**

```
oral_absorption_likely = lipinski_ok AND veber_ok
```

Where:
- `lipinski_ok` = MW ≤ 500, LogP ≤ 5, HBD ≤ 5, HBA ≤ 10
- `veber_ok` = rotatable bonds ≤ 10, TPSA ≤ 140

**CNS penetration:**

```
cns_penetration_likely = TPSA ≤ 90 AND MW ≤ 450 AND HBD ≤ 3 AND 1 ≤ LogP ≤ 4
```

### Pfizer 3/75 Rule

A simple toxicity risk indicator based on LogP and TPSA.

```
at_risk = LogP > 3 AND TPSA < 75
```

Compounds meeting both criteria have statistically higher rates of toxicity and off-target promiscuity.

**Reference:** Hughes et al. (2008). Physiochemical drug properties associated with in vivo toxicological outcomes. *Bioorganic & Medicinal Chemistry Letters*, 18(17), 4872–4875.

### GSK 4/400 Rule

GlaxoSmithKline's compound quality guideline for favorable ADMET properties.

```
favorable = MW ≤ 400 AND LogP ≤ 4
```

**Reference:** Gleeson (2008). Generation of a set of simple, interpretable ADMET rules of thumb. *Journal of Medicinal Chemistry*, 51(4), 817–834.

### Golden Triangle

Abbott's optimal property space for balanced permeability and metabolic stability.

```
in_triangle = 200 ≤ MW ≤ 450 AND −0.5 ≤ LogP ≤ 5
```

Compounds within this space tend to have favorable permeability-metabolism balance.

**Reference:** Johnson et al. (2009). Using the Golden Triangle to optimize clearance and oral absorption. *Bioorganic & Medicinal Chemistry Letters*, 19(17), 5560–5564.

---

## NP-Likeness

Scores how closely a molecule resembles known natural products using fragment-based analysis.

**Range:** −5 to +5 (typical: −3 to +3)
- Positive scores → natural product-like
- Negative scores → synthetic-like
- Near zero → mixed character

**Primary method:** RDKit's NP-likeness model based on Morgan fingerprint fragment contributions from natural product training sets.

**Heuristic fallback** (when RDKit NP model is unavailable):

| Feature | Condition | Score |
|---------|-----------|-------|
| Ring fraction | 0.2–0.6 | +0.5 |
| Ring fraction | > 0.6 | −0.3 |
| O/C ratio | 0.1–0.5 | +0.5 |
| O/C ratio | > 0.5 | +0.3 |
| N/C ratio | > 0.3 | −0.3 |
| Halogens | Per atom (max 3) | −0.3 each |
| Fsp3 | > 0.4 | +0.5 |
| Fsp3 | 0.25–0.4 | +0.3 |
| Fsp3 | ≤ 0.25 | −0.2 |
| Chiral centers | ≥ 3 | +0.5 |
| Chiral centers | 1–2 | +0.2 |
| MW | 200–800 Da | +0.3 |
| MW | < 150 Da | −0.2 |

Final heuristic scaling: `score = max(−3.0, min(3.0, base × 1.5))`

**Interpretation:**

| Score | Category |
|-------|----------|
| ≥ 2.0 | Strong NP-like |
| 1.0–2.0 | NP-like |
| 0.3–1.0 | Moderate NP-like |
| −0.3 to 0.3 | Mixed character |
| −1.0 to −0.3 | Moderate synthetic |
| −2.0 to −1.0 | Synthetic-like |
| < −2.0 | Strong synthetic |

**Reference:** Ertl et al. (2008). Natural product-likeness score and its application for prioritization of compound libraries. *Journal of Chemical Information and Modeling*, 48(1), 68–74.

---

## Aggregator Likelihood

Predicts whether a molecule is likely to form colloidal aggregates in biological assays, which cause false-positive readouts.

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
| Known scaffold match | SMARTS pattern match | +0.20 per pattern |

### Known Aggregator Scaffolds

The following SMARTS patterns are matched against the molecule:

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

### Scoring

```
risk_score = min(1.0, sum of all triggered increments)
confidence = triggered_indicators / total_indicators
```

| Risk Score | Classification |
|------------|---------------|
| ≥ 0.6 | High aggregation risk |
| 0.3–0.59 | Moderate risk |
| < 0.3 | Low risk |

**Reference:** Irwin et al. (2015). An aggregation advisor for ligand discovery. *Journal of Medicinal Chemistry*, 58(17), 7076–7087.

---

## Safety Filters

Screens molecules against curated catalogs of substructural patterns known to cause assay interference, toxicity, or poor pharmacokinetic properties.

### Available Catalogs

| Catalog | Source | Approx. Patterns | Purpose |
|---------|--------|-------------------|---------|
| **PAINS A/B/C** | Baell & Holloway (2010) | ~480 | Pan-Assay Interference Compounds — frequent hitters in HTS |
| **Brenk** | Brenk et al. (2008) | ~105 | Unfavorable chemical moieties for drug development |
| **NIH** | NIH MLSMR Program | ~180 | Molecular Libraries screening exclusion filters |
| **ZINC** | ZINC Database | ~95 | Drug-likeness and reactivity filters |
| **ChEMBL** | Pharma companies | ~700+ | Combined from 7 sub-catalogs (see below) |

### ChEMBL Sub-Catalogs

| Sub-Catalog | Source | Focus |
|-------------|--------|-------|
| BMS | Bristol-Myers Squibb | Reactive functional groups |
| Dundee | University of Dundee | Promiscuous compound filters |
| Glaxo | GlaxoSmithKline | Undesirable moieties |
| Inpharmatica | Inpharmatica Ltd. | Chemical liabilities |
| LINT | — | Lead identification noise filters |
| MLSMR | NIH/MLSMR | Molecular Libraries screening filters |
| SureChEMBL | EMBL-EBI | Patent literature structural alerts |

### Matching Logic

All pattern matching uses RDKit's `FilterCatalog` module, which performs SMARTS substructure matching against the molecule.

**Important context:** 87 FDA-approved drugs trigger PAINS patterns. Alerts are warnings for investigation, not automatic rejection criteria.

**References:**
- Baell & Holloway (2010). New substructure filters for removal of pan assay interference compounds (PAINS). *Journal of Medicinal Chemistry*, 53(7), 2719–2740.
- Brenk et al. (2008). Lessons learnt from assembling screening libraries for drug discovery for neglected diseases. *ChemMedChem*, 3(3), 435–444.

---

## Standardization Pipeline

Normalizes molecular representations using a ChEMBL-compatible pipeline powered by the `chembl_structure_pipeline` and MolVS libraries.

### Pipeline Steps

| Step | Description | Always Runs |
|------|-------------|-------------|
| **Checker** | Detects structural issues (valence errors, stereo problems, size limits) before processing | Yes |
| **Standardizer** | Normalizes functional groups (e.g., nitro groups, metal disconnection, charge neutralization) | Yes |
| **Get Parent** | Extracts the parent molecule by removing salts, solvents, and counterions | Yes |
| **Tautomer** | Canonicalizes tautomeric forms using RDKit's tautomer enumerator | No (opt-in) |

### Options

| Option | Default | Description |
|--------|---------|-------------|
| `include_tautomer` | `false` | Enable tautomer canonicalization. May alter E/Z stereochemistry. |
| `preserve_stereo` | `true` | Attempt to preserve stereochemistry during standardization. |

### Output

The standardization response includes:
- Original and standardized SMILES
- Steps applied with per-step change descriptions
- Checker issues found (pre-standardization)
- Excluded fragments (salts, solvents with classification)
- Stereochemistry comparison (if stereo changed)
- Structure comparison (atom count, formula, mass change percentage)

---

## References

1. Baell, J. B. & Holloway, G. A. (2010). New substructure filters for removal of pan assay interference compounds (PAINS) from screening libraries and for their exclusion in bioassays. *Journal of Medicinal Chemistry*, 53(7), 2719–2740.

2. Bickerton, G. R., Paolini, G. V., Besnard, J., Muresan, S. & Hopkins, A. L. (2012). Quantifying the chemical beauty of drugs. *Nature Chemistry*, 4(2), 90–98.

3. Brenk, R. et al. (2008). Lessons learnt from assembling screening libraries for drug discovery for neglected diseases. *ChemMedChem*, 3(3), 435–444.

4. Congreve, M., Carr, R., Murray, C. & Jhoti, H. (2003). A 'rule of three' for fragment-based lead discovery? *Drug Discovery Today*, 8(19), 876–877.

5. Delaney, J. S. (2004). ESOL: Estimating aqueous solubility directly from molecular structure. *Journal of Chemical Information and Computer Sciences*, 44(3), 1000–1005.

6. Egan, W. J., Merz, K. M. & Baldwin, J. J. (2000). Prediction of drug absorption using multivariate statistics. *Journal of Medicinal Chemistry*, 43(21), 3867–3877.

7. Ertl, P. & Schuffenhauer, A. (2009). Estimation of synthetic accessibility score of drug-like molecules based on molecular complexity and fragment contributions. *Journal of Cheminformatics*, 1(1), 8.

8. Ertl, P., Roggo, S. & Schuffenhauer, A. (2008). Natural product-likeness score and its application for prioritization of compound libraries. *Journal of Chemical Information and Modeling*, 48(1), 68–74.

9. Ghose, A. K., Viswanadhan, V. N. & Wendoloski, J. J. (1999). A knowledge-based approach in designing combinatorial or medicinal chemistry libraries for drug discovery. *Journal of Combinatorial Chemistry*, 1(1), 55–68.

10. Gleeson, M. P. (2008). Generation of a set of simple, interpretable ADMET rules of thumb. *Journal of Medicinal Chemistry*, 51(4), 817–834.

11. Hughes, J. D. et al. (2008). Physiochemical drug properties associated with in vivo toxicological outcomes. *Bioorganic & Medicinal Chemistry Letters*, 18(17), 4872–4875.

12. Irwin, J. J. et al. (2015). An aggregation advisor for ligand discovery. *Journal of Medicinal Chemistry*, 58(17), 7076–7087.

13. Johnson, T. W. et al. (2009). Using the Golden Triangle to optimize clearance and oral absorption. *Bioorganic & Medicinal Chemistry Letters*, 19(17), 5560–5564.

14. Lipinski, C. A. et al. (2001). Experimental and computational approaches to estimate solubility and permeability in drug discovery and development settings. *Advanced Drug Delivery Reviews*, 46(1-3), 3–26.

15. Lovering, F., Bikker, J. & Humblet, C. (2009). Escape from flatland: increasing saturation as an approach to improving clinical success. *Journal of Medicinal Chemistry*, 52(21), 6752–6756.

16. Muegge, I., Heald, S. L. & Brittelli, D. (2001). Simple selection criteria for drug-like chemical matter. *Journal of Medicinal Chemistry*, 44(12), 1841–1846.

17. Veber, D. F. et al. (2002). Molecular properties that influence the oral bioavailability of drug candidates. *Journal of Medicinal Chemistry*, 45(12), 2615–2623.

18. Wager, T. T. et al. (2010). Moving beyond rules: the development of a central nervous system multiparameter optimization (CNS MPO) approach to enable alignment of druglike properties. *ACS Chemical Neuroscience*, 1(6), 435–449.

---

<div align="center">

**Need more help?** Check the [User Guide](USER_GUIDE.md) or [API Reference](API_REFERENCE.md)

</div>
