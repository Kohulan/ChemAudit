# ML Readiness Scorer Redesign

**Date**: 2026-02-25
**Status**: Approved
**Branch**: features_v1

## Problem

The current ML Readiness score is a computation health check — "can RDKit compute descriptors?" For any valid molecule, it's almost always 100/100, making it scientifically meaningless.

## Solution

Redesign as an orchestrator that calls existing scoring/validation services and synthesizes results into a scientifically meaningful 0-100 score across 4 dimensions.

## Score Architecture: 0-100, Four Dimensions

### Dimension 1: Structural Quality (20 pts)

Assesses whether the molecule is structurally clean for ML pipelines.

| Check | Points | Existing Code |
|-------|--------|---------------|
| Single connected component | 5 | `MixtureDetectionCheck` (DVAL-06) |
| Standard organic elements (C,N,O,S,P,F,Cl,Br,I,H) | 5 | `InorganicFilterCheck` (DVAL-08) + `METAL_ATOMIC_NUMS` |
| No radicals | 3 | `RadicalDetectionCheck` (DVAL-09) |
| Reasonable formal charge (\|net\| <= 2) | 3 | `ChargedSpeciesCheck` |
| No dummy atoms (not polymer/fragment) | 4 | `PolymerDetectionCheck` |

Binary per check — pass gets full points, fail gets 0. New code: zero, only a thin wrapper.

### Dimension 2: Property Profile (35 pts)

Physicochemical properties scored against QED-derived ideal ranges (Bickerton et al. 2012).

| Property | Ideal Range | Max Pts |
|----------|-------------|---------|
| Molecular Weight | 200-500 Da | 6 |
| LogP (Crippen) | 0.5-4.0 | 6 |
| TPSA | 40-120 A² | 5 |
| HBD | 0-3 | 4 |
| HBA | 2-8 | 4 |
| Rotatable Bonds | 1-8 | 5 |
| Aromatic Rings | 1-3 | 5 |

**Scoring**: Uses `profile_scoring.desirability(value, min_val, max_val)` (already implemented) for smooth 0-1 mapping. Score = desirability * max_points.

**Supplementary labels** (not scored): Lipinski violations (0-4), Veber pass/fail — from `DrugLikenessScorer`.

**Existing code**: All property values already computed by `druglikeness.py` and `bioavailability_radar.py`. Only need range definitions + desirability calls.

**Citation**: Bickerton et al. (2012) "Quantifying the chemical beauty of drugs." Nature Chemistry 4:90-98. Lipinski et al. (1997) Adv. Drug Deliv. Rev. 23:3-25. Veber et al. (2002) J. Med. Chem. 45:2615-2623.

### Dimension 3: Complexity & Feasibility (25 pts)

| Metric | Max Pts | Existing Code |
|--------|---------|---------------|
| QED (0-1) | 8 | `druglikeness.py` → `_calculate_qed()` |
| SA Score (1-10) | 8 | `admet.py` → `_calculate_synthetic_accessibility()` |
| Fsp3 | 4 | `admet.py` → `_calculate_complexity()` |
| Stereocenter manageability | 5 | `StereoisomerEnumerationCheck` + `UndefinedStereoCentersCheck` |

**Scoring logic**:
- QED: `qed_value * 8` (linear)
- SA Score: inversely mapped via desirability. SA 1-3=8pts, 3-5=6pts, 5-7=3pts, >7=0pts (smooth ramp)
- Fsp3: desirability(fsp3, 0.2, 0.6) * 4
- Stereocenters: 0-4=5pts, 5-8=3pts (linear), >8=0pts. If >50% undefined, apply 0.5x penalty.

**Citation**: QED: Bickerton et al. (2012) Nature Chemistry 4:90-98. SA Score: Ertl & Schuffenhauer (2009) J. Cheminf. 1:8. Fsp3: Lovering et al. (2009) J. Med. Chem. 52:6752-6756.

### Dimension 4: Representation Quality (20 pts)

| Check | Max Pts | Source |
|-------|---------|--------|
| Descriptor completeness (217 RDKit) | 5 | Current `_score_descriptors()` |
| Fingerprint generation (7 types) | 5 | Current `_score_fingerprints()` |
| Fingerprint informativeness (bit density) | 5 | **New** — `fp.GetNumOnBits()/fp.GetNumBits()` |
| Conformer generation (3D feasibility) | 5 | **New** — ETKDGv3 + randomCoords fallback |

**Fingerprint informativeness**: Average bit density across FP types. Ideal 1-30% → 5pts. Too sparse (<1%) or saturated (>45%) → reduced.

**Conformer generation**: ETKDGv3 with `useRandomCoords=True` fallback, 30s timeout. Success=5pts, failure=0pts.

**Citation**: Rogers & Hahn (2010) J. Chem. Inf. Model. 50:742-754. Riniker & Landrum (2015) J. Chem. Inf. Model. 55:2562-2574.

## Optional ChEMBL Context (pill toggle, not scored)

When enabled, live ChEMBL API lookup:
- Known compound? → ChEMBL ID + preferred name
- Clinical phase (max_phase 0-4)
- Bioactivity count + distinct targets
- Data richness: Rich (>100) / Moderate (10-100) / Sparse (1-10) / No data

Purely informational, displayed as supplementary panel.

## Interpretation Thresholds

| Score | Label | Color |
|-------|-------|-------|
| 85-100 | Excellent | Green |
| 70-84 | Good | Blue/Teal |
| 50-69 | Moderate | Amber |
| 30-49 | Limited | Orange |
| 0-29 | Poor | Red |

## Structural Caveats (always displayed when relevant)

Pulled from existing validation checks — no new detection logic:
- Multi-component → `MixtureDetectionCheck`
- Organometallic/Inorganic → `InorganicFilterCheck`
- Radicals → `RadicalDetectionCheck`
- Isotope labels → `IsotopeLabelDetectionCheck`
- Undefined stereocenters → `StereoisomerEnumerationCheck`
- Trivially small → `TrivialMoleculeCheck`

## Frontend Tooltips

Every dimension card gets an InfoTooltip with:
1. **What** it measures (plain English)
2. **How** it's calculated (formula/method)
3. **Why** it matters for ML (practical consequence)
4. **Citation** (author, year, journal)

## Files to Modify

**Backend** (rewrite):
- `backend/app/services/scoring/ml_readiness.py`
- `backend/app/schemas/scoring.py` (ML readiness schemas only)

**Frontend** (redesign):
- `frontend/src/components/scoring/MLReadinessScore.tsx`
- `frontend/src/types/scoring.ts` (ML readiness types)

No new files needed.

## Key Reuse Map

| What We Need | Already Exists In |
|-------------|-------------------|
| Fragment detection | `deep_composition.py` → `MixtureDetectionCheck` |
| Metal/inorganic detection | `deep_composition.py` → `InorganicFilterCheck` + `METAL_ATOMIC_NUMS` |
| Radical detection | `deep_composition.py` → `RadicalDetectionCheck` |
| Isotope detection | `deep_composition.py` → `IsotopeLabelDetectionCheck` |
| Polymer/dummy atoms | `deep_complexity.py` → `PolymerDetectionCheck` |
| Formal charge analysis | `deep_complexity.py` → `ChargedSpeciesCheck` |
| Trivial molecule check | `deep_composition.py` → `TrivialMoleculeCheck` |
| Stereocenter analysis | `deep_stereo_tautomer.py` → `StereoisomerEnumerationCheck` |
| QED score | `druglikeness.py` → `_calculate_qed()` |
| Lipinski/Veber | `druglikeness.py` → multiple methods |
| SA Score | `admet.py` → `_calculate_synthetic_accessibility()` |
| Fsp3, Bertz, stereocenters | `admet.py` → `_calculate_complexity()` |
| Desirability function | `profile_scoring.py` → `desirability()` |
| Descriptor calculability | Current `ml_readiness.py` → `_score_descriptors()` |
| Fingerprint generation | Current `ml_readiness.py` → `_score_fingerprints()` |
| Property values (MW, LogP, TPSA, etc.) | Computed inline via RDKit one-liners |
