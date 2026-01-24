# Understanding Validation Results

This guide helps you interpret validation results and take appropriate action.

## Result Overview

After validating a molecule, you receive:

1. **Quality Score** (0-100): Overall structural quality
2. **ML-Readiness Score** (0-100%): Suitability for machine learning
3. **NP-Likeness Score** (-∞ to +∞): Natural product similarity
4. **Issues List**: Specific validation problems
5. **Structural Alerts**: PAINS and BRENK patterns detected
6. **Molecule Visualization**: 2D structure with highlighting

## Quality Score Deep Dive

### Score Components

The quality score aggregates all validation checks:

```
Base score: 100 points

Deductions per issue:
- CRITICAL (parse failure):     -50 points
- ERROR (chemical invalidity):  -20 points
- WARNING (quality concerns):    -5 points
- INFO (informational):           0 points

Minimum score: 0 (clamped, no negative scores)
```

### Severity Levels Explained

#### CRITICAL
**Meaning:** Structure cannot be processed at all.

**Examples:**
- Parse failure (invalid SMILES/InChI syntax)
- Completely invalid input

**Action required:** Fix input format.

**Impact on workflow:** Cannot proceed with this structure.

#### ERROR
**Meaning:** Chemically invalid or highly problematic.

**Examples:**
- Valence errors (carbon with 5 bonds)
- Sanitization failure
- Conflicting stereochemistry

**Action required:** Review and fix structure.

**Impact on workflow:** Unusable for most applications.

#### WARNING
**Meaning:** Potentially problematic, requires review.

**Examples:**
- Undefined stereocenters
- Multiple fragments (salts)
- Structural alerts (PAINS/BRENK)

**Action required:** Review context, may be acceptable.

**Impact on workflow:** Usable but may need attention.

#### INFO
**Meaning:** Informational, not a quality issue.

**Examples:**
- Roundtrip checks
- Standardization notes
- Descriptor calculations

**Action required:** None, for awareness only.

**Impact on workflow:** No impact.

### Score Ranges and Actions

| Score | Quality | Interpretation | Action |
|-------|---------|----------------|--------|
| 95-100 | Excellent | Perfect structure | Use as-is |
| 85-94 | Very Good | Minor warnings only | Review warnings |
| 75-84 | Good | Some concerns | Investigate issues |
| 60-74 | Fair | Multiple warnings | Consider fixes |
| 40-59 | Poor | Significant issues | Likely needs fixes |
| 20-39 | Very Poor | Major problems | Probably unusable |
| 0-19 | Critical | Invalid structure | Must fix or discard |

### Real-World Examples

#### Example 1: Aspirin (Perfect)
```
Input: CC(=O)Oc1ccccc1C(=O)O

Results:
- Quality Score: 100
- Issues: None
- ML-Readiness: 99.2%
- NP-Likeness: -0.8 (synthetic-like)

Interpretation: Publication-ready, suitable for any use.
```

#### Example 2: Molecule with Undefined Stereo
```
Input: CC(O)C(O)C(O)CO

Results:
- Quality Score: 85 (100 - 3×5 = 85)
- Issues: 3 undefined stereocenters (WARNING)
- ML-Readiness: 98.5%
- Affected atoms: [1, 2, 3]

Interpretation:
- Structurally valid
- 8 stereoisomers possible (2³)
- For racemate: acceptable as-is
- For specific isomer: add stereochemistry

Action:
- If stereochemistry known: add @/@@
- If racemate intended: document and accept
- For ML: may want to enumerate isomers
```

#### Example 3: Salt with PAINS
```
Input: Cl.CN1C=NC2=C1C(=O)N(C(=O)N2C)C

Results:
- Quality Score: 90 (100 - 5 - 5 = 90)
- Issues:
  - Multiple fragments (WARNING): 2 components
  - PAINS_B alert (WARNING): imidazole pattern
- ML-Readiness: 95.0%

Interpretation:
- Caffeine hydrochloride salt
- PAINS alert is false positive (caffeine is safe)
- Salt may need stripping for ML

Action:
- For database: keep salt
- For ML: strip to neutral form
- PAINS: ignore (known safe compound)
```

#### Example 4: Valence Error
```
Input: C(C)(C)(C)(C)C

Results:
- Quality Score: 80 (100 - 20 = 80)
- Issues: Valence error (ERROR)
- Cannot calculate ML-Readiness

Interpretation:
- Carbon with 5 bonds (impossible)
- Likely input error
- Structure is invalid

Action:
- Fix SMILES (probably meant CC(C)(C)C)
- Verify original structure
```

## ML-Readiness Score Deep Dive

### What It Measures

ML-Readiness indicates if a molecule will work well in machine learning pipelines.

### Components (Weighted)

#### 1. Descriptors (40%)
**Checks:** Can RDKit calculate 217 molecular descriptors?

**Measured:**
```
Success rate = (calculated descriptors) / 217
Descriptor score = success_rate × 40
```

**Example:**
- 217/217 calculated: 100% × 40 = 40.0 points
- 200/217 calculated: 92% × 40 = 36.8 points

**Common failures:**
- Very large molecules (>1000 atoms)
- Unusual elements (beyond common organic)
- Exotic valence states

#### 2. Fingerprints (40%)
**Checks:** Can Morgan (ECFP) fingerprints be generated?

**Measured:**
```
Success = can generate fingerprint
Fingerprint score = 40 if success, else 0
```

**Example:**
- Success: 40.0 points
- Failure: 0.0 points

**Common failures:**
- Invalid molecule (failed sanitization)
- Edge cases in fingerprint algorithm

#### 3. Size (20%)
**Checks:** Is molecular weight in typical ML training range?

**Measured:**
```
Ideal range: 180-500 Da (gets 100% of points)
Acceptable range: 50-1000 Da (gets partial points)
Outside range: Gets reduced points
```

**Example:**
- MW = 300 Da: 100% × 20 = 20.0 points (ideal)
- MW = 150 Da: 80% × 20 = 16.0 points (small)
- MW = 800 Da: 60% × 20 = 12.0 points (large)
- MW = 1500 Da: 20% × 20 = 4.0 points (very large)

### Total ML-Readiness Score

```
ML-Readiness = Descriptor_Score + Fingerprint_Score + Size_Score

Maximum: 40 + 40 + 20 = 100
```

### Interpreting ML-Readiness

| Score | Suitability | Action |
|-------|-------------|--------|
| 95-100% | Excellent | Perfect for ML |
| 85-94% | Good | Minor issues, usually OK |
| 70-84% | Fair | May work with preprocessing |
| 50-69% | Poor | Significant limitations |
| <50% | Unsuitable | Avoid for ML |

### When ML-Readiness Matters

**Critical for:**
- Building QSAR models
- Virtual screening
- Property prediction
- Similarity searches

**Less important for:**
- Database storage
- Visualization
- Literature references
- Early discovery research

### Improving Low ML-Readiness

**If low descriptor score:**
- Check for unusual atoms
- Verify structure validity
- Consider simplifying structure

**If low fingerprint score:**
- Structure is likely invalid
- Fix sanitization errors
- Check valence issues

**If low size score:**
- For large molecules: consider fragmenting
- For small molecules: acceptable for some ML tasks
- Document and proceed if intentional

## NP-Likeness Score Deep Dive

### What It Measures

Natural Product Likeness indicates similarity to natural products vs synthetic molecules.

### Scale

```
-∞ ←───────── 0 ─────────→ +∞
Synthetic              Natural Product
```

Typical range: -5 to +5

### Interpretation Table

| Score Range | Classification | Example Compounds |
|-------------|----------------|-------------------|
| +3 to +5 | Very NP-like | Morphine, taxol, strychnine |
| +1 to +3 | NP-like | Simple alkaloids, terpenes |
| -1 to +1 | Borderline | Small heterocycles |
| -3 to -1 | Synthetic-like | Most drugs, fluorinated |
| -5 to -3 | Very synthetic | PEGs, dendrimers, polymers |

### Algorithm (Heuristic Version)

ChemStructVal uses a fragment-based heuristic:

```
Score factors:
+ Fragments common in natural products
+ Structural complexity (rings, stereo)
+ Sp³ carbons (natural products are more saturated)
- Fragments rare in NPs (e.g., fluorine)
- Very simple structures
- Common synthetic motifs
```

**Note:** This is a simplified heuristic. For research-grade scoring, use the full NP model with trained parameters.

### Use Cases

#### Natural Product Research
**Goal:** Find NP-like compounds in library

**Filter:** Score > +1

**Example:**
```
Compound A: NP-score = +2.5 → Keep
Compound B: NP-score = -1.2 → Discard
```

#### Library Design
**Goal:** Balance synthetic accessibility with NP-like properties

**Target:** Score around 0 to +1

**Rationale:** Slight NP-likeness, but synthetically accessible

#### Chemical Space Analysis
**Goal:** Understand library composition

**Analysis:**
```
- Mean NP-score: -0.5 (synthetic-leaning)
- Range: -3.2 to +2.8 (good diversity)
- 30% with score > 0 (NP-like)
```

### NP-Likeness Limitations

**Not a quality metric:**
- Low score ≠ bad molecule
- High score ≠ good molecule
- Contextual interpretation required

**Heuristic version:**
- Less accurate than full model
- Good for relative comparisons
- Not absolute predictions

## Structural Alerts Deep Dive

### PAINS (Pan-Assay Interference Compounds)

#### What They Are
Substructures that frequently cause false positives in biochemical assays through non-specific mechanisms.

#### Three Categories

**PAINS_A (16 patterns - highest concern):**
- Rhodanines
- Hydroxyphenyl hydrazones
- Alkylidene barbiturates

**PAINS_B (55 patterns - moderate concern):**
- Catechols
- Quinones
- Michael acceptors

**PAINS_C (409 patterns - lower concern):**
- Various heterocycles
- Reactive groups

#### Important Context

**87 FDA-approved drugs contain PAINS patterns.**

This means PAINS are warnings, not disqualifications.

#### Decision Framework

**Filter out if:**
- High-throughput screening campaign
- No prior knowledge of compound
- Many alternatives available

**Review carefully if:**
- Lead optimization stage
- Known mechanism of action
- Specific target context

**Likely OK if:**
- Approved drug or analog
- Well-characterized compound
- Non-assay application (e.g., imaging)

#### Example: Rhodanine

```
Structure: O=C1NC(=S)SC1
PAINS: Category A (high concern)

Known issues:
- Aggregation
- Metal chelation
- Redox cycling

But also found in:
- Some anti-diabetic compounds
- Historical use in medicinal chemistry

Decision: Context-dependent
```

### BRENK (Structural Alerts)

#### What They Are
105 patterns for known problematic functional groups.

#### Categories

**Reactive groups:**
- Acid halides (rapid hydrolysis)
- Isocyanates (protein modification)
- Epoxides (DNA alkylation)

**Toxicophores:**
- Nitro aromatics (mutagenicity)
- Hydrazines (hepatotoxicity)

**Metabolically unstable:**
- Nitriles (CYP metabolism)
- Thiols (oxidation)

**Chelating groups:**
- β-diketones
- Hydroxamic acids

#### Decision Framework

**Acceptable if:**
- Prodrug design (intentional reactivity)
- Known safe analog
- Controlled context (e.g., injection)

**Concern if:**
- Oral dosing
- Chronic administration
- Lack of safety data

**Likely problematic if:**
- Multiple alerts
- Early discovery
- No specific justification

## Taking Action

### Decision Tree

```
Is quality score < 50?
├─ YES → Review ERROR issues
│         └─ Can fix? → Fix and re-validate
│         └─ Cannot fix? → Discard
│
└─ NO → Is ML-Readiness < 70% and you need ML?
       ├─ YES → Investigate descriptor failures
       │         └─ Structure too large/complex?
       │
       └─ NO → Are there PAINS/BRENK alerts?
              ├─ YES → Review context
              │         └─ Acceptable? → Document and proceed
              │         └─ Not acceptable? → Discard or modify
              │
              └─ NO → Proceed with confidence!
```

### Action Matrix

| Issue | Severity | For Database | For ML | For Assay Screening |
|-------|----------|--------------|--------|---------------------|
| Parse failure | CRITICAL | Fix/discard | Fix/discard | Fix/discard |
| Valence error | ERROR | Fix | Fix | Fix |
| Undefined stereo | WARNING | Document | Enumerate? | OK if racemate |
| Multiple fragments | WARNING | Keep or strip | Strip salts | Strip salts |
| PAINS alert | WARNING | Keep | Review context | Filter out |
| BRENK alert | WARNING | Keep | Review context | Flag for review |
| Low ML-readiness | INFO | OK | Investigate | OK |

### Workflow Examples

#### Workflow 1: Preparing ML Training Set
```
1. Filter quality_score >= 80
2. Filter ml_readiness >= 90
3. Review PAINS alerts (may keep if known)
4. Strip salts (multiple fragments)
5. For undefined stereo:
   - Enumerate stereoisomers, OR
   - Keep racemate (document)
6. Export clean dataset
```

#### Workflow 2: Database Curation
```
1. Fix CRITICAL and ERROR issues
2. Keep WARNING issues (document)
3. Preserve salts and stereochemistry as-is
4. Add validation metadata
5. Flag low-quality entries for review
```

#### Workflow 3: Screening Library Design
```
1. Filter quality_score >= 70
2. Remove PAINS_A and PAINS_B
3. Review BRENK alerts (chemist decision)
4. Ensure ml_readiness >= 80 (for later modeling)
5. Document acceptable warnings
```

## Batch Results Analysis

When validating thousands of molecules:

### Summary Statistics

```
Total: 10,000 molecules
Success: 9,500 (95%)
Failed: 500 (5%)

Quality score distribution:
- Excellent (90-100): 7,200 (72%)
- Good (80-89): 1,800 (18%)
- Fair (70-79): 400 (4%)
- Poor (<70): 100 (1%)

Common issues:
- Undefined stereo: 3,500 (35%)
- PAINS alerts: 850 (8.5%)
- Multiple fragments: 1,200 (12%)
```

### Identifying Systematic Issues

**High parse failure rate (>10%):**
- Check input format consistency
- Verify encoding (UTF-8)
- Review data source

**Many undefined stereo warnings:**
- Expected for racemates
- Consider if stereochemistry matters
- May enumerate if needed

**High PAINS rate (>20%):**
- May indicate library bias
- Review sourcing
- Consider if acceptable for use case

### Quality Control Thresholds

Recommended filters for different applications:

**For ML modeling:**
```
quality_score >= 80
ml_readiness >= 90
no CRITICAL or ERROR
```

**For database storage:**
```
quality_score >= 60
no CRITICAL or ERROR
```

**For screening:**
```
quality_score >= 70
filter PAINS_A and PAINS_B
review BRENK alerts
```

## Exporting and Sharing Results

### What to Include in Reports

**For collaborators:**
- Quality score
- Key issues
- Recommended actions

**For records:**
- Full validation results
- Input structure
- Timestamp
- ChemStructVal version

**For publications:**
- Validation criteria used
- Filtering thresholds
- Number of molecules processed
- Quality distribution

### Documentation Best Practices

1. **Track validation runs** with timestamps
2. **Document filter criteria** and rationale
3. **Keep original inputs** and results
4. **Version control** datasets
5. **Include metadata** in exports

## Getting Help

If results are unexpected:

1. **Check FAQ** for common issues
2. **Review input format** for errors
3. **Verify RDKit version** compatibility
4. **Compare with reference** molecules
5. **Consult documentation** for specific checks

See [FAQ](../troubleshooting/faq.md) for troubleshooting.
