# Validation Guide

This guide explains each validation check in detail and how to interpret results.

## Validation Checks

ChemStructVal performs comprehensive validation across four categories:

### Basic Checks

#### Parsability
**What it checks:** Whether the input can be parsed into a molecular structure.

**Severity:** CRITICAL (if fails)

**Common causes of failure:**
- Invalid SMILES syntax (e.g., unmatched brackets, invalid atoms)
- Malformed InChI string
- Corrupt MOL block (missing required fields)
- Unsupported format

**Example failures:**
```
C(C                    # Unmatched parenthesis
InChI=invalid          # Invalid InChI format
Zz                     # Invalid element symbol
```

**How to fix:** Check input format, verify syntax, ensure proper encoding.

#### Sanitization
**What it checks:** Whether implicit hydrogens can be added correctly and the molecule is chemically valid.

**Severity:** ERROR

**Common causes:**
- Invalid valence states
- Impossible ring systems
- Conflicting bond orders

**Example:**
```
C1CC2CCC1CC2           # Invalid bridged ring system
c1cccc1                # 5-membered aromatic ring (invalid)
```

**How to fix:** Review molecular structure, check ring systems, verify bond orders.

#### Valence Check
**What it checks:** Whether atoms have chemically reasonable valences.

**Severity:** ERROR

**Valid valences (RDKit defaults):**
- Carbon: 4
- Nitrogen: 3, 5
- Oxygen: 2
- Sulfur: 2, 4, 6
- Phosphorus: 3, 5

**Example failure:**
```
C(C)(C)(C)(C)C         # Carbon with 5 bonds (invalid)
```

**How to fix:** Verify atom connectivity, check for typos in SMILES.

#### Aromaticity Check
**What it checks:** Whether aromatic rings are correctly assigned using Hückel's rule.

**Severity:** WARNING

**Common issues:**
- Non-aromatic rings marked as aromatic
- Aromatic rings not detected

**Example:**
```
c1cccc1                # 5-membered aromatic (should be 6)
C1=CC=CC=C1            # Benzene not recognized as aromatic
```

**How to fix:** Use lowercase for aromatic atoms, verify ring size.

#### Connectivity Check
**What it checks:** Whether the molecule is a single connected fragment.

**Severity:** WARNING

**Why it matters:** Multiple fragments often indicate salts, mixtures, or disconnected structures.

**Example:**
```
CCO.Cl                 # Ethanol + chloride (2 fragments)
[Na+].CC(=O)[O-]       # Sodium acetate salt
```

**When it's OK:** Intentional salts, ionic compounds, coordination complexes.

**How to handle:**
- For ML workflows: use salt stripping
- For database storage: may keep as-is
- Document intention

### Stereochemistry Checks

#### Undefined Stereocenters
**What it checks:** Identifies chiral centers without defined stereochemistry.

**Severity:** WARNING

**Why it matters:** Different stereoisomers can have vastly different biological activity. The drug thalidomide is a famous example where one enantiomer is therapeutic and the other is teratogenic.

**Example:**
```
CC(O)N                 # Chiral center at carbon 2, undefined stereo
```

**With stereochemistry:**
```
C[C@H](O)N             # S-configuration
C[C@@H](O)N            # R-configuration
```

**When it's OK:**
- Racemates (intentional mixture)
- Early-stage research compounds
- When stereochemistry is unknown

**How to fix:** Add `@` or `@@` if stereochemistry is known.

#### Conflicting Stereochemistry
**What it checks:** Stereo annotations that contradict each other or are impossible.

**Severity:** ERROR

**Example issues:**
- Stereo assigned to non-chiral center
- Contradictory wedge bonds in MOL files
- Invalid E/Z assignment

**How to fix:** Review stereo assignments, check for structural errors.

#### E/Z Double Bond Stereochemistry
**What it checks:** Double bonds with defined stereochemistry.

**Notation:**
```
C/C=C/C                # E-configuration
C/C=C\C                # Z-configuration
```

**Warning:** Tautomer canonicalization can remove E/Z stereo. This is why it's OFF by default in standardization.

### Representation Checks

#### SMILES Roundtrip
**What it checks:** Input → canonical SMILES → back matches original structure.

**Severity:** INFO

**What it validates:**
- SMILES is valid
- Structure can be canonicalized
- No information loss

**Example:**
```
Input:  CC(C)C
Output: CC(C)C         # Success
```

**When it fails:**
- Input uses non-canonical notation
- Aromaticity perception differences
- Stereochemistry ambiguity

**Impact:** Usually informational, not a quality issue.

#### InChI Roundtrip
**What it checks:** Structure generates consistent InChI.

**Severity:** INFO

**What it validates:**
- InChI can be generated
- Structure is well-defined
- Unique identifier creation

**Example:**
```
SMILES: CCO
InChI:  InChI=1S/C2H6O/c1-2-3/h3H,2H2,1H3
```

**When it fails:**
- Complex stereochemistry
- Unusual valence states
- Large ring systems

**Note:** InChI differences are common with complex structures and don't necessarily indicate problems.

### Structural Alerts

#### PAINS Alerts
**What it checks:** Pan-Assay Interference Compounds patterns.

**Severity:** WARNING

**Categories:**
- **PAINS_A**: 16 patterns (most problematic)
- **PAINS_B**: 55 patterns (moderate concern)
- **PAINS_C**: 409 patterns (potential issues)

**Total patterns:** 480

**What they mean:**
PAINS are substructures that often give false positives in biochemical assays through non-specific mechanisms (aggregation, redox cycling, fluorescence interference).

**Important context:**
- 87 FDA-approved drugs contain PAINS patterns
- Context matters - not all PAINS are bad
- Consider your specific assay and use case

**Example PAINS:**
```
O=C1NC(=S)SC1          # Rhodanine (PAINS_A)
c1ccc2c(c1)oc(=O)cc2   # Coumarin (PAINS_B)
```

**How to interpret:**
- **Hit screening**: May want to filter out
- **Lead optimization**: Review carefully
- **Approved drugs**: Context-dependent

#### BRENK Alerts
**What it checks:** 105 patterns for known problematic functional groups.

**Severity:** WARNING

**Categories:**
- Reactive groups (e.g., acid halides, isocyanates)
- Known toxicophores (e.g., nitro aromatics)
- Metabolically unstable groups
- Chelating groups

**Example patterns:**
```
CC(=O)Cl               # Acid chloride (reactive)
c1ccc([N+](=O)[O-])cc1 # Nitro aromatic
```

**How to interpret:**
- Flag for chemist review
- May be acceptable in specific contexts
- Consider ADMET implications

## Quality Score Calculation

The overall score (0-100) reflects validation results.

### Scoring Formula

```
Base score: 100

Deductions:
- CRITICAL issue: -50 points each
- ERROR: -20 points each
- WARNING: -5 points each
- INFO: 0 points (no penalty)

Final score: max(0, base - deductions)
```

### Score Interpretation

| Score Range | Quality Level | Interpretation |
|-------------|---------------|----------------|
| 90-100 | Excellent | Publication-ready, ML-ready |
| 80-89 | Good | Minor issues, generally suitable |
| 70-79 | Acceptable | Review warnings, may need fixes |
| 50-69 | Fair | Multiple issues, review required |
| 30-49 | Poor | Significant problems, likely unusable |
| 0-29 | Critical | Structure invalid or severely flawed |

### Example Scoring

**Aspirin (high quality):**
```
SMILES: CC(=O)Oc1ccccc1C(=O)O
Issues: None
Score: 100
```

**Molecule with undefined stereo:**
```
SMILES: CC(O)C(O)C
Issues: 2 undefined stereocenters (WARNING × 2)
Score: 100 - 10 = 90
```

**Molecule with valence error:**
```
SMILES: C(C)(C)(C)(C)C
Issues: Valence error (ERROR)
Score: 100 - 20 = 80
```

**Multiple fragments with PAINS:**
```
SMILES: CCO.Cl contains rhodanine
Issues: Connectivity (WARNING), PAINS (WARNING)
Score: 100 - 10 = 90
```

## ML-Readiness Score

Separate from quality score, indicates suitability for machine learning.

### Components (weighted)

1. **Descriptors (40%)**: Can RDKit calculate 200+ molecular descriptors?
2. **Fingerprints (40%)**: Can Morgan/ECFP fingerprints be generated?
3. **Size (20%)**: Is molecular weight in typical range (50-1000)?

### Scoring

```
Descriptor success rate: 215/217 descriptors = 99% × 40 = 39.6
Fingerprint generation: Success = 100% × 40 = 40.0
Size score: MW=300, target=180-500 = 100% × 20 = 20.0

ML-Readiness: 39.6 + 40.0 + 20.0 = 99.6%
```

### What high score means
Your structure will work well in most ML models and cheminformatics pipelines.

### What low score means
The molecule may cause issues in ML workflows. Common causes:
- Large molecular weight (>1000 Da)
- Unusual elements (beyond common organic)
- Failed descriptor calculation
- No valid fingerprint

### When to care
- **Building ML models**: High score essential
- **Virtual screening**: High score recommended
- **Database curation**: Moderate score acceptable
- **Novel chemistry**: Low score may be expected

## NP-Likeness Score

Indicates similarity to natural products vs synthetic molecules.

### Interpretation

| Score | Meaning | Example Compounds |
|-------|---------|-------------------|
| > +2 | Very natural product-like | Morphine, taxol |
| 0 to +2 | Somewhat NP-like | Simple alkaloids |
| 0 | Borderline | Aromatic heterocycles |
| -2 to 0 | Somewhat synthetic | Fluorinated drugs |
| < -2 | Very synthetic-like | PEGs, dendrimers |

### Algorithm

Uses fragment-based scoring comparing:
- Fragments present vs NP databases
- Structural complexity
- Element composition
- Ring systems

**Note:** ChemStructVal uses a heuristic version. For research-grade NP scoring, consider the full NP model.

### Use Cases

- **Natural product screening**: Filter for positive scores
- **Library design**: Balance synthetic accessibility with NP-like features
- **Chemical space analysis**: Understand library composition

## Taking Action on Results

| Issue Type | Recommended Action |
|------------|-------------------|
| **Parsability error** | Fix syntax, check format |
| **Valence error** | Review structure, likely incorrect |
| **Sanitization failure** | Check ring systems and bonds |
| **Undefined stereo** | Add if known, document if intentional |
| **Conflicting stereo** | Fix annotations |
| **PAINS alert** | Review context, may keep if appropriate |
| **BRENK alert** | Flag for chemist review |
| **Multiple fragments** | Consider salt stripping |
| **Low ML-readiness** | Check size and complexity |

## Advanced: Affected Atoms

Many issues identify specific atoms:

```json
{
  "check_name": "undefined_stereocenters",
  "severity": "WARNING",
  "affected_atoms": [2, 5]
}
```

**How to use:**
- Hover over issue in UI to highlight atoms
- Atom numbering is 0-indexed
- Helps locate problems in large molecules

## Batch Validation Considerations

When validating thousands of molecules:

1. **Set thresholds**: Filter by minimum quality score
2. **Review distributions**: Check for systematic issues
3. **Export results**: Keep audit trail
4. **Handle failures**: Individual failures don't stop batch
5. **Partial results**: Better than nothing

See [Batch Processing Guide](batch-processing.md) for details.
