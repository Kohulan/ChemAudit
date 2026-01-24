# Getting Started with ChemStructVal

ChemStructVal is a web-based chemical structure validation platform. This guide will help you validate your first molecule in under 5 minutes.

## Quick Start

1. **Open the application** at `http://localhost:3000` (or your deployment URL)
2. **Enter a molecule** in the input field:
   - SMILES: `CCO` (ethanol)
   - InChI: `InChI=1S/C2H6O/c1-2-3/h3H,2H2,1H3`
   - Or paste a MOL block
3. **Click "Validate"**
4. **Review results** - score, issues, and molecule visualization

## Input Formats

ChemStructVal automatically detects the input format:

| Format | Example | Use Case |
|--------|---------|----------|
| SMILES | `CCO`, `c1ccccc1` | Quick input, common format |
| InChI | `InChI=1S/...` | Unique identifier, standardized |
| MOL/SDF | V2000/V3000 block | 3D coordinates, explicit H |

### SMILES Examples

```
CCO                     # Ethanol
c1ccccc1               # Benzene
C[C@H](O)N             # L-alaninol (with stereochemistry)
CC(=O)O                # Acetic acid
```

### InChI Examples

```
InChI=1S/C2H6O/c1-2-3/h3H,2H2,1H3           # Ethanol
InChI=1S/C6H6/c1-2-4-6-5-3-1/h1-6H          # Benzene
```

## What Gets Validated?

ChemStructVal performs 11+ validation checks:

### Basic Structure Checks
- **Parsability** - Can the molecule be read?
- **Sanitization** - Are implicit hydrogens valid?
- **Valence** - Are atom valences chemically reasonable?
- **Aromaticity** - Is aromaticity correctly assigned?
- **Connectivity** - Is the molecule a single fragment?

### Stereochemistry Checks
- **Undefined Stereocenters** - Are chiral centers defined?
- **Conflicting Stereochemistry** - Are stereo annotations consistent?
- **Tetrahedral Centers** - Proper chiral center detection

### Quality Checks
- **SMILES Roundtrip** - Does structure convert cleanly?
- **InChI Roundtrip** - Consistent InChI generation
- **Structural Alerts** - PAINS and BRENK pattern detection

### Advanced Scoring
- **ML-Readiness Score** - Suitability for machine learning
- **NP-Likeness Score** - Natural product similarity

## Understanding the Results

After validation, you'll see:

### Quality Score (0-100)
- **80-100**: High quality, suitable for ML/databases
- **50-79**: Acceptable with minor issues
- **<50**: Significant problems, review recommended

### Issues List
Each issue shows:
- **Severity**: CRITICAL, ERROR, WARNING, or INFO
- **Description**: What was found
- **Affected Atoms**: Highlighted in the molecule viewer
- **Recommendation**: How to fix

### Molecule Visualization
- 2D structure rendered with RDKit.js
- Hover over issues to highlight affected atoms
- Interactive viewer with zoom/pan

## Example: Validating Aspirin

1. Enter SMILES: `CC(=O)Oc1ccccc1C(=O)O`
2. Click **Validate**
3. Expected results:
   - Quality Score: ~90-95
   - No critical issues
   - High ML-readiness score
   - Successful roundtrip validation

## Example: Detecting Issues

Try this invalid SMILES: `C(C)(C)(C)(C)C`

Expected results:
- Quality Score: <50
- ERROR: Valence error (carbon with 5 bonds)
- Will not pass validation

## Next Steps

- [Validation Guide](validation.md) - Deep dive into validation checks
- [Batch Processing](batch-processing.md) - Process thousands of molecules
- [Understanding Results](understanding-results.md) - Interpret scores and alerts

## Common Workflows

### Research Chemist
1. Draw structure in ChemDraw
2. Copy SMILES
3. Paste into ChemStructVal
4. Review validation results
5. Fix any issues before database submission

### Computational Chemist
1. Export compound library as SDF
2. Use batch processing
3. Export results as CSV
4. Filter by quality score
5. Use clean dataset for modeling

### Database Curator
1. Upload existing database
2. Batch validate all structures
3. Identify problematic entries
4. Standardize using ChEMBL pipeline
5. Re-export cleaned dataset

## Tips for Best Results

1. **Use canonical SMILES** when possible for consistency
2. **Include stereochemistry** if known (use `@` and `@@`)
3. **Check salt stripping** before ML workflows
4. **Review warnings** - they may indicate intentional features
5. **Export batch results** for record-keeping

## Getting Help

If you encounter issues:
- Check the [FAQ](../troubleshooting/faq.md) for common problems
- Review your input format
- Check browser console for errors
- Verify RDKit.js loaded successfully
