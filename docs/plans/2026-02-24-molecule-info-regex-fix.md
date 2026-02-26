# Fix Regex-Based Molecule Info Bugs — Implementation Plan

> **For Claude:** REQUIRED SUB-SKILL: Use superpowers:executing-plans to implement this plan task-by-task.

**Goal:** Replace all fragile regex-based SMILES parsing in the frontend with proper RDKit.js `get_descriptors()` for preview and authoritative backend data post-validation.

**Architecture:** Backend `MoleculeInfo` schema is enriched with `num_bonds`, `num_rings`, `num_aromatic_rings`, `num_stereocenters`, `has_stereochemistry`. Frontend `useMoleculeInfo` hook is rewritten to use RDKit.js `get_descriptors()` instead of regex. `SingleValidation.tsx` uses backend data exclusively after validation.

**Tech Stack:** Python/FastAPI/RDKit (backend), TypeScript/React/RDKit.js (frontend)

**Branch:** `fix/molecule-info-regex-bugs` from `development`

---

### Task 1: Create branch

**Step 1: Create branch from development**

```bash
git checkout development
git pull origin development
git checkout -b fix/molecule-info-regex-bugs
```

**Step 2: Verify branch**

```bash
git branch --show-current
```

Expected: `fix/molecule-info-regex-bugs`

---

### Task 2: Backend — Write failing tests for enriched MoleculeInfo

**Files:**
- Modify: `backend/tests/test_api/test_validate_endpoint.py`

**Step 1: Write failing tests**

Add a new test class at the end of the file:

```python
class TestMoleculeInfoEnrichedFields:
    """Test that validation response includes enriched molecule info fields."""

    @pytest.mark.asyncio
    async def test_enriched_fields_simple_molecule(self, client):
        """Ethanol should return correct enriched molecular properties."""
        response = await client.post(
            "/api/v1/validate", json={"molecule": "CCO", "format": "smiles"}
        )
        assert response.status_code == 200
        info = response.json()["molecule_info"]

        assert info["num_atoms"] == 3
        assert info["num_bonds"] == 2
        assert info["num_rings"] == 0
        assert info["num_aromatic_rings"] == 0
        assert info["num_stereocenters"] == 0
        assert info["has_stereochemistry"] is False

    @pytest.mark.asyncio
    async def test_enriched_fields_benzene(self, client):
        """Benzene should have 1 ring, 1 aromatic ring."""
        response = await client.post(
            "/api/v1/validate", json={"molecule": "c1ccccc1"}
        )
        assert response.status_code == 200
        info = response.json()["molecule_info"]

        assert info["num_atoms"] == 6
        assert info["num_bonds"] == 6
        assert info["num_rings"] == 1
        assert info["num_aromatic_rings"] == 1
        assert info["num_stereocenters"] == 0
        assert info["has_stereochemistry"] is False

    @pytest.mark.asyncio
    async def test_enriched_fields_isotope_no_false_rings(self, client):
        """Isotope labels must NOT inflate ring count."""
        response = await client.post(
            "/api/v1/validate",
            json={"molecule": "O=[CH][63Ni]([CH]=O)([CH]=O)[13CH]=O"},
        )
        assert response.status_code == 200
        info = response.json()["molecule_info"]

        assert info["num_rings"] == 0, "Isotope digits must not be counted as rings"
        assert info["num_aromatic_rings"] == 0

    @pytest.mark.asyncio
    async def test_enriched_fields_isotope_with_ring(self, client):
        """C-13 cyclohexane: 1 real ring, isotope digit must not inflate."""
        response = await client.post(
            "/api/v1/validate", json={"molecule": "[13C]1CCCCC1"}
        )
        assert response.status_code == 200
        info = response.json()["molecule_info"]

        assert info["num_atoms"] == 6
        assert info["num_rings"] == 1
        assert info["num_aromatic_rings"] == 0

    @pytest.mark.asyncio
    async def test_enriched_fields_stereocenters(self, client):
        """Alanine has 1 stereocenter."""
        response = await client.post(
            "/api/v1/validate", json={"molecule": "N[C@@H](C)C(=O)O"}
        )
        assert response.status_code == 200
        info = response.json()["molecule_info"]

        assert info["num_stereocenters"] == 1
        assert info["has_stereochemistry"] is True

    @pytest.mark.asyncio
    async def test_enriched_fields_ez_stereo(self, client):
        """cis-2-butene has E/Z stereochemistry but no tetrahedral stereocenters."""
        response = await client.post(
            "/api/v1/validate", json={"molecule": r"C/C=C\C"}
        )
        assert response.status_code == 200
        info = response.json()["molecule_info"]

        assert info["num_stereocenters"] == 0
        assert info["has_stereochemistry"] is True

    @pytest.mark.asyncio
    async def test_enriched_fields_aspirin(self, client):
        """Aspirin: 1 ring, 1 aromatic ring, no stereo."""
        response = await client.post(
            "/api/v1/validate", json={"molecule": "CC(=O)Oc1ccccc1C(=O)O"}
        )
        assert response.status_code == 200
        info = response.json()["molecule_info"]

        assert info["num_atoms"] == 13
        assert info["num_bonds"] == 13
        assert info["num_rings"] == 1
        assert info["num_aromatic_rings"] == 1
        assert info["num_stereocenters"] == 0
        assert info["has_stereochemistry"] is False

    @pytest.mark.asyncio
    async def test_enriched_fields_heavy_water(self, client):
        """Heavy water [2H]O[2H]: 0 rings despite isotope digits."""
        response = await client.post(
            "/api/v1/validate", json={"molecule": "[2H]O[2H]"}
        )
        assert response.status_code == 200
        info = response.json()["molecule_info"]

        assert info["num_atoms"] == 3
        assert info["num_bonds"] == 2
        assert info["num_rings"] == 0
```

**Step 2: Run tests to verify they fail**

```bash
cd backend
python -m pytest tests/test_api/test_validate_endpoint.py::TestMoleculeInfoEnrichedFields -v
```

Expected: FAIL — `KeyError: 'num_bonds'` (fields don't exist yet)

**Step 3: Commit failing tests**

```bash
git add tests/test_api/test_validate_endpoint.py
git commit -m "test: add failing tests for enriched MoleculeInfo fields"
```

---

### Task 3: Backend — Enrich MoleculeInfo schema and extract_molecule_info()

**Files:**
- Modify: `backend/app/schemas/validation.py:25-34` (add fields to MoleculeInfo)
- Modify: `backend/app/api/routes/validation.py:263-315` (compute new fields)

**Step 1: Add fields to MoleculeInfo schema**

In `backend/app/schemas/validation.py`, add 5 fields after `num_atoms`:

```python
class MoleculeInfo(BaseModel):
    """Schema for molecule information."""

    input_smiles: str
    canonical_smiles: Optional[str] = None
    inchi: Optional[str] = None
    inchikey: Optional[str] = None
    molecular_formula: Optional[str] = None
    molecular_weight: Optional[float] = None
    num_atoms: Optional[int] = None
    num_bonds: Optional[int] = None
    num_rings: Optional[int] = None
    num_aromatic_rings: Optional[int] = None
    num_stereocenters: Optional[int] = None
    has_stereochemistry: Optional[bool] = None
```

**Step 2: Update extract_molecule_info() to compute new fields**

In `backend/app/api/routes/validation.py`, update `extract_molecule_info()`:

```python
def extract_molecule_info(
    mol: Chem.Mol, input_smiles: str, preserve_aromatic: bool = False
) -> MoleculeInfo:
    """
    Extract molecule properties.

    Args:
        mol: RDKit molecule object
        input_smiles: Original input string
        preserve_aromatic: If True, output aromatic SMILES notation (lowercase atoms)

    Returns:
        MoleculeInfo with molecular properties
    """
    try:
        if preserve_aromatic:
            canonical = Chem.MolToSmiles(mol)
        else:
            try:
                Chem.Kekulize(mol, clearAromaticFlags=False)
                canonical = Chem.MolToSmiles(mol, kekuleSmiles=True)
            except Exception:
                canonical = Chem.MolToSmiles(mol)

        mol_inchi = rdkit_inchi.MolToInchi(mol)
        mol_inchikey = rdkit_inchi.MolToInchiKey(mol) if mol_inchi else None
        formula = rdMolDescriptors.CalcMolFormula(mol)
        mw = Descriptors.MolWt(mol)
        num_atoms = mol.GetNumAtoms()
        num_bonds = mol.GetNumBonds()
        num_rings = rdMolDescriptors.CalcNumRings(mol)
        num_aromatic_rings = Descriptors.NumAromaticRings(mol)

        # Stereocenter detection
        chiral_centers = Chem.FindMolChiralCenters(mol, includeUnassigned=True)
        num_stereocenters = len(chiral_centers)

        # E/Z bond detection
        has_ez = any(
            bond.GetStereo() != Chem.BondStereo.STEREONONE
            for bond in mol.GetBonds()
            if bond.GetBondTypeAsDouble() == 2.0
        )
        has_stereochemistry = num_stereocenters > 0 or has_ez

    except Exception:
        canonical = None
        mol_inchi = None
        mol_inchikey = None
        formula = None
        mw = None
        num_atoms = None
        num_bonds = None
        num_rings = None
        num_aromatic_rings = None
        num_stereocenters = None
        has_stereochemistry = None

    return MoleculeInfo(
        input_smiles=input_smiles,
        canonical_smiles=canonical,
        inchi=mol_inchi,
        inchikey=mol_inchikey,
        molecular_formula=formula,
        molecular_weight=mw,
        num_atoms=num_atoms,
        num_bonds=num_bonds,
        num_rings=num_rings,
        num_aromatic_rings=num_aromatic_rings,
        num_stereocenters=num_stereocenters,
        has_stereochemistry=has_stereochemistry,
    )
```

Note: `Chem.FindMolChiralCenters` and `Descriptors.NumAromaticRings` are already imported at the top of the file. Verify that `Chem.FindMolChiralCenters` is available — if not, add `from rdkit.Chem import FindMolChiralCenters` or use `Chem.FindMolChiralCenters`.

**Step 3: Run tests to verify they pass**

```bash
cd backend
python -m pytest tests/test_api/test_validate_endpoint.py::TestMoleculeInfoEnrichedFields -v
```

Expected: ALL PASS

**Step 4: Run full validation test suite to check for regressions**

```bash
cd backend
python -m pytest tests/test_api/test_validate_endpoint.py -v
```

Expected: ALL PASS

**Step 5: Commit**

```bash
git add backend/app/schemas/validation.py backend/app/api/routes/validation.py
git commit -m "feat: enrich MoleculeInfo with bonds, rings, stereo fields"
```

---

### Task 4: Frontend — Update TypeScript types to match backend

**Files:**
- Modify: `frontend/src/types/validation.ts:12-20`

**Step 1: Add new fields to MoleculeInfo interface**

```typescript
export interface MoleculeInfo {
  input_smiles: string;
  canonical_smiles: string | null;
  inchi: string | null;
  inchikey: string | null;
  molecular_formula: string | null;
  molecular_weight: number | null;
  num_atoms: number | null;
  num_bonds: number | null;
  num_rings: number | null;
  num_aromatic_rings: number | null;
  num_stereocenters: number | null;
  has_stereochemistry: boolean | null;
}
```

**Step 2: Commit**

```bash
git add frontend/src/types/validation.ts
git commit -m "feat: add enriched fields to MoleculeInfo TypeScript interface"
```

---

### Task 5: Frontend — Rewrite useMoleculeInfo.ts with get_descriptors()

**Files:**
- Modify: `frontend/src/hooks/useMoleculeInfo.ts` (full rewrite)
- Modify: `frontend/src/hooks/useRDKit.ts:9-15` (add `get_descriptors` to interface)

**Step 1: Add `get_descriptors` to RDKitMol interface**

In `frontend/src/hooks/useRDKit.ts`, add `get_descriptors` to the `RDKitMol` interface:

```typescript
interface RDKitMol {
  get_svg: (width?: number, height?: number) => string;
  get_svg_with_highlights: (details: string) => string;
  get_smiles: (opts?: string) => string;
  get_molblock: () => string;
  get_descriptors: () => string;
  get_json: () => string;
  delete: () => void;  // CRITICAL: Must call to free WASM memory
}
```

**Step 2: Rewrite useMoleculeInfo.ts**

Replace the entire file with:

```typescript
import { useState, useEffect, useRef } from 'react';
import { getRDKit, RDKitMol } from './useRDKit';

interface MoleculeInfo {
  canonicalSmiles: string;
  kekulizedSmiles: string | null;
  numAtoms: number;
  numBonds: number;
  numRings: number;
  numAromaticRings: number;
  numStereocenters: number;
  hasEZStereo: boolean;
  hasStereochemistry: boolean;
  isValid: boolean;
}

interface UseMoleculeInfoResult {
  info: MoleculeInfo | null;
  isLoading: boolean;
  error: string | null;
}

/**
 * Extract molecular properties from an RDKit.js mol object using get_descriptors().
 * No regex parsing of SMILES — all values come from RDKit's cheminformatics engine.
 */
function extractProperties(mol: RDKitMol): {
  numAtoms: number;
  numBonds: number;
  numRings: number;
  numAromaticRings: number;
  numStereocenters: number;
  hasEZStereo: boolean;
} {
  const descriptorsJson = mol.get_descriptors();
  const desc = JSON.parse(descriptorsJson);

  const numAtoms: number = desc.NumHeavyAtoms ?? 0;
  const numRings: number = desc.NumRings ?? 0;
  const numAromaticRings: number = desc.NumAromaticRings ?? 0;

  // get_descriptors doesn't include bond count or stereo — get from JSON
  let numBonds = 0;
  let numStereocenters = 0;
  let hasEZStereo = false;

  try {
    const molJson = mol.get_json();
    const parsed = JSON.parse(molJson);
    const molData = parsed.molecules?.[0];
    numBonds = molData?.bonds?.length ?? 0;

    // Count stereocenters from atom stereo annotations
    const atoms = molData?.atoms ?? [];
    for (const atom of atoms) {
      if (atom.stereo === 'CCW' || atom.stereo === 'CW') {
        numStereocenters++;
      }
    }

    // Check E/Z stereochemistry from bond stereo annotations
    const bonds = molData?.bonds ?? [];
    for (const bond of bonds) {
      if (
        bond.stereo === 'either' ||
        bond.stereo === 'cis' ||
        bond.stereo === 'trans'
      ) {
        hasEZStereo = true;
        break;
      }
    }
  } catch {
    // get_json not available — bonds stay 0
  }

  return { numAtoms, numBonds, numRings, numAromaticRings, numStereocenters, hasEZStereo };
}

/**
 * Hook to extract basic molecule information using RDKit.js.
 * Shows info immediately when a valid molecule is entered.
 *
 * All molecular properties are computed by RDKit — no regex parsing of SMILES.
 */
export function useMoleculeInfo(smiles: string | null): UseMoleculeInfoResult {
  const [info, setInfo] = useState<MoleculeInfo | null>(null);
  const [isLoading, setIsLoading] = useState(false);
  const [error, setError] = useState<string | null>(null);
  const molRef = useRef<RDKitMol | null>(null);

  useEffect(() => {
    if (!smiles || smiles.trim() === '') {
      setInfo(null);
      setError(null);
      return;
    }

    let cancelled = false;
    setIsLoading(true);
    setError(null);

    getRDKit()
      .then((rdkit) => {
        if (cancelled) return;

        // Clean up previous molecule
        if (molRef.current) {
          molRef.current.delete();
          molRef.current = null;
        }

        try {
          const mol = rdkit.get_mol(smiles);

          if (!mol) {
            setError('Invalid molecule structure');
            setInfo(null);
            setIsLoading(false);
            return;
          }

          molRef.current = mol;

          const canonicalSmiles = mol.get_smiles();

          // Try RDKit.js kekulization — no regex fallback
          let kekulizedSmiles: string | null = null;
          try {
            const result = mol.get_smiles(JSON.stringify({ kekulize: true }));
            if (result && result !== canonicalSmiles) {
              kekulizedSmiles = result;
            }
          } catch {
            // Kekulize option not available in this RDKit.js build
          }

          // Extract all properties via RDKit — no regex
          const props = extractProperties(mol);

          const hasStereochemistry = props.numStereocenters > 0 || props.hasEZStereo;

          setInfo({
            canonicalSmiles,
            kekulizedSmiles,
            numAtoms: props.numAtoms,
            numBonds: props.numBonds,
            numRings: props.numRings,
            numAromaticRings: props.numAromaticRings,
            numStereocenters: props.numStereocenters,
            hasEZStereo: props.hasEZStereo,
            hasStereochemistry,
            isValid: true,
          });
          setError(null);
        } catch (e) {
          setError(e instanceof Error ? e.message : 'Failed to parse molecule');
          setInfo(null);
        } finally {
          setIsLoading(false);
        }
      })
      .catch((err) => {
        if (!cancelled) {
          setError(err.message);
          setIsLoading(false);
        }
      });

    return () => {
      cancelled = true;
      if (molRef.current) {
        molRef.current.delete();
        molRef.current = null;
      }
    };
  }, [smiles]);

  return { info, isLoading, error };
}
```

**Step 3: Verify frontend compiles**

```bash
cd frontend
npx tsc --noEmit 2>&1 | head -30
```

Expected: No errors related to useMoleculeInfo or MoleculeInfo types. If there are type errors from SingleValidation.tsx (due to the new `numAromaticRings` and `hasEZStereo` fields), those will be fixed in Task 6.

**Step 4: Commit**

```bash
git add frontend/src/hooks/useMoleculeInfo.ts frontend/src/hooks/useRDKit.ts
git commit -m "fix: replace regex SMILES parsing with RDKit.js get_descriptors()"
```

---

### Task 6: Frontend — Update SingleValidation.tsx to use backend data post-validation

**Files:**
- Modify: `frontend/src/pages/SingleValidation.tsx`

**Step 1: Update post-validation info tiles (lines 570-598)**

Replace the mixed-source stats grid with backend-only data. The section starting at line 570 `{result ? (` should use `result.molecule_info` for ALL tiles:

```tsx
{result ? (
  <div className="space-y-3">
    {/* Stats row — all from backend */}
    <div className="grid grid-cols-5 gap-2 text-center">
      {result.molecule_info.num_atoms != null && (
        <div className="bg-[var(--color-surface-sunken)] rounded-lg p-2">
          <div className="text-lg font-bold text-[var(--color-text-primary)]">{result.molecule_info.num_atoms}</div>
          <div className="text-[10px] text-[var(--color-text-muted)] uppercase tracking-wider">Atoms</div>
        </div>
      )}
      {result.molecule_info.num_bonds != null && (
        <div className="bg-[var(--color-surface-sunken)] rounded-lg p-2">
          <div className="text-lg font-bold text-[var(--color-text-primary)]">{result.molecule_info.num_bonds}</div>
          <div className="text-[10px] text-[var(--color-text-muted)] uppercase tracking-wider">Bonds</div>
        </div>
      )}
      {result.molecule_info.num_rings != null && (
        <div className="bg-[var(--color-surface-sunken)] rounded-lg p-2">
          <div className="text-lg font-bold text-[var(--color-text-primary)]">{result.molecule_info.num_rings}</div>
          <div className="text-[10px] text-[var(--color-text-muted)] uppercase tracking-wider">Rings</div>
        </div>
      )}
      {result.molecule_info.num_aromatic_rings != null && result.molecule_info.num_aromatic_rings > 0 && (
        <div className="bg-[var(--color-surface-sunken)] rounded-lg p-2">
          <div className="text-lg font-bold text-[var(--color-text-primary)]">{result.molecule_info.num_aromatic_rings}</div>
          <div className="text-[10px] text-[var(--color-text-muted)] uppercase tracking-wider">Arom.</div>
        </div>
      )}
      {result.molecule_info.molecular_weight && (
        <div className="bg-[var(--color-surface-sunken)] rounded-lg p-2">
          <div className="text-lg font-bold text-[var(--color-text-primary)]">{result.molecule_info.molecular_weight.toFixed(1)}</div>
          <div className="text-[10px] text-[var(--color-text-muted)] uppercase tracking-wider">MW</div>
        </div>
      )}
    </div>
    ...rest of detailed info section stays the same...
```

**Step 2: Update stereochemistry indicator (lines 1143-1166)**

Replace the regex-based E/Z detection with `moleculeInfo.hasEZStereo`:

```tsx
{moleculeInfo?.hasStereochemistry && (
  <motion.div
    initial={{ opacity: 0 }}
    animate={{ opacity: 1 }}
    className="mt-3 flex items-center justify-center gap-2"
  >
    <div className="flex items-center gap-1.5 text-xs bg-purple-500/10 text-purple-600 dark:text-purple-400 px-2 py-1 rounded-lg">
      <svg className="w-3 h-3" viewBox="0 0 24 24" fill="none" stroke="currentColor" strokeWidth="2">
        <path d="M12 2L2 7l10 5 10-5-10-5zM2 17l10 5 10-5M2 12l10 5 10-5" />
      </svg>
      <span>
        {moleculeInfo.numStereocenters > 0 && `${moleculeInfo.numStereocenters} stereocenter${moleculeInfo.numStereocenters > 1 ? 's' : ''}`}
        {moleculeInfo.numStereocenters > 0 && moleculeInfo.hasEZStereo && ' + '}
        {moleculeInfo.hasEZStereo && 'E/Z bonds'}
      </span>
    </div>
    {showCIP && (
      <span className="text-xs text-[var(--color-text-muted)]">
        (CIP labels shown)
      </span>
    )}
  </motion.div>
)}
```

**Step 3: Update pre-validation grid (lines 666-681)**

Update the pre-validation view to show 4 tiles (add aromatic rings when > 0):

```tsx
) : moleculeInfo ? (
  <>
    <div className="grid grid-cols-3 gap-3 text-center">
      <div className="bg-[var(--color-surface-sunken)] rounded-lg p-2">
        <div className="text-lg font-bold text-[var(--color-text-primary)]">{moleculeInfo.numAtoms}</div>
        <div className="text-[10px] text-[var(--color-text-muted)] uppercase tracking-wider">Atoms</div>
      </div>
      <div className="bg-[var(--color-surface-sunken)] rounded-lg p-2">
        <div className="text-lg font-bold text-[var(--color-text-primary)]">{moleculeInfo.numBonds}</div>
        <div className="text-[10px] text-[var(--color-text-muted)] uppercase tracking-wider">Bonds</div>
      </div>
      <div className="bg-[var(--color-surface-sunken)] rounded-lg p-2">
        <div className="text-lg font-bold text-[var(--color-text-primary)]">{moleculeInfo.numRings}</div>
        <div className="text-[10px] text-[var(--color-text-muted)] uppercase tracking-wider">Rings</div>
      </div>
    </div>
```

This stays the same — no changes needed for pre-validation grid.

**Step 4: Verify frontend compiles**

```bash
cd frontend
npx tsc --noEmit 2>&1 | head -30
```

Expected: No type errors.

**Step 5: Commit**

```bash
git add frontend/src/pages/SingleValidation.tsx
git commit -m "fix: use backend molecule_info for post-validation display"
```

---

### Task 7: Backend — Run full test suite for regressions

**Step 1: Run all backend tests**

```bash
cd backend
python -m pytest tests/ -v --timeout=60 -x 2>&1 | tail -40
```

Expected: ALL PASS. The schema change adds optional fields, so existing tests should not break.

**Step 2: Specifically run validation tests**

```bash
cd backend
python -m pytest tests/test_api/test_validate_endpoint.py -v
```

Expected: ALL PASS including the new `TestMoleculeInfoEnrichedFields` tests.

---

### Task 8: Frontend — Build verification

**Step 1: Run TypeScript type check**

```bash
cd frontend
npx tsc --noEmit
```

Expected: No errors.

**Step 2: Run frontend build**

```bash
cd frontend
npm run build 2>&1 | tail -20
```

Expected: Build succeeds.

**Step 3: Run frontend tests**

```bash
cd frontend
npm test -- --run 2>&1 | tail -20
```

Expected: Tests pass (or pre-existing failures only — no new failures from our changes).

---

### Task 9: Final commit and verification summary

**Step 1: Review all changes**

```bash
git diff development --stat
```

Expected files changed:
- `backend/app/schemas/validation.py` — 5 new fields
- `backend/app/api/routes/validation.py` — enriched extract_molecule_info()
- `backend/tests/test_api/test_validate_endpoint.py` — new test class
- `frontend/src/types/validation.ts` — 5 new fields
- `frontend/src/hooks/useRDKit.ts` — added get_descriptors/get_json to interface
- `frontend/src/hooks/useMoleculeInfo.ts` — full rewrite (regex → get_descriptors)
- `frontend/src/pages/SingleValidation.tsx` — use backend data, remove regex E/Z

**Step 2: Verify the original bug is fixed**

```bash
cd backend
python -m pytest tests/test_api/test_validate_endpoint.py::TestMoleculeInfoEnrichedFields::test_enriched_fields_isotope_no_false_rings -v
```

Expected: PASS — `num_rings == 0` for the Ni-63 complex.
