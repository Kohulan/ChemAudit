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
