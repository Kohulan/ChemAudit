import { useState, useEffect, useRef } from 'react';
import { getRDKit, RDKitMol } from './useRDKit';

interface UseMoleculeResult {
  svg: string | null;
  isLoading: boolean;
  error: string | null;
  isValid: boolean;
}

interface UseMoleculeOptions {
  width?: number;
  height?: number;
  highlightAtoms?: number[];
}

export function useMolecule(
  smiles: string | null,
  options: UseMoleculeOptions = {}
): UseMoleculeResult {
  const { width = 300, height = 200, highlightAtoms = [] } = options;

  const [svg, setSvg] = useState<string | null>(null);
  const [isLoading, setIsLoading] = useState(false);
  const [error, setError] = useState<string | null>(null);

  // CRITICAL: Track mol reference for cleanup (PITFALL 2)
  const molRef = useRef<RDKitMol | null>(null);

  useEffect(() => {
    // If no SMILES, clear state
    if (!smiles || smiles.trim() === '') {
      setSvg(null);
      setError(null);
      return;
    }

    let cancelled = false;
    setIsLoading(true);
    setError(null);

    getRDKit()
      .then((rdkit) => {
        if (cancelled) return;

        // CRITICAL: Delete previous molecule to prevent memory leak
        if (molRef.current) {
          molRef.current.delete();
          molRef.current = null;
        }

        try {
          const mol = rdkit.get_mol(smiles);

          if (!mol) {
            setError('Invalid molecule');
            setSvg(null);
            setIsLoading(false);
            return;
          }

          // Store reference for cleanup
          molRef.current = mol;

          // Generate SVG
          let svgContent: string;
          if (highlightAtoms.length > 0) {
            const highlightDetails = JSON.stringify({
              atoms: highlightAtoms,
              highlightColour: [1, 0.5, 0.5]  // Red-ish highlight
            });
            svgContent = mol.get_svg_with_highlights(highlightDetails);
          } else {
            svgContent = mol.get_svg(width, height);
          }

          setSvg(svgContent);
          setError(null);
        } catch (e) {
          setError(e instanceof Error ? e.message : 'Unknown error');
          setSvg(null);
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

    // CRITICAL: Cleanup on unmount or smiles change
    return () => {
      cancelled = true;
      if (molRef.current) {
        molRef.current.delete();
        molRef.current = null;
      }
    };
  }, [smiles, width, height, JSON.stringify(highlightAtoms)]);

  return {
    svg,
    isLoading,
    error,
    isValid: svg !== null && error === null
  };
}
