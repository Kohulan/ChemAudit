import { useState, useEffect, useRef } from 'react';
import { getRDKit, RDKitMol } from './useRDKit';

/**
 * Add atom index labels to highlighted atoms in SVG.
 * RDKit draws highlight ellipses in the same order as the atoms array passed to get_svg_with_highlights.
 * We find these ellipses and add text labels showing the atom indices.
 */
function addAtomLabelsToSvg(svgContent: string, atomIndices: number[]): string {
  if (atomIndices.length === 0) return svgContent;

  // RDKit draws atom highlights as ellipses at the beginning of the SVG.
  // The ellipses are drawn in the same order as atomIndices array.
  const ellipsePattern = /<ellipse[^>]*cx=["']([^"']+)["'][^>]*cy=["']([^"']+)["'][^>]*\/>/gi;
  const positions: Array<{ cx: number; cy: number; atomIdx: number }> = [];

  let matchResult;
  let ellipseCount = 0;
  while ((matchResult = ellipsePattern.exec(svgContent)) !== null) {
    // Only process ellipses up to the number of highlighted atoms
    if (ellipseCount < atomIndices.length) {
      positions.push({
        cx: parseFloat(matchResult[1]),
        cy: parseFloat(matchResult[2]),
        atomIdx: atomIndices[ellipseCount],
      });
    }
    ellipseCount++;
  }

  // If no positions found, return original
  if (positions.length === 0) return svgContent;

  // Create text labels positioned near the highlighted atoms
  // Offset slightly up and right, with white stroke for readability
  const labelElements = positions
    .map(({ cx, cy, atomIdx }) => {
      const x = cx + 12;
      const y = cy - 12;
      return `<text x="${x}" y="${y}" font-family="Arial,sans-serif" font-size="12" font-weight="bold" fill="#B45309" stroke="white" stroke-width="3" paint-order="stroke">${atomIdx}</text>`;
    })
    .join('');

  // Insert labels before closing svg tag
  return svgContent.replace(/<\/svg>/, `${labelElements}</svg>`);
}

/**
 * Expand SVG viewBox by adding padding to prevent molecule cutoff.
 * Parses the viewBox, expands it, and shifts content to center.
 */
function expandSvgViewBox(svg: string, padding: number): string {
  // Match viewBox attribute: viewBox="minX minY width height"
  const viewBoxMatch = svg.match(/viewBox=["']([^"']+)["']/);
  if (!viewBoxMatch) {
    // No viewBox found, try to add one based on width/height attributes
    const widthMatch = svg.match(/width=["'](\d+)(?:px)?["']/);
    const heightMatch = svg.match(/height=["'](\d+)(?:px)?["']/);
    if (widthMatch && heightMatch) {
      const w = parseInt(widthMatch[1]);
      const h = parseInt(heightMatch[1]);
      // Add viewBox with padding
      return svg.replace(
        /<svg([^>]*)>/,
        `<svg$1 viewBox="${-padding} ${-padding} ${w + padding * 2} ${h + padding * 2}">`
      );
    }
    return svg;
  }

  const [minX, minY, vbWidth, vbHeight] = viewBoxMatch[1].split(/\s+/).map(Number);

  // Expand viewBox by padding on all sides
  const newMinX = minX - padding;
  const newMinY = minY - padding;
  const newWidth = vbWidth + padding * 2;
  const newHeight = vbHeight + padding * 2;

  return svg.replace(
    /viewBox=["'][^"']+["']/,
    `viewBox="${newMinX} ${newMinY} ${newWidth} ${newHeight}"`
  );
}

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
  showCIP?: boolean;  // Show R/S and E/Z stereochemistry labels
}

export function useMolecule(
  smiles: string | null,
  options: UseMoleculeOptions = {}
): UseMoleculeResult {
  const { width = 300, height = 200, highlightAtoms = [], showCIP = false } = options;

  const [svg, setSvg] = useState<string | null>(null);
  const [isLoading, setIsLoading] = useState(false);
  const [error, setError] = useState<string | null>(null);

  // CRITICAL: Track mol reference for cleanup (PITFALL 2)
  const molRef = useRef<RDKitMol | null>(null);

  // Counter that increments when a new mol is created, triggering SVG re-render
  const [molVersion, setMolVersion] = useState(0);

  // Serialize highlightAtoms for stable dependency comparison
  const highlightAtomsKey = JSON.stringify(highlightAtoms);

  // Effect 1: Create/destroy molecule when SMILES changes (async, shows loading)
  useEffect(() => {
    if (!smiles || smiles.trim() === '') {
      if (molRef.current) {
        molRef.current.delete();
        molRef.current = null;
      }
      setSvg(null);
      setError(null);
      setMolVersion(0);
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
          // Signal that a new mol is ready for SVG rendering
          setMolVersion((v) => v + 1);
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
  }, [smiles]);

  // Effect 2: Render SVG when mol is ready or rendering params change (sync, no loading flash)
  useEffect(() => {
    const mol = molRef.current;
    if (!mol || molVersion === 0) return;

    try {
      let svgContent: string;

      const drawingOptions: Record<string, unknown> = { width, height };
      if (showCIP) {
        drawingOptions.addStereoAnnotation = true;
      }

      if (highlightAtoms.length > 0) {
        const highlightDetails = JSON.stringify({
          atoms: highlightAtoms,
          highlightColour: [1, 0.3, 0],  // Bright orange for better visibility
          highlightRadius: 0.4,          // Slightly larger highlight radius
          ...drawingOptions,
        });

        // Render with highlights
        svgContent = mol.get_svg_with_highlights(highlightDetails);

        // Add atom index labels to highlighted atoms via post-processing
        svgContent = addAtomLabelsToSvg(svgContent, highlightAtoms);
      } else if (showCIP) {
        try {
          svgContent = mol.get_svg_with_highlights(JSON.stringify(drawingOptions));
        } catch {
          // Fallback to basic SVG if options not supported
        }
        svgContent ??= mol.get_svg(width, height);
      } else {
        svgContent = mol.get_svg(width, height);
      }

      svgContent = expandSvgViewBox(svgContent, 20);

      setSvg(svgContent);
      setError(null);
    } catch (e) {
      setError(e instanceof Error ? e.message : 'Unknown error');
      setSvg(null);
    }
    // Use highlightAtomsKey (stable string) instead of highlightAtoms (unstable reference)
    // eslint-disable-next-line react-hooks/exhaustive-deps
  }, [molVersion, width, height, highlightAtomsKey, showCIP]);

  return {
    svg,
    isLoading,
    error,
    isValid: svg !== null && error === null
  };
}
