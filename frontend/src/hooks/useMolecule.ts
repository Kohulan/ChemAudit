import { useState, useEffect, useRef } from 'react';
import { getRDKit, RDKitMol } from './useRDKit';

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

          // Generate SVG with optional CIP stereochemistry labels
          let svgContent: string;

          // Build drawing options for RDKit.js
          const drawingOptions: Record<string, unknown> = {
            width,
            height,
          };

          // Add stereo annotation if showCIP is enabled
          if (showCIP) {
            drawingOptions.addStereoAnnotation = true;
          }

          if (highlightAtoms.length > 0) {
            const highlightDetails = JSON.stringify({
              atoms: highlightAtoms,
              highlightColour: [1, 0.5, 0.5],
              ...drawingOptions
            });
            svgContent = mol.get_svg_with_highlights(highlightDetails);
          } else if (showCIP) {
            // Try RDKit.js get_svg_with_highlights for CIP labels
            try {
              svgContent = mol.get_svg_with_highlights(JSON.stringify(drawingOptions));
            } catch {
              // Fallback to basic SVG if options not supported
            }
            svgContent ??= mol.get_svg(width, height);
          } else {
            svgContent = mol.get_svg(width, height);
          }

          // Post-process SVG to add padding by expanding viewBox
          // This prevents molecules from being cut off at edges
          svgContent = expandSvgViewBox(svgContent, 20);

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
  }, [smiles, width, height, JSON.stringify(highlightAtoms), showCIP]);

  return {
    svg,
    isLoading,
    error,
    isValid: svg !== null && error === null
  };
}
