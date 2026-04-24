import type { ReactElement } from 'react';
import { useMolecule } from '../../hooks/useMolecule';
import { sanitizeSvg } from '../../lib/sanitize';

interface MoleculeViewerProps {
  smiles: string | null;
  highlightAtoms?: number[];
  highlightBonds?: number[];
  /** RGB triple in 0..1 range for highlight fill. Defaults to orange. */
  highlightColor?: [number, number, number];
  /** When false, suppress atom-index number labels around highlighted atoms. */
  showAtomLabels?: boolean;
  width?: number;
  height?: number;
  className?: string;
  showCIP?: boolean;
  /**
   * How the SVG sizes itself inside its container.
   * - 'width'   (default): SVG fills container width, height computed from viewBox.
   *                        Can overflow if the parent has a fixed height.
   * - 'contain': SVG fills container both width and height, letterboxed via
   *              preserveAspectRatio (xMidYMid meet). Use when the parent
   *              is a fixed-size box.
   */
  fit?: 'width' | 'contain';
}

const PLACEHOLDER_BASE =
  'flex items-center justify-center rounded-lg w-full';

/**
 * Strip the fixed width/height attributes and the opaque background rect that
 * RDKit inserts. The SVG keeps its viewBox so the browser can compute the
 * correct aspect ratio; actual sizing is handled via CSS on the container.
 */
function makeSvgResponsive(svgStr: string): string {
  return svgStr
    .replace(/(<svg[^>]*)\s+width=["']\d+(?:px)?["']/, '$1')
    .replace(/(<svg[^>]*)\s+height=["']\d+(?:px)?["']/, '$1')
    .replace(/<rect[^>]*style=['"]opacity:\s*1\.0;fill:#FFFFFF[^"']*['"][^>]*\/>/, '');
}

export function MoleculeViewer({
  smiles,
  highlightAtoms = [],
  highlightBonds = [],
  highlightColor,
  showAtomLabels = true,
  width = 300,
  height = 200,
  className = '',
  showCIP = false,
  fit = 'width',
}: MoleculeViewerProps): ReactElement {
  const { svg, isLoading, error, isValid } = useMolecule(smiles, {
    width,
    height,
    highlightAtoms,
    highlightBonds,
    highlightColor,
    showAtomLabels,
    showCIP,
  });

  // Placeholder sizing matches the expected render height so it doesn't
  // overflow small containers (e.g. 60px cluster thumbnails) or collapse
  // in large containers (e.g. 200px SingleValidation viewer).
  const placeholderStyle = { minHeight: height };

  if (!smiles) {
    return (
      <div className={`${PLACEHOLDER_BASE} bg-gray-100 dark:bg-gray-800 ${className}`} style={placeholderStyle}>
        <span className="text-gray-400 text-sm">Enter a molecule</span>
      </div>
    );
  }

  if (isLoading) {
    return (
      <div className={`${PLACEHOLDER_BASE} bg-gray-100 dark:bg-gray-800 ${className}`} style={placeholderStyle}>
        <div className="animate-spin rounded-full h-8 w-8 border-b-2 border-blue-600"></div>
      </div>
    );
  }

  if (error || !isValid) {
    return (
      <div className={`${PLACEHOLDER_BASE} bg-red-50 dark:bg-red-900/20 border border-red-200 dark:border-red-800 ${className}`} style={placeholderStyle}>
        <span className="text-red-500 text-sm px-4 text-center">{error || 'Invalid molecule'}</span>
      </div>
    );
  }

  const responsiveSvg = makeSvgResponsive(sanitizeSvg(svg));

  // 'contain' letterboxes the SVG inside the parent's box (both w/h 100%,
  // SVG preserveAspectRatio handles centering). 'width' keeps the original
  // width-driven behavior where height scales to aspect ratio — fine for
  // containers with no fixed height, but overflows fixed-height parents.
  const fitClasses =
    fit === 'contain'
      ? '[&>svg]:block [&>svg]:w-full [&>svg]:h-full'
      : '[&>svg]:block [&>svg]:w-full [&>svg]:h-auto';

  return (
    <div
      className={`rounded-lg w-full ${fitClasses} ${className}`}
      dangerouslySetInnerHTML={{ __html: responsiveSvg }}
    />
  );
}
