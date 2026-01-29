import { useMolecule } from '../../hooks/useMolecule';

interface MoleculeViewerProps {
  smiles: string | null;
  highlightAtoms?: number[];
  width?: number;
  height?: number;
  className?: string;
  showCIP?: boolean;  // Show R/S and E/Z stereochemistry labels
}

export function MoleculeViewer({
  smiles,
  highlightAtoms = [],
  width = 300,
  height = 200,
  className = '',
  showCIP = false
}: MoleculeViewerProps) {
  const { svg, isLoading, error, isValid } = useMolecule(smiles, {
    width,
    height,
    highlightAtoms,
    showCIP
  });

  if (!smiles) {
    return (
      <div className={`flex items-center justify-center bg-gray-100 rounded-lg ${className}`}
           style={{ width, height }}>
        <span className="text-gray-400 text-sm">Enter a molecule</span>
      </div>
    );
  }

  if (isLoading) {
    return (
      <div className={`flex items-center justify-center bg-gray-100 rounded-lg ${className}`}
           style={{ width, height }}>
        <div className="animate-spin rounded-full h-8 w-8 border-b-2 border-blue-600"></div>
      </div>
    );
  }

  if (error || !isValid) {
    return (
      <div className={`flex items-center justify-center bg-red-50 rounded-lg border border-red-200 ${className}`}
           style={{ width, height }}>
        <span className="text-red-500 text-sm px-4 text-center">{error || 'Invalid molecule'}</span>
      </div>
    );
  }

  return (
    <div
      className={`bg-white dark:bg-gray-900 rounded-lg flex items-center justify-center ${className}`}
      style={{ width, height }}
      dangerouslySetInnerHTML={{ __html: svg || '' }}
    />
  );
}
