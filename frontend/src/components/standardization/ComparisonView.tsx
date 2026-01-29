/**
 * ComparisonView Component
 *
 * Shows side-by-side comparison of original and standardized molecules.
 */
import { MoleculeViewer } from '../molecules/MoleculeViewer';

interface ComparisonViewProps {
  originalSmiles: string;
  standardizedSmiles: string | null;
  className?: string;
}

export function ComparisonView({
  originalSmiles,
  standardizedSmiles,
  className = ''
}: ComparisonViewProps) {
  const isIdentical = originalSmiles === standardizedSmiles;

  // If identical, show single viewer with badge
  if (isIdentical || !standardizedSmiles) {
    return (
      <div className={`flex flex-col items-center ${className}`}>
        <div className="relative">
          <MoleculeViewer
            smiles={originalSmiles}
            width={350}
            height={250}
            className="border-2 border-gray-200 rounded-lg"
          />
          {isIdentical && standardizedSmiles && (
            <div className="absolute top-2 right-2 bg-yellow-100 dark:bg-yellow-900/30 text-amber-800 dark:text-yellow-400 text-xs font-medium px-2 py-1 rounded">
              No changes needed
            </div>
          )}
        </div>
        <p className="text-sm text-gray-500 mt-2">
          {isIdentical ? 'Structure already standardized' : 'Original Structure'}
        </p>
      </div>
    );
  }

  // Side-by-side comparison
  return (
    <div className={`flex flex-col lg:flex-row gap-6 ${className}`}>
      {/* Original */}
      <div className="flex-1 flex flex-col items-center">
        <div className="border-2 border-gray-300 rounded-lg p-2 bg-gray-50">
          <MoleculeViewer
            smiles={originalSmiles}
            width={280}
            height={200}
          />
        </div>
        <div className="mt-2 text-center">
          <span className="text-sm font-medium text-gray-700">Original</span>
          <p className="text-xs text-gray-500 mt-1 max-w-[280px] truncate font-mono">
            {originalSmiles}
          </p>
        </div>
      </div>

      {/* Arrow */}
      <div className="flex items-center justify-center">
        <div className="hidden lg:block text-3xl text-blue-500">→</div>
        <div className="lg:hidden text-3xl text-blue-500">↓</div>
      </div>

      {/* Standardized */}
      <div className="flex-1 flex flex-col items-center">
        <div className="border-2 border-blue-400 rounded-lg p-2 bg-blue-50">
          <MoleculeViewer
            smiles={standardizedSmiles}
            width={280}
            height={200}
          />
        </div>
        <div className="mt-2 text-center">
          <span className="text-sm font-medium text-blue-700">Standardized</span>
          <p className="text-xs text-gray-500 mt-1 max-w-[280px] truncate font-mono">
            {standardizedSmiles}
          </p>
        </div>
      </div>
    </div>
  );
}
