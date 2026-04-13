import { Microscope } from 'lucide-react';
import { Link } from 'react-router-dom';

interface ProfilerBadgeProps {
  smiles: string;
}

/**
 * Cross-link badge from SingleValidation to CompoundProfiler page.
 * Follows SafetyBadge/DiagnosticsBadge pattern (D-24).
 */
export function ProfilerBadge({ smiles }: ProfilerBadgeProps) {
  return (
    <div className="flex items-center gap-2 p-3 rounded-lg border border-[var(--color-border)] bg-[var(--color-surface-sunken)]">
      <Microscope className="w-4 h-4 text-chem-primary-600 shrink-0" />
      <div className="flex-1 min-w-0">
        <span className="text-sm text-text-secondary">
          PFI, #stars, MPO, ligand efficiency, SA comparison
        </span>
      </div>
      <Link
        to={`/profiler?smiles=${encodeURIComponent(smiles)}`}
        className="text-sm text-chem-primary-600 hover:underline shrink-0"
      >
        Open full profiler &rarr;
      </Link>
    </div>
  );
}
