import { Shield } from 'lucide-react';
import { Link } from 'react-router-dom';
import { Skeleton } from '../ui/Skeleton';

interface SafetyBadgeProps {
  /** Molecule SMILES — passed as /safety?smiles=... query param */
  smiles: string;
  /** Total structural alert count (optional) */
  totalAlerts?: number;
  /** Count of non-pass safety flags (optional) */
  safetyFlags?: number;
  /** True while loading safety data */
  isLoading?: boolean;
}

/**
 * Compact cross-link badge for SingleValidation and CompoundProfiler pages.
 *
 * Renders a shield icon, a summary of alert/safety flag counts, and a link
 * to the full safety report at /safety?smiles=... for SPA navigation.
 *
 * Per D-02 (SV cross-link) and D-03 (Profiler cross-link).
 * The badge does NOT call the safety API — counts are passed in as optional props.
 */
export function SafetyBadge({
  smiles,
  totalAlerts,
  safetyFlags,
  isLoading = false,
}: SafetyBadgeProps) {
  const hasAlerts = (totalAlerts ?? 0) > 0 || (safetyFlags ?? 0) > 0;

  return (
    <div className="flex items-center gap-2 p-3 rounded-lg border border-[var(--color-border)] bg-[var(--color-surface-sunken)]">
      {/* Shield icon */}
      <Shield className="w-4 h-4 text-chem-primary-600 shrink-0" />

      {/* Summary text */}
      <div className="flex-1 min-w-0">
        {isLoading ? (
          <Skeleton variant="rounded" height={16} width={180} />
        ) : hasAlerts ? (
          <span className="text-sm text-text-secondary">
            {totalAlerts ?? 0} alerts &middot; {safetyFlags ?? 0} safety flags
          </span>
        ) : (
          <span className="text-sm text-text-muted">No alerts detected</span>
        )}
      </div>

      {/* Link to full safety report */}
      <Link
        to={`/safety?smiles=${encodeURIComponent(smiles)}`}
        className="text-sm text-chem-primary-600 hover:underline shrink-0"
      >
        View full safety report &rarr;
      </Link>
    </div>
  );
}
