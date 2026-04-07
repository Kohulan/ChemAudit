/**
 * RegistrationTab
 *
 * Registration hash tab showing:
 * - Summary card with RDKit version, hash version, unique count, collision count
 * - Collision groups section with expandable RegHashCollisionGroup cards
 * - No-collisions success state
 * - Computing and error states
 */

// React components use JSX transform (no explicit React import needed)
import { Loader2, AlertTriangle } from 'lucide-react';
import { RegHashCollisionGroup } from './RegHashCollisionGroup';
import type { BatchAnalyticsResponse, RegistrationHashResult } from '../../types/analytics';
import type { BatchResult } from '../../types/batch';

interface RegistrationTabProps {
  /** Full analytics response */
  analyticsData: BatchAnalyticsResponse | null;
  /** Batch results for molecule SMILES lookup */
  results: BatchResult[];
}

export function RegistrationTab({ analyticsData, results }: RegistrationTabProps) {
  const registrationResult: RegistrationHashResult | undefined = analyticsData?.registration;
  const isComputing = analyticsData?.status?.registration?.status === 'computing';
  const errorMessage = analyticsData?.status?.registration?.error;

  // Error state
  if (errorMessage) {
    return (
      <div className="rounded-2xl p-5 bg-amber-500/5 border border-amber-500/20">
        <div className="flex items-center gap-3">
          <AlertTriangle className="w-5 h-5 text-amber-500 flex-shrink-0" />
          <p className="text-sm text-amber-700 dark:text-amber-300">{errorMessage}</p>
        </div>
      </div>
    );
  }

  // Computing state
  if (isComputing && !registrationResult) {
    return (
      <div className="space-y-3">
        <div className="flex items-center gap-3 p-5 rounded-2xl bg-[var(--color-surface-sunken)]/50 border border-[var(--color-border)]">
          <Loader2 className="w-5 h-5 text-[var(--color-primary)] animate-spin" />
          <div>
            <p className="text-sm text-[var(--color-text-primary)]">
              Computing registration hashes...
            </p>
            <p className="text-xs text-[var(--color-text-muted)] mt-0.5">
              Registration hashes are computed automatically after batch processing.
            </p>
          </div>
        </div>
      </div>
    );
  }

  // No data yet
  if (!registrationResult) {
    return (
      <div className="text-center py-12 space-y-2">
        <p className="text-sm font-medium text-[var(--color-text-primary)]">
          No registration hash data
        </p>
        <p className="text-xs text-[var(--color-text-muted)]">
          Registration hashes are computed automatically after batch processing.
        </p>
      </div>
    );
  }

  const hasCollisions = registrationResult.collision_groups.length > 0;

  return (
    <div className="space-y-6">
      {/* Summary card */}
      <div className="rounded-2xl p-5 bg-gradient-to-br from-[var(--color-surface-elevated)] to-[var(--color-surface)] border border-[var(--color-border)]">
        <h4 className="text-base font-semibold text-[var(--color-text-primary)] font-display mb-4">
          Registration Hashes
        </h4>
        <div className="grid grid-cols-2 md:grid-cols-4 gap-4">
          <div>
            <span className="text-xs text-[var(--color-text-muted)] uppercase font-medium">
              RDKit Version
            </span>
            <p className="text-sm text-[var(--color-text-primary)] mt-0.5">
              {registrationResult.rdkit_version}
            </p>
          </div>
          <div>
            <span className="text-xs text-[var(--color-text-muted)] uppercase font-medium">
              Tautomer Hash Version
            </span>
            <p className="text-sm text-[var(--color-text-primary)] mt-0.5">v2</p>
          </div>
          <div>
            <span className="text-xs text-[var(--color-text-muted)] uppercase font-medium">
              Unique Hashes
            </span>
            <p className="text-sm font-semibold text-[var(--color-text-primary)] mt-0.5">
              {registrationResult.unique_count}/{registrationResult.total_count}
            </p>
          </div>
          <div>
            <span className="text-xs text-[var(--color-text-muted)] uppercase font-medium">
              Collision Groups
            </span>
            <p className="text-sm font-semibold text-[var(--color-text-primary)] mt-0.5">
              {registrationResult.collision_groups.length}
            </p>
          </div>
        </div>
      </div>

      {/* Collision groups section */}
      {hasCollisions ? (
        <div className="space-y-3">
          <h4 className="text-base font-semibold text-[var(--color-text-primary)] font-display">
            Hash Collisions
          </h4>
          <p className="text-sm text-[var(--color-text-muted)]">
            These molecules produce the same RegistrationHash despite different SMILES representations.
          </p>
          <div className="space-y-2">
            {registrationResult.collision_groups.map((group) => (
              <RegHashCollisionGroup
                key={group.hash}
                group={group}
                results={results}
              />
            ))}
          </div>
        </div>
      ) : (
        /* No collisions success state */
        <div className="rounded-2xl p-5 bg-green-100 text-green-700 dark:bg-green-900/30 dark:text-green-400 border border-green-500/20">
          <p className="text-sm font-medium">
            No hash collisions detected. All molecules have unique registration hashes.
          </p>
        </div>
      )}
    </div>
  );
}
