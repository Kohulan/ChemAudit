import { motion } from 'framer-motion';
import { AlertTriangle, ShieldX } from 'lucide-react';
import { useSafety } from '../hooks/useSafety';
import { SafetyInput } from '../components/safety/SafetyInput';
import { SafetySummaryStrip } from '../components/safety/SafetySummaryStrip';
import { ClayCard } from '../components/ui/ClayCard';
import { Badge } from '../components/ui/Badge';
import { Skeleton } from '../components/ui/Skeleton';
import { ClayButton } from '../components/ui/ClayButton';

/**
 * Empty state shown before any molecule is screened.
 */
function EmptyState() {
  return (
    <motion.div
      initial={{ opacity: 0, y: 8 }}
      animate={{ opacity: 1, y: 0 }}
      transition={{ duration: 0.3 }}
    >
      <ClayCard variant="flat" size="md" className="text-center py-12">
        <div className="flex flex-col items-center gap-3">
          <div className="w-12 h-12 rounded-2xl bg-[var(--color-surface-sunken)] flex items-center justify-center">
            <ShieldX className="w-6 h-6 text-[var(--color-text-muted)]" />
          </div>
          <h2 className="text-base font-semibold font-display text-[var(--color-text-primary)]">
            No molecule screened yet
          </h2>
          <p className="text-sm text-[var(--color-text-secondary)] max-w-sm">
            Enter a SMILES, InChI, CAS number, or ChEMBL ID above to run a full safety screening.
          </p>
        </div>
      </ClayCard>
    </motion.div>
  );
}

/**
 * Error card for alert screening failures.
 */
function AlertErrorCard({ error, onRetry }: { error: string; onRetry: () => void }) {
  return (
    <motion.div
      initial={{ opacity: 0, y: 8 }}
      animate={{ opacity: 1, y: 0 }}
      transition={{ duration: 0.25 }}
    >
      <ClayCard variant="default" size="sm" className="border border-status-error/20">
        <div className="flex items-start gap-3">
          <AlertTriangle className="w-4 h-4 text-status-error mt-0.5 shrink-0" />
          <div className="flex-1 min-w-0">
            <p className="text-sm font-semibold text-[var(--color-text-primary)] font-display">
              Alert screening failed
            </p>
            <p className="text-sm text-[var(--color-text-secondary)] mt-0.5">{error}</p>
          </div>
          <ClayButton variant="ghost" size="sm" onClick={onRetry} className="shrink-0">
            Retry
          </ClayButton>
        </div>
      </ClayCard>
    </motion.div>
  );
}

/**
 * Error card for safety assessment failures.
 */
function SafetyErrorCard({ error, onRetry }: { error: string; onRetry: () => void }) {
  return (
    <motion.div
      initial={{ opacity: 0, y: 8 }}
      animate={{ opacity: 1, y: 0 }}
      transition={{ duration: 0.25 }}
    >
      <ClayCard variant="default" size="sm" className="border border-status-error/20">
        <div className="flex items-start gap-3">
          <AlertTriangle className="w-4 h-4 text-status-error mt-0.5 shrink-0" />
          <div className="flex-1 min-w-0">
            <p className="text-sm font-semibold text-[var(--color-text-primary)] font-display">
              Safety assessment failed
            </p>
            <p className="text-sm text-[var(--color-text-secondary)] mt-0.5">{error}</p>
          </div>
          <ClayButton variant="ghost" size="sm" onClick={onRetry} className="shrink-0">
            Retry
          </ClayButton>
        </div>
      </ClayCard>
    </motion.div>
  );
}

/**
 * Loading skeleton for the summary strip area.
 */
function SummaryStripSkeleton() {
  return (
    <div className="flex items-center justify-between gap-4 px-5 py-3 rounded-2xl border border-[var(--color-border)] mb-6">
      <Skeleton variant="rounded" height={24} width={240} />
      <div className="flex gap-4">
        {[0, 1, 2, 3].map(i => (
          <Skeleton key={i} variant="rounded" height={16} width={48} />
        ))}
      </div>
      <Skeleton variant="rounded" height={16} width={180} />
    </div>
  );
}

/**
 * Lilly Demerits coming-soon placeholder (D-16).
 * Greyed out with a "Planned" badge to indicate future availability.
 */
function LillyDemeritsPlaceholder() {
  return (
    <div className="opacity-50 pointer-events-none" aria-hidden="true">
      <ClayCard variant="flat" size="md">
        <div className="flex items-center justify-between">
          <div>
            <h3 className="text-base font-semibold font-display text-[var(--color-text-muted)]">
              Lilly Demerits
            </h3>
            <p className="text-sm text-[var(--color-text-muted)] mt-0.5">
              Coming soon — 275-rule graduated demerit scoring
            </p>
          </div>
          <Badge variant="default">Planned</Badge>
        </div>
      </ClayCard>
    </div>
  );
}

/**
 * Safety page — /safety
 *
 * Full-page safety screening layout (D-07):
 * 1. Input (SafetyInput)
 * 2. Summary strip (SafetySummaryStrip) — sticky after results
 * 3. Alert Dashboard slot — Plan 05
 * 4. Safety Flags Grid slot — Plan 06
 * 5. Complexity Radar slot — Plan 06
 * 6. Lilly Demerits placeholder (D-16)
 */
export function Safety() {
  const {
    alertResult,
    safetyResult,
    isLoading,
    alertError,
    safetyError,
    screenMolecule,
  } = useSafety();

  const hasResults = !!(alertResult || safetyResult);

  // Derive SafetySummaryStrip props from both results
  const totalAlerts = alertResult?.total_raw ?? 0;
  const concernGroupCount = Object.keys(alertResult?.concern_groups ?? {}).length;
  const hasCritical = alertResult?.has_critical ?? false;

  // Traffic-light statuses derived from safetyResult
  const cypStatus = 'default' as const; // CYP is always "default" (atom-level, no binary flag)
  const hergStatus: 'success' | 'warning' | 'error' =
    safetyResult?.herg.herg_risk === 'high'
      ? 'error'
      : safetyResult?.herg.herg_risk === 'moderate'
        ? 'warning'
        : 'success';
  const bro5Status: 'success' | 'error' | 'default' =
    !safetyResult?.bro5.applicable
      ? 'default'
      : safetyResult.bro5.passed
        ? 'success'
        : 'error';
  const reosStatus: 'success' | 'warning' | 'error' =
    (safetyResult?.reos.n_violations ?? 0) === 0
      ? 'success'
      : (safetyResult?.reos.n_violations ?? 0) === 1
        ? 'warning'
        : 'error';
  const complexityOutliers = safetyResult?.complexity.n_outliers ?? 0;

  return (
    <div className="max-w-[1200px] mx-auto px-4 pt-16 pb-16">
      {/* Page heading */}
      <div className="mb-6">
        <h1 className="text-2xl font-semibold font-display text-[var(--color-text-primary)]">
          Structural Alerts &amp; Safety
        </h1>
        <p className="text-sm text-[var(--color-text-secondary)] mt-1">
          Screen any molecule against expanded alert libraries and rule-based safety assessments
        </p>
      </div>

      {/* Input */}
      <div className="mb-6">
        <SafetyInput onScreen={screenMolecule} isLoading={isLoading} />
      </div>

      {/* Loading skeleton for summary strip */}
      {isLoading && <SummaryStripSkeleton />}

      {/* Summary strip — shown after any result is available */}
      {!isLoading && hasResults && (
        <SafetySummaryStrip
          totalAlerts={totalAlerts}
          concernGroupCount={concernGroupCount}
          hasCritical={hasCritical}
          cypStatus={cypStatus}
          hergStatus={hergStatus}
          bro5Status={bro5Status}
          reosStatus={reosStatus}
          complexityOutliers={complexityOutliers}
        />
      )}

      {/* Empty state */}
      {!isLoading && !hasResults && !alertError && !safetyError && <EmptyState />}

      {/* Alert Dashboard slot — populated in Plan 05 */}
      {alertResult && (
        <div id="alert-dashboard-slot" className="mb-6" />
      )}

      {/* Alert error card */}
      {alertError && (
        <div className="mb-6">
          <AlertErrorCard
            error={alertError}
            onRetry={() => {
              const lastSmiles = safetyResult?.molecule_info.smiles ?? '';
              if (lastSmiles) screenMolecule(lastSmiles);
            }}
          />
        </div>
      )}

      {/* Safety Flags Grid slot — populated in Plan 06 */}
      {safetyResult && (
        <div id="safety-flags-slot" className="mb-6" />
      )}

      {/* Safety error card */}
      {safetyError && (
        <div className="mb-6">
          <SafetyErrorCard
            error={safetyError}
            onRetry={() => {
              const lastSmiles = alertResult?.molecule_info.smiles ?? safetyResult?.molecule_info.smiles ?? '';
              if (lastSmiles) screenMolecule(lastSmiles);
            }}
          />
        </div>
      )}

      {/* Complexity Radar slot — populated in Plan 06 */}
      {safetyResult && (
        <div id="complexity-radar-slot" className="mb-6" />
      )}

      {/* Lilly Demerits coming-soon placeholder (D-16) */}
      <div className="mt-8">
        <LillyDemeritsPlaceholder />
      </div>
    </div>
  );
}
