import { motion } from 'framer-motion';
import { AlertTriangle } from 'lucide-react';
import { ClayCard } from '../ui/ClayCard';
import { ClayButton } from '../ui/ClayButton';
import { Badge } from '../ui/Badge';
import { Skeleton } from '../ui/Skeleton';
import { SMILESAnnotated } from './SMILESAnnotated';
import { FixChip } from './FixChip';
import type { SMILESDiagnosticsResponse } from '../../types/diagnostics';

const MAX_FIX_CHIPS = 3;

interface SMILESDiagnosticsPanelProps {
  /** Result from the SMILES diagnostics API (DIAG-01). */
  result: SMILESDiagnosticsResponse | null;
  /** True while the API call is in-flight. */
  isLoading: boolean;
  /** Error message if the API call failed. */
  error: string | null;
  /** The original SMILES string that was analyzed. */
  originalSmiles: string;
  /** Called when the user clicks a fix chip; triggers re-validation. */
  onFixApplied: (correctedSmiles: string) => void;
  /** Called when the user clicks the Retry button. */
  onRetry: () => void;
}

/**
 * SMILES Diagnostics panel (DIAG-01).
 *
 * Renders three states:
 * 1. Loading — 3 Skeleton lines.
 * 2. Error — ClayCard error card with retry button.
 * 3. Valid SMILES — SMILESAnnotated with optional warnings + Valid badge.
 * 4. Invalid SMILES — Parse Error heading, SMILESAnnotated with caret, fix chips.
 *
 * Wrapped in a Framer Motion fade-in-up entry animation (0.6s ease-out).
 */
export function SMILESDiagnosticsPanel({
  result,
  isLoading,
  error,
  originalSmiles,
  onFixApplied,
  onRetry,
}: SMILESDiagnosticsPanelProps) {
  // ── Loading state ──
  if (isLoading) {
    return (
      <div className="space-y-3 py-2">
        <Skeleton variant="text" height={16} width="90%" />
        <Skeleton variant="text" height={16} width="75%" />
        <Skeleton variant="text" height={16} width="50%" />
      </div>
    );
  }

  // ── Error state ──
  if (error) {
    return (
      <motion.div
        initial={{ opacity: 0, y: 8 }}
        animate={{ opacity: 1, y: 0 }}
        transition={{ duration: 0.25 }}
      >
        <ClayCard variant="flat" size="sm" className="border border-status-error/20">
          <div className="flex items-start gap-3">
            <AlertTriangle className="w-4 h-4 text-status-error mt-0.5 shrink-0" />
            <div className="flex-1 min-w-0">
              <p className="text-sm font-semibold text-[var(--color-text-primary)] font-display">
                SMILES diagnostics failed
              </p>
              <p className="text-sm text-[var(--color-text-secondary)] mt-0.5">
                {error}
              </p>
            </div>
            <ClayButton variant="ghost" size="sm" onClick={onRetry} className="shrink-0">
              Retry
            </ClayButton>
          </div>
        </ClayCard>
      </motion.div>
    );
  }

  // ── No result yet ──
  if (!result) {
    return null;
  }

  // ── Valid SMILES ──
  if (result.valid) {
    const displaySmiles = result.canonical_smiles ?? originalSmiles;
    // Derive warning positions from atom indices (best-effort; atom_index is zero-indexed)
    const warningPositions = result.warnings
      .filter(w => w.atom_index !== null)
      .map(w => w.atom_index as number);

    return (
      <motion.div
        initial={{ opacity: 0, y: 20 }}
        animate={{ opacity: 1, y: 0 }}
        transition={{ duration: 0.6, ease: 'easeOut' }}
        className="space-y-4"
      >
        {/* Valid status badge row */}
        <div className="flex items-center gap-3">
          <Badge variant="success" size="md">Valid</Badge>
          <span className="text-sm text-[var(--color-text-secondary)]">
            {result.warnings.length === 0
              ? 'Valid — no issues detected'
              : 'Parsed with warnings — see details below'}
          </span>
        </div>

        {/* Annotated SMILES string */}
        <ClayCard variant="flat" size="sm">
          <SMILESAnnotated
            smiles={displaySmiles}
            warningPositions={warningPositions}
          />
        </ClayCard>

        {/* Warnings list */}
        {result.warnings.length > 0 && (
          <div className="space-y-2">
            {result.warnings.map((w, i) => (
              <div
                key={i}
                className="flex items-start gap-2 text-sm text-[var(--color-text-secondary)]"
              >
                <span className="text-status-warning font-semibold shrink-0">⚠</span>
                <span>{w.message}</span>
              </div>
            ))}
          </div>
        )}
      </motion.div>
    );
  }

  // ── Invalid SMILES ──
  const firstError = result.errors[0];
  const errorPosition = firstError?.position ?? null;
  const suggestions = (firstError?.suggestions ?? []).slice(0, MAX_FIX_CHIPS);

  return (
    <motion.div
      initial={{ opacity: 0, y: 20 }}
      animate={{ opacity: 1, y: 0 }}
      transition={{ duration: 0.6, ease: 'easeOut' }}
      className="space-y-4"
    >
      {/* Error heading */}
      <div className="flex items-center gap-3 flex-wrap">
        <span className="text-sm font-semibold text-[var(--color-text-primary)]">
          {errorPosition !== null
            ? `Parse Error at position ${errorPosition}`
            : 'Parse Error'}
        </span>
        {firstError?.error_type && (
          <Badge variant="error" size="sm">
            {firstError.error_type}
          </Badge>
        )}
      </div>

      {/* Annotated SMILES with caret */}
      <ClayCard variant="flat" size="sm">
        <SMILESAnnotated
          smiles={originalSmiles}
          errorPosition={errorPosition}
        />
      </ClayCard>

      {/* Error message */}
      {firstError?.message && (
        <p className="text-sm text-[var(--color-text-secondary)]">
          {firstError.message}
        </p>
      )}

      {/* Fix suggestion chips */}
      {suggestions.length > 0 && (
        <div className="space-y-2">
          <p className="text-xs font-semibold text-[var(--color-text-muted)] uppercase tracking-wide">
            Suggested fixes
          </p>
          <div className="flex flex-wrap gap-2">
            {suggestions.map((suggestion, i) => (
              <FixChip
                key={i}
                suggestion={suggestion}
                onApply={onFixApplied}
              />
            ))}
          </div>
        </div>
      )}
    </motion.div>
  );
}
