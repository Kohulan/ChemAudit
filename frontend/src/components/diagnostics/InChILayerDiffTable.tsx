import { useState } from 'react';
import { motion } from 'framer-motion';
import { AlertTriangle } from 'lucide-react';
import { ClayCard } from '../ui/ClayCard';
import { ClayButton } from '../ui/ClayButton';
import { Badge } from '../ui/Badge';
import { Skeleton } from '../ui/Skeleton';
import { InChILayerRow } from './InChILayerRow';
import type { InChIDiffResponse } from '../../types/diagnostics';

interface InChILayerDiffTableProps {
  /** Result from the InChI diff API (DIAG-02). */
  result: InChIDiffResponse | null;
  /** True while the comparison API call is in-flight. */
  isLoading: boolean;
  /** Error message if the comparison API call failed. */
  error: string | null;
  /** Called when the user submits two InChI strings for comparison. */
  onCompare: (inchiA: string, inchiB: string) => void;
  /** Called when the user clicks the Retry button. */
  onRetry: () => void;
  /** Auto-populated InChI-A from the current molecule analysis. */
  initialInchiA?: string;
}

/**
 * InChI layer-by-layer diff table (DIAG-02).
 *
 * Renders two text inputs (InChI-A and InChI-B), a Compare button, and
 * a table of per-layer match/mismatch indicators using InChILayerRow.
 *
 * Loading: 6 skeleton rows.
 * Error: error card with retry.
 * Result: sticky header + InChILayerRow per layer + summary badge.
 */
export function InChILayerDiffTable({
  result,
  isLoading,
  error,
  onCompare,
  onRetry,
  initialInchiA = '',
}: InChILayerDiffTableProps) {
  const [inchiA, setInchiA] = useState(initialInchiA);
  const [inchiB, setInchiB] = useState('');

  const handleCompare = () => {
    onCompare(inchiA.trim(), inchiB.trim());
  };

  const canCompare = inchiA.trim().length > 0 && inchiB.trim().length > 0;

  // Count mismatching layers
  const mismatchCount = result?.layer_rows.filter(r => !r.match).length ?? 0;

  return (
    <div className="space-y-4">
      {/* ── Input pair ── */}
      <div className="grid grid-cols-1 sm:grid-cols-2 gap-3">
        <div className="space-y-1">
          <label
            htmlFor="inchi-a-input"
            className="text-xs font-semibold text-[var(--color-text-secondary)] uppercase tracking-wide"
          >
            InChI A
          </label>
          <input
            id="inchi-a-input"
            type="text"
            value={inchiA}
            onChange={e => setInchiA(e.target.value)}
            placeholder="InChI=1S/..."
            className={[
              'w-full px-3 py-2 rounded-xl',
              'text-sm font-mono text-[var(--color-text-primary)]',
              'bg-[var(--color-surface-sunken)]',
              'border border-[var(--color-border)]',
              'focus:outline-none focus:ring-2 focus:ring-[var(--color-primary)]/30',
              'placeholder:text-[var(--color-text-muted)]',
            ].join(' ')}
            autoComplete="off"
            spellCheck={false}
          />
        </div>
        <div className="space-y-1">
          <label
            htmlFor="inchi-b-input"
            className="text-xs font-semibold text-[var(--color-text-secondary)] uppercase tracking-wide"
          >
            InChI B
          </label>
          <input
            id="inchi-b-input"
            type="text"
            value={inchiB}
            onChange={e => setInchiB(e.target.value)}
            placeholder="InChI=1S/..."
            className={[
              'w-full px-3 py-2 rounded-xl',
              'text-sm font-mono text-[var(--color-text-primary)]',
              'bg-[var(--color-surface-sunken)]',
              'border border-[var(--color-border)]',
              'focus:outline-none focus:ring-2 focus:ring-[var(--color-primary)]/30',
              'placeholder:text-[var(--color-text-muted)]',
            ].join(' ')}
            autoComplete="off"
            spellCheck={false}
          />
        </div>
      </div>

      {/* Compare button */}
      <div>
        <ClayButton
          variant="primary"
          size="sm"
          onClick={handleCompare}
          disabled={!canCompare || isLoading}
          loading={isLoading}
        >
          Compare Layers
        </ClayButton>
      </div>

      {/* ── Loading state — 6 skeleton rows ── */}
      {isLoading && (
        <div className="space-y-2 pt-2">
          {Array.from({ length: 6 }).map((_, i) => (
            <Skeleton key={i} variant="text" height={32} width="100%" />
          ))}
        </div>
      )}

      {/* ── Error state ── */}
      {!isLoading && error && (
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
                  InChI comparison failed
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
      )}

      {/* ── Result table ── */}
      {!isLoading && result && (
        <motion.div
          initial={{ opacity: 0, y: 12 }}
          animate={{ opacity: 1, y: 0 }}
          transition={{ duration: 0.4, ease: 'easeOut' }}
        >
          <ClayCard variant="default" size="sm" className="overflow-hidden p-0">
            {/* Verdict badge row */}
            <div className="px-4 py-3 border-b border-[var(--color-border)] flex items-center justify-between">
              <span className="text-sm font-semibold font-display text-[var(--color-text-primary)]">
                Layer comparison
              </span>
              {result.identical ? (
                <Badge variant="success" size="md">All layers match</Badge>
              ) : (
                <Badge variant="error" size="md">
                  {mismatchCount} layer{mismatchCount !== 1 ? 's' : ''} differ
                </Badge>
              )}
            </div>

            {/* Scrollable table body */}
            <div className={result.layer_rows.length > 6 ? 'overflow-y-auto max-h-64' : ''}>
              <table className="w-full">
                {/* Sticky header */}
                <thead className="sticky top-0 bg-[var(--color-surface-elevated)] z-10">
                  <tr className="border-b border-[var(--color-border)]">
                    <th className="px-4 py-2 text-left text-xs font-semibold text-[var(--color-text-muted)] uppercase tracking-wide">
                      Layer
                    </th>
                    <th className="px-4 py-2 text-left text-xs font-semibold text-[var(--color-text-muted)] uppercase tracking-wide">
                      InChI A
                    </th>
                    <th className="px-4 py-2 text-left text-xs font-semibold text-[var(--color-text-muted)] uppercase tracking-wide">
                      InChI B
                    </th>
                    <th className="px-4 py-2 text-center text-xs font-semibold text-[var(--color-text-muted)] uppercase tracking-wide">
                      Match
                    </th>
                  </tr>
                </thead>
                {/* InChILayerRow renders <tr> elements — use directly in <tbody> */}
                <tbody>
                  {result.layer_rows.map((row, i) => (
                    <InChILayerRow
                      key={i}
                      layer={row.layer}
                      valueA={row.value_a}
                      valueB={row.value_b}
                      match={row.match}
                    />
                  ))}
                </tbody>
              </table>
            </div>
          </ClayCard>
        </motion.div>
      )}
    </div>
  );
}
