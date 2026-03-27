import { useEffect, useRef } from 'react';
import { motion, AnimatePresence } from 'framer-motion';
import { AlertTriangle } from 'lucide-react';
import { ClayCard } from '../ui/ClayCard';
import { ClayButton } from '../ui/ClayButton';
import { Badge } from '../ui/Badge';
import { Skeleton } from '../ui/Skeleton';
import { Tooltip } from '../ui/Tooltip';
import { MoleculeViewer } from '../molecules/MoleculeViewer';
import type { CrossPipelineResponse, PropertyComparison } from '../../types/diagnostics';

interface CrossPipelinePanelProps {
  result: CrossPipelineResponse | null;
  isLoading: boolean;
  error: string | null;
  currentSmiles: string | null;
  onComparePipelines: (molecule: string) => void;
  onRetry: () => void;
}

const PIPELINE_LABELS = ['RDKit MolStandardize', 'ChEMBL-style', 'Minimal Sanitize'];

/** Get row highlight classes based on disagreement type. */
function getRowHighlightClass(row: PropertyComparison): string {
  if (row.agrees) return '';
  if (row.structural) {
    return 'bg-[rgba(239,68,68,0.15)] border-l-2 border-status-error';
  }
  return 'bg-[rgba(245,158,11,0.15)] border-l-2 border-status-warning';
}

/** Format a property value for display. */
function formatValue(value: string | number, property: string): string {
  if (typeof value === 'number') {
    if (property.toLowerCase().includes('mw') || property === 'MW') {
      return value.toFixed(2);
    }
    return String(value);
  }
  return String(value);
}

/** Truncate long strings for table display. */
function truncate(str: string, maxLen = 28): string {
  if (str.length <= maxLen) return str;
  return str.slice(0, maxLen) + '...';
}

/**
 * Cross-pipeline standardization comparison panel (DIAG-04).
 *
 * Per UI-SPEC (D-09):
 * - 3 MoleculeViewer thumbnails in a 3-column grid
 * - Property comparison table with disagreement cell highlighting
 * - Verdict row: all-agree (amber), warning (property diff), error (structural diff)
 * - Auto-triggers on mount if currentSmiles is set and result is null
 */
export function CrossPipelinePanel({
  result,
  isLoading,
  error,
  currentSmiles,
  onComparePipelines,
  onRetry,
}: CrossPipelinePanelProps) {
  const hasTriggered = useRef(false);

  // Auto-trigger when currentSmiles is available and no result yet
  useEffect(() => {
    if (currentSmiles && !result && !isLoading && !error && !hasTriggered.current) {
      hasTriggered.current = true;
      onComparePipelines(currentSmiles);
    }
  }, [currentSmiles, result, isLoading, error, onComparePipelines]);

  // Network error state
  if (error) {
    return (
      <motion.div
        initial={{ opacity: 0, y: 8 }}
        animate={{ opacity: 1, y: 0 }}
        transition={{ duration: 0.25 }}
      >
        <ClayCard variant="flat" size="sm" className="border border-[var(--color-border)]">
          <div className="flex items-start gap-3">
            <AlertTriangle className="w-4 h-4 text-status-error mt-0.5 shrink-0" />
            <div className="flex-1 min-w-0">
              <p className="text-sm font-semibold text-[var(--color-text-primary)] font-display">
                Pipeline comparison failed
              </p>
              <p className="text-sm text-[var(--color-text-secondary)] mt-0.5">
                Pipeline comparison failed. Verify the SMILES is valid and retry.
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

  // Loading state: 3 skeleton thumbnails + 7 skeleton rows
  if (isLoading) {
    return (
      <div className="space-y-4">
        {/* 3 skeleton MoleculeViewer placeholders */}
        <div className="grid grid-cols-1 sm:grid-cols-3 gap-2">
          {[0, 1, 2].map((i) => (
            <div key={i}>
              <Skeleton variant="rounded" height={160} className="w-full mb-2" />
              <Skeleton variant="text" height={14} className="w-3/4 mx-auto" />
            </div>
          ))}
        </div>
        {/* 7 skeleton table rows */}
        <ClayCard variant="flat" size="sm" className="p-0 overflow-hidden">
          <div className="space-y-0">
            {Array.from({ length: 7 }).map((_, i) => (
              <div key={i} className="flex gap-3 px-4 py-3 border-b border-[var(--color-border)] last:border-0">
                <Skeleton variant="text" height={14} className="w-1/5" />
                <Skeleton variant="text" height={14} className="flex-1" />
                <Skeleton variant="text" height={14} className="flex-1" />
                <Skeleton variant="text" height={14} className="flex-1" />
              </div>
            ))}
          </div>
        </ClayCard>
      </div>
    );
  }

  if (!result) return null;

  // Extract SMILES for the 3 pipelines (in order)
  const pipelineSmiles = result.pipelines.map((p) => p.smiles);

  return (
    <AnimatePresence>
      <motion.div
        initial={{ opacity: 0, y: 12 }}
        animate={{ opacity: 1, y: 0 }}
        transition={{ duration: 0.4, ease: 'easeOut' }}
        className="space-y-4"
      >
        {/* 3 MoleculeViewer thumbnails */}
        <div className="grid grid-cols-1 sm:grid-cols-3 gap-2">
          {PIPELINE_LABELS.map((label, i) => (
            <div key={label} className="flex flex-col items-center gap-2">
              <div
                className="w-full"
                style={{ height: '160px' }}
                aria-label={`${label} — canonical structure rendering`}
              >
                <MoleculeViewer
                  smiles={pipelineSmiles[i] ?? null}
                  height={160}
                  className="h-[160px]"
                />
              </div>
              <span className="text-xs font-semibold text-[var(--color-text-muted)] text-center">
                {label}
              </span>
            </div>
          ))}
        </div>

        {/* Comparison table */}
        <ClayCard variant="flat" size="sm" className="p-0 overflow-hidden">
          <div className="overflow-x-auto">
            <table className="w-full text-sm">
              <thead>
                <tr className="border-b border-[var(--color-border)] bg-[var(--color-surface-sunken)]">
                  <th className="text-left px-4 py-2 text-xs font-semibold text-[var(--color-text-muted)] w-[120px]">
                    Property
                  </th>
                  {PIPELINE_LABELS.map((label) => (
                    <th
                      key={label}
                      className="text-left px-4 py-2 text-xs font-semibold text-[var(--color-text-muted)]"
                    >
                      {label.split(' ')[0]}
                    </th>
                  ))}
                </tr>
              </thead>
              <tbody>
                {result.property_comparison.map((row, rowIdx) => {
                  const highlightClass = getRowHighlightClass(row);
                  const isMonoProperty =
                    row.property === 'SMILES' || row.property === 'InChIKey';
                  return (
                    <tr
                      key={rowIdx}
                      className={`border-b border-[var(--color-border)] last:border-0 ${highlightClass}`}
                      aria-describedby={!row.agrees ? `disagree-summary-${rowIdx}` : undefined}
                    >
                      <td className="px-4 py-2.5 text-xs font-semibold text-[var(--color-text-secondary)]">
                        {row.property}
                        {!row.agrees && (
                          <span id={`disagree-summary-${rowIdx}`} className="sr-only">
                            {row.structural ? 'Structural disagreement' : 'Property disagreement'} in {row.property}
                          </span>
                        )}
                      </td>
                      {row.values.map((val, colIdx) => (
                        <td
                          key={colIdx}
                          className={`px-4 py-2.5 ${isMonoProperty ? 'font-mono text-xs' : 'text-sm'} text-[var(--color-text-primary)]`}
                        >
                          {isMonoProperty ? (
                            <Tooltip content={String(val)}>
                              <span>{truncate(formatValue(val, row.property))}</span>
                            </Tooltip>
                          ) : (
                            formatValue(val, row.property)
                          )}
                        </td>
                      ))}
                    </tr>
                  );
                })}
              </tbody>
            </table>
          </div>

          {/* Verdict row */}
          <div className="px-4 py-3 border-t border-[var(--color-border)] flex items-center gap-2">
            {result.all_agree ? (
              <Badge variant="success">All pipelines agree</Badge>
            ) : result.structural_disagreements > 0 ? (
              <Badge variant="error">
                {result.structural_disagreements} structural disagreement{result.structural_disagreements > 1 ? 's' : ''} — SMILES or InChIKey differ
              </Badge>
            ) : (
              <Badge variant="warning">
                {result.disagreements} disagreement{result.disagreements > 1 ? 's' : ''} found
              </Badge>
            )}
          </div>
        </ClayCard>
      </motion.div>
    </AnimatePresence>
  );
}
