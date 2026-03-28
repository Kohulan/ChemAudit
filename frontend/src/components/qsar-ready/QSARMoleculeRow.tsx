import { useState } from 'react';
import { motion, AnimatePresence } from 'framer-motion';
import { ChevronDown } from 'lucide-react';
import { cn } from '../../lib/utils';
import type { QSARReadyResult, QSARStepResult } from '../../types/qsar_ready';

interface QSARMoleculeRowProps {
  result: QSARReadyResult;
  index: number;
}

// ─── Status helpers ──────────────────────────────────────────────────────────

function getStatusColor(status: QSARReadyResult['status']): string {
  switch (status) {
    case 'ok':
      return '#16a34a';
    case 'rejected':
      return '#dc2626';
    case 'duplicate':
      return '#d97706';
    case 'error':
      return '#7c3aed';
    default:
      return 'var(--color-text-muted)';
  }
}

function getStatusBgClass(status: QSARReadyResult['status']): string {
  switch (status) {
    case 'ok':
      return 'bg-green-50 text-green-700';
    case 'rejected':
      return 'bg-red-50 text-red-700';
    case 'duplicate':
      return 'bg-amber-50 text-amber-700';
    case 'error':
      return 'bg-violet-50 text-violet-700';
    default:
      return 'bg-[var(--color-surface-sunken)] text-[var(--color-text-muted)]';
  }
}

function getStepChipClass(status: QSARStepResult['status']): string {
  switch (status) {
    case 'applied':
      return 'bg-green-50 text-green-700';
    case 'no_change':
      return 'bg-blue-50 text-blue-700';
    case 'skipped':
      return 'bg-[var(--color-surface-sunken)] text-[var(--color-text-muted)]';
    case 'error':
      return 'bg-red-50 text-red-600';
    default:
      return 'bg-[var(--color-surface-sunken)] text-[var(--color-text-muted)]';
  }
}

function getStepDotColor(status: QSARStepResult['status']): string {
  switch (status) {
    case 'applied':
      return '#16a34a';
    case 'no_change':
      return '#2563eb';
    case 'skipped':
      return 'var(--color-text-muted)';
    case 'error':
      return '#dc2626';
    default:
      return 'var(--color-text-muted)';
  }
}

// ─── Component ───────────────────────────────────────────────────────────────

/**
 * A single expandable row in the QSAR batch results table.
 *
 * Per UI-SPEC interaction contract:
 * - Click anywhere on row to expand/collapse
 * - Chevron rotates 180° on expand
 * - aria-expanded on row trigger for accessibility
 * - Expanded content: horizontal mini step chips (flex-wrap, not vertical timeline)
 * - Each chip: step name + status color dot
 * - Expand animation: height 0→auto, 250ms ease-out
 */
export function QSARMoleculeRow({ result, index }: QSARMoleculeRowProps) {
  const [expanded, setExpanded] = useState(false);

  const statusColor = getStatusColor(result.status);
  const statusBgClass = getStatusBgClass(result.status);

  // Truncate long SMILES for display
  const truncate = (s: string | null, maxLen = 24): string => {
    if (!s) return '—';
    return s.length > maxLen ? `${s.slice(0, maxLen)}…` : s;
  };

  return (
    <>
      {/* ── Main row ── */}
      <tr
        className="border-b border-[var(--color-border)] hover:bg-[var(--color-surface-sunken)]/40 cursor-pointer transition-colors"
        onClick={() => setExpanded((v) => !v)}
        aria-expanded={expanded}
      >
        {/* Index */}
        <td className="px-4 py-3 text-xs text-[var(--color-text-muted)] w-12">
          {index + 1}
        </td>

        {/* Original SMILES */}
        <td className="px-4 py-3 max-w-[160px]">
          <span
            className="font-mono text-xs text-[var(--color-text-primary)] truncate block"
            title={result.original_smiles}
          >
            {truncate(result.original_smiles)}
          </span>
        </td>

        {/* Curated SMILES */}
        <td className="px-4 py-3 max-w-[160px]">
          <span
            className="font-mono text-xs text-[var(--color-text-primary)] truncate block"
            title={result.curated_smiles ?? undefined}
          >
            {truncate(result.curated_smiles)}
          </span>
        </td>

        {/* InChIKey change indicator */}
        <td className="px-4 py-3">
          {result.inchikey_changed ? (
            <span className="inline-flex items-center gap-1 px-2 py-0.5 rounded-full text-xs font-medium bg-amber-100 text-amber-700">
              <span
                className="w-1.5 h-1.5 rounded-full flex-shrink-0"
                style={{ background: '#d97706' }}
              />
              InChIKey changed
            </span>
          ) : (
            <span className="inline-flex items-center gap-1 px-2 py-0.5 rounded-full text-xs font-medium bg-[var(--color-surface-sunken)] text-[var(--color-text-muted)]">
              <span
                className="w-1.5 h-1.5 rounded-full flex-shrink-0"
                style={{ background: 'var(--color-text-muted)' }}
              />
              InChIKey preserved
            </span>
          )}
        </td>

        {/* Status badge */}
        <td className="px-4 py-3">
          <span
            className={cn(
              'inline-flex items-center px-2 py-0.5 rounded-full text-xs font-semibold capitalize',
              statusBgClass,
            )}
            style={{ color: statusColor }}
          >
            {result.status}
          </span>
        </td>

        {/* Expand chevron */}
        <td className="px-3 py-3 w-8">
          <ChevronDown
            className={cn(
              'w-4 h-4 text-[var(--color-text-muted)] transition-transform duration-200',
              expanded && 'rotate-180',
            )}
          />
        </td>
      </tr>

      {/* ── Expanded detail row ── */}
      <AnimatePresence>
        {expanded && (
          <tr>
            <td colSpan={6} className="px-4 pb-4 pt-0 bg-[var(--color-surface-sunken)]/30">
              <motion.div
                initial={{ opacity: 0, height: 0 }}
                animate={{ opacity: 1, height: 'auto' }}
                exit={{ opacity: 0, height: 0 }}
                transition={{ duration: 0.25, ease: 'easeOut' }}
                className="overflow-hidden"
              >
                {/* Rejection reason (if any) */}
                {result.rejection_reason && (
                  <p className="text-xs text-[var(--color-text-secondary)] mb-2 mt-2">
                    <span className="font-medium">Reason:</span> {result.rejection_reason}
                  </p>
                )}

                {/* InChIKey values */}
                {(result.original_inchikey || result.standardized_inchikey) && (
                  <div className="flex flex-wrap gap-4 text-xs text-[var(--color-text-muted)] mb-2 mt-2">
                    {result.original_inchikey && (
                      <span>
                        <span className="font-medium">Original InChIKey:</span>{' '}
                        <span
                          className="font-mono"
                          aria-label={`InChIKey: ${result.original_inchikey}`}
                        >
                          {result.original_inchikey}
                        </span>
                      </span>
                    )}
                    {result.standardized_inchikey && result.inchikey_changed && (
                      <span>
                        <span className="font-medium">Curated InChIKey:</span>{' '}
                        <span
                          className="font-mono"
                          aria-label={`InChIKey: ${result.standardized_inchikey}`}
                        >
                          {result.standardized_inchikey}
                        </span>
                      </span>
                    )}
                  </div>
                )}

                {/* Mini step chips */}
                {result.steps.length > 0 && (
                  <div className="flex flex-wrap gap-2 mt-3">
                    {result.steps.map((step) => (
                      <span
                        key={step.step_index}
                        className={cn(
                          'inline-flex items-center gap-1.5 px-2 py-1 rounded-full text-xs',
                          getStepChipClass(step.status),
                        )}
                        title={step.detail ?? step.step_name}
                      >
                        <span
                          className="w-1.5 h-1.5 rounded-full flex-shrink-0"
                          style={{ background: getStepDotColor(step.status) }}
                        />
                        {step.step_name}
                      </span>
                    ))}
                  </div>
                )}
              </motion.div>
            </td>
          </tr>
        )}
      </AnimatePresence>
    </>
  );
}
