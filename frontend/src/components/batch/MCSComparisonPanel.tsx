/**
 * MCSComparisonPanel
 *
 * 480px slide-out drawer showing Maximum Common Substructure (MCS) comparison
 * between two molecules. Displays both molecules side-by-side, Tanimoto
 * similarity gauge, MCS SMARTS as 2D structure, and property delta table.
 */

import React from 'react';
import { motion, AnimatePresence } from 'framer-motion';
import { X, Loader2, AlertTriangle } from 'lucide-react';
import { MoleculeViewer } from '../molecules/MoleculeViewer';
import { MCSPropertyDelta as MCSPropertyDeltaTable } from './MCSPropertyDelta';
import type { BatchResult } from '../../types/batch';
import type { MCSComparisonResult, PropertyStats } from '../../types/analytics';
import { cn } from '../../lib/utils';

interface MCSComparisonPanelProps {
  molecules: BatchResult[];
  mcsResult: MCSComparisonResult | null;
  mcsLoading: boolean;
  mcsError: string | null;
  onClose: () => void;
  onRemoveMolecule: (index: number) => void;
  datasetStats: PropertyStats[] | null;
}

function getSimilarityConfig(score: number): {
  label: string;
  color: string;
  bgColor: string;
  barColor: string;
} {
  if (score >= 0.85) {
    return {
      label: 'Very Similar',
      color: 'text-emerald-600 dark:text-emerald-400',
      bgColor: 'bg-emerald-500/10',
      barColor: 'bg-emerald-500',
    };
  }
  if (score >= 0.7) {
    return {
      label: 'Similar',
      color: 'text-amber-600 dark:text-amber-400',
      bgColor: 'bg-amber-500/10',
      barColor: 'bg-amber-500',
    };
  }
  if (score >= 0.5) {
    return {
      label: 'Moderate',
      color: 'text-orange-600 dark:text-orange-400',
      bgColor: 'bg-orange-500/10',
      barColor: 'bg-orange-500',
    };
  }
  return {
    label: 'Dissimilar',
    color: 'text-red-600 dark:text-red-400',
    bgColor: 'bg-red-500/10',
    barColor: 'bg-red-500',
  };
}

function TanimotoGauge({ tanimoto }: { tanimoto: number }) {
  const config = getSimilarityConfig(tanimoto);
  const pct = Math.round(tanimoto * 100);

  return (
    <div
      className={cn('rounded-xl border border-[var(--color-border)] p-3', config.bgColor)}
      title={`MCS Tanimoto Similarity: ${tanimoto.toFixed(3)}`}
    >
      <div className="flex items-center justify-between mb-1.5">
        <span className="text-xs font-medium text-[var(--color-text-secondary)]">
          MCS Tanimoto Similarity
        </span>
        <span className={cn('text-xs font-semibold', config.color)}>{config.label}</span>
      </div>
      <div className="flex items-center gap-3">
        <div className="flex-1 h-2 rounded-full bg-[var(--color-surface-sunken)] overflow-hidden">
          <motion.div
            initial={{ width: 0 }}
            animate={{ width: `${pct}%` }}
            transition={{ duration: 0.6, ease: 'easeOut' }}
            className={cn('h-full rounded-full', config.barColor)}
          />
        </div>
        <span
          className={cn(
            'text-sm font-bold font-mono tabular-nums min-w-[3ch] text-right',
            config.color
          )}
        >
          {pct}%
        </span>
      </div>
    </div>
  );
}

export const MCSComparisonPanel = React.memo(function MCSComparisonPanel({
  molecules,
  mcsResult,
  mcsLoading,
  mcsError,
  onClose,
  onRemoveMolecule,
}: MCSComparisonPanelProps) {
  return (
    <AnimatePresence>
      {molecules.length === 2 && (
        <>
          {/* Backdrop */}
          <motion.div
            initial={{ opacity: 0 }}
            animate={{ opacity: 1 }}
            exit={{ opacity: 0 }}
            transition={{ duration: 0.2 }}
            className="fixed inset-0 bg-black/30 z-40"
            onClick={onClose}
          />

          {/* Drawer */}
          <motion.div
            initial={{ x: 480 }}
            animate={{ x: 0 }}
            exit={{ x: 480 }}
            transition={{ type: 'spring', damping: 25, stiffness: 300 }}
            role="dialog"
            aria-modal="true"
            aria-labelledby="mcs-panel-title"
            className="fixed right-0 top-0 h-full w-[480px] max-w-[90vw] bg-[var(--color-surface-elevated)] border-l border-[var(--color-border)] shadow-2xl z-50 overflow-y-auto"
          >
            {/* Header */}
            <div className="sticky top-0 z-10 bg-[var(--color-surface-elevated)] border-b border-[var(--color-border)] px-5 py-4 flex items-center justify-between">
              <h2
                id="mcs-panel-title"
                className="text-lg font-semibold text-[var(--color-text-primary)] font-display"
              >
                MCS Comparison
              </h2>
              <button
                onClick={onClose}
                aria-label="Close MCS comparison panel"
                className="p-2 rounded-lg hover:bg-[var(--color-surface-sunken)] transition-colors text-[var(--color-text-muted)]"
              >
                <X className="w-5 h-5" />
              </button>
            </div>

            <div className="p-5 space-y-6">
              {/* Molecule A and B side by side */}
              <div className="grid grid-cols-2 gap-4">
                {molecules.map((mol, i) => (
                  <div
                    key={mol.index}
                    className="rounded-xl border border-[var(--color-border)] bg-[var(--color-surface)] p-3"
                  >
                    <div className="flex items-center justify-between mb-2">
                      <span className="text-xs font-semibold text-[var(--color-text-primary)]">
                        Molecule {String.fromCharCode(65 + i)}
                      </span>
                      <button
                        onClick={() => onRemoveMolecule(i)}
                        className="p-1 rounded hover:bg-red-500/10 transition-colors text-[var(--color-text-muted)] hover:text-red-500"
                        title="Remove from comparison"
                      >
                        <X className="w-3.5 h-3.5" />
                      </button>
                    </div>
                    <div className="flex items-center justify-center rounded-lg bg-white dark:bg-gray-900/50 border border-[var(--color-border-subtle)] min-h-[140px]">
                      <MoleculeViewer
                        smiles={mol.standardization?.standardized_smiles || mol.smiles}
                        width={180}
                        height={140}
                      />
                    </div>
                    <p
                      className="mt-2 text-xs font-mono text-[var(--color-text-muted)] truncate"
                      title={mol.smiles}
                    >
                      {mol.smiles}
                    </p>
                    {mol.name && (
                      <p className="text-xs text-[var(--color-text-secondary)] truncate">
                        {mol.name}
                      </p>
                    )}
                  </div>
                ))}
              </div>

              {/* Tanimoto Similarity Gauge */}
              {mcsResult && <TanimotoGauge tanimoto={mcsResult.tanimoto} />}

              {/* MCS loading state */}
              {mcsLoading && (
                <div className="flex items-center justify-center gap-2 py-6 text-sm text-[var(--color-text-muted)]">
                  <Loader2 className="w-4 h-4 animate-spin text-[var(--color-primary)]" />
                  Computing maximum common substructure...
                </div>
              )}

              {/* MCS error state */}
              {mcsError && !mcsLoading && (
                <div className="rounded-xl p-3 bg-amber-500/5 border border-amber-500/20">
                  <div className="flex items-start gap-2">
                    <AlertTriangle className="w-4 h-4 text-amber-500 flex-shrink-0 mt-0.5" />
                    <p className="text-xs text-amber-600 dark:text-amber-400">
                      MCS computation timed out or failed. The property comparison is still
                      available below.
                    </p>
                  </div>
                </div>
              )}

              {/* MCS Result section */}
              {mcsResult && !mcsLoading && (
                <div className="space-y-3">
                  <h3 className="text-base font-semibold text-[var(--color-text-primary)] font-display">
                    Maximum Common Substructure
                  </h3>

                  {/* MCS SMARTS as 2D structure */}
                  <div
                    className="flex items-center justify-center rounded-xl border border-[var(--color-border)] bg-white dark:bg-gray-900/50 p-3 min-h-[120px]"
                    aria-label={`Maximum common substructure: ${mcsResult.mcs_smarts}`}
                  >
                    <MoleculeViewer smiles={mcsResult.mcs_smarts} width={200} height={120} />
                  </div>

                  {/* Matched count */}
                  <p className="text-sm text-[var(--color-text-secondary)]">
                    Matched: {mcsResult.num_atoms} atoms, {mcsResult.num_bonds} bonds
                  </p>

                  {/* Timeout warning */}
                  {mcsResult.timed_out && (
                    <div className="rounded-xl p-3 bg-amber-500/5 border border-amber-500/20">
                      <p className="text-xs text-amber-600 dark:text-amber-400">
                        MCS search is limited to 10 seconds. Complex molecule pairs may not produce
                        a result.
                      </p>
                    </div>
                  )}
                </div>
              )}

              {/* Property delta table */}
              {mcsResult?.property_deltas && mcsResult.property_deltas.length > 0 && (
                <div className="space-y-2">
                  <h3 className="text-sm font-semibold text-[var(--color-text-primary)] font-display">
                    Property Comparison
                  </h3>
                  <MCSPropertyDeltaTable deltas={mcsResult.property_deltas} />
                </div>
              )}
            </div>
          </motion.div>
        </>
      )}
    </AnimatePresence>
  );
});
