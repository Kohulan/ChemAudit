/**
 * MoleculeComparisonPanel (VIZ-07)
 *
 * Side-by-side molecule comparison drawer with 2D structures,
 * properties comparison table, and overlaid property radar.
 * Slides in from the right.
 */

import React, { useMemo } from 'react';
import { motion, AnimatePresence } from 'framer-motion';
import { X } from 'lucide-react';
import { MoleculeViewer } from '../molecules/MoleculeViewer';
import { MoleculePropertyRadar } from './MoleculePropertyRadar';
import type { BatchResult } from '../../types/batch';
import type { PropertyStats } from '../../types/analytics';
import { cn } from '../../lib/utils';

interface MoleculeComparisonPanelProps {
  molecules: BatchResult[];
  datasetStats: PropertyStats[] | null;
  onClose: () => void;
  onRemoveMolecule: (index: number) => void;
}

interface ComparisonRow {
  label: string;
  key: string;
  higherIsBetter: boolean;
  format: (val: number | null) => string;
}

const COMPARISON_ROWS: ComparisonRow[] = [
  {
    label: 'Overall Score',
    key: 'overall_score',
    higherIsBetter: true,
    format: (v) => (v !== null ? v.toFixed(0) : '-'),
  },
  {
    label: 'QED',
    key: 'qed',
    higherIsBetter: true,
    format: (v) => (v !== null ? v.toFixed(3) : '-'),
  },
  {
    label: 'SA Score',
    key: 'sa_score',
    higherIsBetter: false,
    format: (v) => (v !== null ? v.toFixed(1) : '-'),
  },
  {
    label: 'Fsp3',
    key: 'fsp3',
    higherIsBetter: true,
    format: (v) => (v !== null ? v.toFixed(3) : '-'),
  },
  {
    label: 'Lipinski Violations',
    key: 'lipinski_violations',
    higherIsBetter: false,
    format: (v) => (v !== null ? v.toFixed(0) : '-'),
  },
  {
    label: 'Alert Count',
    key: 'alert_count',
    higherIsBetter: false,
    format: (v) => (v !== null ? v.toFixed(0) : '-'),
  },
];

function getValue(r: BatchResult, key: string): number | null {
  switch (key) {
    case 'overall_score':
      return r.validation?.overall_score ?? null;
    case 'qed':
      return r.scoring?.druglikeness?.qed_score ?? null;
    case 'sa_score':
      return r.scoring?.admet?.sa_score ?? null;
    case 'fsp3':
      return r.scoring?.admet?.fsp3 ?? null;
    case 'lipinski_violations':
      return r.scoring?.druglikeness?.lipinski_violations ?? null;
    case 'alert_count':
      return r.alerts?.alert_count ?? null;
    default:
      return null;
  }
}

function getCellColor(
  val: number | null,
  otherVal: number | null,
  higherIsBetter: boolean
): string {
  if (val === null || otherVal === null) return '';
  if (val === otherVal) return '';
  const isBetter = higherIsBetter ? val > otherVal : val < otherVal;
  return isBetter
    ? 'text-emerald-600 dark:text-emerald-400 bg-emerald-500/5'
    : 'text-red-600 dark:text-red-400 bg-red-500/5';
}

export const MoleculeComparisonPanel = React.memo(function MoleculeComparisonPanel({
  molecules,
  datasetStats,
  onClose,
  onRemoveMolecule,
}: MoleculeComparisonPanelProps) {
  const hasTwoMolecules = molecules.length === 2;

  const values = useMemo(() => {
    return molecules.map((mol) => {
      const row: Record<string, number | null> = {};
      for (const r of COMPARISON_ROWS) {
        row[r.key] = getValue(mol, r.key);
      }
      return row;
    });
  }, [molecules]);

  return (
    <AnimatePresence>
      {molecules.length > 0 && (
        <>
          {/* Backdrop */}
          <motion.div
            initial={{ opacity: 0 }}
            animate={{ opacity: 1 }}
            exit={{ opacity: 0 }}
            className="fixed inset-0 bg-black/30 z-40"
            onClick={onClose}
          />

          {/* Drawer */}
          <motion.div
            initial={{ x: 480 }}
            animate={{ x: 0 }}
            exit={{ x: 480 }}
            transition={{ type: 'spring', damping: 25, stiffness: 300 }}
            className="fixed right-0 top-0 h-full w-[480px] max-w-[90vw] bg-[var(--color-surface-elevated)] border-l border-[var(--color-border)] shadow-2xl z-50 overflow-y-auto"
          >
            {/* Header */}
            <div className="sticky top-0 z-10 bg-[var(--color-surface-elevated)] border-b border-[var(--color-border)] px-5 py-4 flex items-center justify-between">
              <h2 className="text-lg font-semibold text-[var(--color-text-primary)] font-display">
                Compare Molecules
              </h2>
              <button
                onClick={onClose}
                className="p-2 rounded-lg hover:bg-[var(--color-surface-sunken)] transition-colors text-[var(--color-text-muted)]"
              >
                <X className="w-5 h-5" />
              </button>
            </div>

            <div className="p-5 space-y-6">
              {/* Structure cards */}
              <div className={cn(
                'grid gap-4',
                hasTwoMolecules ? 'grid-cols-2' : 'grid-cols-1'
              )}>
                {molecules.map((mol, i) => (
                  <div
                    key={mol.index}
                    className="rounded-xl border border-[var(--color-border)] bg-[var(--color-surface)] p-3"
                  >
                    <div className="flex items-center justify-between mb-2">
                      <span className="text-xs font-semibold text-[var(--color-text-primary)]">
                        {hasTwoMolecules ? `Molecule ${String.fromCharCode(65 + i)}` : 'Molecule'}
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
                      className="mt-2 text-[10px] font-mono text-[var(--color-text-muted)] truncate"
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

              {/* Properties comparison table */}
              <div className="rounded-xl border border-[var(--color-border)] overflow-hidden">
                <table className="w-full text-xs">
                  <thead>
                    <tr className="bg-[var(--color-surface-sunken)]">
                      <th className="px-3 py-2 text-left text-[var(--color-text-muted)] font-medium">
                        Property
                      </th>
                      {molecules.map((_, i) => (
                        <th
                          key={i}
                          className="px-3 py-2 text-center text-[var(--color-text-muted)] font-medium"
                        >
                          {hasTwoMolecules ? `Mol ${String.fromCharCode(65 + i)}` : 'Value'}
                        </th>
                      ))}
                    </tr>
                  </thead>
                  <tbody className="divide-y divide-[var(--color-border)]">
                    {COMPARISON_ROWS.map((row) => (
                      <tr key={row.key}>
                        <td className="px-3 py-2 text-[var(--color-text-secondary)] font-medium">
                          {row.label}
                        </td>
                        {molecules.map((_, i) => {
                          const val = values[i]?.[row.key] ?? null;
                          const otherVal = hasTwoMolecules
                            ? (values[1 - i]?.[row.key] ?? null)
                            : null;
                          const cellColor = hasTwoMolecules
                            ? getCellColor(val, otherVal, row.higherIsBetter)
                            : '';

                          return (
                            <td
                              key={i}
                              className={cn(
                                'px-3 py-2 text-center font-mono',
                                cellColor || 'text-[var(--color-text-primary)]'
                              )}
                            >
                              {row.format(val)}
                            </td>
                          );
                        })}
                      </tr>
                    ))}
                  </tbody>
                </table>
              </div>

              {/* Radar Chart */}
              <div>
                <h3 className="text-sm font-semibold text-[var(--color-text-primary)] mb-2 font-display">
                  Property Radar
                </h3>
                <MoleculePropertyRadar
                  molecules={molecules}
                  datasetStats={datasetStats}
                />
              </div>
            </div>
          </motion.div>
        </>
      )}
    </AnimatePresence>
  );
});
