/**
 * MoleculeComparisonPanel (VIZ-07)
 *
 * Side-by-side molecule comparison drawer with 2D structures,
 * properties comparison table, and overlaid property radar.
 * Slides in from the right.
 */

import React, { useEffect, useMemo, useState } from 'react';
import { motion, AnimatePresence } from 'framer-motion';
import { X } from 'lucide-react';
import { MoleculeViewer } from '../molecules/MoleculeViewer';
import { MoleculePropertyRadar } from './MoleculePropertyRadar';
import { validationApi } from '../../services/api';
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

function getSimilarityConfig(score: number): { label: string; color: string; bgColor: string; barColor: string } {
  if (score >= 0.85) return { label: 'Very Similar', color: 'text-emerald-600 dark:text-emerald-400', bgColor: 'bg-emerald-500/10', barColor: 'bg-emerald-500' };
  if (score >= 0.70) return { label: 'Similar', color: 'text-amber-600 dark:text-amber-400', bgColor: 'bg-amber-500/10', barColor: 'bg-amber-500' };
  if (score >= 0.50) return { label: 'Moderate', color: 'text-orange-600 dark:text-orange-400', bgColor: 'bg-orange-500/10', barColor: 'bg-orange-500' };
  return { label: 'Dissimilar', color: 'text-red-600 dark:text-red-400', bgColor: 'bg-red-500/10', barColor: 'bg-red-500' };
}

function SimilarityGauge({ similarity, loading }: {
  similarity: { score: number; common_bits: number; bits_a: number; bits_b: number } | null;
  loading: boolean;
}) {
  if (loading) {
    return (
      <div className="rounded-xl border border-[var(--color-border)] bg-[var(--color-surface)] p-3 flex items-center justify-center gap-2 text-xs text-[var(--color-text-muted)]">
        <div className="w-3.5 h-3.5 border-2 border-current border-t-transparent rounded-full animate-spin" />
        Computing similarity...
      </div>
    );
  }

  if (!similarity) return null;

  const { score, common_bits, bits_a, bits_b } = similarity;
  const config = getSimilarityConfig(score);
  const pct = Math.round(score * 100);

  return (
    <div className={cn('rounded-xl border border-[var(--color-border)] p-3', config.bgColor)} title={`ECFP4 Tanimoto | Common bits: ${common_bits} | Bits A: ${bits_a} | Bits B: ${bits_b}`}>
      <div className="flex items-center justify-between mb-1.5">
        <span className="text-xs font-medium text-[var(--color-text-secondary)]">ECFP4 Tanimoto Similarity</span>
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
        <span className={cn('text-sm font-bold font-mono tabular-nums min-w-[3ch] text-right', config.color)}>
          {pct}%
        </span>
      </div>
      <div className="flex items-center justify-between mt-1.5 text-[10px] text-[var(--color-text-muted)]">
        <span>Common: {common_bits} bits</span>
        <span>A: {bits_a} | B: {bits_b}</span>
      </div>
    </div>
  );
}

export const MoleculeComparisonPanel = React.memo(function MoleculeComparisonPanel({
  molecules,
  datasetStats,
  onClose,
  onRemoveMolecule,
}: MoleculeComparisonPanelProps) {
  const hasTwoMolecules = molecules.length === 2;

  // ECFP4 Tanimoto similarity
  const [similarity, setSimilarity] = useState<{
    score: number;
    common_bits: number;
    bits_a: number;
    bits_b: number;
  } | null>(null);
  const [simLoading, setSimLoading] = useState(false);

  useEffect(() => {
    if (!hasTwoMolecules) {
      setSimilarity(null);
      return;
    }
    const smilesA = molecules[0].smiles;
    const smilesB = molecules[1].smiles;
    if (!smilesA || !smilesB) return;

    let cancelled = false;
    setSimLoading(true);
    validationApi
      .getSimilarity(smilesA, smilesB)
      .then((res) => {
        if (!cancelled) {
          setSimilarity({
            score: res.tanimoto_similarity,
            common_bits: res.common_bits,
            bits_a: res.bits_a,
            bits_b: res.bits_b,
          });
        }
      })
      .catch(() => {
        if (!cancelled) setSimilarity(null);
      })
      .finally(() => {
        if (!cancelled) setSimLoading(false);
      });
    return () => { cancelled = true; };
  }, [hasTwoMolecules, molecules]);

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

              {/* Tanimoto Similarity Gauge */}
              {hasTwoMolecules && (
                <SimilarityGauge similarity={similarity} loading={simLoading} />
              )}

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
