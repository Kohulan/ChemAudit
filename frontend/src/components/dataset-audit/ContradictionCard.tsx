import { useState } from 'react';
import { motion, AnimatePresence } from 'framer-motion';
import { ChevronDown, AlertTriangle, Copy, Check } from 'lucide-react';
import type { ContradictoryLabelResult } from '../../types/dataset_intelligence';
import { MoleculeViewer } from '../molecules/MoleculeViewer';

// =============================================================================
// Types
// =============================================================================

interface ContradictionCardProps {
  contradiction: ContradictoryLabelResult;
  isExpanded: boolean;
  onToggle: () => void;
  index: number;
}

// =============================================================================
// Severity — uses project theme vars (primary = #c41e3a, accent = #d97706)
// =============================================================================

// All severity levels use the same neutral elevated background.
// Color is conveyed via left border, text, badge, and shadow tint only.
const CLAY_BASE = {
  clayBg: 'var(--color-surface-elevated)',
  clayShadow: `
    2px 4px 8px 0 rgba(0, 0, 0, 0.06),
    4px 8px 16px 0 rgba(0, 0, 0, 0.04),
    inset 2px 2px 6px rgba(255, 255, 255, 0.55),
    inset -2px -2px 6px rgba(0, 0, 0, 0.04)
  `,
  clayHoverShadow: `
    2px 4px 12px 0 rgba(0, 0, 0, 0.1),
    6px 12px 24px 0 rgba(0, 0, 0, 0.06),
    inset 2px 2px 6px rgba(255, 255, 255, 0.65),
    inset -2px -2px 6px rgba(0, 0, 0, 0.06)
  `,
};

function getSeverity(fold: number) {
  if (fold >= 100) {
    return {
      label: 'Critical',
      cssColor: 'var(--color-primary)',
      textClass: 'text-[var(--color-primary)]',
      badgeClass: 'bg-[var(--color-primary-lighter)] text-[var(--color-primary-dark)]',
      borderColor: 'var(--color-primary)',
      ...CLAY_BASE,
    };
  }
  if (fold >= 10) {
    return {
      label: 'Significant',
      cssColor: 'var(--color-accent)',
      textClass: 'text-[var(--color-accent-dark)]',
      badgeClass: 'bg-[var(--color-accent-lighter)] text-[var(--color-accent-dark)]',
      borderColor: 'var(--color-accent)',
      ...CLAY_BASE,
    };
  }
  return {
    label: 'Moderate',
    cssColor: 'var(--color-accent-light)',
    textClass: 'text-[var(--color-accent-dark)]',
    badgeClass: 'bg-[var(--color-accent-lighter)] text-[var(--color-accent-dark)]',
    borderColor: 'var(--color-accent-light)',
    ...CLAY_BASE,
  };
}

// =============================================================================
// Activity Range Bar — dots rendered outside overflow area to prevent clipping
// =============================================================================

function ActivityRangeBar({ entries, severity }: {
  entries: ContradictoryLabelResult['entries'];
  severity: ReturnType<typeof getSeverity>;
}) {
  const activities = entries
    .map((e) => (typeof e.activity === 'number' ? e.activity : null))
    .filter((v): v is number => v !== null);

  if (activities.length === 0) return null;

  const min = Math.min(...activities);
  const max = Math.max(...activities);

  // Use log scale for large ranges, linear for small
  const useLog = max / Math.max(min, 1e-10) > 10;

  const toPosition = (val: number): number => {
    if (!useLog) {
      const range = max - min || 1;
      return ((val - min) / range) * 100;
    }
    const logMin = Math.log10(Math.max(min, 1e-10));
    const logMax = Math.log10(Math.max(max, 1e-10));
    const logRange = logMax - logMin || 1;
    const logVal = Math.log10(Math.max(val, 1e-10));
    return ((logVal - logMin) / logRange) * 100;
  };

  return (
    <div className="mt-3">
      <div className="flex items-center justify-between text-[10px] text-[var(--color-text-muted)] mb-1.5">
        <span className="tabular-nums font-mono">{min.toPrecision(3)}</span>
        <span className="text-[var(--color-text-secondary)] font-medium">
          Activity Range{useLog ? ' (log scale)' : ''}
        </span>
        <span className="tabular-nums font-mono">{max.toPrecision(3)}</span>
      </div>
      {/* Track + dots share same relative parent. Dots use % left within 4-96% to avoid edge clipping. */}
      <div className="relative h-6">
        {/* Track background */}
        <div className="absolute left-0 right-0 top-1/2 -translate-y-1/2 h-2.5 rounded-full bg-[var(--color-surface-sunken)]">
          <motion.div
            className="absolute inset-y-0 left-0 rounded-full"
            style={{ backgroundColor: severity.cssColor, opacity: 0.15 }}
            initial={{ width: 0 }}
            animate={{ width: '100%' }}
            transition={{ duration: 0.6, ease: 'easeOut' }}
          />
        </div>
        {/* Dots — clamped to 4%-96% so they stay fully within the track */}
        {activities.map((val, i) => {
          const rawPct = toPosition(val);
          const clampedPct = 4 + (rawPct / 100) * 92; // maps 0-100% → 4-96%
          // Use negative margins for centering instead of transform:translate,
          // because Framer Motion's scale animation overrides the transform property.
          return (
            <motion.div
              key={i}
              className="absolute w-3.5 h-3.5 rounded-full"
              style={{
                top: '50%',
                left: `${clampedPct}%`,
                marginTop: -7,  // half of 14px (h-3.5 = 0.875rem ≈ 14px)
                marginLeft: -7, // half of 14px
                backgroundColor: severity.cssColor,
                border: '2.5px solid var(--color-surface-elevated)',
                boxShadow: 'var(--shadow-sm)',
              }}
              initial={{ scale: 0, opacity: 0 }}
              animate={{ scale: 1, opacity: 1 }}
              transition={{ delay: 0.3 + i * 0.12, duration: 0.35, type: 'spring', stiffness: 300 }}
              title={`Activity: ${val.toPrecision(4)}`}
            />
          );
        })}
      </div>
    </div>
  );
}

// =============================================================================
// Copy InChIKey button
// =============================================================================

function CopyInChIKey({ inchikey }: { inchikey: string }) {
  const [copied, setCopied] = useState(false);

  const handleCopy = (e: React.MouseEvent) => {
    e.stopPropagation();
    navigator.clipboard.writeText(inchikey);
    setCopied(true);
    setTimeout(() => setCopied(false), 1500);
  };

  return (
    <button
      onClick={handleCopy}
      className="p-0.5 rounded hover:bg-[var(--color-surface-sunken)] transition-colors"
      title="Copy InChIKey"
    >
      {copied ? (
        <Check className="w-3 h-3 text-emerald-500" />
      ) : (
        <Copy className="w-3 h-3 text-[var(--color-text-muted)]" />
      )}
    </button>
  );
}

// =============================================================================
// Component
// =============================================================================

export function ContradictionCard({
  contradiction,
  isExpanded,
  onToggle,
  index,
}: ContradictionCardProps) {
  const { inchikey, fold_difference, entry_count, smiles, entries } = contradiction;
  const severity = getSeverity(fold_difference);

  const activities = entries
    .map((e) => (typeof e.activity === 'number' ? e.activity : null))
    .filter((v): v is number => v !== null);
  const minActivity = activities.length > 0 ? Math.min(...activities) : null;
  const maxActivity = activities.length > 0 ? Math.max(...activities) : null;

  return (
    <motion.div
      initial={{ opacity: 0, y: 16 }}
      animate={{ opacity: 1, y: 0 }}
      transition={{ delay: index * 0.06, duration: 0.35, ease: [0.4, 0, 0.2, 1] }}
    >
      <div
        className="relative rounded-2xl cursor-pointer transition-all duration-300 hover:-translate-y-1 overflow-hidden"
        style={{
          backgroundColor: severity.clayBg,
          borderLeft: `4px solid ${severity.borderColor}`,
          boxShadow: severity.clayShadow,
        }}
        onMouseEnter={(e) => { e.currentTarget.style.boxShadow = severity.clayHoverShadow; }}
        onMouseLeave={(e) => { e.currentTarget.style.boxShadow = severity.clayShadow; }}
        onClick={onToggle}
        role="button"
        tabIndex={0}
        aria-expanded={isExpanded}
        aria-label={`${inchikey}: ${fold_difference.toFixed(1)}x fold-difference, ${entry_count} entries`}
        onKeyDown={(e) => { if (e.key === 'Enter' || e.key === ' ') { e.preventDefault(); onToggle(); } }}
      >
        <div className="p-5">
          <div className="flex items-start gap-5">
            {/* 2D structure — uses the proven MoleculeViewer component */}
            <div className="shrink-0 w-[130px] h-[100px] rounded-2xl overflow-hidden bg-[var(--color-surface-elevated)] shadow-[var(--shadow-sm)] border border-[var(--color-border)]/60">
              <MoleculeViewer
                smiles={smiles}
                width={130}
                height={100}
                className="[&>svg]:bg-transparent"
              />
            </div>

            <div className="flex-1 min-w-0">
              {/* Severity badge + fold value */}
              <div className="flex items-center gap-2.5 flex-wrap">
                <span className={`inline-flex items-center gap-1 px-2.5 py-1 text-xs font-semibold rounded-full shadow-sm ${severity.badgeClass}`}>
                  <AlertTriangle className="w-3 h-3" />
                  {severity.label}
                </span>
                <span className={`text-xl font-bold tabular-nums font-display ${severity.textClass}`}>
                  {fold_difference.toFixed(1)}x
                </span>
                <span className="text-xs text-[var(--color-text-muted)]">fold-difference</span>
              </div>

              {/* InChIKey with copy */}
              <div className="flex items-center gap-1.5 mt-2.5">
                <p className="text-xs font-mono text-[var(--color-text-secondary)] truncate" title={inchikey}>
                  {inchikey}
                </p>
                <CopyInChIKey inchikey={inchikey} />
              </div>

              {/* Stats row */}
              <div className="flex items-center gap-4 mt-3">
                <div className="flex items-center gap-1.5">
                  <span className="text-xs text-[var(--color-text-muted)]">Entries:</span>
                  <span className="text-sm font-semibold text-[var(--color-text-primary)] tabular-nums">{entry_count}</span>
                </div>
                {minActivity !== null && maxActivity !== null && (
                  <>
                    <div className="w-px h-4 bg-[var(--color-border)]" />
                    <div className="flex items-center gap-1.5">
                      <span className="text-xs text-[var(--color-text-muted)]">Min:</span>
                      <span className="text-xs font-mono font-semibold text-[var(--color-text-primary)] tabular-nums">{minActivity.toPrecision(3)}</span>
                    </div>
                    <div className="flex items-center gap-1.5">
                      <span className="text-xs text-[var(--color-text-muted)]">Max:</span>
                      <span className="text-xs font-mono font-semibold text-[var(--color-text-primary)] tabular-nums">{maxActivity.toPrecision(3)}</span>
                    </div>
                  </>
                )}
              </div>

              {/* Activity range bar */}
              <ActivityRangeBar entries={entries} severity={severity} />
            </div>
          </div>

          {/* Expand indicator */}
          <div className="flex items-center justify-center gap-1.5 mt-4 pt-3 border-t border-[var(--color-border)]/40">
            <motion.span
              animate={{ rotate: isExpanded ? 180 : 0 }}
              transition={{ duration: 0.25, ease: 'easeOut' }}
            >
              <ChevronDown className="w-4 h-4 text-[var(--color-text-muted)]" />
            </motion.span>
            <span className="text-xs font-medium text-[var(--color-text-muted)]">
              {isExpanded ? 'Hide entries' : 'Show all entries'}
            </span>
          </div>
        </div>

        {/* Expandable entries */}
        <AnimatePresence initial={false}>
          {isExpanded && (
            <motion.div
              initial={{ height: 0, opacity: 0 }}
              animate={{ height: 'auto', opacity: 1 }}
              exit={{ height: 0, opacity: 0 }}
              transition={{ duration: 0.3, ease: [0.4, 0, 0.2, 1] }}
              className="overflow-hidden"
            >
              <div className="px-5 pb-5">
                <table className="w-full text-xs">
                  <thead>
                    <tr className="text-[var(--color-text-secondary)]">
                      <th className="text-left py-2 pr-3 font-semibold text-[10px] uppercase tracking-wider">Row</th>
                      <th className="text-left py-2 pr-3 font-semibold text-[10px] uppercase tracking-wider">SMILES</th>
                      <th className="text-right py-2 pr-3 font-semibold text-[10px] uppercase tracking-wider">Activity</th>
                      <th className="text-right py-2 font-semibold text-[10px] uppercase tracking-wider w-[90px]">Relative</th>
                    </tr>
                  </thead>
                  <tbody>
                    {entries.map((entry, i) => {
                      const val = typeof entry.activity === 'number' ? entry.activity : 0;
                      const barPct = maxActivity && maxActivity > 0 ? (val / maxActivity) * 100 : 0;
                      return (
                        <motion.tr
                          key={entry.row_index}
                          className="border-t border-[var(--color-border)]/30"
                          initial={{ opacity: 0, x: -8 }}
                          animate={{ opacity: 1, x: 0 }}
                          transition={{ delay: i * 0.06, duration: 0.25 }}
                        >
                          <td className="py-2.5 pr-3 text-[var(--color-text-muted)] tabular-nums">
                            {entry.row_index}
                          </td>
                          <td className="py-2.5 pr-3 font-mono text-[var(--color-text-primary)] truncate max-w-[280px]" title={entry.smiles}>
                            {entry.smiles}
                          </td>
                          <td className="py-2.5 pr-3 text-right font-mono text-[var(--color-text-primary)] tabular-nums font-semibold">
                            {typeof entry.activity === 'number' ? entry.activity.toPrecision(4) : String(entry.activity)}
                          </td>
                          <td className="py-2.5">
                            <div className="h-2.5 rounded-full bg-[var(--color-surface-sunken)] overflow-hidden">
                              <motion.div
                                className="h-full rounded-full"
                                style={{ backgroundColor: severity.cssColor }}
                                initial={{ width: 0 }}
                                animate={{ width: `${barPct}%` }}
                                transition={{ delay: 0.15 + i * 0.06, duration: 0.4, ease: 'easeOut' }}
                              />
                            </div>
                          </td>
                        </motion.tr>
                      );
                    })}
                  </tbody>
                </table>
              </div>
            </motion.div>
          )}
        </AnimatePresence>
      </div>
    </motion.div>
  );
}
