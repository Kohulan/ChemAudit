/**
 * Slide-over panel for inspecting and scoring selected molecules inline.
 *
 * Two tabs:
 * 1. Validation — shows existing batch results for selected molecules (instant)
 * 2. Score — apply a profile to selected molecules via inline API (no Celery)
 *
 * Results appear directly in the panel with a download option.
 */

import { useState, useEffect, useCallback, useMemo } from 'react';
import { motion, AnimatePresence } from 'framer-motion';
import {
  X, Eye, BarChart3, Download, ChevronDown, CheckCircle2,
  XCircle, AlertTriangle, Beaker,
} from 'lucide-react';
import { ClayButton } from '../ui/ClayButton';
import { subsetApi, profilesApi, batchApi } from '../../services/api';
import type { InlineScoreResponse, InlineScoredMolecule } from '../../services/api';
import { cn } from '../../lib/utils';
import type { ScoringProfile } from '../../types/workflow';
import type { BatchResult } from '../../types/batch';

/* ------------------------------------------------------------------ */
/*  Types & helpers                                                    */
/* ------------------------------------------------------------------ */

type TabId = 'validation' | 'scoring';

interface SubsetActionPanelProps {
  jobId: string;
  selectedIndices: Set<number>;
  isOpen: boolean;
  onClose: () => void;
}

const PROPERTY_LABELS: Record<string, string> = {
  mw: 'MW',
  logp: 'LogP',
  hbd: 'HBD',
  hba: 'HBA',
  tpsa: 'TPSA',
  rotatable_bonds: 'RotB',
  aromatic_rings: 'AroR',
  fsp3: 'Fsp3',
};

function scoreColor(score: number): string {
  if (score >= 80) return 'text-emerald-600 dark:text-emerald-400';
  if (score >= 50) return 'text-amber-600 dark:text-amber-400';
  return 'text-red-500 dark:text-red-400';
}

function scoreBg(score: number): string {
  if (score >= 80) return 'bg-emerald-500/15 border-emerald-500/25';
  if (score >= 50) return 'bg-amber-500/15 border-amber-500/25';
  return 'bg-red-500/15 border-red-500/25';
}

function barColor(score: number): string {
  if (score >= 80) return 'bg-emerald-500';
  if (score >= 50) return 'bg-amber-500';
  return 'bg-red-500';
}

/** Generate a CSV blob from scored molecules for download. */
function buildScoreCsv(data: InlineScoreResponse): Blob {
  const propKeys = data.molecules.length > 0
    ? Object.keys(data.molecules[0].profile.properties)
    : [];
  const header = ['Index', 'Name', 'SMILES', 'Profile Score',
    ...propKeys.map(k => `${PROPERTY_LABELS[k] || k} Value`),
    ...propKeys.map(k => `${PROPERTY_LABELS[k] || k} In Range`),
  ].join(',');

  const rows = data.molecules.map(m => {
    const p = m.profile;
    return [
      m.index,
      `"${(m.name || '').replace(/"/g, '""')}"`,
      `"${m.smiles.replace(/"/g, '""')}"`,
      p.score ?? '',
      ...propKeys.map(k => p.properties[k]?.value ?? ''),
      ...propKeys.map(k => p.properties[k]?.in_range ?? ''),
    ].join(',');
  });

  return new Blob([header + '\n' + rows.join('\n')], { type: 'text/csv' });
}

/** Generate a CSV blob from validation batch results for download. */
function buildValidationCsv(results: BatchResult[]): Blob {
  const header = 'Index,Name,SMILES,Status,Validation Score,QED,Issues';
  const rows = results.map(r => [
    r.index,
    `"${(r.name || '').replace(/"/g, '""')}"`,
    `"${r.smiles.replace(/"/g, '""')}"`,
    r.status,
    r.validation?.overall_score ?? '',
    r.scoring?.druglikeness?.qed_score ?? '',
    (r.validation?.issues?.filter(i => !i.passed).length ?? 0),
  ].join(','));

  return new Blob([header + '\n' + rows.join('\n')], { type: 'text/csv' });
}

function triggerDownload(blob: Blob, filename: string) {
  const url = URL.createObjectURL(blob);
  const a = document.createElement('a');
  a.href = url;
  a.download = filename;
  a.click();
  URL.revokeObjectURL(url);
}

/* ------------------------------------------------------------------ */
/*  Stagger animation                                                  */
/* ------------------------------------------------------------------ */

const stagger = {
  container: {
    hidden: { opacity: 0 },
    show: { opacity: 1, transition: { staggerChildren: 0.04 } },
  },
  item: {
    hidden: { opacity: 0, y: 8 },
    show: { opacity: 1, y: 0, transition: { duration: 0.25, ease: 'easeOut' } },
  },
};

/* ------------------------------------------------------------------ */
/*  Mini Score Distribution Bar                                        */
/* ------------------------------------------------------------------ */

function MiniDistribution({ molecules }: { molecules: InlineScoredMolecule[] }) {
  const buckets = useMemo(() => {
    const b = { high: 0, med: 0, low: 0, none: 0 };
    for (const m of molecules) {
      const s = m.profile.score;
      if (s == null) b.none++;
      else if (s >= 80) b.high++;
      else if (s >= 50) b.med++;
      else b.low++;
    }
    return b;
  }, [molecules]);

  const total = molecules.length || 1;

  return (
    <div className="space-y-1.5">
      <div className="flex items-center gap-1.5 h-5 rounded-full overflow-hidden bg-[var(--color-surface-sunken)]">
        {buckets.high > 0 && (
          <motion.div
            initial={{ width: 0 }}
            animate={{ width: `${(buckets.high / total) * 100}%` }}
            transition={{ duration: 0.5, ease: 'easeOut' }}
            className="h-full bg-emerald-500 rounded-l-full first:rounded-l-full"
          />
        )}
        {buckets.med > 0 && (
          <motion.div
            initial={{ width: 0 }}
            animate={{ width: `${(buckets.med / total) * 100}%` }}
            transition={{ duration: 0.5, delay: 0.1, ease: 'easeOut' }}
            className="h-full bg-amber-500"
          />
        )}
        {buckets.low > 0 && (
          <motion.div
            initial={{ width: 0 }}
            animate={{ width: `${(buckets.low / total) * 100}%` }}
            transition={{ duration: 0.5, delay: 0.2, ease: 'easeOut' }}
            className="h-full bg-red-500 last:rounded-r-full"
          />
        )}
      </div>
      <div className="flex items-center gap-3 text-[10px] text-[var(--color-text-muted)]">
        {buckets.high > 0 && (
          <span className="flex items-center gap-1">
            <span className="w-2 h-2 rounded-full bg-emerald-500" />
            {buckets.high} pass
          </span>
        )}
        {buckets.med > 0 && (
          <span className="flex items-center gap-1">
            <span className="w-2 h-2 rounded-full bg-amber-500" />
            {buckets.med} moderate
          </span>
        )}
        {buckets.low > 0 && (
          <span className="flex items-center gap-1">
            <span className="w-2 h-2 rounded-full bg-red-500" />
            {buckets.low} fail
          </span>
        )}
      </div>
    </div>
  );
}

/* ------------------------------------------------------------------ */
/*  Validation Result Card                                             */
/* ------------------------------------------------------------------ */

function ValidationCard({ result }: { result: BatchResult }) {
  const score = result.validation?.overall_score ?? 0;
  const issues = result.validation?.issues?.filter(i => !i.passed) ?? [];
  const warnings = result.validation?.issues?.filter(i => i.passed && i.severity === 'warning') ?? [];
  const isError = result.status === 'error';

  return (
    <motion.div
      variants={stagger.item}
      className={cn(
        'rounded-xl border p-3 transition-colors',
        isError
          ? 'bg-red-500/5 border-red-500/15'
          : issues.length > 0
            ? 'bg-amber-500/5 border-amber-500/15'
            : 'bg-emerald-500/5 border-emerald-500/15',
      )}
    >
      <div className="flex items-start justify-between gap-2 mb-1.5">
        <div className="flex items-center gap-1.5 min-w-0">
          {isError ? (
            <XCircle className="w-3.5 h-3.5 text-red-500 flex-shrink-0" />
          ) : issues.length > 0 ? (
            <AlertTriangle className="w-3.5 h-3.5 text-amber-500 flex-shrink-0" />
          ) : (
            <CheckCircle2 className="w-3.5 h-3.5 text-emerald-500 flex-shrink-0" />
          )}
          <span className="text-xs font-medium text-[var(--color-text-primary)] truncate">
            {result.name || `Molecule #${result.index}`}
          </span>
        </div>
        <span className={cn('text-xs font-semibold tabular-nums flex-shrink-0', scoreColor(score))}>
          {score.toFixed(0)}
        </span>
      </div>

      {/* Key metrics row */}
      <div className="flex items-center gap-2 text-[10px] text-[var(--color-text-muted)] mb-1">
        {result.scoring?.druglikeness && (
          <>
            <span>MW {result.scoring.druglikeness.mw?.toFixed(0)}</span>
            <span className="opacity-30">|</span>
            <span>LogP {result.scoring.druglikeness.logp?.toFixed(1)}</span>
            <span className="opacity-30">|</span>
            <span>QED {result.scoring.druglikeness.qed_score?.toFixed(2)}</span>
          </>
        )}
      </div>

      {/* Issues */}
      {issues.length > 0 && (
        <div className="flex flex-wrap gap-1 mt-1.5">
          {issues.slice(0, 3).map((issue, i) => (
            <span
              key={i}
              className="inline-flex items-center gap-0.5 px-1.5 py-0.5 rounded-md text-[9px] font-medium bg-red-500/10 text-red-600 dark:text-red-400"
            >
              {issue.check_name}
            </span>
          ))}
          {issues.length > 3 && (
            <span className="text-[9px] text-[var(--color-text-muted)]">
              +{issues.length - 3} more
            </span>
          )}
        </div>
      )}
      {warnings.length > 0 && issues.length === 0 && (
        <span className="text-[9px] text-amber-600 dark:text-amber-400">
          {warnings.length} warning{warnings.length > 1 ? 's' : ''}
        </span>
      )}
      {isError && result.error && (
        <p className="text-[9px] text-red-500 mt-1 line-clamp-2">{result.error}</p>
      )}
    </motion.div>
  );
}

/* ------------------------------------------------------------------ */
/*  Scored Molecule Card                                               */
/* ------------------------------------------------------------------ */

function ScoredCard({ mol }: { mol: InlineScoredMolecule }) {
  const score = mol.profile.score;
  const props = mol.profile.properties;
  const propEntries = Object.entries(props);

  return (
    <motion.div
      variants={stagger.item}
      className={cn(
        'rounded-xl border p-3',
        score != null ? scoreBg(score) : 'bg-[var(--color-surface-sunken)] border-[var(--color-border)]',
      )}
    >
      <div className="flex items-center justify-between gap-2 mb-2">
        <span className="text-xs font-medium text-[var(--color-text-primary)] truncate">
          {mol.name || `Molecule #${mol.index}`}
        </span>
        {score != null ? (
          <span className={cn('text-sm font-bold tabular-nums', scoreColor(score))}>
            {score.toFixed(0)}
          </span>
        ) : (
          <span className="text-xs text-[var(--color-text-muted)]">N/A</span>
        )}
      </div>

      {/* Score bar */}
      {score != null && (
        <div className="h-1.5 rounded-full bg-black/5 dark:bg-white/5 overflow-hidden mb-2">
          <motion.div
            initial={{ width: 0 }}
            animate={{ width: `${Math.min(score, 100)}%` }}
            transition={{ duration: 0.6, ease: 'easeOut' }}
            className={cn('h-full rounded-full', barColor(score))}
          />
        </div>
      )}

      {/* Property pills */}
      <div className="flex flex-wrap gap-1">
        {propEntries.map(([key, p]) => (
          <span
            key={key}
            className={cn(
              'inline-flex items-center gap-0.5 px-1.5 py-0.5 rounded-md text-[9px] font-medium',
              p.in_range
                ? 'bg-emerald-500/10 text-emerald-700 dark:text-emerald-400'
                : 'bg-red-500/10 text-red-600 dark:text-red-400',
            )}
          >
            {p.in_range ? (
              <CheckCircle2 className="w-2.5 h-2.5" />
            ) : (
              <XCircle className="w-2.5 h-2.5" />
            )}
            {PROPERTY_LABELS[key] || key}
          </span>
        ))}
      </div>
    </motion.div>
  );
}

/* ------------------------------------------------------------------ */
/*  Shimmer placeholder                                                */
/* ------------------------------------------------------------------ */

function ScoringShimmer() {
  return (
    <div className="space-y-3 animate-pulse">
      {[...Array(3)].map((_, i) => (
        <div key={i} className="rounded-xl border border-[var(--color-border)] p-3">
          <div className="flex justify-between mb-2">
            <div className="h-3 w-24 rounded bg-[var(--color-surface-sunken)]" />
            <div className="h-4 w-8 rounded bg-[var(--color-surface-sunken)]" />
          </div>
          <div className="h-1.5 rounded-full bg-[var(--color-surface-sunken)] mb-2" />
          <div className="flex gap-1">
            {[...Array(4)].map((_, j) => (
              <div key={j} className="h-4 w-10 rounded-md bg-[var(--color-surface-sunken)]" />
            ))}
          </div>
        </div>
      ))}
    </div>
  );
}

/* ------------------------------------------------------------------ */
/*  Main Component                                                     */
/* ------------------------------------------------------------------ */

export function SubsetActionPanel({
  jobId,
  selectedIndices,
  isOpen,
  onClose,
}: SubsetActionPanelProps) {
  const [tab, setTab] = useState<TabId>('validation');
  const [profiles, setProfiles] = useState<ScoringProfile[]>([]);
  const [selectedProfileId, setSelectedProfileId] = useState<number | undefined>();
  const [isScoring, setIsScoring] = useState(false);
  const [scoreResult, setScoreResult] = useState<InlineScoreResponse | null>(null);
  const [scoreError, setScoreError] = useState<string | null>(null);
  const [selectedResults, setSelectedResults] = useState<BatchResult[]>([]);
  const [loadingResults, setLoadingResults] = useState(false);

  // Fetch results for selected indices across all pages when panel opens
  useEffect(() => {
    if (!isOpen || selectedIndices.size === 0) return;
    let cancelled = false;
    setLoadingResults(true);
    setSelectedResults([]);
    const idxSet = selectedIndices;
    const found: BatchResult[] = [];

    (async () => {
      try {
        // Paginate through results (max page_size=100) until we find all selected
        let page = 1;
        const pageSize = 100;
        // eslint-disable-next-line no-constant-condition
        while (true) {
          const res = await batchApi.getBatchResults(jobId, page, pageSize);
          if (cancelled) return;
          for (const r of res.results) {
            if (idxSet.has(r.index)) found.push(r);
          }
          // Stop if we found all or reached the last page
          if (found.length >= idxSet.size || page >= res.total_pages) break;
          page++;
        }
        if (!cancelled) setSelectedResults(found);
      } catch {
        // ignore
      } finally {
        if (!cancelled) setLoadingResults(false);
      }
    })();

    return () => { cancelled = true; };
  }, [isOpen, jobId, selectedIndices]);

  // Load profiles when panel opens
  useEffect(() => {
    if (isOpen) {
      profilesApi.getProfiles().then(setProfiles).catch(() => {});
    }
  }, [isOpen]);

  // Reset scoring state when selection changes
  useEffect(() => {
    setScoreResult(null);
    setScoreError(null);
  }, [selectedIndices]);

  // Close on Escape
  useEffect(() => {
    const handleKey = (e: KeyboardEvent) => {
      if (e.key === 'Escape') onClose();
    };
    if (isOpen) {
      document.addEventListener('keydown', handleKey);
      return () => document.removeEventListener('keydown', handleKey);
    }
  }, [isOpen, onClose]);

  // Summary stats for validation tab
  const validationSummary = useMemo(() => {
    let pass = 0, warn = 0, fail = 0;
    for (const r of selectedResults) {
      if (r.status === 'error') { fail++; continue; }
      const issues = r.validation?.issues?.filter(i => !i.passed) ?? [];
      if (issues.length > 0) warn++;
      else pass++;
    }
    return { pass, warn, fail, total: selectedResults.length };
  }, [selectedResults]);

  const handleScore = useCallback(async () => {
    if (!selectedProfileId) return;
    setIsScoring(true);
    setScoreError(null);
    setScoreResult(null);
    try {
      const indices = Array.from(selectedIndices);
      const res = await subsetApi.scoreInline(jobId, indices, selectedProfileId);
      setScoreResult(res);
    } catch (err) {
      setScoreError(err instanceof Error ? err.message : 'Scoring failed');
    } finally {
      setIsScoring(false);
    }
  }, [jobId, selectedIndices, selectedProfileId]);

  const handleDownloadValidation = useCallback(() => {
    const blob = buildValidationCsv(selectedResults);
    triggerDownload(blob, `validation_subset_${jobId.slice(0, 8)}.csv`);
  }, [selectedResults, jobId]);

  const handleDownloadScores = useCallback(() => {
    if (!scoreResult) return;
    const blob = buildScoreCsv(scoreResult);
    triggerDownload(blob, `profile_scores_${jobId.slice(0, 8)}.csv`);
  }, [scoreResult, jobId]);

  if (!isOpen) return null;

  return (
    <div className="fixed inset-0 z-50 flex justify-end pt-14" onClick={onClose}>
      {/* Backdrop */}
      <motion.div
        initial={{ opacity: 0 }}
        animate={{ opacity: 1 }}
        exit={{ opacity: 0 }}
        className="absolute inset-0 bg-black/30 backdrop-blur-[2px]"
      />

      {/* Panel */}
      <motion.div
        initial={{ x: '100%' }}
        animate={{ x: 0 }}
        exit={{ x: '100%' }}
        transition={{ type: 'spring', damping: 30, stiffness: 300 }}
        className="relative w-full max-w-md bg-[var(--color-surface-elevated)] border-l border-[var(--color-border)] shadow-2xl flex flex-col h-full"
        onClick={(e) => e.stopPropagation()}
      >
        {/* ── Header ── */}
        <div className="flex items-center justify-between px-5 py-4 border-b border-[var(--color-border)]">
          <div className="flex items-center gap-2.5">
            <div className="w-8 h-8 rounded-lg bg-gradient-to-br from-rose-500/20 to-amber-500/20 flex items-center justify-center">
              <Beaker className="w-4 h-4 text-[var(--color-primary)]" />
            </div>
            <div>
              <h3 className="text-sm font-semibold text-[var(--color-text-primary)] font-display">
                Inspect Selection
              </h3>
              <p className="text-[10px] text-[var(--color-text-muted)] mt-0.5">
                {selectedIndices.size} molecule{selectedIndices.size !== 1 ? 's' : ''} selected
              </p>
            </div>
          </div>
          <button
            onClick={onClose}
            className="p-1.5 rounded-lg hover:bg-[var(--color-surface-sunken)] transition-colors"
          >
            <X className="w-4 h-4 text-[var(--color-text-muted)]" />
          </button>
        </div>

        {/* ── Tab Switcher ── */}
        <div className="flex border-b border-[var(--color-border)]">
          {([
            { id: 'validation' as TabId, label: 'Validation', icon: Eye },
            { id: 'scoring' as TabId, label: 'Profile Score', icon: BarChart3 },
          ]).map(({ id, label, icon: Icon }) => (
            <button
              key={id}
              onClick={() => setTab(id)}
              className={cn(
                'flex-1 flex items-center justify-center gap-1.5 px-4 py-2.5 text-xs font-medium transition-all relative',
                tab === id
                  ? 'text-[var(--color-text-primary)]'
                  : 'text-[var(--color-text-muted)] hover:text-[var(--color-text-secondary)]',
              )}
            >
              <Icon className="w-3.5 h-3.5" />
              {label}
              {tab === id && (
                <motion.div
                  layoutId="subset-tab-indicator"
                  className="absolute bottom-0 left-2 right-2 h-0.5 rounded-full bg-[var(--color-primary)]"
                  transition={{ type: 'spring', stiffness: 400, damping: 30 }}
                />
              )}
            </button>
          ))}
        </div>

        {/* ── Content ── */}
        <div className="relative flex-1 overflow-y-auto">
          <AnimatePresence mode="wait">
            {tab === 'validation' ? (
              <motion.div
                key="validation"
                initial={{ opacity: 0, x: -10 }}
                animate={{ opacity: 1, x: 0 }}
                exit={{ opacity: 0, x: -10 }}
                transition={{ duration: 0.15 }}
                className="p-4 space-y-3"
              >
                {/* Summary strip */}
                <div className="flex items-center gap-2 p-2.5 rounded-xl bg-[var(--color-surface-sunken)]">
                  <div className="flex items-center gap-3 text-[10px] font-medium">
                    {validationSummary.pass > 0 && (
                      <span className="flex items-center gap-1 text-emerald-600 dark:text-emerald-400">
                        <CheckCircle2 className="w-3 h-3" />
                        {validationSummary.pass} pass
                      </span>
                    )}
                    {validationSummary.warn > 0 && (
                      <span className="flex items-center gap-1 text-amber-600 dark:text-amber-400">
                        <AlertTriangle className="w-3 h-3" />
                        {validationSummary.warn} issues
                      </span>
                    )}
                    {validationSummary.fail > 0 && (
                      <span className="flex items-center gap-1 text-red-500">
                        <XCircle className="w-3 h-3" />
                        {validationSummary.fail} errors
                      </span>
                    )}
                  </div>
                </div>

                {loadingResults && (
                  <div className="text-xs text-[var(--color-text-muted)] text-center py-4 animate-pulse">
                    Loading molecule data...
                  </div>
                )}

                {/* Molecule cards */}
                <motion.div
                  variants={stagger.container}
                  initial="hidden"
                  animate="show"
                  className="space-y-2"
                >
                  {selectedResults.map(r => (
                    <ValidationCard key={r.index} result={r} />
                  ))}
                </motion.div>

                {selectedResults.length === 0 && !loadingResults && (
                  <p className="text-xs text-[var(--color-text-muted)] text-center py-6">
                    No results found for selected molecules.
                  </p>
                )}

                {/* Download */}
                {selectedResults.length > 0 && (
                  <div className="pt-2">
                    <ClayButton
                      size="sm"
                      variant="default"
                      onClick={handleDownloadValidation}
                      leftIcon={<Download className="w-3.5 h-3.5" />}
                      className="w-full"
                    >
                      Download Validation CSV
                    </ClayButton>
                  </div>
                )}
              </motion.div>
            ) : (
              <motion.div
                key="scoring"
                initial={{ opacity: 0, x: 10 }}
                animate={{ opacity: 1, x: 0 }}
                exit={{ opacity: 0, x: 10 }}
                transition={{ duration: 0.15 }}
                className="p-4 space-y-4"
              >
                {/* Profile picker */}
                <div className="space-y-2">
                  <label className="text-[10px] font-medium uppercase tracking-wider text-[var(--color-text-muted)]">
                    Scoring Profile
                  </label>
                  <div className="relative">
                    <select
                      value={selectedProfileId ?? ''}
                      onChange={(e) => {
                        setSelectedProfileId(e.target.value ? Number(e.target.value) : undefined);
                        setScoreResult(null);
                        setScoreError(null);
                      }}
                      className={cn(
                        'w-full px-3 py-2.5 rounded-xl text-sm appearance-none',
                        'bg-[var(--color-surface-sunken)] border border-[var(--color-border)]',
                        'text-[var(--color-text-primary)]',
                        'focus:outline-none focus:ring-2 focus:ring-[var(--color-primary)]/30',
                      )}
                    >
                      <option value="">Select a profile...</option>
                      {profiles.map((p) => (
                        <option key={p.id} value={p.id}>
                          {p.name} {p.is_preset ? '(Preset)' : ''}
                        </option>
                      ))}
                    </select>
                    <ChevronDown className="absolute right-3 top-1/2 -translate-y-1/2 w-4 h-4 text-[var(--color-text-muted)] pointer-events-none" />
                  </div>
                  <ClayButton
                    variant="accent"
                    size="sm"
                    onClick={handleScore}
                    disabled={!selectedProfileId || isScoring}
                    loading={isScoring}
                    leftIcon={<BarChart3 className="w-3.5 h-3.5" />}
                    className="w-full"
                  >
                    {isScoring ? 'Scoring...' : `Score ${selectedIndices.size} Molecules`}
                  </ClayButton>
                </div>

                {/* Loading shimmer */}
                {isScoring && <ScoringShimmer />}

                {/* Error */}
                {scoreError && (
                  <motion.div
                    initial={{ opacity: 0, y: 8 }}
                    animate={{ opacity: 1, y: 0 }}
                    className="p-3 rounded-xl bg-red-500/10 border border-red-500/20 text-xs text-red-600 dark:text-red-400"
                  >
                    {scoreError}
                  </motion.div>
                )}

                {/* Results */}
                <AnimatePresence>
                  {scoreResult && (
                    <motion.div
                      initial={{ opacity: 0, y: 12 }}
                      animate={{ opacity: 1, y: 0 }}
                      exit={{ opacity: 0, y: -8 }}
                      transition={{ duration: 0.3 }}
                      className="space-y-3"
                    >
                      {/* Summary header */}
                      <div className="flex items-center justify-between">
                        <div>
                          <h4 className="text-xs font-semibold text-[var(--color-text-primary)]">
                            {scoreResult.profile_name}
                          </h4>
                          <p className="text-[10px] text-[var(--color-text-muted)]">
                            {scoreResult.molecules.length} molecules scored
                          </p>
                        </div>
                        {(() => {
                          const scores = scoreResult.molecules
                            .map(m => m.profile.score)
                            .filter((s): s is number => s != null);
                          if (scores.length === 0) return null;
                          const avg = scores.reduce((a, b) => a + b, 0) / scores.length;
                          return (
                            <div className="text-right">
                              <div className={cn('text-lg font-bold tabular-nums', scoreColor(avg))}>
                                {avg.toFixed(0)}
                              </div>
                              <div className="text-[9px] text-[var(--color-text-muted)]">avg score</div>
                            </div>
                          );
                        })()}
                      </div>

                      {/* Distribution bar */}
                      <MiniDistribution molecules={scoreResult.molecules} />

                      {/* Scored cards */}
                      <motion.div
                        variants={stagger.container}
                        initial="hidden"
                        animate="show"
                        className="space-y-2"
                      >
                        {scoreResult.molecules
                          .sort((a, b) => (b.profile.score ?? 0) - (a.profile.score ?? 0))
                          .map(mol => (
                            <ScoredCard key={mol.index} mol={mol} />
                          ))}
                      </motion.div>

                      {/* Download */}
                      <div className="pt-1">
                        <ClayButton
                          size="sm"
                          variant="default"
                          onClick={handleDownloadScores}
                          leftIcon={<Download className="w-3.5 h-3.5" />}
                          className="w-full"
                        >
                          Download Scores CSV
                        </ClayButton>
                      </div>
                    </motion.div>
                  )}
                </AnimatePresence>

                {/* Empty state */}
                {!isScoring && !scoreResult && !scoreError && (
                  <div className="text-center py-8">
                    <BarChart3 className="w-8 h-8 mx-auto text-[var(--color-text-muted)] opacity-30 mb-2" />
                    <p className="text-xs text-[var(--color-text-muted)]">
                      Select a profile and click Score to see results.
                    </p>
                  </div>
                )}
              </motion.div>
            )}
          </AnimatePresence>
        </div>
      </motion.div>
    </div>
  );
}
