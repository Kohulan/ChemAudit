import { motion, AnimatePresence } from 'framer-motion';
import { CheckCircle2, AlertTriangle, XCircle, HelpCircle } from 'lucide-react';
import { safeHref } from '../../lib/sanitize';
import { MoleculeViewer } from '../molecules/MoleculeViewer';
import { CopyButton } from '../ui/CopyButton';
import type { ConsistencyResult, DatabaseEntry, PropertyComparison } from '../../types/integrations';

interface DatabaseComparisonPanelProps {
  result: ConsistencyResult;
}

/* ── Verdict config ────────────────────────────────── */

const VERDICT_CONFIG: Record<string, {
  gradient: string;
  ring: string;
  label: string;
  Icon: typeof CheckCircle2;
}> = {
  consistent: {
    gradient: 'from-emerald-500 via-teal-500 to-cyan-500',
    ring: 'ring-emerald-500/30',
    label: 'Consistent',
    Icon: CheckCircle2,
  },
  minor_differences: {
    gradient: 'from-amber-500 via-orange-500 to-yellow-500',
    ring: 'ring-amber-500/30',
    label: 'Minor Differences',
    Icon: AlertTriangle,
  },
  major_discrepancies: {
    gradient: 'from-rose-500 via-red-500 to-pink-500',
    ring: 'ring-rose-500/30',
    label: 'Major Discrepancies',
    Icon: XCircle,
  },
  no_data: {
    gradient: 'from-slate-400 via-gray-400 to-zinc-400',
    ring: 'ring-gray-400/30',
    label: 'No Data',
    Icon: HelpCircle,
  },
};

/* ── Status styling ────────────────────────────────── */

const STATUS_STYLES: Record<string, {
  badge: string;
  badgeText: string;
  accent: string;
  detailBg: string;
  icon: string;
  glow: string;
}> = {
  match: {
    badge: 'bg-emerald-500/10 text-emerald-700 border border-emerald-500/20',
    badgeText: 'Match',
    accent: 'bg-emerald-500',
    detailBg: 'bg-emerald-50/60 border-emerald-200/50',
    icon: '●',
    glow: 'shadow-emerald-500/5',
  },
  minor_difference: {
    badge: 'bg-amber-500/10 text-amber-700 border border-amber-500/20',
    badgeText: 'Minor Difference',
    accent: 'bg-amber-500',
    detailBg: 'bg-amber-50/60 border-amber-200/50',
    icon: '◐',
    glow: 'shadow-amber-500/5',
  },
  mismatch: {
    badge: 'bg-rose-500/10 text-rose-700 border border-rose-500/20',
    badgeText: 'Mismatch',
    accent: 'bg-rose-500',
    detailBg: 'bg-rose-50/60 border-rose-200/50',
    icon: '○',
    glow: 'shadow-rose-500/5',
  },
  missing: {
    badge: 'bg-gray-500/10 text-gray-500 border border-gray-300/30',
    badgeText: 'Missing',
    accent: 'bg-gray-400',
    detailBg: 'bg-gray-50/60 border-gray-200/50',
    icon: '—',
    glow: '',
  },
};

const PROP_LABELS: Record<string, string> = {
  canonical_smiles: 'SMILES',
  inchikey: 'InChIKey',
  inchi: 'InChI',
};

/* ── Database colors ─────────────────────────────── */

const DB_COLORS: Record<string, { dot: string; text: string; bg: string; border: string }> = {
  PubChem:  { dot: 'bg-blue-500',    text: 'text-blue-700',    bg: 'bg-blue-50',    border: 'border-blue-200' },
  ChEMBL:   { dot: 'bg-violet-500',  text: 'text-violet-700',  bg: 'bg-violet-50',  border: 'border-violet-200' },
  COCONUT:  { dot: 'bg-amber-500',   text: 'text-amber-700',   bg: 'bg-amber-50',   border: 'border-amber-200' },
  Wikidata: { dot: 'bg-sky-500',     text: 'text-sky-700',     bg: 'bg-sky-50',     border: 'border-sky-200' },
  Resolved: { dot: 'bg-emerald-500', text: 'text-emerald-700', bg: 'bg-emerald-50', border: 'border-emerald-200' },
};

const fallbackDB = { dot: 'bg-indigo-500', text: 'text-indigo-700', bg: 'bg-indigo-50', border: 'border-indigo-200' };

/* ── Sub-components ────────────────────────────────── */

function StructureCard({ entry }: { entry: DatabaseEntry & { canonical_smiles: string } }) {
  const c = DB_COLORS[entry.database] || fallbackDB;
  // For Resolved, prefer kekulized SMILES for display
  const displaySmiles = (entry.database === 'Resolved' && entry.kekulized_smiles) || entry.canonical_smiles;

  return (
    <motion.div
      initial={{ opacity: 0, scale: 0.95 }}
      animate={{ opacity: 1, scale: 1 }}
      whileHover={{ y: -2 }}
      className="group relative border border-gray-200/80 rounded-xl overflow-hidden bg-white hover:shadow-xl hover:shadow-gray-200/50 transition-all duration-300 cursor-pointer"
    >
      <div className={`h-1 ${c.dot}`} />
      <div className="p-2.5">
        <div className="flex items-center justify-between mb-1.5">
          <div className="flex items-center gap-1.5">
            <span className={`w-1.5 h-1.5 rounded-full ${c.dot}`} />
            <span className={`text-[10px] font-bold ${c.text} uppercase tracking-widest`}>
              {entry.database}
            </span>
          </div>
          {entry.url && (
            <a href={safeHref(entry.url)} target="_blank" rel="noopener noreferrer"
              className="text-gray-400 hover:text-gray-700 transition-colors">
              <svg className="w-3 h-3" fill="none" stroke="currentColor" viewBox="0 0 24 24">
                <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={2} d="M10 6H6a2 2 0 00-2 2v10a2 2 0 002 2h10a2 2 0 002-2v-4M14 4h6m0 0v6m0-6L10 14" />
              </svg>
            </a>
          )}
        </div>
        <MoleculeViewer smiles={displaySmiles} width={180} height={130} className="bg-white rounded-lg" />
      </div>
    </motion.div>
  );
}

/** Resolved SMILES row — shows canonical kekulized form */
function ResolvedSmilesRow({ kekulized, canonical }: { kekulized?: string | null; canonical: string }) {
  const display = kekulized || canonical;
  return (
    <div className="flex items-start gap-2">
      <code className="text-[12px] font-mono text-gray-800 break-all leading-relaxed select-all flex-1">{display}</code>
      <div className="shrink-0">
        <CopyButton text={display} size={13} />
      </div>
    </div>
  );
}

/** Full-width card for a single identifier property */
function IdentifierCard({ comparison, databases, entries }: { comparison: PropertyComparison; databases: string[]; entries: DatabaseEntry[] }) {
  const style = STATUS_STYLES[comparison.status] || STATUS_STYLES.missing;
  const propLabel = PROP_LABELS[comparison.property_name] || comparison.property_name;
  // Databases first, Resolved last
  const sortedDbs = [...databases.filter(d => d !== 'Resolved'), ...databases.filter(d => d === 'Resolved')];

  return (
    <motion.div
      initial={{ opacity: 0, y: 8 }}
      animate={{ opacity: 1, y: 0 }}
      className={`rounded-xl border border-gray-200/80 overflow-hidden bg-white shadow-sm ${style.glow}`}
    >
      {/* Card header */}
      <div className="flex items-center justify-between px-4 py-2.5 bg-gradient-to-r from-gray-50 to-slate-50 border-b border-gray-100">
        <span className="text-[11px] font-bold text-slate-700 uppercase tracking-widest">
          {propLabel}
        </span>
        <span className={`inline-flex items-center gap-1.5 px-2.5 py-0.5 text-[10px] font-bold rounded-full ${style.badge}`}>
          <span className="text-[8px]">{style.icon}</span>
          {style.badgeText}
        </span>
      </div>

      {/* Per-database values */}
      <div className="divide-y divide-gray-100/80">
        {sortedDbs.map((db) => {
          const value = comparison.values[db];
          const c = DB_COLORS[db] || fallbackDB;
          const isResolved = db === 'Resolved';
          const isSmilesProperty = comparison.property_name === 'canonical_smiles';

          return (
            <div key={db} className={`flex items-start gap-3 px-4 py-2.5 group hover:bg-slate-50/50 transition-colors ${isResolved ? 'bg-emerald-50/20' : ''}`}>
              <div className="flex items-center gap-2 pt-0.5 shrink-0 w-[90px]">
                <span className={`w-2 h-2 rounded-full shrink-0 ${c.dot}`} />
                <span className={`text-[10px] font-bold ${c.text} uppercase tracking-wider`}>{db}</span>
              </div>
              <div className="flex-1 min-w-0">
                {value != null ? (
                  isResolved && isSmilesProperty ? (
                    <ResolvedSmilesRow canonical={value} kekulized={entries.find(e => e.database === 'Resolved')?.kekulized_smiles} />
                  ) : (
                    <div className="flex items-start gap-2">
                      <code className="text-[12px] font-mono text-gray-800 break-all leading-relaxed select-all flex-1">{value}</code>
                      <div className="shrink-0 opacity-0 group-hover:opacity-100 transition-opacity">
                        <CopyButton text={value} size={13} />
                      </div>
                    </div>
                  )
                ) : (
                  <span className="text-[11px] text-gray-300 italic">not available</span>
                )}
              </div>
            </div>
          );
        })}
      </div>

      {/* Detail explanation */}
      {comparison.detail && (
        <div className={`flex gap-3 px-4 py-2.5 border-t ${style.detailBg}`}>
          <div className={`w-0.5 shrink-0 rounded-full self-stretch ${style.accent}`} />
          <p className="text-[11px] text-gray-600 leading-relaxed">{comparison.detail}</p>
        </div>
      )}
    </motion.div>
  );
}

/* ── Main component ────────────────────────────────── */

export function DatabaseComparisonPanel({ result }: DatabaseComparisonPanelProps) {
  const config = VERDICT_CONFIG[result.overall_verdict] || VERDICT_CONFIG.no_data;
  const VerdictIcon = config.Icon;
  const foundEntries = result.entries.filter(e => e.found);
  const databases = foundEntries.map(e => e.database);
  // Structure cards: external databases first, Resolved last
  const externalWithSmiles = foundEntries.filter(
    (e): e is DatabaseEntry & { canonical_smiles: string } => !!e.canonical_smiles && e.database !== 'Resolved'
  );
  const resolvedWithSmiles = foundEntries.filter(
    (e): e is DatabaseEntry & { canonical_smiles: string } => !!e.canonical_smiles && e.database === 'Resolved'
  );
  const allStructures = [...externalWithSmiles, ...resolvedWithSmiles];
  const filtered = result.comparisons.filter(c => c.property_name in PROP_LABELS);

  return (
    <motion.div
      initial={{ opacity: 0, y: 16 }}
      animate={{ opacity: 1, y: 0 }}
      transition={{ duration: 0.5, ease: [0.22, 1, 0.36, 1] }}
      className={`rounded-2xl border border-gray-200/60 bg-white overflow-hidden shadow-lg shadow-gray-200/40 ring-1 ${config.ring}`}
    >
      {/* ── Gradient header ── */}
      <div className={`relative bg-gradient-to-r ${config.gradient} px-6 py-4 overflow-hidden`}>
        {/* Decorative mesh overlay */}
        <div className="absolute inset-0 opacity-10" style={{
          backgroundImage: 'radial-gradient(circle at 20% 50%, white 1px, transparent 1px), radial-gradient(circle at 80% 20%, white 1px, transparent 1px)',
          backgroundSize: '30px 30px',
        }} />
        <div className="relative flex items-center justify-between">
          <div className="flex items-center gap-3">
            <div className="w-10 h-10 rounded-xl bg-white/20 backdrop-blur-sm flex items-center justify-center ring-1 ring-white/30">
              <VerdictIcon className="w-5 h-5 text-white" />
            </div>
            <div>
              <h4 className="text-base font-bold text-white tracking-tight">
                Cross-Database Comparison
              </h4>
              <p className="text-[11px] text-white/75 mt-0.5 max-w-md">
                {result.summary}
              </p>
            </div>
          </div>
          <div className="flex items-center gap-2">
            {/* Database count pill */}
            <span className="px-2.5 py-1 text-[10px] font-bold rounded-lg bg-white/15 backdrop-blur-sm text-white/90 border border-white/20">
              {foundEntries.filter(e => e.database !== 'Resolved').length} databases
            </span>
            <span className="px-3 py-1 text-[11px] font-bold rounded-lg bg-white/25 backdrop-blur-sm text-white border border-white/30">
              {config.label}
            </span>
          </div>
        </div>
      </div>

      <div className="p-5 space-y-5">
        {/* ── Database status badges ── */}
        <div className="flex flex-wrap gap-2">
          <AnimatePresence>
            {result.entries.map((entry, i) => {
              const c = DB_COLORS[entry.database] || fallbackDB;
              return (
                <motion.div
                  key={entry.database}
                  initial={{ opacity: 0, scale: 0.9 }}
                  animate={{ opacity: 1, scale: 1 }}
                  transition={{ delay: i * 0.05 }}
                  className={`flex items-center gap-2 px-3 py-1.5 rounded-lg text-[11px] font-semibold border transition-all hover:shadow-md cursor-default ${
                    entry.found
                      ? `${c.bg} ${c.text} ${c.border}`
                      : 'bg-gray-50 text-gray-400 border-gray-200'
                  }`}
                >
                  <span className={`w-2 h-2 rounded-full ${entry.found ? c.dot : 'bg-gray-300'}`} />
                  {entry.database}
                  {entry.found && entry.url && (
                    <a href={safeHref(entry.url)} target="_blank" rel="noopener noreferrer"
                      className="ml-1 opacity-50 hover:opacity-100 transition-opacity">
                      <svg className="w-3 h-3" fill="none" stroke="currentColor" viewBox="0 0 24 24">
                        <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={2} d="M10 6H6a2 2 0 00-2 2v10a2 2 0 002 2h10a2 2 0 002-2v-4M14 4h6m0 0v6m0-6L10 14" />
                      </svg>
                    </a>
                  )}
                  {!entry.found && <span className="text-[9px] font-normal opacity-60 ml-1">N/A</span>}
                </motion.div>
              );
            })}
          </AnimatePresence>
        </div>

        {/* ── Structure comparison — all in one horizontal row ── */}
        {allStructures.length >= 2 && (
          <div>
            <div className="flex items-center gap-2 mb-3">
              <div className="h-px flex-1 bg-gradient-to-r from-gray-200 to-transparent" />
              <span className="text-[10px] font-bold text-gray-400 uppercase tracking-widest px-2">
                Structures
              </span>
              <div className="h-px flex-1 bg-gradient-to-l from-gray-200 to-transparent" />
            </div>
            <div className="flex gap-2.5 overflow-x-auto pb-1">
              {allStructures.map((entry, i) => (
                <motion.div
                  key={entry.database}
                  initial={{ opacity: 0, x: -10 }}
                  animate={{ opacity: 1, x: 0 }}
                  transition={{ delay: i * 0.06 }}
                  className="shrink-0"
                  style={{ width: `${Math.min(200, Math.floor(100 / allStructures.length))}%`, minWidth: 160 }}
                >
                  <StructureCard entry={entry} />
                </motion.div>
              ))}
            </div>
          </div>
        )}

        {/* ── Identifier comparison — one card per property ── */}
        {filtered.length > 0 && databases.length >= 2 && (
          <div>
            <div className="flex items-center gap-2 mb-3">
              <div className="h-px flex-1 bg-gradient-to-r from-gray-200 to-transparent" />
              <span className="text-[10px] font-bold text-gray-400 uppercase tracking-widest px-2">
                Identifiers
              </span>
              <div className="h-px flex-1 bg-gradient-to-l from-gray-200 to-transparent" />
            </div>
            <div className="space-y-3">
              {filtered.map((comp, i) => (
                <motion.div
                  key={comp.property_name}
                  initial={{ opacity: 0, y: 10 }}
                  animate={{ opacity: 1, y: 0 }}
                  transition={{ delay: 0.15 + i * 0.08 }}
                >
                  <IdentifierCard comparison={comp} databases={databases} entries={foundEntries} />
                </motion.div>
              ))}
            </div>
          </div>
        )}
      </div>
    </motion.div>
  );
}
