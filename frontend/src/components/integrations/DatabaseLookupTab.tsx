import { useId, type ComponentProps } from 'react';
import { AnimatePresence, motion } from 'framer-motion';
import { ChevronDown, Database, GitCompareArrows, Info, Search } from 'lucide-react';
import { ClayButton } from '../ui/ClayButton';
import { DatabaseLookupResults } from './DatabaseLookupResults';

interface DatabaseLookupControlsProps {
  /** Run the four-database lookup for the current molecule. */
  onLookup: () => void;
  /** Disable the Look Up button (no molecule, or another request in flight). */
  lookupDisabled: boolean;
  /** A database lookup is currently running. */
  databaseLoading: boolean;
  /** A cross-database comparison is currently running. */
  isComparing: boolean;
  /** Whether comparison runs automatically after a lookup. */
  autoCompare: boolean;
  /** Toggle the auto-compare preference. */
  onToggleAutoCompare: () => void;
  /** Show the manual Compare button (auto-compare off, results present, not yet compared). */
  showManualCompare: boolean;
  /** Run a manual cross-database comparison. */
  onCompare: () => void;
}

/**
 * Database Lookup tab — left-column input controls. Lists the queried
 * databases, the Look Up button, the auto-compare toggle, and the manual
 * Compare button. Purely presentational; all state and handlers live in the
 * parent page (they are coupled to result caching, snapshots, and export).
 */
export function DatabaseLookupControls({
  onLookup,
  lookupDisabled,
  databaseLoading,
  isComparing,
  autoCompare,
  onToggleAutoCompare,
  showManualCompare,
  onCompare,
}: DatabaseLookupControlsProps) {
  const autoCompareLabelId = useId();
  return (
    <div className="space-y-4">
      <div className="flex items-start gap-3">
        <div className="w-8 h-8 rounded-lg bg-[var(--color-primary)]/10 flex items-center justify-center text-[var(--color-primary)] flex-shrink-0">
          <Info className="w-4 h-4" />
        </div>
        <div>
          <p className="text-[var(--color-text-secondary)] text-sm mb-3">
            Click Look Up to query four databases for this molecule:
          </p>
          <ul className="list-none space-y-1 text-sm text-[var(--color-text-secondary)]">
            <li className="flex items-center gap-2">
              <span className="w-1.5 h-1.5 rounded-full bg-blue-500"></span>
              <strong className="text-[var(--color-text-primary)]">PubChem</strong> — NIH compound database (100M+ compounds)
            </li>
            <li className="flex items-center gap-2">
              <span className="w-1.5 h-1.5 rounded-full bg-purple-500"></span>
              <strong className="text-[var(--color-text-primary)]">ChEMBL</strong> — Bioactivity and drug data
            </li>
            <li className="flex items-center gap-2">
              <span className="w-1.5 h-1.5 rounded-full bg-amber-500"></span>
              <strong className="text-[var(--color-text-primary)]">COCONUT</strong> — Natural products database
            </li>
            <li className="flex items-center gap-2">
              <span className="w-1.5 h-1.5 rounded-full bg-teal-500"></span>
              <strong className="text-[var(--color-text-primary)]">Wikidata</strong> — Open knowledge base
            </li>
            <li className="flex items-center gap-2">
              <span className="w-1.5 h-1.5 rounded-full bg-rose-500"></span>
              <strong className="text-[var(--color-text-primary)]">SureChEMBL</strong> — Patent literature presence
            </li>
          </ul>
        </div>
      </div>
      <div className="flex items-center gap-3 flex-wrap">
        <ClayButton
          variant="primary"
          onClick={onLookup}
          disabled={lookupDisabled}
          loading={databaseLoading || isComparing}
          leftIcon={<Search className="w-4 h-4" />}
        >
          {isComparing ? 'Comparing...' : databaseLoading ? 'Looking up...' : 'Look Up'}
        </ClayButton>
        {/* Auto-compare toggle */}
        <label className="flex items-center gap-2 cursor-pointer select-none">
          <button
            type="button"
            role="switch"
            aria-checked={autoCompare}
            aria-labelledby={autoCompareLabelId}
            onClick={onToggleAutoCompare}
            className={`relative w-8 h-[18px] rounded-full transition-colors focus-visible:outline-none focus-visible:ring-2 focus-visible:ring-[var(--color-primary)] focus-visible:ring-offset-2 ${autoCompare ? 'bg-[var(--color-primary)]' : 'bg-chem-dark-300 dark:bg-chem-dark-600'}`}
          >
            <span className={`absolute top-[2px] w-[14px] h-[14px] rounded-full bg-chem-dark-50 shadow-sm transition-transform ${autoCompare ? 'translate-x-[16px]' : 'translate-x-[2px]'}`} />
          </button>
          <span id={autoCompareLabelId} className="text-[11px] text-[var(--color-text-muted)] font-medium">Auto-compare</span>
        </label>
        {/* Manual compare button — shown when toggle is off, lookup is done, no comparison yet */}
        {showManualCompare && (
          <ClayButton
            variant="outline"
            size="sm"
            onClick={onCompare}
            disabled={isComparing}
            loading={isComparing}
            leftIcon={<GitCompareArrows className="w-3.5 h-3.5" />}
          >
            Compare
          </ClayButton>
        )}
      </div>
    </div>
  );
}

type DatabaseResults = ComponentProps<typeof DatabaseLookupResults>['results'];

interface DatabaseResultsBoxProps {
  /** Per-database results to render inside the collapsible body. */
  results: DatabaseResults;
  /** Whether the box is expanded. */
  expanded: boolean;
  /** Toggle the expanded state. */
  onToggleExpanded: () => void;
}

/**
 * Database Lookup tab — collapsible "Database Results" box shown below the
 * grid once a lookup completes. Header toggles the per-database detail cards.
 * Presentational; expand state is owned by the parent so it persists across
 * tab switches.
 */
export function DatabaseResultsBox({ results, expanded, onToggleExpanded }: DatabaseResultsBoxProps) {
  return (
    <motion.div
      initial={{ opacity: 0, y: 20 }}
      animate={{ opacity: 1, y: 0 }}
      exit={{ opacity: 0, y: -20 }}
      className="card overflow-hidden"
    >
      <button
        onClick={onToggleExpanded}
        className="w-full flex items-center gap-3 p-4 sm:p-5 hover:bg-[var(--color-surface-hover)] transition-colors cursor-pointer"
      >
        <div className="w-8 h-8 rounded-lg bg-gradient-to-br from-[var(--color-primary)]/10 to-[var(--color-accent)]/10 flex items-center justify-center text-[var(--color-primary)]">
          <Database className="w-4 h-4" />
        </div>
        <div className="flex-1 text-left">
          <h4 className="font-semibold text-[var(--color-text-primary)] text-sm tracking-tight">
            Database Results
          </h4>
          <p className="text-[11px] text-[var(--color-text-muted)]">Individual PubChem, ChEMBL, COCONUT, Wikidata, SureChEMBL details</p>
        </div>
        <ChevronDown className={`w-4 h-4 text-[var(--color-text-muted)] transition-transform duration-200 ${expanded ? 'rotate-180' : ''}`} />
      </button>
      <AnimatePresence>
        {expanded && (
          <motion.div
            initial={{ height: 0, opacity: 0 }}
            animate={{ height: 'auto', opacity: 1 }}
            exit={{ height: 0, opacity: 0 }}
            transition={{ duration: 0.25 }}
            className="overflow-hidden"
          >
            <div className="px-4 pb-4 sm:px-5 sm:pb-5 border-t border-[var(--color-border)]">
              <div className="pt-4">
                <DatabaseLookupResults results={results} />
              </div>
            </div>
          </motion.div>
        )}
      </AnimatePresence>
    </motion.div>
  );
}
