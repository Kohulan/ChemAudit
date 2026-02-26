import { useState } from 'react';
import { motion } from 'framer-motion';
import {
  Shield,
  ShieldCheck,
  ShieldAlert,
  AlertTriangle,
  ChevronDown,
  Info
} from 'lucide-react';
import type { SafetyFilterResult, FilterAlertResult, ChEMBLAlertsResult } from '../../types/scoring';
import { InfoTooltip } from '../ui/Tooltip';
import { cn } from '../../lib/utils';

const CATEGORY_BADGE: Record<string, string> = {
  'Reactive Group': 'bg-red-100 text-red-700',
  'Toxicophore': 'bg-rose-100 text-rose-700',
  'Metabolic Liability': 'bg-amber-100 text-amber-700',
  'Assay Interference': 'bg-purple-100 text-purple-700',
  'Physicochemical': 'bg-slate-100 text-slate-700',
  'Unwanted Functionality': 'bg-gray-100 text-gray-600',
};

function formatPatternName(name: string): string {
  return name
    .replace(/\(\d+\)$/, '')
    .trim()
    .split('_')
    .map((w) => w.charAt(0).toUpperCase() + w.slice(1))
    .join(' ');
}

interface SafetyFiltersScoreProps {
  result: SafetyFilterResult;
}

/**
 * Catalog metadata for display — sourced from AVAILABLE_CATALOGS (backend).
 * Keeps frontend descriptions consistent with the enriched alert data.
 */
const FILTER_INFO: Record<string, { label: string; description: string; citation: string }> = {
  PAINS: {
    label: 'PAINS (Pan-Assay Interference Compounds)',
    description: 'Identifies compounds that appear active in multiple assay types due to non-specific mechanisms such as aggregation, redox cycling, or fluorescence interference.',
    citation: 'Baell JB, Holloway GA. J Med Chem 53 (2010) 2719-2740.',
  },
  Brenk: {
    label: 'Brenk Structural Alerts',
    description: 'Flags compounds with known problematic functional groups associated with toxicity or unfavourable pharmacokinetic properties.',
    citation: 'Brenk R et al. ChemMedChem 3 (2008) 435-444.',
  },
  NIH: {
    label: 'NIH MLPCN Exclusion Filters',
    description: 'Reactive and interference compounds excluded from NIH Molecular Libraries screening collection.',
    citation: 'Jadhav A et al. J Med Chem 53 (2010) 37-51.',
  },
  ZINC: {
    label: 'ZINC Druglike Filters',
    description: 'Reactivity and drug-likeness filters used by the ZINC purchasable compound database.',
    citation: 'Irwin JJ, Shoichet BK. J Chem Inf Model 45 (2005) 177-182.',
  },
};

const CHEMBL_FILTER_INFO: Record<string, { label: string; description: string; citation: string }> = {
  bms: {
    label: 'BMS HTS Desirability Filters',
    description: 'Functional group filters for HTS compound selection at Bristol-Myers Squibb.',
    citation: 'Pearce BC et al. J Chem Inf Model 46 (2006) 1060-1068.',
  },
  dundee: {
    label: 'Dundee NTD Screening Filters',
    description: 'Alerts developed for neglected tropical disease screening at the University of Dundee.',
    citation: 'Brenk R et al. ChemMedChem 3 (2008) 435-444.',
  },
  glaxo: {
    label: 'Glaxo Hard Filters',
    description: 'Reactive and toxic group filters for lead optimisation at GlaxoSmithKline.',
    citation: 'Hann M et al. J Chem Inf Comput Sci 39 (1999) 897-902.',
  },
  inpharmatica: {
    label: 'Inpharmatica Unwanted Fragments',
    description: 'Fragments flagged as undesirable in drug-like compounds by Inpharmatica.',
    citation: 'ChEMBL structural alerts collection.',
  },
  lint: {
    label: 'Lilly MedChem Rules (LINT)',
    description: 'Medicinal chemistry structural alerts developed at Eli Lilly to flag undesirable compounds.',
    citation: 'Bruns RF, Watson IA. J Med Chem 55 (2012) 9763-9772.',
  },
  mlsmr: {
    label: 'NIH MLSMR Excluded Functionality',
    description: 'Functional groups excluded from the NIH Molecular Libraries Small Molecule Repository.',
    citation: 'NIH Chemical Genomics Center guidelines.',
  },
  schembl: {
    label: 'SureChEMBL Non-chemical Patterns',
    description: 'Filters for non-chemical entities and patent noise in the SureChEMBL patent-mined database.',
    citation: 'Papadatos G et al. Nucleic Acids Res 44 (2016) D1220-D1228.',
  },
};

/**
 * Individual filter result display
 */
function FilterCard({
  name,
  description,
  citation,
  result,
  delay = 0
}: {
  name: string;
  description: string;
  citation?: string;
  result: FilterAlertResult;
  delay?: number;
}) {
  const [showAlerts, setShowAlerts] = useState(false);

  return (
    <motion.div
      initial={{ opacity: 0, y: 10 }}
      animate={{ opacity: 1, y: 0 }}
      transition={{ delay }}
      className={cn(
        'rounded-xl p-3',
        'bg-[var(--color-surface-sunken)]',
        'border',
        result.passed
          ? 'border-emerald-500/20'
          : 'border-red-500/20'
      )}
    >
      <div className="flex items-center justify-between mb-2">
        <div className="flex items-center gap-2">
          {result.passed ? (
            <ShieldCheck className="w-4 h-4 text-emerald-500" />
          ) : (
            <ShieldAlert className="w-4 h-4 text-red-500" />
          )}
          <span className="text-sm font-medium text-[var(--color-text-primary)]">{name}</span>
          <InfoTooltip title={name} content={
            <div className="text-xs space-y-1">
              <p>{description}</p>
              {citation && <p className="text-[var(--color-text-muted)] italic">{citation}</p>}
            </div>
          } />
        </div>
        <span
          className={cn(
            'text-xs font-medium px-2 py-0.5 rounded-full',
            result.passed
              ? 'bg-emerald-500/10 text-emerald-600 dark:text-emerald-400'
              : 'bg-red-500/10 text-red-600 dark:text-red-400'
          )}
        >
          {result.passed ? 'Clear' : `${result.alert_count} alert${result.alert_count !== 1 ? 's' : ''}`}
        </span>
      </div>

      {/* Show alerts if any */}
      {!result.passed && result.alerts.length > 0 && (
        <>
          <button
            onClick={() => setShowAlerts(!showAlerts)}
            className={cn(
              'w-full flex items-center justify-between px-2 py-1 mt-1 rounded-lg',
              'text-xs text-red-600 dark:text-red-400',
              'hover:bg-red-500/5 transition-colors'
            )}
          >
            <span>View alerts</span>
            <ChevronDown className={cn(
              'w-3 h-3 transition-transform',
              showAlerts && 'rotate-180'
            )} />
          </button>
          {showAlerts && (
            <motion.div
              initial={{ opacity: 0, height: 0 }}
              animate={{ opacity: 1, height: 'auto' }}
              className="mt-2 p-2 rounded-lg bg-red-500/5 text-xs max-h-32 overflow-y-auto space-y-1"
            >
              {(result.alert_details && result.alert_details.length > 0
                ? result.alert_details
                : result.alerts.map((a) => ({ name: a, category: 'Unwanted Functionality' }))
              ).map((detail, i) => (
                <div key={i} className="flex items-center gap-1.5 py-0.5">
                  <AlertTriangle className="w-3 h-3 flex-shrink-0 text-red-500" />
                  <span className="text-[var(--color-text-secondary)]">{formatPatternName(detail.name)}</span>
                  <span className={cn('px-1.5 py-0 rounded-full text-[10px] font-medium', CATEGORY_BADGE[detail.category] || 'bg-gray-100 text-gray-600')}>
                    {detail.category}
                  </span>
                </div>
              ))}
            </motion.div>
          )}
        </>
      )}
    </motion.div>
  );
}

/**
 * ChEMBL alerts section
 */
function ChEMBLAlertsSection({ chembl }: { chembl: ChEMBLAlertsResult }) {
  const [expanded, setExpanded] = useState(false);

  const alertFilterKeys = ['bms', 'dundee', 'glaxo', 'inpharmatica', 'lint', 'mlsmr', 'schembl'] as const;

  return (
    <motion.div
      initial={{ opacity: 0, y: 10 }}
      animate={{ opacity: 1, y: 0 }}
      transition={{ delay: 0.3 }}
      className={cn(
        'rounded-xl p-3',
        'bg-[var(--color-surface-sunken)]',
        'border',
        chembl.passed ? 'border-emerald-500/20' : 'border-amber-500/20'
      )}
    >
      <button
        onClick={() => setExpanded(!expanded)}
        className="w-full flex items-center justify-between"
      >
        <div className="flex items-center gap-2">
          {chembl.passed ? (
            <ShieldCheck className="w-4 h-4 text-emerald-500" />
          ) : (
            <ShieldAlert className="w-4 h-4 text-amber-500" />
          )}
          <span className="text-sm font-medium text-[var(--color-text-primary)]">ChEMBL Alerts</span>
          <InfoTooltip
            asSpan
            title="ChEMBL Structural Alerts"
            content={
              <span className="text-xs">
                Curated structural alert filter sets from pharmaceutical companies and screening centres, aggregated in ChEMBL.
              </span>
            }
          />
        </div>
        <div className="flex items-center gap-2">
          <span
            className={cn(
              'text-xs font-medium px-2 py-0.5 rounded-full',
              chembl.passed
                ? 'bg-emerald-500/10 text-emerald-600 dark:text-emerald-400'
                : 'bg-amber-500/10 text-amber-600 dark:text-amber-400'
            )}
          >
            {chembl.passed ? 'Clear' : `${chembl.total_alerts} alert${chembl.total_alerts !== 1 ? 's' : ''}`}
          </span>
          <ChevronDown className={cn(
            'w-4 h-4 text-[var(--color-text-muted)] transition-transform',
            expanded && 'rotate-180'
          )} />
        </div>
      </button>

      {expanded && (
        <motion.div
          initial={{ opacity: 0, height: 0 }}
          animate={{ opacity: 1, height: 'auto' }}
          className="mt-3 space-y-1"
        >
          {alertFilterKeys.map(key => {
            const filterResult = chembl[key as keyof ChEMBLAlertsResult] as FilterAlertResult | null;
            const info = CHEMBL_FILTER_INFO[key];
            if (!filterResult || !info) return null;

            return (
              <div key={key} className="py-1.5">
                <div className="flex items-center justify-between text-xs">
                  <div className="flex items-center gap-2">
                    {filterResult.passed ? (
                      <ShieldCheck className="w-3 h-3 text-emerald-500" />
                    ) : (
                      <ShieldAlert className="w-3 h-3 text-red-500" />
                    )}
                    <span className="text-[var(--color-text-secondary)]">{info.label}</span>
                    <InfoTooltip
                      asSpan
                      title={info.label}
                      content={
                        <div className="text-xs space-y-1">
                          <p>{info.description}</p>
                          <p className="text-[var(--color-text-muted)] italic">{info.citation}</p>
                        </div>
                      }
                    />
                  </div>
                  <span className={cn(
                    'text-xs',
                    filterResult.passed ? 'text-emerald-500' : 'text-red-500'
                  )}>
                    {filterResult.passed ? 'Pass' : `${filterResult.alert_count} alert${filterResult.alert_count !== 1 ? 's' : ''}`}
                  </span>
                </div>
                {/* Show enriched alert details when filter has alerts */}
                {!filterResult.passed && filterResult.alert_details && filterResult.alert_details.length > 0 && (
                  <div className="ml-5 mt-1 space-y-0.5">
                    {filterResult.alert_details.map((detail, i) => (
                      <div key={i} className="flex items-center gap-1.5 text-[11px]">
                        <AlertTriangle className="w-2.5 h-2.5 flex-shrink-0 text-red-400" />
                        <span className="text-[var(--color-text-muted)]">{formatPatternName(detail.name)}</span>
                        <span className={cn('px-1.5 rounded-full text-[10px] font-medium', CATEGORY_BADGE[detail.category] || 'bg-gray-100 text-gray-600')}>
                          {detail.category}
                        </span>
                      </div>
                    ))}
                  </div>
                )}
              </div>
            );
          })}
        </motion.div>
      )}
    </motion.div>
  );
}

/**
 * Displays safety filter results including PAINS, Brenk, NIH, ZINC, and ChEMBL.
 */
export function SafetyFiltersScore({ result }: SafetyFiltersScoreProps) {
  const { pains, brenk, nih, zinc, chembl, all_passed, total_alerts, interpretation } = result;

  return (
    <motion.div
      initial={{ opacity: 0, y: 20 }}
      animate={{ opacity: 1, y: 0 }}
      className={cn(
        'rounded-2xl p-5',
        'bg-gradient-to-br from-[var(--color-surface-elevated)] to-[var(--color-surface)]',
        'border',
        all_passed
          ? 'border-emerald-500/20'
          : 'border-red-500/20'
      )}
    >
      {/* Header */}
      <div className="flex items-center justify-between mb-4">
        <div className="flex items-center gap-3">
          <div className={cn(
            'w-10 h-10 rounded-xl flex items-center justify-center',
            all_passed
              ? 'bg-emerald-500/10 text-emerald-500'
              : 'bg-red-500/10 text-red-500'
          )}>
            <Shield className="w-5 h-5" />
          </div>
          <div>
            <h3 className="font-semibold text-[var(--color-text-primary)]">Safety Filters</h3>
            <p className="text-xs text-[var(--color-text-muted)]">Structural alert screening</p>
          </div>
        </div>

        {/* Overall status badge */}
        <div className={cn(
          'flex items-center gap-2 px-3 py-1.5 rounded-full',
          all_passed
            ? 'bg-emerald-500/10 text-emerald-600 dark:text-emerald-400'
            : 'bg-red-500/10 text-red-600 dark:text-red-400'
        )}>
          {all_passed ? (
            <>
              <ShieldCheck className="w-4 h-4" />
              <span className="text-sm font-medium">All Clear</span>
            </>
          ) : (
            <>
              <AlertTriangle className="w-4 h-4" />
              <span className="text-sm font-medium">{total_alerts} Alert{total_alerts !== 1 ? 's' : ''}</span>
            </>
          )}
        </div>
      </div>

      {/* Filter cards */}
      <div className="space-y-2">
        {([
          { key: 'PAINS', result: pains, delay: 0.1 },
          { key: 'Brenk', result: brenk, delay: 0.15 },
          { key: 'NIH', result: nih, delay: 0.2 },
          { key: 'ZINC', result: zinc, delay: 0.25 },
        ] as const).map(({ key, result: filterResult, delay }) => {
          if (!filterResult) return null;
          const info = FILTER_INFO[key];
          return (
            <FilterCard
              key={key}
              name={info.label}
              description={info.description}
              citation={info.citation}
              result={filterResult}
              delay={delay}
            />
          );
        })}
        {chembl && <ChEMBLAlertsSection chembl={chembl} />}
      </div>

      {/* Interpretation */}
      <div className={cn(
        'mt-4 p-3 rounded-xl text-xs',
        all_passed
          ? 'bg-emerald-500/5 border border-emerald-500/10 text-emerald-700 dark:text-emerald-300'
          : 'bg-red-500/5 border border-red-500/10 text-red-700 dark:text-red-300'
      )}>
        <div className="flex items-start gap-2">
          <Info className="w-4 h-4 flex-shrink-0 mt-0.5" />
          <span>{interpretation}</span>
        </div>
      </div>
    </motion.div>
  );
}
