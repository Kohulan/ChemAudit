/**
 * DeepValidationTab Component
 *
 * Main tab container for the "Deep Validation" section of SingleValidation.
 * Features:
 * - Segmented control: Category view | Severity view
 * - Gear icon to open SeverityConfigPanel
 * - Category view: 3 collapsible domain sections (Stereo & Tautomers, Chemical Composition, Structural Complexity)
 * - Severity view: checks grouped by effective severity level (ERROR → WARNING → INFO → PASS)
 * - Dynamic verdict based on user severity overrides
 * - Summary header showing passed/issue counts
 * - Loading skeleton, empty, and error states
 */
import { useState } from 'react';
import { motion, AnimatePresence } from 'framer-motion';
import { Settings, ChevronDown, Microscope } from 'lucide-react';
import { Badge } from '../ui/Badge';
import { DeepCheckCard } from './DeepCheckCard';
import { SeverityConfigPanel } from './SeverityConfigPanel';
import { useDeepValidationConfig } from '../../hooks/useDeepValidationConfig';
import { DEEP_CHECK_DOMAINS } from '../../types/validation';
import { cn } from '../../lib/utils';
import type { CheckResult } from '../../types/validation';

type ViewMode = 'category' | 'severity';

// All deep check names (union of all domain check arrays)
const ALL_DEEP_CHECK_NAMES = new Set<string>(
  Object.values(DEEP_CHECK_DOMAINS).flatMap((d) => d.checks)
);

interface DeepValidationTabProps {
  checks: CheckResult[];
  onHighlightAtoms: (atoms: number[]) => void;
}

interface DomainSectionProps {
  label: string;
  checks: CheckResult[];
  effectiveSeverityFn: (checkName: string, defaultSeverity: string) => string;
  onHighlightAtoms: (atoms: number[]) => void;
}

/** Computes the worst severity from a list of failing checks */
function worstSeverity(checks: CheckResult[], getSeverity: (c: CheckResult) => string): string {
  const order = ['error', 'critical', 'warning', 'info', 'pass'];
  let worst = 'pass';
  for (const check of checks) {
    if (!check.passed) {
      const sev = getSeverity(check);
      if (order.indexOf(sev) < order.indexOf(worst)) {
        worst = sev;
      }
    }
  }
  return worst;
}

function getDomainVariant(
  severity: string
): 'error' | 'warning' | 'info' | 'success' | 'default' {
  switch (severity) {
    case 'error':
    case 'critical':
      return 'error';
    case 'warning':
      return 'warning';
    case 'info':
      return 'info';
    default:
      return 'success';
  }
}

/** Collapsible domain section */
function DomainSection({ label, checks, effectiveSeverityFn, onHighlightAtoms }: DomainSectionProps) {
  const [open, setOpen] = useState(true);
  const failCount = checks.filter((c) => !c.passed).length;
  const worst = worstSeverity(checks, (c) => effectiveSeverityFn(c.check_name, c.severity));

  return (
    <div className="rounded-xl border border-[var(--color-border)] overflow-hidden">
      {/* Section header */}
      <button
        onClick={() => setOpen((prev) => !prev)}
        className={cn(
          'w-full flex items-center gap-3 px-4 py-3',
          'bg-[var(--color-surface-sunken)] hover:bg-[var(--color-surface-elevated)]',
          'transition-colors text-left'
        )}
      >
        <ChevronDown
          className={cn(
            'w-4 h-4 text-[var(--color-text-muted)] transition-transform duration-200 flex-shrink-0',
            !open && '-rotate-90'
          )}
        />
        <span className="flex-1 text-sm font-semibold text-[var(--color-text-primary)]">
          {label}
        </span>
        <div className="flex items-center gap-2">
          <span className="text-xs text-[var(--color-text-muted)]">
            {checks.length} checks
          </span>
          {failCount > 0 ? (
            <Badge variant={getDomainVariant(worst)} size="sm">
              {failCount} issue{failCount !== 1 ? 's' : ''}
            </Badge>
          ) : (
            <Badge variant="success" size="sm">All pass</Badge>
          )}
        </div>
      </button>

      {/* Section content */}
      <AnimatePresence initial={false}>
        {open && (
          <motion.div
            initial={{ opacity: 0, height: 0 }}
            animate={{ opacity: 1, height: 'auto' }}
            exit={{ opacity: 0, height: 0 }}
            transition={{ duration: 0.2 }}
            className="overflow-hidden"
          >
            <div className="p-3 space-y-2">
              {checks.map((check) => (
                <DeepCheckCard
                  key={check.check_name}
                  check={check}
                  effectiveSeverity={effectiveSeverityFn(check.check_name, check.severity)}
                  onHighlightAtoms={onHighlightAtoms}
                />
              ))}
            </div>
          </motion.div>
        )}
      </AnimatePresence>
    </div>
  );
}

/**
 * Renders the Deep Validation tab with category/severity views and config panel.
 */
export function DeepValidationTab({ checks, onHighlightAtoms }: DeepValidationTabProps) {
  const [viewMode, setViewMode] = useState<ViewMode>('category');
  const [configOpen, setConfigOpen] = useState(false);

  const { config, setSeverityOverride, removeSeverityOverride, resetAllOverrides, getEffectiveSeverity } =
    useDeepValidationConfig();

  // Filter to only deep validation checks
  const deepChecks = checks.filter((c) => ALL_DEEP_CHECK_NAMES.has(c.check_name));

  // Empty state
  if (deepChecks.length === 0) {
    return (
      <div className="flex flex-col items-center justify-center py-12 text-center">
        <div className="w-12 h-12 rounded-xl bg-[var(--color-primary)]/10 flex items-center justify-center text-[var(--color-primary)] mb-4">
          <Microscope className="w-6 h-6" />
        </div>
        <h3 className="text-base font-semibold text-[var(--color-text-primary)] mb-2">
          No Deep Validation Results
        </h3>
        <p className="text-sm text-[var(--color-text-secondary)] max-w-sm">
          Deep validation results are not available. Run a validation to see advanced
          structure checks including stereo, tautomers, composition, and complexity analysis.
        </p>
      </div>
    );
  }

  // Compute summary counts
  const totalChecks = deepChecks.length;
  const passedChecks = deepChecks.filter((c) => c.passed).length;

  // Compute dynamic verdict based on effective severities
  const failingChecks = deepChecks.filter((c) => !c.passed);
  const hasEffectiveError = failingChecks.some(
    (c) => getEffectiveSeverity(c.check_name, c.severity) === 'error' ||
           getEffectiveSeverity(c.check_name, c.severity) === 'critical'
  );
  const hasEffectiveWarning = failingChecks.some(
    (c) => getEffectiveSeverity(c.check_name, c.severity) === 'warning'
  );

  const verdict = hasEffectiveError ? 'FAIL' : hasEffectiveWarning ? 'WARNING' : 'PASS';
  const verdictVariant =
    verdict === 'FAIL' ? 'error' : verdict === 'WARNING' ? 'warning' : 'success';

  // Build check map by domain for category view
  const checksByDomain: Record<string, CheckResult[]> = {};
  for (const [domainKey, domain] of Object.entries(DEEP_CHECK_DOMAINS)) {
    checksByDomain[domainKey] = domain.checks
      .map((name) => deepChecks.find((c) => c.check_name === name))
      .filter((c): c is CheckResult => Boolean(c));
  }

  // Build check groups by severity for severity view
  const checksBySeverity: Record<string, CheckResult[]> = {
    error: [],
    warning: [],
    info: [],
    pass: [],
  };
  for (const check of deepChecks) {
    const sev = getEffectiveSeverity(check.check_name, check.severity);
    if (check.passed) {
      checksBySeverity.pass.push(check);
    } else if (sev === 'error' || sev === 'critical') {
      checksBySeverity.error.push(check);
    } else if (sev === 'warning') {
      checksBySeverity.warning.push(check);
    } else {
      checksBySeverity.info.push(check);
    }
  }

  const severityGroups = ([
    { key: 'error', label: 'Error', variant: 'error' as const, checks: checksBySeverity.error },
    { key: 'warning', label: 'Warning', variant: 'warning' as const, checks: checksBySeverity.warning },
    { key: 'info', label: 'Info', variant: 'info' as const, checks: checksBySeverity.info },
    { key: 'pass', label: 'Pass', variant: 'success' as const, checks: checksBySeverity.pass },
  ]).filter((g) => g.checks.length > 0);

  return (
    <div className="space-y-4">
      {/* Summary header + verdict */}
      <div className="flex items-center justify-between flex-wrap gap-3">
        <div className="flex items-center gap-3 flex-wrap">
          <div className="flex items-center gap-2">
            <span className="text-sm font-medium text-[var(--color-text-secondary)]">
              {passedChecks}/{totalChecks} passed
            </span>
          </div>

          {/* Verdict badge */}
          <Badge variant={verdictVariant} size="md">
            Deep Validation: {verdict}
          </Badge>

          {/* Issue count badges */}
          {checksBySeverity.error.length > 0 && (
            <Badge variant="error" size="sm">
              {checksBySeverity.error.length} error{checksBySeverity.error.length !== 1 ? 's' : ''}
            </Badge>
          )}
          {checksBySeverity.warning.length > 0 && (
            <Badge variant="warning" size="sm">
              {checksBySeverity.warning.length} warning{checksBySeverity.warning.length !== 1 ? 's' : ''}
            </Badge>
          )}
        </div>

        {/* View toggle + gear icon */}
        <div className="flex items-center gap-2">
          {/* Segmented control */}
          <div className="flex items-center rounded-lg border border-[var(--color-border)] bg-[var(--color-surface-sunken)] p-0.5">
            <button
              onClick={() => setViewMode('category')}
              className={cn(
                'px-3 py-1.5 text-xs font-medium rounded-md transition-all',
                viewMode === 'category'
                  ? 'bg-[var(--color-surface-elevated)] text-[var(--color-text-primary)] shadow-sm'
                  : 'text-[var(--color-text-muted)] hover:text-[var(--color-text-secondary)]'
              )}
            >
              Category
            </button>
            <button
              onClick={() => setViewMode('severity')}
              className={cn(
                'px-3 py-1.5 text-xs font-medium rounded-md transition-all',
                viewMode === 'severity'
                  ? 'bg-[var(--color-surface-elevated)] text-[var(--color-text-primary)] shadow-sm'
                  : 'text-[var(--color-text-muted)] hover:text-[var(--color-text-secondary)]'
              )}
            >
              Severity
            </button>
          </div>

          {/* Gear icon — opens SeverityConfigPanel */}
          <button
            onClick={() => setConfigOpen(true)}
            className={cn(
              'p-2 rounded-lg border border-[var(--color-border)] transition-all',
              'text-[var(--color-text-muted)] hover:text-[var(--color-text-primary)]',
              'hover:bg-[var(--color-surface-sunken)]',
              Object.keys(config.severityOverrides).length > 0 &&
                'border-[var(--color-primary)]/30 text-[var(--color-primary)]'
            )}
            title="Configure severity levels"
            aria-label="Open severity configuration"
          >
            <Settings className="w-4 h-4" />
          </button>
        </div>
      </div>

      {/* Content area */}
      <AnimatePresence mode="wait">
        <motion.div
          key={viewMode}
          initial={{ opacity: 0, y: 8 }}
          animate={{ opacity: 1, y: 0 }}
          exit={{ opacity: 0, y: -8 }}
          transition={{ duration: 0.15 }}
        >
          {viewMode === 'category' ? (
            /* Category view: 3 domain sections */
            <div className="space-y-3">
              {Object.entries(DEEP_CHECK_DOMAINS).map(([domainKey, domain]) => {
                const domainChecks = checksByDomain[domainKey] ?? [];
                if (domainChecks.length === 0) return null;
                return (
                  <DomainSection
                    key={domainKey}
                    label={domain.label}
                    checks={domainChecks}
                    effectiveSeverityFn={getEffectiveSeverity}
                    onHighlightAtoms={onHighlightAtoms}
                  />
                );
              })}
            </div>
          ) : (
            /* Severity view: grouped by effective severity */
            <div className="space-y-3">
              {severityGroups.map((group) => (
                <div key={group.key} className="rounded-xl border border-[var(--color-border)] overflow-hidden">
                  <div className="flex items-center gap-3 px-4 py-3 bg-[var(--color-surface-sunken)]">
                    <Badge variant={group.variant} size="sm">
                      {group.label.toUpperCase()}
                    </Badge>
                    <span className="text-xs text-[var(--color-text-muted)]">
                      {group.checks.length} check{group.checks.length !== 1 ? 's' : ''}
                    </span>
                  </div>
                  <div className="p-3 space-y-2">
                    {group.checks.map((check) => (
                      <DeepCheckCard
                        key={check.check_name}
                        check={check}
                        effectiveSeverity={getEffectiveSeverity(check.check_name, check.severity)}
                        onHighlightAtoms={onHighlightAtoms}
                      />
                    ))}
                  </div>
                </div>
              ))}
            </div>
          )}
        </motion.div>
      </AnimatePresence>

      {/* Severity Config Panel */}
      <SeverityConfigPanel
        isOpen={configOpen}
        onClose={() => setConfigOpen(false)}
        checks={deepChecks}
        config={config}
        onSetSeverity={setSeverityOverride}
        onRemoveOverride={removeSeverityOverride}
        onResetAll={resetAllOverrides}
      />
    </div>
  );
}
