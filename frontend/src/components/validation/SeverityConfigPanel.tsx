/**
 * SeverityConfigPanel Component
 *
 * A modal panel for configuring per-check severity level overrides.
 * Lists all 16 deep validation checks grouped by domain, allowing users
 * to override each check's severity (ERROR | WARNING | INFO).
 * Overridden checks are highlighted with a visual indicator dot.
 * Enters/exits with Framer Motion slide animation.
 */
import { motion, AnimatePresence } from 'framer-motion';
import { X, RotateCcw } from 'lucide-react';
import { Badge } from '../ui/Badge';
import { cn } from '../../lib/utils';
import {
  DEEP_CHECK_DOMAINS,
} from '../../types/validation';
import type { CheckResult, DeepValidationConfig, SeverityOverride } from '../../types/validation';

interface SeverityConfigPanelProps {
  isOpen: boolean;
  onClose: () => void;
  checks: CheckResult[];
  config: DeepValidationConfig;
  onSetSeverity: (checkName: string, severity: SeverityOverride) => void;
  onRemoveOverride: (checkName: string) => void;
  onResetAll: () => void;
}

const SEVERITY_OPTIONS: SeverityOverride[] = ['error', 'warning', 'info'];

function getSeverityVariant(
  severity: string
): 'error' | 'warning' | 'info' | 'success' | 'default' {
  switch (severity) {
    case 'error':
      return 'error';
    case 'warning':
      return 'warning';
    case 'info':
      return 'info';
    default:
      return 'default';
  }
}

/** Formats snake_case check name to Title Case */
function formatCheckName(name: string): string {
  return name
    .split('_')
    .map((word) => word.charAt(0).toUpperCase() + word.slice(1))
    .join(' ');
}

/**
 * Modal panel for configuring per-check severity overrides.
 */
export function SeverityConfigPanel({
  isOpen,
  onClose,
  checks,
  config,
  onSetSeverity,
  onRemoveOverride,
  onResetAll,
}: SeverityConfigPanelProps) {
  // Build a lookup for current checks by name
  const checkLookup: Record<string, CheckResult> = {};
  for (const check of checks) {
    checkLookup[check.check_name] = check;
  }

  const hasOverrides = Object.keys(config.severityOverrides).length > 0;

  return (
    <AnimatePresence>
      {isOpen && (
        <>
          {/* Backdrop */}
          <motion.div
            key="backdrop"
            initial={{ opacity: 0 }}
            animate={{ opacity: 1 }}
            exit={{ opacity: 0 }}
            transition={{ duration: 0.15 }}
            className="fixed inset-0 z-40 bg-black/40 backdrop-blur-sm"
            onClick={onClose}
          />

          {/* Panel - slides in from right */}
          <motion.div
            key="panel"
            initial={{ opacity: 0, x: '100%' }}
            animate={{ opacity: 1, x: 0 }}
            exit={{ opacity: 0, x: '100%' }}
            transition={{ type: 'spring', damping: 25, stiffness: 300 }}
            className={cn(
              'fixed top-0 right-0 bottom-0 z-50',
              'w-full max-w-md',
              'bg-[var(--color-surface-elevated)] border-l border-[var(--color-border)]',
              'flex flex-col shadow-2xl'
            )}
          >
            {/* Header */}
            <div className="flex items-center justify-between px-5 py-4 border-b border-[var(--color-border)]">
              <div>
                <h3 className="font-semibold text-[var(--color-text-primary)] text-base">
                  Severity Configuration
                </h3>
                <p className="text-xs text-[var(--color-text-muted)] mt-0.5">
                  Override per-check severity. Changes affect the verdict dynamically.
                </p>
              </div>
              <button
                onClick={onClose}
                className={cn(
                  'p-2 rounded-lg transition-colors flex-shrink-0',
                  'text-[var(--color-text-muted)] hover:text-[var(--color-text-primary)]',
                  'hover:bg-[var(--color-surface-sunken)]'
                )}
                aria-label="Close settings"
              >
                <X className="w-5 h-5" />
              </button>
            </div>

            {/* Scrollable check list */}
            <div className="flex-1 overflow-y-auto px-5 py-4 space-y-6">
              {Object.entries(DEEP_CHECK_DOMAINS).map(([domainKey, domain]) => (
                <div key={domainKey}>
                  <h4 className="text-xs font-semibold text-[var(--color-text-muted)] uppercase tracking-wider mb-3">
                    {domain.label}
                  </h4>
                  <div className="space-y-2">
                    {domain.checks.map((checkName) => {
                      const check = checkLookup[checkName];
                      const override = config.severityOverrides[checkName];
                      const defaultSeverity = check?.severity ?? 'info';
                      const effectiveSeverity = override ?? defaultSeverity;
                      const isOverridden = Boolean(override);

                      return (
                        <div
                          key={checkName}
                          className={cn(
                            'p-3 rounded-xl border transition-all',
                            isOverridden
                              ? 'border-[var(--color-primary)]/30 bg-[var(--color-primary)]/5'
                              : 'border-[var(--color-border)] bg-[var(--color-surface-sunken)]'
                          )}
                        >
                          <div className="flex items-start gap-2 mb-2">
                            {/* Override indicator */}
                            <div
                              className={cn(
                                'w-1.5 h-1.5 rounded-full mt-1.5 flex-shrink-0 transition-colors',
                                isOverridden
                                  ? 'bg-[var(--color-primary)]'
                                  : 'bg-[var(--color-text-muted)]/30'
                              )}
                            />
                            <div className="flex-1 min-w-0">
                              <div className="flex items-center gap-2 flex-wrap">
                                <span className="text-xs font-medium text-[var(--color-text-primary)]">
                                  {formatCheckName(checkName)}
                                </span>
                                <Badge
                                  variant={getSeverityVariant(effectiveSeverity)}
                                  size="sm"
                                >
                                  {effectiveSeverity.toUpperCase()}
                                </Badge>
                                {isOverridden && (
                                  <span className="text-[10px] text-[var(--color-text-muted)] italic">
                                    overridden
                                  </span>
                                )}
                              </div>
                            </div>
                          </div>

                          {/* Severity toggle buttons */}
                          <div className="flex items-center gap-1 ml-3.5">
                            {SEVERITY_OPTIONS.map((sev) => {
                              const isActive = effectiveSeverity === sev;
                              const isDefault = defaultSeverity === sev && !isOverridden;
                              return (
                                <button
                                  key={sev}
                                  onClick={() => {
                                    if (isActive && isOverridden) {
                                      onRemoveOverride(checkName);
                                    } else {
                                      onSetSeverity(checkName, sev);
                                    }
                                  }}
                                  className={cn(
                                    'px-2 py-1 text-[10px] font-medium rounded-md transition-all',
                                    'border',
                                    isActive && !isDefault
                                      ? sev === 'error'
                                        ? 'bg-red-500/15 text-red-600 dark:text-red-400 border-red-500/30'
                                        : sev === 'warning'
                                        ? 'bg-amber-500/15 text-amber-600 dark:text-amber-400 border-amber-500/30'
                                        : 'bg-sky-500/15 text-sky-600 dark:text-sky-400 border-sky-500/30'
                                      : isActive && isDefault
                                      ? 'bg-[var(--color-surface-elevated)] text-[var(--color-text-primary)] border-[var(--color-border-strong)]'
                                      : 'bg-transparent text-[var(--color-text-muted)] border-[var(--color-border)] hover:border-[var(--color-border-strong)] hover:text-[var(--color-text-secondary)]'
                                  )}
                                >
                                  {sev.toUpperCase()}
                                </button>
                              );
                            })}
                            {isOverridden && (
                              <button
                                onClick={() => onRemoveOverride(checkName)}
                                className={cn(
                                  'ml-1 px-2 py-1 text-[10px] font-medium rounded-md transition-all',
                                  'text-[var(--color-text-muted)] hover:text-[var(--color-text-secondary)]',
                                  'border border-[var(--color-border)] hover:border-[var(--color-border-strong)]'
                                )}
                              >
                                Reset
                              </button>
                            )}
                          </div>
                        </div>
                      );
                    })}
                  </div>
                </div>
              ))}
            </div>

            {/* Footer */}
            <div className="px-5 py-4 border-t border-[var(--color-border)]">
              <button
                onClick={onResetAll}
                disabled={!hasOverrides}
                className={cn(
                  'w-full flex items-center justify-center gap-2',
                  'px-4 py-2.5 rounded-xl text-sm font-medium transition-all',
                  'border',
                  hasOverrides
                    ? 'border-[var(--color-border)] text-[var(--color-text-secondary)] hover:bg-[var(--color-surface-sunken)] hover:text-[var(--color-text-primary)]'
                    : 'border-[var(--color-border)] text-[var(--color-text-muted)] opacity-50 cursor-not-allowed'
                )}
              >
                <RotateCcw className="w-4 h-4" />
                Reset All Overrides
              </button>
              <p className="text-[10px] text-[var(--color-text-muted)] text-center mt-2">
                {hasOverrides
                  ? `${Object.keys(config.severityOverrides).length} override(s) active — saved in localStorage`
                  : 'No overrides active — using backend defaults'}
              </p>
            </div>
          </motion.div>
        </>
      )}
    </AnimatePresence>
  );
}
