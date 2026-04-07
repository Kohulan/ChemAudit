import { CheckCircle2, XCircle } from 'lucide-react';

/**
 * Drug-likeness rule result shape from backend druglikeness scorer.
 */
interface DruglikenessRuleResult {
  passed: boolean;
  violations?: number;
  violations_count?: number;
  num_violations?: number;
  [key: string]: unknown;
}

interface DrugLikenessGridProps {
  druglikeness: Record<string, unknown>;
}

/**
 * The 5 rules shown in the compact grid (per D-03).
 * Backend keys match the scoring module output from druglikeness.py.
 */
const RULES = [
  { key: 'lipinski', label: 'Lipinski' },
  { key: 'veber', label: 'Veber' },
  { key: 'egan', label: 'Egan' },
  { key: 'muegge', label: 'Muegge' },
  { key: 'ghose', label: 'Ghose' },
] as const;

/**
 * Extract violation count from a rule result object.
 * Backend may use different field names across scorers.
 */
function getViolationCount(result: DruglikenessRuleResult): number | null {
  if (typeof result.violations === 'number') return result.violations;
  if (typeof result.violations_count === 'number') return result.violations_count;
  if (typeof result.num_violations === 'number') return result.num_violations;
  return null;
}

/**
 * DrugLikenessGrid — compact pass/fail grid for 5 drug-likeness rules.
 *
 * Per D-03: inline on profiler page, no expand/collapse needed.
 * Rules: Lipinski, Veber, Egan, Muegge, Ghose.
 */
export function DrugLikenessGrid({ druglikeness }: DrugLikenessGridProps) {
  return (
    <div className="grid grid-cols-2 md:grid-cols-3 lg:grid-cols-5 gap-2">
      {RULES.map(({ key, label }) => {
        const raw = druglikeness[key];

        // Handle missing / unavailable rule data
        if (!raw || typeof raw !== 'object') {
          return (
            <div
              key={key}
              className="clay-card-sm p-3 flex flex-col items-center gap-1 text-center"
            >
              <span className="text-xs font-semibold font-display text-text-secondary">
                {label}
              </span>
              <span className="text-[10px] text-text-muted mt-0.5">N/A</span>
            </div>
          );
        }

        const result = raw as DruglikenessRuleResult;
        const passed = !!result.passed;
        const violationCount = getViolationCount(result);

        return (
          <div
            key={key}
            className={[
              'clay-card-sm p-3 flex flex-col items-center gap-1 text-center',
              'border',
              passed
                ? 'border-emerald-500/20 bg-emerald-500/5 dark:bg-emerald-500/10'
                : 'border-rose-500/20 bg-rose-500/5 dark:bg-rose-500/10',
            ].join(' ')}
          >
            {/* Rule name */}
            <span className="text-xs font-semibold font-display text-text-primary">
              {label}
            </span>

            {/* Pass / Fail icon */}
            {passed ? (
              <CheckCircle2 className="w-4 h-4 text-emerald-500" aria-label="Pass" />
            ) : (
              <XCircle className="w-4 h-4 text-rose-500" aria-label="Fail" />
            )}

            {/* Violation count */}
            <span
              className={[
                'text-[10px] leading-tight',
                passed ? 'text-emerald-600 dark:text-emerald-400' : 'text-rose-600 dark:text-rose-400',
              ].join(' ')}
            >
              {violationCount !== null
                ? violationCount === 0
                  ? '0 violations'
                  : `${violationCount} violation${violationCount === 1 ? '' : 's'}`
                : passed
                ? 'Pass'
                : 'Fail'}
            </span>
          </div>
        );
      })}
    </div>
  );
}
