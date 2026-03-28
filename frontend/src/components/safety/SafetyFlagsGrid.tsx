import { MetricCard } from '../profiler/MetricCard';
import { CypSoftSpotsPanel } from './CypSoftSpotsPanel';
import type { CypResult, HergResult, Bro5Result, ReosResult, Bro5Violation, ReosViolation } from '../../types/safety';

interface SafetyFlagsGridProps {
  cyp: CypResult;
  herg: HergResult;
  bro5: Bro5Result;
  reos: ReosResult;
  onHighlightAtoms: (atoms: number[], alertName?: string) => void;
  highlightedAlert: string | null;
}

/**
 * Formats a descriptor key from snake_case or camelCase to readable title.
 */
function formatKey(key: string): string {
  return key
    .replace(/_/g, ' ')
    .replace(/([a-z])([A-Z])/g, '$1 $2')
    .replace(/\b\w/g, (c) => c.toUpperCase());
}

/**
 * Capitalize first letter of string.
 */
function capitalize(str: string): string {
  return str.charAt(0).toUpperCase() + str.slice(1);
}

/**
 * hERG risk detail: flags + descriptor key-value table.
 */
function HergDetail({ flags, descriptors }: { flags: string[]; descriptors: Record<string, number | boolean> }) {
  return (
    <div className="space-y-3">
      {flags.length > 0 && (
        <div>
          <p className="text-xs font-semibold text-text-secondary mb-1.5">Risk Flags</p>
          <ul className="space-y-1">
            {flags.map((flag, i) => (
              <li key={i} className="text-sm text-text-primary flex items-start gap-1.5">
                <span className="text-status-warning mt-0.5 shrink-0">&#8226;</span>
                <span>{flag}</span>
              </li>
            ))}
          </ul>
        </div>
      )}
      {Object.keys(descriptors).length > 0 && (
        <div>
          <p className="text-xs font-semibold text-text-secondary mb-1.5">Descriptors</p>
          <div className="grid grid-cols-2 gap-x-4 gap-y-1">
            {Object.entries(descriptors).map(([key, value]) => (
              <>
                <span key={`${key}-label`} className="text-xs text-text-muted">{formatKey(key)}</span>
                <span key={`${key}-value`} className="text-xs text-text-primary font-medium tabular-nums">
                  {typeof value === 'boolean' ? (value ? 'Yes' : 'No') : typeof value === 'number' ? value.toFixed(2) : String(value)}
                </span>
              </>
            ))}
          </div>
        </div>
      )}
    </div>
  );
}

/**
 * bRo5 detail: violations table + property values.
 */
function Bro5Detail({ violations, values }: { violations: Bro5Violation[]; values: Record<string, number> }) {
  return (
    <div className="space-y-3">
      {violations.length > 0 && (
        <div>
          <p className="text-xs font-semibold text-text-secondary mb-1.5">Violations</p>
          <div className="space-y-1">
            {violations.map((v, i) => (
              <div key={i} className="grid grid-cols-3 gap-x-2 text-xs">
                <span className="text-text-muted">{formatKey(v.property)}</span>
                <span className="text-text-primary font-medium tabular-nums">{v.value.toFixed(2)}</span>
                <span className="text-text-muted">(threshold: {v.threshold})</span>
              </div>
            ))}
          </div>
        </div>
      )}
      {Object.keys(values).length > 0 && (
        <div>
          <p className="text-xs font-semibold text-text-secondary mb-1.5">Values</p>
          <div className="grid grid-cols-2 gap-x-4 gap-y-1">
            {Object.entries(values).map(([key, value]) => (
              <>
                <span key={`${key}-label`} className="text-xs text-text-muted">{formatKey(key)}</span>
                <span key={`${key}-value`} className="text-xs text-text-primary font-medium tabular-nums">
                  {typeof value === 'number' ? value.toFixed(2) : String(value)}
                </span>
              </>
            ))}
          </div>
        </div>
      )}
    </div>
  );
}

/**
 * REOS detail: violations table + descriptor key-value.
 */
function ReosDetail({ violations, descriptors }: { violations: ReosViolation[]; descriptors: Record<string, number> }) {
  return (
    <div className="space-y-3">
      {violations.length > 0 && (
        <div>
          <p className="text-xs font-semibold text-text-secondary mb-1.5">Violations</p>
          <div className="space-y-1">
            {violations.map((v, i) => (
              <div key={i} className="grid grid-cols-3 gap-x-2 text-xs">
                <span className="text-text-muted">{formatKey(v.property)}</span>
                <span className="text-text-primary font-medium tabular-nums">{v.value.toFixed(2)}</span>
                <span className="text-text-muted">
                  [{v.range[0]}, {v.range[1]}]{v.exceeded ? ' (exceeded)' : ' (below min)'}
                </span>
              </div>
            ))}
          </div>
        </div>
      )}
      {Object.keys(descriptors).length > 0 && (
        <div>
          <p className="text-xs font-semibold text-text-secondary mb-1.5">Descriptors</p>
          <div className="grid grid-cols-2 gap-x-4 gap-y-1">
            {Object.entries(descriptors).map(([key, value]) => (
              <>
                <span key={`${key}-label`} className="text-xs text-text-muted">{formatKey(key)}</span>
                <span key={`${key}-value`} className="text-xs text-text-primary font-medium tabular-nums">
                  {typeof value === 'number' ? value.toFixed(2) : String(value)}
                </span>
              </>
            ))}
          </div>
        </div>
      )}
    </div>
  );
}

/**
 * Safety Flags Grid — 2x2 MetricCard grid (D-13).
 *
 * Renders four MetricCard instances for CYP, hERG, bRo5, and REOS.
 * CYP is always 'default' variant (informational — no binary pass/fail per Pitfall 5).
 * bRo5 always renders, showing N/A when MW <= 500 (D-15).
 */
export function SafetyFlagsGrid({
  cyp,
  herg,
  bro5,
  reos,
  onHighlightAtoms,
  highlightedAlert,
}: SafetyFlagsGridProps) {
  // hERG traffic-light
  const hergVariant =
    herg.herg_risk === 'low' ? 'success' : herg.herg_risk === 'moderate' ? 'warning' : 'error';

  // bRo5 traffic-light
  const bro5Variant = !bro5.applicable ? 'default' : bro5.passed ? 'success' : 'error';
  const bro5Score = !bro5.applicable
    ? 'N/A'
    : bro5.passed
      ? 'Passed'
      : `${bro5.violations.length} violation(s)`;
  const bro5Classification = !bro5.applicable
    ? 'Not applicable'
    : bro5.passed
      ? 'Pass'
      : `${bro5.violations.length} violation(s)`;

  // REOS traffic-light
  const reosVariant =
    reos.n_violations === 0 ? 'success' : reos.n_violations === 1 ? 'warning' : 'error';

  return (
    <div className="space-y-4">
      <h2 className="text-lg font-semibold font-display">Safety Flags</h2>

      <div className="grid grid-cols-1 sm:grid-cols-2 gap-6">
        {/* CYP Soft-Spots — always 'default', never warning/error (Pitfall 5) */}
        <MetricCard
          title="CYP Soft-Spots (Informational)"
          score={cyp.n_sites}
          classification={cyp.n_sites > 0 ? `${cyp.n_sites} site(s)` : 'None detected'}
          classificationVariant="default"
          defaultExpanded={true}
        >
          <CypSoftSpotsPanel
            sites={cyp.sites}
            onHighlightAtoms={onHighlightAtoms}
            highlightedAlert={highlightedAlert}
          />
        </MetricCard>

        {/* hERG Risk */}
        <MetricCard
          title="hERG Risk"
          score={`${herg.risk_score}/${herg.max_score}`}
          classification={capitalize(herg.herg_risk)}
          classificationVariant={hergVariant}
          defaultExpanded={true}
        >
          <HergDetail flags={herg.flags} descriptors={herg.descriptors} />
        </MetricCard>

        {/* bRo5 — ALWAYS rendered, never hidden (D-15) */}
        <MetricCard
          title="bRo5 (Beyond Rule of 5)"
          score={bro5Score}
          classification={bro5Classification}
          classificationVariant={bro5Variant}
          defaultExpanded={true}
        >
          {!bro5.applicable ? (
            <p className="text-text-muted text-sm">N/A — MW &le; 500, use standard Ro5 for this molecule</p>
          ) : (
            <Bro5Detail violations={bro5.violations} values={bro5.values} />
          )}
        </MetricCard>

        {/* REOS Filter */}
        <MetricCard
          title="REOS Filter"
          score={reos.n_violations === 0 ? 'Passed' : `${reos.n_violations} violation(s)`}
          classification={reos.passed ? 'Pass' : `${reos.n_violations} violation(s)`}
          classificationVariant={reosVariant}
          defaultExpanded={true}
        >
          <ReosDetail violations={reos.violations} descriptors={reos.descriptors} />
        </MetricCard>
      </div>
    </div>
  );
}
