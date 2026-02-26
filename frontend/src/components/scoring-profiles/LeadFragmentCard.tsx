import { FlaskConical, Beaker } from 'lucide-react';
import type { LeadLikeness, RuleOfThreeResult, SaltInventory, LigandEfficiency } from '../../types/scoring';
import { cn } from '../../lib/utils';

interface LeadFragmentCardProps {
  leadLikeness: LeadLikeness | null;
  ro3: RuleOfThreeResult | null;
  saltInventory: SaltInventory | null;
  ligandEfficiency: LigandEfficiency | null;
}

function PropertyRow({
  label,
  value,
  threshold,
  passed
}: {
  label: string;
  value: string | number;
  threshold?: string;
  passed?: boolean;
}) {
  return (
    <div className="flex items-center justify-between py-1.5 border-b border-[var(--color-border)]/50 last:border-0">
      <span className="text-xs text-[var(--color-text-muted)]">{label}</span>
      <div className="flex items-center gap-2">
        <span className={cn(
          'text-sm font-medium',
          passed === undefined
            ? 'text-[var(--color-text-primary)]'
            : passed
              ? 'text-emerald-600 dark:text-emerald-400'
              : 'text-red-600 dark:text-red-400'
        )}>
          {typeof value === 'number' ? value.toFixed(2) : value}
        </span>
        {threshold && (
          <span className="text-xs text-[var(--color-text-muted)]">({threshold})</span>
        )}
      </div>
    </div>
  );
}

function CategoryBadge({ category }: { category: string }) {
  const colors: Record<string, string> = {
    counterion: 'bg-blue-500/10 text-blue-600 dark:text-blue-400',
    salt: 'bg-purple-500/10 text-purple-600 dark:text-purple-400',
    solvent: 'bg-amber-500/10 text-amber-600 dark:text-amber-400',
    drug: 'bg-emerald-500/10 text-emerald-600 dark:text-emerald-400',
    unknown: 'bg-gray-500/10 text-gray-600 dark:text-gray-400',
  };

  return (
    <span className={cn(
      'text-xs px-2 py-0.5 rounded-full font-medium',
      colors[category] || colors.unknown
    )}>
      {category}
    </span>
  );
}

export function LeadFragmentCard({ leadLikeness, ro3, saltInventory, ligandEfficiency }: LeadFragmentCardProps) {
  return (
    <div className="space-y-5">
      {/* Lead-Likeness */}
      {leadLikeness && (
        <div>
          <div className="flex items-center gap-2 mb-3">
            <FlaskConical className="w-4 h-4 text-[var(--color-primary)]" />
            <h4 className="text-sm font-semibold text-[var(--color-text-primary)]">Lead-Likeness</h4>
            <span className={cn(
              'text-xs px-2 py-0.5 rounded-full font-medium',
              leadLikeness.passed
                ? 'bg-emerald-500/10 text-emerald-600 dark:text-emerald-400'
                : 'bg-red-500/10 text-red-600 dark:text-red-400'
            )}>
              {leadLikeness.passed ? 'PASS' : `${leadLikeness.violations} violation${leadLikeness.violations !== 1 ? 's' : ''}`}
            </span>
          </div>
          <div className="rounded-xl bg-[var(--color-surface-sunken)] p-3">
            {leadLikeness.violation_details.map((v) => (
              <PropertyRow
                key={v.property}
                label={v.property}
                value={v.value}
                threshold={v.threshold}
                passed={v.result === 'pass'}
              />
            ))}
          </div>
        </div>
      )}

      {/* Fragment-Likeness (Ro3) */}
      {ro3 && (
        <div>
          <div className="flex items-center gap-2 mb-3">
            <Beaker className="w-4 h-4 text-[var(--color-accent)]" />
            <h4 className="text-sm font-semibold text-[var(--color-text-primary)]">Fragment-Likeness (Ro3)</h4>
            <span className={cn(
              'text-xs px-2 py-0.5 rounded-full font-medium',
              ro3.passed
                ? 'bg-emerald-500/10 text-emerald-600 dark:text-emerald-400'
                : 'bg-red-500/10 text-red-600 dark:text-red-400'
            )}>
              {ro3.passed ? 'PASS' : `${ro3.violations} violation${ro3.violations !== 1 ? 's' : ''}`}
            </span>
          </div>
          <div className="rounded-xl bg-[var(--color-surface-sunken)] p-3">
            <PropertyRow label="MW" value={ro3.mw} threshold="<=300" passed={ro3.mw <= 300} />
            <PropertyRow label="LogP" value={ro3.logp} threshold="<=3" passed={ro3.logp <= 3} />
            <PropertyRow label="HBD" value={ro3.hbd} threshold="<=3" passed={ro3.hbd <= 3} />
            <PropertyRow label="HBA" value={ro3.hba} threshold="<=3" passed={ro3.hba <= 3} />
            <PropertyRow label="RotBonds" value={ro3.rotatable_bonds} threshold="<=3" passed={ro3.rotatable_bonds <= 3} />
            <PropertyRow label="TPSA" value={ro3.tpsa} threshold="<=60" passed={ro3.tpsa <= 60} />
          </div>
        </div>
      )}

      {/* Salt Inventory */}
      {saltInventory && saltInventory.has_salts && (
        <div>
          <div className="flex items-center gap-2 mb-3">
            <h4 className="text-sm font-semibold text-[var(--color-text-primary)]">Salt Inventory</h4>
            <span className="text-xs px-2 py-0.5 rounded-full font-medium bg-blue-500/10 text-blue-600 dark:text-blue-400">
              {saltInventory.total_fragments} fragment{saltInventory.total_fragments !== 1 ? 's' : ''}
            </span>
          </div>
          <div className="rounded-xl overflow-hidden border border-[var(--color-border)]">
            <table className="w-full text-xs">
              <thead>
                <tr className="bg-[var(--color-surface-sunken)]">
                  <th className="text-left px-3 py-2 text-[var(--color-text-muted)] font-medium">SMILES</th>
                  <th className="text-left px-3 py-2 text-[var(--color-text-muted)] font-medium">Name</th>
                  <th className="text-left px-3 py-2 text-[var(--color-text-muted)] font-medium">Category</th>
                  <th className="text-right px-3 py-2 text-[var(--color-text-muted)] font-medium">MW</th>
                </tr>
              </thead>
              <tbody>
                {saltInventory.fragments.map((f, i) => (
                  <tr key={i} className="border-t border-[var(--color-border)]/50">
                    <td className="px-3 py-2 font-mono text-[var(--color-text-secondary)]">{f.smiles}</td>
                    <td className="px-3 py-2 text-[var(--color-text-primary)]">{f.name}</td>
                    <td className="px-3 py-2"><CategoryBadge category={f.category} /></td>
                    <td className="px-3 py-2 text-right text-[var(--color-text-secondary)]">{f.mw.toFixed(1)}</td>
                  </tr>
                ))}
              </tbody>
            </table>
          </div>
          <p className="text-xs text-[var(--color-text-muted)] mt-2">{saltInventory.interpretation}</p>
        </div>
      )}

      {/* Ligand Efficiency */}
      {ligandEfficiency && (
        <div>
          <div className="flex items-center gap-2 mb-3">
            <h4 className="text-sm font-semibold text-[var(--color-text-primary)]">Ligand Efficiency</h4>
            {ligandEfficiency.proxy_used && (
              <span className="text-xs px-2 py-0.5 rounded-full font-medium bg-amber-500/10 text-amber-600 dark:text-amber-400">
                BEI proxy
              </span>
            )}
          </div>
          <div className="rounded-xl bg-[var(--color-surface-sunken)] p-3">
            {ligandEfficiency.le !== null && (
              <div className="flex items-center justify-between py-1.5">
                <span className="text-xs text-[var(--color-text-muted)]">LE value</span>
                <span className="text-sm font-bold text-[var(--color-text-primary)]">
                  {ligandEfficiency.le.toFixed(3)}
                </span>
              </div>
            )}
            <div className="flex items-center justify-between py-1.5 border-t border-[var(--color-border)]/50">
              <span className="text-xs text-[var(--color-text-muted)]">Heavy atoms</span>
              <span className="text-sm text-[var(--color-text-primary)]">{ligandEfficiency.heavy_atom_count}</span>
            </div>
          </div>
          <p className="text-xs text-[var(--color-text-muted)] mt-2">{ligandEfficiency.interpretation}</p>
        </div>
      )}
    </div>
  );
}
