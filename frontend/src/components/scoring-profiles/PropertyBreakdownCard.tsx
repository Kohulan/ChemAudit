import { useState } from 'react';
import type { TPSABreakdown, LogPBreakdown, BertzDetail, Fsp3Detail } from '../../types/scoring';
import { cn } from '../../lib/utils';

interface PropertyBreakdownCardProps {
  tpsa: TPSABreakdown | null;
  logp: LogPBreakdown | null;
  bertz: BertzDetail | null;
  fsp3: Fsp3Detail | null;
}

type SubTab = 'tpsa' | 'logp' | 'bertz' | 'fsp3';

function PropertyRow({ label, value, unit }: { label: string; value: string | number; unit?: string }) {
  return (
    <div className="flex items-center justify-between py-1.5 border-b border-[var(--color-border)]/50 last:border-0">
      <span className="text-xs text-[var(--color-text-muted)]">{label}</span>
      <span className="text-sm font-medium text-[var(--color-text-primary)]">
        {typeof value === 'number' ? (Number.isInteger(value) ? String(value) : value.toFixed(2)) : value}
        {unit && <span className="text-xs text-[var(--color-text-muted)] ml-1">{unit}</span>}
      </span>
    </div>
  );
}

function GroupTable({ groups }: { groups: { group_name: string; contribution: number; atom_indices: number[] }[] }) {
  if (groups.length === 0) return null;

  return (
    <div className="rounded-xl overflow-hidden border border-[var(--color-border)] mt-3">
      <table className="w-full text-xs">
        <thead>
          <tr className="bg-[var(--color-surface-sunken)]">
            <th className="text-left px-3 py-2 text-[var(--color-text-muted)] font-medium">Group</th>
            <th className="text-right px-3 py-2 text-[var(--color-text-muted)] font-medium">Contribution</th>
            <th className="text-right px-3 py-2 text-[var(--color-text-muted)] font-medium">Atoms</th>
          </tr>
        </thead>
        <tbody>
          {groups.map((g, i) => (
            <tr key={i} className="border-t border-[var(--color-border)]/50">
              <td className="px-3 py-2 text-[var(--color-text-primary)]">{g.group_name}</td>
              <td className={cn(
                'px-3 py-2 text-right font-mono',
                g.contribution > 0 ? 'text-emerald-600 dark:text-emerald-400' : 'text-[var(--color-text-secondary)]'
              )}>
                {g.contribution > 0 ? '+' : ''}{g.contribution.toFixed(2)}
              </td>
              <td className="px-3 py-2 text-right text-[var(--color-text-muted)]">{g.atom_indices.length}</td>
            </tr>
          ))}
        </tbody>
      </table>
    </div>
  );
}

export function PropertyBreakdownCard({ tpsa, logp, bertz, fsp3 }: PropertyBreakdownCardProps) {
  const tabs: { id: SubTab; label: string; available: boolean }[] = [
    { id: 'tpsa', label: 'TPSA', available: tpsa !== null },
    { id: 'logp', label: 'LogP', available: logp !== null },
    { id: 'bertz', label: 'Bertz', available: bertz !== null },
    { id: 'fsp3', label: 'Fsp3', available: fsp3 !== null },
  ];

  const firstAvailable = tabs.find(t => t.available)?.id || 'tpsa';
  const [activeTab, setActiveTab] = useState<SubTab>(firstAvailable);

  return (
    <div className="space-y-4">
      {/* Sub-Tab Bar */}
      <div className="flex gap-1 p-1 rounded-xl bg-[var(--color-surface-sunken)]">
        {tabs.filter(t => t.available).map((tab) => (
          <button
            key={tab.id}
            onClick={() => setActiveTab(tab.id)}
            className={cn(
              'flex-1 px-3 py-1.5 text-xs font-medium rounded-lg transition-all',
              activeTab === tab.id
                ? 'bg-[var(--color-surface-elevated)] text-[var(--color-primary)] shadow-sm'
                : 'text-[var(--color-text-muted)] hover:text-[var(--color-text-secondary)]'
            )}
          >
            {tab.label}
          </button>
        ))}
      </div>

      {/* TPSA Tab */}
      {activeTab === 'tpsa' && tpsa && (
        <div>
          <div className="rounded-xl bg-[var(--color-surface-sunken)] p-3">
            <PropertyRow label="Total TPSA" value={tpsa.total} unit="A^2" />
            <PropertyRow
              label="Contributing atoms"
              value={tpsa.atom_contributions.filter(a => a.contribution > 0).length}
            />
          </div>
          <GroupTable groups={tpsa.functional_group_summary} />
        </div>
      )}

      {/* LogP Tab */}
      {activeTab === 'logp' && logp && (
        <div>
          <div className="rounded-xl bg-[var(--color-surface-sunken)] p-3">
            <PropertyRow label="Total LogP" value={logp.total} />
            <PropertyRow
              label="Total atoms"
              value={logp.atom_contributions.length}
            />
          </div>
          <GroupTable groups={logp.functional_group_summary} />
        </div>
      )}

      {/* Bertz Tab */}
      {activeTab === 'bertz' && bertz && (
        <div className="rounded-xl bg-[var(--color-surface-sunken)] p-3">
          <PropertyRow label="Bertz CT" value={bertz.bertz_ct} />
          <PropertyRow label="Bonds" value={bertz.num_bonds} />
          <PropertyRow label="Heavy atoms" value={bertz.num_atoms} />
          <PropertyRow label="Rings" value={bertz.num_rings} />
          <PropertyRow label="Aromatic rings" value={bertz.num_aromatic_rings} />
          <PropertyRow label="Ring complexity" value={bertz.ring_complexity} />
          <p className="text-xs text-[var(--color-text-muted)] mt-2 pt-2 border-t border-[var(--color-border)]/50">
            {bertz.interpretation}
          </p>
        </div>
      )}

      {/* Fsp3 Tab */}
      {activeTab === 'fsp3' && fsp3 && (
        <div>
          <div className="rounded-xl bg-[var(--color-surface-sunken)] p-3">
            <PropertyRow label="Fsp3" value={fsp3.fsp3} />
            <PropertyRow label="Total carbons" value={fsp3.total_carbons} />
            <PropertyRow label="sp3 carbons" value={fsp3.sp3_count} />
            <PropertyRow label="sp2 carbons" value={fsp3.sp2_count} />
            <PropertyRow label="sp carbons" value={fsp3.sp_count} />
            <p className="text-xs text-[var(--color-text-muted)] mt-2 pt-2 border-t border-[var(--color-border)]/50">
              {fsp3.interpretation}
            </p>
          </div>

          {/* Carbon hybridization mini-table */}
          {fsp3.per_carbon.length > 0 && fsp3.per_carbon.length <= 30 && (
            <div className="mt-3">
              <p className="text-xs font-medium text-[var(--color-text-muted)] mb-2">Carbon hybridization map</p>
              <div className="flex flex-wrap gap-1">
                {fsp3.per_carbon.map((c) => (
                  <span
                    key={c.atom_index}
                    className={cn(
                      'inline-flex items-center justify-center w-7 h-7 rounded-lg text-xs font-mono font-medium',
                      c.hybridization === 'sp3'
                        ? 'bg-emerald-500/15 text-emerald-600 dark:text-emerald-400'
                        : c.hybridization === 'sp2'
                          ? 'bg-amber-500/15 text-amber-600 dark:text-amber-400'
                          : c.hybridization === 'sp'
                            ? 'bg-blue-500/15 text-blue-600 dark:text-blue-400'
                            : 'bg-gray-500/10 text-gray-500'
                    )}
                    title={`C${c.atom_index}: ${c.hybridization}`}
                  >
                    {c.atom_index}
                  </span>
                ))}
              </div>
              <div className="flex gap-3 mt-2 text-xs text-[var(--color-text-muted)]">
                <span className="flex items-center gap-1">
                  <span className="w-2.5 h-2.5 rounded bg-emerald-500/30" /> sp3
                </span>
                <span className="flex items-center gap-1">
                  <span className="w-2.5 h-2.5 rounded bg-amber-500/30" /> sp2
                </span>
                <span className="flex items-center gap-1">
                  <span className="w-2.5 h-2.5 rounded bg-blue-500/30" /> sp
                </span>
              </div>
            </div>
          )}
        </div>
      )}
    </div>
  );
}
