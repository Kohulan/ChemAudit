import { useState, useEffect } from 'react';
import { motion } from 'framer-motion';
import { BarChart3, AlertTriangle } from 'lucide-react';
import { scoringApi } from '../../services/api';
import { ConsensusScoreCard } from './ConsensusScoreCard';
import { LeadFragmentCard } from './LeadFragmentCard';
import { PropertyBreakdownCard } from './PropertyBreakdownCard';
import { BioavailabilityCard } from './BioavailabilityCard';
import { AtomContributionViewer } from './AtomContributionViewer';
import { cn } from '../../lib/utils';
import type { ScoringResponse, ScoringError } from '../../types/scoring';

interface ScoringProfilesTabProps {
  smiles: string;
}

function LoadingSkeleton() {
  return (
    <div className="space-y-4 animate-pulse">
      {[1, 2, 3, 4].map((i) => (
        <div key={i} className="rounded-2xl bg-[var(--color-surface-sunken)] p-5">
          <div className="h-4 w-32 rounded bg-[var(--color-border)] mb-3" />
          <div className="h-20 rounded-xl bg-[var(--color-border)]/50" />
        </div>
      ))}
    </div>
  );
}

export function ScoringProfilesTab({ smiles }: ScoringProfilesTabProps) {
  const [data, setData] = useState<ScoringResponse | null>(null);
  const [loading, setLoading] = useState(false);
  const [error, setError] = useState<string | null>(null);

  useEffect(() => {
    if (!smiles.trim()) {
      setData(null);
      setError(null);
      return;
    }

    let cancelled = false;

    async function fetchData() {
      setLoading(true);
      setError(null);
      try {
        const response = await scoringApi.getScoringProfiles(smiles.trim());
        if (!cancelled) {
          setData(response);
        }
      } catch (err) {
        if (!cancelled) {
          const scoringErr = err as ScoringError;
          setError(scoringErr?.error || 'Failed to load scoring profiles');
        }
      } finally {
        if (!cancelled) {
          setLoading(false);
        }
      }
    }

    fetchData();

    return () => {
      cancelled = true;
    };
  }, [smiles]);

  if (!smiles.trim()) {
    return (
      <div className="text-center py-12">
        <BarChart3 className="w-12 h-12 mx-auto text-[var(--color-text-muted)] mb-3 opacity-50" />
        <p className="text-sm text-[var(--color-text-muted)]">
          Enter a molecule above to see scoring profiles
        </p>
      </div>
    );
  }

  if (loading) {
    return <LoadingSkeleton />;
  }

  if (error) {
    return (
      <div className="flex items-start gap-3 p-4 rounded-xl bg-red-500/10 border border-red-500/20">
        <AlertTriangle className="w-5 h-5 text-red-500 flex-shrink-0 mt-0.5" />
        <div>
          <p className="text-sm font-medium text-red-600 dark:text-red-400">Scoring Error</p>
          <p className="text-xs text-[var(--color-text-secondary)] mt-1">{error}</p>
        </div>
      </div>
    );
  }

  if (!data) return null;

  const cardVariants = {
    hidden: { opacity: 0, y: 20 },
    visible: (i: number) => ({
      opacity: 1,
      y: 0,
      transition: { delay: i * 0.1, duration: 0.4, ease: [0.25, 0.46, 0.45, 0.94] },
    }),
  };

  return (
    <div className="space-y-5">
      {/* Consensus Score Card */}
      {data.consensus && (
        <motion.div
          custom={0}
          initial="hidden"
          animate="visible"
          variants={cardVariants}
          className={cn(
            'rounded-2xl p-5',
            'bg-gradient-to-br from-[var(--color-surface-elevated)] to-[var(--color-surface)]',
            'border border-[var(--color-border)]'
          )}
        >
          <div className="flex items-center gap-3 mb-4">
            <div className="w-10 h-10 rounded-xl bg-gradient-to-br from-[var(--color-primary)]/10 to-[var(--color-accent)]/10 flex items-center justify-center text-[var(--color-primary)]">
              <BarChart3 className="w-5 h-5" />
            </div>
            <div>
              <h3 className="font-semibold text-[var(--color-text-primary)] text-sm tracking-tight">
                Consensus Score
              </h3>
              <p className="text-xs text-[var(--color-text-muted)]">
                Drug-likeness across 5 rule sets
              </p>
            </div>
          </div>
          <ConsensusScoreCard data={data.consensus} />
        </motion.div>
      )}

      {/* Lead & Fragment Card */}
      {(data.lead_likeness || data.salt_inventory || data.ligand_efficiency) && (
        <motion.div
          custom={1}
          initial="hidden"
          animate="visible"
          variants={cardVariants}
          className={cn(
            'rounded-2xl p-5',
            'bg-gradient-to-br from-[var(--color-surface-elevated)] to-[var(--color-surface)]',
            'border border-[var(--color-border)]'
          )}
        >
          <div className="flex items-center gap-3 mb-4">
            <div className="w-10 h-10 rounded-xl bg-gradient-to-br from-teal-500/10 to-emerald-500/10 flex items-center justify-center text-teal-500">
              <svg className="w-5 h-5" viewBox="0 0 24 24" fill="none" stroke="currentColor" strokeWidth="2">
                <path d="M10 2v7.527a2 2 0 0 1-.211.896L4.72 20.55a1 1 0 0 0 .9 1.45h12.76a1 1 0 0 0 .9-1.45l-5.069-10.127A2 2 0 0 1 14 9.527V2" />
                <path d="M8.5 2h7" />
              </svg>
            </div>
            <div>
              <h3 className="font-semibold text-[var(--color-text-primary)] text-sm tracking-tight">
                Lead & Fragment Profile
              </h3>
              <p className="text-xs text-[var(--color-text-muted)]">
                Lead-likeness, Ro3, salts, and ligand efficiency
              </p>
            </div>
          </div>
          <LeadFragmentCard
            leadLikeness={data.lead_likeness}
            ro3={data.druglikeness?.ro3 || null}
            saltInventory={data.salt_inventory}
            ligandEfficiency={data.ligand_efficiency}
          />
        </motion.div>
      )}

      {/* Property Breakdown Card */}
      {(data.tpsa_breakdown || data.logp_breakdown || data.bertz_detail || data.fsp3_detail) && (
        <motion.div
          custom={2}
          initial="hidden"
          animate="visible"
          variants={cardVariants}
          className={cn(
            'rounded-2xl p-5',
            'bg-gradient-to-br from-[var(--color-surface-elevated)] to-[var(--color-surface)]',
            'border border-[var(--color-border)]'
          )}
        >
          <div className="flex items-center gap-3 mb-4">
            <div className="w-10 h-10 rounded-xl bg-gradient-to-br from-purple-500/10 to-pink-500/10 flex items-center justify-center text-purple-500">
              <svg className="w-5 h-5" viewBox="0 0 24 24" fill="none" stroke="currentColor" strokeWidth="2">
                <path d="M3 3v18h18" />
                <path d="M7 16l4-8 4 4 4-10" />
              </svg>
            </div>
            <div>
              <h3 className="font-semibold text-[var(--color-text-primary)] text-sm tracking-tight">
                Property Breakdown
              </h3>
              <p className="text-xs text-[var(--color-text-muted)]">
                TPSA, LogP, complexity, and 3D character
              </p>
            </div>
          </div>
          <PropertyBreakdownCard
            tpsa={data.tpsa_breakdown}
            logp={data.logp_breakdown}
            bertz={data.bertz_detail}
            fsp3={data.fsp3_detail}
          />

          {/* Atom Contribution Viewers */}
          {data.tpsa_breakdown && data.tpsa_breakdown.atom_contributions.length > 0 && (
            <div className="mt-5 pt-4 border-t border-[var(--color-border)]">
              <AtomContributionViewer
                contributions={data.tpsa_breakdown.atom_contributions}
                label="TPSA"
                unit="A^2"
              />
            </div>
          )}
          {data.logp_breakdown && data.logp_breakdown.atom_contributions.length > 0 && (
            <div className="mt-4 pt-4 border-t border-[var(--color-border)]">
              <AtomContributionViewer
                contributions={data.logp_breakdown.atom_contributions}
                label="LogP"
              />
            </div>
          )}
        </motion.div>
      )}

      {/* Bioavailability Card */}
      {(data.bioavailability_radar || data.boiled_egg) && (
        <motion.div
          custom={3}
          initial="hidden"
          animate="visible"
          variants={cardVariants}
          className={cn(
            'rounded-2xl p-5',
            'bg-gradient-to-br from-[var(--color-surface-elevated)] to-[var(--color-surface)]',
            'border border-[var(--color-border)]'
          )}
        >
          <div className="flex items-center gap-3 mb-4">
            <div className="w-10 h-10 rounded-xl bg-gradient-to-br from-indigo-500/10 to-blue-500/10 flex items-center justify-center text-indigo-500">
              <svg className="w-5 h-5" viewBox="0 0 24 24" fill="none" stroke="currentColor" strokeWidth="2">
                <circle cx="12" cy="12" r="10" />
                <path d="M12 2a15.3 15.3 0 0 1 4 10 15.3 15.3 0 0 1-4 10 15.3 15.3 0 0 1-4-10 15.3 15.3 0 0 1 4-10z" />
                <path d="M2 12h20" />
              </svg>
            </div>
            <div>
              <h3 className="font-semibold text-[var(--color-text-primary)] text-sm tracking-tight">
                Bioavailability & Permeation
              </h3>
              <p className="text-xs text-[var(--color-text-muted)]">
                Radar profile and BOILED-Egg classification
              </p>
            </div>
          </div>
          <BioavailabilityCard
            radar={data.bioavailability_radar}
            boiledEgg={data.boiled_egg}
          />
        </motion.div>
      )}

      {/* Execution time */}
      {data.execution_time_ms !== undefined && (
        <p className="text-xs text-[var(--color-text-muted)] text-right">
          Computed in {data.execution_time_ms}ms
        </p>
      )}
    </div>
  );
}
