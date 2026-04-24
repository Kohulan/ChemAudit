import { useState, useCallback } from 'react';
import { AnimatePresence } from 'framer-motion';
import { api } from '../../services/api';
import { useProfiler } from '../../hooks/useProfiler';
import type { ProfileResponse } from '../../types/profiler';
import { HeroSection } from './HeroSection';
import { PFIPanel } from './PFIPanel';
import { StarsPanel } from './StarsPanel';
import { BioavailabilityPanel } from './BioavailabilityPanel';
import { ConsensusLogPPanel } from './ConsensusLogPPanel';
import { SkinPermeationPanel } from './SkinPermeationPanel';
import { DrugLikenessGrid } from './DrugLikenessGrid';
import { SAComparisonPanel } from './SAComparisonPanel';
import { Shape3DPanel } from './Shape3DPanel';
import { LigandEfficiencyPanel } from './LigandEfficiencyPanel';
import { MPOPanel } from './MPOPanel';
import { ComparisonBar, type PinnedMolecule } from './ComparisonBar';
import { ComparisonView } from './ComparisonView';

interface ProfilerAccordionProps {
  /** Current molecule SMILES */
  smiles: string;
  /** Profile result from /api/v1/profiler/full */
  profile: ProfileResponse | null;
  /** Whether the profile is currently loading */
  isLoading: boolean;
  /** Error message if profile fetch failed */
  error: string | null;
}

/**
 * ProfilerAccordion — wraps all profiler standalone components into a single
 * panel suitable for rendering inside a DrillDownSection accordion on
 * SingleValidation.
 *
 * Layout:
 * 1. HeroSection (2D structure + radar)
 * 2. Core metrics grid (PFI, Stars, Bioavailability, ConsensusLogP, SkinPermeation)
 * 3. DrugLikenessGrid (5-card rule pass/fail)
 * 4. SAComparisonPanel (SA/SCScore/SYBA 4-card)
 * 5. Collapsible lazy-compute sections (Shape3D, LigandEfficiency, MPO)
 * 6. ComparisonBar (sticky bottom, pin up to 5)
 */
export function ProfilerAccordion({
  smiles,
  profile,
  isLoading,
  error,
}: ProfilerAccordionProps) {
  const {
    compute3DShape,
    computeEfficiency,
    computeCustomMPO,
  } = useProfiler();

  // Comparison state (multi-molecule pinning)
  const [pinnedMolecules, setPinnedMolecules] = useState<PinnedMolecule[]>([]);
  const [showComparison, setShowComparison] = useState(false);

  const handlePin = useCallback(() => {
    if (!profile || !smiles) return;
    if (pinnedMolecules.length >= 5) return;
    if (pinnedMolecules.some((m) => m.smiles === smiles)) return;
    setPinnedMolecules((prev) => [
      ...prev,
      {
        smiles,
        label: smiles.substring(0, 20),
        profile,
      },
    ]);
  }, [profile, smiles, pinnedMolecules]);

  const handleRemove = useCallback((smilesStr: string) => {
    setPinnedMolecules((prev) => prev.filter((m) => m.smiles !== smilesStr));
  }, []);

  const handleAddMultiple = useCallback(
    async (smilesList: string[]) => {
      const remaining = 5 - pinnedMolecules.length;
      const toAdd = smilesList.slice(0, remaining);

      const results = await Promise.allSettled(
        toAdd.map(async (smi) => {
          const { data } = await api.post<ProfileResponse>('/profiler/full', { smiles: smi });
          return { smiles: smi, label: smi.substring(0, 20), profile: data };
        })
      );

      const newPinned = results
        .filter((r): r is PromiseFulfilledResult<PinnedMolecule> => r.status === 'fulfilled')
        .map((r) => r.value)
        .filter((m) => !pinnedMolecules.some((existing) => existing.smiles === m.smiles));

      if (newPinned.length > 0) {
        setPinnedMolecules((prev) => [...prev, ...newPinned]);
      }
    },
    [pinnedMolecules]
  );

  // Empty state: no SMILES entered
  if (!smiles && !isLoading && !error) {
    return (
      <div className="text-center py-8 text-text-muted">
        <p className="text-sm">Enter a molecule to see its profile</p>
      </div>
    );
  }

  // Loading state: skeleton placeholders
  if (isLoading) {
    return (
      <div className="space-y-6">
        {/* Hero skeleton */}
        <div className="animate-pulse rounded-xl bg-surface-sunken h-48 w-full" />
        {/* Metrics grid skeleton */}
        <div className="grid grid-cols-1 md:grid-cols-2 lg:grid-cols-3 gap-4">
          {Array.from({ length: 5 }).map((_, i) => (
            <div
              key={i}
              className="animate-pulse rounded-xl bg-surface-sunken h-32"
            />
          ))}
        </div>
        {/* Drug-likeness skeleton */}
        <div className="animate-pulse rounded-xl bg-surface-sunken h-24 w-full" />
      </div>
    );
  }

  // Error state
  if (error) {
    return (
      <div className="p-4 rounded-xl bg-status-error/10 border border-status-error/20 text-status-error text-sm">
        {error}
      </div>
    );
  }

  // No profile data yet (but smiles is set — waiting for data)
  if (!profile) {
    return null;
  }

  // Data loaded: full profiler content
  return (
    <div className="space-y-8">
      {/* Property radar chart */}
      <HeroSection smiles={smiles} profile={profile} onPin={handlePin} />

      {/* Core Metrics grid */}
      <section>
        <h3 className="text-lg font-semibold font-display mb-5">Core Metrics</h3>
        <div className="grid grid-cols-1 md:grid-cols-2 gap-5">
          <PFIPanel data={profile.pfi} />
          <StarsPanel data={profile.stars} />
          <BioavailabilityPanel data={profile.abbott} />
          <ConsensusLogPPanel data={profile.consensus_logp} />
          <SkinPermeationPanel data={profile.skin_permeation} className="md:col-span-2" />
        </div>
      </section>

      {/* 3. Drug-Likeness Rules Grid */}
      {profile.druglikeness && (
        <section>
          <h3 className="text-lg font-semibold font-display mb-4">Drug-Likeness Rules</h3>
          <DrugLikenessGrid druglikeness={profile.druglikeness} />
        </section>
      )}

      {/* 4. SA Comparison */}
      <section>
        <SAComparisonPanel data={profile.sa_comparison} />
      </section>

      {/* 5. Collapsible lazy-compute sections */}
      <section>
        <Shape3DPanel smiles={smiles} compute3DShape={compute3DShape} />
      </section>

      <section>
        <LigandEfficiencyPanel smiles={smiles} computeEfficiency={computeEfficiency} />
      </section>

      <section>
        <MPOPanel
          smiles={smiles}
          cnsMPO={profile.cns_mpo}
          computeCustomMPO={computeCustomMPO}
        />
      </section>

      {/* ComparisonBar — fixed bottom strip when >= 1 molecule pinned */}
      {pinnedMolecules.length > 0 && (
        <ComparisonBar
          molecules={pinnedMolecules}
          onRemove={handleRemove}
          onCompare={() => setShowComparison(true)}
          onAddMultiple={handleAddMultiple}
        />
      )}

      {/* Comparison sidebar — slides in from right */}
      <AnimatePresence>
        {showComparison && pinnedMolecules.length >= 2 && (
          <ComparisonView
            molecules={pinnedMolecules}
            onClose={() => setShowComparison(false)}
          />
        )}
      </AnimatePresence>
    </div>
  );
}
