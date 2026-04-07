import { useState, useCallback } from 'react';
import { useProfiler } from '../../hooks/useProfiler';
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
import type { ProfileResponse } from '../../types/profiler';

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
 * Layout mirrors CompoundProfiler page:
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
      // Stub for comparison bar batch paste — limited use in accordion context
      void smilesList;
    },
    []
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
      {/* 1. Hero section: 2D structure + 6-axis property radar */}
      <HeroSection smiles={smiles} profile={profile} onPin={handlePin} />

      {/* 2. Core Metrics grid */}
      <section className="space-y-2">
        <h3 className="text-lg font-semibold font-display mb-4">Core Metrics</h3>
        <div className="grid grid-cols-1 md:grid-cols-2 lg:grid-cols-3 gap-4">
          <PFIPanel data={profile.pfi} />
          <StarsPanel data={profile.stars} />
          <BioavailabilityPanel data={profile.abbott} />
        </div>
        <div className="grid grid-cols-1 md:grid-cols-2 gap-4 mt-4">
          <ConsensusLogPPanel data={profile.consensus_logp} />
          <SkinPermeationPanel data={profile.skin_permeation} />
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

      {/* 6. ComparisonBar — sticky bottom strip when >= 1 molecule pinned */}
      {pinnedMolecules.length > 0 && (
        <ComparisonBar
          molecules={pinnedMolecules}
          onRemove={handleRemove}
          onCompare={() => {}}
          onAddMultiple={handleAddMultiple}
        />
      )}
    </div>
  );
}
