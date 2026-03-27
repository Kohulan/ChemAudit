import { useState, useCallback } from 'react';
import { motion } from 'framer-motion';
import { useProfiler } from '../hooks/useProfiler';
import { ProfilerInput } from '../components/profiler/ProfilerInput';
import { HeroSection } from '../components/profiler/HeroSection';
import { PFIPanel } from '../components/profiler/PFIPanel';
import { StarsPanel } from '../components/profiler/StarsPanel';
import { BioavailabilityPanel } from '../components/profiler/BioavailabilityPanel';
import { ConsensusLogPPanel } from '../components/profiler/ConsensusLogPPanel';
import { SkinPermeationPanel } from '../components/profiler/SkinPermeationPanel';
import { Shape3DPanel } from '../components/profiler/Shape3DPanel';
import { MPOPanel } from '../components/profiler/MPOPanel';
import { LigandEfficiencyPanel } from '../components/profiler/LigandEfficiencyPanel';
import { SAComparisonPanel } from '../components/profiler/SAComparisonPanel';

/**
 * Compound Profiler page — /profiler
 *
 * Scrolling single-page layout per UI-SPEC D-01.
 * Plans 05 and 06 will fill in the section panels (HeroSection, Core Metrics,
 * Shape 3D, Custom MPO, Ligand Efficiency, SA Comparison).
 */
export function CompoundProfilerPage() {
  const {
    profile,
    isLoading,
    error,
    profileCompound,
    compute3DShape,
    computeEfficiency,
    computeCustomMPO,
  } = useProfiler();
  const [currentSmiles, setCurrentSmiles] = useState<string>('');

  const handleProfile = useCallback(
    (smiles: string) => {
      setCurrentSmiles(smiles);
      profileCompound(smiles);
    },
    [profileCompound]
  );

  return (
    <motion.div
      className="max-w-7xl mx-auto px-4 py-16"
      initial={{ opacity: 0, y: 20 }}
      animate={{ opacity: 1, y: 0 }}
      transition={{ duration: 0.6, ease: 'easeOut' }}
    >
      {/* Page heading */}
      <div className="mb-8">
        <h1 className="text-3xl font-semibold font-display text-text-primary mb-1">
          Compound Profiler
        </h1>
        <p className="text-sm text-text-secondary">
          Comprehensive multi-parameter profiling: PFI, #stars, bioavailability, LogP,
          skin permeation, 3D shape, MPO, ligand efficiency, and synthetic accessibility.
        </p>
      </div>

      {/* Input section */}
      <ProfilerInput onSubmit={handleProfile} isLoading={isLoading} />

      {/* API error */}
      {error && (
        <motion.div
          initial={{ opacity: 0, y: -4 }}
          animate={{ opacity: 1, y: 0 }}
          className="mt-4 p-4 rounded-xl bg-status-error/10 border border-status-error/20 text-status-error text-sm"
        >
          {error === 'Profile failed'
            ? 'Unable to compute profile. Check that the structure is valid and try again.'
            : error}
        </motion.div>
      )}

      {/* Profile result sections */}
      {profile && (
        <div className="mt-8 space-y-8">
          {/* Hero section: 2D structure + 6-axis property radar */}
          <HeroSection smiles={currentSmiles} profile={profile} />

          {/* Core Metrics section */}
          <section className="mt-12 space-y-2">
            <h2 className="text-2xl font-semibold font-display mb-6">Core Metrics</h2>
            <div className="grid grid-cols-1 md:grid-cols-2 lg:grid-cols-3 gap-6">
              <PFIPanel data={profile.pfi} />
              <StarsPanel data={profile.stars} />
              <BioavailabilityPanel data={profile.abbott} />
            </div>
            <div className="grid grid-cols-1 md:grid-cols-2 gap-6 mt-6">
              <ConsensusLogPPanel data={profile.consensus_logp} />
              <SkinPermeationPanel data={profile.skin_permeation} />
            </div>
          </section>

          {/* Shape & 3D — per D-26 collapsed, lazy compute */}
          <section className="mt-12">
            <Shape3DPanel smiles={currentSmiles} compute3DShape={compute3DShape} />
          </section>

          {/* Custom MPO */}
          <section className="mt-12">
            <MPOPanel
              smiles={currentSmiles}
              cnsMPO={profile.cns_mpo}
              computeCustomMPO={computeCustomMPO}
            />
          </section>

          {/* Ligand Efficiency — per D-09 collapsed */}
          <section className="mt-12">
            <LigandEfficiencyPanel smiles={currentSmiles} computeEfficiency={computeEfficiency} />
          </section>

          {/* SA Comparison */}
          <section className="mt-12">
            <SAComparisonPanel data={profile.sa_comparison} />
          </section>
        </div>
      )}

      {/* Empty state */}
      {!profile && !isLoading && !error && (
        <div className="mt-16 text-center text-text-muted">
          <h2 className="text-2xl font-semibold font-display mb-2 text-text-primary">
            Enter a molecule to begin profiling
          </h2>
          <p className="text-sm">
            Paste a SMILES, InChI, CAS number, ChEMBL ID, PubChem CID, or DrugBank ID above.
            Try aspirin:{' '}
            <code className="font-mono text-xs bg-surface-sunken px-1 py-0.5 rounded">
              CC(=O)Oc1ccccc1C(=O)O
            </code>
          </p>
        </div>
      )}
    </motion.div>
  );
}
