import { useState, useCallback } from 'react';
import { motion } from 'framer-motion';
import { useProfiler } from '../hooks/useProfiler';
import { ProfilerInput } from '../components/profiler/ProfilerInput';

/**
 * Compound Profiler page — /profiler
 *
 * Scrolling single-page layout per UI-SPEC D-01.
 * Plans 05 and 06 will fill in the section panels (HeroSection, Core Metrics,
 * Shape 3D, Custom MPO, Ligand Efficiency, SA Comparison).
 */
export function CompoundProfilerPage() {
  const { profile, isLoading, error, profileCompound } = useProfiler();
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
          {/* HeroSection — Plan 05 will implement */}
          {/* Core Metrics section (PFI, Stars, Bioavailability, Consensus LogP, Skin Permeation) — Plan 05 will implement */}
          {/* Shape 3D section — Plan 06 will implement */}
          {/* Custom MPO section — Plan 06 will implement */}
          {/* Ligand Efficiency section — Plan 06 will implement */}
          {/* SA Comparison section — Plan 06 will implement */}

          {/* Temporary debug display — will be replaced by section panels */}
          <div className="text-xs font-mono text-text-muted p-4 bg-surface-sunken rounded-xl">
            Profile loaded for: <span className="text-text-primary">{currentSmiles.slice(0, 60)}{currentSmiles.length > 60 ? '...' : ''}</span>
          </div>
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
