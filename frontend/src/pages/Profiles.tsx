/**
 * Profiles page — hosts PresetPicker and ProfileBuilder.
 * Accessible at /profiles via lazy-loaded route.
 */

import { useState, useEffect, useCallback } from 'react';
import { motion } from 'framer-motion';
import { SlidersHorizontal } from 'lucide-react';
import { PresetPicker } from '../components/profiles/PresetPicker';
import { ProfileBuilder } from '../components/profiles/ProfileBuilder';
import { profilesApi } from '../services/api';
import type { ScoringProfile, ScoringProfileCreate } from '../types/workflow';

type ActiveSection = 'presets' | 'builder';

/** Scoring Profiles page with preset picker and custom builder. */
export function ProfilesPage() {
  const [activeSection, setActiveSection] = useState<ActiveSection>('presets');
  const [editingProfile, setEditingProfile] = useState<ScoringProfile | null>(null);
  const [presets, setPresets] = useState<ScoringProfile[]>([]);
  const [activeProfileId, setActiveProfileId] = useState<number | null>(null);
  const [isLoading, setIsLoading] = useState(true);
  const [isSaving, setIsSaving] = useState(false);

  const fetchProfiles = useCallback(async () => {
    setIsLoading(true);
    try {
      const profiles = await profilesApi.getProfiles();
      setPresets(profiles);
    } catch {
      // handle silently
    } finally {
      setIsLoading(false);
    }
  }, []);

  useEffect(() => {
    fetchProfiles();
  }, [fetchProfiles]);

  /** Apply a preset as the active profile. */
  const handleApply = useCallback((profile: ScoringProfile) => {
    setActiveProfileId(profile.id);
  }, []);

  /** Duplicate a preset and switch to builder to customize it. */
  const handleDuplicate = useCallback((profile: ScoringProfile) => {
    // Create a copy without ID so it becomes a new profile in builder
    const copy: ScoringProfile = {
      ...profile,
      id: 0,
      name: `${profile.name} (Copy)`,
      is_preset: false,
    };
    setEditingProfile(copy);
    setActiveSection('builder');
  }, []);

  /** Save a new or edited profile. */
  const handleSave = useCallback(async (data: ScoringProfileCreate) => {
    setIsSaving(true);
    try {
      if (editingProfile && editingProfile.id > 0) {
        await profilesApi.updateProfile(editingProfile.id, data);
      } else {
        await profilesApi.createProfile(data);
      }
      await fetchProfiles();
      setActiveSection('presets');
      setEditingProfile(null);
    } finally {
      setIsSaving(false);
    }
  }, [editingProfile, fetchProfiles]);

  /** Cancel builder and return to preset picker. */
  const handleCancel = useCallback(() => {
    setActiveSection('presets');
    setEditingProfile(null);
  }, []);

  return (
    <div className="mx-auto max-w-5xl px-4 sm:px-6 space-y-8 pb-16">
      {/* Page header */}
      <motion.div
        className="text-center pt-4"
        initial={{ opacity: 0, y: -20 }}
        animate={{ opacity: 1, y: 0 }}
        transition={{ duration: 0.6, ease: [0.25, 0.46, 0.45, 0.94] }}
      >
        <div className="inline-flex items-center gap-2 px-4 py-1.5 rounded-full bg-[var(--color-primary)]/10 border border-[var(--color-primary)]/20 mb-4">
          <SlidersHorizontal className="w-4 h-4 text-[var(--color-primary)]" />
          <span className="text-sm font-medium text-[var(--color-primary)]">Customisable Scoring</span>
        </div>
        <h1 className="text-3xl sm:text-4xl font-bold text-gradient tracking-tight font-display">
          Scoring Profiles
        </h1>
        <p className="text-[var(--color-text-secondary)] mt-3 text-base sm:text-lg max-w-2xl mx-auto">
          Apply preset scoring rules or build a custom profile with your own thresholds and weights.
        </p>
      </motion.div>

      {/* Tab bar */}
      <motion.div
        initial={{ opacity: 0, y: 10 }}
        animate={{ opacity: 1, y: 0 }}
        transition={{ delay: 0.1 }}
        className="flex gap-2 border-b border-[var(--color-border)] pb-0"
      >
        <button
          onClick={() => setActiveSection('presets')}
          className={`px-5 py-2.5 text-sm font-medium rounded-t-lg transition-all duration-200 border-b-2 ${
            activeSection === 'presets'
              ? 'text-[var(--color-primary)] border-[var(--color-primary)] bg-[var(--color-primary)]/5'
              : 'text-[var(--color-text-secondary)] border-transparent hover:text-[var(--color-text-primary)]'
          }`}
        >
          Presets
        </button>
        <button
          onClick={() => setActiveSection('builder')}
          className={`px-5 py-2.5 text-sm font-medium rounded-t-lg transition-all duration-200 border-b-2 ${
            activeSection === 'builder'
              ? 'text-[var(--color-primary)] border-[var(--color-primary)] bg-[var(--color-primary)]/5'
              : 'text-[var(--color-text-secondary)] border-transparent hover:text-[var(--color-text-primary)]'
          }`}
        >
          Custom Builder
        </button>
      </motion.div>

      {/* Content */}
      <motion.div
        key={activeSection}
        initial={{ opacity: 0, y: 10 }}
        animate={{ opacity: 1, y: 0 }}
        transition={{ duration: 0.2 }}
      >
        {activeSection === 'presets' && (
          <div className="space-y-4">
            {isLoading ? (
              <div className="flex items-center justify-center py-16">
                <div className="animate-spin rounded-full h-8 w-8 border-b-2 border-[var(--color-primary)]" />
              </div>
            ) : (
              <PresetPicker
                presets={presets}
                activeProfileId={activeProfileId}
                onApply={handleApply}
                onDuplicate={handleDuplicate}
              />
            )}
          </div>
        )}

        {activeSection === 'builder' && (
          <ProfileBuilder
            profile={editingProfile}
            onSave={handleSave}
            onCancel={handleCancel}
            isSaving={isSaving}
          />
        )}
      </motion.div>
    </div>
  );
}
