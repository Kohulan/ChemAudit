/**
 * ProvenanceTimeline Component
 *
 * Vertical timeline rendering per-stage standardization provenance.
 * Shows each pipeline stage as an expandable ProvenanceStageCard with
 * status icons and change counts.
 *
 * Falls back to StepsList when provenance is null (handled in parent).
 * Shows tautomer summary and stereo warning banners when applicable.
 */
import { useState } from 'react';
import { ProvenanceStageCard } from './ProvenanceStageCard';
import type { StandardizationProvenance } from '../../types/standardization';

interface ProvenanceTimelineProps {
  provenance: StandardizationProvenance;
}

export function ProvenanceTimeline({ provenance }: ProvenanceTimelineProps) {
  // Track which stages are expanded (all collapsed by default, except first with changes)
  const [expandedStages, setExpandedStages] = useState<Set<number>>(() => {
    // Auto-expand any stage that has changes
    const initial = new Set<number>();
    provenance.stages.forEach((stage, idx) => {
      const hasChanges =
        stage.charge_changes.length > 0 ||
        stage.bond_changes.length > 0 ||
        stage.ring_changes.length > 0 ||
        stage.radical_changes.length > 0 ||
        stage.fragment_removals.length > 0;
      if (hasChanges) initial.add(idx);
    });
    return initial;
  });

  const toggleStage = (idx: number) => {
    setExpandedStages((prev) => {
      const next = new Set(prev);
      if (next.has(idx)) {
        next.delete(idx);
      } else {
        next.add(idx);
      }
      return next;
    });
  };

  const hasTautomer = Boolean(provenance.tautomer);
  const hasStereoWarning =
    provenance.stereo_summary?.stereo_stripped ||
    (provenance.stereo_summary && provenance.stereo_summary.centers_lost > 0) ||
    (provenance.stereo_summary && provenance.stereo_summary.bonds_lost > 0);

  return (
    <div className="space-y-3">
      {/* Stereo warning banner (when stereocenters were lost) */}
      {hasStereoWarning && provenance.stereo_summary && (
        <div className="bg-amber-50 border-l-4 border-amber-400 p-3 rounded-r-lg">
          <div className="flex items-start gap-2">
            <span className="text-amber-500 font-bold flex-shrink-0">!</span>
            <div>
              <p className="text-sm font-semibold text-amber-800">
                Stereochemistry Changed During Standardization
              </p>
              <div className="flex flex-wrap gap-3 mt-1 text-xs text-amber-700">
                {provenance.stereo_summary.centers_lost > 0 && (
                  <span>{provenance.stereo_summary.centers_lost} stereocenter(s) lost</span>
                )}
                {provenance.stereo_summary.bonds_lost > 0 && (
                  <span>{provenance.stereo_summary.bonds_lost} E/Z bond(s) lost</span>
                )}
                {provenance.stereo_summary.per_center.length > 0 && (
                  <span>{provenance.stereo_summary.per_center.length} center detail(s) available</span>
                )}
              </div>
            </div>
          </div>
        </div>
      )}

      {/* Tautomer summary banner */}
      {hasTautomer && provenance.tautomer && (
        <div className="bg-blue-50 border border-blue-200 rounded-lg p-3">
          <div className="flex items-center justify-between flex-wrap gap-2">
            <div>
              <p className="text-xs font-semibold text-blue-800">Tautomer Canonicalization</p>
              <p className="text-xs text-blue-600 mt-0.5">
                {provenance.tautomer.num_tautomers_enumerated} tautomer
                {provenance.tautomer.num_tautomers_enumerated !== 1 ? 's' : ''} enumerated
                {provenance.tautomer.modified_atoms.length > 0 && (
                  <span> · {provenance.tautomer.modified_atoms.length} atom(s) modified</span>
                )}
              </p>
            </div>
            <div className="flex items-center gap-2">
              {provenance.tautomer.stereo_stripped && (
                <span className="inline-flex items-center px-2 py-0.5 rounded-full text-xs font-medium bg-amber-100 text-amber-700">
                  Stereo stripped
                </span>
              )}
              {provenance.tautomer.complexity_flag && (
                <span className="inline-flex items-center px-2 py-0.5 rounded-full text-xs font-medium bg-orange-100 text-orange-700">
                  High complexity (&gt;100 tautomers)
                </span>
              )}
            </div>
          </div>
        </div>
      )}

      {/* Vertical timeline of stages */}
      <div className="relative">
        {/* Vertical connector line */}
        <div
          className="absolute left-[22px] top-5 bottom-5 w-0.5 bg-gray-200 z-0"
          aria-hidden="true"
        />

        <div className="space-y-2 relative z-10">
          {provenance.stages.map((stage, idx) => (
            <div key={`${stage.stage_name}-${idx}`} className="flex items-start gap-3">
              {/* Timeline dot */}
              <div className="flex-shrink-0 w-[18px] h-[18px] mt-[13px] rounded-full border-2 border-gray-300 bg-white z-10" />

              {/* Stage card */}
              <div className="flex-1 min-w-0">
                <ProvenanceStageCard
                  stage={stage}
                  isExpanded={expandedStages.has(idx)}
                  onToggle={() => toggleStage(idx)}
                />
              </div>
            </div>
          ))}
        </div>
      </div>

      {/* Stage count summary */}
      <p className="text-xs text-gray-400 text-right">
        {provenance.stages.filter((s) => s.applied).length} of {provenance.stages.length} stages applied
      </p>
    </div>
  );
}
