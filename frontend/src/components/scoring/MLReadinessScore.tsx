import { useState } from 'react';
import type { MLReadinessResult } from '../../types/scoring';
import { ScoreChart, ScoreBreakdownBar } from './ScoreChart';
import { InfoTooltip } from '../ui/Tooltip';

interface MLReadinessScoreProps {
  result: MLReadinessResult;
}

/**
 * Displays ML-readiness score with radial chart, breakdown bars, and informative tooltips.
 */
export function MLReadinessScore({ result }: MLReadinessScoreProps) {
  const [showFailedDescriptors, setShowFailedDescriptors] = useState(false);
  const { score, breakdown, interpretation, failed_descriptors } = result;

  // Calculation explanation
  const calculation = `Score = Descriptors (${breakdown.descriptors_max}pts) + Fingerprints (${breakdown.fingerprints_max}pts) + Size (${breakdown.size_max}pts)
= ${breakdown.descriptors_score.toFixed(0)} + ${breakdown.fingerprints_score.toFixed(0)} + ${breakdown.size_score.toFixed(0)} = ${score}`;

  return (
    <div className="card-chem p-6">
      {/* Header */}
      <div className="flex items-start justify-between mb-6">
        <div className="flex items-center gap-3">
          <div className="section-header-icon">
            <svg className="w-5 h-5" viewBox="0 0 24 24" fill="none" stroke="currentColor" strokeWidth="2">
              <path d="M12 2v4M12 18v4M4.93 4.93l2.83 2.83M16.24 16.24l2.83 2.83M2 12h4M18 12h4M4.93 19.07l2.83-2.83M16.24 7.76l2.83-2.83" />
            </svg>
          </div>
          <div>
            <h3 className="text-lg font-semibold text-chem-dark">ML-Readiness Score</h3>
            <p className="text-sm text-chem-dark/50">Suitability for machine learning models</p>
          </div>
        </div>

        {/* Main Score */}
        <ScoreChart
          score={score}
          label="ML-Readiness"
          size={120}
          calculation={calculation}
          interpretation={interpretation}
          compact
        />
      </div>

      {/* Interpretation */}
      <div className="bg-chem-primary/5 rounded-xl p-4 mb-6">
        <p className="text-sm text-chem-dark/80">{interpretation}</p>
      </div>

      {/* Score Breakdown */}
      <div className="space-y-5">
        <h4 className="text-sm font-semibold text-chem-dark/70 uppercase tracking-wide flex items-center gap-2">
          Score Breakdown
          <InfoTooltip
            title="How ML-Readiness is Calculated"
            content={
              <div className="space-y-2 text-xs">
                <p>The ML-Readiness score measures how suitable a molecule is for machine learning models.</p>
                <ul className="list-disc list-inside space-y-1 text-white/70">
                  <li>Descriptors (40pts): % of RDKit descriptors calculable</li>
                  <li>Fingerprints (40pts): Ability to generate Morgan, MACCS, AtomPair</li>
                  <li>Size (20pts): Molecular weight and atom count constraints</li>
                </ul>
              </div>
            }
          />
        </h4>

        {/* Descriptors */}
        <ScoreBreakdownBar
          label="Molecular Descriptors"
          score={breakdown.descriptors_score}
          maxScore={breakdown.descriptors_max}
          detail={`${breakdown.descriptors_successful}/${breakdown.descriptors_total} calculated`}
          calculation={`Score = ${breakdown.descriptors_max} x (successful / total descriptors)
= ${breakdown.descriptors_max} x (${breakdown.descriptors_successful} / ${breakdown.descriptors_total})
= ${breakdown.descriptors_score.toFixed(1)}`}
          interpretation={breakdown.descriptors_score >= 35
            ? 'Most molecular descriptors can be calculated successfully. This molecule is well-suited for descriptor-based ML models.'
            : breakdown.descriptors_score >= 20
            ? 'Some descriptors failed to calculate. This may limit compatibility with certain ML models.'
            : 'Many descriptors cannot be calculated. Consider structural issues that may affect ML model performance.'}
        />

        {/* Fingerprints */}
        <ScoreBreakdownBar
          label="Molecular Fingerprints"
          score={breakdown.fingerprints_score}
          maxScore={breakdown.fingerprints_max}
          detail={breakdown.fingerprints_successful.join(', ') || 'none'}
          calculation={`Score based on fingerprint generation:
- Morgan (radius=2): 15 pts
- MACCS keys: 15 pts
- Atom pair: 10 pts
Successful: ${breakdown.fingerprints_successful.join(', ') || 'none'}`}
          interpretation={breakdown.fingerprints_score >= 35
            ? 'All major fingerprint types can be generated. Compatible with fingerprint-based similarity and ML methods.'
            : breakdown.fingerprints_score >= 20
            ? 'Some fingerprints types failed. May affect compatibility with certain similarity methods.'
            : 'Major fingerprint generation issues. Check for unusual atoms or bonds.'}
        />
        {breakdown.fingerprints_failed.length > 0 && (
          <p className="text-xs text-status-error ml-4 -mt-3">
            Failed: {breakdown.fingerprints_failed.join(', ')}
          </p>
        )}

        {/* Size Constraints */}
        <ScoreBreakdownBar
          label="Size Constraints"
          score={breakdown.size_score}
          maxScore={breakdown.size_max}
          detail={breakdown.size_category}
          calculation={`Based on molecular weight and atom count:
- Optimal (100-900 Da, 3-100 atoms): 20 pts
- Acceptable (50-1200 Da, 1-150 atoms): 10 pts
- Out of range: 0 pts
Category: ${breakdown.size_category}`}
          interpretation={breakdown.size_score >= 15
            ? 'Molecule is within optimal size range for most ML models and drug-likeness criteria.'
            : breakdown.size_score >= 5
            ? 'Molecule is at the edge of typical size ranges. Some models may handle it differently.'
            : 'Molecule is outside typical size ranges for drug-like compounds. May require specialized models.'}
        />
        {breakdown.molecular_weight !== null && breakdown.num_atoms !== null && (
          <p className="text-xs text-chem-dark/50 ml-4 -mt-3">
            MW: {breakdown.molecular_weight.toFixed(1)} Da | Atoms: {breakdown.num_atoms}
          </p>
        )}
      </div>

      {/* Failed Descriptors (collapsible) */}
      {failed_descriptors.length > 0 && (
        <div className="mt-6 pt-4 border-t border-chem-dark/10">
          <button
            onClick={() => setShowFailedDescriptors(!showFailedDescriptors)}
            className="flex items-center gap-2 text-sm text-chem-dark/60 hover:text-chem-dark transition-colors"
          >
            <svg
              className={`w-4 h-4 transition-transform ${showFailedDescriptors ? 'rotate-90' : ''}`}
              viewBox="0 0 24 24"
              fill="none"
              stroke="currentColor"
              strokeWidth="2"
            >
              <path d="M9 18l6-6-6-6" />
            </svg>
            <span>{failed_descriptors.length} descriptors could not be calculated</span>
          </button>
          {showFailedDescriptors && (
            <div className="mt-3 p-4 bg-chem-dark/5 rounded-lg text-xs text-chem-dark/60 font-mono max-h-32 overflow-y-auto scrollbar-thin">
              {failed_descriptors.join(', ')}
            </div>
          )}
        </div>
      )}
    </div>
  );
}
