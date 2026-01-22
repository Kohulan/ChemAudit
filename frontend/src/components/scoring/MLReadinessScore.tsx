import { useState } from 'react';
import type { MLReadinessResult } from '../../types/scoring';

interface MLReadinessScoreProps {
  result: MLReadinessResult;
}

/**
 * Displays ML-readiness score with breakdown and interpretation.
 */
export function MLReadinessScore({ result }: MLReadinessScoreProps) {
  const [showFailedDescriptors, setShowFailedDescriptors] = useState(false);
  const { score, breakdown, interpretation, failed_descriptors } = result;

  // Color based on score
  const getScoreColor = (score: number) => {
    if (score >= 80) return { text: 'text-green-600', bg: 'bg-green-500' };
    if (score >= 50) return { text: 'text-yellow-600', bg: 'bg-yellow-500' };
    return { text: 'text-red-600', bg: 'bg-red-500' };
  };

  const scoreColor = getScoreColor(score);

  // Calculate bar widths as percentages
  const descriptorsWidth = (breakdown.descriptors_score / breakdown.descriptors_max) * 100;
  const fingerprintsWidth = (breakdown.fingerprints_score / breakdown.fingerprints_max) * 100;
  const sizeWidth = (breakdown.size_score / breakdown.size_max) * 100;

  return (
    <div className="bg-white rounded-lg border border-gray-200 p-6">
      <div className="flex items-center justify-between mb-4">
        <h3 className="text-lg font-semibold text-gray-900">ML-Readiness Score</h3>
        <div className={`text-3xl font-bold ${scoreColor.text}`}>
          {score}<span className="text-lg font-normal text-gray-500">/100</span>
        </div>
      </div>

      <p className="text-sm text-gray-600 mb-6">{interpretation}</p>

      {/* Score Breakdown */}
      <div className="space-y-4">
        {/* Descriptors */}
        <div>
          <div className="flex justify-between text-sm mb-1">
            <span className="text-gray-700">Descriptors</span>
            <span className="text-gray-600">
              {breakdown.descriptors_score.toFixed(0)}/{breakdown.descriptors_max}
              <span className="text-gray-400 ml-1">
                ({breakdown.descriptors_successful}/{breakdown.descriptors_total} calculated)
              </span>
            </span>
          </div>
          <div className="h-2 bg-gray-200 rounded-full overflow-hidden">
            <div
              className={`h-full ${scoreColor.bg} transition-all duration-500`}
              style={{ width: `${descriptorsWidth}%` }}
            />
          </div>
        </div>

        {/* Fingerprints */}
        <div>
          <div className="flex justify-between text-sm mb-1">
            <span className="text-gray-700">Fingerprints</span>
            <span className="text-gray-600">
              {breakdown.fingerprints_score.toFixed(0)}/{breakdown.fingerprints_max}
              <span className="text-gray-400 ml-1">
                ({breakdown.fingerprints_successful.join(', ') || 'none'})
              </span>
            </span>
          </div>
          <div className="h-2 bg-gray-200 rounded-full overflow-hidden">
            <div
              className={`h-full ${scoreColor.bg} transition-all duration-500`}
              style={{ width: `${fingerprintsWidth}%` }}
            />
          </div>
          {breakdown.fingerprints_failed.length > 0 && (
            <p className="text-xs text-red-500 mt-1">
              Failed: {breakdown.fingerprints_failed.join(', ')}
            </p>
          )}
        </div>

        {/* Size */}
        <div>
          <div className="flex justify-between text-sm mb-1">
            <span className="text-gray-700">Size Constraints</span>
            <span className="text-gray-600">
              {breakdown.size_score.toFixed(0)}/{breakdown.size_max}
              <span className="text-gray-400 ml-1">
                ({breakdown.size_category})
              </span>
            </span>
          </div>
          <div className="h-2 bg-gray-200 rounded-full overflow-hidden">
            <div
              className={`h-full ${scoreColor.bg} transition-all duration-500`}
              style={{ width: `${sizeWidth}%` }}
            />
          </div>
          {breakdown.molecular_weight !== null && breakdown.num_atoms !== null && (
            <p className="text-xs text-gray-500 mt-1">
              MW: {breakdown.molecular_weight.toFixed(1)} Da | Atoms: {breakdown.num_atoms}
            </p>
          )}
        </div>
      </div>

      {/* Failed Descriptors (collapsible) */}
      {failed_descriptors.length > 0 && (
        <div className="mt-4 pt-4 border-t border-gray-200">
          <button
            onClick={() => setShowFailedDescriptors(!showFailedDescriptors)}
            className="flex items-center text-sm text-gray-600 hover:text-gray-800"
          >
            <span className={`mr-2 transform transition-transform ${showFailedDescriptors ? 'rotate-90' : ''}`}>
              &#9654;
            </span>
            {failed_descriptors.length} descriptors could not be calculated
          </button>
          {showFailedDescriptors && (
            <div className="mt-2 p-3 bg-gray-50 rounded text-xs text-gray-600 max-h-32 overflow-y-auto">
              {failed_descriptors.join(', ')}
            </div>
          )}
        </div>
      )}
    </div>
  );
}
