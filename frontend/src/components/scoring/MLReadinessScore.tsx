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
  const [showCalculationInfo, setShowCalculationInfo] = useState(false);
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

      <p className="text-sm text-gray-600 mb-4">{interpretation}</p>

      {/* How Calculated Info */}
      <div className="mb-6">
        <button
          onClick={() => setShowCalculationInfo(!showCalculationInfo)}
          className="flex items-center text-sm text-blue-600 hover:text-blue-800"
        >
          <svg className="w-4 h-4 mr-1" fill="currentColor" viewBox="0 0 20 20">
            <path fillRule="evenodd" d="M18 10a8 8 0 11-16 0 8 8 0 0116 0zm-7-4a1 1 0 11-2 0 1 1 0 012 0zM9 9a1 1 0 000 2v3a1 1 0 001 1h1a1 1 0 100-2v-3a1 1 0 00-1-1H9z" clipRule="evenodd" />
          </svg>
          How is this calculated?
        </button>
        {showCalculationInfo && (
          <div className="mt-2 p-4 bg-blue-50 rounded-lg text-sm text-gray-700 border border-blue-100">
            <p className="font-medium text-gray-900 mb-2">ML-Readiness Score Formula:</p>
            <ul className="space-y-2 text-xs">
              <li>
                <span className="font-semibold">Descriptors (40 pts max):</span> Percentage of RDKit molecular descriptors that can be successfully calculated. Score = 40 × (successful / total descriptors).
              </li>
              <li>
                <span className="font-semibold">Fingerprints (40 pts max):</span> Ability to generate common fingerprint types:
                <ul className="ml-4 mt-1 text-gray-600">
                  <li>• Morgan fingerprint (radius=2, 2048 bits): 15 pts</li>
                  <li>• MACCS keys (166 bits): 15 pts</li>
                  <li>• Atom pair fingerprint: 10 pts</li>
                </ul>
              </li>
              <li>
                <span className="font-semibold">Size Constraints (20 pts max):</span> Based on molecular weight and atom count:
                <ul className="ml-4 mt-1 text-gray-600">
                  <li>• Optimal (100-900 Da, 3-100 atoms): 20 pts</li>
                  <li>• Acceptable (50-1200 Da, 1-150 atoms): 10 pts</li>
                  <li>• Out of range: 0 pts</li>
                </ul>
              </li>
            </ul>
            <p className="mt-3 text-xs text-gray-500">
              Higher scores indicate better suitability for machine learning models. Most drug-like molecules score 80+.
            </p>
          </div>
        )}
      </div>

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
