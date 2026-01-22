import { MLReadinessScore } from './MLReadinessScore';
import { NPLikenessScore } from './NPLikenessScore';
import { ScaffoldDisplay } from './ScaffoldDisplay';
import type { ScoringResponse } from '../../types/scoring';

interface ScoringResultsProps {
  scoringResponse: ScoringResponse;
}

/**
 * Combined display for all scoring results.
 */
export function ScoringResults({ scoringResponse }: ScoringResultsProps) {
  const { ml_readiness, np_likeness, scaffold, execution_time_ms } = scoringResponse;

  return (
    <div className="space-y-6">
      {/* Header */}
      <div className="flex items-center justify-between">
        <h2 className="text-xl font-semibold text-gray-900">Scoring Results</h2>
        <span className="text-sm text-gray-500">
          Completed in {execution_time_ms}ms
        </span>
      </div>

      {/* Results grid */}
      <div className="grid grid-cols-1 lg:grid-cols-2 gap-6">
        {/* ML-Readiness (full width on larger screens for detail) */}
        {ml_readiness && (
          <div className="lg:col-span-2">
            <MLReadinessScore result={ml_readiness} />
          </div>
        )}

        {/* NP-Likeness */}
        {np_likeness && (
          <NPLikenessScore result={np_likeness} />
        )}

        {/* Scaffold */}
        {scaffold && (
          <ScaffoldDisplay result={scaffold} />
        )}
      </div>
    </div>
  );
}
