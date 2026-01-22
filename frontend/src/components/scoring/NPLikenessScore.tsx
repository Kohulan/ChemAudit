import type { NPLikenessResult } from '../../types/scoring';

interface NPLikenessScoreProps {
  result: NPLikenessResult;
}

/**
 * Displays NP-likeness score with visual scale and interpretation.
 */
export function NPLikenessScore({ result }: NPLikenessScoreProps) {
  const { score, interpretation, caveats } = result;

  // Scale is -5 to +5, position as percentage (0-100)
  // -5 = 0%, 0 = 50%, +5 = 100%
  const markerPosition = ((score + 5) / 10) * 100;
  const clampedPosition = Math.max(0, Math.min(100, markerPosition));

  // Color based on score
  const getScoreColor = () => {
    if (score >= 1.0) return 'text-green-600';
    if (score >= -0.3) return 'text-gray-700';
    return 'text-orange-600';
  };

  return (
    <div className="bg-white rounded-lg border border-gray-200 p-6">
      <div className="flex items-center justify-between mb-4">
        <h3 className="text-lg font-semibold text-gray-900">NP-Likeness Score</h3>
        <div className={`text-2xl font-bold ${getScoreColor()}`}>
          {score >= 0 ? '+' : ''}{score.toFixed(2)}
        </div>
      </div>

      {/* Visual Scale */}
      <div className="mb-4">
        <div className="relative h-8">
          {/* Background gradient bar */}
          <div
            className="absolute inset-x-0 top-3 h-2 rounded-full"
            style={{
              background: 'linear-gradient(to right, #f97316, #fbbf24, #a3a3a3, #86efac, #22c55e)'
            }}
          />

          {/* Marker */}
          <div
            className="absolute top-0 transform -translate-x-1/2 transition-all duration-500"
            style={{ left: `${clampedPosition}%` }}
          >
            <div className="w-0 h-0 border-l-[6px] border-r-[6px] border-t-[8px] border-l-transparent border-r-transparent border-t-gray-800" />
            <div className="w-3 h-3 bg-gray-800 rounded-full mt-[-2px] ml-[-1.5px]" />
          </div>
        </div>

        {/* Scale labels */}
        <div className="flex justify-between text-xs text-gray-500 mt-2">
          <span>-5 (Synthetic)</span>
          <span>0 (Mixed)</span>
          <span>+5 (Natural Product)</span>
        </div>
      </div>

      {/* Interpretation */}
      <p className="text-sm text-gray-600 mb-4">{interpretation}</p>

      {/* Caveats */}
      {caveats.length > 0 && (
        <div className="space-y-1">
          {caveats.map((caveat, index) => (
            <div key={index} className="flex items-start gap-2 text-xs text-amber-700 bg-amber-50 rounded px-2 py-1">
              <span className="text-amber-500">&#9888;</span>
              <span>{caveat}</span>
            </div>
          ))}
        </div>
      )}
    </div>
  );
}
