import type { NPLikenessResult } from '../../types/scoring';
import { CalculationTooltip, InfoTooltip } from '../ui/Tooltip';

interface NPLikenessScoreProps {
  result: NPLikenessResult;
}

/**
 * Displays NP-likeness score with visual scale, marker, and informative tooltips.
 */
export function NPLikenessScore({ result }: NPLikenessScoreProps) {
  const { score, interpretation, caveats } = result;

  // Scale is -5 to +5, position as percentage (0-100)
  // -5 = 0%, 0 = 50%, +5 = 100%
  const markerPosition = ((score + 5) / 10) * 100;
  const clampedPosition = Math.max(0, Math.min(100, markerPosition));

  // Color and label based on score
  const getScoreInfo = () => {
    if (score >= 1.0) return {
      text: 'text-score-excellent',
      bg: 'bg-score-excellent/10',
      label: 'Natural Product-like',
      color: '#059669',
    };
    if (score >= -0.3) return {
      text: 'text-chem-dark/70',
      bg: 'bg-chem-dark/5',
      label: 'Mixed Character',
      color: '#64748b',
    };
    return {
      text: 'text-status-warning',
      bg: 'bg-status-warning/10',
      label: 'Synthetic-like',
      color: '#d97706',
    };
  };

  const scoreInfo = getScoreInfo();

  const calculation = `NP-likeness score based on fragment analysis.
Range: -5 (synthetic) to +5 (natural product)
Score: ${score.toFixed(2)}`;

  const fullInterpretation = score >= 1.0
    ? 'This molecule has structural features commonly found in natural products. It may have evolved drug-like properties and good bioavailability.'
    : score >= -0.3
    ? 'This molecule has mixed characteristics of both synthetic and natural compounds. It combines features from both chemical spaces.'
    : 'This molecule has features more common in synthetic compounds. This is typical for many successful drugs and is not inherently negative.';

  return (
    <div className="card-chem p-6">
      {/* Header */}
      <div className="flex items-start justify-between mb-6">
        <div className="flex items-center gap-3">
          <div className="section-header-icon bg-emerald-100 text-emerald-600">
            <svg className="w-5 h-5" viewBox="0 0 24 24" fill="none" stroke="currentColor" strokeWidth="2">
              <path d="M12 2a10 10 0 0110 10c0 5.523-4.477 10-10 10S2 17.523 2 12 6.477 2 12 2z" />
              <path d="M12 2c-2.761 0-5 4.477-5 10s2.239 10 5 10 5-4.477 5-10-2.239-10-5-10z" />
              <path d="M2 12h20M12 2v20" />
            </svg>
          </div>
          <div>
            <h3 className="text-lg font-semibold text-chem-dark flex items-center gap-2">
              NP-Likeness Score
              <InfoTooltip
                title="What is NP-Likeness?"
                content={
                  <div className="space-y-2 text-xs">
                    <p>NP-likeness measures how similar a molecule is to natural products.</p>
                    <ul className="list-disc list-inside space-y-1 text-white/70">
                      <li>Positive scores: Natural product-like</li>
                      <li>Around zero: Mixed character</li>
                      <li>Negative scores: Synthetic-like</li>
                    </ul>
                    <p className="text-white/50 text-xs mt-2">Based on fragment analysis comparing to known natural products.</p>
                  </div>
                }
              />
            </h3>
            <p className="text-sm text-chem-dark/50">Natural product similarity</p>
          </div>
        </div>

        {/* Score Display */}
        <CalculationTooltip
          calculation={calculation}
          interpretation={fullInterpretation}
          title="NP-Likeness Score"
          value={`${score >= 0 ? '+' : ''}${score.toFixed(2)}`}
          position="left"
        >
          <div className={`text-center px-4 py-3 rounded-xl ${scoreInfo.bg}`}>
            <div className={`text-3xl font-bold ${scoreInfo.text}`}>
              {score >= 0 ? '+' : ''}{score.toFixed(2)}
            </div>
            <div className={`text-xs mt-1 ${scoreInfo.text}`}>{scoreInfo.label}</div>
          </div>
        </CalculationTooltip>
      </div>

      {/* Visual Scale */}
      <div className="mb-6">
        <div className="relative h-12 pt-2">
          {/* Background gradient bar */}
          <div
            className="absolute inset-x-0 top-4 h-4 rounded-full shadow-inner"
            style={{
              background: 'linear-gradient(to right, #d97706 0%, #fbbf24 25%, #94a3b8 50%, #86efac 75%, #059669 100%)'
            }}
          />

          {/* Tick marks */}
          {[-5, -2.5, 0, 2.5, 5].map((tick) => {
            const pos = ((tick + 5) / 10) * 100;
            return (
              <div
                key={tick}
                className="absolute top-4 w-0.5 h-4 bg-white/50"
                style={{ left: `${pos}%` }}
              />
            );
          })}

          {/* Marker */}
          <div
            className="absolute top-0 transform -translate-x-1/2 transition-all duration-700 ease-out"
            style={{ left: `${clampedPosition}%` }}
          >
            {/* Arrow */}
            <div
              className="w-0 h-0 border-l-[8px] border-r-[8px] border-t-[10px] border-l-transparent border-r-transparent"
              style={{ borderTopColor: scoreInfo.color }}
            />
            {/* Dot */}
            <div
              className="w-4 h-4 rounded-full -mt-1 shadow-lg"
              style={{ backgroundColor: scoreInfo.color, marginLeft: '-2px' }}
            />
          </div>
        </div>

        {/* Scale labels */}
        <div className="flex justify-between text-xs mt-3">
          <span className="text-status-warning font-medium">-5 Synthetic</span>
          <span className="text-chem-dark/50">0 Mixed</span>
          <span className="text-score-excellent font-medium">+5 Natural</span>
        </div>
      </div>

      {/* Interpretation */}
      <div className="bg-chem-primary/5 rounded-xl p-4 mb-4">
        <p className="text-sm text-chem-dark/80">{interpretation}</p>
      </div>

      {/* Caveats */}
      {caveats.length > 0 && (
        <div className="space-y-2">
          <h4 className="text-xs font-semibold text-chem-dark/50 uppercase tracking-wide">
            Important Notes
          </h4>
          {caveats.map((caveat, index) => (
            <div
              key={index}
              className="flex items-start gap-2 text-xs text-status-warning bg-status-warning-light/50 rounded-lg px-3 py-2"
            >
              <svg className="w-4 h-4 flex-shrink-0 mt-0.5" viewBox="0 0 24 24" fill="currentColor">
                <path d="M12 2L1 21h22L12 2zm0 3.17L20.12 19H3.88L12 5.17zM11 10v4h2v-4h-2zm0 6v2h2v-2h-2z" />
              </svg>
              <span>{caveat}</span>
            </div>
          ))}
        </div>
      )}
    </div>
  );
}
