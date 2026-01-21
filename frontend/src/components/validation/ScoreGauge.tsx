interface ScoreGaugeProps {
  score: number;
  size?: number;
  className?: string;
}

export function ScoreGauge({ score, size = 120, className = '' }: ScoreGaugeProps) {
  // Clamp score to 0-100
  const clampedScore = Math.max(0, Math.min(100, score));

  // Determine color based on score
  const getColor = (score: number) => {
    if (score >= 80) return { stroke: '#10b981', text: 'text-green-600', bg: 'bg-green-50' };
    if (score >= 50) return { stroke: '#f59e0b', text: 'text-yellow-600', bg: 'bg-yellow-50' };
    return { stroke: '#ef4444', text: 'text-red-600', bg: 'bg-red-50' };
  };

  const color = getColor(clampedScore);

  // Calculate circle properties
  const strokeWidth = 10;
  const radius = (size - strokeWidth) / 2;
  const circumference = 2 * Math.PI * radius;
  const offset = circumference - (clampedScore / 100) * circumference;

  return (
    <div className={`flex flex-col items-center ${className}`}>
      <div className="relative" style={{ width: size, height: size }}>
        <svg
          width={size}
          height={size}
          className="transform -rotate-90"
        >
          {/* Background circle */}
          <circle
            cx={size / 2}
            cy={size / 2}
            r={radius}
            stroke="#e5e7eb"
            strokeWidth={strokeWidth}
            fill="none"
          />
          {/* Progress circle */}
          <circle
            cx={size / 2}
            cy={size / 2}
            r={radius}
            stroke={color.stroke}
            strokeWidth={strokeWidth}
            fill="none"
            strokeDasharray={circumference}
            strokeDashoffset={offset}
            strokeLinecap="round"
            className="transition-all duration-500 ease-out"
          />
        </svg>
        {/* Score text in center */}
        <div className="absolute inset-0 flex flex-col items-center justify-center">
          <span className={`text-3xl font-bold ${color.text}`}>
            {Math.round(clampedScore)}
          </span>
          <span className="text-xs text-gray-500 uppercase tracking-wide">Score</span>
        </div>
      </div>
      <div className={`mt-3 px-4 py-1.5 rounded-full ${color.bg}`}>
        <span className={`text-sm font-medium ${color.text}`}>
          {clampedScore >= 80 && 'Excellent'}
          {clampedScore >= 50 && clampedScore < 80 && 'Fair'}
          {clampedScore < 50 && 'Poor'}
        </span>
      </div>
    </div>
  );
}
