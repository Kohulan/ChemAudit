import {
  BarChart,
  Bar,
  XAxis,
  YAxis,
  CartesianGrid,
  ReferenceLine,
  ResponsiveContainer,
  Tooltip,
} from 'recharts';
import type { PropertyDistributions, HistogramData } from '../../types/dataset_intelligence';

// =============================================================================
// Types
// =============================================================================

interface PropertyDistOverlayProps {
  /** Property distribution histograms (MW, LogP, TPSA). */
  distributions: PropertyDistributions;
}

// =============================================================================
// Drug-space reference ranges (Lipinski / Veber)
// =============================================================================

interface PropertyConfig {
  key: keyof PropertyDistributions;
  label: string;
  unit: string;
  refLow: number;
  refHigh: number;
  refLabel: string;
}

const PROPERTY_CONFIGS: PropertyConfig[] = [
  {
    key: 'mw',
    label: 'Molecular Weight',
    unit: 'Da',
    refLow: 150,
    refHigh: 500,
    refLabel: 'Lipinski (150-500)',
  },
  {
    key: 'logp',
    label: 'LogP',
    unit: '',
    refLow: -0.4,
    refHigh: 5,
    refLabel: 'Lipinski (-0.4 to 5)',
  },
  {
    key: 'tpsa',
    label: 'TPSA',
    unit: '\u00C5\u00B2',
    refLow: 20,
    refHigh: 130,
    refLabel: 'Veber (20-130)',
  },
];

// =============================================================================
// Helpers
// =============================================================================

/** Convert histogram bins + counts to bar chart data. */
function histogramToBarData(hist: HistogramData): Array<{ bin: string; binCenter: number; count: number }> {
  const result: Array<{ bin: string; binCenter: number; count: number }> = [];
  for (let i = 0; i < hist.counts.length; i++) {
    const lo = hist.bins[i];
    const hi = hist.bins[i + 1] ?? lo;
    const center = (lo + hi) / 2;
    result.push({
      bin: `${lo.toFixed(0)}-${hi.toFixed(0)}`,
      binCenter: center,
      count: hist.counts[i],
    });
  }
  return result;
}

// =============================================================================
// Single histogram chart
// =============================================================================

function PropertyHistogram({ config, data }: { config: PropertyConfig; data: HistogramData }) {
  const barData = histogramToBarData(data);

  if (barData.length === 0) {
    return (
      <div className="flex items-center justify-center h-[200px] text-sm text-[var(--color-text-muted)]">
        No data
      </div>
    );
  }

  return (
    <div>
      <h4 className="text-xs font-semibold text-[var(--color-text-secondary)] mb-2 text-center">
        {config.label} {config.unit && `(${config.unit})`}
      </h4>
      <ResponsiveContainer width="100%" height={200}>
        <BarChart data={barData} margin={{ top: 5, right: 10, left: 0, bottom: 20 }}>
          <CartesianGrid strokeDasharray="3 3" stroke="var(--color-border)" opacity={0.5} />
          <XAxis
            dataKey="bin"
            tick={{ fontSize: 9, fill: 'var(--color-text-muted)' }}
            angle={-45}
            textAnchor="end"
            height={40}
          />
          <YAxis
            tick={{ fontSize: 10, fill: 'var(--color-text-muted)' }}
            label={{
              value: 'Count',
              angle: -90,
              position: 'insideLeft',
              style: { fontSize: 10, fill: 'var(--color-text-muted)' },
            }}
          />
          <Tooltip
            contentStyle={{
              background: 'var(--color-surface-elevated)',
              border: '1px solid var(--color-border)',
              borderRadius: '8px',
              fontSize: '12px',
            }}
            formatter={(value: number) => [value, 'Count']}
            labelFormatter={(label: string) => `Range: ${label}`}
          />
          {/* Dataset distribution bars */}
          <Bar
            dataKey="count"
            fill="var(--color-primary)"
            fillOpacity={0.6}
            radius={[2, 2, 0, 0]}
          />
          {/* Drug-space reference boundary lines (dashed) */}
          <ReferenceLine
            x={findClosestBin(barData, config.refLow)}
            stroke="#059669"
            strokeDasharray="4 4"
            strokeWidth={2}
            label={{
              value: String(config.refLow),
              position: 'top',
              style: { fontSize: 9, fill: '#059669' },
            }}
          />
          <ReferenceLine
            x={findClosestBin(barData, config.refHigh)}
            stroke="#059669"
            strokeDasharray="4 4"
            strokeWidth={2}
            label={{
              value: String(config.refHigh),
              position: 'top',
              style: { fontSize: 9, fill: '#059669' },
            }}
          />
        </BarChart>
      </ResponsiveContainer>
      <p className="text-[10px] text-[var(--color-text-muted)] text-center mt-1">
        Dashed lines: {config.refLabel}
      </p>
    </div>
  );
}

/** Find the bin label closest to a given numeric value. */
function findClosestBin(
  barData: Array<{ bin: string; binCenter: number }>,
  target: number,
): string {
  let closest = barData[0];
  let minDist = Math.abs(closest.binCenter - target);
  for (const item of barData) {
    const dist = Math.abs(item.binCenter - target);
    if (dist < minDist) {
      closest = item;
      minDist = dist;
    }
  }
  return closest?.bin ?? '';
}

// =============================================================================
// Main Component
// =============================================================================

/**
 * 3 property distribution histograms (MW, LogP, TPSA) with drug-space reference overlays.
 *
 * Per UI-SPEC:
 * - Responsive grid: lg:grid-cols-3, md:grid-cols-2, mobile: stack
 * - Blue bars for dataset distribution
 * - Dashed green reference lines at drug-space boundaries
 * - MW: Lipinski 150-500
 * - LogP: Lipinski -0.4 to 5
 * - TPSA: Veber 20-130
 */
export function PropertyDistOverlay({ distributions }: PropertyDistOverlayProps) {
  return (
    <div className="space-y-3">
      <h3 className="text-base font-semibold font-display text-[var(--color-text-primary)]">
        Property Distributions
      </h3>
      <div className="grid grid-cols-1 md:grid-cols-2 lg:grid-cols-3 gap-6">
        {PROPERTY_CONFIGS.map((config) => (
          <PropertyHistogram
            key={config.key}
            config={config}
            data={distributions[config.key]}
          />
        ))}
      </div>
    </div>
  );
}
