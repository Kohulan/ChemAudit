import { motion } from 'framer-motion';
import {
  BarChart,
  Bar,
  XAxis,
  YAxis,
  Tooltip,
  ResponsiveContainer,
} from 'recharts';
import { ClayCard } from '../ui/ClayCard';
import type { QSARBatchSummary as QSARBatchSummaryType } from '../../types/qsar_ready';

interface QSARBatchSummaryProps {
  summary: QSARBatchSummaryType;
}

// ─── Metric card data ───────────────────────────────────────────────────────

interface MetricCard {
  label: string;
  value: number;
  color: string;
}

function getMetricCards(summary: QSARBatchSummaryType): MetricCard[] {
  return [
    { label: 'OK', value: summary.ok, color: '#16a34a' },
    { label: 'Rejected', value: summary.rejected, color: '#dc2626' },
    { label: 'Duplicate', value: summary.duplicate, color: '#d97706' },
    { label: 'Error', value: summary.error, color: '#7c3aed' },
  ];
}

// ─── Component ───────────────────────────────────────────────────────────────

/**
 * Batch summary panel for QSAR-Ready results.
 *
 * Per UI-SPEC batch result layout:
 * - 4 metric cards in grid grid-cols-2 md:grid-cols-4 gap-4
 * - Each card: count in text-2xl font-semibold + label below
 * - Steps distribution bar chart: Recharts BarChart, full width, 80px height
 * - Reveal animation: fade-in 300ms ease-out
 */
export function QSARBatchSummary({ summary }: QSARBatchSummaryProps) {
  const cards = getMetricCards(summary);

  // Build chart data from steps_applied_counts
  const chartData = Object.entries(summary.steps_applied_counts).map(
    ([step, count]) => ({ step, count }),
  );

  return (
    <motion.div
      initial={{ opacity: 0 }}
      animate={{ opacity: 1 }}
      transition={{ duration: 0.3, ease: 'easeOut' }}
      className="space-y-6"
    >
      {/* ── 4 Metric cards ── */}
      <div className="grid grid-cols-2 md:grid-cols-4 gap-4">
        {cards.map((card) => (
          <ClayCard key={card.label} variant="default" size="sm">
            <div className="flex flex-col items-center py-2">
              <span
                className="text-2xl font-semibold font-display"
                style={{ color: card.color }}
              >
                {card.value.toLocaleString()}
              </span>
              <span className="text-xs text-[var(--color-text-secondary)] mt-1">
                {card.label}
              </span>
            </div>
          </ClayCard>
        ))}
      </div>

      {/* ── Steps distribution bar chart ── */}
      {chartData.length > 0 && (
        <div>
          <p className="text-xs font-medium text-[var(--color-text-secondary)] mb-2">
            Steps Applied
          </p>
          <ResponsiveContainer width="100%" height={80}>
            <BarChart
              data={chartData}
              margin={{ top: 0, right: 0, left: -20, bottom: 0 }}
            >
              <XAxis
                dataKey="step"
                tick={{ fontSize: 10, fill: 'var(--color-text-muted)' }}
                axisLine={false}
                tickLine={false}
              />
              <YAxis
                tick={{ fontSize: 10, fill: 'var(--color-text-muted)' }}
                axisLine={false}
                tickLine={false}
              />
              <Tooltip
                contentStyle={{
                  background: 'var(--color-surface-elevated)',
                  border: '1px solid var(--color-border)',
                  borderRadius: '6px',
                  fontSize: '12px',
                }}
                cursor={{ fill: 'var(--color-surface-sunken)' }}
              />
              <Bar
                dataKey="count"
                fill="var(--color-primary)"
                radius={[2, 2, 0, 0]}
              />
            </BarChart>
          </ResponsiveContainer>
        </div>
      )}
    </motion.div>
  );
}
