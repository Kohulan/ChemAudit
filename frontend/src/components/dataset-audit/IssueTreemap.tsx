import React, { useMemo, useCallback } from 'react';
import { Treemap, ResponsiveContainer, Tooltip } from 'recharts';
import type { IssueEntry } from '../../types/dataset_intelligence';
import { InfoTooltip } from '../ui/Tooltip';

// =============================================================================
// Types
// =============================================================================

interface IssueTreemapProps {
  /** Array of issue entries from health audit. */
  issues: IssueEntry[];
  /** Callback when a treemap cell is clicked. */
  onCellClick: (issueType: string) => void;
  /** Currently selected issue type (for ring highlight). */
  selectedIssueType: string | null;
}

// =============================================================================
// Category colors per UI-SPEC
// =============================================================================

const CATEGORY_COLORS: Record<string, string> = {
  'Parse failures': '#c41e3a',
  'Undefined stereo': '#d97706',
  'Duplicates': '#ea580c',
  'Alert hits': '#059669',
  'Std. inconsistency': '#2563eb',
};

const DEFAULT_COLOR = '#6b7280';

function getCategoryColor(issueType: string): string {
  return CATEGORY_COLORS[issueType] ?? DEFAULT_COLOR;
}

// =============================================================================
// Custom content renderer
// =============================================================================

function CustomContent(props: Record<string, unknown>) {
  const {
    x,
    y,
    width,
    height,
    name,
    value,
    color,
    selectedIssueType,
    onCellClick,
  } = props as {
    x: number;
    y: number;
    width: number;
    height: number;
    name: string;
    value: number;
    color: string;
    selectedIssueType: string | null;
    onCellClick: (name: string) => void;
  };

  if (width < 4 || height < 4) return null;

  const isSelected = selectedIssueType === name;
  const isDimmed = selectedIssueType && !isSelected;

  return (
    <g
      role="button"
      tabIndex={0}
      aria-label={`${name}: ${value} molecules affected`}
      onClick={(e) => {
        e.stopPropagation();
        onCellClick(name);
      }}
      onKeyDown={(e) => {
        if (e.key === 'Enter' || e.key === ' ') {
          e.preventDefault();
          onCellClick(name);
        }
      }}
      style={{ cursor: 'pointer', outline: 'none' }}
    >
      <rect
        x={x}
        y={y}
        width={Math.max(width, 44)} // Minimum 44px touch target
        height={Math.max(height, 44)}
        fill={color || DEFAULT_COLOR}
        fillOpacity={isDimmed ? 0.4 : 0.85}
        stroke={isSelected ? '#fff' : 'var(--color-surface)'}
        strokeWidth={isSelected ? 3 : 2}
        rx={4}
        style={{ transition: 'opacity 0.15s ease' }}
      />
      {/* Category name */}
      {width > 60 && height > 30 && (
        <text
          x={x + width / 2}
          y={y + height / 2 - (height > 50 ? 6 : 0)}
          textAnchor="middle"
          dominantBaseline="middle"
          fill="#fff"
          fontSize={11}
          fontWeight={600}
          style={{ pointerEvents: 'none' }}
        >
          {name.length > Math.floor(width / 7)
            ? name.slice(0, Math.floor(width / 7)) + '...'
            : name}
        </text>
      )}
      {/* Count */}
      {width > 60 && height > 50 && (
        <text
          x={x + width / 2}
          y={y + height / 2 + 12}
          textAnchor="middle"
          dominantBaseline="middle"
          fill="rgba(255,255,255,0.8)"
          fontSize={10}
          style={{ pointerEvents: 'none' }}
        >
          {value}
        </text>
      )}
    </g>
  );
}

// =============================================================================
// Component
// =============================================================================

/**
 * Recharts Treemap showing issue categories with clickable cells.
 *
 * Follows the ValidationTreemap.tsx pattern from Phase 3 VIZ-04.
 *
 * Per UI-SPEC:
 * - 6 category colors for specific issue types
 * - Clickable cells with role="button" + keyboard Enter/Space
 * - Selected cell gets white ring highlight
 * - Dimmed non-selected cells on selection
 * - Minimum 44px touch target
 */
export const IssueTreemap = React.memo(function IssueTreemap({
  issues,
  onCellClick,
  selectedIssueType,
}: IssueTreemapProps) {
  // Group issues by issue_type and count
  const treemapData = useMemo(() => {
    const counts: Record<string, number> = {};
    for (const issue of issues) {
      counts[issue.issue_type] = (counts[issue.issue_type] ?? 0) + 1;
    }
    return Object.entries(counts)
      .sort(([, a], [, b]) => b - a)
      .map(([name, value]) => ({
        name,
        value,
        color: getCategoryColor(name),
      }));
  }, [issues]);

  const handleCellClick = useCallback(
    (name: string) => {
      // Toggle: clicking active filter clears it
      onCellClick(selectedIssueType === name ? '' : name);
    },
    [onCellClick, selectedIssueType],
  );

  if (treemapData.length === 0) {
    return (
      <div className="space-y-3">
        <h3 className="text-base font-semibold font-display text-[var(--color-text-primary)]">
          Issue Breakdown
        </h3>
        <div className="flex items-center justify-center h-[240px] text-sm text-[var(--color-text-muted)]">
          No issues detected in this dataset.
        </div>
      </div>
    );
  }

  return (
    <div className="space-y-3">
      <div className="flex items-center gap-1.5">
        <h3 className="text-base font-semibold font-display text-[var(--color-text-primary)]">
          Issue Breakdown
        </h3>
        <InfoTooltip
          title="Issue Breakdown"
          content={
            <div className="text-xs space-y-1">
              <p>Treemap visualization showing the relative prevalence of each issue category in your dataset.</p>
              <ul className="mt-1 text-white/70 space-y-0.5">
                <li><strong>Parse failures</strong> &mdash; SMILES that RDKit cannot interpret</li>
                <li><strong>Undefined stereo</strong> &mdash; Chiral centers or double-bond geometry not specified</li>
                <li><strong>Duplicates</strong> &mdash; Multiple rows with the same InChIKey</li>
                <li><strong>Alert hits</strong> &mdash; Molecules flagging structural alerts (PAINS, etc.)</li>
                <li><strong>Std. inconsistency</strong> &mdash; Pipelines disagree on standardized form</li>
              </ul>
              <p className="mt-1 text-white/60">Click a cell to drill down into affected molecules. Larger cells represent more prevalent issues.</p>
            </div>
          }
          size="small"
        />
      </div>
      <ResponsiveContainer width="100%" height={280}>
        <Treemap
          data={treemapData}
          dataKey="value"
          aspectRatio={4 / 3}
          content={
            <CustomContent
              selectedIssueType={selectedIssueType}
              onCellClick={handleCellClick}
            />
          }
        >
          <Tooltip
            content={({ active, payload }) => {
              if (!active || !payload?.[0]) return null;
              const entry = payload[0].payload as { name: string; value: number };
              return (
                <div className="rounded-lg bg-[var(--color-surface-elevated)] border border-[var(--color-border)] shadow-lg p-3 text-xs max-w-[280px]">
                  <p className="font-semibold text-[var(--color-text-primary)] mb-1">
                    {entry.name}
                  </p>
                  <p className="text-[var(--color-text-muted)]">
                    {entry.value} molecule{entry.value !== 1 ? 's' : ''} affected
                  </p>
                  <p className="text-[var(--color-primary)] mt-1.5 font-medium">
                    Click to view affected molecules
                  </p>
                </div>
              );
            }}
          />
        </Treemap>
      </ResponsiveContainer>
    </div>
  );
});
