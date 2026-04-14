import { Badge } from '../ui/Badge';
import type { FileIssue } from '../../types/diagnostics';

interface FileIssueListProps {
  issues: FileIssue[];
}

/**
 * Map backend issue_type to human-readable group header.
 *
 * Per UI-SPEC copywriting contract (D-10):
 * - "missing_m_end" -> "Missing M END block"
 * - "malformed_count_line" -> "Malformed atom count line"
 * - "encoding_fallback" | "encoding" -> "Encoding issue"
 * - "missing_smiles_column" -> "Missing SMILES column"
 * - default: capitalize issue_type replacing underscores with spaces
 */
function formatIssueType(issueType: string): string {
  switch (issueType) {
    case 'missing_m_end':
      return 'Missing M END block';
    case 'malformed_count_line':
      return 'Malformed atom count line';
    case 'encoding_fallback':
    case 'encoding':
      return 'Encoding issue';
    case 'missing_smiles_column':
      return 'Missing SMILES column';
    case 'duplicate_columns':
      return 'Duplicate column headers';
    default:
      // Capitalize words, replace underscores with spaces
      return issueType
        .split('_')
        .map((w) => w.charAt(0).toUpperCase() + w.slice(1))
        .join(' ');
  }
}

/** Get Badge variant for a severity level. */
function getSeverityVariant(severity: string): 'error' | 'warning' | 'default' {
  switch (severity.toLowerCase()) {
    case 'error':
      return 'error';
    case 'warning':
      return 'warning';
    default:
      return 'default';
  }
}

/** Get Badge label for a severity level (always text, not icon-only per accessibility). */
function getSeverityLabel(severity: string): string {
  switch (severity.toLowerCase()) {
    case 'error':
      return 'Error';
    case 'warning':
      return 'Warning';
    case 'info':
      return 'Info';
    default:
      return severity.charAt(0).toUpperCase() + severity.slice(1);
  }
}

/**
 * Scrollable grouped file issue list (DIAG-05, D-10).
 *
 * Groups issues by issue_type, shows:
 * - Group header with human-readable type label
 * - Each issue row: line number, severity badge (text, not icon-only), description
 * - Min-height 44px per row (WCAG 2.5.5 touch target)
 * - Max height 256px (max-h-64) with overflow scroll
 *
 * Renders null when issues.length === 0 (parent handles empty state).
 */
export function FileIssueList({ issues }: FileIssueListProps) {
  if (issues.length === 0) return null;

  // Group issues by issue_type
  const groups = issues.reduce<Record<string, FileIssue[]>>((acc, issue) => {
    const key = issue.issue_type;
    if (!acc[key]) acc[key] = [];
    acc[key].push(issue);
    return acc;
  }, {});

  return (
    <div className="max-h-64 overflow-y-auto rounded-lg border border-[var(--color-border)]">
      {Object.entries(groups).map(([issueType, groupIssues]) => (
        <div key={issueType}>
          {/* Group header */}
          <div className="px-4 py-2 bg-[var(--color-surface-sunken)] border-b border-[var(--color-border)]">
            <span className="text-xs font-semibold text-[var(--color-text-muted)] uppercase tracking-wide">
              {formatIssueType(issueType)}
            </span>
            <span className="ml-2 text-xs text-[var(--color-text-muted)]">
              ({groupIssues.length})
            </span>
          </div>

          {/* Issue rows */}
          {groupIssues.map((issue, idx) => (
            <div
              key={idx}
              className="flex items-start gap-3 px-4 py-2 border-b border-[var(--color-border)] last:border-0 min-h-[44px]"
            >
              {/* Line number */}
              <span className="text-xs font-mono text-[var(--color-text-muted)] shrink-0 mt-0.5 w-16">
                Line {issue.line}
              </span>

              {/* Severity badge — always text label, never icon-only (accessibility) */}
              <div className="shrink-0 mt-0.5">
                <Badge variant={getSeverityVariant(issue.severity)} size="sm">
                  {getSeverityLabel(issue.severity)}
                </Badge>
              </div>

              {/* Description */}
              <span className="text-sm text-[var(--color-text-secondary)] flex-1 min-w-0">
                {issue.description}
              </span>
            </div>
          ))}
        </div>
      ))}
    </div>
  );
}
