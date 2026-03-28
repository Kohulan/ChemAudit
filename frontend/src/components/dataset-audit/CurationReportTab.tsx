import { useState, useCallback } from 'react';
import { motion, AnimatePresence } from 'framer-motion';
import { ChevronDown, FileJson, FileSpreadsheet } from 'lucide-react';
import { ClayButton } from '../ui/ClayButton';
import { datasetApi } from '../../services/api';

// =============================================================================
// Types
// =============================================================================

interface CurationReportTabProps {
  /** Curation report data (null if not yet available). */
  report: Record<string, unknown> | null;
  /** Job ID for download URL construction. */
  jobId: string | null;
  /** Whether a curated CSV is available for download. */
  curatedCsvAvailable: boolean;
}

// =============================================================================
// Section config
// =============================================================================

interface ReportSection {
  key: string;
  title: string;
  /** Extract the collapsed summary line from the report data. */
  summary: (report: Record<string, unknown>) => string;
  /** Extract the expanded detail content. */
  detail: (report: Record<string, unknown>) => React.ReactNode;
}

const REPORT_SECTIONS: ReportSection[] = [
  {
    key: 'file_metadata',
    title: 'File Metadata',
    summary: (r) => {
      const meta = (r.file_metadata ?? r.metadata ?? {}) as Record<string, unknown>;
      return `${meta.filename ?? 'Unknown file'} - ${meta.format ?? 'Unknown format'} - ${meta.row_count ?? '?'} rows`;
    },
    detail: (r) => {
      const meta = (r.file_metadata ?? r.metadata ?? {}) as Record<string, unknown>;
      return (
        <dl className="grid grid-cols-2 gap-x-4 gap-y-2 text-xs">
          <dt className="text-[var(--color-text-secondary)]">Filename</dt>
          <dd className="font-mono">{String(meta.filename ?? '-')}</dd>
          <dt className="text-[var(--color-text-secondary)]">Format</dt>
          <dd className="font-mono">{String(meta.format ?? '-')}</dd>
          <dt className="text-[var(--color-text-secondary)]">Row Count</dt>
          <dd className="font-mono tabular-nums">{String(meta.row_count ?? '-')}</dd>
          <dt className="text-[var(--color-text-secondary)]">SHA-256</dt>
          <dd className="font-mono truncate" title={String(meta.sha256 ?? '')}>
            {String(meta.sha256 ?? '-')}
          </dd>
        </dl>
      );
    },
  },
  {
    key: 'health_audit',
    title: 'Health Audit Summary',
    summary: (r) => {
      const health = (r.health_audit ?? {}) as Record<string, unknown>;
      return `Overall score: ${health.overall_score ?? '?'}/100`;
    },
    detail: (r) => {
      const health = (r.health_audit ?? {}) as Record<string, unknown>;
      const subScores = (health.sub_scores ?? []) as Array<Record<string, unknown>>;
      return (
        <div className="space-y-2 text-xs">
          <p>
            <span className="text-[var(--color-text-secondary)]">Overall: </span>
            <span className="font-semibold tabular-nums">
              {String(health.overall_score ?? '-')}/100
            </span>
          </p>
          {subScores.length > 0 && (
            <ul className="space-y-1">
              {subScores.map((s) => (
                <li key={String(s.name)} className="flex justify-between">
                  <span className="text-[var(--color-text-secondary)] capitalize">
                    {String(s.name)}
                  </span>
                  <span className="font-mono tabular-nums">
                    {typeof s.score === 'number' ? (s.score * 100).toFixed(0) : '-'}%
                  </span>
                </li>
              ))}
            </ul>
          )}
        </div>
      );
    },
  },
  {
    key: 'contradictory_labels',
    title: 'Contradictory Labels',
    summary: (r) => {
      const count = (r.contradictory_labels_count ?? (r.contradictory_labels as unknown[])?.length) ?? 0;
      return `${count} found`;
    },
    detail: (r) => {
      const count = (r.contradictory_labels_count ?? (r.contradictory_labels as unknown[])?.length) ?? 0;
      return (
        <p className="text-xs">
          <span className="font-semibold tabular-nums">{String(count)}</span>{' '}
          contradictory label groups detected across the dataset.
        </p>
      );
    },
  },
  {
    key: 'standardization',
    title: 'Standardization Decisions',
    summary: (r) => {
      const count = (r.standardization_changes_count ?? (r.standardization_decisions as unknown[])?.length) ?? 0;
      return `${count} changed`;
    },
    detail: (r) => {
      const count = (r.standardization_changes_count ?? (r.standardization_decisions as unknown[])?.length) ?? 0;
      return (
        <p className="text-xs">
          <span className="font-semibold tabular-nums">{String(count)}</span>{' '}
          molecules were modified during standardization.
        </p>
      );
    },
  },
  {
    key: 'deduplication',
    title: 'Deduplication Summary',
    summary: (r) => {
      const count = (r.duplicate_count ?? (r.dedup_groups as unknown[])?.length) ?? 0;
      return `${count} duplicates`;
    },
    detail: (r) => {
      const count = (r.duplicate_count ?? (r.dedup_groups as unknown[])?.length) ?? 0;
      return (
        <p className="text-xs">
          <span className="font-semibold tabular-nums">{String(count)}</span>{' '}
          duplicate groups identified by InChIKey.
        </p>
      );
    },
  },
  {
    key: 'issues',
    title: 'Issues Flagged',
    summary: (r) => {
      const count = (r.total_issues ?? (r.issues as unknown[])?.length) ?? 0;
      return `${count} total`;
    },
    detail: (r) => {
      const count = (r.total_issues ?? (r.issues as unknown[])?.length) ?? 0;
      return (
        <p className="text-xs">
          <span className="font-semibold tabular-nums">{String(count)}</span>{' '}
          total issues flagged across all categories.
        </p>
      );
    },
  },
];

// =============================================================================
// Component
// =============================================================================

/**
 * Curation Report tab content for the Dataset Audit page.
 *
 * Per UI-SPEC D-14, D-15:
 * - Report preview with 6 collapsible accordion sections (Framer Motion)
 * - Two download buttons: JSON Report + Curated CSV
 * - Both disabled when report is null
 */
export function CurationReportTab({
  report,
  jobId,
  curatedCsvAvailable,
}: CurationReportTabProps) {
  const [expandedSections, setExpandedSections] = useState<Set<string>>(
    new Set(),
  );

  const toggleSection = useCallback((key: string) => {
    setExpandedSections((prev) => {
      const next = new Set(prev);
      if (next.has(key)) {
        next.delete(key);
      } else {
        next.add(key);
      }
      return next;
    });
  }, []);

  const handleDownloadReport = useCallback(() => {
    if (!jobId) return;
    const url = datasetApi.downloadReport(jobId);
    window.open(url, '_blank');
  }, [jobId]);

  const handleDownloadCsv = useCallback(() => {
    if (!jobId) return;
    const url = datasetApi.downloadCsv(jobId);
    window.open(url, '_blank');
  }, [jobId]);

  if (!report) {
    return (
      <div className="min-h-[200px] flex items-center justify-center">
        <p className="text-sm text-[var(--color-text-muted)]">
          Curation report will be available once the audit is complete.
        </p>
      </div>
    );
  }

  return (
    <div className="min-h-[200px] space-y-6">
      {/* Heading */}
      <h2 className="text-base font-semibold font-display text-[var(--color-text-primary)]">
        Curation Report
      </h2>

      {/* Report preview */}
      <div>
        <h3 className="text-sm font-medium text-[var(--color-text-secondary)] mb-3">
          Report Preview
        </h3>
        <div className="border border-[var(--color-border)] rounded-xl overflow-hidden divide-y divide-[var(--color-border)]">
          {REPORT_SECTIONS.map((section) => {
            const isOpen = expandedSections.has(section.key);
            return (
              <div key={section.key}>
                {/* Section header */}
                <button
                  className="w-full flex items-center justify-between px-4 py-3 text-left hover:bg-[var(--color-surface-sunken)]/50 transition-colors"
                  onClick={() => toggleSection(section.key)}
                  aria-expanded={isOpen}
                >
                  <div className="flex items-center gap-2">
                    <motion.span
                      animate={{ rotate: isOpen ? 180 : 0 }}
                      transition={{ duration: 0.2, ease: 'easeOut' }}
                    >
                      <ChevronDown className="w-4 h-4 text-[var(--color-text-muted)]" />
                    </motion.span>
                    <span className="text-sm font-medium text-[var(--color-text-primary)]">
                      {section.title}
                    </span>
                  </div>
                  <span className="text-xs text-[var(--color-text-muted)]">
                    {section.summary(report)}
                  </span>
                </button>

                {/* Section content */}
                <AnimatePresence initial={false}>
                  {isOpen && (
                    <motion.div
                      initial={{ height: 0, opacity: 0 }}
                      animate={{ height: 'auto', opacity: 1 }}
                      exit={{ height: 0, opacity: 0 }}
                      transition={{ duration: 0.2, ease: 'easeOut' }}
                      className="overflow-hidden"
                    >
                      <div className="px-4 pb-4 pt-1">
                        {section.detail(report)}
                      </div>
                    </motion.div>
                  )}
                </AnimatePresence>
              </div>
            );
          })}
        </div>
      </div>

      {/* Download buttons */}
      <div className="flex gap-4">
        <ClayButton
          variant="primary"
          onClick={handleDownloadReport}
          disabled={!report || !jobId}
          aria-label="Download full curation report as JSON file"
          leftIcon={<FileJson className="w-4 h-4" />}
        >
          Download JSON Report
        </ClayButton>

        <ClayButton
          variant="primary"
          onClick={handleDownloadCsv}
          disabled={!report || !jobId || !curatedCsvAvailable}
          aria-label="Download curated dataset with issue flags as CSV"
          leftIcon={<FileSpreadsheet className="w-4 h-4" />}
        >
          Download Curated CSV
        </ClayButton>
      </div>
    </div>
  );
}
