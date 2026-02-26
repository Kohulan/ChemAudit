import { useState } from 'react';
import { ChevronDown, ChevronUp } from 'lucide-react';
import type { AlertResult, AlertSeverity } from '../../types/alerts';

interface AlertCardProps {
  alert: AlertResult;
  onAtomHover?: (atoms: number[]) => void;
  className?: string;
}

const CATEGORY_STYLES: Record<string, { bg: string; text: string }> = {
  'Reactive Group': { bg: 'bg-red-100', text: 'text-red-700' },
  'Toxicophore': { bg: 'bg-rose-100', text: 'text-rose-700' },
  'Metabolic Liability': { bg: 'bg-amber-100', text: 'text-amber-700' },
  'Assay Interference': { bg: 'bg-purple-100', text: 'text-purple-700' },
  'Physicochemical': { bg: 'bg-slate-100', text: 'text-slate-700' },
  'Unwanted Functionality': { bg: 'bg-gray-100', text: 'text-gray-600' },
};

/**
 * Known PAINS patterns that appear in FDA-approved drugs.
 * Provides educational context to help users understand alerts are warnings.
 */
const APPROVED_DRUG_EXAMPLES: Record<string, string[]> = {
  rhodanine: ['Methotrexate', 'Epalrestat'],
  catechol: ['Dopamine', 'Epinephrine', 'Levodopa'],
  quinone: ['Doxorubicin', 'Vitamin K'],
  michael_acceptor: ['Ibrutinib', 'Afatinib'],
  azo: ['Sulfasalazine', 'Phenazopyridine'],
  thiourea: ['Methimazole', 'Thiouracil'],
};

function getApprovedDrugNote(patternName: string): string | null {
  const patternLower = patternName.toLowerCase();
  for (const [key, drugs] of Object.entries(APPROVED_DRUG_EXAMPLES)) {
    if (patternLower.includes(key)) {
      return `Found in approved drugs: ${drugs.slice(0, 2).join(', ')}`;
    }
  }
  return null;
}

export function AlertCard({ alert, onAtomHover, className = '' }: AlertCardProps) {
  const [expanded, setExpanded] = useState(false);

  const getSeverityStyles = (severity: AlertSeverity) => {
    switch (severity) {
      case 'critical':
        return {
          bg: 'bg-red-50',
          border: 'border-red-200',
          badge: 'bg-red-100 text-red-800',
          icon: '!!!',
        };
      case 'warning':
        return {
          bg: 'bg-amber-50',
          border: 'border-amber-200',
          badge: 'bg-amber-100 text-amber-800',
          icon: '!',
        };
      case 'info':
        return {
          bg: 'bg-blue-50',
          border: 'border-blue-200',
          badge: 'bg-blue-100 text-blue-800',
          icon: 'i',
        };
      default:
        return {
          bg: 'bg-gray-50',
          border: 'border-gray-200',
          badge: 'bg-gray-100 text-gray-800',
          icon: '?',
        };
    }
  };

  const styles = getSeverityStyles(alert.severity);
  const categoryName = alert.category ?? null;
  const categoryStyle = categoryName
    ? (CATEGORY_STYLES[categoryName] || CATEGORY_STYLES['Unwanted Functionality'])
    : null;

  const formatPatternName = (name: string) => {
    const cleaned = name.replace(/\(\d+\)$/, '').trim();
    return cleaned
      .split('_')
      .map((word) => word.charAt(0).toUpperCase() + word.slice(1))
      .join(' ');
  };

  const approvedDrugNote = getApprovedDrugNote(alert.pattern_name);

  return (
    <div
      className={`${styles.bg} border ${styles.border} rounded-lg p-4 transition-all hover:shadow-md ${className}`}
      onMouseEnter={() => onAtomHover?.(alert.matched_atoms)}
      onMouseLeave={() => onAtomHover?.([])}
    >
      <div className="flex items-start gap-3">
        <div
          className={`flex items-center justify-center w-8 h-8 rounded-full ${styles.badge} font-bold text-sm flex-shrink-0`}
        >
          {styles.icon}
        </div>
        <div className="flex-1 min-w-0">
          {/* Row 1: Pattern name + badges */}
          <div className="flex items-center gap-2 mb-1.5 flex-wrap">
            <h4 className="font-medium text-gray-900">
              {formatPatternName(alert.pattern_name)}
            </h4>
            <span className={`px-2 py-0.5 text-xs font-medium rounded-full ${styles.badge}`}>
              {alert.severity.toUpperCase()}
            </span>
            {categoryStyle && categoryName && (
              <span
                className={`px-2 py-0.5 text-xs font-medium rounded-full ${categoryStyle.bg} ${categoryStyle.text}`}
              >
                {categoryName}
              </span>
            )}
          </div>

          {/* Row 2: Human-readable catalog name */}
          <p className="text-sm text-gray-600 mb-1.5">
            {alert.catalog_description || alert.catalog_source}
          </p>

          {/* Row 3: Scope — what this filter screens for */}
          {alert.scope && <p className="text-sm text-gray-500 italic mb-2">{alert.scope}</p>}

          {/* Matched atoms */}
          {alert.matched_atoms.length > 0 && (
            <div className="text-xs text-gray-500 mb-2">
              <span className="font-medium">Matched atoms:</span>{' '}
              {alert.matched_atoms.join(', ')}
              <span className="ml-2 text-amber-600">(hover to highlight)</span>
            </div>
          )}

          {/* Approved drug note */}
          {approvedDrugNote && (
            <div className="text-xs text-amber-700 bg-yellow-50 rounded px-2 py-1 inline-block mb-2">
              {approvedDrugNote}
            </div>
          )}

          {/* Expandable reference section */}
          {alert.reference && (
            <button
              onClick={(e) => {
                e.stopPropagation();
                setExpanded(!expanded);
              }}
              className="flex items-center gap-1 text-xs text-gray-400 hover:text-gray-600 transition-colors"
            >
              {expanded ? (
                <ChevronUp className="w-3 h-3" />
              ) : (
                <ChevronDown className="w-3 h-3" />
              )}
              {expanded ? 'Hide reference' : 'Show reference'}
            </button>
          )}
          {expanded && alert.reference && (
            <div className="mt-2 text-xs text-gray-500 bg-white/50 rounded p-2 border border-gray-100">
              <span className="font-medium">Reference:</span> {alert.reference}
            </div>
          )}
        </div>
      </div>
    </div>
  );
}
