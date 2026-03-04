import { useState } from 'react';
import { safeHref } from '../../lib/sanitize';
import type { ConsistencyResult, PropertyComparison } from '../../types/integrations';

interface DatabaseComparisonPanelProps {
  result: ConsistencyResult;
}

const VERDICT_STYLES: Record<string, { bg: string; text: string; label: string }> = {
  consistent: { bg: 'bg-green-100', text: 'text-green-800', label: 'Consistent' },
  minor_differences: { bg: 'bg-yellow-100', text: 'text-yellow-800', label: 'Minor Differences' },
  major_discrepancies: { bg: 'bg-red-100', text: 'text-red-800', label: 'Major Discrepancies' },
  no_data: { bg: 'bg-gray-100', text: 'text-gray-600', label: 'No Data' },
};

const STATUS_COLORS: Record<string, string> = {
  match: 'bg-green-50 text-green-800',
  minor_difference: 'bg-yellow-50 text-yellow-800',
  mismatch: 'bg-red-50 text-red-800',
  missing: 'bg-gray-50 text-gray-500',
};

const PROP_LABELS: Record<string, string> = {
  canonical_smiles: 'SMILES',
  molecular_formula: 'Formula',
  molecular_weight: 'MW (Da)',
  inchikey: 'InChIKey',
  inchi: 'InChI',
};

function ComparisonRow({ comparison }: { comparison: PropertyComparison }) {
  const [expanded, setExpanded] = useState(false);

  return (
    <tr className="border-t border-gray-100">
      <td className="py-2 px-3 text-sm font-medium text-gray-700">
        <button
          onClick={() => setExpanded(!expanded)}
          className="flex items-center gap-1 hover:text-indigo-600"
        >
          <svg
            className={`w-3 h-3 transition-transform ${expanded ? 'rotate-90' : ''}`}
            fill="none" stroke="currentColor" viewBox="0 0 24 24"
          >
            <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={2} d="M9 5l7 7-7 7" />
          </svg>
          {PROP_LABELS[comparison.property_name] || comparison.property_name}
        </button>
        {expanded && comparison.detail && (
          <div className="text-xs text-gray-500 mt-1 ml-4">{comparison.detail}</div>
        )}
      </td>
      {Object.entries(comparison.values).map(([db, value]) => (
        <td key={db} className={`py-2 px-3 text-xs font-mono truncate max-w-[200px] ${STATUS_COLORS[comparison.status] || ''}`}>
          {value ?? <span className="text-gray-400 italic">—</span>}
        </td>
      ))}
      <td className="py-2 px-3">
        <span className={`px-2 py-0.5 text-xs rounded-full ${STATUS_COLORS[comparison.status] || 'bg-gray-100 text-gray-600'}`}>
          {comparison.status}
        </span>
      </td>
    </tr>
  );
}

export function DatabaseComparisonPanel({ result }: DatabaseComparisonPanelProps) {
  const verdict = VERDICT_STYLES[result.overall_verdict] || VERDICT_STYLES.no_data;
  const foundEntries = result.entries.filter(e => e.found);
  const databases = foundEntries.map(e => e.database);

  return (
    <div className="bg-white border border-gray-200 rounded-lg p-4 mt-4">
      {/* Verdict badge */}
      <div className="flex items-center justify-between mb-4">
        <h4 className="text-sm font-semibold text-gray-900">Cross-Database Comparison</h4>
        <span className={`px-3 py-1 text-xs font-medium rounded-full ${verdict.bg} ${verdict.text}`}>
          {verdict.label}
        </span>
      </div>

      <p className="text-sm text-gray-600 mb-4">{result.summary}</p>

      {/* Database presence */}
      <div className="flex gap-2 mb-4">
        {result.entries.map((entry) => (
          <div
            key={entry.database}
            className={`px-3 py-1.5 rounded-lg text-xs font-medium ${
              entry.found
                ? 'bg-blue-50 text-blue-800 border border-blue-200'
                : 'bg-gray-50 text-gray-400 border border-gray-200'
            }`}
          >
            {entry.database}
            {entry.found ? ' ✓' : ' —'}
            {entry.found && entry.url && (
              <a
                href={safeHref(entry.url)}
                target="_blank"
                rel="noopener noreferrer"
                className="ml-1 text-blue-500 hover:text-blue-700"
              >
                ↗
              </a>
            )}
          </div>
        ))}
      </div>

      {/* Comparison table */}
      {result.comparisons.length > 0 && databases.length >= 2 && (
        <div className="overflow-x-auto">
          <table className="w-full text-left">
            <thead>
              <tr className="border-b border-gray-200">
                <th className="py-2 px-3 text-xs font-medium text-gray-500 uppercase">Property</th>
                {databases.map((db) => (
                  <th key={db} className="py-2 px-3 text-xs font-medium text-gray-500 uppercase">{db}</th>
                ))}
                <th className="py-2 px-3 text-xs font-medium text-gray-500 uppercase">Status</th>
              </tr>
            </thead>
            <tbody>
              {result.comparisons.map((comp) => (
                <ComparisonRow key={comp.property_name} comparison={comp} />
              ))}
            </tbody>
          </table>
        </div>
      )}
    </div>
  );
}
