import { useState } from 'react';
import type { ValidationResponse } from '../../types/validation';
import { ScoreGauge } from './ScoreGauge';
import { IssueCard } from './IssueCard';

interface ValidationResultsProps {
  result: ValidationResponse;
  onHighlightAtoms?: (atoms: number[]) => void;
  className?: string;
}

export function ValidationResults({
  result,
  onHighlightAtoms,
  className = '',
}: ValidationResultsProps) {
  const [showAllChecks, setShowAllChecks] = useState(false);

  const { molecule_info, overall_score, issues, all_checks, execution_time_ms } = result;

  return (
    <div className={`space-y-6 ${className}`}>
      {/* Score Summary */}
      <div className="bg-white rounded-lg shadow-md p-6">
        <h3 className="text-lg font-semibold text-gray-900 mb-4 text-center">
          Validation Score
        </h3>
        <ScoreGauge score={overall_score} size={140} className="mx-auto" />
        <div className="mt-4 text-center text-sm text-gray-500">
          Completed in {execution_time_ms.toFixed(0)}ms
        </div>
      </div>

      {/* Molecule Information */}
      <div className="bg-white rounded-lg shadow-md p-6">
        <h3 className="text-lg font-semibold text-gray-900 mb-4">
          Molecule Information
        </h3>
        <div className="space-y-2 text-sm">
          {molecule_info.canonical_smiles && (
            <div className="flex">
              <span className="font-medium text-gray-700 w-32">SMILES:</span>
              <code className="flex-1 text-gray-600 font-mono text-xs break-all">
                {molecule_info.canonical_smiles}
              </code>
            </div>
          )}
          {molecule_info.molecular_formula && (
            <div className="flex">
              <span className="font-medium text-gray-700 w-32">Formula:</span>
              <span className="text-gray-600">{molecule_info.molecular_formula}</span>
            </div>
          )}
          {molecule_info.molecular_weight && (
            <div className="flex">
              <span className="font-medium text-gray-700 w-32">Mol. Weight:</span>
              <span className="text-gray-600">
                {molecule_info.molecular_weight.toFixed(2)} g/mol
              </span>
            </div>
          )}
          {molecule_info.inchikey && (
            <div className="flex">
              <span className="font-medium text-gray-700 w-32">InChIKey:</span>
              <code className="flex-1 text-gray-600 font-mono text-xs break-all">
                {molecule_info.inchikey}
              </code>
            </div>
          )}
          {molecule_info.num_atoms !== null && (
            <div className="flex">
              <span className="font-medium text-gray-700 w-32">Atoms:</span>
              <span className="text-gray-600">{molecule_info.num_atoms}</span>
            </div>
          )}
        </div>
      </div>

      {/* Issues */}
      {issues.length > 0 ? (
        <div className="bg-white rounded-lg shadow-md p-6">
          <h3 className="text-lg font-semibold text-gray-900 mb-4">
            Issues Found ({issues.length})
          </h3>
          <div className="space-y-3">
            {issues.map((issue, index) => (
              <IssueCard
                key={`${issue.check_name}-${index}`}
                issue={issue}
                onAtomHover={onHighlightAtoms}
              />
            ))}
          </div>
        </div>
      ) : (
        <div className="bg-green-50 border border-green-200 rounded-lg p-6 text-center">
          <div className="text-4xl mb-2">✓</div>
          <h3 className="text-lg font-semibold text-green-900 mb-1">
            No Issues Found
          </h3>
          <p className="text-sm text-green-700">
            All validation checks passed successfully
          </p>
        </div>
      )}

      {/* All Checks (collapsible) */}
      <div className="bg-white rounded-lg shadow-md p-6">
        <button
          onClick={() => setShowAllChecks(!showAllChecks)}
          className="w-full flex items-center justify-between text-left"
        >
          <h3 className="text-lg font-semibold text-gray-900">
            All Checks ({all_checks.length})
          </h3>
          <svg
            className={`w-5 h-5 text-gray-500 transition-transform ${
              showAllChecks ? 'transform rotate-180' : ''
            }`}
            fill="none"
            stroke="currentColor"
            viewBox="0 0 24 24"
          >
            <path
              strokeLinecap="round"
              strokeLinejoin="round"
              strokeWidth={2}
              d="M19 9l-7 7-7-7"
            />
          </svg>
        </button>

        {showAllChecks && (
          <div className="mt-4 space-y-2">
            {all_checks.map((check, index) => (
              <div
                key={`${check.check_name}-${index}`}
                className="flex items-center justify-between py-2 px-3 bg-gray-50 rounded"
              >
                <div className="flex items-center gap-2">
                  <span>{check.passed ? '✓' : '✗'}</span>
                  <span className="text-sm font-medium text-gray-900">
                    {check.check_name.replace(/_/g, ' ')}
                  </span>
                </div>
                <span
                  className={`text-xs px-2 py-1 rounded ${
                    check.passed
                      ? 'bg-green-100 text-green-800'
                      : check.severity === 'critical'
                      ? 'bg-red-100 text-red-800'
                      : check.severity === 'error'
                      ? 'bg-orange-100 text-orange-800'
                      : check.severity === 'warning'
                      ? 'bg-yellow-100 text-yellow-800'
                      : 'bg-blue-100 text-blue-800'
                  }`}
                >
                  {check.passed ? 'PASS' : check.severity.toUpperCase()}
                </span>
              </div>
            ))}
          </div>
        )}
      </div>
    </div>
  );
}
