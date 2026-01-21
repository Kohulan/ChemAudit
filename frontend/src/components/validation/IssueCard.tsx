import type { CheckResult, Severity } from '../../types/validation';

interface IssueCardProps {
  issue: CheckResult;
  onAtomHover?: (atoms: number[]) => void;
  className?: string;
}

export function IssueCard({ issue, onAtomHover, className = '' }: IssueCardProps) {
  const getSeverityStyles = (severity: Severity) => {
    switch (severity) {
      case 'critical':
        return {
          bg: 'bg-red-50',
          border: 'border-red-200',
          badge: 'bg-red-100 text-red-800',
          icon: '🚫',
        };
      case 'error':
        return {
          bg: 'bg-orange-50',
          border: 'border-orange-200',
          badge: 'bg-orange-100 text-orange-800',
          icon: '⚠️',
        };
      case 'warning':
        return {
          bg: 'bg-yellow-50',
          border: 'border-yellow-200',
          badge: 'bg-yellow-100 text-yellow-800',
          icon: '⚡',
        };
      case 'info':
        return {
          bg: 'bg-blue-50',
          border: 'border-blue-200',
          badge: 'bg-blue-100 text-blue-800',
          icon: 'ℹ️',
        };
      default:
        return {
          bg: 'bg-gray-50',
          border: 'border-gray-200',
          badge: 'bg-gray-100 text-gray-800',
          icon: '✓',
        };
    }
  };

  const styles = getSeverityStyles(issue.severity);

  const formatCheckName = (name: string) => {
    return name
      .split('_')
      .map((word) => word.charAt(0).toUpperCase() + word.slice(1))
      .join(' ');
  };

  return (
    <div
      className={`${styles.bg} border ${styles.border} rounded-lg p-4 ${className}`}
      onMouseEnter={() => onAtomHover?.(issue.affected_atoms)}
      onMouseLeave={() => onAtomHover?.([])}
    >
      <div className="flex items-start gap-3">
        <span className="text-2xl mt-0.5">{styles.icon}</span>
        <div className="flex-1 min-w-0">
          <div className="flex items-center gap-2 mb-2">
            <h4 className="font-medium text-gray-900">
              {formatCheckName(issue.check_name)}
            </h4>
            <span
              className={`px-2 py-0.5 text-xs font-medium rounded-full ${styles.badge}`}
            >
              {issue.severity.toUpperCase()}
            </span>
          </div>
          <p className="text-sm text-gray-700">{issue.message}</p>
          {issue.affected_atoms.length > 0 && (
            <div className="mt-2 text-xs text-gray-500">
              Affected atoms: {issue.affected_atoms.join(', ')}
            </div>
          )}
          {Object.keys(issue.details).length > 0 && (
            <div className="mt-2 text-xs text-gray-500">
              <details className="cursor-pointer">
                <summary className="hover:text-gray-700">Additional details</summary>
                <pre className="mt-1 p-2 bg-white rounded border border-gray-200 overflow-x-auto">
                  {JSON.stringify(issue.details, null, 2)}
                </pre>
              </details>
            </div>
          )}
        </div>
      </div>
    </div>
  );
}
