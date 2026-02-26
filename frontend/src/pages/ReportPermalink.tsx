import { useEffect, useState } from 'react';
import { useParams, useNavigate } from 'react-router-dom';
import { permalinksApi } from '../services/api';
import { MoleculeLoader } from '../components/ui/MoleculeLoader';

/**
 * Resolves a /report/:shortId permalink and navigates to /batch
 * with the resolved job_id passed via location state.
 */
export function ReportPermalinkPage() {
  const { shortId } = useParams<{ shortId: string }>();
  const navigate = useNavigate();
  const [error, setError] = useState<string | null>(null);

  useEffect(() => {
    if (!shortId) {
      setError('No report ID provided');
      return;
    }

    let cancelled = false;

    permalinksApi.resolvePermalink(shortId).then((data) => {
      if (cancelled) return;
      navigate('/batch', {
        replace: true,
        state: { permalinkJobId: data.job_id, permalinkSnapshot: data.snapshot_data },
      });
    }).catch((err) => {
      if (cancelled) return;
      const status = err?.response?.status;
      if (status === 410) {
        setError('This report link has expired.');
      } else if (status === 404) {
        setError('Report not found. The link may be invalid or the report may have been removed.');
      } else {
        setError('Failed to load report. Please try again later.');
      }
    });

    return () => { cancelled = true; };
  }, [shortId, navigate]);

  if (error) {
    return (
      <div className="flex flex-col items-center justify-center h-64 gap-4">
        <div className="w-14 h-14 rounded-2xl flex items-center justify-center bg-red-500/10">
          <svg className="w-7 h-7 text-red-500" viewBox="0 0 24 24" fill="none" stroke="currentColor" strokeWidth="2">
            <circle cx="12" cy="12" r="10" />
            <path d="M15 9l-6 6M9 9l6 6" />
          </svg>
        </div>
        <p className="text-text-secondary text-center max-w-md">{error}</p>
        <button
          onClick={() => navigate('/')}
          className="text-sm text-accent-primary hover:underline"
        >
          Go to home page
        </button>
      </div>
    );
  }

  return (
    <div className="flex items-center justify-center h-64">
      <MoleculeLoader size="lg" text="Loading shared report..." />
    </div>
  );
}
