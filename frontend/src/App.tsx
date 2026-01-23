import { lazy, Suspense } from 'react';
import { BrowserRouter as Router, Routes, Route, NavLink } from 'react-router-dom';
import { Layout } from './components/layout/Layout';
import { useRDKit } from './hooks/useRDKit';

// Lazy-loaded page components for code splitting
const SingleValidationPage = lazy(() =>
  import('./pages/SingleValidation').then(module => ({ default: module.SingleValidationPage }))
);
const BatchValidationPage = lazy(() =>
  import('./pages/BatchValidation').then(module => ({ default: module.BatchValidationPage }))
);

/**
 * Page loading fallback with teal chemistry-themed spinner
 */
function PageLoader() {
  return (
    <div className="flex items-center justify-center h-64">
      <div className="text-center">
        <div className="relative">
          {/* Outer ring */}
          <div className="w-16 h-16 rounded-full border-4 border-chem-primary/20"></div>
          {/* Spinning ring */}
          <div className="absolute top-0 left-0 w-16 h-16 rounded-full border-4 border-transparent border-t-chem-primary animate-spin"></div>
          {/* Inner molecule icon */}
          <div className="absolute top-1/2 left-1/2 -translate-x-1/2 -translate-y-1/2">
            <svg className="w-6 h-6 text-chem-primary" viewBox="0 0 24 24" fill="currentColor">
              <circle cx="12" cy="12" r="3" />
              <circle cx="12" cy="4" r="2" />
              <circle cx="20" cy="12" r="2" />
              <circle cx="12" cy="20" r="2" />
              <circle cx="4" cy="12" r="2" />
              <line x1="12" y1="9" x2="12" y2="6" stroke="currentColor" strokeWidth="1.5" />
              <line x1="15" y1="12" x2="18" y2="12" stroke="currentColor" strokeWidth="1.5" />
              <line x1="12" y1="15" x2="12" y2="18" stroke="currentColor" strokeWidth="1.5" />
              <line x1="9" y1="12" x2="6" y2="12" stroke="currentColor" strokeWidth="1.5" />
            </svg>
          </div>
        </div>
        <p className="mt-4 text-chem-primary font-medium">Loading...</p>
      </div>
    </div>
  );
}

/**
 * Navigation tabs with teal chemistry theme
 */
function NavigationTabs() {
  return (
    <div className="flex justify-center mb-6">
      <nav className="flex bg-chem-primary/5 rounded-xl p-1.5 border border-chem-primary/10">
        <NavLink
          to="/"
          end
          className={({ isActive }) =>
            `px-5 py-2.5 rounded-lg text-sm font-medium transition-all duration-200 ${
              isActive
                ? 'bg-white text-chem-primary shadow-md shadow-chem-primary/10'
                : 'text-chem-dark/60 hover:text-chem-primary hover:bg-white/50'
            }`
          }
        >
          <span className="flex items-center gap-2">
            <svg className="w-4 h-4" viewBox="0 0 24 24" fill="none" stroke="currentColor" strokeWidth="2">
              <circle cx="12" cy="12" r="3" />
              <path d="M12 2v4M12 18v4M2 12h4M18 12h4" />
            </svg>
            Single Validation
          </span>
        </NavLink>
        <NavLink
          to="/batch"
          className={({ isActive }) =>
            `px-5 py-2.5 rounded-lg text-sm font-medium transition-all duration-200 ${
              isActive
                ? 'bg-white text-chem-primary shadow-md shadow-chem-primary/10'
                : 'text-chem-dark/60 hover:text-chem-primary hover:bg-white/50'
            }`
          }
        >
          <span className="flex items-center gap-2">
            <svg className="w-4 h-4" viewBox="0 0 24 24" fill="none" stroke="currentColor" strokeWidth="2">
              <rect x="3" y="3" width="7" height="7" rx="1" />
              <rect x="14" y="3" width="7" height="7" rx="1" />
              <rect x="3" y="14" width="7" height="7" rx="1" />
              <rect x="14" y="14" width="7" height="7" rx="1" />
            </svg>
            Batch Processing
          </span>
        </NavLink>
      </nav>
    </div>
  );
}

/**
 * Main application component with lazy-loaded routes
 */
function App() {
  const { rdkit, loading: rdkitLoading, error: rdkitError } = useRDKit();

  if (rdkitLoading) {
    return (
      <Layout>
        <div className="flex items-center justify-center h-64">
          <div className="text-center">
            <div className="relative">
              <div className="w-16 h-16 rounded-full border-4 border-chem-primary/20"></div>
              <div className="absolute top-0 left-0 w-16 h-16 rounded-full border-4 border-transparent border-t-chem-primary animate-spin"></div>
            </div>
            <p className="mt-4 text-chem-primary font-medium">Loading RDKit.js...</p>
          </div>
        </div>
      </Layout>
    );
  }

  if (rdkitError) {
    return (
      <Layout>
        <div className="bg-status-error/10 border border-status-error/30 rounded-xl p-6 text-center">
          <div className="w-12 h-12 mx-auto mb-3 rounded-full bg-status-error/20 flex items-center justify-center">
            <svg className="w-6 h-6 text-status-error" viewBox="0 0 24 24" fill="none" stroke="currentColor" strokeWidth="2">
              <circle cx="12" cy="12" r="10" />
              <path d="M15 9l-6 6M9 9l6 6" />
            </svg>
          </div>
          <h2 className="text-status-error font-semibold mb-2">Failed to load RDKit.js</h2>
          <p className="text-status-error/80 text-sm">{rdkitError}</p>
        </div>
      </Layout>
    );
  }

  return (
    <Router>
      <Layout>
        {rdkit && (
          <div className="text-center mb-4">
            <p className="text-xs text-chem-dark/40">
              RDKit.js {rdkit.version()}
            </p>
          </div>
        )}
        <NavigationTabs />
        <Suspense fallback={<PageLoader />}>
          <Routes>
            <Route path="/" element={<SingleValidationPage />} />
            <Route path="/batch" element={<BatchValidationPage />} />
          </Routes>
        </Suspense>
      </Layout>
    </Router>
  );
}

export default App;
