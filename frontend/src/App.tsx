import { lazy, Suspense } from 'react';
import { BrowserRouter as Router, Routes, Route, NavLink } from 'react-router-dom';
import { AnimatePresence, motion } from 'framer-motion';
import { Atom, Grid3X3 } from 'lucide-react';
import { Layout } from './components/layout/Layout';
import { MoleculeLoader } from './components/ui/MoleculeLoader';
import { useRDKit } from './hooks/useRDKit';
import { cn } from './lib/utils';

// Lazy-loaded page components for code splitting
const SingleValidationPage = lazy(() =>
  import('./pages/SingleValidation').then(module => ({ default: module.SingleValidationPage }))
);
const BatchValidationPage = lazy(() =>
  import('./pages/BatchValidation').then(module => ({ default: module.BatchValidationPage }))
);

/**
 * Page loading fallback with molecule loader
 */
function PageLoaderFallback() {
  return (
    <div className="flex items-center justify-center h-64">
      <MoleculeLoader size="lg" text="Loading..." />
    </div>
  );
}

/**
 * Navigation tabs with clay styling
 */
function NavigationTabs() {
  return (
    <div className="flex justify-center mb-8">
      <nav className="clay-card-sm p-1.5 flex gap-1">
        <NavLink
          to="/"
          end
          className={({ isActive }) =>
            cn(
              'px-5 py-2.5 rounded-xl text-sm font-medium transition-all duration-200',
              'flex items-center gap-2',
              isActive
                ? 'bg-white dark:bg-white/10 text-chem-primary-600 dark:text-chem-primary-400 shadow-md'
                : 'text-text-secondary hover:text-text-primary hover:bg-white/50 dark:hover:bg-white/5'
            )
          }
        >
          <Atom className="w-4 h-4" />
          Single Validation
        </NavLink>
        <NavLink
          to="/batch"
          className={({ isActive }) =>
            cn(
              'px-5 py-2.5 rounded-xl text-sm font-medium transition-all duration-200',
              'flex items-center gap-2',
              isActive
                ? 'bg-white dark:bg-white/10 text-chem-primary-600 dark:text-chem-primary-400 shadow-md'
                : 'text-text-secondary hover:text-text-primary hover:bg-white/50 dark:hover:bg-white/5'
            )
          }
        >
          <Grid3X3 className="w-4 h-4" />
          Batch Processing
        </NavLink>
      </nav>
    </div>
  );
}

/**
 * RDKit loading/error states
 */
function RDKitStatus() {
  const { rdkit, loading, error } = useRDKit();

  if (loading) {
    return (
      <div className="flex items-center justify-center h-64">
        <MoleculeLoader size="lg" text="Loading RDKit.js..." />
      </div>
    );
  }

  if (error) {
    return (
      <motion.div
        initial={{ opacity: 0, y: 10 }}
        animate={{ opacity: 1, y: 0 }}
        className="clay-card p-8 text-center max-w-md mx-auto"
      >
        <div className={cn(
          'w-14 h-14 mx-auto mb-4 rounded-2xl flex items-center justify-center',
          'bg-red-500/10 dark:bg-red-400/10'
        )}>
          <svg
            className="w-7 h-7 text-red-500 dark:text-red-400"
            viewBox="0 0 24 24"
            fill="none"
            stroke="currentColor"
            strokeWidth="2"
          >
            <circle cx="12" cy="12" r="10" />
            <path d="M15 9l-6 6M9 9l6 6" />
          </svg>
        </div>
        <h2 className="text-lg font-semibold text-red-600 dark:text-red-400 mb-2">
          Failed to load RDKit.js
        </h2>
        <p className="text-sm text-text-secondary">{error}</p>
      </motion.div>
    );
  }

  return rdkit ? (
    <p className="text-center text-xs text-text-muted mb-2">
      RDKit.js {rdkit.version()}
    </p>
  ) : null;
}

/**
 * Main application component
 */
function App() {
  const { loading: rdkitLoading, error: rdkitError } = useRDKit();

  // Show loading/error states before router
  if (rdkitLoading || rdkitError) {
    return (
      <Layout>
        <RDKitStatus />
      </Layout>
    );
  }

  return (
    <Router>
      <Layout>
        <RDKitStatus />
        <NavigationTabs />
        <AnimatePresence mode="wait">
          <Suspense fallback={<PageLoaderFallback />}>
            <Routes>
              <Route
                path="/"
                element={
                  <motion.div
                    key="single"
                    initial={{ opacity: 0, y: 10 }}
                    animate={{ opacity: 1, y: 0 }}
                    exit={{ opacity: 0, y: -10 }}
                    transition={{ duration: 0.2 }}
                  >
                    <SingleValidationPage />
                  </motion.div>
                }
              />
              <Route
                path="/batch"
                element={
                  <motion.div
                    key="batch"
                    initial={{ opacity: 0, y: 10 }}
                    animate={{ opacity: 1, y: 0 }}
                    exit={{ opacity: 0, y: -10 }}
                    transition={{ duration: 0.2 }}
                  >
                    <BatchValidationPage />
                  </motion.div>
                }
              />
            </Routes>
          </Suspense>
        </AnimatePresence>
      </Layout>
    </Router>
  );
}

export default App;
