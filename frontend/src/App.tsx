import { BrowserRouter as Router, Routes, Route, NavLink } from 'react-router-dom';
import { Layout } from './components/layout/Layout';
import { SingleValidationPage } from './pages/SingleValidation';
import { BatchValidationPage } from './pages/BatchValidation';
import { useRDKit } from './hooks/useRDKit';

function NavigationTabs() {
  return (
    <div className="flex justify-center mb-6">
      <nav className="flex bg-gray-100 rounded-lg p-1">
        <NavLink
          to="/"
          end
          className={({ isActive }) =>
            `px-4 py-2 rounded-md text-sm font-medium transition-colors ${
              isActive
                ? 'bg-white text-blue-600 shadow-sm'
                : 'text-gray-600 hover:text-gray-900'
            }`
          }
        >
          Single Validation
        </NavLink>
        <NavLink
          to="/batch"
          className={({ isActive }) =>
            `px-4 py-2 rounded-md text-sm font-medium transition-colors ${
              isActive
                ? 'bg-white text-blue-600 shadow-sm'
                : 'text-gray-600 hover:text-gray-900'
            }`
          }
        >
          Batch Processing
        </NavLink>
      </nav>
    </div>
  );
}

function App() {
  const { rdkit, loading: rdkitLoading, error: rdkitError } = useRDKit();

  if (rdkitLoading) {
    return (
      <Layout>
        <div className="flex items-center justify-center h-64">
          <div className="text-center">
            <div className="animate-spin rounded-full h-12 w-12 border-b-2 border-blue-600 mx-auto mb-4"></div>
            <p className="text-gray-500">Loading RDKit.js...</p>
          </div>
        </div>
      </Layout>
    );
  }

  if (rdkitError) {
    return (
      <Layout>
        <div className="bg-red-50 border border-red-200 rounded-lg p-6 text-center">
          <h2 className="text-red-700 font-semibold mb-2">Failed to load RDKit.js</h2>
          <p className="text-red-600 text-sm">{rdkitError}</p>
        </div>
      </Layout>
    );
  }

  return (
    <Router>
      <Layout>
        {rdkit && (
          <div className="text-center mb-4">
            <p className="text-xs text-gray-400">
              RDKit.js {rdkit.version()}
            </p>
          </div>
        )}
        <NavigationTabs />
        <Routes>
          <Route path="/" element={<SingleValidationPage />} />
          <Route path="/batch" element={<BatchValidationPage />} />
        </Routes>
      </Layout>
    </Router>
  );
}

export default App;
