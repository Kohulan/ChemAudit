import { useState } from 'react';
import { Layout } from './components/layout/Layout';
import { StructureInput } from './components/molecules/StructureInput';
import { MoleculeViewer } from './components/molecules/MoleculeViewer';
import { useRDKit } from './hooks/useRDKit';

function App() {
  const [molecule, setMolecule] = useState('');
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
    <Layout>
      <div className="max-w-4xl mx-auto space-y-8">
        <div className="text-center">
          <h2 className="text-2xl font-bold text-gray-900">Chemical Structure Validation</h2>
          <p className="text-gray-500 mt-2">
            Enter a SMILES, InChI, or MOL block to validate
          </p>
          {rdkit && (
            <p className="text-xs text-gray-400 mt-1">
              RDKit.js {rdkit.version()}
            </p>
          )}
        </div>

        <div className="grid grid-cols-1 md:grid-cols-2 gap-6">
          <div className="space-y-4">
            <h3 className="font-medium text-gray-900">Input</h3>
            <StructureInput
              value={molecule}
              onChange={setMolecule}
            />
          </div>

          <div className="space-y-4">
            <h3 className="font-medium text-gray-900">Preview</h3>
            <MoleculeViewer
              smiles={molecule}
              width={400}
              height={300}
              className="mx-auto"
            />
          </div>
        </div>
      </div>
    </Layout>
  );
}

export default App;
