import { useState } from 'react';
import { StructureInput } from '../components/molecules/StructureInput';
import { MoleculeViewer } from '../components/molecules/MoleculeViewer';
import { ValidationResults } from '../components/validation/ValidationResults';
import { AlertResults } from '../components/alerts/AlertResults';
import { useValidation } from '../hooks/useValidation';
import { alertsApi } from '../services/api';
import type { AlertScreenResponse, AlertError } from '../types/alerts';

const EXAMPLE_MOLECULES = [
  {
    name: 'Aspirin',
    smiles: 'CC(=O)Oc1ccccc1C(=O)O',
  },
  {
    name: 'Caffeine',
    smiles: 'CN1C=NC2=C1C(=O)N(C(=O)N2C)C',
  },
  {
    name: 'Ethanol',
    smiles: 'CCO',
  },
  {
    name: 'Rhodanine (PAINS)',
    smiles: 'O=C1NC(=S)SC1',
  },
  {
    name: 'Quinone (PAINS)',
    smiles: 'O=C1C=CC(=O)C=C1',
  },
];

export function SingleValidationPage() {
  const [molecule, setMolecule] = useState('');
  const [highlightedAtoms, setHighlightedAtoms] = useState<number[]>([]);
  const { validate, result, error, isLoading, reset } = useValidation();

  // Alert screening state
  const [alertResult, setAlertResult] = useState<AlertScreenResponse | null>(null);
  const [alertError, setAlertError] = useState<AlertError | null>(null);
  const [alertsLoading, setAlertsLoading] = useState(false);
  const [selectedCatalogs, setSelectedCatalogs] = useState<string[]>(['PAINS', 'BRENK']);

  const handleValidate = async () => {
    if (!molecule.trim()) return;

    await validate({
      molecule: molecule.trim(),
      format: 'auto',
    });
  };

  const handleScreenAlerts = async () => {
    if (!molecule.trim()) return;

    setAlertsLoading(true);
    setAlertError(null);

    try {
      const response = await alertsApi.screenAlerts({
        molecule: molecule.trim(),
        format: 'auto',
        catalogs: selectedCatalogs,
      });
      setAlertResult(response);
    } catch (err) {
      setAlertError(err as AlertError);
      setAlertResult(null);
    } finally {
      setAlertsLoading(false);
    }
  };

  const handleExampleClick = (smiles: string) => {
    setMolecule(smiles);
    reset();
    setAlertResult(null);
    setAlertError(null);
    setHighlightedAtoms([]);
  };

  const handleReset = () => {
    setMolecule('');
    reset();
    setAlertResult(null);
    setAlertError(null);
    setHighlightedAtoms([]);
  };

  const toggleCatalog = (catalog: string) => {
    setSelectedCatalogs((prev) =>
      prev.includes(catalog)
        ? prev.filter((c) => c !== catalog)
        : [...prev, catalog]
    );
  };

  return (
    <div className="max-w-7xl mx-auto space-y-6">
      {/* Header */}
      <div className="text-center">
        <h2 className="text-3xl font-bold text-gray-900">
          Single Molecule Validation
        </h2>
        <p className="text-gray-500 mt-2">
          Validate chemical structures with comprehensive checks
        </p>
      </div>

      {/* Example molecules */}
      <div className="flex flex-wrap gap-2 justify-center">
        <span className="text-sm text-gray-600 self-center">Try examples:</span>
        {EXAMPLE_MOLECULES.map((example) => (
          <button
            key={example.name}
            onClick={() => handleExampleClick(example.smiles)}
            className="px-3 py-1 text-sm bg-blue-50 hover:bg-blue-100 text-blue-700 rounded-md transition-colors"
          >
            {example.name}
          </button>
        ))}
      </div>

      {/* Input and Preview Section */}
      <div className="grid grid-cols-1 lg:grid-cols-2 gap-6">
        {/* Input */}
        <div className="space-y-4">
          <div className="bg-white rounded-lg shadow-md p-6">
            <h3 className="font-semibold text-gray-900 mb-4">Input Structure</h3>
            <StructureInput value={molecule} onChange={setMolecule} />
            <div className="mt-4 flex gap-2">
              <button
                onClick={handleValidate}
                disabled={!molecule.trim() || isLoading}
                className="flex-1 px-6 py-2.5 bg-blue-600 hover:bg-blue-700 disabled:bg-gray-300 disabled:cursor-not-allowed text-white font-medium rounded-lg transition-colors"
              >
                {isLoading ? 'Validating...' : 'Validate'}
              </button>
              <button
                onClick={handleReset}
                disabled={!molecule && !result && !error && !alertResult}
                className="px-6 py-2.5 bg-gray-200 hover:bg-gray-300 disabled:bg-gray-100 disabled:cursor-not-allowed text-gray-700 font-medium rounded-lg transition-colors"
              >
                Reset
              </button>
            </div>

            {/* Alert Screening Section */}
            <div className="mt-4 pt-4 border-t border-gray-100">
              <h4 className="font-medium text-gray-700 mb-2">Structural Alerts</h4>
              <div className="flex flex-wrap gap-2 mb-3">
                {['PAINS', 'BRENK', 'NIH', 'ZINC'].map((catalog) => (
                  <label
                    key={catalog}
                    className={`px-3 py-1 text-sm rounded-md cursor-pointer transition-colors ${
                      selectedCatalogs.includes(catalog)
                        ? 'bg-amber-100 text-amber-800 border border-amber-300'
                        : 'bg-gray-100 text-gray-600 border border-gray-200'
                    }`}
                  >
                    <input
                      type="checkbox"
                      checked={selectedCatalogs.includes(catalog)}
                      onChange={() => toggleCatalog(catalog)}
                      className="sr-only"
                    />
                    {catalog}
                  </label>
                ))}
              </div>
              <button
                onClick={handleScreenAlerts}
                disabled={!molecule.trim() || alertsLoading || selectedCatalogs.length === 0}
                className="w-full px-4 py-2 bg-amber-500 hover:bg-amber-600 disabled:bg-gray-300 disabled:cursor-not-allowed text-white font-medium rounded-lg transition-colors"
              >
                {alertsLoading ? 'Screening...' : 'Screen for Alerts'}
              </button>
            </div>
          </div>
        </div>

        {/* Preview */}
        <div className="space-y-4">
          <div className="bg-white rounded-lg shadow-md p-6">
            <h3 className="font-semibold text-gray-900 mb-4">Structure Preview</h3>
            <MoleculeViewer
              smiles={result?.molecule_info.canonical_smiles || alertResult?.molecule_info.canonical_smiles || molecule}
              highlightAtoms={highlightedAtoms}
              width={400}
              height={300}
              className="mx-auto"
            />
            {highlightedAtoms.length > 0 && (
              <div className="mt-2 text-xs text-center text-amber-600 font-medium">
                Highlighting atoms: {highlightedAtoms.join(', ')}
              </div>
            )}
          </div>
        </div>
      </div>

      {/* Loading State */}
      {(isLoading || alertsLoading) && (
        <div className="bg-blue-50 border border-blue-200 rounded-lg p-8 text-center">
          <div className="animate-spin rounded-full h-12 w-12 border-b-2 border-blue-600 mx-auto mb-4"></div>
          <p className="text-blue-700 font-medium">
            {isLoading ? 'Running validation checks...' : 'Screening for structural alerts...'}
          </p>
        </div>
      )}

      {/* Error State */}
      {(error || alertError) && !isLoading && !alertsLoading && (
        <div className="bg-red-50 border border-red-200 rounded-lg p-6">
          <div className="flex items-start gap-3">
            <span className="text-2xl">!</span>
            <div className="flex-1">
              <h3 className="font-semibold text-red-900 mb-2">
                {(error?.error || alertError?.error)?.includes('parse') ? 'Parse Error' : 'Error'}
              </h3>
              <p className="text-red-700 text-sm">{error?.error || alertError?.error}</p>
              {(error?.details || alertError?.details) && (
                <div className="mt-3 text-xs text-red-600">
                  {(error?.details?.errors || alertError?.details?.errors) && (
                    <div>
                      <strong>Errors:</strong>
                      <ul className="list-disc list-inside mt-1">
                        {(error?.details?.errors || alertError?.details?.errors)?.map((err, i) => (
                          <li key={i}>{err}</li>
                        ))}
                      </ul>
                    </div>
                  )}
                </div>
              )}
            </div>
          </div>
        </div>
      )}

      {/* Validation Results */}
      {result && !isLoading && (
        <ValidationResults
          result={result}
          onHighlightAtoms={setHighlightedAtoms}
        />
      )}

      {/* Alert Results */}
      {alertResult && !alertsLoading && (
        <div className="bg-white rounded-lg shadow-md p-6">
          <h3 className="text-lg font-semibold text-gray-900 mb-4">
            Structural Alert Screening Results
          </h3>
          <AlertResults
            alerts={alertResult.alerts}
            screenedCatalogs={alertResult.screened_catalogs}
            educationalNote={alertResult.educational_note}
            onHighlightAtoms={setHighlightedAtoms}
          />
          <div className="mt-4 text-xs text-gray-500 text-right">
            Completed in {alertResult.execution_time_ms}ms
          </div>
        </div>
      )}
    </div>
  );
}
