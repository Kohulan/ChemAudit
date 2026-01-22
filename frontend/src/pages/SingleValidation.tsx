import { useState } from 'react';
import { StructureInput } from '../components/molecules/StructureInput';
import { MoleculeViewer } from '../components/molecules/MoleculeViewer';
import { ValidationResults } from '../components/validation/ValidationResults';
import { AlertResults } from '../components/alerts/AlertResults';
import { ScoringResults } from '../components/scoring/ScoringResults';
import { StandardizationResults } from '../components/standardization/StandardizationResults';
import { useValidation } from '../hooks/useValidation';
import { alertsApi, scoringApi, standardizationApi } from '../services/api';
import type { AlertScreenResponse, AlertError } from '../types/alerts';
import type { ScoringResponse, ScoringError } from '../types/scoring';
import type { StandardizeResponse, StandardizeError } from '../types/standardization';

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
    name: 'Morphine',
    smiles: 'CN1CCC23C4=C5C=CC(O)=C4OC2C(O)C=CC3C1C5',
  },
  {
    name: 'Amine HCl (salt)',
    smiles: 'CCN.Cl',
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

  // Scoring state
  const [scoringResult, setScoringResult] = useState<ScoringResponse | null>(null);
  const [scoringError, setScoringError] = useState<ScoringError | null>(null);
  const [scoringLoading, setScoringLoading] = useState(false);

  // Standardization state
  const [standardizationResult, setStandardizationResult] = useState<StandardizeResponse | null>(null);
  const [standardizationError, setStandardizationError] = useState<StandardizeError | null>(null);
  const [standardizationLoading, setStandardizationLoading] = useState(false);
  const [includeTautomer, setIncludeTautomer] = useState(false);

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

  const handleCalculateScores = async () => {
    if (!molecule.trim()) return;

    setScoringLoading(true);
    setScoringError(null);

    try {
      const response = await scoringApi.getScoring(molecule.trim(), 'auto');
      setScoringResult(response);
    } catch (err) {
      setScoringError(err as ScoringError);
      setScoringResult(null);
    } finally {
      setScoringLoading(false);
    }
  };

  const handleStandardize = async () => {
    if (!molecule.trim()) return;

    setStandardizationLoading(true);
    setStandardizationError(null);

    try {
      const response = await standardizationApi.standardize({
        molecule: molecule.trim(),
        format: 'auto',
        options: {
          include_tautomer: includeTautomer,
          preserve_stereo: true,
        },
      });
      setStandardizationResult(response);
    } catch (err) {
      setStandardizationError(err as StandardizeError);
      setStandardizationResult(null);
    } finally {
      setStandardizationLoading(false);
    }
  };

  const handleExampleClick = (smiles: string) => {
    setMolecule(smiles);
    reset();
    setAlertResult(null);
    setAlertError(null);
    setScoringResult(null);
    setScoringError(null);
    setStandardizationResult(null);
    setStandardizationError(null);
    setHighlightedAtoms([]);
  };

  const handleReset = () => {
    setMolecule('');
    reset();
    setAlertResult(null);
    setAlertError(null);
    setScoringResult(null);
    setScoringError(null);
    setStandardizationResult(null);
    setStandardizationError(null);
    setHighlightedAtoms([]);
  };

  const toggleCatalog = (catalog: string) => {
    setSelectedCatalogs((prev) =>
      prev.includes(catalog)
        ? prev.filter((c) => c !== catalog)
        : [...prev, catalog]
    );
  };

  const isAnyLoading = isLoading || alertsLoading || scoringLoading || standardizationLoading;

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
                disabled={!molecule.trim() || isAnyLoading}
                className="flex-1 px-6 py-2.5 bg-blue-600 hover:bg-blue-700 disabled:bg-gray-300 disabled:cursor-not-allowed text-white font-medium rounded-lg transition-colors"
              >
                {isLoading ? 'Validating...' : 'Validate'}
              </button>
              <button
                onClick={handleReset}
                disabled={!molecule && !result && !error && !alertResult && !scoringResult && !standardizationResult}
                className="px-6 py-2.5 bg-gray-200 hover:bg-gray-300 disabled:bg-gray-100 disabled:cursor-not-allowed text-gray-700 font-medium rounded-lg transition-colors"
              >
                Reset
              </button>
            </div>

            {/* Scoring Section */}
            <div className="mt-4 pt-4 border-t border-gray-100">
              <h4 className="font-medium text-gray-700 mb-2">ML-Readiness Scoring</h4>
              <p className="text-xs text-gray-500 mb-3">
                Calculate ML-readiness, NP-likeness, and extract Murcko scaffold
              </p>
              <button
                onClick={handleCalculateScores}
                disabled={!molecule.trim() || isAnyLoading}
                className="w-full px-4 py-2 bg-purple-600 hover:bg-purple-700 disabled:bg-gray-300 disabled:cursor-not-allowed text-white font-medium rounded-lg transition-colors"
              >
                {scoringLoading ? 'Calculating...' : 'Calculate Scores'}
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
                disabled={!molecule.trim() || isAnyLoading || selectedCatalogs.length === 0}
                className="w-full px-4 py-2 bg-amber-500 hover:bg-amber-600 disabled:bg-gray-300 disabled:cursor-not-allowed text-white font-medium rounded-lg transition-colors"
              >
                {alertsLoading ? 'Screening...' : 'Screen for Alerts'}
              </button>
            </div>

            {/* Standardization Section */}
            <div className="mt-4 pt-4 border-t border-gray-100">
              <h4 className="font-medium text-gray-700 mb-2">Standardization</h4>
              <p className="text-xs text-gray-500 mb-3">
                Standardize structure using ChEMBL pipeline (salt stripping, normalization)
              </p>
              <label className="flex items-center gap-2 mb-3 cursor-pointer group">
                <input
                  type="checkbox"
                  checked={includeTautomer}
                  onChange={(e) => setIncludeTautomer(e.target.checked)}
                  className="rounded border-gray-300 text-green-600 focus:ring-green-500"
                />
                <span className="text-sm text-gray-600">Include tautomer canonicalization</span>
                <span className="relative">
                  <span className="text-amber-500 cursor-help" title="May lose E/Z stereochemistry">!</span>
                  <span className="absolute left-6 top-0 w-48 p-2 bg-amber-50 border border-amber-200 rounded text-xs text-amber-700 hidden group-hover:block z-10">
                    Warning: Tautomer canonicalization may remove E/Z double bond stereochemistry
                  </span>
                </span>
              </label>
              <button
                onClick={handleStandardize}
                disabled={!molecule.trim() || isAnyLoading}
                className="w-full px-4 py-2 bg-green-600 hover:bg-green-700 disabled:bg-gray-300 disabled:cursor-not-allowed text-white font-medium rounded-lg transition-colors"
              >
                {standardizationLoading ? 'Standardizing...' : 'Standardize'}
              </button>
            </div>
          </div>
        </div>

        {/* Preview */}
        <div className="space-y-4">
          <div className="bg-white rounded-lg shadow-md p-6">
            <h3 className="font-semibold text-gray-900 mb-4">Structure Preview</h3>
            <MoleculeViewer
              smiles={result?.molecule_info.canonical_smiles || alertResult?.molecule_info.canonical_smiles || scoringResult?.molecule_info.canonical_smiles || molecule}
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
      {isAnyLoading && (
        <div className="bg-blue-50 border border-blue-200 rounded-lg p-8 text-center">
          <div className="animate-spin rounded-full h-12 w-12 border-b-2 border-blue-600 mx-auto mb-4"></div>
          <p className="text-blue-700 font-medium">
            {isLoading ? 'Running validation checks...' :
             scoringLoading ? 'Calculating scores...' :
             standardizationLoading ? 'Running standardization pipeline...' :
             'Screening for structural alerts...'}
          </p>
        </div>
      )}

      {/* Error State */}
      {(error || alertError || scoringError || standardizationError) && !isAnyLoading && (
        <div className="bg-red-50 border border-red-200 rounded-lg p-6">
          <div className="flex items-start gap-3">
            <span className="text-2xl">!</span>
            <div className="flex-1">
              <h3 className="font-semibold text-red-900 mb-2">
                {(error?.error || alertError?.error || scoringError?.error || standardizationError?.error)?.includes('parse') ? 'Parse Error' : 'Error'}
              </h3>
              <p className="text-red-700 text-sm">{error?.error || alertError?.error || scoringError?.error || standardizationError?.error}</p>
              {(error?.details || alertError?.details || scoringError?.details || standardizationError?.details) && (
                <div className="mt-3 text-xs text-red-600">
                  {(error?.details?.errors || alertError?.details?.errors || scoringError?.details?.errors || standardizationError?.details?.errors) && (
                    <div>
                      <strong>Errors:</strong>
                      <ul className="list-disc list-inside mt-1">
                        {(error?.details?.errors || alertError?.details?.errors || scoringError?.details?.errors || standardizationError?.details?.errors)?.map((err, i) => (
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

      {/* Scoring Results */}
      {scoringResult && !scoringLoading && (
        <div className="bg-white rounded-lg shadow-md p-6">
          <ScoringResults scoringResponse={scoringResult} />
        </div>
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

      {/* Standardization Results */}
      {standardizationResult && !standardizationLoading && (
        <div className="bg-white rounded-lg shadow-md p-6">
          <h3 className="text-lg font-semibold text-gray-900 mb-4">
            Standardization Results
          </h3>
          <StandardizationResults result={standardizationResult.result} />
          <div className="mt-4 text-xs text-gray-500 text-right">
            Completed in {standardizationResult.execution_time_ms}ms
          </div>
        </div>
      )}
    </div>
  );
}
