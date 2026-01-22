import { useState } from 'react';
import { StructureInput } from '../components/molecules/StructureInput';
import { MoleculeViewer } from '../components/molecules/MoleculeViewer';
import { ValidationResults } from '../components/validation/ValidationResults';
import { useValidation } from '../hooks/useValidation';

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
    name: 'Chiral Center Test',
    smiles: 'CC(O)CC',
  },
];

export function SingleValidationPage() {
  const [molecule, setMolecule] = useState('');
  const [highlightedAtoms, setHighlightedAtoms] = useState<number[]>([]);
  const { validate, result, error, isLoading, reset } = useValidation();

  const handleValidate = async () => {
    if (!molecule.trim()) return;

    await validate({
      molecule: molecule.trim(),
      format: 'auto',
    });
  };

  const handleExampleClick = (smiles: string) => {
    setMolecule(smiles);
    reset();
    setHighlightedAtoms([]);
  };

  const handleReset = () => {
    setMolecule('');
    reset();
    setHighlightedAtoms([]);
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
                disabled={!molecule && !result && !error}
                className="px-6 py-2.5 bg-gray-200 hover:bg-gray-300 disabled:bg-gray-100 disabled:cursor-not-allowed text-gray-700 font-medium rounded-lg transition-colors"
              >
                Reset
              </button>
            </div>
          </div>
        </div>

        {/* Preview */}
        <div className="space-y-4">
          <div className="bg-white rounded-lg shadow-md p-6">
            <h3 className="font-semibold text-gray-900 mb-4">Structure Preview</h3>
            <MoleculeViewer
              smiles={result?.molecule_info.canonical_smiles || molecule}
              highlightAtoms={highlightedAtoms}
              width={400}
              height={300}
              className="mx-auto"
            />
            {highlightedAtoms.length > 0 && (
              <div className="mt-2 text-xs text-center text-gray-500">
                Highlighting atoms: {highlightedAtoms.join(', ')}
              </div>
            )}
          </div>
        </div>
      </div>

      {/* Loading State */}
      {isLoading && (
        <div className="bg-blue-50 border border-blue-200 rounded-lg p-8 text-center">
          <div className="animate-spin rounded-full h-12 w-12 border-b-2 border-blue-600 mx-auto mb-4"></div>
          <p className="text-blue-700 font-medium">Running validation checks...</p>
        </div>
      )}

      {/* Error State */}
      {error && !isLoading && (
        <div className="bg-red-50 border border-red-200 rounded-lg p-6">
          <div className="flex items-start gap-3">
            <span className="text-2xl">⚠️</span>
            <div className="flex-1">
              <h3 className="font-semibold text-red-900 mb-2">
                {error.error?.includes('parse') ? 'Parse Error' : 'Validation Error'}
              </h3>
              <p className="text-red-700 text-sm">{error.error}</p>
              {error.details && (
                <div className="mt-3 text-xs text-red-600">
                  {error.details.errors && (
                    <div>
                      <strong>Errors:</strong>
                      <ul className="list-disc list-inside mt-1">
                        {error.details.errors.map((err, i) => (
                          <li key={i}>{err}</li>
                        ))}
                      </ul>
                    </div>
                  )}
                  {error.details.warnings && error.details.warnings.length > 0 && (
                    <div className="mt-2">
                      <strong>Warnings:</strong>
                      <ul className="list-disc list-inside mt-1">
                        {error.details.warnings.map((warn, i) => (
                          <li key={i}>{warn}</li>
                        ))}
                      </ul>
                    </div>
                  )}
                  {error.details.format_detected && (
                    <div className="mt-2">
                      <strong>Format detected:</strong> {error.details.format_detected}
                    </div>
                  )}
                </div>
              )}
            </div>
          </div>
        </div>
      )}

      {/* Results */}
      {result && !isLoading && (
        <ValidationResults
          result={result}
          onHighlightAtoms={setHighlightedAtoms}
        />
      )}
    </div>
  );
}
