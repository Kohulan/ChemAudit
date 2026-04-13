import { useState } from 'react';
import { integrationsApi } from '../../services/api';
import { logger } from '../../lib/logger';
import type { PubChemResult, ChEMBLResult, COCONUTResult, WikidataResult, ConsistencyResult } from '../../types/integrations';
import { DatabaseComparisonPanel } from './DatabaseComparisonPanel';
import { DatabaseLookupResults } from './DatabaseLookupResults';
import { ClayButton } from '../ui/ClayButton';

interface DatabaseLookupProps {
  inchikey?: string;
  smiles?: string;
}

export function DatabaseLookup({ inchikey, smiles }: DatabaseLookupProps) {
  const [isLoading, setIsLoading] = useState(false);
  const [results, setResults] = useState<{
    pubchem: PubChemResult | null;
    chembl: ChEMBLResult | null;
    coconut: COCONUTResult | null;
    wikidata: WikidataResult | null;
  } | null>(null);
  const [error, setError] = useState<string | null>(null);
  const [comparisonResult, setComparisonResult] = useState<ConsistencyResult | null>(null);
  const [isComparing, setIsComparing] = useState(false);

  const handleLookup = async () => {
    if (!inchikey && !smiles) return;

    setIsLoading(true);
    setError(null);

    try {
      const result = await integrationsApi.lookupAll({ inchikey, smiles });
      setResults(result);
    } catch (err) {
      setError('Failed to look up databases');
      logger.error('Database lookup error:', err);
    } finally {
      setIsLoading(false);
    }
  };

  const handleCompare = async () => {
    setIsComparing(true);
    try {
      const result = await integrationsApi.compareAcrossDatabases({ smiles, inchikey });
      setComparisonResult(result);
    } catch (err) {
      logger.error('Comparison error:', err);
    } finally {
      setIsComparing(false);
    }
  };

  if (!inchikey && !smiles) {
    return null;
  }

  return (
    <div className="bg-white rounded-lg shadow-md p-6">
      <div className="flex items-center justify-between mb-4">
        <h3 className="text-lg font-semibold text-gray-900">
          Database Cross-Reference
        </h3>
        {!results && (
          <ClayButton
            variant="primary"
            size="sm"
            onClick={handleLookup}
            disabled={isLoading}
            loading={isLoading}
          >
            {isLoading ? 'Looking up...' : 'Look Up in Databases'}
          </ClayButton>
        )}
      </div>

      {!results && !isLoading && (
        <p className="text-sm text-gray-500">
          Check if this molecule exists in PubChem, ChEMBL, COCONUT, or Wikidata.
        </p>
      )}

      {isLoading && (
        <div className="flex items-center justify-center py-8">
          <div className="animate-spin rounded-full h-8 w-8 border-b-2 border-indigo-600"></div>
          <span className="ml-3 text-gray-600">Querying databases...</span>
        </div>
      )}

      {error && (
        <div className="bg-red-50 border border-red-200 rounded-lg p-4 text-red-700 text-sm">
          {error}
        </div>
      )}

      {results && (
        <div className="space-y-4">
          <DatabaseLookupResults results={results} />

          {!comparisonResult && (
            <div className="flex justify-end">
              <ClayButton
                variant="primary"
                size="sm"
                onClick={handleCompare}
                disabled={isComparing}
                loading={isComparing}
              >
                {isComparing ? 'Comparing...' : 'Compare Across Databases'}
              </ClayButton>
            </div>
          )}
          {comparisonResult && <DatabaseComparisonPanel result={comparisonResult} />}
        </div>
      )}
    </div>
  );
}
