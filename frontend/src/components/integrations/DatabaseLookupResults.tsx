import type { PubChemResult, ChEMBLResult, COCONUTResult, WikidataResult, SureChEMBLResult } from '../../types/integrations';
import { safeHref } from '../../lib/sanitize';
const pubchemLogo = '/assets/logos/pubchem.png';
const chemblLogo = '/assets/logos/chembl.png';
const coconutLogo = '/assets/logos/coconut.png';
const wikidataLogo = '/assets/logos/wikidata.svg';
const surechemblLogo = '/assets/logos/surechembl.png';

interface DatabaseLookupResultsProps {
  results: {
    pubchem: PubChemResult | null;
    chembl: ChEMBLResult | null;
    coconut: COCONUTResult | null;
    wikidata: WikidataResult | null;
    surechembl?: SureChEMBLResult | null;
  };
}

export function DatabaseLookupResults({ results }: DatabaseLookupResultsProps) {
  return (
    <div className="space-y-4">
      {/* PubChem Result */}
      <PubChemCard result={results.pubchem} />

      {/* ChEMBL Result */}
      <ChEMBLCard result={results.chembl} />

      {/* COCONUT Result */}
      <COCONUTCard result={results.coconut} />

      {/* Wikidata Result */}
      <WikidataCard result={results.wikidata} />

      {/* SureChEMBL Result */}
      {results.surechembl !== undefined && (
        <SureChEMBLCard result={results.surechembl ?? null} />
      )}
    </div>
  );
}

function PubChemCard({ result }: { result: PubChemResult | null }) {
  if (!result) {
    return (
      <div className="border border-[var(--color-border)] rounded-lg p-4">
        <div className="flex items-center gap-2">
          <img src={pubchemLogo} alt="PubChem" className="w-5 h-5 rounded-sm object-contain" />
          <span className="font-medium text-[var(--color-text-primary)]">PubChem</span>
          <span className="text-xs text-[var(--color-text-muted)]">Lookup unavailable</span>
        </div>
      </div>
    );
  }

  return (
    <div className={`border rounded-lg p-4 ${result.found ? 'border-blue-200 bg-blue-50 dark:border-blue-800 dark:bg-blue-900/20' : 'border-[var(--color-border)]'}`}>
      <div className="flex items-center gap-2 mb-2">
        <img src={pubchemLogo} alt="PubChem" className="w-5 h-5 rounded-sm object-contain" />
        <span className="font-medium text-[var(--color-text-primary)]">PubChem</span>
        {result.found ? (
          <span className="px-2 py-0.5 bg-blue-100 text-blue-800 dark:bg-blue-900/30 dark:text-blue-400 text-xs rounded-full">Found</span>
        ) : (
          <span className="px-2 py-0.5 bg-[var(--color-surface-sunken)] text-[var(--color-text-secondary)] text-xs rounded-full">Not Found</span>
        )}
      </div>

      {result.found && (
        <div className="mt-3 space-y-2 text-sm">
          <div className="grid grid-cols-2 gap-2">
            <div>
              <span className="text-[var(--color-text-secondary)]">CID:</span>{' '}
              <span className="font-mono text-blue-700 dark:text-blue-400">{result.cid}</span>
            </div>
            {result.iupac_name && (
              <div className="col-span-2">
                <span className="text-[var(--color-text-secondary)]">IUPAC:</span>{' '}
                <span className="text-[var(--color-text-primary)]">{result.iupac_name}</span>
              </div>
            )}
            {result.molecular_formula && (
              <div>
                <span className="text-[var(--color-text-secondary)]">Formula:</span>{' '}
                <span className="font-mono">{result.molecular_formula}</span>
              </div>
            )}
            {result.molecular_weight && (
              <div>
                <span className="text-[var(--color-text-secondary)]">MW:</span>{' '}
                <span>{result.molecular_weight.toFixed(2)}</span>
              </div>
            )}
          </div>
          {result.synonyms && result.synonyms.length > 0 && (
            <div>
              <span className="text-[var(--color-text-secondary)]">Synonyms:</span>{' '}
              <span className="text-[var(--color-text-primary)]">{result.synonyms.slice(0, 5).join(', ')}</span>
            </div>
          )}
          {result.url && (
            <a
              href={safeHref(result.url)}
              target="_blank"
              rel="noopener noreferrer"
              className="inline-flex items-center gap-1 text-blue-600 hover:text-blue-800 dark:text-blue-400 dark:hover:text-blue-300"
            >
              View on PubChem
              <svg className="w-3 h-3" fill="none" stroke="currentColor" viewBox="0 0 24 24">
                <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={2} d="M10 6H6a2 2 0 00-2 2v10a2 2 0 002 2h10a2 2 0 002-2v-4M14 4h6m0 0v6m0-6L10 14" />
              </svg>
            </a>
          )}
        </div>
      )}
    </div>
  );
}

function ChEMBLCard({ result }: { result: ChEMBLResult | null }) {
  if (!result) {
    return (
      <div className="border border-[var(--color-border)] rounded-lg p-4">
        <div className="flex items-center gap-2">
          <img src={chemblLogo} alt="ChEMBL" className="w-5 h-5 rounded-sm object-contain" />
          <span className="font-medium text-[var(--color-text-primary)]">ChEMBL</span>
          <span className="text-xs text-[var(--color-text-muted)]">Lookup unavailable</span>
        </div>
      </div>
    );
  }

  const phaseLabels: Record<number, string> = {
    0: 'Preclinical',
    1: 'Phase I',
    2: 'Phase II',
    3: 'Phase III',
    4: 'Approved',
  };

  return (
    <div className={`border rounded-lg p-4 ${result.found ? 'border-purple-200 bg-purple-50 dark:border-purple-800 dark:bg-purple-900/20' : 'border-[var(--color-border)]'}`}>
      <div className="flex items-center gap-2 mb-2">
        <img src={chemblLogo} alt="ChEMBL" className="w-5 h-5 rounded-sm object-contain" />
        <span className="font-medium text-[var(--color-text-primary)]">ChEMBL</span>
        {result.found ? (
          <span className="px-2 py-0.5 bg-purple-100 text-purple-800 dark:bg-purple-900/30 dark:text-purple-400 text-xs rounded-full">Found</span>
        ) : (
          <span className="px-2 py-0.5 bg-[var(--color-surface-sunken)] text-[var(--color-text-secondary)] text-xs rounded-full">Not Found</span>
        )}
        {result.found && result.max_phase !== undefined && result.max_phase > 0 && (
          <span className="px-2 py-0.5 bg-yellow-100 text-amber-800 dark:bg-yellow-900/30 dark:text-yellow-400 text-xs rounded-full">
            {phaseLabels[result.max_phase] || `Phase ${result.max_phase}`}
          </span>
        )}
      </div>

      {result.found && (
        <div className="mt-3 space-y-2 text-sm">
          <div className="grid grid-cols-2 gap-2">
            <div>
              <span className="text-[var(--color-text-secondary)]">ChEMBL ID:</span>{' '}
              <span className="font-mono text-purple-700 dark:text-purple-400">{result.chembl_id}</span>
            </div>
            {result.pref_name && (
              <div>
                <span className="text-[var(--color-text-secondary)]">Name:</span>{' '}
                <span className="text-[var(--color-text-primary)]">{result.pref_name}</span>
              </div>
            )}
            {result.molecule_type && (
              <div>
                <span className="text-[var(--color-text-secondary)]">Type:</span>{' '}
                <span>{result.molecule_type}</span>
              </div>
            )}
            <div>
              <span className="text-[var(--color-text-secondary)]">Bioactivities:</span>{' '}
              <span className="font-medium">{result.bioactivity_count}</span>
            </div>
          </div>

          {result.bioactivities.length > 0 && (
            <div className="mt-2">
              <span className="text-[var(--color-text-secondary)] text-xs">Top bioactivities:</span>
              <div className="mt-1 space-y-1">
                {result.bioactivities.slice(0, 3).map((act, i) => (
                  <div key={i} className="text-xs bg-[var(--color-surface-elevated)] rounded px-2 py-1 border border-purple-100 dark:border-purple-800">
                    <span className="font-medium">{act.target_name || act.target_chembl_id}</span>
                    {act.activity_value && (
                      <span className="ml-2 text-[var(--color-text-secondary)]">
                        {act.activity_type}: {act.activity_value} {act.activity_unit}
                      </span>
                    )}
                  </div>
                ))}
              </div>
            </div>
          )}

          {result.url && (
            <a
              href={safeHref(result.url)}
              target="_blank"
              rel="noopener noreferrer"
              className="inline-flex items-center gap-1 text-purple-600 hover:text-purple-800 dark:text-purple-400 dark:hover:text-purple-300"
            >
              View on ChEMBL
              <svg className="w-3 h-3" fill="none" stroke="currentColor" viewBox="0 0 24 24">
                <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={2} d="M10 6H6a2 2 0 00-2 2v10a2 2 0 002 2h10a2 2 0 002-2v-4M14 4h6m0 0v6m0-6L10 14" />
              </svg>
            </a>
          )}
        </div>
      )}
    </div>
  );
}

function COCONUTCard({ result }: { result: COCONUTResult | null }) {
  if (!result) {
    return (
      <div className="border border-[var(--color-border)] rounded-lg p-4">
        <div className="flex items-center gap-2">
          <img src={coconutLogo} alt="COCONUT" className="w-5 h-5 rounded-sm object-contain" />
          <span className="font-medium text-[var(--color-text-primary)]">COCONUT</span>
          <span className="text-xs text-[var(--color-text-muted)]">Lookup unavailable</span>
        </div>
      </div>
    );
  }

  return (
    <div className={`border rounded-lg p-4 ${result.found ? 'border-amber-200 bg-amber-50 dark:border-amber-800 dark:bg-amber-900/20' : 'border-[var(--color-border)]'}`}>
      <div className="flex items-center gap-2 mb-2">
        <img src={coconutLogo} alt="COCONUT" className="w-5 h-5 rounded-sm object-contain" />
        <span className="font-medium text-[var(--color-text-primary)]">COCONUT</span>
        {result.found ? (
          <span className="px-2 py-0.5 bg-amber-100 text-amber-800 dark:bg-amber-900/30 dark:text-amber-400 text-xs rounded-full">Natural Product</span>
        ) : (
          <span className="px-2 py-0.5 bg-[var(--color-surface-sunken)] text-[var(--color-text-secondary)] text-xs rounded-full">Not Found</span>
        )}
      </div>

      {result.found && (
        <div className="mt-3 space-y-2 text-sm">
          <div className="grid grid-cols-2 gap-2">
            <div>
              <span className="text-[var(--color-text-secondary)]">COCONUT ID:</span>{' '}
              <span className="font-mono text-amber-700 dark:text-amber-400">{result.coconut_id}</span>
            </div>
            {result.name && (
              <div>
                <span className="text-[var(--color-text-secondary)]">Name:</span>{' '}
                <span className="text-[var(--color-text-primary)]">{result.name}</span>
              </div>
            )}
            {result.organism && (
              <div className="col-span-2">
                <span className="text-[var(--color-text-secondary)]">Organism:</span>{' '}
                <span className="italic text-[var(--color-text-primary)]">{result.organism}</span>
                {result.organism_type && (
                  <span className="text-[var(--color-text-secondary)] ml-1">({result.organism_type})</span>
                )}
              </div>
            )}
            {result.nplikeness != null && (
              <div>
                <span className="text-[var(--color-text-secondary)]">NP-likeness:</span>{' '}
                <span className={result.nplikeness > 0 ? 'text-amber-700 dark:text-amber-400' : 'text-[var(--color-text-secondary)]'}>
                  {result.nplikeness.toFixed(2)}
                </span>
              </div>
            )}
          </div>

          {result.url && (
            <a
              href={safeHref(result.url)}
              target="_blank"
              rel="noopener noreferrer"
              className="inline-flex items-center gap-1 text-amber-600 hover:text-amber-800 dark:text-amber-400 dark:hover:text-amber-300"
            >
              View on COCONUT
              <svg className="w-3 h-3" fill="none" stroke="currentColor" viewBox="0 0 24 24">
                <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={2} d="M10 6H6a2 2 0 00-2 2v10a2 2 0 002 2h10a2 2 0 002-2v-4M14 4h6m0 0v6m0-6L10 14" />
              </svg>
            </a>
          )}
        </div>
      )}

      {!result.found && (
        <p className="text-xs text-[var(--color-text-secondary)] mt-2">
          Not found in the COCONUT natural products database.
        </p>
      )}
    </div>
  );
}

function WikidataCard({ result }: { result: WikidataResult | null }) {
  if (!result) {
    return (
      <div className="border border-[var(--color-border)] rounded-lg p-4">
        <div className="flex items-center gap-2">
          <img src={wikidataLogo} alt="Wikidata" className="w-5 h-5 rounded-sm object-contain" />
          <span className="font-medium text-[var(--color-text-primary)]">Wikidata</span>
          <span className="text-xs text-[var(--color-text-muted)]">Lookup unavailable</span>
        </div>
      </div>
    );
  }

  return (
    <div className={`border rounded-lg p-4 ${result.found ? 'border-teal-200 bg-teal-50 dark:border-teal-800 dark:bg-teal-900/20' : 'border-[var(--color-border)]'}`}>
      <div className="flex items-center gap-2 mb-2">
        <img src={wikidataLogo} alt="Wikidata" className="w-5 h-5 rounded-sm object-contain" />
        <span className="font-medium text-[var(--color-text-primary)]">Wikidata</span>
        {result.found ? (
          <span className="px-2 py-0.5 bg-teal-100 text-teal-800 dark:bg-teal-900/30 dark:text-teal-400 text-xs rounded-full">Found</span>
        ) : (
          <span className="px-2 py-0.5 bg-[var(--color-surface-sunken)] text-[var(--color-text-secondary)] text-xs rounded-full">Not Found</span>
        )}
      </div>

      {result.found && (
        <div className="mt-3 space-y-2 text-sm">
          <div className="grid grid-cols-2 gap-2">
            {result.label && (
              <div className="col-span-2">
                <span className="text-[var(--color-text-secondary)]">Name:</span>{' '}
                <span className="text-[var(--color-text-primary)]">{result.label}</span>
              </div>
            )}
            {result.cas && (
              <div>
                <span className="text-[var(--color-text-secondary)]">CAS:</span>{' '}
                <span className="font-mono text-teal-700 dark:text-teal-400">{result.cas}</span>
              </div>
            )}
            {result.molecular_formula && (
              <div>
                <span className="text-[var(--color-text-secondary)]">Formula:</span>{' '}
                <span className="font-mono">{result.molecular_formula}</span>
              </div>
            )}
            {result.molecular_weight != null && (
              <div>
                <span className="text-[var(--color-text-secondary)]">MW:</span>{' '}
                <span>{result.molecular_weight.toFixed(2)}</span>
              </div>
            )}
          </div>

          {result.url && (
            <a
              href={safeHref(result.url)}
              target="_blank"
              rel="noopener noreferrer"
              className="inline-flex items-center gap-1 text-teal-600 hover:text-teal-800 dark:text-teal-400 dark:hover:text-teal-300"
            >
              View on Wikidata
              <svg className="w-3 h-3" fill="none" stroke="currentColor" viewBox="0 0 24 24">
                <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={2} d="M10 6H6a2 2 0 00-2 2v10a2 2 0 002 2h10a2 2 0 002-2v-4M14 4h6m0 0v6m0-6L10 14" />
              </svg>
            </a>
          )}
        </div>
      )}
    </div>
  );
}

function SureChEMBLCard({ result }: { result: SureChEMBLResult | null }) {
  if (!result) {
    return (
      <div className="border border-[var(--color-border)] rounded-lg p-4">
        <div className="flex items-center gap-2">
          <img src={surechemblLogo} alt="SureChEMBL" className="w-5 h-5 rounded-sm object-contain" />
          <span className="font-medium text-[var(--color-text-primary)]">SureChEMBL</span>
          <span className="text-xs text-[var(--color-text-muted)]">Lookup unavailable</span>
        </div>
      </div>
    );
  }

  return (
    <div className={`border rounded-lg p-4 ${result.found ? 'border-rose-200 bg-rose-50 dark:border-rose-800 dark:bg-rose-900/20' : 'border-[var(--color-border)]'}`}>
      <div className="flex items-center gap-2 mb-2">
        <img src={surechemblLogo} alt="SureChEMBL" className="w-5 h-5 rounded-sm object-contain" />
        <span className="font-medium text-[var(--color-text-primary)]">SureChEMBL</span>
        {result.found ? (
          <span className="px-2 py-0.5 bg-rose-100 text-rose-800 dark:bg-rose-900/30 dark:text-rose-400 text-xs rounded-full">In Patents</span>
        ) : (
          <span className="px-2 py-0.5 bg-[var(--color-surface-sunken)] text-[var(--color-text-secondary)] text-xs rounded-full">Not Found</span>
        )}
      </div>

      {result.found && (
        <div className="mt-3 space-y-2 text-sm">
          <div className="grid grid-cols-2 gap-2">
            {result.schembl_id && (
              <div>
                <span className="text-[var(--color-text-secondary)]">SCHEMBL ID:</span>{' '}
                <span className="font-mono text-rose-700 dark:text-rose-400">{result.schembl_id}</span>
              </div>
            )}
            {result.patent_count != null && (
              <div>
                <span className="text-[var(--color-text-secondary)]">Patents:</span>{' '}
                <span className="font-medium">{result.patent_count.toLocaleString()}</span>
              </div>
            )}
            {result.source && (
              <div>
                <span className="text-[var(--color-text-secondary)]">Source:</span>{' '}
                <span className="text-[var(--color-text-primary)]">{result.source === 'surechembl_api' ? 'SureChEMBL API' : result.source ?? 'UniChem'}</span>
              </div>
            )}
          </div>

          {result.url && (
            <a
              href={safeHref(result.url)}
              target="_blank"
              rel="noopener noreferrer"
              className="inline-flex items-center gap-1 text-rose-600 hover:text-rose-800 dark:text-rose-400 dark:hover:text-rose-300"
            >
              View on SureChEMBL
              <svg className="w-3 h-3" fill="none" stroke="currentColor" viewBox="0 0 24 24">
                <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={2} d="M10 6H6a2 2 0 00-2 2v10a2 2 0 002 2h10a2 2 0 002-2v-4M14 4h6m0 0v6m0-6L10 14" />
              </svg>
            </a>
          )}
        </div>
      )}

      {!result.found && (
        <p className="text-xs text-[var(--color-text-secondary)] mt-2">
          Not found in patent literature via SureChEMBL.
        </p>
      )}
    </div>
  );
}
