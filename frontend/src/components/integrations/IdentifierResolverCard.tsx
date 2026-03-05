import { safeHref } from '../../lib/sanitize';
import type { ResolvedCompound } from '../../types/integrations';

interface IdentifierResolverCardProps {
  result: ResolvedCompound;
}

const TYPE_LABELS: Record<string, string> = {
  smiles: 'SMILES',
  inchi: 'InChI',
  inchikey: 'InChIKey',
  pubchem_cid: 'PubChem CID',
  chembl_id: 'ChEMBL ID',
  cas: 'CAS Number',
  drugbank_id: 'DrugBank ID',
  chebi_id: 'ChEBI ID',
  unii: 'UNII',
  wikipedia: 'Wikipedia',
  name: 'Compound Name',
};

const CONFIDENCE_COLORS: Record<string, string> = {
  high: 'bg-green-100 text-green-800',
  medium: 'bg-yellow-100 text-yellow-800',
  low: 'bg-red-100 text-red-800',
};

export function IdentifierResolverCard({ result }: IdentifierResolverCardProps) {
  if (!result.resolved) return null;

  const xrefs = result.cross_references;
  const dbLinks: { label: string; id: string; url: string }[] = [];

  if (xrefs.pubchem_cid) {
    dbLinks.push({
      label: 'PubChem',
      id: `CID:${xrefs.pubchem_cid}`,
      url: `https://pubchem.ncbi.nlm.nih.gov/compound/${xrefs.pubchem_cid}`,
    });
  }
  if (xrefs.chembl_id) {
    dbLinks.push({
      label: 'ChEMBL',
      id: xrefs.chembl_id,
      url: `https://www.ebi.ac.uk/chembl/compound_report_card/${xrefs.chembl_id}`,
    });
  }
  if (xrefs.drugbank_id) {
    dbLinks.push({
      label: 'DrugBank',
      id: xrefs.drugbank_id,
      url: `https://go.drugbank.com/drugs/${xrefs.drugbank_id}`,
    });
  }
  if (xrefs.chebi_id) {
    dbLinks.push({
      label: 'ChEBI',
      id: xrefs.chebi_id,
      url: `https://www.ebi.ac.uk/chebi/searchId.do?chebiId=${xrefs.chebi_id}`,
    });
  }
  if (xrefs.kegg_id) {
    dbLinks.push({
      label: 'KEGG',
      id: xrefs.kegg_id,
      url: `https://www.genome.jp/entry/${xrefs.kegg_id}`,
    });
  }

  return (
    <div className="bg-indigo-50 border border-indigo-200 rounded-lg p-4 mb-4">
      <div className="flex items-center gap-2 mb-3">
        <span className="px-2 py-0.5 bg-indigo-100 text-indigo-800 text-xs font-medium rounded-full">
          {TYPE_LABELS[result.identifier_type_detected] || result.identifier_type_detected}
        </span>
        <span className={`px-2 py-0.5 text-xs font-medium rounded-full ${CONFIDENCE_COLORS[result.confidence] || 'bg-gray-100 text-gray-600'}`}>
          {result.confidence} confidence
        </span>
        <span className="text-xs text-gray-500 ml-auto">
          via {result.resolution_source}
        </span>
      </div>

      {/* Resolution chain */}
      {result.resolution_chain.length > 0 && (
        <div className="text-xs text-gray-600 mb-3 font-mono bg-white rounded px-2 py-1">
          {result.resolution_chain.join(' → ')}
        </div>
      )}

      {/* Cross-references */}
      {dbLinks.length > 0 && (
        <div className="mt-2">
          <span className="text-xs text-gray-500 font-medium">Cross-references:</span>
          <div className="flex flex-wrap gap-2 mt-1">
            {dbLinks.map((link) => (
              <a
                key={link.label}
                href={safeHref(link.url)}
                target="_blank"
                rel="noopener noreferrer"
                className="inline-flex items-center gap-1 px-2 py-0.5 bg-white border border-gray-200 rounded text-xs text-blue-600 hover:text-blue-800 hover:border-blue-300 transition-colors"
              >
                <span className="font-medium text-gray-700">{link.label}:</span>
                <span className="font-mono">{link.id}</span>
                <svg className="w-2.5 h-2.5" fill="none" stroke="currentColor" viewBox="0 0 24 24">
                  <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={2} d="M10 6H6a2 2 0 00-2 2v10a2 2 0 002 2h10a2 2 0 002-2v-4M14 4h6m0 0v6m0-6L10 14" />
                </svg>
              </a>
            ))}
          </div>
        </div>
      )}
    </div>
  );
}
