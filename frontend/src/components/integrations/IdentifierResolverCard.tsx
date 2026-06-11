import { ExternalLink } from 'lucide-react';
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
  high: 'bg-status-success-light dark:bg-status-success/15 text-status-success-dark dark:text-status-success',
  medium: 'bg-yellow-100 text-yellow-800 dark:bg-yellow-900/30 dark:text-yellow-400',
  low: 'bg-red-100 text-red-800 dark:bg-red-900/30 dark:text-red-400',
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
    <div className="bg-[var(--color-surface-sunken)] border border-[var(--color-border)] rounded-lg p-4 mb-4">
      <div className="flex items-center gap-2 mb-3">
        <span className="px-2 py-0.5 bg-chem-primary-100 text-chem-primary-800 dark:bg-chem-primary-900/30 dark:text-chem-primary-300 text-xs font-medium rounded-full">
          {TYPE_LABELS[result.identifier_type_detected] || result.identifier_type_detected}
        </span>
        <span className={`px-2 py-0.5 text-xs font-medium rounded-full ${CONFIDENCE_COLORS[result.confidence] || 'bg-[var(--color-surface-sunken)] text-[var(--color-text-secondary)]'}`}>
          {result.confidence} confidence
        </span>
        <span className="text-xs text-[var(--color-text-secondary)] ml-auto">
          via {result.resolution_source}
        </span>
      </div>

      {/* Resolution chain */}
      {result.resolution_chain.length > 0 && (
        <div className="text-xs text-[var(--color-text-secondary)] mb-3 font-mono bg-[var(--color-surface-elevated)] rounded px-2 py-1">
          {result.resolution_chain.join(' → ')}
        </div>
      )}

      {/* Cross-references */}
      {dbLinks.length > 0 && (
        <div className="mt-2">
          <span className="text-xs text-[var(--color-text-secondary)] font-medium">Cross-references:</span>
          <div className="flex flex-wrap gap-2 mt-1">
            {dbLinks.map((link) => (
              <a
                key={link.label}
                href={safeHref(link.url)}
                target="_blank"
                rel="noopener noreferrer"
                className="inline-flex items-center gap-1 px-2 py-0.5 bg-[var(--color-surface-elevated)] border border-[var(--color-border)] rounded text-xs text-blue-600 hover:text-blue-800 hover:border-blue-300 transition-colors"
              >
                <span className="font-medium text-[var(--color-text-primary)]">{link.label}:</span>
                <span className="font-mono">{link.id}</span>
                <ExternalLink className="w-2.5 h-2.5" aria-hidden="true" />
              </a>
            ))}
          </div>
        </div>
      )}
    </div>
  );
}
