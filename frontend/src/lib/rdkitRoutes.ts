/**
 * Routes that must not BLOCK first paint on RDKit WASM init (they skip the
 * app-level splash gate). This is not "never renders structures": /diagnostics
 * mounts CrossPipelinePanel, which renders MoleculeViewer. Structures on these
 * routes still work because getRDKit() is a lazy singleton that initializes on
 * demand when a MoleculeViewer mounts — these routes merely defer that cost
 * instead of paying it up front.
 *
 * History and Bookmarks are deliberately NOT listed: both render
 * MoleculeViewer thumbnails above the fold, so eager init is worth it there.
 */
const RDKIT_FREE_ROUTES = new Set(['/about', '/privacy', '/diagnostics']);

export function routeNeedsRDKit(pathname: string): boolean {
  const normalized = pathname.replace(/\/+$/, '') || '/';
  return !RDKIT_FREE_ROUTES.has(normalized);
}
