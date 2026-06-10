/**
 * Routes that never render molecular structures and therefore must not
 * pay the RDKit WASM init cost (or the splash wait) on first load.
 *
 * History and Bookmarks are deliberately NOT listed: both render
 * MoleculeViewer thumbnails.
 */
const RDKIT_FREE_ROUTES = new Set(['/about', '/privacy', '/diagnostics']);

export function routeNeedsRDKit(pathname: string): boolean {
  const normalized = pathname.replace(/\/+$/, '') || '/';
  return !RDKIT_FREE_ROUTES.has(normalized);
}
