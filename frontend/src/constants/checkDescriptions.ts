/**
 * Human-readable descriptions for all validation checks.
 *
 * Used by the treemap tooltip, DeepCheckCard, and any UI that needs
 * to explain what a check does. Keys are the snake_case check names
 * returned by the backend.
 */
export const CHECK_DESCRIPTIONS: Record<string, string> = {
  // Basic checks
  parsability: 'Verifies the molecule can be parsed from the input SMILES or other format.',
  sanitization: 'Verifies the molecule passes RDKit chemical sanitization rules.',
  valence: 'Checks that all atoms have chemically valid valences.',
  aromaticity: 'Verifies aromatic rings can be correctly kekulized (alternating single/double bonds).',
  connectivity: 'Detects disconnected molecular fragments (e.g. salts, mixtures encoded as dot-separated SMILES).',

  // Stereochemistry checks
  undefined_stereocenters: 'Detects chiral centers (sp3 atoms with 4 different substituents) where R/S configuration is not specified.',
  undefined_doublebond_stereo: 'Detects double bonds where E/Z (cis/trans) geometry is not specified.',
  conflicting_stereo: 'Detects contradictory stereochemistry annotations on the same atom or bond.',

  // Representation checks
  smiles_roundtrip: 'Converts SMILES to molecule and back; flags if the output differs from the input.',
  inchi_generation: 'Checks whether a valid InChI identifier can be generated for the molecule.',
  inchi_roundtrip: 'Converts to InChI and back; flags if the reconstructed molecule differs.',

  // Stereo & Tautomer (deep)
  stereoisomer_enumeration: 'Enumerates possible stereoisomers for undefined stereocenters (up to 128). Helps identify ambiguous chirality.',
  tautomer_detection: 'Detects tautomeric forms and identifies the canonical tautomer. Reports whether the input matches the canonical form.',
  aromatic_system_validation: 'Checks for unusual aromatic ring sizes (not 5- or 6-membered) and charged aromatic atoms that may indicate errors.',
  coordinate_dimension: 'Reports whether the molecule has 2D coordinates, 3D coordinates, or no coordinate information.',

  // Chemical Composition (deep)
  mixture_detection: 'Detects multi-fragment inputs and classifies each component as drug, salt, solvent, or unknown.',
  solvent_contamination: 'Screens for common lab solvents (water, DMSO, DMF, methanol, etc.) that may contaminate the structure.',
  inorganic_filter: 'Flags molecules lacking carbon (inorganic) or containing metal atoms (organometallic).',
  radical_detection: 'Identifies atoms with unpaired radical electrons indicating unstable or reactive species.',
  isotope_label_detection: 'Detects isotope-labeled atoms (deuterium, carbon-13, etc.) often used in pharmacokinetic studies.',
  trivial_molecule: 'Flags molecules with 3 or fewer heavy atoms as too small for meaningful chemical validation.',

  // Structural Complexity (deep)
  hypervalent_atoms: 'Detects atoms exceeding their normal valence limits, which may indicate unusual bonding or input errors.',
  polymer_detection: 'Identifies possible polymers via SGroup markers, molecular weight above 1500 Da, or dummy atom attachment points.',
  ring_strain: 'Flags 3-membered (cyclopropane) and 4-membered (cyclobutane) rings that have significant ring strain.',
  macrocycle_detection: 'Identifies macrocyclic rings with more than 12 atoms, common in natural products and cyclic peptides.',
  charged_species: 'Reports formal charges and identifies zwitterions (net charge zero with both positive and negative atoms).',
  explicit_hydrogen_audit: 'Audits explicit hydrogen specifications and detects H atom objects from AddHs() processing.',
};

/** Format a snake_case check name to a readable title. */
export function formatCheckName(name: string): string {
  return name.replace(/_/g, ' ').replace(/\b\w/g, (c) => c.toUpperCase());
}
