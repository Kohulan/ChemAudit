"""
Pattern Category Classifier for Structural Alerts

Maps structural alert pattern names to broad concern categories
using keyword-based rules. Categories help chemists quickly assess
what type of concern an alert represents.

Categories:
- Reactive Group: Electrophilic/nucleophilic reactivity
- Metabolic Liability: Metabolically labile groups
- Toxicophore: Known toxic substructures
- Assay Interference: HTS assay artifacts (PAINS default)
- Physicochemical: Property-based filters
- Unwanted Functionality: General undesirable groups
"""

from typing import List, Tuple

# Ordered list of (keywords, category) — first match wins.
# Keywords are matched case-insensitively against the pattern name.
_CATEGORY_RULES: List[Tuple[List[str], str]] = [
    # --- Reactive Groups ---
    (
        [
            "michael acceptor",
            "acyl halide",
            "acyl_halide",
            "acid halide",
            "acid_halide",
            "acyl chloride",
            "acyl fluoride",
            "acyl cyanide",
            "acyl_cyanide",
            "epoxide",
            "aziridine",
            "isocyanate",
            "isothiocyanate",
            "ketene",
            "acyl_phosph",
            "sulfonyl_halide",
            "sulfonyl halide",
            "allyl_halide",
            "alkyl halide",
            "alkyl_halide",
            "vinyl_halide",
            "activated_acetylene",
            "activated_vinyl",
            "activated_4mem",
            "activated_diazo",
            "activated_S#O",
            "beta_lactam",
            "beta-keto",
            "anhydride",
            "peroxide",
            "acetal",
            "halo_ether",
            "halo ether",
            "N-halo",
            "N_halo",
            "triflate",
            "pentafluorophenyl ester",
            "phosphorane",
            "phosphor",
            "silicon halogen",
            "reactive",
            "Reactive",
        ],
        "Reactive Group",
    ),
    # --- Toxicophores ---
    (
        [
            "nitro",
            "nitroso",
            "N-nitroso",
            "N_nitroso",
            "azide",
            "azido",
            "diazo",
            "hydrazine",
            "hydrazide",
            "acyl hydrazine",
            "acyl_hydrazine",
            "nitrogen_mustard",
            "nitrogen mustard",
            "aniline",
            "diaminobenzene",
            "heavy metal",
            "benzidine",
            "polycyclic aromatic",
            "cyanamide",
            "cyanohydrin",
            "hydroxamic acid",
            "hydroxamic_acid",
        ],
        "Toxicophore",
    ),
    # --- Metabolic Liability ---
    (
        [
            "isolated alkene",
            "isolated_alkene",
            "acyclic C=C",
            "acyclic_C=C",
            "stilbene",
            "polyene",
            "aldehyde",
            "thiol",
            "disulphide",
            "disulfide",
            "enamine",
            "imine",
            "oxime",
            "thioester",
            "thiocarbonyl",
            "thio_ester",
            "ester_of_HOBT",
            "ester of HOBT",
            "phenol ester",
            "phenol_ester",
            "phenyl carbonate",
            "phenyl_carbonate",
            "four member lactone",
            "four_member_lactone",
            "catechol",
            "hydroquinone",
            "quinone",
            "chinone",
            "N oxide",
            "N_oxide",
            "crown ether",
            "crown_ether",
            "Oxygen-nitrogen single bond",
            "oxygen_nitrogen",
            "sulfinic acid",
            "sulfonic acid",
            "Sulfonic_acid",
            "sulphate",
            "sulfur oxygen",
            "sulfur_oxygen",
            "Sulphur-nitrogen",
            "sulphur_nitrogen",
            "triple bond",
            "triple_bond",
        ],
        "Metabolic Liability",
    ),
    # --- Physicochemical ---
    (
        [
            "Non-Hydrogen_atoms",
            "Non-Hydrogen atoms",
            "carbons",
            "N,O,S",
            "halogens",
            "rings",
            "rotatable",
            "chiral",
            "molecular_weight",
            "Molecular weight",
            "Aliphatic long chain",
            "aliphatic_long_chain",
            "perfluorinated",
        ],
        "Physicochemical",
    ),
    # --- Assay Interference (PAINS-like) ---
    (
        [
            "azo",
            "Azo",
            "rhodanine",
            "thiazole_amine",
            "thiaz_ene",
            "pyrrole_A",
            "catechol_A",
            "ene_six_het",
            "hzone",
            "anil_di_alk",
            "mannich",
            "het_5_pyrazole",
            "quinone_A",
            "imine_one",
            "keto_phenone",
            "azo_A",
            "sulfonamide_A",
            "indol_3yl_alk",
            "phenol_sulfite",
            "biotin",
        ],
        "Assay Interference",
    ),
]

# Catalog-level defaults: if catalog is PAINS*, default to Assay Interference
_CATALOG_DEFAULTS: dict[str, str] = {
    "PAINS": "Assay Interference",
    "PAINS_A": "Assay Interference",
    "PAINS_B": "Assay Interference",
    "PAINS_C": "Assay Interference",
}


def classify_pattern(pattern_name: str, catalog_source: str) -> str:
    """
    Classify a structural alert pattern into a concern category.

    Args:
        pattern_name: The pattern name from RDKit FilterCatalog
        catalog_source: The catalog source (PAINS, BRENK, CHEMBL_DUNDEE, etc.)

    Returns:
        Category string: one of "Reactive Group", "Metabolic Liability",
        "Toxicophore", "Assay Interference", "Physicochemical",
        or "Unwanted Functionality" (default)
    """
    name_lower = pattern_name.lower()

    for keywords, category in _CATEGORY_RULES:
        for kw in keywords:
            if kw.lower() in name_lower:
                return category

    # Catalog-level default
    catalog_upper = catalog_source.upper()
    for prefix, default_cat in _CATALOG_DEFAULTS.items():
        if catalog_upper.startswith(prefix):
            return default_cat

    return "Unwanted Functionality"
