"""
FilterCatalog Singleton for Structural Alert Screening

Provides cached access to RDKit FilterCatalog for PAINS, BRENK, and other
structural alert pattern sets. Uses lru_cache for singleton behavior.

PAINS (Pan-Assay INterference compoundS):
- PAINS_A: 16 patterns (most severe)
- PAINS_B: 55 patterns
- PAINS_C: 409 patterns
- Total: 480 patterns

BRENK: 105 patterns for known problematic functional groups

Note: Structural alerts are warnings, not automatic rejections.
87 FDA-approved drugs contain PAINS patterns.
"""

from dataclasses import dataclass
from functools import lru_cache
from typing import Dict

from rdkit.Chem.FilterCatalog import FilterCatalog, FilterCatalogParams

# Available catalog types and their descriptions
AVAILABLE_CATALOGS: Dict[str, Dict[str, str]] = {
    "PAINS": {
        "name": "PAINS (Pan-Assay Interference Compounds)",
        "description": "Substructure filters for removal of frequent hitters in high-throughput screens",
        "scope": "Identifies compounds that appear active in multiple assay types due to non-specific mechanisms such as aggregation, redox cycling, fluorescence interference, or chemical reactivity rather than genuine target binding",
        "reference": "Baell JB, Holloway GA. New Substructure Filters for Removal of Pan Assay Interference Compounds (PAINS) from Screening Libraries and for Their Exclusion in Bioassays. J Med Chem 53 (2010) 2719-2740.",
        "doi": "10.1021/jm901137j",
        "pmid": "20131845",
        "pattern_count": "480",
        "severity": "warning",
        "note": "87 FDA-approved drugs contain PAINS patterns — these are warnings, not rejections",
    },
    "PAINS_A": {
        "name": "PAINS Class A (Most Frequent Hitters)",
        "description": "16 most promiscuous PAINS patterns — highest hit rates across multiple assay technologies",
        "scope": "PAINS filters (family A)",
        "reference": "Baell JB, Holloway GA. J Med Chem 53 (2010) 2719-2740.",
        "doi": "10.1021/jm901137j",
        "pmid": "20131845",
        "pattern_count": "16",
        "severity": "warning",
    },
    "PAINS_B": {
        "name": "PAINS Class B (Moderate Hitters)",
        "description": "55 moderately promiscuous PAINS patterns",
        "scope": "PAINS filters (family B)",
        "reference": "Baell JB, Holloway GA. J Med Chem 53 (2010) 2719-2740.",
        "doi": "10.1021/jm901137j",
        "pmid": "20131845",
        "pattern_count": "55",
        "severity": "warning",
    },
    "PAINS_C": {
        "name": "PAINS Class C (Lower Frequency Hitters)",
        "description": "409 less frequent but still notable PAINS patterns",
        "scope": "PAINS filters (family C)",
        "reference": "Baell JB, Holloway GA. J Med Chem 53 (2010) 2719-2740.",
        "doi": "10.1021/jm901137j",
        "pmid": "20131845",
        "pattern_count": "409",
        "severity": "warning",
    },
    "BRENK": {
        "name": "Brenk Structural Alerts",
        "description": "Filters for potentially problematic functional groups identified during assembly of screening libraries for neglected tropical disease drug discovery",
        "scope": "Unwanted functionality due to potential toxicity reasons or unfavourable pharmacokinetic properties",
        "reference": "Brenk R, Schipani A, James D, et al. Lessons Learnt from Assembling Screening Libraries for Drug Discovery for Neglected Diseases. ChemMedChem 3 (2008) 435-444.",
        "doi": "10.1002/cmdc.200700139",
        "pmid": "18064617",
        "pattern_count": "105",
        "severity": "warning",
    },
    "NIH": {
        "name": "NIH/NCGC Alerts",
        "description": "Structural alert patterns from the NIH Chemical Genomics Center to annotate compounds with problematic functional groups",
        "scope": "Annotate compounds with problematic functional groups causing aggregation, autofluorescence, and reactivity artifacts",
        "reference": "Jadhav A, Ferreira RS, Klumpp C, et al. Quantitative Analyses of Aggregation, Autofluorescence, and Reactivity Artifacts in a Screen for Inhibitors of a Thiol Protease. J Med Chem 53 (2010) 37-51.",
        "doi": "10.1021/jm901070c",
        "pmid": "19908840",
        "pattern_count": "180",
        "severity": "warning",
    },
    "ZINC": {
        "name": "ZINC Drug-likeness Filters",
        "description": "Property and functional group filters from the ZINC database used to ensure drug-likeness of screening compounds",
        "scope": "Drug-likeness and unwanted functional group filters",
        "reference": "Irwin JJ, Shoichet BK. ZINC — A Free Database of Commercially Available Compounds for Virtual Screening. J Chem Inf Model 45 (2005) 177-182.",
        "doi": "10.1021/ci049714+",
        "pmid": "15667143",
        "pattern_count": "50",
        "severity": "warning",
    },
    "CHEMBL_BMS": {
        "name": "Bristol-Myers Squibb HTS Deck Filters",
        "description": "Empirically derived filters to remove promiscuous compounds from BMS high-throughput screening decks",
        "scope": "Identifies reactive and promiscuous functional groups that cause high false-positive rates in HTS",
        "reference": "Pearce BC, Sofia MJ, Good AC, et al. An Empirical Process for the Design of High-Throughput Screening Deck Filters. J Chem Inf Model 46 (2006) 1060-1068.",
        "doi": "10.1021/ci050504m",
        "pmid": "16711725",
        "pattern_count": "180",
        "severity": "warning",
    },
    "CHEMBL_DUNDEE": {
        "name": "University of Dundee NTD Screening Filters",
        "description": "Compound filters developed at the University of Dundee for neglected tropical disease screening library curation",
        "scope": "Screens for unwanted functionality due to potential toxicity or unfavourable pharmacokinetic properties",
        "reference": "Brenk R, Schipani A, James D, et al. Lessons Learnt from Assembling Screening Libraries for Drug Discovery for Neglected Diseases. ChemMedChem 3 (2008) 435-444.",
        "doi": "10.1002/cmdc.200700139",
        "pmid": "18064617",
        "pattern_count": "105",
        "severity": "warning",
    },
    "CHEMBL_GLAXO": {
        "name": "Glaxo Wellcome Hard Filters",
        "description": "Compound quality filters developed at GlaxoSmithKline to remove reactive and undesirable compounds before HTS pooling",
        "scope": "Identifies reactive, unstable, or otherwise undesirable compounds to exclude from screening pools",
        "reference": "Hann M, Hudson B, Lewell X, et al. Strategic Pooling of Compounds for High-Throughput Screening. J Chem Inf Comput Sci 39 (1999) 897-902.",
        "doi": "10.1021/ci990423o",
        "pmid": "10529988",
        "pattern_count": "55",
        "severity": "warning",
    },
    "CHEMBL_INPHARMATICA": {
        "name": "Inpharmatica Unwanted Fragment Filters",
        "description": "Filters developed at Inpharmatica Ltd to identify unwanted fragments including mutagenic groups, groups with unfavourable pharmacokinetic properties, and reactive functionalities",
        "scope": "Identifies unwanted fragments: mutagenic groups (e.g. nitro), poor PK groups (e.g. sulfates), and reactive groups (e.g. thiols, 2-halopyridines)",
        "reference": "Inpharmatica Ltd. (no peer-reviewed publication; included in ChEMBL structural alerts since ChEMBL 23)",
        "pattern_count": "91",
        "severity": "warning",
    },
    "CHEMBL_LINT": {
        "name": "Lilly MedChem Rules (LINT)",
        "description": "275 rules developed over 18 years at Eli Lilly to identify reactive or promiscuous compounds that interfere with biological assays",
        "scope": "Identifies reactive compounds, assay interference (fluorescence/absorbance), protein-damaging agents, unstable compounds, and non-drug-like molecules",
        "reference": "Bruns RF, Watson IA. Rules for Identifying Potentially Reactive or Promiscuous Compounds. J Med Chem 55 (2012) 9763-9772.",
        "doi": "10.1021/jm301008n",
        "pmid": "23061697",
        "pattern_count": "57",
        "severity": "warning",
    },
    "CHEMBL_MLSMR": {
        "name": "NIH MLSMR Excluded Functionality Filters",
        "description": "Structural filters from the NIH Molecular Libraries Small Molecule Repository to exclude compounds with problematic functionalities",
        "scope": "Excludes compounds with reactive, toxic, or otherwise undesirable functional groups from the NIH screening collection",
        "reference": "NIH Molecular Libraries Program. MLSMR Excluded Functionality Filters (included in ChEMBL structural alerts since ChEMBL 23).",
        "pattern_count": "116",
        "severity": "warning",
    },
    "CHEMBL_SURECHEMBL": {
        "name": "SureChEMBL Patent-derived Alerts",
        "description": "Structural alerts used in the SureChEMBL patent chemistry database for compound quality assessment",
        "scope": "Identifies known problematic substructures in patent-extracted chemical entities",
        "reference": "Papadatos G, Davies M, Dedman N, et al. SureChEMBL: a large-scale, chemically annotated patent document database. Nucleic Acids Res 44 (2016) D1220-D1228.",
        "doi": "10.1093/nar/gkv1253",
        "pmid": "26582922",
        "pattern_count": "166",
        "severity": "warning",
    },
    "ALL": {
        "name": "All Catalogs Combined",
        "description": "Screen against all available structural alert patterns",
        "scope": "Combined screening across all filter sets",
        "pattern_count": "variable",
        "severity": "warning",
    },
}


@dataclass
class CatalogInfo:
    """Information about a loaded FilterCatalog."""

    catalog: FilterCatalog
    catalog_type: str
    num_entries: int


def _build_catalog_params(catalog_type: str) -> FilterCatalogParams:
    """
    Build FilterCatalogParams for the specified catalog type.

    Args:
        catalog_type: One of PAINS, PAINS_A, PAINS_B, PAINS_C, BRENK, NIH, ZINC, ALL

    Returns:
        Configured FilterCatalogParams
    """
    params = FilterCatalogParams()
    catalog_type_upper = catalog_type.upper()

    if catalog_type_upper == "PAINS":
        # All PAINS classes (A + B + C)
        params.AddCatalog(FilterCatalogParams.FilterCatalogs.PAINS_A)
        params.AddCatalog(FilterCatalogParams.FilterCatalogs.PAINS_B)
        params.AddCatalog(FilterCatalogParams.FilterCatalogs.PAINS_C)
    elif catalog_type_upper == "PAINS_A":
        params.AddCatalog(FilterCatalogParams.FilterCatalogs.PAINS_A)
    elif catalog_type_upper == "PAINS_B":
        params.AddCatalog(FilterCatalogParams.FilterCatalogs.PAINS_B)
    elif catalog_type_upper == "PAINS_C":
        params.AddCatalog(FilterCatalogParams.FilterCatalogs.PAINS_C)
    elif catalog_type_upper == "BRENK":
        params.AddCatalog(FilterCatalogParams.FilterCatalogs.BRENK)
    elif catalog_type_upper == "NIH":
        params.AddCatalog(FilterCatalogParams.FilterCatalogs.NIH)
    elif catalog_type_upper == "ZINC":
        params.AddCatalog(FilterCatalogParams.FilterCatalogs.ZINC)
    elif catalog_type_upper == "CHEMBL_BMS":
        params.AddCatalog(FilterCatalogParams.FilterCatalogs.CHEMBL_BMS)
    elif catalog_type_upper == "CHEMBL_DUNDEE":
        params.AddCatalog(FilterCatalogParams.FilterCatalogs.CHEMBL_Dundee)
    elif catalog_type_upper == "CHEMBL_GLAXO":
        params.AddCatalog(FilterCatalogParams.FilterCatalogs.CHEMBL_Glaxo)
    elif catalog_type_upper == "CHEMBL_INPHARMATICA":
        params.AddCatalog(FilterCatalogParams.FilterCatalogs.CHEMBL_Inpharmatica)
    elif catalog_type_upper == "CHEMBL_LINT":
        params.AddCatalog(FilterCatalogParams.FilterCatalogs.CHEMBL_LINT)
    elif catalog_type_upper == "CHEMBL_MLSMR":
        params.AddCatalog(FilterCatalogParams.FilterCatalogs.CHEMBL_MLSMR)
    elif catalog_type_upper == "CHEMBL_SURECHEMBL":
        params.AddCatalog(FilterCatalogParams.FilterCatalogs.CHEMBL_SureChEMBL)
    elif catalog_type_upper == "ALL":
        # Add all available catalogs
        params.AddCatalog(FilterCatalogParams.FilterCatalogs.PAINS_A)
        params.AddCatalog(FilterCatalogParams.FilterCatalogs.PAINS_B)
        params.AddCatalog(FilterCatalogParams.FilterCatalogs.PAINS_C)
        params.AddCatalog(FilterCatalogParams.FilterCatalogs.BRENK)
        params.AddCatalog(FilterCatalogParams.FilterCatalogs.NIH)
        params.AddCatalog(FilterCatalogParams.FilterCatalogs.ZINC)
        params.AddCatalog(FilterCatalogParams.FilterCatalogs.CHEMBL_BMS)
        params.AddCatalog(FilterCatalogParams.FilterCatalogs.CHEMBL_Dundee)
        params.AddCatalog(FilterCatalogParams.FilterCatalogs.CHEMBL_Glaxo)
        params.AddCatalog(FilterCatalogParams.FilterCatalogs.CHEMBL_Inpharmatica)
        params.AddCatalog(FilterCatalogParams.FilterCatalogs.CHEMBL_LINT)
        params.AddCatalog(FilterCatalogParams.FilterCatalogs.CHEMBL_MLSMR)
        params.AddCatalog(FilterCatalogParams.FilterCatalogs.CHEMBL_SureChEMBL)
    else:
        raise ValueError(
            f"Unknown catalog type: {catalog_type}. Available: {list(AVAILABLE_CATALOGS.keys())}"
        )

    return params


@lru_cache(maxsize=16)
def get_filter_catalog(catalog_type: str = "PAINS") -> CatalogInfo:
    """
    Get a cached FilterCatalog for the specified catalog type.

    Uses lru_cache for singleton-like behavior per catalog type.
    First call initializes the catalog; subsequent calls return cached instance.

    Args:
        catalog_type: Catalog type (PAINS, PAINS_A, PAINS_B, PAINS_C, BRENK, NIH, ZINC, ALL)

    Returns:
        CatalogInfo with initialized FilterCatalog

    Raises:
        ValueError: If catalog_type is unknown
    """
    catalog_type_upper = catalog_type.upper()

    if catalog_type_upper not in AVAILABLE_CATALOGS:
        raise ValueError(
            f"Unknown catalog type: {catalog_type}. Available: {list(AVAILABLE_CATALOGS.keys())}"
        )

    params = _build_catalog_params(catalog_type_upper)
    catalog = FilterCatalog(params)

    return CatalogInfo(
        catalog=catalog,
        catalog_type=catalog_type_upper,
        num_entries=catalog.GetNumEntries(),
    )


def list_available_catalogs() -> Dict[str, Dict[str, str]]:
    """
    List all available filter catalogs with their descriptions.

    Returns:
        Dictionary mapping catalog type to description info
    """
    return AVAILABLE_CATALOGS.copy()


def get_catalog_stats() -> Dict[str, int]:
    """
    Get pattern counts for each catalog type.

    Note: This loads each catalog which may be slow on first call.

    Returns:
        Dictionary mapping catalog type to pattern count
    """
    stats = {}
    for catalog_type in AVAILABLE_CATALOGS:
        if catalog_type != "ALL":  # Skip ALL to avoid counting duplicates
            info = get_filter_catalog(catalog_type)
            stats[catalog_type] = info.num_entries
    return stats
