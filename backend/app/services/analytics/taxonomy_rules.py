"""
Chemical Taxonomy SMARTS Rules

Curated set of ~50 chemotype SMARTS rules for ClassyFire-style chemical
taxonomy classification. Rules are organized into three groups:

- Ring Systems (~20 rules)
- Functional Groups (~15 rules)
- Pharmacophoric / Drug-class (~15 rules)

Each rule has a ``name``, ``smarts``, ``category``, and ``description``.
Rules are compiled to RDKit mol objects via ``get_compiled_rules()`` which
returns a module-level singleton (same pattern as custom_smarts.py).
"""

from __future__ import annotations

import logging
from typing import Optional

from rdkit import Chem

logger = logging.getLogger(__name__)

# ---------------------------------------------------------------------------
# Chemotype Rules (~50 entries)
# ---------------------------------------------------------------------------

CHEMOTYPE_RULES: list[dict] = [
    # ===== Ring Systems (~20) =====
    {
        "name": "beta_lactam",
        "smarts": "[#7]1~[#6](~[#8])~[#6]~[#6]1",
        "category": "Lactams",
        "description": "Beta-lactam ring system (penicillins, cephalosporins)",
    },
    {
        "name": "benzodiazepine",
        "smarts": "[#6]1=[#7][#6][#6](=[#8])[#7]c2ccccc21",
        "category": "Benzodiazepines",
        "description": "Benzodiazepine fused ring system (anxiolytics, sedatives)",
    },
    {
        "name": "indole",
        "smarts": "c1ccc2[nH]ccc2c1",
        "category": "Indoles",
        "description": "Indole heterocycle (tryptophan derivatives, serotonin)",
    },
    {
        "name": "quinoline",
        "smarts": "c1ccc2ncccc2c1",
        "category": "Quinolines",
        "description": "Quinoline heterocycle (antimalarials, fluoroquinolones)",
    },
    {
        "name": "isoquinoline",
        "smarts": "c1ccc2cnccc2c1",
        "category": "Isoquinolines",
        "description": "Isoquinoline heterocycle (opioids, papaverine)",
    },
    {
        "name": "pyridine",
        "smarts": "c1ccncc1",
        "category": "Pyridines",
        "description": "Pyridine ring (niacin, isoniazid)",
    },
    {
        "name": "pyrimidine",
        "smarts": "c1ncncc1",
        "category": "Pyrimidines",
        "description": "Pyrimidine ring (nucleic acid bases, trimethoprim)",
    },
    {
        "name": "imidazole",
        "smarts": "c1cnc[nH]1",
        "category": "Imidazoles",
        "description": "Imidazole ring (histidine, antifungals)",
    },
    {
        "name": "benzimidazole",
        "smarts": "c1ccc2[nH]cnc2c1",
        "category": "Benzimidazoles",
        "description": "Benzimidazole fused ring (omeprazole, albendazole)",
    },
    {
        "name": "thiazole",
        "smarts": "c1cscn1",
        "category": "Thiazoles",
        "description": "Thiazole ring (penicillin side chain, ritonavir)",
    },
    {
        "name": "oxazole",
        "smarts": "c1cocn1",
        "category": "Oxazoles",
        "description": "Oxazole ring (oxacillin, sulfamethoxazole metabolites)",
    },
    {
        "name": "furan",
        "smarts": "c1ccoc1",
        "category": "Furans",
        "description": "Furan ring (nitrofurantoin, furazolidone)",
    },
    {
        "name": "thiophene",
        "smarts": "c1ccsc1",
        "category": "Thiophenes",
        "description": "Thiophene ring (ticlopidine, duloxetine)",
    },
    {
        "name": "pyrrole",
        "smarts": "c1cc[nH]c1",
        "category": "Pyrroles",
        "description": "Pyrrole ring (porphyrins, atorvastatin)",
    },
    {
        "name": "piperidine",
        "smarts": "C1CCNCC1",
        "category": "Piperidines",
        "description": "Piperidine ring (fentanyl, methylphenidate)",
    },
    {
        "name": "piperazine",
        "smarts": "C1CNCCN1",
        "category": "Piperazines",
        "description": "Piperazine ring (ciprofloxacin, aripiprazole)",
    },
    {
        "name": "morpholine",
        "smarts": "C1COCCN1",
        "category": "Morpholines",
        "description": "Morpholine ring (linezolid, gefitinib)",
    },
    {
        "name": "triazole_123",
        "smarts": "c1nn[nH]c1",
        "category": "Triazoles",
        "description": "1,2,3-Triazole ring (click chemistry products)",
    },
    {
        "name": "triazole_124",
        "smarts": "c1ncn[nH]1",
        "category": "Triazoles",
        "description": "1,2,4-Triazole ring (fluconazole, letrozole)",
    },
    {
        "name": "tetrazole",
        "smarts": "c1nnn[nH]1",
        "category": "Tetrazoles",
        "description": "Tetrazole ring (losartan, valsartan bioisostere)",
    },
    # ===== Functional Groups (~15) =====
    {
        "name": "phenol",
        "smarts": "[OX2H]c1ccccc1",
        "category": "Phenols",
        "description": "Phenol group (tyrosine, estradiol)",
    },
    {
        "name": "aniline",
        "smarts": "[NX3H2]c1ccccc1",
        "category": "Anilines",
        "description": "Aniline group (sulfonamide precursors)",
    },
    {
        "name": "sulfonamide",
        "smarts": "[NX3][SX4](=[OX1])(=[OX1])",
        "category": "Sulfonamides",
        "description": "Sulfonamide group (antibiotics, diuretics)",
    },
    {
        "name": "urea",
        "smarts": "[NX3][CX3](=[OX1])[NX3]",
        "category": "Ureas",
        "description": "Urea group (sorafenib, hydroxyurea)",
    },
    {
        "name": "carbamate",
        "smarts": "[NX3][CX3](=[OX1])[OX2]",
        "category": "Carbamates",
        "description": "Carbamate group (rivastigmine, carbaryl)",
    },
    {
        "name": "amide",
        "smarts": "[NX3][CX3](=[OX1])[#6]",
        "category": "Amides",
        "description": "Amide bond (peptides, lidocaine)",
    },
    {
        "name": "ester",
        "smarts": "[#6][CX3](=[OX1])[OX2][#6]",
        "category": "Esters",
        "description": "Ester group (aspirin, prodrugs)",
    },
    {
        "name": "hydroxamic_acid",
        "smarts": "[OX2H][NX3][CX3](=[OX1])",
        "category": "Hydroxamic Acids",
        "description": "Hydroxamic acid group (HDAC inhibitors, vorinostat)",
    },
    {
        "name": "nitro",
        "smarts": "[NX3+](=O)[O-]",
        "category": "Nitro Compounds",
        "description": "Nitro group (nitrofurantoin, nitroglycerin)",
    },
    {
        "name": "nitrile",
        "smarts": "[CX2]#[NX1]",
        "category": "Nitriles",
        "description": "Nitrile group (cyanide, letrozole)",
    },
    {
        "name": "thiol",
        "smarts": "[SX2H]",
        "category": "Thiols",
        "description": "Thiol group (cysteine, captopril)",
    },
    {
        "name": "disulfide",
        "smarts": "[SX2][SX2]",
        "category": "Disulfides",
        "description": "Disulfide bridge (cystine, glutathione oxidized)",
    },
    {
        "name": "aldehyde",
        "smarts": "[CX3H1](=O)",
        "category": "Aldehydes",
        "description": "Aldehyde group (formaldehyde, retinal)",
    },
    {
        "name": "guanidine",
        "smarts": "[NX3][CX3](=[NX2])[NX3]",
        "category": "Guanidines",
        "description": "Guanidine group (arginine, metformin)",
    },
    {
        "name": "lactone",
        "smarts": "[#6]1~[#6]~[#6]~[#8]~[#6](~[#8])1",
        "category": "Lactones",
        "description": "Lactone ring (statins, macrolides)",
    },
    # ===== Pharmacophoric / Drug-class (~15) =====
    {
        "name": "barbiturate",
        "smarts": "O=C1CC(=O)NC(=O)N1",
        "category": "Barbiturates",
        "description": "Barbiturate ring (phenobarbital, thiopental)",
    },
    {
        "name": "hydantoin",
        "smarts": "O=C1NC(=O)NC1",
        "category": "Hydantoins",
        "description": "Hydantoin ring (phenytoin, fosphenytoin)",
    },
    {
        "name": "steroid",
        "smarts": "[#6]~1~[#6]~[#6]~[#6]~2~[#6]~[#6]~[#6]~[#6]~3~[#6]~[#6]~[#6]~[#6]~4~[#6]~[#6]~[#6]~[#6]~[#6]~4~[#6]~3~[#6]~2~[#6]~1",
        "category": "Steroids",
        "description": "Steroid four-fused-ring skeleton (cholesterol, dexamethasone)",
    },
    {
        "name": "peptide_bond",
        "smarts": "[NX3][CX3](=[OX1])[CX4][NX3]",
        "category": "Peptide-like",
        "description": "Peptide bond motif (peptides, peptidomimetics)",
    },
    {
        "name": "crown_ether",
        "smarts": "[OX2]1~[#6]~[#6]~[OX2]~[#6]~[#6]~[OX2]~[#6]~[#6]1",
        "category": "Crown Ethers",
        "description": "Crown ether macrocycle (18-crown-6)",
    },
    {
        "name": "azide",
        "smarts": "[NX1]=[NX2]=[NX1]",
        "category": "Azides",
        "description": "Azide group (click chemistry, azidothymidine)",
    },
    {
        "name": "epoxide",
        "smarts": "[OX2r3]1[#6r3][#6r3]1",
        "category": "Epoxides",
        "description": "Epoxide ring (carbamazepine-10,11-epoxide)",
    },
    {
        "name": "aziridine",
        "smarts": "[NX3r3]1[#6r3][#6r3]1",
        "category": "Aziridines",
        "description": "Aziridine ring (mitomycin C, nitrogen mustards)",
    },
    {
        "name": "acetal",
        "smarts": "[OX2][CX4][OX2]",
        "category": "Acetals",
        "description": "Acetal group (carbohydrate anomeric center)",
    },
    {
        "name": "phosphoester",
        "smarts": "[PX4](=[OX1])([OX2])([OX2])[OX2]",
        "category": "Phosphoesters",
        "description": "Phosphoester group (nucleotides, phospholipids)",
    },
    {
        "name": "coumarin",
        "smarts": "[#8]1[#6](=[OX1])[#6]~[#6]c2ccccc12",
        "category": "Coumarins",
        "description": "Coumarin lactone (warfarin, umbelliferone)",
    },
    {
        "name": "biphenyl",
        "smarts": "c1ccc(-c2ccccc2)cc1",
        "category": "Biphenyls",
        "description": "Biphenyl scaffold (valsartan, flurbiprofen)",
    },
    {
        "name": "adamantane",
        "smarts": "C1(CC2CC3CC(C1)CC(C2)C3)",
        "category": "Adamantanes",
        "description": "Adamantane cage (amantadine, memantine)",
    },
    {
        "name": "spiro",
        "smarts": "[R]1~[R]~[R]([R]2~[R]~[R]~[R]~[R]2)~[R]~[R]1",
        "category": "Spiro Compounds",
        "description": "Spiro junction (spironolactone, buspirone)",
    },
]

# ---------------------------------------------------------------------------
# Compiled rules singleton
# ---------------------------------------------------------------------------

_COMPILED: Optional[list[tuple[dict, Chem.Mol]]] = None


def get_compiled_rules() -> list[tuple[dict, Chem.Mol]]:
    """
    Compile all CHEMOTYPE_RULES SMARTS to RDKit mol objects.

    Returns a module-level singleton list of (rule_dict, compiled_mol) tuples.
    Rules whose SMARTS fail to compile are logged as warnings and skipped.

    Returns:
        List of (rule_dict, compiled_mol) tuples.
    """
    global _COMPILED
    if _COMPILED is not None:
        return _COMPILED

    compiled: list[tuple[dict, Chem.Mol]] = []
    for rule in CHEMOTYPE_RULES:
        mol = Chem.MolFromSmarts(rule["smarts"])
        if mol is None:
            logger.warning(
                "Taxonomy rule '%s' SMARTS failed to compile: %s",
                rule["name"],
                rule["smarts"],
            )
            continue
        compiled.append((rule, mol))

    _COMPILED = compiled
    logger.info("Compiled %d/%d taxonomy SMARTS rules", len(compiled), len(CHEMOTYPE_RULES))
    return _COMPILED
