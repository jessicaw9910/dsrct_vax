"""Constants for dsrct_vax."""

from pathlib import Path

DICT_GPS_PEPTIDES = {
    "WT1-A1": {
        "seq_wt": "RMFPNAPYL",  # WT1 126-134
        "mhc_class": "I",  # CD8
        "dict_mods": {0: "Y"},  # R126Y
    },
    "WT1-122A1-long": {
        "seq_wt": "SGQARMFPNAPYLPSCLES",  # WT1 122-140
        "mhc_class": "I/II",  # CD4 + CD8; embeds WT1-A1
        "dict_mods": {4: "Y"},  # R126Y
    },
    "WT1-331-long": {
        "seq_wt": "PGCNKRYFKLSHLQMHSRKHTG",  # WT1 331-352
        "mhc_class": "II",  # CD4
        "dict_mods": {},  # endogenous
    },
    "WT1-427-long": {
        "seq_wt": "RSDELVRHHNMHQRNMTKL",  # WT1 427-445
        "mhc_class": "II",  # CD4
        "dict_mods": {},  # endogenous
    },
}
"""dict[str, dict]: galinpepimut-S (GPS) peptides keyed by name; ``dict_mods`` maps
0-indexed peptide positions to the substituted residue (heteroclitic changes).
Sequences from Koehne et al., "WT1 heteroclitic epitope immunization following
autologous stem cell transplantation in patients with high-risk multiple myeloma
(MM)", J Clin Oncol 35(15_suppl):8016 (2017); see :data:`STR_GPS_PEPTIDES_SOURCE`.
The residue numbers in the comments are as given there; the canonical UniProt
positions are computed (WT1-331-long maps to P19544 328-349)."""

STR_GPS_PEPTIDES_SOURCE = "https://doi.org/10.1200/JCO.2017.35.15_suppl.8016"
"""str: DOI of the source of the GPS peptide sequences (Koehne et al., JCO 2017)."""

STR_GENE = "WT1"
"""str: HGNC symbol of the antigen's parent gene."""

PATH_DATA = Path(__file__).parent / "data"
"""Path: package data directory (gene snapshots and AlphaFold structures)."""

STR_ALPHAFOLD_DIR = "AlphaFold"
"""str: directory name of the AlphaFold mmCIFs inside the archive."""

PATH_ALPHAFOLD_TAR = PATH_DATA / f"{STR_ALPHAFOLD_DIR}.tar.gz"
"""Path: committed archive of AlphaFold DB mmCIF files, one per isoform."""

STR_PEPTIDE_DIR = "PeptideInfo"
"""str: directory name of the PeptideInfo JSONs inside the archive."""

PATH_PEPTIDE_TAR = PATH_DATA / f"{STR_PEPTIDE_DIR}.tar.gz"
"""Path: committed archive of serialized PeptideInfo objects, one JSON per peptide."""
