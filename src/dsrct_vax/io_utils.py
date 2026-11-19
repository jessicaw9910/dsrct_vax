"""I/O utilities for dsrct_vax.

Committed package data (``src/dsrct_vax/data``), read without network access:

- ``PeptideInfo.tar.gz``: one serialized :class:`~dsrct_vax.schema.PeptideInfo` JSON
  per peptide, read one at a time in memory (as mkt does for ``KinaseInfo.tar.gz``).
- ``AlphaFold.tar.gz``: one AlphaFold DB mmCIF per isoform.
- ``<gene>.json``: the :class:`~dsrct_vax.schema.GeneInfo` snapshot.

:func:`return_peptide_dict` mirrors ``mkt.databases.io_utils.return_kinase_dict``.
"""

import io
import json
import logging
import tarfile
from functools import lru_cache
from pathlib import Path

from mkt.schema.io_utils import untar_files_in_memory

from dsrct_vax.constants import (
    DICT_GPS_PEPTIDES,
    PATH_ALPHAFOLD_TAR,
    PATH_DATA,
    PATH_PEPTIDE_TAR,
    STR_ALPHAFOLD_DIR,
    STR_PEPTIDE_DIR,
)
from dsrct_vax.schema import GeneInfo, PeptideInfo

logger = logging.getLogger(__name__)


def return_gene_info_path(str_gene: str, path_dir: Path = PATH_DATA) -> Path:
    """Return the path of a gene's snapshot JSON."""
    return Path(path_dir) / f"{str_gene}.json"


@lru_cache
def load_gene_info(str_gene: str) -> GeneInfo:
    """Load a gene's committed snapshot (cached per gene).

    Raises
    ------
    FileNotFoundError
        If no snapshot exists; run ``generate_data --gene <gene>`` to create one.
    """
    path = return_gene_info_path(str_gene)
    if not path.exists():
        raise FileNotFoundError(
            f"No snapshot for {str_gene} at {path}; run `generate_data --gene {str_gene}`."
        )
    return GeneInfo.model_validate_json(path.read_text())


def save_gene_info(obj_gene: GeneInfo, path_dir: Path = PATH_DATA) -> Path:
    """Write a gene snapshot to ``<path_dir>/<gene>.json``."""
    path = return_gene_info_path(obj_gene.gene, path_dir)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(obj_gene.model_dump(mode="json"), indent=2) + "\n")
    load_gene_info.cache_clear()
    return path


def write_tar(dict_files: dict[str, str], path_tar: Path, str_dir: str) -> Path:
    """Add text files to a ``.tar.gz`` under ``str_dir/``, keeping other entries.

    Parameters
    ----------
    dict_files : dict[str, str]
        File name -> contents to add or replace.
    path_tar : Path
        Archive path; created if absent.
    str_dir : str
        Directory the entries sit under inside the archive.

    Returns
    -------
    Path
        The archive path.
    """
    dict_all = {}
    if Path(path_tar).exists():
        dict_all = {
            Path(k).name: v for k, v in untar_files_in_memory(str(path_tar))[1].items()
        }
    dict_all.update(dict_files)
    Path(path_tar).parent.mkdir(parents=True, exist_ok=True)
    with tarfile.open(path_tar, "w:gz") as tar:
        for name in sorted(dict_all):
            bytes_file = dict_all[name].encode("utf-8")
            info = tarfile.TarInfo(name=f"{str_dir}/{name}")
            info.size = len(bytes_file)
            tar.addfile(info, io.BytesIO(bytes_file))
    return Path(path_tar)


def _read_one(path_tar: Path, str_id: str) -> str | None:
    """Return one entry (by file stem) from an archive, read in memory."""
    if not Path(path_tar).exists():
        return None
    _, dict_files = untar_files_in_memory(str(path_tar), list_ids=[str_id])
    return next(iter(dict_files.values()), None)


def save_alphafold_tar(dict_cif: dict[str, str]) -> Path:
    """Add mmCIF files (name -> contents) to ``AlphaFold.tar.gz``."""
    return write_tar(dict_cif, PATH_ALPHAFOLD_TAR, STR_ALPHAFOLD_DIR)


@lru_cache
def load_alphafold_cif(entry_id: str) -> str | None:
    """Return an AlphaFold mmCIF (e.g. ``"AF-P19544-2-F1"``) or None if absent."""
    return _read_one(PATH_ALPHAFOLD_TAR, entry_id)


def save_peptide_tar(dict_peptide: dict[str, PeptideInfo]) -> Path:
    """Serialize PeptideInfo objects into ``PeptideInfo.tar.gz`` (one JSON each)."""
    dict_json = {
        f"{name}.json": obj.model_dump_json(indent=2) for name, obj in dict_peptide.items()
    }
    return write_tar(dict_json, PATH_PEPTIDE_TAR, STR_PEPTIDE_DIR)


@lru_cache
def load_peptide_info(str_name: str) -> PeptideInfo:
    """Load one PeptideInfo from ``PeptideInfo.tar.gz`` without unpacking the rest."""
    str_json = _read_one(PATH_PEPTIDE_TAR, str_name)
    if str_json is None:
        raise FileNotFoundError(f"{str_name} not in {PATH_PEPTIDE_TAR}")
    return PeptideInfo.model_validate_json(str_json)


def list_peptide_names() -> list[str]:
    """Return the peptide names in ``PeptideInfo.tar.gz`` without extracting them."""
    list_ids, _ = untar_files_in_memory(str(PATH_PEPTIDE_TAR), bool_extract=False)
    return sorted(list_ids, key=lambda x: list(DICT_GPS_PEPTIDES).index(x)
                  if x in DICT_GPS_PEPTIDES else len(DICT_GPS_PEPTIDES))


def generate_peptide_dict(
    dict_in: dict[str, dict] = DICT_GPS_PEPTIDES,
) -> dict[str, PeptideInfo]:
    """Build PeptideInfo objects from their definitions (needs the gene snapshot)."""
    return {str_name: PeptideInfo(name=str_name) for str_name in dict_in}


def return_peptide_dict(list_names: list[str] | None = None) -> dict[str, PeptideInfo]:
    """Return a dictionary of PeptideInfo objects from ``PeptideInfo.tar.gz``.

    Parameters
    ----------
    list_names : list[str] | None, optional
        Peptides to load, by default None (all in the archive).

    Returns
    -------
    dict[str, PeptideInfo]
        Dictionary of PeptideInfo objects with peptide names as keys

    """
    list_names = list_names or list_peptide_names()
    return {str_name: load_peptide_info(str_name) for str_name in list_names}
