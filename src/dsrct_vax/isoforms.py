"""Query a gene's UniProt isoforms and AlphaFold models and assemble its snapshot.

Sources:

- UniProt: every isoform sequence and its Ensembl cross-references, via
  :mod:`dsrct_vax.uniprot`.
- GENCODE transcripts on each build, via :mod:`dsrct_vax.gencode`.
- AlphaFold DB: one model per isoform (``AF-<isoform>-F1``), via
  :class:`mkt.databases.alphafold.AlphaFoldStructure`.
- InterPro: entries of type ``"domain"`` on the canonical isoform (every hit), via
  :class:`dsrct_vax.interpro.InterPro`.

:func:`generate_gene_info` assembles them into a :class:`~dsrct_vax.schema.GeneInfo`.
"""

import datetime
import logging

from mkt.databases.alphafold import AlphaFoldStructure
from mkt.databases.hgnc import HGNC
from mkt.databases.uniprot import UniProtJSON

from dsrct_vax.gencode import DICT_BUILD2UCSC, fetch_gencode_transcripts
from dsrct_vax.interpro import InterPro, return_domains
from dsrct_vax.schema import AlphaFoldModel, GeneInfo
from dsrct_vax.uniprot import UniProtIsoformFASTA, return_isoform_xrefs

logger = logging.getLogger(__name__)


def get_uniprot_id(str_gene: str) -> str:
    """Return the canonical UniProt accession of a gene from HGNC."""
    str_field = "uniprot_ids"
    obj_hgnc = HGNC(input_symbol_or_id=str_gene)
    dict_info = obj_hgnc.maybe_get_info_from_hgnc_fetch(list_to_extract=[str_field])
    return dict_info[str_field][0][0]


def fetch_alphafold_models(
    dict_seq: dict[str, str],
) -> tuple[dict[str, AlphaFoldModel], dict[str, str]]:
    """Fetch the AlphaFold DB model of every isoform.

    Parameters
    ----------
    dict_seq : dict[str, str]
        Isoform accession -> sequence.

    Returns
    -------
    tuple[dict[str, AlphaFoldModel], dict[str, str]]
        Isoform accession -> model metadata, and mmCIF file name -> contents.
    """
    dict_model, dict_cif = {}, {}
    for iso, seq in dict_seq.items():
        obj_af = AlphaFoldStructure(uniprot_id=iso)
        if obj_af._json is None or obj_af._cif is None:
            logger.warning(f"No AlphaFold DB model for {iso}")
            continue
        seq_af = obj_af._json.get("uniprotSequence") or obj_af._json.get("sequence")
        entry_id = obj_af._json["entryId"]
        dict_model[iso] = AlphaFoldModel(
            entry_id=entry_id,
            version=obj_af._json.get("latestVersion"),
            model_created=obj_af._json.get("modelCreatedDate"),
            file_name=f"{entry_id}.cif",
            bool_seq_match=seq_af == seq,
        )
        if seq_af != seq:
            logger.warning(f"{entry_id} sequence differs from UniProt {iso}")
        dict_cif[f"{entry_id}.cif"] = obj_af._cif
    return dict_model, dict_cif


def generate_gene_info(str_gene: str) -> tuple[GeneInfo, dict[str, str]]:
    """Query every source and assemble a gene snapshot.

    Parameters
    ----------
    str_gene : str
        HGNC symbol.

    Returns
    -------
    tuple[GeneInfo, dict[str, str]]
        The snapshot and the AlphaFold mmCIF file name -> contents.
    """
    uniprot_id = get_uniprot_id(str_gene)
    dict_seq = UniProtIsoformFASTA(uniprot_id=uniprot_id)._dict_seq
    if not dict_seq:
        raise ValueError(f"No UniProt isoforms for {uniprot_id}")
    dict_xref = return_isoform_xrefs(UniProtJSON(uniprot_id=uniprot_id))
    dict_tx = {b: fetch_gencode_transcripts(str_gene, b) for b in DICT_BUILD2UCSC}
    dict_af, dict_cif = fetch_alphafold_models(dict_seq)
    list_domains = return_domains(InterPro(uniprot_id=uniprot_id)._pfam, "domain")
    obj_gene = GeneInfo(
        gene=str_gene,
        uniprot_id=uniprot_id,
        query_date=datetime.date.today().isoformat(),
        dict_isoform_seq=dict_seq,
        dict_isoform_xref=dict_xref,
        dict_transcripts=dict_tx,
        dict_alphafold=dict_af,
        list_domains=list_domains,
    )
    return obj_gene, dict_cif
