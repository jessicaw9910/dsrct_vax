"""Regenerate the committed package data for a gene.

Queries UniProt, UCSC (GENCODE) and AlphaFold DB, then writes ``data/<gene>.json``,
adds the isoform models to ``data/AlphaFold.tar.gz`` and serializes the peptides into
``data/PeptideInfo.tar.gz``.

Example
-------
``generate_data --gene WT1``
"""

import logging

import typer

from dsrct_vax.constants import DICT_GPS_PEPTIDES, STR_GENE
from dsrct_vax.io_utils import (
    generate_peptide_dict,
    save_alphafold_tar,
    save_gene_info,
    save_peptide_tar,
)
from dsrct_vax.isoforms import generate_gene_info
from dsrct_vax.schema import PeptideStatus

app = typer.Typer(add_completion=False)
logger = logging.getLogger(__name__)


@app.command()
def main(gene: str = typer.Option(STR_GENE, help="HGNC symbol of the parent gene.")):
    """Query every source for GENE and rewrite the committed package data."""
    logging.basicConfig(level=logging.INFO, format="%(levelname)s %(name)s: %(message)s")
    obj_gene, dict_cif = generate_gene_info(gene)
    typer.echo(f"Wrote {save_gene_info(obj_gene)}")
    typer.echo(f"Wrote {save_alphafold_tar(dict_cif)} ({len(dict_cif)} models)")

    dict_def = {k: v for k, v in DICT_GPS_PEPTIDES.items() if v.get("gene", STR_GENE) == gene}
    dict_peptide = generate_peptide_dict(dict_def)
    typer.echo(f"Wrote {save_peptide_tar(dict_peptide)} ({len(dict_peptide)} peptides)")

    for iso, obj_iso in obj_gene.dict_isoforms.items():
        str_map = ", ".join(
            f"{b}: {m.method} {m.transcript_id}" for b, m in obj_iso.dict_build.items()
        )
        typer.echo(f"  {iso:10s} {len(obj_iso.seq):4d} aa | {str_map}")
    for name, obj in dict_peptide.items():
        list_flag = [
            f"{iso}:{m.status}"
            + (f"({m.float_coverage:.0%})" if m.status == PeptideStatus.partial else "")
            for iso, m in obj.dict_isoforms.items()
            if m.status != PeptideStatus.present
        ]
        typer.echo(f"  {name}: not fully present in {list_flag or 'none'}")


if __name__ == "__main__":
    app()
