"""AlphaFold structures of an isoform with its vaccine peptides or exons highlighted.

:class:`PeptideStructure` stands in for mkt's ``StructureVisualizer`` (which loads
structures from a ``KinaseInfo``): it reads the isoform's model from the committed
``AlphaFold.tar.gz`` and exposes the same ``pdb_text``/``residues``/
``get_highlight_data()``/``config.str_attr`` interface that
:class:`~dsrct_vax.app.visualizers.StructureVisualizerGenerator` consumes.
:func:`build_highlight_config` colors the peptides, :func:`build_exon_config` the
exons and :func:`build_domain_config` the InterPro domains.
"""

import io
from dataclasses import dataclass, field

from Bio.PDB import MMCIFParser
from mkt.databases.utils import convert_structure2string

from dsrct_vax.io_utils import load_alphafold_cif
from dsrct_vax.schema import GeneInfo, IsoformInfo, PeptideInfo


@dataclass
class HighlightConfig:
    """Residues to highlight: 1-indexed isoform positions with color/style/label."""

    str_attr: str = "Peptides"
    """Highlight category (``"Peptides"``, ``"Exons"`` or ``"Domains"``), keys the
    viewer defaults."""
    list_idx: list[int] = field(default_factory=list)
    list_color: list[str] = field(default_factory=list)
    list_style: list[str] = field(default_factory=list)
    list_label: list[str | None] = field(default_factory=list)


def build_highlight_config(
    obj_iso: IsoformInfo,
    list_peptide: list[PeptideInfo],
    dict_color: dict[str, str],
    str_mod_color: str,
) -> HighlightConfig:
    """Highlight every selected peptide residue present in an isoform.

    Peptides are painted longest first, so a nested peptide (e.g. WT1-A1 inside
    WT1-122A1-long) shows over its host; heteroclitic positions become labelled
    sticks (e.g. ``"R126Y"``, numbered on this isoform).
    """
    dict_hl: dict[int, tuple[str, str, str | None]] = {}
    for obj_pep in sorted(list_peptide, key=lambda p: -len(p.seq_wt)):
        match = obj_pep.dict_isoforms[obj_iso.isoform_id]
        for idx_pep, idx_iso in zip(match.list_idx_present, match.list_idx_isoform):
            aa_mod = obj_pep.dict_mods.get(idx_pep - 1)
            if aa_mod is not None:
                str_label = f"{obj_iso.seq[idx_iso - 1]}{idx_iso}{aa_mod}"
                dict_hl[idx_iso] = (str_mod_color, "stick", str_label)
            elif idx_iso not in dict_hl or dict_hl[idx_iso][1] != "stick":
                dict_hl[idx_iso] = (dict_color[obj_pep.name], "cartoon", None)
    list_idx = sorted(dict_hl)
    return HighlightConfig(
        list_idx=list_idx,
        list_color=[dict_hl[i][0] for i in list_idx],
        list_style=[dict_hl[i][1] for i in list_idx],
        list_label=[dict_hl[i][2] for i in list_idx],
    )


def build_exon_config(
    obj_iso: IsoformInfo, build: str, list_colors: list[str]
) -> HighlightConfig:
    """Color an isoform's residues by exon, labelling each exon at its middle residue.

    Parameters:
    -----------
    obj_iso : IsoformInfo
        Isoform whose residues are colored.
    build : str
        Build whose exon labels are used.
    list_colors : list[str]
        Colors cycled over consecutive exons, so neighbours always differ.

    Returns:
    --------
    HighlightConfig
        Every residue with an exon, as cartoon; ``"Exon <label>"`` on each exon's
        middle residue.
    """
    dict_exon = obj_iso.dict_build[build].dict_idx2exon
    # an exon is a contiguous run of residues; split runs where the label changes
    list_run: list[tuple[str, list[int]]] = []
    for idx in sorted(dict_exon):
        if list_run and list_run[-1][0] == dict_exon[idx] and list_run[-1][1][-1] == idx - 1:
            list_run[-1][1].append(idx)
        else:
            list_run.append((dict_exon[idx], [idx]))

    config = HighlightConfig(str_attr="Exons")
    for i, (str_exon, list_idx) in enumerate(list_run):
        idx_mid = list_idx[len(list_idx) // 2]
        for idx in list_idx:
            config.list_idx.append(idx)
            config.list_color.append(list_colors[i % len(list_colors)])
            config.list_style.append("cartoon")
            config.list_label.append(f"Exon {str_exon}" if idx == idx_mid else None)
    return config


def build_domain_config(
    obj_gene: GeneInfo, str_isoform: str, dict_color: dict[str, str]
) -> HighlightConfig:
    """Color an isoform's residues by the canonical InterPro domain they fall in.

    Domains are in canonical positions; each residue is carried to the isoform by its
    genomic codon, so a domain residue the isoform lacks is skipped. Each domain is
    labelled at its middle residue present in the isoform.

    Parameters:
    -----------
    obj_gene : GeneInfo
        Gene snapshot with ``list_domains``.
    str_isoform : str
        Isoform whose residues are colored.
    dict_color : dict[str, str]
        Domain accession -> color.

    Returns:
    --------
    HighlightConfig
        Domain residues as cartoon, labelled ``"<accession> <start>-<end>"``.
    """
    build = obj_gene.list_builds[0]
    dict_can = obj_gene.dict_isoforms[obj_gene.str_canonical].dict_build[build].dict_idx2key
    dict_key2idx = {
        k: i for i, k in obj_gene.dict_isoforms[str_isoform].dict_build[build].dict_idx2key.items()
    }
    config = HighlightConfig(str_attr="Domains")
    for obj_dom in obj_gene.list_domains:
        list_idx = [
            dict_key2idx[dict_can[i]]
            for i in range(obj_dom.start, obj_dom.end + 1)
            if i in dict_can and dict_can[i] in dict_key2idx
        ]
        if not list_idx:
            continue
        idx_mid = list_idx[len(list_idx) // 2]
        # short labels: neighbouring domains (e.g. zinc fingers) sit close in 3D
        str_label = f"{obj_dom.accession} {obj_dom.start}–{obj_dom.end}"
        for idx in list_idx:
            config.list_idx.append(idx)
            config.list_color.append(dict_color[obj_dom.accession])
            config.list_style.append("cartoon")
            config.list_label.append(str_label if idx == idx_mid else None)
    return config


class PeptideStructure:
    """An isoform's AlphaFold model plus the peptide residues to highlight.

    Parameters
    ----------
    obj_iso : IsoformInfo
        Isoform whose AlphaFold model is loaded.
    config : HighlightConfig
        Residues to highlight (see :func:`build_highlight_config`,
        :func:`build_exon_config` and :func:`build_domain_config`).
    """

    def __init__(self, obj_iso: IsoformInfo, config: HighlightConfig):
        if obj_iso.alphafold is None:
            raise ValueError(f"No AlphaFold model for {obj_iso.isoform_id}")
        str_cif = load_alphafold_cif(obj_iso.alphafold.entry_id)
        if str_cif is None:
            raise ValueError(f"{obj_iso.alphafold.entry_id} not in AlphaFold.tar.gz")
        self.config = config
        self.structure = MMCIFParser(QUIET=True).get_structure(
            obj_iso.isoform_id, io.StringIO(str_cif)
        )
        self.pdb_text = convert_structure2string(self.structure)
        self.residues = list(self.structure.get_residues())

    def get_highlight_data(
        self,
    ) -> tuple[list[int], dict[int, str], dict[int, str], dict[int, str | None]]:
        """Return highlight positions and color/style/label maps (mkt interface)."""
        cfg = self.config
        return (
            cfg.list_idx,
            dict(zip(cfg.list_idx, cfg.list_color)),
            dict(zip(cfg.list_idx, cfg.list_style)),
            dict(zip(cfg.list_idx, cfg.list_label)),
        )
