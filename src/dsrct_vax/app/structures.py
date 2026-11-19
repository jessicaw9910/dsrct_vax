"""AlphaFold structures of an isoform with its vaccine peptides highlighted.

:class:`PeptideStructure` stands in for mkt's ``StructureVisualizer`` (which loads
structures from a ``KinaseInfo``): it reads the isoform's model from the committed
``AlphaFold.tar.gz`` and exposes the same ``pdb_text``/``residues``/
``get_highlight_data()``/``config.str_attr`` interface that
:class:`~dsrct_vax.app.visualizers.StructureVisualizerGenerator` consumes.
"""

import io
from dataclasses import dataclass, field

from Bio.PDB import MMCIFParser
from mkt.databases.utils import convert_structure2string

from dsrct_vax.io_utils import load_alphafold_cif
from dsrct_vax.schema import IsoformInfo, PeptideInfo


@dataclass
class PeptideHighlightConfig:
    """Residues to highlight: 1-indexed isoform positions with color/style/label."""

    str_attr: str = "Peptides"
    """Highlight category, keys the viewer's opacity/style defaults."""
    list_idx: list[int] = field(default_factory=list)
    list_color: list[str] = field(default_factory=list)
    list_style: list[str] = field(default_factory=list)
    list_label: list[str | None] = field(default_factory=list)


def build_highlight_config(
    obj_iso: IsoformInfo,
    list_peptide: list[PeptideInfo],
    dict_color: dict[str, str],
    str_mod_color: str,
) -> PeptideHighlightConfig:
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
    return PeptideHighlightConfig(
        list_idx=list_idx,
        list_color=[dict_hl[i][0] for i in list_idx],
        list_style=[dict_hl[i][1] for i in list_idx],
        list_label=[dict_hl[i][2] for i in list_idx],
    )


class PeptideStructure:
    """An isoform's AlphaFold model plus the peptide residues to highlight.

    Parameters
    ----------
    obj_iso : IsoformInfo
        Isoform whose AlphaFold model is loaded.
    config : PeptideHighlightConfig
        Residues to highlight (see :func:`build_highlight_config`).
    """

    def __init__(self, obj_iso: IsoformInfo, config: PeptideHighlightConfig):
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
