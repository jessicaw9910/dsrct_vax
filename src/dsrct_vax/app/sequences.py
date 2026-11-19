"""Build the isoform/peptide sequence tracks for :func:`render_alignment_plot`.

Columns are every isoform residue's genomic codon, in transcription order, so
isoforms line up exon by exon (a skipped exon is a gap) without a sequence
alignment. Peptide tracks sit on the columns of an isoform that contains them, and
the domain bars (:func:`build_domain_spans`) on the canonical isoform's columns.
"""

from dataclasses import dataclass, field

from dsrct_vax.app.constants import (
    STR_GAP_COLOR,
    STR_MOD_COLOR,
    STR_RESIDUE_COLOR,
    return_domain_colors,
)
from dsrct_vax.app.visualizers import TrackSpan
from dsrct_vax.schema import GeneInfo, PeptideInfo, PeptideStatus


@dataclass
class SequenceTracks:
    """Inputs for :func:`dsrct_vax.app.visualizers.render_alignment_plot`."""

    list_sequences: list[str]
    list_ids: list[str]
    list_colors: list[list[str]]
    exon_map: dict[int, str]
    """1-based column -> exon label on the chosen build."""
    dict_tick: dict[int, str]
    """1-based column -> canonical residue number ("" where canonical has none)."""
    int_first_peptide: int = 1
    """1-based column of the first peptide residue (opens the plot there)."""
    int_n_peptide: int = 0
    """Number of peptide tracks (the last entries, i.e. drawn on top)."""
    list_spans: list[TrackSpan] = field(default_factory=list)
    """Labelled bars (InterPro domains) on the isoform-section tracks."""

    def return_section(self, bool_peptide: bool) -> dict[str, list]:
        """Return the peptide or isoform tracks as ``render_alignment_plot`` kwargs.

        Parameters:
        -----------
        bool_peptide : bool
            True for the peptide tracks, False for the isoform tracks.

        Returns:
        --------
        dict[str, list]
            ``list_sequences``, ``list_ids``, ``list_colors`` and ``list_spans`` of
            that section.
        """
        int_split = len(self.list_ids) - self.int_n_peptide
        sl = slice(int_split, None) if bool_peptide else slice(0, int_split)
        return {
            "list_sequences": self.list_sequences[sl],
            "list_ids": self.list_ids[sl],
            "list_colors": self.list_colors[sl],
            "list_spans": [] if bool_peptide else self.list_spans,
        }

    @property
    def formatter_code(self) -> str:
        """JavaScript tick formatter that labels columns with canonical positions."""
        return f"var m = {self.dict_tick}; return m[tick] || '';"


def build_domain_spans(
    obj_gene: GeneInfo, dict_key2col: dict[float, int], str_track: str
) -> list[TrackSpan]:
    """Return one labelled bar per InterPro domain hit on the alignment columns.

    A bar runs from the column of the domain's first canonical residue to that of its
    last, so it stays continuous across columns the canonical isoform lacks (e.g.
    exon 1b); it is colored by family (:func:`return_domain_colors`).

    Parameters:
    -----------
    obj_gene : GeneInfo
        Gene snapshot with ``list_domains`` in canonical positions.
    dict_key2col : dict[float, int]
        Residue key on the reference build -> 0-based column.
    str_track : str
        Track id the bars sit on.

    Returns:
    --------
    list[TrackSpan]
        Bars with ``"<accession> <name>"`` labels.
    """
    build_ref = obj_gene.list_builds[0]
    obj_can = obj_gene.dict_isoforms[obj_gene.str_canonical]
    dict_key = obj_can.dict_build[build_ref].dict_idx2key
    dict_color = return_domain_colors([d.accession for d in obj_gene.list_domains])
    list_spans = []
    for obj_dom in obj_gene.list_domains:
        list_col = [
            dict_key2col[dict_key[idx]] + 1
            for idx in range(obj_dom.start, obj_dom.end + 1)
            if idx in dict_key
        ]
        if list_col:
            list_spans.append(
                TrackSpan(
                    str_track=str_track,
                    int_start=min(list_col),
                    int_end=max(list_col),
                    str_color=dict_color[obj_dom.accession],
                    str_label=f"{obj_dom.accession} {obj_dom.name}",
                )
            )
    return list_spans


def build_tracks(
    obj_gene: GeneInfo,
    list_peptide: list[PeptideInfo],
    list_isoforms: list[str],
    build: str,
    dict_color: dict[str, str],
    bool_domains: bool = True,
) -> SequenceTracks:
    """Return peptide tracks (top), the domain track, then isoform tracks.

    Isoform residues inside a selected peptide take the peptide's color (the shortest
    peptide wins where they nest); peptide tracks show the modified sequence with
    heteroclitic residues in :data:`STR_MOD_COLOR`. The domain track heads the isoform
    section when ``bool_domains`` and the snapshot has domains.

    Parameters:
    -----------
    obj_gene : GeneInfo
        Gene snapshot.
    list_peptide : list[PeptideInfo]
        Peptides to draw.
    list_isoforms : list[str]
        Isoforms to draw, top to bottom.
    build : str
        Build whose exon labels are used.
    dict_color : dict[str, str]
        Peptide name -> color.
    bool_domains : bool
        Add the InterPro domain track, by default True.

    Returns:
    --------
    SequenceTracks
        Tracks in bokeh category order (bottom first).
    """
    build_ref = obj_gene.list_builds[0]
    list_key = obj_gene.return_column_keys()
    dict_key2col = {k: i for i, k in enumerate(list_key)}
    n_col = len(list_key)

    exon_map: dict[int, str] = {}
    dict_tick: dict[int, str] = {}
    for iso, obj_iso in obj_gene.dict_isoforms.items():
        dict_key = obj_iso.dict_build[build_ref].dict_idx2key
        dict_exon = obj_iso.dict_build[build].dict_idx2exon
        for idx, key in dict_key.items():
            col = dict_key2col[key] + 1
            if idx in dict_exon:
                exon_map.setdefault(col, dict_exon[idx])
            if obj_iso.bool_canonical:
                dict_tick[col] = str(idx)

    list_seq, list_id, list_col = [], [], []
    list_pep_col: list[int] = []
    for obj_pep in list_peptide:
        list_present = [
            (iso, m) for iso, m in obj_pep.dict_isoforms.items()
            if m.status == PeptideStatus.present
        ]
        chars, colors = ["-"] * n_col, [STR_GAP_COLOR] * n_col
        if list_present:
            iso, match = list_present[0]
            dict_key = obj_gene.dict_isoforms[iso].dict_build[build_ref].dict_idx2key
            for i, idx_iso in enumerate(match.list_idx_isoform):
                col = dict_key2col[dict_key[idx_iso]]
                list_pep_col.append(col + 1)
                chars[col] = obj_pep.seq_mod[i]
                bool_mod = i in obj_pep.dict_mods
                colors[col] = STR_MOD_COLOR if bool_mod else dict_color[obj_pep.name]
        list_seq.append("".join(chars))
        list_id.append(obj_pep.name)
        list_col.append(colors)

    list_spans: list[TrackSpan] = []
    if bool_domains and obj_gene.list_domains:
        str_track = f"Domains ({obj_gene.str_canonical})"
        list_spans = build_domain_spans(obj_gene, dict_key2col, str_track)
        list_seq.append(" " * n_col)
        list_id.append(str_track)
        list_col.append([STR_GAP_COLOR] * n_col)

    for iso in list_isoforms:
        obj_iso = obj_gene.dict_isoforms[iso]
        dict_key = obj_iso.dict_build[build_ref].dict_idx2key
        chars, colors = ["-"] * n_col, [STR_GAP_COLOR] * n_col
        for idx, key in dict_key.items():
            chars[dict_key2col[key]] = obj_iso.seq[idx - 1]
            colors[dict_key2col[key]] = STR_RESIDUE_COLOR
        for obj_pep in sorted(list_peptide, key=lambda p: -len(p.seq_wt)):
            for idx in obj_pep.dict_isoforms[iso].list_idx_isoform:
                colors[dict_key2col[dict_key[idx]]] = dict_color[obj_pep.name]
        str_label = f"{iso} (canonical)" if obj_iso.bool_canonical else iso
        list_seq.append("".join(chars))
        list_id.append(str_label)
        list_col.append(colors)

    # bokeh draws the first category at the bottom; put peptides on top
    return SequenceTracks(
        list_sequences=list_seq[::-1],
        list_ids=list_id[::-1],
        list_colors=list_col[::-1],
        exon_map=exon_map,
        dict_tick=dict_tick,
        int_first_peptide=min(list_pep_col, default=1),
        int_n_peptide=len(list_peptide),
        list_spans=list_spans,
    )
