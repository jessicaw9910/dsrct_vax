"""Streamlit dashboard of the GPS vaccine peptides across WT1 isoforms.

Reads only committed package data (``PeptideInfo.tar.gz``, ``AlphaFold.tar.gz`` and
the gene snapshot), so it runs hosted without network access to the sources.

Run locally with ``streamlit run src/dsrct_vax/app/app.py``.
"""

import logging
from collections import defaultdict

import pandas as pd
import streamlit as st
from bokeh.embed import file_html
from bokeh.layouts import column
from bokeh.resources import CDN, INLINE, Resources

from dsrct_vax.app.constants import (
    DICT_BUILD_LABEL,
    DICT_PEPTIDE_COLOR,
    DICT_RESOURCE_URLS,
    FLOAT_PX_PER_RESIDUE,
    INT_CDN_CHECK_TTL,
    INT_IFRAME_PAD,
    INT_PANEL_ISOFORM_EXTRA,
    INT_PANEL_PEPTIDE_EXTRA,
    INT_PLOT_BORDER_LEFT,
    LIST_EXON_COLORS,
    STR_MOD_COLOR,
    STR_RESIDUE_COLOR,
    STR_TITLE_INK,
    return_domain_colors,
    return_text_color,
)
from dsrct_vax.app.offline import inline_3dmol, return_3dmol_url, return_is_reachable
from dsrct_vax.app.sequences import build_tracks
from dsrct_vax.app.structures import (
    PeptideStructure,
    build_domain_config,
    build_exon_config,
    build_highlight_config,
)
from dsrct_vax.app.visualizers import StructureVisualizerGenerator, render_alignment_plot
from dsrct_vax.constants import STR_GENE
from dsrct_vax.io_utils import list_peptide_names, load_gene_info, load_peptide_info
from dsrct_vax.schema import GeneInfo, PeptideInfo, PeptideStatus

logger = logging.getLogger(__name__)
# module constants live in dsrct_vax.app.constants: Streamlit renders any bare string
# in this script (e.g. a constant's trailing docstring) as markdown


@st.cache_resource
def get_gene() -> GeneInfo:
    """Load the gene snapshot once per server process."""
    return load_gene_info(STR_GENE)


@st.cache_resource
def get_peptide(str_name: str) -> PeptideInfo:
    """Load one peptide from the archive once per server process."""
    return load_peptide_info(str_name)


@st.cache_resource(ttl=INT_CDN_CHECK_TTL)
def get_is_reachable(str_url: str) -> bool:
    """Return whether a CDN URL is reachable, cached for :data:`INT_CDN_CHECK_TTL` s.

    Checked from the server, which is also the browser's machine when the app runs
    locally offline; the cache expiry picks up a reconnect.
    """
    bool_ok = return_is_reachable(str_url)
    if not bool_ok:
        logger.warning(f"{str_url} unreachable; using the offline fallback")
    return bool_ok


def get_bokeh_resources() -> Resources:
    """Return Bokeh's CDN resources if reachable, else inline (offline) resources."""
    return CDN if get_is_reachable(CDN.js_files[0]) else INLINE


def return_structure_html(str_html: str) -> str:
    """Return a py3Dmol viewer's HTML, inlining 3Dmol.js if its CDN is unreachable.

    Parameters:
    -----------
    str_html : str
        Output of ``py3Dmol.view._make_html()``.

    Returns:
    --------
    str
        HTML that loads 3Dmol.js from the CDN, or embeds the committed copy.
    """
    str_url = return_3dmol_url(str_html)
    if str_url is None or get_is_reachable(str_url):
        return str_html
    return inline_3dmol(str_html)


def return_coverage_table(list_peptide: list[PeptideInfo], obj_gene: GeneInfo) -> pd.DataFrame:
    """Peptide x isoform table of where each peptide sits (or how much is missing)."""
    dict_rows = {}
    for obj_pep in list_peptide:
        row = {}
        for iso, m in obj_pep.dict_isoforms.items():
            if m.status == PeptideStatus.present:
                row[iso] = f"{m.int_start}–{m.int_end}"
            elif m.status == PeptideStatus.partial:
                row[iso] = f"partial ({m.float_coverage:.0%})"
            else:
                row[iso] = "absent"
        dict_rows[obj_pep.name] = row
    return pd.DataFrame.from_dict(dict_rows, orient="index")[list(obj_gene.dict_isoforms)]


def return_peptide_table(
    list_peptide: list[PeptideInfo], obj_gene: GeneInfo, build: str
) -> pd.DataFrame:
    """One row per peptide: sequences, length, MHC class, position, genomic span and
    domains."""
    list_rows = []
    for obj_pep in list_peptide:
        str_mods = ", ".join(
            f"{obj_pep.seq_wt[i]}{obj_pep.int_start + i if obj_pep.int_start else i + 1}{aa}"
            for i, aa in obj_pep.dict_mods.items()
        )
        span = None
        if obj_pep.int_start is not None:
            span = obj_gene.return_genomic_span(
                obj_gene.str_canonical, build, obj_pep.int_start, obj_pep.int_end
            )
        list_rows.append(
            {
                "Peptide": obj_pep.name,
                "MHC class": obj_pep.mhc_class,
                "Wild type": obj_pep.seq_wt,
                "Vaccine": obj_pep.seq_mod,
                "Length (aa)": len(obj_pep.seq_wt),
                "Modifications": str_mods or "none",
                "Canonical": f"{obj_pep.int_start}–{obj_pep.int_end}",
                f"Genomic ({build})": span,
                "Domains": "; ".join(
                    f"{d.accession} {d.name} ({d.start}–{d.end})"
                    for d in obj_pep.list_domains
                )
                or "none",
            }
        )
    return pd.DataFrame(list_rows).set_index("Peptide")


def return_isoform_table(obj_gene: GeneInfo, build: str) -> pd.DataFrame:
    """One row per isoform: length, encoding transcript and AlphaFold model."""
    list_rows = []
    for iso, obj_iso in obj_gene.dict_isoforms.items():
        mapping = obj_iso.dict_build[build]
        str_method = mapping.method.value
        if mapping.int_offset:
            str_method += f" (+{mapping.int_offset} aa)"
        af = obj_iso.alphafold
        list_rows.append(
            {
                "Isoform": iso + (" (canonical)" if obj_iso.bool_canonical else ""),
                "Length": len(obj_iso.seq),
                "Transcript": mapping.transcript_id,
                "Mapping": str_method,
                "AlphaFold": f"{af.entry_id} v{af.version}" if af else None,
            }
        )
    return pd.DataFrame(list_rows).set_index("Isoform")


def return_exon_table(obj_gene: GeneInfo, build: str) -> pd.DataFrame:
    """One row per coding exon: genomic coding segments and each isoform's residues."""
    dict_range: dict[str, dict[str, str]] = defaultdict(dict)
    for iso, obj_iso in obj_gene.dict_isoforms.items():
        dict_idx: dict[str, list[int]] = defaultdict(list)
        for idx, label in obj_iso.dict_build[build].dict_idx2exon.items():
            dict_idx[label].append(idx)
        for label, list_idx in dict_idx.items():
            dict_range[label][iso] = f"{min(list_idx)}–{max(list_idx)}"
    chrom = obj_gene.dict_transcripts[build][0].chrom
    dict_span = obj_gene.return_exon_spans(build)
    bool_rev = obj_gene.strand == "-"
    # transcription order
    list_label = sorted(dict_span, key=lambda k: dict_span[k][0][0], reverse=bool_rev)
    list_rows = [
        {
            "Exon": label,
            f"Coding span ({build})": "; ".join(
                f"{chrom}:{lo:,}-{hi:,}" for lo, hi in dict_span[label]
            ),
        }
        | {iso: dict_range[label].get(iso, "—") for iso in obj_gene.dict_isoforms}
        for label in list_label
    ]
    return pd.DataFrame(list_rows).set_index("Exon")


def return_domain_table(obj_gene: GeneInfo, build: str) -> pd.DataFrame:
    """One row per InterPro domain hit on the canonical isoform.

    Parameters:
    -----------
    obj_gene : GeneInfo
        Gene snapshot.
    build : str
        Build of the genomic span column.

    Returns:
    --------
    pd.DataFrame
        Accession, name, canonical residues, genomic span and member signatures,
        in residue order.
    """
    return pd.DataFrame(
        [
            {
                "Accession": dom.accession,
                "Name": dom.name,
                "Residues": f"{dom.start}–{dom.end}",
                f"Genomic ({build})": obj_gene.return_genomic_span(
                    obj_gene.str_canonical, build, dom.start, dom.end
                ),
                "Member signatures": "; ".join(
                    f"{db}: {', '.join(acc)}" for db, acc in dom.dict_members.items()
                ),
            }
            for dom in obj_gene.list_domains
        ]
    )


def return_tag_css() -> str:
    """Return CSS coloring the sidebar's selected-value tags.

    Peptide tags take their legend color (text in whichever of white or near-black
    contrasts more); isoform tags are the neutral residue grey. Tags are matched by
    the multiselect's key class and the tag's ``aria-label``.

    Returns:
    --------
    str
        A ``<style>`` block for :func:`streamlit.html`.
    """
    list_rule = [
        f".st-key-sel_isoforms [data-tag] {{ background-color: {STR_RESIDUE_COLOR} "
        f"!important; color: {STR_TITLE_INK} !important; }}"
    ]
    for str_name, str_color in DICT_PEPTIDE_COLOR.items():
        list_rule.append(
            f'.st-key-sel_peptides [data-tag][aria-label="{str_name}"] '
            f"{{ background-color: {str_color} !important; "
            f"color: {return_text_color(str_color)} !important; }}"
        )
    return "<style>\n" + "\n".join(list_rule) + "\n</style>"


def setup_sidebar(obj_gene: GeneInfo) -> tuple[list[str], list[str], str]:
    """Render sidebar controls; return selected peptides, isoforms and build."""
    st.sidebar.title("Peptide options")
    list_names = list_peptide_names()
    list_pep = st.sidebar.multiselect(
        "Peptides", list_names, default=list_names, key="sel_peptides"
    )
    list_iso_all = list(obj_gene.dict_isoforms)
    list_iso = st.sidebar.multiselect(
        "Isoforms", list_iso_all, default=list_iso_all, key="sel_isoforms"
    )
    st.html(return_tag_css())
    bool_37 = st.sidebar.checkbox(
        "Use GRCh37 (hg19) coordinates",
        value=False,
        help="Genomic spans and exon labels from GRCh37 instead of GRCh38. "
        "Exon labels follow the MANE Select transcript's numbering on that build.",
    )
    build = "GRCh37" if bool_37 else "GRCh38"
    st.sidebar.markdown("## Legend")
    for name in list_names:
        st.sidebar.markdown(
            f"<span style='color:{DICT_PEPTIDE_COLOR[name]}'>■</span> {name}",
            unsafe_allow_html=True,
        )
    st.sidebar.markdown(
        f"<span style='color:{STR_MOD_COLOR}'>■</span> heteroclitic residue",
        unsafe_allow_html=True,
    )
    if obj_gene.list_domains:
        st.sidebar.markdown("**Domains (InterPro)**")
        dict_name = {d.accession: d.name for d in obj_gene.list_domains}
        for acc, color in return_domain_colors(list(dict_name)).items():
            st.sidebar.markdown(
                f"<span style='color:{color}'>■</span> {acc} {dict_name[acc]}",
                unsafe_allow_html=True,
            )
    st.sidebar.markdown("## Database resources")
    for text, url in DICT_RESOURCE_URLS.items():
        st.sidebar.link_button(text, url)
    st.sidebar.caption(f"Snapshot queried {obj_gene.query_date}.")
    return list_pep, list_iso, build


def main() -> None:
    """Render the dashboard."""
    st.set_page_config(page_title=f"{STR_GENE} vaccine peptides", layout="wide")
    obj_gene = get_gene()
    list_pep, list_iso, build = setup_sidebar(obj_gene)
    list_peptide = [get_peptide(name) for name in list_pep]

    st.title(f"{STR_GENE} vaccine peptides across isoforms")
    if not list_peptide or not list_iso:
        st.info("Select at least one peptide and one isoform in the sidebar.")
        return

    # one download per selected peptide, sized to its label and packed left; the row
    # wraps on a narrow screen
    row_download = st.container(
        horizontal=True, horizontal_alignment="left", gap="small", wrap=True
    )
    for obj_pep in list_peptide:
        row_download.download_button(
            label=f"{obj_pep.name} JSON",
            data=obj_pep.model_dump_json(indent=2),
            file_name=f"{obj_pep.name}.json",
            mime="application/json",
            icon=":material/download:",
            width="content",
        )

    with st.expander("Sequences", expanded=True):
        st.markdown(
            "The vaccine peptides (top) sit over the wild-type isoforms (bottom); both "
            "panels pan and zoom together. Isoforms are lined up by genomic codon, so a "
            "skipped exon shows as a gap, and their residues are shaded by the peptide "
            "covering them (a nested peptide shows over its host). Heteroclitic residues "
            "in the vaccine sequences are red. The domain track above the isoforms shows "
            f"the domains of the canonical isoform ({obj_gene.str_canonical}), colored by "
            "InterPro entry as in the legend. The x-axis gives canonical "
            f"({obj_gene.str_canonical}) positions; exon labels are from {DICT_BUILD_LABEL[build]}."
        )
        tracks = build_tracks(obj_gene, list_peptide, list_iso, build, DICT_PEPTIDE_COLOR)
        dict_shared = {
            "formatter_code": tracks.formatter_code,
            "exon_map": tracks.exon_map,
            "int_min_border_left": INT_PLOT_BORDER_LEFT,
            "bool_stretch_width": True,
        }
        # two stacked panels sharing one x-range so they pan/zoom together
        plot_pep = render_alignment_plot(
            **tracks.return_section(bool_peptide=True),
            **dict_shared,
            int_view_start=tracks.int_first_peptide,
            str_title="Vaccine peptides",
            bool_xaxis=False,
            int_extra_height=INT_PANEL_PEPTIDE_EXTRA,
            float_px_per_col=FLOAT_PX_PER_RESIDUE,
        )
        plot_iso = render_alignment_plot(
            **tracks.return_section(bool_peptide=False),
            **dict_shared,
            str_title=f"{STR_GENE} isoforms (wild type)",
            bool_exon_labels=False,
            x_range=plot_pep.x_range,
            int_extra_height=INT_PANEL_ISOFORM_EXTRA,
        )
        # embedded as a standalone page rather than via streamlit_bokeh, which scales
        # the height with the container width and distorts stacked panels; here the
        # panels stretch to the page width at a fixed total height. BokehJS comes
        # from the CDN, or is inlined when the CDN is unreachable (offline)
        layout = column(plot_pep, plot_iso, sizing_mode="stretch_width")
        st.iframe(
            file_html(layout, get_bokeh_resources(), title="Sequences"),
            height=plot_pep.height + plot_iso.height + INT_IFRAME_PAD,
        )

    with st.expander("Peptide coverage by isoform", expanded=True):
        st.markdown(
            "1-based position of each peptide in each UniProt isoform. *Partial* means "
            "some peptide residues sit at the same genomic codons but the rest are "
            "spliced out; *absent* means the isoform has none of them."
        )
        st.dataframe(return_coverage_table(list_peptide, obj_gene)[list_iso])
        st.dataframe(return_peptide_table(list_peptide, obj_gene, build))

    col1, col2 = st.columns(2)
    with col1:
        with st.expander("Structure", expanded=True):
            str_iso = st.selectbox("AlphaFold model of isoform", list_iso, index=0)
            str_view = st.radio(
                "Color structure by", ["Peptides", "Exons", "Domains"], horizontal=True
            )
            obj_iso = obj_gene.dict_isoforms[str_iso]
            dict_domain_color = return_domain_colors(
                [d.accession for d in obj_gene.list_domains]
            )
            if obj_iso.alphafold is None:
                st.error(f"No AlphaFold DB model for {str_iso}.", icon="⚠️")
            else:
                try:
                    if str_view == "Domains":
                        config = build_domain_config(
                            obj_gene, str_iso, dict_domain_color
                        )
                        str_caption = (
                            "Residues in an InterPro domain of the canonical isoform are "
                            "colored by entry (see the legend) and each domain is "
                            "labelled with its accession and canonical residues; residues "
                            "this isoform "
                            "lacks are skipped."
                        )
                    elif str_view == "Exons":
                        config = build_exon_config(obj_iso, build, LIST_EXON_COLORS)
                        str_caption = (
                            "Consecutive exons alternate in color and each is labelled "
                            "at its middle residue; numbering follows the MANE Select "
                            f"transcript on {DICT_BUILD_LABEL[build]}."
                        )
                    else:
                        config = build_highlight_config(
                            obj_iso, list_peptide, DICT_PEPTIDE_COLOR, STR_MOD_COLOR
                        )
                        str_caption = (
                            "Peptide residues present in this isoform are colored; "
                            "heteroclitic residues are labelled sticks."
                        )
                    viz = StructureVisualizerGenerator(PeptideStructure(obj_iso, config))
                    st.markdown(f"#### {obj_iso.alphafold.entry_id}")
                    st.iframe(return_structure_html(viz.html), height=600)
                    st.caption(str_caption)
                except Exception as e:
                    logger.exception(f"Error rendering {str_iso}: {e}")
                    st.error(f"Could not render the structure of {str_iso}.", icon="⚠️")
    with col2:
        with st.expander("Isoforms", expanded=True):
            st.dataframe(return_isoform_table(obj_gene, build).loc[
                [i + (" (canonical)" if obj_gene.dict_isoforms[i].bool_canonical else "")
                 for i in list_iso]
            ])
        with st.expander("Exons", expanded=True):
            st.markdown(
                "Exons are numbered on the MANE Select transcript; an exon outside it is "
                "labelled by the MANE exon upstream plus a letter (e.g. *1b*). Ranges are "
                "1-based isoform residues."
            )
            st.dataframe(return_exon_table(obj_gene, build)[
                [f"Coding span ({build})"] + list_iso
            ])
        if obj_gene.list_domains:
            with st.expander("Domains", expanded=True):
                st.markdown(
                    f"InterPro domains of the canonical isoform ({obj_gene.str_canonical}), "
                    "with the member-database signatures each one integrates (Pfam, PROSITE, "
                    "SMART, …); residues are 1-based canonical positions."
                )
                st.dataframe(return_domain_table(obj_gene, build), hide_index=True)


main()
