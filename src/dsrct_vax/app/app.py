"""Streamlit dashboard of the GPS vaccine peptides across WT1 isoforms.

Reads only committed package data (``PeptideInfo.tar.gz``, ``AlphaFold.tar.gz`` and
the gene snapshot), so it runs hosted without network access to the sources.

Run locally with ``streamlit run src/dsrct_vax/app/app.py``.
"""

import logging
from collections import defaultdict

import pandas as pd
import streamlit as st
from bokeh.layouts import column
from streamlit_bokeh import streamlit_bokeh

from dsrct_vax.app.constants import (
    DICT_BUILD_LABEL,
    DICT_PEPTIDE_COLOR,
    DICT_RESOURCE_URLS,
    INT_PANEL_ISOFORM_EXTRA,
    INT_PANEL_PEPTIDE_EXTRA,
    INT_PLOT_BORDER_LEFT,
    INT_PLOT_PAD,
    INT_PLOT_WIDTH,
    STR_MOD_COLOR,
)
from dsrct_vax.app.sequences import build_tracks
from dsrct_vax.app.structures import PeptideStructure, build_highlight_config
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
    """One row per peptide: sequences, MHC class, canonical position and genomic span."""
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
                "Modifications": str_mods or "none",
                "Canonical": f"{obj_pep.int_start}–{obj_pep.int_end}",
                f"Genomic ({build})": span,
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


def setup_sidebar(obj_gene: GeneInfo) -> tuple[list[str], list[str], str]:
    """Render sidebar controls; return selected peptides, isoforms and build."""
    st.sidebar.title("Peptide options")
    list_names = list_peptide_names()
    list_pep = st.sidebar.multiselect("Peptides", list_names, default=list_names)
    list_iso_all = list(obj_gene.dict_isoforms)
    list_iso = st.sidebar.multiselect("Isoforms", list_iso_all, default=list_iso_all)
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

    with st.expander("Sequences", expanded=True):
        st.markdown(
            "The vaccine peptides (top) sit over the wild-type isoforms (bottom); both "
            "panels pan and zoom together. Isoforms are lined up by genomic codon, so a "
            "skipped exon shows as a gap, and their residues are shaded by the peptide "
            "covering them (a nested peptide shows over its host). Heteroclitic residues "
            "in the vaccine sequences are red. The x-axis gives canonical "
            f"({obj_gene.str_canonical}) positions; exon labels are from {DICT_BUILD_LABEL[build]}."
        )
        tracks = build_tracks(obj_gene, list_peptide, list_iso, build, DICT_PEPTIDE_COLOR)
        dict_shared = {
            "formatter_code": tracks.formatter_code,
            "exon_map": tracks.exon_map,
            "plot_width": INT_PLOT_WIDTH,
            "int_min_border_left": INT_PLOT_BORDER_LEFT,
        }
        # two stacked panels sharing one x-range so they pan/zoom together
        plot_pep = render_alignment_plot(
            **tracks.return_section(bool_peptide=True),
            **dict_shared,
            int_view_start=tracks.int_first_peptide,
            str_title="Vaccine peptides",
            bool_xaxis=False,
            int_extra_height=INT_PANEL_PEPTIDE_EXTRA,
        )
        plot_iso = render_alignment_plot(
            **tracks.return_section(bool_peptide=False),
            **dict_shared,
            str_title=f"{STR_GENE} isoforms (wild type)",
            bool_exon_labels=False,
            x_range=plot_pep.x_range,
            int_extra_height=INT_PANEL_ISOFORM_EXTRA,
        )
        # streamlit_bokeh sizes the root from its own width/height (a small default
        # when unset, scaled with the container width), and bokeh hands any mismatch
        # to the first panel; each panel has a fixed total height, so the stack is
        # exactly their sum
        layout = column(
            plot_pep,
            plot_iso,
            width=INT_PLOT_BORDER_LEFT + INT_PLOT_WIDTH + INT_PLOT_PAD,
            height=plot_pep.height + plot_iso.height,
        )
        streamlit_bokeh(layout, use_container_width=False, key=f"plot_alignment_{build}")

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
            obj_iso = obj_gene.dict_isoforms[str_iso]
            if obj_iso.alphafold is None:
                st.error(f"No AlphaFold DB model for {str_iso}.", icon="⚠️")
            else:
                try:
                    config = build_highlight_config(
                        obj_iso, list_peptide, DICT_PEPTIDE_COLOR, STR_MOD_COLOR
                    )
                    viz = StructureVisualizerGenerator(PeptideStructure(obj_iso, config))
                    st.markdown(f"#### {obj_iso.alphafold.entry_id}")
                    st.iframe(viz.html, height=600)
                    st.caption(
                        "Peptide residues present in this isoform are colored; "
                        "heteroclitic residues are labelled sticks."
                    )
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

    for obj_pep in list_peptide:
        st.download_button(
            label=f"Download {obj_pep.name} JSON",
            data=obj_pep.model_dump_json(indent=2),
            file_name=f"{obj_pep.name}.json",
            mime="application/json",
            icon=":material/download:",
        )


main()
