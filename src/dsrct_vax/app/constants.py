"""Dashboard constants: colors, captions and resource links."""

from dsrct_vax.constants import DICT_GPS_PEPTIDES, STR_GPS_PEPTIDES_SOURCE

LIST_SERIES_COLORS = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4", "#008300", "#6250d6"]
"""list[str]: Categorical slots in fixed order (dataviz reference palette, light). The
first four pass the CVD and normal-vision checks as adjacent series; two sit below
3:1 contrast, so every track is direct-labelled and a table view is provided."""

DICT_PEPTIDE_COLOR = {
    name: LIST_SERIES_COLORS[i % len(LIST_SERIES_COLORS)]
    for i, name in enumerate(DICT_GPS_PEPTIDES)
}
"""dict[str, str]: Peptide -> color, fixed by definition order so a filter never
repaints the remaining peptides."""

STR_MOD_COLOR = "#e34948"
"""str: Heteroclitic (modified) residues."""

STR_RESIDUE_COLOR = "#c3c2b7"
"""str: Isoform residues outside every selected peptide."""

STR_GAP_COLOR = "#ffffff"
"""str: Alignment gaps."""

DICT_RESOURCE_URLS = {
    "Peptide source (Koehne et al., JCO 2017)": STR_GPS_PEPTIDES_SOURCE,
    "UniProt": "https://www.uniprot.org/uniprotkb/P19544/entry",
    "GENCODE (UCSC)": "https://genome.ucsc.edu/",
    "AlphaFold DB": "https://alphafold.ebi.ac.uk/",
    "missense-kinase-toolkit": "https://github.com/choderalab/missense-kinase-toolkit",
}
"""dict[str, str]: Resources linked from the sidebar."""

DICT_BUILD_LABEL = {"GRCh38": "GRCh38 (hg38)", "GRCh37": "GRCh37 (hg19)"}
"""dict[str, str]: Display name of each genome build."""

INT_PLOT_WIDTH = 920
"""int: Frame width of the sequence panels in pixels."""

INT_PLOT_BORDER_LEFT = 140
"""int: Left border of both sequence panels, wide enough for the longest track name,
so the stacked frames line up column for column."""

INT_PLOT_PAD = 20
"""int: Extra width in pixels to the right of the sequence frames."""

INT_PANEL_PEPTIDE_EXTRA = 70
"""int: Height in pixels of the peptide panel beyond its tracks (title, exon labels)."""

INT_PANEL_ISOFORM_EXTRA = 120
"""int: Height in pixels of the isoform panel beyond its tracks (title, rotated
position labels, axis label, toolbar)."""
