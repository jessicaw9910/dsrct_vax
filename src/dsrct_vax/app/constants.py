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


def return_text_color(str_hex: str) -> str:
    """Return white or near-black, whichever contrasts more with a background color.

    Parameters:
    -----------
    str_hex : str
        Background color as ``"#rrggbb"``.

    Returns:
    --------
    str
        ``"#ffffff"`` or ``"#0b0b0b"`` (WCAG relative-luminance contrast).
    """

    def _channel(int_c: int) -> float:
        float_c = int_c / 255
        return float_c / 12.92 if float_c <= 0.04045 else ((float_c + 0.055) / 1.055) ** 2.4

    r, g, b = (int(str_hex[i : i + 2], 16) for i in (1, 3, 5))
    float_lum = 0.2126 * _channel(r) + 0.7152 * _channel(g) + 0.0722 * _channel(b)
    float_lum_ink = 0.0033  # relative luminance of #0b0b0b
    float_white = 1.05 / (float_lum + 0.05)
    float_ink = (float_lum + 0.05) / (float_lum_ink + 0.05)
    return "#ffffff" if float_white >= float_ink else "#0b0b0b"


STR_MOD_COLOR = "#e34948"
"""str: Heteroclitic (modified) residues."""

STR_RESIDUE_COLOR = "#c3c2b7"
"""str: Isoform residues outside every selected peptide."""

STR_TITLE_INK = "#0b0b0b"
"""str: Primary text ink on light fills (e.g. the grey isoform tags)."""

STR_GAP_COLOR = "#ffffff"
"""str: Alignment gaps."""

LIST_EXON_COLORS = [LIST_SERIES_COLORS[4], "#8f8d84"]
"""list[str]: Alternating structure colors for consecutive exons (magenta and a
neutral, kept apart from the peptide and domain hues); each exon is also labelled, so
color only marks the boundaries."""

LIST_DOMAIN_COLORS = [LIST_SERIES_COLORS[5], LIST_SERIES_COLORS[6]]
"""list[str]: Domain entry colors in order of first occurrence (green, violet);
validated with the four peptide colors as one six-slot palette."""


def return_domain_colors(list_accession: list[str]) -> dict[str, str]:
    """Return domain accession -> color, fixed by order of first occurrence.

    Parameters:
    -----------
    list_accession : list[str]
        Accessions in residue order (repeats allowed).

    Returns:
    --------
    dict[str, str]
        One color per family, cycling :data:`LIST_DOMAIN_COLORS`.
    """
    list_unique = list(dict.fromkeys(list_accession))
    return {
        acc: LIST_DOMAIN_COLORS[i % len(LIST_DOMAIN_COLORS)]
        for i, acc in enumerate(list_unique)
    }


DICT_RESOURCE_URLS = {
    "Peptide source (Koehne et al., JCO 2017)": STR_GPS_PEPTIDES_SOURCE,
    "UniProt": "https://www.uniprot.org/uniprotkb/P19544/entry",
    "GENCODE (UCSC)": "https://genome.ucsc.edu/",
    "AlphaFold DB": "https://alphafold.ebi.ac.uk/",
    "InterPro": "https://www.ebi.ac.uk/interpro/protein/UniProt/P19544/",
    "missense-kinase-toolkit": "https://github.com/choderalab/missense-kinase-toolkit",
}
"""dict[str, str]: Resources linked from the sidebar."""

DICT_BUILD_LABEL = {"GRCh38": "GRCh38 (hg38)", "GRCh37": "GRCh37 (hg19)"}
"""dict[str, str]: Display name of each genome build."""

FLOAT_PX_PER_RESIDUE = 11.5
"""float: Width of one sequence column in pixels; the panels fill the page width and
show as many residues as fit at this size."""

INT_PLOT_BORDER_LEFT = 140
"""int: Left border of both sequence panels, wide enough for the longest track name,
so the stacked frames line up column for column."""

INT_PANEL_PEPTIDE_EXTRA = 70
"""int: Height in pixels of the peptide panel beyond its tracks (title, exon labels)."""

INT_PANEL_ISOFORM_EXTRA = 120
"""int: Height in pixels of the isoform panel beyond its tracks (title, rotated
position labels, axis label, toolbar)."""

INT_IFRAME_PAD = 16
"""int: Extra iframe height in pixels around the embedded sequence panels."""

INT_CDN_CHECK_TTL = 300
"""int: Seconds to cache the Bokeh CDN reachability check before probing again."""
