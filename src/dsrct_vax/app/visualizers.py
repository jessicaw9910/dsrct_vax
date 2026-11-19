"""Bokeh sequence-track and py3Dmol structure viewers for the dashboard.

Copied from the mkt Streamlit app
(``missense-kinase-toolkit/missense_kinase_toolkit/app/visualizers.py``), which is
not pip-installable. Changes from upstream:

- :func:`render_alignment_plot` takes exon labels as strings (e.g. ``"1b"``), sizes
  its height by track count, can open at a given column, and can be stacked as
  titled panels sharing an x-range (axis/exon labels optional per panel; titles sit
  below the exon labels in a fixed dark ink); the
  x-axis formatter is supplied by the caller (canonical positions, not KLIFS).
- ``SequenceAlignmentGenerator`` (KinaseInfo/KLIFS-specific) is not copied.
- :class:`StructureVisualizerGenerator` is unchanged apart from its docstring; it is
  duck-typed and takes :class:`dsrct_vax.app.structures.PeptideStructure`.
"""

import logging
from itertools import chain
from typing import Any

import numpy as np
import py3Dmol
from bokeh.models import (
    ColumnDataSource,
    CustomJSTickFormatter,
    FixedTicker,
    LinearAxis,
    Range1d,
    Span,
    Title,
)
from bokeh.models.glyphs import Rect, Text
from bokeh.plotting import figure

logger = logging.getLogger(__name__)

STR_TITLE_COLOR = "#0b0b0b"
"""str: Panel-title ink, fixed so it stays legible on the white plot under any theme."""


def render_alignment_plot(
    list_sequences: list[str],
    list_ids: list[str],
    list_colors: list[list[str]],
    font_size: int = 9,
    plot_width: int = 1200,
    formatter_code: str | None = None,
    exon_map: dict[int, str] | None = None,
    int_view_start: int = 1,
    int_track_height: int = 22,
    str_title: str | None = None,
    bool_xaxis: bool = True,
    bool_exon_labels: bool = True,
    int_min_border_left: int = 0,
    x_range: Range1d | None = None,
    int_extra_height: int | None = None,
) -> figure:
    """Render a Bokeh sequence alignment plot from pre-computed alignment data.

    Adapted from https://dmnfarrell.github.io/bioinformatics/bokeh-sequence-aligner.

    Parameters:
    -----------
    list_sequences : list[str]
        List of aligned sequence strings (all same length).
    list_ids : list[str]
        List of sequence labels for the y-axis.
    list_colors : list[list[str]]
        List of per-residue color lists, one per sequence.
    font_size : int
        Font size for residue characters, by default 9.
    plot_width : int
        Width of the plot in pixels, by default 1200.
    formatter_code : str | None
        Optional JavaScript code string for a CustomJSTickFormatter on the
        x-axis. If None, default numeric tick labels are used.
    exon_map : dict[int, str] | None
        Column -> exon-label map; when given, an exon ruler (a divider at
        each exon boundary + one centered exon number above the tracks) is drawn.
        Columns are in transcription order, so exons are contiguous runs. By default None.

    int_view_start : int
        First column in the initial (pannable) view, by default 1.
    int_track_height : int
        Pixel height per track (upstream fixes the frame at 100 px), by default 22.
    str_title : str | None
        Panel title, by default None.
    bool_xaxis : bool
        Draw the residue-position axis and toolbar, by default True.
    bool_exon_labels : bool
        Label exons on a top axis (boundaries are drawn either way), by default True.
    int_min_border_left : int
        Left border in pixels; set equal across stacked panels so frames align.
    x_range : Range1d | None
        Shared x-range to link stacked panels, by default None (a new range).
    int_extra_height : int | None
        If given, fix the figure's total height to the tracks plus this many pixels
        (title, axes, toolbar) instead of fixing the frame height, so a stack of
        panels has an exact total height, by default None.

    Returns:
    --------
    figure
        Bokeh figure with the sequence alignment visualization.
    """
    list_text = [i for s in list_sequences for i in s]

    # a missing track (all "-") gets crimson label and residue text so it reads apart from a
    # misaligned one; its cell background keeps the default gap color
    str_missing_color = "#DC143C"
    empty_sequences = [all(c == "-" for c in seq) for seq in list_sequences]
    text_colors = [
        str_missing_color if empty else "black"
        for seq, empty in zip(list_sequences, empty_sequences)
        for _ in seq
    ]
    colors = list(chain(*list_colors))

    N = len(list_sequences[0])
    S = len(list_sequences)

    x = np.arange(1, N + 1)
    y = np.arange(0, S, 1)
    # creates a 2D grid of coords from the 1D arrays
    xx, yy = np.meshgrid(x, y)
    # flattens the arrays
    gx = xx.ravel()
    gy = yy.flatten()
    # use recty for rect coords with an offset
    recty = gy + 0.5
    # now we can create the ColumnDataSource with all the arrays
    source = ColumnDataSource(
        dict(
            x=gx,
            y=gy,
            recty=recty,
            text=list_text,
            colors=colors,
            text_colors=text_colors,
        )
    )

    viewlen = min(N, 80)

    # sequence text view with ability to scroll along x axis
    # view_range is for the close up view
    if x_range is None:
        x_range = Range1d(int_view_start - 9, int_view_start - 1 + viewlen)
    plot_height = max(int_track_height * 2, int_track_height * S)
    dict_height = (
        {"frame_height": plot_height}
        if int_extra_height is None
        else {"height": plot_height + int_extra_height}
    )
    p1 = figure(
        title=None,
        frame_width=plot_width,
        **dict_height,
        x_range=x_range,
        y_range=list_ids,
        tools="xpan, xwheel_zoom, reset, save",
        min_border=0,
        min_border_left=int_min_border_left,
        toolbar_location="below" if bool_xaxis else None,
        background_fill_color="white",
        border_fill_color="white",
    )
    glyph = Text(
        x="x",
        y="y",
        text="text",
        text_align="center",
        text_baseline="bottom",
        text_color="text_colors",
        text_font_size=f"{str(font_size)}pt",
    )
    rects = Rect(
        x="x",
        y="recty",
        width=1,
        height=1,
        fill_color="colors",
        line_color=None,
        fill_alpha=0.4,
    )
    p1.add_glyph(source, glyph)
    p1.add_glyph(source, rects)
    p1.grid.visible = False

    p1.xaxis.ticker = FixedTicker(ticks=list(range(1, N + 1)))

    if formatter_code is not None:
        p1.xaxis.formatter = CustomJSTickFormatter(code=formatter_code)

    # rotate labels 90 degrees
    p1.xaxis.major_label_orientation = np.pi / 2
    # add some space between axis and labels
    p1.xaxis.major_label_standoff = 2
    p1.xaxis.axis_label = "Residue Position"
    p1.xaxis.visible = bool_xaxis

    # show the categorical y-axis labels; upstream hides them and draws Label
    # annotations at data x=0, which are clipped once the view opens past column 1
    p1.yaxis.major_label_text_font_size = "9pt"
    p1.yaxis.major_label_text_font_style = "bold"

    p1.yaxis.minor_tick_line_width = 0
    p1.yaxis.major_tick_line_width = 0
    p1.yaxis.axis_label_text_color = "black"
    p1.yaxis.major_label_text_color = "black"
    p1.xaxis.axis_label_text_color = "black"
    p1.xaxis.major_label_text_color = "black"

    # panel title in a fixed dark ink so the streamlit theme can't fade it against the
    # white plot background; added before the exon axis because bokeh stacks later
    # "above" layouts further out, which puts the title between exon labels and tracks
    if str_title is not None:
        p1.add_layout(
            Title(
                text=str_title,
                text_color=STR_TITLE_COLOR,
                text_font_size="11pt",
                text_font_style="bold",
            ),
            "above",
        )

    # exon ruler: a dotted divider at each exon boundary + one centered exon number above the
    # tracks (only for entries with an exon map). exons are contiguous UniProt-index runs, so we
    # collapse to per-exon segments -- no per-residue clutter.
    if exon_map:
        segments = []  # [exon_no, start_x, end_x]
        for idx, exon_no in sorted((int(k), str(v)) for k, v in exon_map.items()):
            if not 1 <= idx <= N:
                continue
            if segments and segments[-1][0] == exon_no and idx == segments[-1][2] + 1:
                segments[-1][2] = idx
            else:
                segments.append([exon_no, idx, idx])
        for _, _, seg_end in segments[:-1]:
            p1.add_layout(
                Span(
                    location=seg_end + 0.5,
                    dimension="height",
                    line_color="#444444",
                    line_width=1.5,
                    line_dash="dashed",
                    line_alpha=0.9,
                )
            )
        # exon numbers on a dedicated top axis (ticks at each exon center); an axis renders
        # reliably above the tracks, avoiding the categorical-y clipping that floating labels hit
        overrides = {
            (seg_start + seg_end) / 2: f"Exon {exon_no}"
            for exon_no, seg_start, seg_end in segments
        }
    if exon_map and bool_exon_labels:
        exon_axis = LinearAxis(
            ticker=FixedTicker(ticks=list(overrides)),
            major_label_overrides=overrides,
            major_label_text_font_size="11pt",
            major_label_text_color="#333333",
            axis_line_color=None,
            major_tick_line_color=None,
            minor_tick_line_color=None,
        )
        p1.add_layout(exon_axis, "above")

    return p1


DICT_VIZ_OPACITY = {
    "None": 1.0,
    "Peptides": 1.0,
    "lowlight": 0.5,
}
"""dict[str, float]: Opacity for the py3Dmol viewer."""

DICT_VIZ_STYLE = {
    "None": "cartoon",
    "lowlight": "cartoon",
}
"""dict[str, str]: Default styles for py3Dmol viewer."""

DICT_VIZ_COLOR = {
    "None": "spectrum",
    "lowlight": "gray",
}
"""dict[str, str]: Default colors for py3Dmol viewer."""


class StructureVisualizerGenerator:
    """Generate py3Dmol structure visualizations for kinase structures.

    Uses composition with a structure object to access structure data and
    highlight information from its config.

    Parameters
    ----------
    viz : PeptideStructure
        Object exposing ``pdb_text``, ``residues``, ``get_highlight_data()`` and
        ``config.str_attr`` (see :class:`dsrct_vax.app.structures.PeptideStructure`).
    dict_opacity : dict[str, float], optional
        Opacity settings for different highlight types.
    bool_show : bool, optional
        Whether to show the structure in the viewer or return HTML.
    dict_dims : dict[str, int], optional
        Dimensions for the py3Dmol viewer.

    Attributes
    ----------
    viz : PeptideStructure
        The underlying structure visualizer.
    html : str | None
        HTML representation of the py3Dmol viewer.
    """

    def __init__(
        self,
        viz,
        dict_opacity: dict[str, float] | None = None,
        bool_show: bool = False,
        dict_dims: dict[str, int] | None = None,
    ):
        self.viz = viz
        self.dict_opacity = dict_opacity or DICT_VIZ_OPACITY.copy()
        self.bool_show = bool_show
        self.dict_dims = dict_dims or {"width": 600, "height": 600}

        # Generate visualization on init
        self.html = self.visualize_structure()

    @property
    def str_attr(self) -> str:
        """Get the attribute being highlighted from the config."""
        return self.viz.config.str_attr

    def _return_style_dict(
        self,
        str_key: str,
        str_color: str | None = None,
        str_style: str | None = None,
        float_opacity: float | None = None,
    ) -> dict[str, Any]:
        """Return the style dictionary for the given key.

        Parameters
        ----------
        str_key : str
            Key for the style dictionary.
        str_color : str, optional
            Color for the style dictionary, by default None.
        str_style : str, optional
            Style for the style dictionary, by default None.
        float_opacity : float, optional
            Opacity for the style dictionary, by default None.

        Returns
        -------
        dict[str, Any]
            Style dictionary for the given key.
        """
        if str_color is None:
            str_color = DICT_VIZ_COLOR.get(str_key, "gray")
        if str_style is None:
            str_style = DICT_VIZ_STYLE.get(str_key, "cartoon")
        if float_opacity is None:
            float_opacity = self.dict_opacity.get(str_key, 1.0)

        dict_style = {
            str_style: {
                "color": str_color,
                "opacity": float_opacity,
            }
        }

        return dict_style

    def visualize_structure(self) -> str | None:
        """Visualize the structure using py3Dmol.

        Returns
        -------
        str | None
            HTML representation of the py3Dmol viewer or None if bool_show=True.
        """
        view = py3Dmol.view(**self.dict_dims)

        view.addModel(self.viz.pdb_text, "pdb")

        # Get highlight data from the visualizer (which gets it from config)
        list_highlight, dict_color, dict_style, dict_label = (
            self.viz.get_highlight_data()
        )

        # Get opacity for this attribute type
        attr_opacity = self.dict_opacity.get(self.str_attr, 1.0)

        for residue in self.viz.residues:
            res_no = residue.get_id()[1]
            # set lowlight background
            view.setStyle(
                {"resi": str(res_no)},
                self._return_style_dict("lowlight"),
            )
            # add highlights for the selected attribute
            if res_no in list_highlight:
                view.addStyle(
                    {"resi": str(res_no)},
                    self._return_style_dict(
                        str_key=self.str_attr,
                        str_color=dict_color[res_no],
                        str_style=dict_style[res_no],
                        float_opacity=attr_opacity,
                    ),
                )
                # add label if present
                label_text = dict_label.get(res_no)
                if label_text is not None:
                    # addLabel signature: (text, options, sel)
                    # options = label styling, sel = atom selector for positioning
                    view.addLabel(
                        label_text,
                        {
                            "backgroundColor": "white",
                            "fontColor": "black",
                            "fontSize": 10,
                            "showBackground": True,
                            "backgroundOpacity": 0.7,
                        },
                        {"resi": str(res_no), "atom": "CA"},
                    )

        view.zoomTo()
        if self.bool_show:
            view.show()
            return None
        else:
            return view._make_html()
