"""Fallbacks that let the dashboard render without internet access.

The sequence panels (BokehJS) and the structure viewer (3Dmol.js) load their
JavaScript from public CDNs. :func:`return_is_reachable` probes a CDN; when it is
unreachable, the app embeds BokehJS inline and :func:`inline_3dmol` swaps py3Dmol's
CDN loader for the committed copy of 3Dmol.js (:data:`PATH_3DMOL_JS`, BSD-3-Clause,
see ``static/3Dmol-LICENSE``).
"""

import re
import urllib.request
from pathlib import Path

STR_3DMOL_VERSION = "2.5.5"
"""str: Version of the committed 3Dmol.js, matching the py3Dmol release in use."""

PATH_3DMOL_JS = Path(__file__).parent / "static" / f"3Dmol-{STR_3DMOL_VERSION}-min.js"
"""Path: Committed 3Dmol.js build used when its CDN is unreachable."""

RE_3DMOL_LOADER = re.compile(r"loadScriptAsync\('(https?://[^']*3Dmol-min\.js)'\)")
"""re.Pattern: py3Dmol's call that fetches 3Dmol.js from the CDN; group 1 is the URL."""


def return_is_reachable(str_url: str, float_timeout: float = 2.0) -> bool:
    """Return whether a URL answers a HEAD request within the timeout.

    Parameters:
    -----------
    str_url : str
        URL to probe.
    float_timeout : float
        Seconds to wait, by default 2.0.

    Returns:
    --------
    bool
        True on any HTTP response below 400, False on error or timeout.
    """
    try:
        req = urllib.request.Request(str_url, method="HEAD")
        with urllib.request.urlopen(req, timeout=float_timeout) as res:
            return res.status < 400
    except Exception:
        return False


def return_3dmol_url(str_html: str) -> str | None:
    """Return the 3Dmol.js CDN URL a py3Dmol viewer's HTML loads, if any.

    Parameters:
    -----------
    str_html : str
        Output of ``py3Dmol.view._make_html()``.

    Returns:
    --------
    str | None
        The CDN URL, or None if the loader call is not found.
    """
    match = RE_3DMOL_LOADER.search(str_html)
    return match.group(1) if match else None


def inline_3dmol(str_html: str) -> str:
    """Replace py3Dmol's CDN loader with the committed 3Dmol.js, inlined.

    Parameters:
    -----------
    str_html : str
        Output of ``py3Dmol.view._make_html()``.

    Returns:
    --------
    str
        The same viewer with 3Dmol.js defined in a leading ``<script>`` and the
        loader resolved immediately; unchanged if the loader call is not found.
    """
    if return_3dmol_url(str_html) is None:
        return str_html
    str_js = PATH_3DMOL_JS.read_text()
    # a lambda keeps re.sub from parsing backslashes in the replacement
    str_body = RE_3DMOL_LOADER.sub(lambda _: "Promise.resolve()", str_html, count=1)
    return f"<script>{str_js}</script>\n{str_body}"
