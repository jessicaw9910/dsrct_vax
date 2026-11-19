"""InterPro API clients that keep every hit: Pfam signatures and InterPro entries.

:class:`Pfam` is our version of :class:`mkt.databases.pfam.Pfam` (same endpoint,
cached session and provenance stamp). The mkt client keeps only the first location
of each Pfam entry, which suits single-domain kinases but drops repeated domains
such as WT1's zinc fingers; this one returns one row per hit. :class:`InterPro`
queries the integrated InterPro entries the website's domain view shows, through
the same parser. :func:`return_domains` turns either into
:class:`~dsrct_vax.schema.ProteinDomain` records, optionally of one entry type.
"""

import json
import logging
from dataclasses import dataclass

import pandas as pd
from mkt.databases import requests_wrapper, utils_requests
from mkt.databases.pfam import Pfam as PfamMKT

from dsrct_vax.schema import ProteinDomain

logger = logging.getLogger(__name__)

STR_INTERPRO_API = "https://www.ebi.ac.uk/interpro/api/entry"
"""str: InterPro API entry endpoint; the member database and protein follow."""


@dataclass
class Pfam(PfamMKT):
    """Pfam signatures of a UniProt accession, one row per hit (every location).

    ``_pfam`` (named as in mkt) holds a DataFrame with ``accession``, ``name``,
    ``type``, ``start`` and ``end`` (1-based inclusive), plus ``member_databases``,
    ``integrated``, ``score``, ``uniprot`` and ``protein_length``; None if the
    request fails or finds nothing.
    """

    def query_api(self) -> pd.DataFrame | None:
        """Query the InterPro API and return every hit as one row.

        Returns:
        --------
        pd.DataFrame | None
            One row per hit (fragment), sorted by start; None on failure or no hits.
        """
        res = requests_wrapper.get_cached_session().get(
            f"{self.url}{self.uniprot_id}", headers={"Accept": "application/json"}
        )
        self._stamp_from_response(res)
        if not res.ok:
            utils_requests.print_status_code_if_res_not_ok(res)
            return None
        if len(res.text) == 0:
            logger.warning(f"No {type(self).__name__} hits found: {self.uniprot_id}")
            return None

        list_rows = []
        for entry in json.loads(res.text)["results"]:
            dict_meta = entry["metadata"]
            dict_prot = entry["proteins"][0]
            # mkt reads only entry_protein_locations[0]; keep every location
            for loc in dict_prot["entry_protein_locations"]:
                for frag in loc["fragments"]:
                    list_rows.append(
                        {
                            "uniprot": dict_prot["accession"],
                            "protein_length": dict_prot["protein_length"],
                            "accession": dict_meta["accession"],
                            "name": dict_meta["name"],
                            "type": dict_meta["type"],
                            "integrated": dict_meta.get("integrated"),
                            "member_databases": dict_meta.get("member_databases"),
                            "start": frag["start"],
                            "end": frag["end"],
                            "score": loc.get("score"),
                        }
                    )
        return pd.DataFrame(list_rows).sort_values("start", ignore_index=True)


@dataclass
class InterPro(Pfam):
    """Integrated InterPro entries of a UniProt accession, one row per hit."""

    url: str = f"{STR_INTERPRO_API}/interpro/protein/UniProt/"
    """InterPro entries endpoint."""


def return_domains(
    df_hits: pd.DataFrame | None, str_type: str | None = None
) -> list[ProteinDomain]:
    """Return hits from :class:`Pfam` or :class:`InterPro` as domain records.

    Parameters:
    -----------
    df_hits : pd.DataFrame | None
        ``_pfam`` of a :class:`Pfam` or :class:`InterPro` client.
    str_type : str | None
        Keep only this entry type (e.g. ``"domain"``), by default None (all).

    Returns:
    --------
    list[ProteinDomain]
        One record per hit, sorted by start.
    """
    if df_hits is None:
        return []
    if str_type is not None:
        df_hits = df_hits.loc[df_hits["type"] == str_type]
    return [
        ProteinDomain(
            accession=row.accession,
            name=row.name,
            type=row.type,
            start=row.start,
            end=row.end,
            dict_members={
                db: sorted(acc) for db, acc in (row.member_databases or {}).items()
            },
        )
        for row in df_hits.itertuples(index=False)
    ]
