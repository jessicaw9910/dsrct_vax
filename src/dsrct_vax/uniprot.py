"""UniProt isoform client.

:class:`UniProtIsoformFASTA` follows ``mkt.databases.uniprot``: it extends mkt's
:class:`~mkt.databases.uniprot.UniProt` base, inherits
:class:`~mkt.databases.api_schema.RESTAPIClient` and queries on instantiation. The
Ensembl cross-references per isoform come from mkt's
:class:`~mkt.databases.uniprot.UniProtJSON` record via :func:`return_isoform_xrefs`.
"""

import io
import logging
from dataclasses import dataclass, field

from Bio import SeqIO
from mkt.databases import requests_wrapper, utils_requests
from mkt.databases.api_schema import RESTAPIClient
from mkt.databases.uniprot import UniProt, UniProtJSON

logger = logging.getLogger(__name__)


@dataclass
class UniProtIsoformFASTA(UniProt, RESTAPIClient):
    """Class to download every isoform sequence of a UniProt entry."""

    url_query: str | None = None
    """Query URL."""
    _dict_seq: dict[str, str] | None = field(default=None, repr=False)
    """Isoform accession -> sequence, canonical first; None if the request failed."""

    def __post_init__(self):
        self.create_query_url()
        self._dict_seq = self.query_api()

    def create_query_url(self) -> None:
        """Build the search URL that returns the entry and all its isoforms."""
        self.url_query = (
            f"{self.url}/search?query=accession:{self.uniprot_id}"
            "&includeIsoform=true&format=fasta"
        )

    def query_api(self) -> dict[str, str] | None:
        """Query UniProt and return isoform sequences; None if the request fails.

        The canonical isoform is keyed by the bare accession (as in the FASTA), the
        others by ``<accession>-<n>``.
        """
        res = requests_wrapper.get_cached_session().get(self.url_query)
        self._stamp_from_response(res)
        if not res.ok:
            utils_requests.print_status_code_if_res_not_ok(res)
            return None
        dict_seq = {
            rec.id.split("|")[1]: str(rec.seq)
            for rec in SeqIO.parse(io.StringIO(res.text), "fasta")
        }
        return dict(
            sorted(
                dict_seq.items(),
                key=lambda kv: (kv[0] != self.uniprot_id, len(kv[0]), kv[0]),
            )
        )


def return_isoform_xrefs(obj_json: UniProtJSON) -> dict[str, list[str]]:
    """Return isoform accession -> unversioned Ensembl transcripts from a UniProt record.

    UniProt names the canonical isoform ``<accession>-1``; it is keyed by the bare
    accession to match :class:`UniProtIsoformFASTA`.
    """
    uniprot_id = obj_json.uniprot_id
    dict_xref: dict[str, list[str]] = {}
    for xref in (obj_json._json or {}).get("uniProtKBCrossReferences", []):
        if xref["database"] != "Ensembl":
            continue
        iso = xref.get("isoformId") or uniprot_id
        iso = uniprot_id if iso == f"{uniprot_id}-1" else iso
        dict_xref.setdefault(iso, []).append(xref["id"].split(".")[0])
    return dict_xref
