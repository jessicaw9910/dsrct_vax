"""GENCODE transcripts of a gene, from the UCSC Genome Browser REST API.

UCSC serves the GENCODE ``knownGene`` track for GRCh38 (hg38) and GRCh37 (hg19, lifted
over), with exon blocks, CDS bounds and tags. The REST clients follow
``mkt.databases`` (each inherits :class:`mkt.databases.api_schema.RESTAPIClient`,
queries on instantiation through the shared cached session and stamps provenance):

- :class:`UCSCGeneRegion`: a gene's HGNC span (``/search``).
- :class:`UCSCTrack`: rows of a track in a region (``/getData/track``).
- :class:`UCSCSequence`: genome sequence of a region (``/getData/sequence``).

:func:`fetch_gencode_transcripts` combines them: it locates the gene, pulls every
coding transcript, assembles each CDS from the genome sequence and translates it, so
isoforms can be tied to the transcript that encodes them by sequence. (Genome Nexus'
exon records are out of date for newer transcripts; this reads the same GENCODE
annotation Ensembl publishes.)

:class:`Transcript` holds one transcript and maps its residues to genomic codons and
exons; :func:`return_mane_transcript` picks the transcript that numbers exons.
"""

import logging
from dataclasses import dataclass, field
from typing import Literal

from Bio.Seq import Seq
from mkt.databases import requests_wrapper, utils_requests
from mkt.databases.api_schema import RESTAPIClient
from pydantic import BaseModel

logger = logging.getLogger(__name__)

DICT_BUILD2UCSC = {"GRCh38": "hg38", "GRCh37": "hg19"}
"""dict[str, str]: Genome build -> UCSC assembly name; the first build is the
reference used to order residues across isoforms."""

STR_UCSC_API = "https://api.genome.ucsc.edu"
"""str: UCSC Genome Browser REST API (serves GENCODE ``knownGene`` on both builds)."""

STR_UCSC_TRACK = "knownGene"
"""str: UCSC track of GENCODE transcripts."""

STR_MANE_TAG = "Ensembl_canonical"
"""str: GENCODE tag of the transcript whose exon numbering labels all isoforms."""

INT_REGION_PAD = 50_000
"""int: bp added either side of the HGNC gene span when querying transcripts."""


class Transcript(BaseModel):
    """A GENCODE protein-coding transcript on one genome build."""

    transcript_id: str
    """Versioned Ensembl transcript id (e.g. ``"ENST00000452863.10"``)."""
    build: str
    """Genome build (e.g. ``"GRCh38"``)."""
    chrom: str
    """Chromosome (UCSC style, e.g. ``"chr11"``)."""
    strand: Literal["+", "-"]
    """Genomic strand."""
    list_exons: list[tuple[int, int]]
    """Exon (start, end), 0-based half-open, in transcription order."""
    cds_start: int
    """Genomic start of the CDS (``thickStart``), 0-based."""
    cds_end: int
    """Genomic end of the CDS incl. stop codon (``thickEnd``), 0-based exclusive."""
    list_tags: list[str] = []
    """GENCODE tags (e.g. ``"CCDS"``, ``"Ensembl_canonical"``, ``"non_ATG_start"``)."""
    seq_protein: str
    """Translation of the CDS, stop codon removed."""

    @property
    def id_unversioned(self) -> str:
        """Ensembl id without its version."""
        return self.transcript_id.split(".")[0]

    def return_cds_positions(self) -> list[tuple[int, int]]:
        """Return (genomic position, exon rank) of each CDS base in coding order."""
        list_out = []
        for rank, (start, end) in enumerate(self.list_exons, start=1):
            lo, hi = max(start, self.cds_start), min(end, self.cds_end)
            if lo >= hi:
                continue
            positions = range(lo, hi) if self.strand == "+" else range(hi - 1, lo - 1, -1)
            list_out.extend((pos, rank) for pos in positions)
        return list_out

    def return_residue_codons(self) -> list[tuple[int, int, tuple[int, int, int]]]:
        """Return (middle-base position, exon rank, codon base positions) per residue.

        A residue is assigned to the exon holding its codon's middle base, so a codon
        split across a junction goes with the exon that holds two of its bases.
        """
        list_cds = self.return_cds_positions()
        list_out = []
        for idx in range(len(self.seq_protein)):
            codon = list_cds[3 * idx : 3 * idx + 3]
            if len(codon) < 3:
                break
            list_out.append((codon[1][0], codon[1][1], tuple(pos for pos, _ in codon)))
        return list_out

    def return_coding_exons(self) -> list[tuple[int, tuple[int, int]]]:
        """Return (exon rank, coding (start, end) 0-based half-open) of coding exons."""
        list_out = []
        for rank, (start, end) in enumerate(self.list_exons, start=1):
            lo, hi = max(start, self.cds_start), min(end, self.cds_end)
            if lo < hi:
                list_out.append((rank, (lo, hi)))
        return list_out


def return_mane_transcript(list_tx: list[Transcript]) -> Transcript:
    """Return the MANE Select (``Ensembl_canonical``) transcript, else the first."""
    return next((t for t in list_tx if STR_MANE_TAG in t.list_tags), list_tx[0])


@dataclass
class UCSC:
    """Class to interact with the UCSC Genome Browser REST API."""

    genome: str
    """UCSC assembly (e.g. ``"hg38"``, ``"hg19"``)."""
    url: str = STR_UCSC_API
    """URL for the UCSC REST API."""


@dataclass
class UCSCGeneRegion(UCSC, RESTAPIClient):
    """Class to query a gene's HGNC span on a UCSC assembly."""

    gene: str | None = None
    """HGNC symbol to search for."""
    url_query: str | None = None
    """Query URL."""
    chrom: str | None = None
    """Chromosome of the gene (UCSC style, e.g. ``"chr11"``)."""
    start: int | None = None
    """Gene start, 0-based."""
    end: int | None = None
    """Gene end, 0-based exclusive."""

    def __post_init__(self):
        if self.gene is None:
            raise ValueError("gene must be provided.")
        self.create_query_url()
        self.query_api()

    def create_query_url(self) -> None:
        """Build the search URL from the gene and assembly."""
        self.url_query = f"{self.url}/search?search={self.gene};genome={self.genome}"

    def query_api(self) -> None:
        """Query UCSC search and store the gene's span from its ``hgnc`` match."""
        res = requests_wrapper.get_cached_session().get(self.url_query)
        self._stamp_from_response(res)
        if not res.ok:
            utils_requests.print_status_code_if_res_not_ok(res)
            return
        for track in res.json().get("positionMatches", []):
            if track.get("trackName") != "hgnc":
                continue
            for match in track["matches"]:
                if match.get("posName") == self.gene:
                    self.chrom, span = match["position"].split(":")
                    start, end = (int(x) for x in span.split("-"))
                    self.start, self.end = start - 1, end
                    return
        logger.error(f"{self.gene} not found on {self.genome}")


@dataclass
class UCSCTrack(UCSC, RESTAPIClient):
    """Class to query the rows of a UCSC track within a region."""

    chrom: str | None = None
    """Chromosome (UCSC style)."""
    start: int | None = None
    """Region start, 0-based."""
    end: int | None = None
    """Region end, 0-based exclusive."""
    track: str = STR_UCSC_TRACK
    """Track name."""
    url_query: str | None = None
    """Query URL."""
    _rows: list[dict] | None = field(default=None, repr=False)
    """Track rows; None if the request failed."""

    def __post_init__(self):
        self.create_query_url()
        self._rows = self.query_api()

    def create_query_url(self) -> None:
        """Build the track URL from the assembly, track and region."""
        self.url_query = (
            f"{self.url}/getData/track?genome={self.genome};track={self.track};"
            f"chrom={self.chrom};start={self.start};end={self.end}"
        )

    def query_api(self) -> list[dict] | None:
        """Query the track and return its rows; None if the request fails."""
        res = requests_wrapper.get_cached_session().get(self.url_query)
        self._stamp_from_response(res)
        if res.ok:
            return res.json().get(self.track, [])
        utils_requests.print_status_code_if_res_not_ok(res)
        return None


@dataclass
class UCSCSequence(UCSC, RESTAPIClient):
    """Class to query the genome sequence of a region (plus strand, upper case)."""

    chrom: str | None = None
    """Chromosome (UCSC style)."""
    start: int | None = None
    """Region start, 0-based."""
    end: int | None = None
    """Region end, 0-based exclusive."""
    url_query: str | None = None
    """Query URL."""
    _sequence: str | None = field(default=None, repr=False)
    """Region sequence; None if the request failed."""

    def __post_init__(self):
        self.create_query_url()
        self._sequence = self.query_api()

    def create_query_url(self) -> None:
        """Build the sequence URL from the assembly and region."""
        self.url_query = (
            f"{self.url}/getData/sequence?genome={self.genome};chrom={self.chrom};"
            f"start={self.start};end={self.end}"
        )

    def query_api(self) -> str | None:
        """Query the sequence and return it upper-cased; None if the request fails."""
        res = requests_wrapper.get_cached_session().get(self.url_query)
        self._stamp_from_response(res)
        if res.ok:
            return res.json()["dna"].upper()
        utils_requests.print_status_code_if_res_not_ok(res)
        return None


def _parse_transcript(row: dict, dna: str, offset: int, build: str) -> Transcript:
    """Build a :class:`Transcript` from a ``knownGene`` row and its region sequence.

    Parameters
    ----------
    row : dict
        ``knownGene`` row (bigGenePred fields).
    dna : str
        Plus-strand genome sequence starting at 0-based ``offset``.
    offset : int
        Genomic position of ``dna[0]``.
    build : str
        Genome build of the row.
    """
    list_start = [int(x) for x in row["chromStarts"].rstrip(",").split(",")]
    list_size = [int(x) for x in row["blockSizes"].rstrip(",").split(",")]
    list_exons = [
        (row["chromStart"] + s, row["chromStart"] + s + z)
        for s, z in zip(list_start, list_size)
    ]
    cds = "".join(
        dna[max(s, row["thickStart"]) - offset : min(e, row["thickEnd"]) - offset]
        for s, e in list_exons
        if min(e, row["thickEnd"]) > max(s, row["thickStart"])
    )
    if row["strand"] == "-":
        cds = str(Seq(cds).reverse_complement())
        list_exons = list_exons[::-1]
    if len(cds) % 3:
        logger.warning(f"{row['name']} CDS length {len(cds)} not a multiple of 3")
        cds = cds[: len(cds) - len(cds) % 3]
    prot = str(Seq(cds).translate())
    if prot.endswith("*"):
        prot = prot[:-1]
    return Transcript(
        # hg19 lift-overs carry a "_N" suffix on the versioned id
        transcript_id=row["name"].split("_")[0],
        build=build,
        chrom=row["chrom"],
        strand=row["strand"],
        list_exons=list_exons,
        cds_start=row["thickStart"],
        cds_end=row["thickEnd"],
        list_tags=[t for t in (row.get("tag") or "").split(",") if t],
        seq_protein=prot,
    )


def fetch_gencode_transcripts(str_gene: str, build: str) -> list[Transcript]:
    """Return a gene's GENCODE protein-coding transcripts on one build.

    Parameters
    ----------
    str_gene : str
        HGNC symbol.
    build : str
        Genome build, a key of :data:`DICT_BUILD2UCSC`.

    Returns
    -------
    list[Transcript]
        Coding transcripts whose ``geneName`` is ``str_gene``.

    Raises
    ------
    ValueError
        If the gene, its transcripts or their sequence cannot be retrieved.
    """
    genome = DICT_BUILD2UCSC[build]
    obj_region = UCSCGeneRegion(genome=genome, gene=str_gene)
    if obj_region.chrom is None:
        raise ValueError(f"{str_gene} not found on {build}")

    obj_track = UCSCTrack(
        genome=genome,
        chrom=obj_region.chrom,
        start=max(0, obj_region.start - INT_REGION_PAD),
        end=obj_region.end + INT_REGION_PAD,
    )
    list_rows = [
        row
        for row in obj_track._rows or []
        if row.get("geneName") == str_gene and row["thickEnd"] > row["thickStart"]
    ]
    if not list_rows:
        raise ValueError(f"No coding {str_gene} transcripts on {build}")

    lo = min(row["chromStart"] for row in list_rows)
    hi = max(row["chromEnd"] for row in list_rows)
    obj_seq = UCSCSequence(genome=genome, chrom=obj_region.chrom, start=lo, end=hi)
    if obj_seq._sequence is None:
        raise ValueError(f"No {build} sequence for {obj_region.chrom}:{lo}-{hi}")
    return [_parse_transcript(row, obj_seq._sequence, lo, build) for row in list_rows]
