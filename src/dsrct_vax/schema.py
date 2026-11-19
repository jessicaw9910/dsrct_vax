"""Schema for vaccine peptides and the isoforms of their parent gene.

:class:`GeneInfo` holds a snapshot of a gene's UniProt isoforms, its GENCODE
transcripts on each genome build and its AlphaFold models, and on instantiation
maps every isoform residue to a genomic codon and exon (:class:`IsoformInfo`);
transcripts are :class:`dsrct_vax.gencode.Transcript`.
:class:`PeptideInfo` is populated from :data:`dsrct_vax.constants.DICT_GPS_PEPTIDES`
by name and locates the peptide in every isoform (:class:`PeptideMatch`).

All public residue positions are 1-based and inclusive; genomic coordinates are
stored 0-based half-open (UCSC convention) and reported 1-based.
"""

import logging
from enum import StrEnum

from mkt.databases.aligners import BL2UniProtAligner, Kincore2UniProtAligner
from mkt.databases.isoform import map_positions_by_alignment
from pydantic import BaseModel, Field, model_validator

from dsrct_vax.constants import DICT_GPS_PEPTIDES
from dsrct_vax.gencode import (
    DICT_BUILD2UCSC,
    STR_MANE_TAG,
    Transcript,
    return_mane_transcript,
)

logger = logging.getLogger(__name__)

FLOAT_INSERT_STEP = 1e-3
"""float: Key step for isoform residues with no genomic codon (projected isoforms),
placed between their aligned neighbours."""


class MatchMethod(StrEnum):
    """How an isoform was tied to a genomic transcript."""

    exact = "exact"
    """A transcript translates to the isoform (a non-AUG start may differ at residue 1)."""
    alt_start = "alt_start"
    """The isoform is a transcript's translation from a downstream in-frame start."""
    projected = "projected"
    """No transcript encodes it; residues are mapped by alignment to another isoform."""
    none = "none"
    """No transcript or alignable isoform."""


class PeptideStatus(StrEnum):
    """Whether an isoform contains a peptide."""

    present = "present"
    """The whole peptide occurs in the isoform."""
    partial = "partial"
    """Some peptide residues occur at the same genomic codons; the rest are spliced out."""
    absent = "absent"
    """No peptide residue occurs in the isoform."""


class AlphaFoldModel(BaseModel):
    """An AlphaFold DB model of one isoform, stored as package data."""

    entry_id: str
    """AlphaFold DB entry (e.g. ``"AF-P19544-2-F1"``)."""
    version: int | None = None
    """AlphaFold DB model version."""
    model_created: str | None = None
    """Model creation date reported by AlphaFold DB."""
    file_name: str
    """mmCIF file name inside :data:`dsrct_vax.constants.PATH_ALPHAFOLD_TAR`."""
    bool_seq_match: bool
    """Whether the model sequence equals the UniProt isoform sequence."""


class ProteinDomain(BaseModel):
    """One hit of an InterPro or Pfam entry on the canonical isoform."""

    accession: str
    """Entry accession (e.g. ``"IPR013087"`` or ``"PF00096"``)."""
    name: str
    """Entry name (e.g. ``"Zinc finger C2H2-type"``)."""
    type: str
    """Entry type (e.g. ``"domain"``, ``"family"``)."""
    start: int
    """First canonical residue, 1-based inclusive."""
    end: int
    """Last canonical residue, 1-based inclusive."""
    dict_members: dict[str, list[str]] = {}
    """Member database -> signatures integrated in the entry (e.g. ``"pfam"``)."""


class IsoformMapping(BaseModel):
    """An isoform's residues mapped to genomic codons and exons on one build."""

    build: str
    """Genome build."""
    transcript_id: str | None
    """Transcript the mapping comes from (for ``projected``, the reference isoform's)."""
    method: MatchMethod
    """How the isoform was tied to the transcript."""
    int_offset: int = 0
    """Transcript residue = isoform residue + offset (``alt_start`` only)."""
    dict_idx2key: dict[int, float] = {}
    """1-based residue -> genomic position of its codon's middle base (fractional for
    residues with no codon); orders residues across isoforms."""
    dict_idx2exon: dict[int, str] = {}
    """1-based residue -> exon label (numbered on the MANE transcript)."""
    dict_idx2codon: dict[int, tuple[int, int, int]] = {}
    """1-based residue -> 1-based genomic positions of its codon's bases, in coding order."""


class IsoformInfo(BaseModel):
    """A UniProt isoform with its genomic mapping on each build."""

    isoform_id: str
    """UniProt isoform accession (the canonical isoform has no ``-N`` suffix)."""
    seq: str
    """Isoform sequence."""
    bool_canonical: bool
    """Whether this is the UniProt canonical (displayed) isoform."""
    list_ensembl: list[str] = []
    """Unversioned Ensembl transcripts UniProt cross-references for this isoform."""
    dict_build: dict[str, IsoformMapping] = {}
    """Build -> genomic mapping."""
    alphafold: AlphaFoldModel | None = None
    """AlphaFold DB model, if one exists."""


class GeneInfo(BaseModel):
    """Snapshot of a gene's isoforms, transcripts and structures.

    The raw inputs are serialized; ``dict_isoforms`` and ``dict_exon_label`` are
    derived on instantiation.
    """

    gene: str
    """HGNC symbol."""
    uniprot_id: str
    """UniProt accession of the canonical isoform."""
    query_date: str
    """ISO date the snapshot was queried."""
    dict_isoform_seq: dict[str, str]
    """Isoform accession -> sequence (canonical first)."""
    dict_isoform_xref: dict[str, list[str]] = {}
    """Isoform accession -> unversioned Ensembl transcripts from UniProt."""
    dict_transcripts: dict[str, list[Transcript]]
    """Build -> GENCODE protein-coding transcripts of the gene."""
    dict_alphafold: dict[str, AlphaFoldModel] = {}
    """Isoform accession -> AlphaFold DB model."""
    list_domains: list[ProteinDomain] = []
    """InterPro domain hits on the canonical isoform, sorted by start."""
    dict_isoforms: dict[str, IsoformInfo] = Field(default_factory=dict, exclude=True)
    """Isoform accession -> derived :class:`IsoformInfo`."""
    dict_exon_label: dict[str, dict[tuple[int, int], str]] = Field(
        default_factory=dict, exclude=True
    )
    """Build -> genomic exon (start, end) -> label."""

    @property
    def list_builds(self) -> list[str]:
        """Builds in snapshot order; the first orders residues across isoforms."""
        return [b for b in DICT_BUILD2UCSC if b in self.dict_transcripts]

    @property
    def str_canonical(self) -> str:
        """Canonical isoform accession."""
        return self.uniprot_id

    @property
    def strand(self) -> str:
        """Strand of the gene on the reference build."""
        return self.dict_transcripts[self.list_builds[0]][0].strand

    @model_validator(mode="after")
    def populate_isoforms(self) -> "GeneInfo":
        """Derive exon labels and every isoform's genomic mapping."""
        for build in self.list_builds:
            self.dict_exon_label[build] = self._label_exons(build)
        self.dict_isoforms = {
            iso: IsoformInfo(
                isoform_id=iso,
                seq=seq,
                bool_canonical=iso == self.uniprot_id,
                list_ensembl=self.dict_isoform_xref.get(iso, []),
                alphafold=self.dict_alphafold.get(iso),
            )
            for iso, seq in self.dict_isoform_seq.items()
        }
        for build in self.list_builds:
            # transcript-backed isoforms first, then project the rest onto them
            for obj_iso in self.dict_isoforms.values():
                mapping = self._map_to_transcript(obj_iso, build)
                if mapping is not None:
                    obj_iso.dict_build[build] = mapping
            for obj_iso in self.dict_isoforms.values():
                if build not in obj_iso.dict_build:
                    obj_iso.dict_build[build] = self._project(obj_iso, build)
        return self

    def _label_exons(self, build: str) -> dict[tuple[int, int], str]:
        """Label every exon on a build by the MANE transcript's exon numbering.

        An exon overlapping MANE exon ``n`` is labelled ``"n"``; one overlapping none
        is labelled by the MANE exons upstream of it plus a letter (``"0a"`` lies
        upstream of exon 1, ``"4a"`` between exons 4 and 5).
        """
        list_tx = self.dict_transcripts[build]
        mane = return_mane_transcript(list_tx)
        dict_label: dict[tuple[int, int], str] = {}
        dict_novel: dict[int, list[tuple[int, int]]] = {}
        for tx in list_tx:
            for exon in tx.list_exons:
                if exon in dict_label or any(exon in v for v in dict_novel.values()):
                    continue
                list_rank = [
                    rank
                    for rank, (s, e) in enumerate(mane.list_exons, start=1)
                    if min(e, exon[1]) > max(s, exon[0])
                ]
                if list_rank:
                    dict_label[exon] = "/".join(str(r) for r in list_rank)
                    continue
                if mane.strand == "+":
                    n_up = sum(1 for _, e in mane.list_exons if e <= exon[0])
                else:
                    n_up = sum(1 for s, _ in mane.list_exons if s >= exon[1])
                dict_novel.setdefault(n_up, []).append(exon)
        for n_up, list_exon in dict_novel.items():
            list_exon.sort(key=lambda x: x[0], reverse=mane.strand == "-")
            for i, exon in enumerate(list_exon):
                dict_label[exon] = f"{n_up}{chr(ord('a') + i)}"
        return dict_label

    def _map_to_transcript(
        self, obj_iso: IsoformInfo, build: str
    ) -> IsoformMapping | None:
        """Tie an isoform to the transcript that encodes it, if any."""

        def priority(tx: Transcript) -> tuple:
            return (
                tx.id_unversioned not in obj_iso.list_ensembl,
                STR_MANE_TAG not in tx.list_tags,
                "CCDS" not in tx.list_tags,
                tx.transcript_id,
            )

        seq = obj_iso.seq
        list_tx = sorted(self.dict_transcripts[build], key=priority)
        for tx in list_tx:
            prot = tx.seq_protein
            if len(prot) == len(seq) and prot[1:] == seq[1:]:
                return self._build_mapping(tx, build, MatchMethod.exact, 0, len(seq))
        for tx in list_tx:
            prot = tx.seq_protein
            if len(prot) > len(seq) and prot.endswith(seq):
                offset = len(prot) - len(seq)
                return self._build_mapping(
                    tx, build, MatchMethod.alt_start, offset, len(seq)
                )
        return None

    def _build_mapping(
        self,
        tx: Transcript,
        build: str,
        method: MatchMethod,
        offset: int,
        n_res: int,
    ) -> IsoformMapping:
        """Build an isoform's mapping from a transcript's residue codons."""
        list_codon = tx.return_residue_codons()
        dict_label = self.dict_exon_label[build]
        dict_key, dict_exon, dict_codon = {}, {}, {}
        for idx in range(1, n_res + 1):
            mid, rank, tuple_pos = list_codon[idx - 1 + offset]
            dict_key[idx] = float(mid)
            dict_exon[idx] = dict_label[tx.list_exons[rank - 1]]
            dict_codon[idx] = tuple(pos + 1 for pos in tuple_pos)
        return IsoformMapping(
            build=build,
            transcript_id=tx.transcript_id,
            method=method,
            int_offset=offset,
            dict_idx2key=dict_key,
            dict_idx2exon=dict_exon,
            dict_idx2codon=dict_codon,
        )

    def _project(self, obj_iso: IsoformInfo, build: str) -> IsoformMapping:
        """Map an isoform with no transcript by aligning it to a mapped isoform.

        The reference is the canonical isoform if it is transcript-backed, else the
        first that is. Aligned residues inherit the reference's codon and exon;
        unaligned ones get fractional keys between their neighbours and no exon.
        """
        list_ref = [
            o
            for o in self.dict_isoforms.values()
            if build in o.dict_build
            and o.dict_build[build].method in (MatchMethod.exact, MatchMethod.alt_start)
        ]
        if not list_ref:
            logger.warning(f"{obj_iso.isoform_id}: no transcript-backed isoform on {build}")
            return IsoformMapping(build=build, transcript_id=None, method=MatchMethod.none)
        ref = next((o for o in list_ref if o.bool_canonical), list_ref[0])
        map_ref = ref.dict_build[build]
        dict_map = map_positions_by_alignment(obj_iso.seq, ref.seq, BL2UniProtAligner())
        logger.warning(
            f"{obj_iso.isoform_id}: no transcript on {build}; projected onto "
            f"{ref.isoform_id} ({len(dict_map)}/{len(obj_iso.seq)} residues aligned)"
        )
        int_dir = 1 if self.strand == "+" else -1
        dict_key, dict_exon, dict_codon = {}, {}, {}
        for idx, idx_ref in dict_map.items():
            dict_key[idx] = map_ref.dict_idx2key[idx_ref]
            dict_exon[idx] = map_ref.dict_idx2exon[idx_ref]
            dict_codon[idx] = map_ref.dict_idx2codon[idx_ref]
        # unaligned residues: step off the previous aligned key (or back from the next)
        prev_key, n_since = None, 0
        list_pending = []
        for idx in range(1, len(obj_iso.seq) + 1):
            if idx in dict_key:
                for j, idx_pend in enumerate(list_pending, start=1):
                    dict_key[idx_pend] = dict_key[idx] - int_dir * FLOAT_INSERT_STEP * (
                        len(list_pending) - j + 1
                    )
                list_pending, prev_key, n_since = [], dict_key[idx], 0
            elif prev_key is not None:
                n_since += 1
                dict_key[idx] = prev_key + int_dir * FLOAT_INSERT_STEP * n_since
            else:
                list_pending.append(idx)
        return IsoformMapping(
            build=build,
            transcript_id=map_ref.transcript_id,
            method=MatchMethod.projected,
            dict_idx2key=dict_key,
            dict_idx2exon=dict_exon,
            dict_idx2codon=dict_codon,
        )

    def return_genomic_span(
        self, isoform_id: str, build: str, int_start: int, int_end: int
    ) -> str | None:
        """Return the 1-based genomic span of isoform residues ``int_start``-``int_end``.

        Returns
        -------
        str | None
            e.g. ``"chr11:32,392,754-32,392,800; chr11:32,396,258-32,396,320"``, or
            None if a residue has no codon.
        """
        dict_codon = self.dict_isoforms[isoform_id].dict_build[build].dict_idx2codon
        list_codon = [dict_codon.get(i) for i in range(int_start, int_end + 1)]
        if any(c is None for c in list_codon):
            return None
        chrom = self.dict_transcripts[build][0].chrom
        return "; ".join(
            f"{chrom}:{lo:,}-{hi:,}"
            for lo, hi in _return_runs({p for c in list_codon for p in c})
        )

    def return_exon_spans(self, build: str) -> dict[str, list[tuple[int, int]]]:
        """Return exon label -> 1-based coding (start, end) segments on a build.

        Taken from the transcripts that encode the isoforms, so an exon with an
        alternative splice site (e.g. WT1 exon 9 +/-KTS) has one segment per variant.
        """
        dict_label = self.dict_exon_label[build]
        dict_tx = {t.transcript_id: t for t in self.dict_transcripts[build]}
        set_tx = {
            o.dict_build[build].transcript_id
            for o in self.dict_isoforms.values()
            if o.dict_build[build].method in (MatchMethod.exact, MatchMethod.alt_start)
        }
        dict_out: dict[str, set[tuple[int, int]]] = {}
        for tx_id in set_tx:
            tx = dict_tx[tx_id]
            for rank, (lo, hi) in tx.return_coding_exons():
                label = dict_label[tx.list_exons[rank - 1]]
                dict_out.setdefault(label, set()).add((lo + 1, hi))
        return {k: sorted(v) for k, v in dict_out.items()}

    def return_column_keys(self, list_isoforms: list[str] | None = None) -> list[float]:
        """Return the union of residue keys in transcription order (alignment columns)."""
        build = self.list_builds[0]
        list_isoforms = list_isoforms or list(self.dict_isoforms)
        set_key = {
            key
            for iso in list_isoforms
            for key in self.dict_isoforms[iso].dict_build[build].dict_idx2key.values()
        }
        return sorted(set_key, reverse=self.strand == "-")


def _return_runs(set_pos: set[int]) -> list[tuple[int, int]]:
    """Collapse genomic positions into contiguous (start, end) runs, ascending."""
    list_run: list[list[int]] = []
    for pos in sorted(set_pos):
        if list_run and pos == list_run[-1][1] + 1:
            list_run[-1][1] = pos
        else:
            list_run.append([pos, pos])
    return [(lo, hi) for lo, hi in list_run]


def _align_peptide(seq_pep: str, seq_prot: str) -> tuple[int, int]:
    """Align a peptide to a protein; return 0-based half-open (start, end)."""
    alignments = Kincore2UniProtAligner().align(seq1=seq_pep, seq2=seq_prot)
    if len(alignments) != 1:
        raise ValueError(f"Expected 1 alignment of {seq_pep}, found {len(alignments)}.")
    start, end = (int(i) for i in alignments[0].aligned[1][0])
    if end - start != len(seq_pep):
        raise ValueError(f"Alignment of {seq_pep} does not span the full peptide.")
    return start, end


def get_start(seq_pep: str, seq_prot: str) -> int:
    """Return the 1-based position of the peptide's first residue in the protein."""
    return _align_peptide(seq_pep, seq_prot)[0] + 1


def get_end(seq_pep: str, seq_prot: str) -> int:
    """Return the 1-based position of the peptide's last residue (inclusive)."""
    return _align_peptide(seq_pep, seq_prot)[1]


def apply_mods(seq_wt: str, dict_mods: dict[int, str]) -> str:
    """Apply 0-indexed residue substitutions to a wild-type sequence."""
    list_seq = list(seq_wt)
    for idx, aa in dict_mods.items():
        if not 0 <= idx < len(list_seq):
            raise ValueError(f"Mod position {idx} out of range for {seq_wt}.")
        list_seq[idx] = aa
    return "".join(list_seq)


class PeptideMatch(BaseModel):
    """Where a peptide occurs in one isoform."""

    isoform_id: str
    """UniProt isoform accession."""
    status: PeptideStatus
    """Whether the isoform contains the peptide."""
    int_start: int | None = None
    """1-based start in the isoform (``present`` only)."""
    int_end: int | None = None
    """1-based inclusive end in the isoform (``present`` only)."""
    int_n_hits: int = 0
    """Number of exact occurrences (the first is reported)."""
    list_idx_present: list[int] = []
    """1-based peptide positions present in the isoform."""
    list_idx_isoform: list[int] = []
    """1-based isoform positions of ``list_idx_present``."""
    float_coverage: float = 0.0
    """Fraction of peptide residues present."""


class PeptideInfo(BaseModel):
    """Vaccine peptide populated from ``DICT_GPS_PEPTIDES`` by ``name``.

    On instantiation it derives the modified sequence, its canonical position and its
    match in every isoform of ``gene`` (from the package snapshot, see
    :func:`dsrct_vax.io_utils.load_gene_info`).

    Example
    -------
    >>> obj = PeptideInfo(name="WT1-A1")
    >>> obj.seq_mod, obj.int_start, obj.int_end
    ('YMFPNAPYL', 126, 134)
    """

    name: str
    """str: key in ``DICT_GPS_PEPTIDES``."""
    gene: str = "WT1"
    """str: HGNC symbol of the parent protein."""
    seq_wt: str | None = None
    """str: wild-type peptide sequence."""
    mhc_class: str | None = None
    """str: MHC class presented on (``"I"``, ``"II"``, or ``"I/II"``)."""
    dict_mods: dict[int, str] = {}
    """dict[int, str]: 0-indexed peptide position -> substituted residue."""
    seq_mod: str | None = None
    """str: peptide sequence with ``dict_mods`` applied."""
    int_start: int | None = None
    """int: 1-based start of ``seq_wt`` in the canonical UniProt sequence."""
    int_end: int | None = None
    """int: 1-based inclusive end of ``seq_wt`` in the canonical UniProt sequence."""
    dict_isoforms: dict[str, PeptideMatch] = {}
    """dict[str, PeptideMatch]: isoform accession -> match."""
    list_domains: list[ProteinDomain] = []
    """list[ProteinDomain]: canonical InterPro domains overlapping
    ``int_start``-``int_end``."""

    @model_validator(mode="before")
    @classmethod
    def populate_from_dict(cls, data: dict) -> dict:
        """Fill fields not explicitly passed from ``DICT_GPS_PEPTIDES[name]``."""
        str_name = data.get("name")
        if str_name not in DICT_GPS_PEPTIDES:
            raise ValueError(
                f"{str_name} not in DICT_GPS_PEPTIDES: {list(DICT_GPS_PEPTIDES)}"
            )
        return {**DICT_GPS_PEPTIDES[str_name], **data}

    @model_validator(mode="after")
    def populate_derived(self) -> "PeptideInfo":
        """Derive ``seq_mod``, canonical position, isoform matches and domains.

        An object deserialized from ``PeptideInfo.tar.gz`` already carries its isoform
        matches and is returned as stored.
        """
        from dsrct_vax.io_utils import load_gene_info

        if self.dict_isoforms:
            return self

        self.seq_mod = apply_mods(self.seq_wt, self.dict_mods)
        obj_gene = load_gene_info(self.gene)
        seq_canonical = obj_gene.dict_isoform_seq[obj_gene.str_canonical]
        if self.seq_wt in seq_canonical:
            self.int_start = get_start(self.seq_wt, seq_canonical)
            self.int_end = get_end(self.seq_wt, seq_canonical)
        else:
            logger.warning(f"{self.name} is not in canonical {obj_gene.str_canonical}")
        self.dict_isoforms = self.match_isoforms(obj_gene)
        self.list_domains = self.return_domain_overlap(obj_gene)
        return self

    def return_domain_overlap(self, obj_gene: GeneInfo) -> list[ProteinDomain]:
        """Return the canonical domains the peptide overlaps.

        Parameters:
        -----------
        obj_gene : GeneInfo
            Gene snapshot with ``list_domains`` in canonical positions.

        Returns:
        --------
        list[ProteinDomain]
            Domains sharing at least one residue with ``int_start``-``int_end``;
            empty if the peptide is not in the canonical isoform.
        """
        if self.int_start is None or self.int_end is None:
            return []
        return [
            dom
            for dom in obj_gene.list_domains
            if dom.start <= self.int_end and dom.end >= self.int_start
        ]

    def match_isoforms(self, obj_gene: GeneInfo) -> dict[str, PeptideMatch]:
        """Locate the peptide in every isoform.

        An exact occurrence is ``present``. Otherwise the peptide's residues are taken
        at their genomic codons in a reference isoform that contains it (canonical
        first) and looked up in this isoform: residues found at the same codon with
        the same amino acid are ``partial`` coverage (e.g. a peptide spanning a
        skipped exon); none found is ``absent``.
        """
        seq = self.seq_wt
        build = obj_gene.list_builds[0]
        list_ref = sorted(
            (o for o in obj_gene.dict_isoforms.values() if seq in o.seq),
            key=lambda o: not o.bool_canonical,
        )
        list_key = None
        if list_ref:
            ref = list_ref[0]
            idx0 = ref.seq.find(seq)
            dict_key = ref.dict_build[build].dict_idx2key
            list_key = [dict_key.get(idx0 + 1 + i) for i in range(len(seq))]
        else:
            logger.warning(f"{self.name} occurs in no {self.gene} isoform")

        dict_out = {}
        for iso, obj_iso in obj_gene.dict_isoforms.items():
            n_hits = obj_iso.seq.count(seq)
            if n_hits:
                start = obj_iso.seq.find(seq) + 1
                dict_out[iso] = PeptideMatch(
                    isoform_id=iso,
                    status=PeptideStatus.present,
                    int_start=start,
                    int_end=start + len(seq) - 1,
                    int_n_hits=n_hits,
                    list_idx_present=list(range(1, len(seq) + 1)),
                    list_idx_isoform=list(range(start, start + len(seq))),
                    float_coverage=1.0,
                )
                if n_hits > 1:
                    logger.warning(f"{self.name} occurs {n_hits}x in {iso}; using first")
                continue
            list_present, list_iso = [], []
            if list_key is not None:
                dict_key2idx = {
                    k: i for i, k in obj_iso.dict_build[build].dict_idx2key.items()
                }
                for i, key in enumerate(list_key, start=1):
                    idx_iso = dict_key2idx.get(key)
                    if idx_iso is not None and obj_iso.seq[idx_iso - 1] == seq[i - 1]:
                        list_present.append(i)
                        list_iso.append(idx_iso)
            dict_out[iso] = PeptideMatch(
                isoform_id=iso,
                status=PeptideStatus.partial if list_present else PeptideStatus.absent,
                list_idx_present=list_present,
                list_idx_isoform=list_iso,
                float_coverage=len(list_present) / len(seq),
            )
        return dict_out
