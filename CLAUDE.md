# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Skills

This repo ships project skills under `.claude/skills/`. **Invoke the relevant
skill before doing that kind of work.**

The skills are thin: they `WebFetch` my canonical baseline conventions from the
central skills repo (`github.com/jessicaw9910/skills`) and layer on only the
repo-specific additions. If a fetch fails, the skill will say so rather than
silently skipping the baseline. Don't copy the central skills into this repo.

- **`formatting`** — before writing or editing any Python.
- **`git`** — before creating commits, branches, or pull requests.
- **`packaging`** — before changing `pyproject.toml`, extras, or entry points.
- **`packaging-data`** — before adding, regenerating, or debugging the committed
  data in `src/dsrct_vax/data/`.
- **`versioning`** — before tagging or touching versioningit config.
- **`tests`** — before adding or running tests.
- **`ci`** / **`codecov`** — before editing `.github/workflows`, pre-commit, or
  `.codecov.yml`.
- **`release`** — before tagging a release or adding a publish workflow.
- **`readme`** — before creating or editing the README.
- **`gitignore`** — before creating or extending the `.gitignore`.

## Project Overview

`dsrct_vax` locates the galinpepimut-S (GPS) WT1 vaccine peptides
(`DICT_GPS_PEPTIDES` in `constants.py`; sequences from
[Koehne et al., JCO 2017](https://doi.org/10.1200/JCO.2017.35.15_suppl.8016)) in every UniProt isoform of WT1, maps each
isoform residue to its genomic codon and exon on GRCh38 and GRCh37, and serves the
results — with AlphaFold DB structures — in a Streamlit dashboard. It builds on
`mkt-databases`/`mkt-schema` from
[missense-kinase-toolkit](https://github.com/choderalab/missense-kinase-toolkit).

## Environment Setup

```bash
uv venv --python 3.12 .venv
uv pip install -p .venv -e ".[app]"   # extras: dev, test, notebooks, app
```

Supported Python versions: 3.10–3.12.

## Key Commands

```bash
generate_data --gene WT1                  # query sources, rewrite src/dsrct_vax/data/*
streamlit run src/dsrct_vax/app/app.py    # dashboard (reads committed data only)
```

Streamlit Community Cloud: entrypoint `src/dsrct_vax/app/app.py`, dependencies from
the root `requirements.txt` (`.[app]`).

## Code Architecture

API clients follow `mkt.databases`: one module per source, each client a
dataclass inheriting `mkt.databases.api_schema.RESTAPIClient` that queries on
instantiation through mkt's cached session and stamps provenance.

```
gencode.py     UCSCGeneRegion / UCSCTrack / UCSCSequence clients, Transcript model,
               fetch_gencode_transcripts (GENCODE knownGene on hg38 + hg19)
uniprot.py     UniProtIsoformFASTA (extends mkt UniProt); Ensembl xrefs via mkt UniProtJSON
interpro.py    Pfam (our version of mkt Pfam: one row per hit, not just the first) and
               InterPro (integrated entries, same parser); return_domains -> ProteinDomain
isoforms.py    assembles a GeneInfo snapshot (UniProt + GENCODE + AlphaFold + InterPro domains)
schema.py      GeneInfo -> IsoformInfo/IsoformMapping (codon + exon per residue);
               PeptideInfo -> PeptideMatch per isoform (present / partial / absent)
               and list_domains (canonical InterPro domains it overlaps)
io_utils.py    load/save the snapshot and the PeptideInfo / AlphaFold archives
cli/           generate_data
app/           Streamlit dashboard; visualizers.py is adapted from the mkt app;
               offline.py inlines BokehJS / the committed static/3Dmol.js when CDNs
               are unreachable
```

### Data flow

1. `generate_data` queries UniProt, UCSC (GENCODE), AlphaFold DB and InterPro, then writes
   `data/<gene>.json`, `data/AlphaFold.tar.gz` and `data/PeptideInfo.tar.gz`.
2. `GeneInfo` ties each isoform to the transcript that translates to it (exact,
   or a downstream in-frame start), else projects it by alignment; exons are
   numbered on the MANE Select transcript (`1b` = an exon outside it).
3. `PeptideInfo(name=...)` derives `seq_mod`, canonical `int_start`/`int_end`
   (1-based inclusive) and per-isoform matches; objects loaded from the archive
   keep their stored matches.
4. The app reads only the committed data — no network at runtime. BokehJS and
   3Dmol.js load from CDNs when reachable, else are embedded (works offline).
