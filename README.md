DSRCT vaccine (`dsrct_vax`)
==============================
[//]: # (Badges)
[//]: # (Add codecov, pre-commit.ci, CI and Streamlit App badges once each is live.)

## Intro

`dsrct_vax` supports peptide vaccine design for desmoplastic small round cell tumor
(DSRCT). It finds the four galinpepimut-S (GPS) WT1 vaccine peptides in every UniProt
isoform of WT1 and maps each isoform residue to its genomic codon and exon on GRCh38 and
GRCh37. A Streamlit dashboard shows the results alongside the InterPro domains of the
canonical protein and AlphaFold DB structures.

Peptide sequences are from
[Koehne et al., J Clin Oncol 35(15_suppl):8016 (2017)](https://doi.org/10.1200/JCO.2017.35.15_suppl.8016).
They are defined in `DICT_GPS_PEPTIDES` in `src/dsrct_vax/constants.py`.

| Peptide          | Sequence (wild type)     | Canonical P19544 | Modification |
| :-               | :-                       | :-:              | :-           |
| `WT1-A1`         | `RMFPNAPYL`              | 126–134          | R126Y        |
| `WT1-122A1-long` | `SGQARMFPNAPYLPSCLES`    | 122–140          | R126Y        |
| `WT1-331-long`   | `PGCNKRYFKLSHLQMHSRKHTG` | 328–349          | none         |
| `WT1-427-long`   | `RSDELVRHHNMHQRNMTKL`    | 427–445          | none         |

Positions are 1-based and inclusive. The paper numbers `WT1-331-long` as 331–352, but
on the canonical UniProt isoform it sits at 328–349.

## Getting started

### Installation

Requires Python 3.10–3.12. Install directly from GitHub:

```
pip install "git+https://github.com/jessicaw9910/dsrct_vax.git"
```

For development, clone the repository and make an editable install with the extras you
need (`dev`, `test`, `notebooks`, `app`):

```
git clone https://github.com/jessicaw9910/dsrct_vax.git
cd dsrct_vax
uv venv --python 3.12 .venv
uv pip install -p .venv -e ".[app]"
```

`mkt-schema` and `mkt-databases` are installed from
[missense-kinase-toolkit](https://github.com/choderalab/missense-kinase-toolkit).

### Running the dashboard

```
streamlit run src/dsrct_vax/app/app.py
```

The app reads only the committed data in `src/dsrct_vax/data/`, so it makes no network
calls for data at runtime. Its plotting libraries (BokehJS, 3Dmol.js) load from public
CDNs when reachable; offline, the app embeds them instead (3Dmol.js from a committed
copy in `src/dsrct_vax/app/static/`), so it also runs without internet access. It shows:

- **Sequences:** the vaccine peptides above the wild-type isoforms, lined up by
  genomic codon with exon boundaries, and a track of the canonical protein's InterPro
  domains. A sidebar checkbox switches between GRCh38 and GRCh37.
- **Peptide coverage by isoform:** where each peptide sits in each isoform. A peptide
  can be present, partial (some residues spliced out) or absent.
- **Structure:** the AlphaFold DB model of an isoform, colored by peptide, exon or
  domain.
- **Isoforms, exons and domains:** the transcript matched to each isoform, the coding
  exon coordinates, and each domain's residues, genomic span and member signatures
  (Pfam, PROSITE, SMART).

To deploy on Streamlit Community Cloud, set the entrypoint to `src/dsrct_vax/app/app.py`.
Dependencies are installed from the root `requirements.txt`, which installs `.[app]`.

### Package data

| File                 | Contents                                                       |
| :-                   | :-                                                             |
| `WT1.json`           | Gene snapshot: UniProt isoforms, GENCODE transcripts on both builds, AlphaFold metadata, canonical InterPro domains |
| `PeptideInfo.tar.gz` | One serialized `PeptideInfo` per peptide with its per-isoform matches and overlapping domains |
| `AlphaFold.tar.gz`   | One AlphaFold DB mmCIF per isoform                             |

Regenerate all three by querying UniProt, UCSC (GENCODE `knownGene` on hg38 and hg19),
AlphaFold DB and InterPro (entries of type `domain`, every hit):

```
generate_data --gene WT1
```

This overwrites the files in `src/dsrct_vax/data/`. Review the diff before committing:
apart from `query_date`, the snapshot should only change if a source has been updated.

### Using the data in Python

```python
from dsrct_vax.io_utils import load_gene_info, return_peptide_dict

obj_gene = load_gene_info("WT1")  # GeneInfo: isoforms, codon and exon mappings
dict_peptides = return_peptide_dict()  # {name: PeptideInfo}, read from the archive

obj_pep = dict_peptides["WT1-A1"]
obj_pep.seq_mod  # vaccine sequence with heteroclitic substitutions: YMFPNAPYL
obj_pep.int_start, obj_pep.int_end  # canonical position, 1-based inclusive: 126, 134
for match in obj_pep.dict_isoforms.values():
    print(match.isoform_id, match.status, match.int_start, match.int_end)
obj_pep.list_domains  # canonical InterPro domains it overlaps (IPR000976 for WT1-A1)
```

A new `PeptideInfo(name=...)` built from `DICT_GPS_PEPTIDES` computes its matches from
the committed gene snapshot. One loaded from the archive keeps its stored matches.

### Copyright
Copyright (c) 2026, Jess White

#### Acknowledgements
+ [missense-kinase-toolkit](https://github.com/choderalab/missense-kinase-toolkit) for the API clients and dashboard components
+ [UniProt](https://www.uniprot.org/), [GENCODE](https://www.gencodegenes.org/) via the [UCSC Genome Browser API](https://genome.ucsc.edu/goldenPath/help/api.html), [AlphaFold DB](https://alphafold.ebi.ac.uk/), and [InterPro](https://www.ebi.ac.uk/interpro/)
