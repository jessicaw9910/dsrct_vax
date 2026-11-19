---
name: formatting
description: >-
  dsrct_vax-specific Python style; extends my central `formatting` skill. Use whenever writing or editing any Python code in this repo.
---

# dsrct_vax formatting

## Baseline — fetch first

Apply my canonical `formatting` conventions (imports, NumPy-style docstrings, trailing-docstring constants, type-prefixed names, lowercase comments, double quotes) before the repo-specific
notes below. WebFetch and follow:

https://raw.githubusercontent.com/jessicaw9910/skills/main/.claude/skills/formatting/SKILL.md

If the fetch fails (no network / non-200), **tell me the central `formatting` skill
could not be retrieved** and confirm how to proceed — do not silently skip the
baseline.

## Repo-specific additions

- Local imports come from the `dsrct_vax.*` namespace
  (e.g. `from dsrct_vax.schema import PeptideInfo`).
- `PeptideInfo.populate_derived` imports `load_gene_info` inside the method on
  purpose: `io_utils` imports `schema`, so a top-level import would be circular.
- Dashboard colors live in `dsrct_vax/app/constants.py`; peptide colors are fixed
  by `DICT_GPS_PEPTIDES` order so a filter never repaints the others.
- Public residue positions are **1-based inclusive**; 0-based indices stay
  internal (genomic coordinates are stored 0-based half-open, UCSC style, and
  reported 1-based).
