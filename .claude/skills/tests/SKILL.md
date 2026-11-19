---
name: tests
description: >-
  dsrct_vax testing; extends my central `tests` skill. Use when adding or running tests.
---

# dsrct_vax tests

## Baseline — fetch first

Apply my canonical `tests` conventions (pytest layout, network marker, session fixtures, graceful skips on API failures) before the repo-specific
notes below. WebFetch and follow:

https://raw.githubusercontent.com/jessicaw9910/skills/main/.claude/skills/tests/SKILL.md

If the fetch fails (no network / non-200), **tell me the central `tests` skill
could not be retrieved** and confirm how to proceed — do not silently skip the
baseline.

## Repo-specific additions

- Tests live in `tests/` (configured in `pyproject.toml`).
- The API clients (`gencode.py`, `uniprot.py`, AlphaFold via mkt) hit UCSC,
  UniProt and AlphaFold DB — mark those `network`.
- Schema/peptide/app logic can be tested offline against the committed snapshot
  in `src/dsrct_vax/data/` (e.g. WT1-A1 at canonical 126–134, absent from
  P19544-6/-9); the dashboard can be smoke-tested with
  `streamlit.testing.v1.AppTest`.
