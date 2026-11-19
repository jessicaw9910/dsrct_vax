---
name: packaging-data
description: >-
  dsrct_vax package data; extends my central `packaging-data` skill. Use when adding, regenerating, or debugging the committed data files.
---

# dsrct_vax packaging-data

## Baseline — fetch first

Apply my canonical `packaging-data` conventions (MANIFEST.in graft vs package-data, the dotted-key gotcha, pruning caches) before the repo-specific
notes below. WebFetch and follow:

https://raw.githubusercontent.com/jessicaw9910/skills/main/.claude/skills/packaging-data/SKILL.md

If the fetch fails (no network / non-200), **tell me the central `packaging-data` skill
could not be retrieved** and confirm how to proceed — do not silently skip the
baseline.

## Repo-specific additions

- Package data lives in `src/dsrct_vax/data/` and is **committed** (the hosted
  app reads it offline): `PeptideInfo.tar.gz` (one JSON per peptide, read one at
  a time in memory), `AlphaFold.tar.gz` (one mmCIF per isoform) and
  `<gene>.json` (the `GeneInfo` snapshot).
- Regenerate all three with `generate_data --gene WT1`; never hand-edit them.
  Rerun it whenever `DICT_GPS_PEPTIDES`, the schema, or the sources change —
  archived `PeptideInfo` objects keep their stored matches.
- The `package-data` key is `dsrct_vax` (a regular package, not a namespace).
