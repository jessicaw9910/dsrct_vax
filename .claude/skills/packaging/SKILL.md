---
name: packaging
description: >-
  dsrct_vax packaging; extends my central `packaging` skill. Use when changing pyproject.toml, extras, entry points, or what gets installed.
---

# dsrct_vax packaging

## Baseline — fetch first

Apply my canonical `packaging` conventions (setuptools + pyproject layout, dynamic version, optional-dependencies, console scripts) before the repo-specific
notes below. WebFetch and follow:

https://raw.githubusercontent.com/jessicaw9910/skills/main/.claude/skills/packaging/SKILL.md

If the fetch fails (no network / non-200), **tell me the central `packaging` skill
could not be retrieved** and confirm how to proceed — do not silently skip the
baseline.

## Repo-specific additions

- Single distribution (not a monorepo), `src/` layout: import package
  `dsrct_vax` under `src/dsrct_vax/`.
- `mkt-databases`/`mkt-schema` are git dependencies on
  `choderalab/missense-kinase-toolkit` (mirrors `tansey-lab/mkt_impact`).
- Extras: `dev`, `test`, `notebooks` (ipykernel), `app` (Streamlit dashboard).
- Console script: `generate_data` (`dsrct_vax.cli.generate_data:app`).
- Streamlit Community Cloud installs the root `requirements.txt` (`.[app]`);
  entrypoint `src/dsrct_vax/app/app.py`.
