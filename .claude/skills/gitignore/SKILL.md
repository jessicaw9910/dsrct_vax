---
name: gitignore
description: >-
  dsrct_vax .gitignore; extends my central `gitignore` skill. Use when creating or extending the .gitignore.
---

# dsrct_vax gitignore

## Baseline — fetch first

Apply my canonical `gitignore` conventions (the MolSSI base plus personal add-ons and the keep-skills .claude rule) before the repo-specific
notes below. WebFetch and follow:

https://raw.githubusercontent.com/jessicaw9910/skills/main/.claude/skills/gitignore/SKILL.md

If the fetch fails (no network / non-200), **tell me the central `gitignore` skill
could not be retrieved** and confirm how to proceed — do not silently skip the
baseline.

## Repo-specific additions

- **Do not ignore `src/dsrct_vax/data/`** — the baseline's `data/` rule would
  drop the committed package data; anchor any analysis-output ignore to the repo
  root (`/data/`) instead.
