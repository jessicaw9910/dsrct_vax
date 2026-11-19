---
name: git
description: >-
  dsrct_vax git conventions; defers to my central `git` skill. Use when creating commits, branches, or pull requests.
---

# dsrct_vax git

## Baseline — fetch first

Apply my canonical `git` conventions (branch first, scoped imperative commits, Claude attribution, PR-template handling) before the repo-specific
notes below. WebFetch and follow:

https://raw.githubusercontent.com/jessicaw9910/skills/main/.claude/skills/git/SKILL.md

If the fetch fails (no network / non-200), **tell me the central `git` skill
could not be retrieved** and confirm how to proceed — do not silently skip the
baseline.

## Repo-specific additions

- Remote is `jessicaw9910/dsrct_vax` (personal GitHub, not tansey-lab); default
  branch `main`.
- Regenerated package data (`src/dsrct_vax/data/*`) goes in its own commit,
  separate from the code change that required it.
