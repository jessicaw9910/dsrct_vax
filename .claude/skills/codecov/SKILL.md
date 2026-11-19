---
name: codecov
description: >-
  dsrct_vax Codecov config; defers to my central `codecov` skill. Use when adding or tuning .codecov.yml.
---

# dsrct_vax codecov

## Baseline — fetch first

Apply my canonical `codecov` conventions (the .codecov.yml carryforward/status/ignore conventions) before the repo-specific
notes below. WebFetch and follow:

https://raw.githubusercontent.com/jessicaw9910/skills/main/.claude/skills/codecov/SKILL.md

If the fetch fails (no network / non-200), **tell me the central `codecov` skill
could not be retrieved** and confirm how to proceed — do not silently skip the
baseline.

## Repo-specific additions

- Single package: one flag (`dsrct_vax`); ignore `src/dsrct_vax/app/` if the
  dashboard is covered only by smoke tests.
