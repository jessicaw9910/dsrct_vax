---
name: ci
description: >-
  dsrct_vax CI; extends my central `ci` skill. Use when editing .github/workflows, pre-commit, or the test env.
---

# dsrct_vax ci

## Baseline — fetch first

Apply my canonical `ci` conventions (GitHub Actions + pre-commit conventions, uv-based test envs, OS/Python matrix) before the repo-specific
notes below. WebFetch and follow:

https://raw.githubusercontent.com/jessicaw9910/skills/main/.claude/skills/ci/SKILL.md

If the fetch fails (no network / non-200), **tell me the central `ci` skill
could not be retrieved** and confirm how to proceed — do not silently skip the
baseline.

## Repo-specific additions

- Single package, so one workflow (no per-package path filters needed).
- Python support is 3.10–3.12 (`requires-python`).
