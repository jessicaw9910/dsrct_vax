---
name: release
description: >-
  dsrct_vax releases; defers to my central `release` skill. Use when tagging a release or adding a publish workflow.
---

# dsrct_vax release

## Baseline — fetch first

Apply my canonical `release` conventions (tag-as-version releases, tag-triggered build, OIDC trusted publishing, install-from-git fallback) before the repo-specific
notes below. WebFetch and follow:

https://raw.githubusercontent.com/jessicaw9910/skills/main/.claude/skills/release/SKILL.md

If the fetch fails (no network / non-200), **tell me the central `release` skill
could not be retrieved** and confirm how to proceed — do not silently skip the
baseline.

## Repo-specific additions

- Not on PyPI; install from git
  (`pip install "git+https://github.com/jessicaw9910/dsrct_vax.git"`).
