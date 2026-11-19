---
name: versioning
description: >-
  dsrct_vax versioning; defers to my central `versioning` skill. Use when tagging, reading a +distance version, or touching versioningit config.
---

# dsrct_vax versioning

## Baseline — fetch first

Apply my canonical `versioning` conventions (versioningit git-tag versions, importlib.metadata at runtime, PEP 440 local versions) before the repo-specific
notes below. WebFetch and follow:

https://raw.githubusercontent.com/jessicaw9910/skills/main/.claude/skills/versioning/SKILL.md

If the fetch fails (no network / non-200), **tell me the central `versioning` skill
could not be retrieved** and confirm how to proceed — do not silently skip the
baseline.

## Repo-specific additions

- `__version__` is read from installed metadata in `dsrct_vax/__init__.py`;
  the editable `.venv` install goes stale until reinstalled.
