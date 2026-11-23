# Handoff: MultiPlex5000 is now a standalone producer

Snapshot: 2026-11-23. Supersedes the "pipeline stays in Senior" statements of the
2026-10-20 plan and the 2026-10-28 coordination handoff for new pipeline changes.

## State

- Repository: `smith6jt-cop/MultiPlex5000`, release `multiplex5000-v0.1.0` (commit `a4d7f6a`).
- Senior `main` at `fd114ff` has no `core/multiplex5000` or `core/panels` change that is
  absent from the standalone repository. Senior PRs #110 and #111 were ported unchanged
  (MultiPlex5000#4 and #2) and closed.
- Parity: the installed 0.1.0 wheel replayed donor 6436's multicell assignment (139 per
  weight at 0.25, 0.50, 0.75), identical to the retained 2026-11-07 replay; spleen `06_typing`
  parity is recorded in `MultiPlex5000/docs/MIGRATION.md` at `f993723`.

## Rules from now on

1. New pipeline, stage, typing and display-exporter changes go to MultiPlex5000. Do not
   change `core/multiplex5000` or `core/panels` here; Senior's copy is frozen at `fd114ff`
   until its callers move.
2. Senior owns the Shiny app, the `senior.display/1` reader and UFRC deployment.
3. The unchanged scientific rules in `CLAUDE.md` still apply (donor exclusions, `object_id`,
   no inferred positivity, missing is not zero).

## Remaining work before Senior's copy can be retired

Each caller needs a tested replacement that uses the installed wheel:

- the editable pipeline install in `.github/workflows/xenium-analysis.yml` (`pip install --no-deps -e ../../core/multiplex5000`)
- `.github/workflows/script-tests.yml`: it installs only third-party dependencies; its checkout dependency is the
  `sys.path` insertion of `core/multiplex5000/scripts/typing` in `scripts/senior/build_tissuumaps_regions.py`
  (lines 63-67), which imports `build_islet_geometry`, `cohort` and `subpanel_io` directly. That import migration is separate from an install change
- `scripts/senior/*` layer builders that import typing helpers from the checkout
- the same-slide lane scripts (`core/multiplex5000/scripts/sameslide`)
- a staged `mplx5k-display plan/export/verify` run read by the Senior reader with hash verification
- rerun HiPerGator spleen parity at 0.1.0

Roll back with the previous Senior environment and launcher settings; never with
`git reset --hard` or `git clean -xfd`.
