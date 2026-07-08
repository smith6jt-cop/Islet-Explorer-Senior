# Superseded / backup files archived 2026-07-08

Moved here during the `Phenocycler_Analysis` submodule migration + 8-class broad-lineage rerun
(plan: `~/.claude/plans/make-a-plan-to-linked-allen.md`). Git history is preserved (files were
`mv`d; the rename is recorded at commit). Nothing here is part of the live pipeline or app.

## app_pre_modularization/
- `app_original.R` — pre-modularization monolithic Shiny app (5,397 lines); backup only. The live
  entrypoint is `app/shiny_app/app.R` + `app/shiny_app/R/`.
- `app_modular.R` — orphan intermediate from the modularization (0 references across `.R`/`.md`;
  verified before archiving).

## notebooks_superseded/
- `rebuild_trajectory.ipynb` — superseded by `scripts/reaggregate_islets.py` (per `scripts/CLAUDE.md`).
- `SCVI_CODEX_copy.ipynb` — a small stub "copy"; the real scVI notebook is
  `single_cell_analysis/SCVI_CODEX_v2.ipynb`.

## misc/
- `FileReview.txt` — stale "Phase 2" implementation-summary scratch note.

## Also deleted (not archived — regenerable or self-flagged)
- `HANDOVER.md` — self-flagged *"Delete after reading; do not commit"* (superseded by
  `docs/HANDOFF_22donor_rebuild.md`).
- `notebooks/senior_phenotyping_redsea_restore_executed.ipynb` — gitignored executed-notebook output.
- `app/shiny_app/Channel_names.junior35.bak` — gitignored `.bak` of the old 35-plex Junior channel panel.

## Deferred (NOT done here)
- Retiring the superseded Senior step-1–5 scripts (`redsea_full.py`, `restore_normalize.py`,
  `apply_hormone_floor.py`, `assign_broad_lineage.py`, `marker_taxonomy.py`, `plot_celltype_markers.py`,
  `redsea_reassess_diag.py`, `build_cells_parquet.py`, `export_broad_class_for_qupath.py`, the two
  `redsea_*_validate.py`) → do this **after** PR #3 is merged and the submodule is committed, then
  repoint downstream imports.
- The Junior-era top-level `scripts/` pipeline + one-off diagnostics/converters → left in place pending
  explicit confirmation.
