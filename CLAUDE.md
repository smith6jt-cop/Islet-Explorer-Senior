# Islet Explorer Senior

## Project scope

This repository contains three related but separate systems:

1. `core/phenocycler/` currently produces UUID-keyed upstream expression
   data from raw QuPath exports. It does not currently produce verified marker
   calls or cell identities.
2. `analysis/xenium/` is the Xenium 5K spatial-transcriptomics lane, absorbed
   2026-08-28 from `Xenium_Analysis@8062e554` (consolidation Phase 3). It carries the
   established two-tier identity-core -> scANVI typer, the panel machinery that reads
   `core/panels`, a dormant spatial-domain benchmark, and a HiPerGator-only `hpc/`
   directory. Runs in `scvi-env`, **not** `xenium_analysis` (broken here). Its numbered
   do-not-regress decisions #1-#17 are cited by number from other documents and must not
   be renumbered. See `analysis/xenium/README.md` and `analysis/xenium/CLAUDE.md`.
3. `apps/islet/` serves aggregated islet, spatial, trajectory, and
   single-cell views. Its shared downstream R code lives once in
   `apps/R/ieCore`; exact display contracts live in
   `apps/islet/config/islet_default.yml` and `apps/lymph/config/lymphoid_default.yml`.
   `apps/lymph` is configuration plus synthetic test data, not another app.

The panel taxonomy has **one** authority: the `core/panels` subtree (= XeniumPanelExplorer
@`076d1192`, verified byte-identical 2026-08-28). `analysis/xenium/scripts/fetch_subpanels.py`
generates a regenerable, gitignored snapshot from it into
`data/development/xenium_panels/subpanels/`, stamped with the source tree id. Never commit a
second copy of those CSVs, and never hand-edit the snapshot.

The source workflow is authoritative for upstream measurement provenance. The
app is a downstream consumer, but its current H5AD is not connected to current
content-addressed source runs. App code and existing H5AD phenotype columns
must not define thresholds, pairs, or train the source workflow.

DOMAIN values such as pixel size, palettes, and phenotype vocabulary are UI
metadata. They are not evidence of cell identity or marker positivity.

## Non-negotiable analysis rules

- Never infer positivity from a desired frequency, cross-donor uniformity, or
  disease ordering.
- There are no verified production target/reference pairs. Do not infer one
  from historical registries, outputs, or expected cell frequencies.
- Existing QuPath classifications, annotations, phenotypes, clusters, disease
  labels, and gold-standard reviews are audit context only. They are not
  training, tuning, or release-gating inputs.
- Preserve the QuPath UUID as `object_id` throughout the upstream workflow.
  Do not substitute centroid matching or donor-name parsing.
- A missing measurement is unavailable, not zero or negative.
- Donors `6457` and `6579` are centrally excluded in
  `core/phenocycler/phenocycler/cohort.py`. Configuration must not bypass
  those exclusions.

## Authoritative Phenocycler workflow

The current source workflow is:

```text
QuPath cohort manifest
  -> ingest
  -> geometry QC
  -> REDSEA boundary-spillover correction
  -> one authoritative intensity per marker
```

Ordinary execution stops after expression. Pair-dependent `controls`,
`calibrate`, `type`, state calling, and assignment export are blocked.

**The content-addressed run framework was removed on 2026-09-16.** Stages write to
stable roots -- `data/cells`, `data/geometry_qc`, `data/redsea`, `data/expression` --
and a donor is done when its partition has content. There is no `run_id`, no stage
manifest, no donor receipt, no fingerprint validation and no `data/runs/<run_id>/`.
A stage whose output is already present is skipped; `--force` reruns it.

Run it **from `core/phenocycler/`**, which owns the only real `config.ini`:

```bash
python -m pip install -e core/phenocycler
cd core/phenocycler
python -m phenocycler.pipeline run --config config.ini --pipelined --through expression
python -m phenocycler.pipeline status --config config.ini
```

There is exactly one `config.ini`, and it is `core/phenocycler/config.ini`. A
repo-root copy existed until 2026-08-10 with no `[restore]` section at all, so
running the command from the repository root silently fell back to the stale
defaults it omitted. If a root `config.ini` ever reappears, it is a bug.

Raw CSV, qptiff, cell/nucleus geometry, pixel calibration, channel mapping,
panel identity, and segmentation version belong in
`core/phenocycler/data/qupath_manifest.json`. That manifest still records a
`content_id` and per-file sha256, but **nothing validates them** -- it is a
description, so it can be hand-edited when a source file moves. Generated results
live under `core/phenocycler/data/<stage>/`. Runs made before 2026-09-16 are still
on disk under `core/phenocycler/data/runs/<run_id>/`; that tree is historical and
nothing writes to it. Active mutable review artifacts live under the repo-root
`data/development/` — a level above `core/`, not inside the package. Only the
ingest/geometry and REDSEA/expression notebooks
describe current stages. Calibration, hierarchical/probabilistic typing,
marker-debug, and Otsu material is historical.

Authoritative documentation:

- `core/phenocycler/README.md`
- `core/phenocycler/docs/WORKFLOW.md`
- `core/phenocycler/DATA_README.md`
- `core/phenocycler/docs/METHOD_HISTORY.md`
- `core/phenocycler/phenocycler/marker_registry.json`
- `core/phenocycler/phenocycler/fundamental_pair_diagnostics.json`
- `docs/WORKFLOW_2026-09-02_MAP.md` — what runs today, lane by lane, with evidence
- `docs/PLAN_2026-09-02_AUXILIARY_PACKAGES.md` — the planned workflow and where every
  auxiliary package sits in it
- `docs/EXTERNAL_PACKAGES.md` — the pin index for external code
- `docs/SOURCE_RECONCILIATION_2026-09-04.md` — source-repository tips, unique-change
  dispositions, and the archive/release gate

The mutable diagnostic roster contains only raw CD3e/Ker8_18 and exact signed
QuPath Iba1_norm/raw Ker8_18. CD79a is excluded as a pair candidate because its
expression is too sparse in normal tissue. These are rank/overlap comparisons,
not production pairs or calls.

Retired, and not to be revived as compatibility fallbacks: RESTORE pair webs,
hierarchical/simultaneous/probabilistic typing, marker-debug/oracle and Otsu
lanes, the Gate 2-5 apparatus, and flat `cells_redsea*`/`restore*` output roots.
Preserve their history in Git.

Note the flat `data/restore_gated_redsea/` and `data/phenotype/broad/` roots at
the repository root are still read by the Shiny app's Card C. They are retired
for the source workflow, not deleted for the app.

## Downstream application boundary

The source workflow currently ends at upstream expression. The existing
application aggregation scripts consume a legacy single-cell H5AD and are not a
supported bridge from a current Phenocycler run.

The next downstream implementation must make one canonical adapter:

```text
cells + expression + marker_evidence + assignments + cell_states
  -> single-cell H5AD
  -> islet, spatial, trajectory, and application aggregates
```

That future adapter must join by donor and `object_id` and embed the run ID plus
relevant manifest content IDs. Until a verified pair/typing workflow and its
adapter exist, do not describe the legacy H5AD rebuild as current, and do not
feed provisional or stale labels into application scripts.

## App architecture

- Entrypoint: `apps/islet/app.R`
- Shared constants and state: `apps/islet/R/00_globals.R`
- Data loading: `apps/islet/R/data_loading.R`
- Plot modules: `apps/islet/R/mod_plot_*.R`
- Drill-down and segmentation helpers:
  `apps/islet/R/drilldown_helpers.R`,
  `apps/islet/R/segmentation_helpers.R`
- Spatial helpers: `apps/islet/R/spatial_helpers.R`
- Viewer helpers: `apps/islet/R/viewer_helpers.R`

## Viewer backend (Avivator | TissUUmaps, Jul 2026)

The Viewer tab can render through **Avivator** (default) or **TissUUmaps** — flip with
`ISLET_VIEWER_BACKEND=tissuumaps`. TissUUmaps adds region drawing, GeoJSON import/export and
per-cell marker overlays, and tiles images itself (no HTTP-Range requirement, so the Viewer
also works under `shiny::runApp`). It consumes **`.tmap` projects** built per donor by
`scripts/senior/build_tissuumaps_project.py`, which names each channel layer from the slide's
own OME-XML — required because the cohort mixes panels (batch-1 SST@35 vs batch-2 SST@1).
The server we run is **our fork, `smith6jt-cop/TissUUmaps`**, built locally by
`deploy/build_tissuumaps_image.sh` (no registry pull, no auto-update).
Findings, costs and the deployment recipe: [`docs/tissuumaps_evaluation.md`](docs/tissuumaps_evaluation.md);
app wiring in `apps/islet/CLAUDE.md` § "Viewer Tab — pluggable backend".

App-specific UI, deployment, browser, plotting, and reticulate guidance lives
in `apps/islet/CLAUDE.md`. Script-level guidance lives in
`scripts/CLAUDE.md`.

Important app conventions:

- Donor status order is `ND < Aab+ < T1D`.
- Donor colors are centralized in `00_globals.R`.
- Use `Sex` in user-facing labels; the stored column remains `gender`.
- Use `scales::pseudo_log_trans(base = 10)` when zeros must remain visible.
- Core parents are `Islet_N`; peri-islet parents are `Islet_N_exp20um`.
- Check `total_cells_peri > 0` before using peri-islet metrics.
- Keep root-level outputs guarded by the active tab and keep module inputs
  correctly namespaced.
- Build downloadable ggplot objects once and share them between rendering and
  `downloadHandler`.

## Running the app

Development:

```bash
cd apps/islet
Rscript -e 'shiny::runApp(".", port = 7777)'
```

Production is served through nginx and shiny-server. See
`deploy/DEPLOY.md` and `apps/islet/CLAUDE.md`; do not treat the development
server as a production deployment test.

## Source-control hygiene

- There are no submodules. `Phenocycler_Analysis` and `XeniumPanelExplorer`
  became the `core/phenocycler` and `core/panels` subtrees in Phase 2, so the
  old "commit the submodule first, then the pointer" rule is retired. One
  commit now covers a change and its tests.
- Do not mix generated `data/`, run manifests, notebooks with outputs, caches,
  or local model/image paths into Git.
- When replacing a workflow, remove obsolete live entrypoints and rely on Git
  history rather than leaving contradictory runnable scripts.
