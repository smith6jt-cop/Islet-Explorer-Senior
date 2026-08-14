# Islet Explorer Senior

## Project scope

This repository contains two related but separate systems:

1. `Phenocycler_Analysis/` currently produces UUID-keyed upstream expression
   data from raw QuPath exports. It does not currently produce verified marker
   calls or cell identities.
2. `apps/islet/` serves aggregated islet, spatial, trajectory, and
   single-cell views.

The source workflow is authoritative for upstream measurement provenance. The
app is a downstream consumer, but its current H5AD is not connected to current
content-addressed source runs. App code and existing H5AD phenotype columns
must not define thresholds, pairs, or train the source workflow.

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
  `Phenocycler_Analysis/phenocycler/cohort.py`. Configuration must not bypass
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
`calibrate`, `type`, state calling, and assignment export are blocked. The
run/stage manifest architecture remains authoritative for upstream stages.

Run it **from `Phenocycler_Analysis/`**, which owns the only real `config.ini`:

```bash
python -m pip install -e Phenocycler_Analysis
cd Phenocycler_Analysis
python -m phenocycler.pipeline run --config config.ini --pipelined --through expression
python -m phenocycler.pipeline status --config config.ini
```

There is exactly one `config.ini`, and it is the submodule's. A repo-root copy
existed until 2026-08-10 with no `[restore]` section at all, so running the
command from the repository root silently fell back to the stale defaults it
omitted. If a root `config.ini` ever reappears, it is a bug.

Raw CSV, qptiff, cell/nucleus geometry, pixel calibration, channel mapping,
panel identity, and segmentation version belong in
`Phenocycler_Analysis/data/qupath_manifest.json`. Generated results live under
`Phenocycler_Analysis/data/runs/<run_id>/`; stages and runs are immutable and
content addressed. Active mutable review artifacts live under
`Phenocycler_Analysis/../data/development/`; they never write a stage manifest
or mutate `data/runs/`. Only the ingest/geometry and REDSEA/expression notebooks
describe current stages. Calibration, hierarchical/probabilistic typing,
marker-debug, and Otsu material is historical.

Authoritative documentation:

- `Phenocycler_Analysis/README.md`
- `Phenocycler_Analysis/docs/WORKFLOW.md`
- `Phenocycler_Analysis/DATA_README.md`
- `Phenocycler_Analysis/docs/METHOD_HISTORY.md`
- `Phenocycler_Analysis/phenocycler/marker_registry.json`
- `Phenocycler_Analysis/phenocycler/fundamental_pair_diagnostics.json`

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

The source workflow currently ends at immutable upstream expression. The
existing application aggregation scripts consume a legacy single-cell H5AD and
are not a supported bridge from a current content-addressed Phenocycler run.

The next downstream implementation must make one canonical adapter:

```text
<future verified run>/cells + expression + marker_evidence + assignments + cell_states
  -> run-aware single-cell H5AD
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

- Commit the `Phenocycler_Analysis` submodule first, push it, then commit the
  parent submodule pointer.
- Do not mix generated `data/`, run manifests, notebooks with outputs, caches,
  or local model/image paths into Git.
- When replacing a workflow, remove obsolete live entrypoints and rely on Git
  history rather than leaving contradictory runnable scripts.
