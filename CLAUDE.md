# Islet Explorer Senior

## Project scope

This repository contains two related but separate systems:

1. `Phenocycler_Analysis/` produces per-cell marker evidence and hierarchical
   cell identities from raw QuPath exports.
2. `app/shiny_app/` serves aggregated islet, spatial, trajectory, and
   single-cell views.

The source workflow is authoritative for marker calibration and cell identity.
The app is a downstream consumer, but its current H5AD is not yet connected to
content-addressed source runs. App code and existing H5AD phenotype columns
must not redefine thresholds or train the source workflow.

## Non-negotiable analysis rules

- Never infer positivity from a desired frequency, cross-donor uniformity, or
  disease ordering.
- Marker calibration is donor-local and uses the marker's registered,
  mutually exclusive reference control.
- Existing QuPath classifications, annotations, phenotypes, clusters, disease
  labels, and gold-standard reviews are audit context only. They are not
  training, tuning, or release-gating inputs.
- Preserve the QuPath UUID as `object_id` from ingest through QuPath import.
  Do not substitute centroid matching or donor-name parsing.
- A missing measurement is unavailable, not zero or negative. An invalid
  donor-marker model is indeterminate and must fail closed.
- Keep identity separate from process/state biology. Ki67 and PCNA are
  orthogonal state outputs, not identity evidence.
- Donors `6457` and `6579` are centrally excluded in
  `Phenocycler_Analysis/phenocycler/cohort.py`. Configuration must not bypass
  those exclusions.

## Authoritative Phenocycler workflow

The current source workflow is:

<<<<<<< HEAD
```text
QuPath cohort manifest
  -> ingest
  -> geometry QC
  -> REDSEA boundary-spillover correction
  -> one authoritative intensity per marker
  -> deterministic donor-local reference controls
  -> marker calibration with uncertainty
  -> simultaneous broad typing
  -> parent-constrained subtype typing
  -> orthogonal state annotation
  -> uncertainty-preserving QuPath export
=======
## App Architecture (Modular - Feb 2026)

The Shiny app lives in `app/shiny_app/` and uses a modular architecture:

```
app/shiny_app/
  app.R                        # Entrypoint (~510 lines): UI layout, CSS, server wiring, root-level outputs
  R/
    00_globals.R               # Package loads, constants, feature flags, paths (h5ad_path, master_path)
    data_loading.R             # load_master_auto(), prep_data(), bin_islet_sizes(), H5AD loading
    utils_stats.R              # summary_stats(), per_bin_anova(), per_bin_kendall(), cohens_d(), eta_squared(), pairwise_wilcox()
    utils_safe_join.R          # safe_left_join(), add_islet_key(), compute_diameter_um()
    segmentation_helpers.R     # GeoJSON cache, build_segmentation_base_plot(), render_islet_segmentation_plot()
    drilldown_helpers.R        # PHENOTYPE_COLORS, load_islet_cells(), render_islet_drilldown_plot(), render_drilldown_summary()
    viewer_helpers.R           # Viewer backend (Avivator | TissUUmaps) URL builders, channel config, environment detection
    ai_helpers.R               # OpenAI credential loading, call_openai_chat()
    mod_plot_ui.R / mod_plot_server.R           # Plot tab (scatter, distribution, outliers, segmentation+drilldown panel)
    mod_trajectory_ui.R / mod_trajectory_server.R  # Trajectory tab (UMAP, heatmap, pseudotime, segmentation+drilldown panel)
    mod_viewer_ui.R / mod_viewer_server.R       # Viewer tab (OME-TIFF iframe; ISLET_VIEWER_BACKEND=avivator|tissuumaps)
    mod_statistics_ui.R / mod_statistics_server.R  # Statistics tab (5-section narrative)
    spatial_helpers.R           # Per-donor tissue CSV loader with new.env() caching
    mod_spatial_ui.R / mod_spatial_server.R     # Spatial tab (3-panel: sidebar controls, tissue scatter, Leiden panel)
    mod_ai_assistant_ui.R / mod_ai_assistant_server.R  # AI chat panel
  www/                         # Static assets (images, logos)
  app_original.R               # Backup of pre-modularization monolithic app (5,397 lines)
>>>>>>> origin/main
```

Run it from the repository root after installing the submodule package:

```bash
conda activate scvi-env
python -m pip install -e Phenocycler_Analysis
python -m phenocycler.pipeline run --config config.ini
python -m phenocycler.pipeline status --config config.ini
```

Raw CSV, qptiff, cell/nucleus geometry, pixel calibration, channel mapping,
panel identity, and segmentation version belong in
`data/qupath_manifest.json`. Generated results live under
`data/runs/<run_id>/`; stages and runs are immutable and content addressed.

The four notebooks are review interfaces over the production APIs:

- `01_qupath_ingest_and_geometry_qc.ipynb`
- `02_redsea_spillover.ipynb`
- `03_marker_calibration.ipynb`
- `04_hierarchical_cell_typing.ipynb`

They must not duplicate algorithms or write around stage manifests.

Authoritative documentation:

- `Phenocycler_Analysis/README.md`
- `Phenocycler_Analysis/docs/WORKFLOW.md`
- `Phenocycler_Analysis/DATA_README.md`
- `Phenocycler_Analysis/phenocycler/marker_registry.json`
- `Phenocycler_Analysis/phenocycler/typing_rules.json`

RESTORE, NNMF/GMM/Otsu gate searches, ordered residual lineage, Gate 5
counterfactuals, flat `cells_redsea*`/`restore*` outputs, and the former
anatomy-check release arms are retired. Preserve their history in Git; do not
revive them as compatibility fallbacks.

## Downstream application boundary

The source workflow currently ends at immutable per-cell assignments and
QuPath export. The existing application aggregation scripts consume a legacy
single-cell H5AD and are not yet a supported bridge from a content-addressed
Phenocycler run.

The next downstream implementation must make one canonical adapter:

```text
<run>/cells + expression + marker_evidence + assignments + cell_states
  -> run-aware single-cell H5AD
  -> islet, spatial, trajectory, and application aggregates
```

That adapter must join by donor and `object_id`, retain ambiguous/unavailable
assignments, and embed the run ID plus relevant manifest content IDs. Until it
exists, do not describe the legacy H5AD rebuild as current, and do not feed
provisional or stale labels into the application scripts.

## App architecture

<<<<<<< HEAD
- Entrypoint: `app/shiny_app/app.R`
- Shared constants and state: `app/shiny_app/R/00_globals.R`
- Data loading: `app/shiny_app/R/data_loading.R`
- Plot modules: `app/shiny_app/R/mod_plot_*.R`
- Drill-down and segmentation helpers:
  `app/shiny_app/R/drilldown_helpers.R`,
  `app/shiny_app/R/segmentation_helpers.R`
- Spatial helpers: `app/shiny_app/R/spatial_helpers.R`
- Viewer helpers: `app/shiny_app/R/viewer_helpers.R`
=======
## Viewer backend (Avivator | TissUUmaps, Jul 2026)

The Viewer tab can render through **Avivator** (default) or **TissUUmaps** — flip with
`ISLET_VIEWER_BACKEND=tissuumaps`. TissUUmaps adds region drawing, GeoJSON import/export and
per-cell marker overlays, and tiles images itself (no HTTP-Range requirement, so the Viewer
also works under `shiny::runApp`). It consumes **`.tmap` projects** built per donor by
`scripts/senior/build_tissuumaps_project.py`, which names each channel layer from the slide's
own OME-XML — required because the cohort mixes panels (batch-1 SST@35 vs batch-2 SST@1).
Findings, costs and the deployment recipe: [`docs/tissuumaps_evaluation.md`](docs/tissuumaps_evaluation.md);
app wiring in `app/shiny_app/CLAUDE.md` § "Viewer Tab — pluggable backend".

## Spatial / Statistics / Trajectory / AI tabs
>>>>>>> origin/main

App-specific UI, deployment, browser, plotting, and reticulate guidance lives
in `app/shiny_app/CLAUDE.md`. Script-level guidance lives in
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
cd app/shiny_app
Rscript -e 'shiny::runApp(".", port = 7777)'
```

Production is served through nginx and shiny-server. See
`deploy/DEPLOY.md` and `app/shiny_app/CLAUDE.md`; do not treat the development
server as a production deployment test.

## Source-control hygiene

- Commit the `Phenocycler_Analysis` submodule first, push it, then commit the
  parent submodule pointer.
- Do not mix generated `data/`, run manifests, notebooks with outputs, caches,
  or local model/image paths into Git.
- When replacing a workflow, remove obsolete live entrypoints and rely on Git
  history rather than leaving contradictory runnable scripts.
