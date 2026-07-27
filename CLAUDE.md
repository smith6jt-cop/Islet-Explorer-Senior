# CLAUDE.md - Islet Explorer

This file provides guidance to Claude Code when working with this repository.

## Project Overview

Interactive Shiny app for exploring pancreatic islet CODEX multiplexed imaging data from nPOD donors. Tracks insulin loss and immune infiltration across ND -> Aab+ -> T1D progression. Features islet-level aggregated views, pseudotime trajectory, spatial neighborhood analysis, and single-cell drill-down.

## ⚠️ Working rule — no frequency/uniformity assumptions

**Never assume marker or cell-type frequencies, and never treat per-donor variation as an error to
normalize toward a cohort norm — donor biology is genuinely variable — BUT validate each call against
that donor's OWN data** (its histograms, negative-control population, spatial coherence of the calls),
not against other donors' values or a clinical label. Positivity is **RESTORE-threshold-defined**, not
a raw-MFI floor. Do not call a marker rare/abundant or a fraction "implausible" from assumption —
measure it per-donor. Justify method changes on methodological grounds (stability/representativeness),
not on "this donor should look like the others."

**Hard donor exclusions — 6457 and 6579 are excluded until further notice from every analysis.**
Donor 6457's poor DAPI staining prevented accurate Phenocycler channel registration; donor 6579 is
excluded by maintainer decision as a pancreatitis outlier. Raw/REDSEA files may remain for
provenance, but no RESTORE, phenotype, clustering, trajectory, figure, or app artifact may include
either donor. Enforced fail-closed in `Phenocycler_Analysis/phenocycler/cohort.py` (automatic
discovery filters both; explicit requests naming either raise). **22 acquired → 20 analysis-eligible.**
Note **donor 6476 is NOT a pancreatitis donor** (older docs said so). Its high CD3e (and high CD11b)
is high **channel background**, not extreme true expression (maintainer, 2026-07-24) — RESTORE's divisor
is designed to account for it, and the v9 gate fix ensures the QC gates do not spuriously abstain on it.

**Data note:** the phenotyping pipeline lives in the `Phenocycler_Analysis` submodule — see the
**Status** block below. `data/` + raw inputs are NOT in git (~21 GB, gitignored).

## Sub-directory CLAUDE.md files

- `app/shiny_app/CLAUDE.md` — per-tab feature details (Drill-Down, Trajectory, Plot, AI Assistant) + conventions including the May 2026 font-standardisation system
  - `app/shiny_app/docs/spatial.md` — Spatial Neighborhood Analysis deep reference (split out May 2026)
  - `app/shiny_app/docs/statistics.md` — Statistics Tab deep reference (split out May 2026)
- `scripts/CLAUDE.md` — pipeline scripts, batch correction, pseudotime modes, **the Senior phenotyping pipeline (REDSEA spillover → RESTORE normalization → broad-lineage assignment → QuPath inspection)**, and the **Senior trajectory batch-correction + robustness fix (Jul 2026)** — scVI `batch_key`, Harmony θ-sweep, robust root, hard gates, combined mode
  - (`islet_analysis/` was removed; its content folded into `scripts/senior/` + `scripts/CLAUDE.md`)

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
    viewer_helpers.R           # Avivator URL builders, channel config, environment detection
    ai_helpers.R               # OpenAI credential loading, call_openai_chat()
    mod_plot_ui.R / mod_plot_server.R           # Plot tab (scatter, distribution, outliers, segmentation+drilldown panel)
    mod_trajectory_ui.R / mod_trajectory_server.R  # Trajectory tab (UMAP, heatmap, pseudotime, segmentation+drilldown panel)
    mod_viewer_ui.R / mod_viewer_server.R       # Viewer tab (OME-TIFF, Avivator iframe)
    mod_statistics_ui.R / mod_statistics_server.R  # Statistics tab (5-section narrative)
    spatial_helpers.R           # Per-donor tissue CSV loader with new.env() caching
    mod_spatial_ui.R / mod_spatial_server.R     # Spatial tab (3-panel: sidebar controls, tissue scatter, Leiden panel)
    mod_ai_assistant_ui.R / mod_ai_assistant_server.R  # AI chat panel
  www/                         # Static assets (images, logos)
  app_original.R               # Backup of pre-modularization monolithic app (5,397 lines)
```

### Data Directory Layout (Senior, relocated Jun 2026)

`data/` is now a **real directory at the repo root** (was `~/IO60panc2nd/senior_data`, one level up; moved in and gitignored — never committed; ~21 GB). It holds the raw/intermediate pipeline files at top level (`singlecell_protein.h5ad`, `cells/` parquet, `islet_explorer_senior.h5ad`, `cluster_*`, `marker_audit_*`, …) **plus** a curated **`data/app_data/`** subdirectory that is the app's view (per-islet drill-down CSVs, GeoJSONs, `islet_spatial_lookup.csv`, and relative symlinks `islet_explorer.h5ad`→`../islet_explorer_senior.h5ad`, `adata_ins_root.h5ad`→`../islets_core_senior.h5ad`, `parquet`, `phenotype_rules.csv`). **The app reads from `data/app_data/`** (R paths use `file.path("..","..","data","app_data", …)`); the **Senior pipeline scripts** read/write `data/` top-level via `Path(__file__).resolve().parents[2] / "data"` (move-proof, no absolute paths).

### Data Loading (H5AD with Excel Fallback)

The app uses `load_master_auto()` which tries H5AD first, falls back to Excel:
- `h5ad_path` -> `data/app_data/islet_explorer.h5ad` (Senior; -> `../islet_explorer_senior.h5ad`)
- `master_path` -> `data/app_data/master_results.xlsx` (Excel fallback; absent in Senior)

Both produce the same `list(markers, targets, comp, lgals3, phenotypes, donor_demographics, neighborhood)` structure consumed by `prep_data()`. The `phenotypes`, `donor_demographics`, and `neighborhood` elements are extracted from H5AD `.obs` and are `NULL` when loading from Excel (graceful fallback).

### Performance: Reticulate Bridge Optimization (Mar 2026)

**Critical**: `ad$uns[[key]]` crosses the Python-R bridge on EVERY call (~221ms each). With ~62 groovy columns × 3 sheets = ~186 calls, this took ~16s. Fix: cache `uns <- ad$uns` once (0.2s) then index the R list. `reconstruct_groovy_df_from_list(uns, sheet)` replaces per-key bridge crossings.

**Deferred loading**: `annotations.tsv` (72 MB, 576K rows) was loaded eagerly at startup but only serves as a fallback when `islet_spatial_lookup` doesn't have the islet (never happens). Now lazy-loaded via `.seg_lazy` environment on first fallback access.

Total startup: **~22s → ~3.5s**.

### Shared Reactive State (wired in app.R server)

- `prepared()` - core data reactive from H5AD or Excel, consumed by Plot, Trajectory, Statistics, Spatial
- `selected_islet` - reactiveVal with `case_id`, `islet_key`, `centroid_x/y`, `donor_status`, `donor_age`, `donor_gender`
- `forced_image` - reactiveVal, written by Trajectory, read by Viewer
- `active_tab` - `reactive(input$tabs)`, passed to Plot and Trajectory to prevent duplicate output IDs
- `active_donor_colors` - reactive from `donor_palette_name()`, synced across 3 non-namespaced inputs: `sidebar_donor_palette` (Plot sidebar), `traj_donor_palette` (Trajectory), `spatial_donor_palette` (Spatial)

### Critical: Root-Level Outputs + Active Tab Guard

The following render outputs are defined in **app.R** at the root level, NOT inside any module:
- `output$islet_segmentation_view` - GeoJSON boundary plot (shared by Plot + Trajectory)
- `output$islet_drilldown_view` - Single-cell spatial overlay plot
- `output$islet_drilldown_summary` - Cell composition bar chart
- `output$islet_drilldown_table` - Core/peri cell counts

Both Plot and Trajectory tabs embed these using non-namespaced `plotOutput("islet_segmentation_view")` etc. **IMPORTANT**: Each module's `segmentation_viewer_panel` renderUI must guard with `if (active_tab() != "Plot") return(NULL)` (or `"Trajectory"` respectively). Without this guard, both renderUIs fire when `selected_islet()` changes, creating duplicate DOM IDs — Shiny can only bind one, causing the other to grey out. Shiny `renderUI` fires based on reactive dependencies, NOT tab visibility.

### Non-Namespaced Inputs for Root-Level Outputs

The segmentation panel `renderUI` in both modules generates non-namespaced inputs:
- `drilldown_view_mode` - "Boundaries" or "Single Cells" toggle
- `drilldown_color_by` - "phenotype" or a marker column name
- `drilldown_show_peri` - Whether to show peri-islet cells

These are read by the root-level `renderPlot` outputs in `app.R`.

### Plotly Click Events in Modules

- Plot module: `source = ns("plot_scatter")` with `event_register("plotly_click")`
- Trajectory module: `source = ns("traj_scatter")` with `event_register("plotly_click")`
- Both use namespaced sources; `event_data()` calls must use matching `ns()` prefix
- **macOS native-click bridge (Jun 2026)**: `plotly_click` is suppressed on macOS, so both scatters also pipe `islet_click_bridge(ns("<plot|traj>_native_click"))` (native mousedown/mouseup → `Shiny.setInputValue`) after `event_register`. The new server observers and the legacy plotly_click handler both call the shared `select_islet_from_key()` (`segmentation_helpers.R`). See `app/shiny_app/CLAUDE.md` § "macOS Native-Click Bridge + Shared Resolver".

## Key Data Files

- `data/islet_explorer.h5ad` - **Primary app data** (~70 MB, 5,214 islets, groovy + trajectory + donor + neighborhood + Leiden)
- `data/master_results.xlsx` - Aggregated islet-level data (composition, markers, targets) -- Excel fallback
- `data/adata_ins_root.h5ad` - Trajectory h5ad (5,214 islets, 31 vars, DPT pseudotime, UMAP)
- `data/neighborhood_metrics.csv` - Per-islet peri-islet metrics (5,214 rows, 62 columns)
- `data/cells/*.csv` - Per-islet single-cell CSVs for drill-down (~5,214 files, ~203 MB total)
- `data/donors/*.csv` - Per-donor tissue-wide cell CSVs for Spatial tab scatter (15 files, ~78 MB total)
- `data/islet_spatial_lookup.csv` - Centroid coordinates for segmentation viewer
- `data/json/*.geojson` / `data/gson/*.geojson.gz` - QuPath segmentation boundaries
- `data/DATA_PROVENANCE.md` - Documents canonical H5AD lineage and data sources

## Data Pipeline

### Canonical H5AD Lineage

```
single_cell_analysis/CODEX_scvi_BioCov_phenotyped_newDuctal.h5ad  (2.6M cells, 31 markers)
  |
  |-- scripts/reaggregate_islets.py  (one-stop: aggregate + trajectory + Leiden)
  |     -> islet_analysis/islets_core_fixed.h5ad  (5,214 islets, min_cells=0, require_paired=True)
  |     -> islet_analysis/islets_peri_fixed.h5ad   (5,214 peri-islet)
  |     -> islet_analysis/islets_merged_fixed.h5ad (5,214 paired core+peri)
  |     -> data/adata_ins_root.h5ad  (+ core + combined pseudotimes, UMAP, X_scVI_harmony_pt, X_scVI_combined_harmony_pt)
  |     -> islet_analysis/islets_core_clustered.h5ad  (+ Leiden at 4 resolutions, Harmony θ=2 on imageid)
  |
  |-- scripts/compute_neighborhood_metrics.py
  |     -> data/neighborhood_metrics.csv  (5,214 rows, ~100 cols — per-phenotype enrich_z and min_dist)
  |
  |-- scripts/extract_per_islet_cells.py
  |     -> data/cells/*.csv  (5,214 files, 37 cols: coords, phenotype, region, 31 markers)
  |
  |-- scripts/extract_per_donor_tissue.py
  |     -> data/donors/*.csv  (15 files, 5 cols: X/Y coords, phenotype, cell_region, islet_name)
  |
  +-- scripts/build_h5ad_for_app.py (trajectory + groovy + donor + neighborhood + Leiden)
        -> data/islet_explorer.h5ad  (all app data in one file, incl. leiden_* + leiden_umap_*)
```

`scripts/reaggregate_islets.py` replaces the manual notebook workflow (aggregation + trajectory + Leiden in one script). `islet_analysis/islets_core_clustered.h5ad` provides Leiden clustering (4 resolutions) + UMAP coords, merged by `build_h5ad_for_app.py` step 4.7.

See `scripts/CLAUDE.md` for per-script details, batch correction (Phase 17), and the two-pseudotime architecture (`dpt_pseudotime` core-only vs `dpt_pseudotime_combined` core+peri+structures).

### Senior phenotyping pipeline (current, Jun 2026) — supersedes the provisional auto-labels

The Senior 59-plex data is being **re-phenotyped from scratch** because global clustering produced a large fake "Unassigned" fraction caused by **upstream artifacts**, fixed in order **REDSEA → RESTORE → broad lineage** (full detail in `scripts/CLAUDE.md`):

```
data/cells/donor_id=*/  (raw 59-marker MFI + object_id UUID; 25,151,028 canonical cells after cell QC)
  |-- QuPath cell detection (scripts/groovy/0_cell_detection.groovy — NEVER selectAnnotations(); see the handoff)
  |-- QuPath per-cell GeoJSON export (Phenocycler_Analysis/scripts/groovy/export_cells_geojson.groovy)
  |-- phenocycler.cell_qc                  Tier-1 exclusions BEFORE any parameter estimation
  |     -> 31,454,041 -> 25,151,028 (20.04% removed); pre-QC copies in data/cells_preqc/
  |-- phenocycler.redsea (v10)             pixel-level lateral-spillover correction
  |     (comp_mode=1 subtract+reinforce, norm_form=donor, α=1, 1px band, nuclear/ECM channels uncorrected)
  |     -> data/cells_redsea/donor_id=*/    (59 spillover-corrected means; the r≈0.998 figure is a RASTERISATION-fidelity check on
  |                                            UNcorrected means, NOT raw-vs-corrected agreement — actual r 0.81-0.98)
  |-- RESTORE pair validation (phenocycler.restore_validation)  NOT YET ACCEPTED — see Status
  |     -> manuscript-faithful screen; emits NO normalized values and NO phenotype calls
  +-- lineage (phenocycler.lineage)  every cell typed, 0 Unassigned (Other = explicit residual)
        -> data/phenotype/broad/donor_id=*/  (5 compartments + Other: Immune/Epithelial/Endothelial/Neural/Mesenchymal; endocrine/exocrine deferred to a validated 2nd pass)
        -> QuPath inspection: export_broad_class_for_qupath.py + scripts/groovy/import_broad_lineage.groovy
```

> ### ⚠ START HERE — [`docs/GATE4_HANDOFF_2026-07-27.md`](docs/GATE4_HANDOFF_2026-07-27.md)
>
> **Gate 4 is COMPLETE and clean on all 20 donors** (2026-07-27 evening): 479/480 donor-pairs VALID,
> zero technical failures, 0% Unresolved, 25,151,028 cells typed. The gate set was expanded 7 → 24
> markers on maintainer sign-off and `Other` fell from ~65% to **13.7%**; islet-core cells — endocrine
> by QuPath annotation, hormone-validated — went from **1–6% Epithelial to 78.9%** (60.5–85.5% per
> donor). Policy frozen at `reference_only`/`p99`; `COMPARTMENT_ORDER` is now
> Immune > **Endothelial** > Epithelial > Neural > Mesenchymal.
>
> **Nothing here is verified against ground truth.** The gold-standard bundle (`data/gold_standard_v1/`)
> is built, blinded and waiting on annotator time; the per-image review of the 14 new markers is
> outstanding. Read the handoff for the open items before changing anything.
>
> Background, still authoritative for their own subjects:
> [`docs/ANATOMY_CHECK_2026-07-27.md`](docs/ANATOMY_CHECK_2026-07-27.md) (anatomy evidence, divisor
> sweep) and [`docs/HANDOFF_2026-07-27_qc_and_v10.md`](docs/HANDOFF_2026-07-27_qc_and_v10.md) (cell QC,
> REDSEA v10, duplicate detections)
>
> **The anatomy check has been RUN.** Its finding: the **single-marker Epithelial gate** is a defect on
> its own terms, and it is at least as large a lever as the divisor. Islet-core cells — endocrine by
> QuPath annotation, hormone-validated, unambiguously epithelial — are called `Epithelial` in only
> **1–6%** of cases under the gate and divisor **production is currently wired to**, and `Other` in
> **82–94%**. Changing only the gate to a 4-marker union (E-cadherin ∪ Pan-CK ∪ Ker8/18 ∪ EpCAM) moves
> that to 16–61%; changing only the divisor to `reference_only`/p99 moves it to 11–44%; both give
> 28–72%. Also: `reference_only`/**max is not a viable policy at all** — `restore_apply` halts on it
> (MODEL_UNSTABLE), and it is the arm production has silently defaulted to since Gate 4 shipped
> (now fixed: the v10 knobs are `config.ini` keys). Divisor still **undecided** — decide it *after* the
> Epithelial gate, not before.
>
> **The Status block below is superseded on four points (2026-07-25 → 27). Read the handoff too.**
>
> ⚠ `HANDOFF*` is in `.gitignore` (under "Temporary files"), so
> `docs/HANDOFF_2026-07-27_qc_and_v10.md` is **local-only and will not survive a fresh clone**, despite
> being the authoritative status document. `docs/ANATOMY_CHECK_2026-07-27.md` IS tracked and carries the
> current findings. Decide whether to keep that ignore rule.
>
> 1. **311,340 duplicate cell detections found and removed** in 7 donors (6450, 6505, 6523, 6534,
>    6566, 6591, 6623) — the same physical cell segmented twice with two UUIDs, 100% inside islets,
>    inflating islet counts ~2×. Cause: detection was run with `selectAnnotations()`, which selects
>    Tissue AND every nested `Islet_N`; InstanSeg segments each selected parent independently.
>    Fixed at source by `scripts/groovy/0_cell_detection.groovy`.
> 2. **New cell-QC stage** (`phenocycler/cell_qc.py`), applied: 31,454,041 → **25,151,028** canonical
>    cells (pre-QC copies in `data/cells_preqc/`). Runs before any parameter estimation.
> 3. **REDSEA v10** — `comp_mode 0→1` (reinforcement restored) and `norm_form recipient→donor` (the
>    published mass-conserving operator); nuclear/ECM channels no longer corrected. So the
>    "subtract-only, 1px band, α=1" description above is **out of date**.
> 4. **RESTORE v10** — the negative control was gated on the target's own Otsu floor, so the divisor
>    *was* that floor (median ratio 1.044 over 49 donor-pairs). Now
>    `control_definition` + `divisor_statistic` are configurable; `reference_only` + `p99` is the
>    leading candidate (E-cadherin 43–54%, still short of the ~75% the accounting implies).
>
> The E-cadherin under-call is **narrowed, not closed**. ~~Next step is the anatomy check~~ — **done,
> 2026-07-27**; see the round report above. Next step is Gate-1→3 validation of Pan-Cytokeratin,
> Ker8/18 and EpCAM as production RESTORE targets, then the divisor decision, then Gate 4→5.

**Status (updated 2026-07-24) — manuscript-faithful rebuild, at Gate 2 of 7.** `RESTORE.pdf` +
`RESTOREsup.pdf` (repo root) are **ground truth and override prior code, handoffs, and docs**. The
gate plan is [`docs/restore_faithful_rebuild_plan.md`](docs/restore_faithful_rebuild_plan.md).
Submodule branch `redsea-restore-curated-5compartment`, METHOD_VERSION `v9`.

- **Divisor = MAXIMUM target intensity in the accepted target-negative controls** (the manuscript
  equations: controls ≤ 1, signal > 1). Released-code `mean+3σ` is **sensitivity analysis only**.
  Each marker's background comes from **its own** negative control — never shared, averaged, or
  medianed across markers or references.
- **Separator = 2-component NNMF** on raw QuPath compartment means (membrane + whole-cell), trained
  on deterministic equal-size arms with robust min-max scaling; GMM kept as comparator. Ordered
  controls: raw double-low QC → NNMF assignment → Step 1 anchor/separator → target-component
  support. Marker-specific raw arm floors, no silent fallback (CD3e/CD20 use max(Otsu, Triangle)).
- **The positive CALL is divisor-gated only: `target_redsea > Step 2 divisor`** (METHOD_VERSION v7,
  2026-07-23). The former `| target_component_supported` rescue — which additionally called every
  NNMF-target-component cell left of Step 1 even below the divisor — is removed from the call. It was
  negligible for clean markers (≤0.3 pp) but **doubled a high-background marker like CD11b** by
  retaining background cells the divisor correctly excludes (6521 CD11b ← PanCK: 14% divisor-gated vs
  28% by component-support, up to 95% via CD3e). The below-divisor component cells stay **tracked**
  (`target_supported_below_divisor`) and are shown in the review bundle as a distinct
  "retained, not called" category (so a marker they help stays visible for a per-marker decision) —
  they are simply not counted positive. Investigation: `CD11b_from_PanCK` is a genuine overcall in the
  *call*, not display (unlike `CD3e ← CD163`); RESTORE's reciprocal divisor (max CD11b among PanCK⁺
  controls) sets the right threshold — the retention was the leak.
- **CD11b ← CD3e RESOLVED + one reference FROZEN per target (Gate-2 reference selection, 2026-07-24,
  METHOD_VERSION v8).** The screen never combines references (each pair emits its own divisor; extra
  references per immune target are comparators), so the `CD11b ← CD3e` overcall was an **isolation
  artifact** — CD3e is the *minimum/outlier* divisor of CD11b's references on 19/20 donors; any single
  epithelial reference (or median/max across them) gives the validated ~10%. No resampling can fix it:
  CD3e's control arm is ≤1,706 cells on every donor (below the 2,000 equal-N cap) so its equal-N max
  equals its full max, and nothing can sample a CD11b value above CD3e's brightest control cell (9,923)
  — the false-positive epithelial cells sit at 9,923–13,026, a background only an epithelial reference
  contains. **Maintainer rule:** one predeclared reference per target; selection condition = the
  reference-positive negative control must be ≥ the target-positive population — a **pure proportion**
  `reference_n/target_n ≥ 1`, **no absolute floor** (images vary widely in total/per-marker counts).
  `restore_validation.py` adds `reference_to_target_ratio` + `reference_undersized` as **reported-only**
  diagnostics (never gate `PairState`, never touch a divisor); the v10 screen is byte-identical to v9
  on every divisor. New module `phenocycler/reference_selection.py` builds the dossier
  (`data/restore_pair_validation/reference_selection_v1/`: `reference_matrix.csv`,
  `frozen_reference_matrix.csv`, 7 per-target figures) — **no divisor/normalized/phenotype output.**
  **Frozen references:** CD31←E_cadherin, Vimentin←E_cadherin, B3TUBB←EpCAM, CD3e←E_cadherin,
  CD68←E_cadherin, **CD11b←EpCAM** (maintainer override of the marginal E_cadherin baseline),
  **E_cadherin←CD31** (accepted despite failing the proportion — advisory for a structural target;
  endothelium is exclusive but sparse). **Rejected by the rule:** `CD11b ← CD3e` (0.04×),
  `CD11b ← Pan_Cytokeratin` (0.78×), `CD68 ← CD3e` (0.20×); the divisor-consistency cross-check agrees
  (the undersized CD3e reference is also the low-divisor outlier). Production wiring
  (`BASELINE_PAIRS`/`config.py`/`lineage.py`) is **unchanged** — that is Gate 4.
- **Gate 3 full-cell stress pilot RUN and PASSED; two high-background gate defects fixed (METHOD v8→v9,
  2026-07-24).** Isolated pilot `data/restore_pair_validation/gate3_pilot_v2/` (noncanonical): the 7
  `ACCEPTED_PAIRS` × 7 donors (6539, 6476 + the weakest-case donor per pair), marked ACCEPTED via
  `pair_reviews_accepted.csv`. Confirmed the divisor policy on full cells — every full-cell divisor is
  **higher** than the 5% screen's (max over the *complete* control population; median ~1.15×) so calls
  are correspondingly lower; **CD11b←EpCAM lands ~7–10%** (validated). The first run surfaced two
  abstentions on **high-background** markers (maintainer: 6476 has high CD3e *background*, 6591 high
  CD11b background — not extreme true expression) where the provisional QC gates misfired though the
  divisor was well-defined. Fix, since RESTORE depends only on the divisor: **seed-stability now gates
  on divisor reproducibility** (`min/max seed divisor ≥ min_divisor_reproducibility` 0.90), not
  control-set membership (`control_jaccard` kept as a reported diagnostic — verified on 6476:
  `divisor_reproducibility=1.0` while `control_jaccard=0.8846`); and a **low target-arm mean SBR is now
  a reported flag `target_low_sbr`, not an abstention** (manuscript: flag weak calls, don't zero them;
  genuine/confirmed absence still abstains). Re-run under v9: **all 49 donor-pairs VALID, divisors
  byte-identical**, invariants hold. 145 tests. Submodule `8c77d05`. The 20-donor 5% screen was **not**
  regenerated (its divisors + frequency ratios are v9-identical; only diagnostic state columns shift).
- **CD20 takes NO part in RESTORE — neither target nor reference — and is deferred to a second pass
  after the broad cell types are settled** (maintainer decision 2026-07-23). As a reference its
  positive cells would BE the negative control, and its exclusive arm ran 41–2,026 cells/donor vs
  2,206–50,000 for the structural markers, producing inverted fits (`6591 CD11b ← CD20`: 50 reference
  cells thresholding 30,507 target cells). As a target its Step-1 anchor bottomed out at 14 cells with
  anchor recovery to 0.217, and reviewer inspection of the CD20 review PDFs found overcalling
  throughout, with both qptiff channels stretched from near-background into block-like tile artifacts
  — consistent with CD20 being the panel's second-weakest channel (display ceiling median ~503 counts
  vs ~6,800 for E-cadherin; only CD163 at ~181 is lower). Enforced at import via
  `RESTORE_EXCLUDED_MARKERS` in `restore_validation.py`, which rejects an excluded marker in **either**
  role, plus a self-pair guard. CD20 is still **extracted** (`DEFERRED_EXTRACTION_MARKERS` feeds
  `QUPATH_MARKERS`) and keeps its `max(Otsu, Triangle)` input-floor policy, so the deferred pass needs
  no QuPath re-export. CD3e is the immune reference in its place. Consequences: CD20 is not in
  `COMPARTMENT_GATES["Immune"]` (B cells reach Immune via CD79a); `lineage.py` keeps its inert
  `pos("CD20")` term, and since CD79a is also unscreened the B/Plasma branch yields nothing until the
  second pass — those cells stay `Immune-other`, an abstention, not a negative call.
- **Nothing is accepted yet.** No pair is approved, no divisor emitted, no normalized/phenotype/app
  output written — `data/` deliberately holds only `cells/`, `cells_redsea/`, `redsea_scratch/`,
  `_redsea_validate/`, `singlecell_protein.h5ad`, and the review bundles. Of the 400 screen rows
  (20 donors × 20 pairs) 384 sit at `PAIR_REVIEW_PENDING`, 10 `MODEL_UNSTABLE`,
  6 `LOW_SIGNAL_UNRESOLVED`.
- **NEXT: Gate 4 — production wiring.** Gate 2 (reference selection, frozen) and Gate 3 (full-cell
  stress pilot, passed under v9) are complete; the maintainer signed off the pairs ("all pairs look
  good enough to continue"). Gate 4: update `config.py` to the approved markers/rules; ordered residual
  **Immune → Epithelial → Endothelial → Neural → Mesenchymal → Other**; remove the silent
  missing-marker fallback in `lineage.py` (`_pos_norm_getters`); implement **tri-state
  `positive|negative|unavailable`** calls + per-donor-marker QC disposition + provenance; `Other` only
  when every required gate is valid-and-negative, blocked cells `Unresolved`. Then Gates 5–7. The
  frozen references are in `reference_selection_v1/frozen_reference_matrix.csv`; the accepted-pair set
  is `restore_validation.ACCEPTED_PAIRS`. Read `data/restore_pair_validation/README.txt` for active vs
  superseded artifacts.
- **CD3e ← CD163 investigated (2026-07-23) — display-only; the calculation is sound; NOT the CD20
  failure mode.** The blocky `CD3e_from_CD163` qptiff panel is a *display* stretch: the bundle windows
  every channel to `[1st, 99.8th pct] + γ0.65`, so a channel whose crop median is ~7 counts (CD163 is
  the panel's weakest, display ceiling 60–425) has its noise floor amplified into a mottled/block-like
  look; a linear render removes it and there is no acquisition damage (no dead tiles in the per-cell
  field). The window is computed in `_display_marker_channel` for the RGB panels + manifest metadata
  only — it never re-enters a computed metric. The calculation is unaffected: across all 20 donors
  CD3e's divisor via CD163 matches its divisor via the five strong references within ±20% (median
  1.08×), and the NNMF control arm is a real CD163⁺/CD3e-low macrophage population of 1,336–10,201
  cells — nothing like CD20's 14-cell anchor. CD163 is therefore **weak-but-valid** (`reference_fold`
  4.4× vs 7.6× CD68 and 19–25× the keratins, above the 2× SBR floor) with a mild +8% divisor bias.
  This **corrects** an earlier "same failure mode as CD20" flag, which the numbers disprove. Note CD163
  is reference-only (never a screen target), so it gets no threshold of its own regardless — whether
  CD3e should lean on the weakest reference when CD68/keratins exist is a Gate-2 *preference*, not a
  correctness fix, and nothing needs regenerating.
- **Superseded**, parked on submodule branch `superseded/multiref-median-8gate`: the multi-reference
  median-threshold RESTORE and the partial 8-gate lineage edit. Do not revive before Gate 4.
- The app still serves the **provisional** un-QC'd labels until the Gate-7 rebuild.

**Reproducible workflow:** the submodule's numbered notebooks **`Phenocycler_Analysis/notebooks/00–06`** run the pipeline end-to-end (thin orchestrators over `python -m phenocycler.pipeline --only <stage>`), ending with a **cell-type × marker dotplot + heatmap** (`03b_restore_diagnostic.ipynb` traces per-image RESTORE thresholds; `06_figures.ipynb` builds the composition/dotplot/heatmap). The old IES orchestrator `notebooks/senior_phenotyping_redsea_restore.ipynb` and the `scripts/senior/*.py` step scripts were retired (scripts archived to **`archive/phenotyping_senior_preport/`**; `run_pipeline.sh` is now a thin wrapper over the submodule). App tab-by-tab capability is tracked in `app/shiny_app/docs/capability_status.md`.

## Spatial / Statistics / Trajectory / AI tabs

Per-tab implementation details — Spatial Neighborhood Analysis (Phase 7), Single-Cell Drill-Down (Phase 8), Statistics Tab pseudoreplication fix (Phase 16), Cell-Count-Weighted Trajectory (Phase 12), Pseudotime Heatmap, AI Assistant (Phase 10) — live in `app/shiny_app/CLAUDE.md`.

## Important Conventions

- `ggplot2::coord_sf()` and `ggplot2::geom_sf()` - these are ggplot2 functions, NOT sf functions
- `PIXEL_SIZE_UM = 0.3774` - micrometers per pixel conversion constant
- UMAP "Selected Feature" uses continuous viridis (inferno) colormap scaled to data min/max
- Donor status colors: ND = steel blue (#4477AA), Aab+ = burnt umber (#CC6633), T1D = forest green (#228833) — Paul Tol bright variant, colorblind-safe, centralized as `DONOR_COLORS` in `00_globals.R`
- Phenotyping uses Rules1 (`data/phenotype_rules.csv`) - 21 phenotypes, 18 markers. CSV headers use adata.var_names (CD8a, HLADR, PDPN — not CD8, HLA-DR, Podoplanin)
- **Donor metadata**: Comes from `islets_core_fixed.h5ad` obs, NOT from `CODEX_Pancreas_Donors.xlsx` (different cohort)
- **User-facing terminology**: Use "Sex" not "Gender" in all UI labels, plot titles, and methods text. The underlying data column remains `gender` for backwards compatibility.
- **H5AD obs index**: `islets_core_fixed.h5ad` index name is `islet_id` -- same as column, use `reset_index(drop=True)`
- **Case ID zero-padding**: GeoJSON files use `0112.geojson` (4-digit padded), data uses `112` (unpadded). `load_case_geojson()` and click handlers use `sprintf("%04d", ...)` fallback
- **Log-scale with zeros**: Use `scales::pseudo_log_trans(base=10)` instead of `scale_y_log10()` -- zeros map to 0 (visible) instead of -Infinity (dropped)
- **Plot defaults**: Point size = 3.0, transparency = 0.6
- **Plot fonts & downloads**: Every plot uses `theme_islet()` (ggplot) and `plotly_islet_fonts()` + `plotly_islet_config()` (plotly) — single source of truth in `app/shiny_app/R/00_globals.R`. Base 16pt / titles 18pt / axes 14-16pt / legends 14-16pt / heatmap rows 12pt min. Plotly modebar PNG downloads ship 2800 × 1800 px (scale=2). ggplot `renderPlot` outputs that should be downloadable are refactored into a `reactive()` consumed by both `renderPlot` and a `downloadHandler` using `ggsave(..., width=10, height=7, dpi=300)` (heatmaps 12 × 8). See `app/shiny_app/CLAUDE.md` § "Plot Fonts & Publication-Quality Downloads".
- **Phenotype marker rules in UI**: `data/phenotype_rules.csv` is parsed once at startup into `PHENOTYPE_RULES_CSV` in `00_globals.R`, with hierarchical composition (sub-phenotypes inherit parent constraints) and an alias map for the four PHENOTYPE_COLORS-vs-CSV naming mismatches (`CD8a Tcell` ↔ `CD8 T cell`, etc.). Surfaced via (1) inline blue hint under every phenotype `selectInput` showing the active rule, (2) global "Phenotype rules — marker definitions" button at the top of every tab opening a reference modal listing all 21 phenotypes + markers, (3) expandable rules list in the drill-down composition panel naming each phenotype present in the clicked islet. See `app/shiny_app/CLAUDE.md` § "Phenotype Rules — Marker Definitions in the UI".
- **Peri-islet data guard**: Always check `total_cells_peri > 0` before using peri metrics
- **Column name sanitization**: Phenotype names use `_` for spaces, `plus` for `+` in peri-islet columns
- **Single-cell Parent column**: `Islet_N` = core cells, `Islet_N_exp20um` = peri-islet cells
- **macOS plotly clicks & Viewer iframe (Jun 2026)**: `plotly_click` is suppressed on macOS — drill-down uses a native-click bridge (`islet_click_bridge()` in `00_globals.R`). The Viewer iframe is mounted once and its `src` updated in place via `session$sendCustomMessage` to avoid recreating the WebGL canvas (macOS WebRender flicker). Details + the `reactiveVal`-dedup gotcha in `app/shiny_app/CLAUDE.md`.

App-specific UI conventions (CSS overflow, font minimums, large-scatter rendering, palette syncing, reticulate bridge, synthetic-union rows, macOS click bridge & Viewer iframe persistence) live in `app/shiny_app/CLAUDE.md`.

## Deployment Architecture

- **Production URL**: `http://<server-ip>:8080/islet-explorer/`
- nginx (port 8080) reverse-proxies `/islet-explorer/` -> shiny-server (port 3838)
- shiny-server serves from symlink: `/srv/shiny-server/islet-explorer` -> `app/shiny_app/`
- Changes take effect when shiny-server spawns a fresh R worker (no restart needed; kill stale workers if needed)
- Dev server (`Rscript -e 'shiny::runApp(".", port=7777)'`) is for local testing only

## Running the App

Production (always use this):
```
http://<server-ip>:8080/islet-explorer/
```

Development only:
```bash
cd app/shiny_app
Rscript -e 'shiny::runApp(".", port = 7777)'
```

## Running Pipeline Scripts

See `scripts/CLAUDE.md` for the full per-script command list. Quick rebuild:

```bash
conda activate scvi-env
python scripts/reaggregate_islets.py          # ~15 min: aggregation + trajectory + Leiden
python scripts/compute_neighborhood_metrics.py # ~3 min
python scripts/build_h5ad_for_app.py          # ~1 min: merge into app H5AD
```

## Analysis Pipeline

See `scripts/CLAUDE.md` for the Python pipeline; phenotyping steps 1–5 (REDSEA → RESTORE → 7-class broad lineage) live in the `Phenocycler_Analysis` submodule. (The old `islet_analysis/` dir was removed; its content folded into `scripts/senior/` + the submodule.)
See `data/DATA_PROVENANCE.md` for full data lineage documentation.
See `docs/user_guide.md` for end-user documentation.

## Cross-Session Knowledge

Cross-session development knowledge for this project lives in Claude's auto-memory (see `MEMORY.md`). The earlier `Skills_Registry/` submodule + `/advise` + `/retrospective` workflow was retired Apr 2026 in favour of the auto-memory mechanism.

## R Dependencies

shiny, shinyjs, plotly, ggplot2, dplyr, tidyr, readxl, sf, jsonlite, RColorBrewer, scales, anndata, reticulate, httr2, stringr, lmerTest, lme4, emmeans
