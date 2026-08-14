# Xenium Panel Explorer

Interactive Shiny app for exploring 49 curated subpanels of the
**10x Genomics Xenium Prime 5K Human Pan Tissue and Pathways** panel
(plus a 100-gene custom T1D-GWAS add-on), and for re-clustering
user-loaded Xenium runs against any chosen subpanel with full control
of the clustering / sub-clustering pipeline and `clustree` for
resolution-stability inspection.

The app is built for a single-user research workflow on a workstation or
on HiPerGator. There is no authentication or multi-user state. See
[`CLAUDE.md`](CLAUDE.md) for the full build spec.

---

## Table of contents

- [Quick start](#quick-start)
- [Installation](#installation)
- [Launching the app](#launching-the-app)
- [Using the app](#using-the-app)
  - [Overview](#overview)
  - [Panel Browser](#panel-browser)
  - [Load Xenium](#load-xenium)
  - [Panel-vs-Data](#panel-vs-data)
  - [Cluster](#cluster)
  - [Subcluster](#subcluster)
  - [Markers](#markers)
  - [Clustree](#clustree)
  - [Export](#export)
- [Typical workflows](#typical-workflows)
- [Repo layout](#repo-layout)
- [Tests](#tests)
- [Running on HiPerGator](#running-on-hipergator)
- [Known environment notes](#known-environment-notes)
- [Out of scope](#out-of-scope)

---

## Quick start

```bash
# 1. system deps (Ubuntu 22.04 / 24.04)
sudo apt-get install -y build-essential cmake gfortran pkg-config \
  libcurl4-openssl-dev libssl-dev libxml2-dev \
  libhdf5-dev hdf5-tools \
  libfontconfig1-dev libfreetype6-dev libpng-dev libtiff5-dev libjpeg-dev \
  libcairo2-dev libxt-dev libharfbuzz-dev libfribidi-dev \
  libgit2-dev libuv1-dev

# 2. R packages (R ≥ 4.4)
R -e 'install.packages("renv"); renv::restore()'

# 3. launch
R -e 'shiny::runApp("app.R", launch.browser = TRUE)'
```

The app comes up at `http://127.0.0.1:<random_port>`. With no Xenium data
loaded the **Overview** and **Panel Browser** tabs are fully usable
against the audit CSVs in `data/panel_audit/`; the analytical tabs gate
behind a successful **Load Xenium**.

---

## Installation

R 4.4 or newer is required. The R *package* layer is reproducible
through `renv`; the R *interpreter* itself is host-managed. If your
workstation is set up by IT or carries multiple R series, the simplest
way to land on a known-good R is `rig`.

### Option A — `rig` (recommended for multi-R-version workstations)

[`rig`](https://github.com/r-lib/rig) is a small R-installation
manager (the R-world analogue of `nvm` / `pyenv`). It puts every R
release side-by-side under `~/.local/share/R/` and lets you pick a
default per-project.

```bash
# install rig (Ubuntu / Debian)
curl -L https://rig.r-pkg.org/deb/rig.gpg \
  | sudo tee /etc/apt/trusted.gpg.d/rig.gpg > /dev/null
sudo curl -L -o /etc/apt/sources.list.d/rig.list https://rig.r-pkg.org/deb/rig.list
sudo apt-get update
sudo apt-get install r-rig

# install + select the R version this project pins
rig add 4.4.1
rig default 4.4.1
```

`rig` doesn't replace any existing system R; the executable just
becomes the first one on `$PATH` for new shells. Verify with
`R --version` from a fresh terminal — you should see 4.4.1.

### Option B — system R, however your host gives it to you

If you already have R ≥ 4.4 from `apt`, conda, brew, or a HiPerGator
module, skip `rig` entirely.

### R packages

From the project root, with R ≥ 4.4 on `$PATH`:

```r
install.packages("renv")
renv::restore()                   # pulls every package pinned in renv.lock
```

If `libuv1-dev` isn't on the host, `fs` will refuse to build. Either
install the dev package via apt, or set `USE_BUNDLED_LIBUV=1` before
`renv::restore` so `fs` builds with its bundled libuv.

`presto` is installed from GitHub (`immunogenomics/presto`); `renv`
records that source in the lockfile so `restore()` reproduces it.

The first run takes ~10–15 minutes for source builds (Seurat, harmony,
arrow, rhdf5, clustree and their transitive deps). Subsequent restores
are seconds when the renv cache is warm.

### 10x pre-designed panel data

The Panel Browser compares your tissue subpanels against 10x Genomics'
public pre-designed Xenium panels. Fetch them once with:

```bash
Rscript scripts/fetch_10x_panels.R
```

This writes 11 CSVs (Human Lung / Brain / Breast / Colon / Skin /
Multi-Tissue / Immuno-Oncology / Prime 5K, plus Mouse Brain / Tissue
Atlas / Prime 5K) into `data/reference_panels/` and a `manifest.csv`
indexing them. The fetched files are checked into the repo so a fresh
clone has the panels available; re-run the script to pick up 10x's
revisions.

### Startup environment check

The app runs `check_environment()` at session start and prints any
issues to the terminal **and** to a banner card on the Overview tab:

- R version below the pinned minimum (4.4),
- any `Imports:` package that fails to load via `requireNamespace`,
- a missing `data/panel_audit/subpanel_summary_v2.csv`,
- (informational) `chromote` couldn't find a Chrome binary — only
  matters for the `shinytest2` end-to-end test.

These are warnings, not errors: the app keeps starting so a
near-miss configuration is obvious rather than blocking.

---

## Launching the app

From the project root:

```bash
R -e 'shiny::runApp("app.R", launch.browser = TRUE)'
```

Or in an interactive R session:

```r
pkgload::load_all(".")
shiny::runApp(xenium_panel_app(), launch.browser = TRUE)
```

`app.R` is a two-line entry point — `pkgload::load_all(".")` followed by
`xenium_panel_app()`, which returns a `shiny.appobj`.

To pin the port (useful behind an SSH tunnel):

```r
shiny::runApp(xenium_panel_app(), port = 4321, host = "127.0.0.1",
              launch.browser = FALSE)
```

---

## Using the app

The navbar has nine tabs, in roughly the order you'd touch them. The
app maintains a small piece of cross-tab state (`app_state`) so loading a
dataset or running the cluster pipeline once propagates to every
downstream tab without re-doing work.

### Overview

The landing page summarising the audit. No inputs to fill; everything
comes from `data/panel_audit/subpanel_summary_v2.csv`.

**What you see**
- Sidebar — counts (subpanels indexed, CSV files loaded, 5K-shared
  genes, custom-100 genes).
- "Genes per subpanel" — horizontal `plotly` bar of `n_genes` for every
  subpanel; hover for the description.
- "Subpanels" `DT::datatable` — the full summary, sortable / filterable
  / searchable.

**Interactions**
- **Click a row** in the summary table to jump to **Panel Browser**
  pre-filtered to that subpanel. The cross-tab handoff resolves the
  long summary names (e.g. `99d_truly_unannotated_subset_of_99c`) to
  the short file stem (`99d_truly_unannotated`) automatically.

When to use it: orientation. If you're new to the audit, sort by
`n_genes` to see which biology is over- or under-represented.

### Panel Browser

The deep-dive browser for any single subpanel.

**Sidebar**
- **Primary panel** — selectize across the active tissue's subpanels,
  the three uncategorized residuals (`99`, `99c`, `99d`), the tissue's
  custom panel (e.g. `T1D-GWAS custom 100`), `xenium5k_in_audit` (the
  joined audit reference), and 10x Genomics' eleven pre-designed Xenium
  panels (Lung, Brain, Breast, Colon, Skin, Multi-Tissue & Cancer,
  Immuno-Oncology, Mouse Brain, Mouse Tissue Atlas, and the Human /
  Mouse Prime 5K Pan-Tissue panels). The 10x panels are grouped under
  a "10x pre-designed panels" optgroup. Their gene lists are fetched
  by `scripts/fetch_10x_panels.R` (see Installation below).
- **Compare with (optional)** — pick a second panel to overlay.
- **Hide genes with `exclude_recommended = yes`** — subpanels carry the
  audit's `exclude_recommended` flag; checking this filters them out
  from both the table and the scatter.

**Main panel**
- **Gene `DT`** — every column from the subpanel CSV, filterable per
  column. Hidden by default until you select a panel.
- **Detection scatter** — `detection_pct_0041323` (x) vs
  `detection_pct_0041326` (y) for every gene in the primary panel.
  Colour is `log2_detection_ratio_326_over_323` on a centred RdBu
  scale (red = higher in 326, blue = higher in 323). The y = x diagonal
  is dashed grey. Hover shows the gene symbol.
- **Overlap with comparison panel** — only populated when a comparison
  panel is selected. Shows `(A only / both / B only / union)` counts
  plus the first 60 shared gene symbols.

When to use it:
- Drill into one biology (e.g. `08_antigen_presentation_MHC` to inspect
  the audit's MHC machinery genes).
- Sanity-check overlap between two related subpanels (e.g. immune
  T-cell vs antigen-presentation MHC).
- Identify high-detection genes you'd want to feature in marker plots
  later.

### Load Xenium

Brings a Xenium dataset into the app. **The analytical tabs (Compare,
Cluster, Subcluster, Markers, Clustree, Export) require a load.**

**Three ingest paths** (CLAUDE.md §2)
1. **Xenium output bundle directory** — must contain
   `cell_feature_matrix.h5` (10x's CSC HDF5) and `cells.parquet`
   (or `cells.csv.gz`). `transcripts.parquet` is optional and is
   recorded in `xen@misc$transcripts_parquet` for downstream use; the
   transcripts themselves aren't loaded into the Seurat object.
2. **Saved Seurat object** as `.rds` or `.qs2`.
3. **`.h5ad`** — wired but not implemented. Pick one of the first two
   for now.

**Sidebar**
- **Pick Xenium bundle directory** — `shinyFiles::shinyDirButton`,
  scoped to your home / cwd / system root.
- **Pick saved Seurat (.rds / .qs2)** — file picker.
- **Or paste a path** — text input + load button. Useful when the
  picker can't reach the volume your data lives on.
- **Bypass cache (force re-ingest)** — checkbox. Cache hits are keyed
  by `rlang::hash(abs_path + mtime)`; tick this if you've replaced a
  file in place without changing its mtime.

**Main panel — after a successful load**
- Status card — path / class / cells / genes / cache state.
- Panel-validation alert — green ≥95%, yellow ≥50%, blue otherwise.
  Reports: data gene count vs. reference (`xenium5k_in_audit ∪
  custom_T1D_GWAS_panel`), intersection size and pct of reference
  covered, missing-from-data, and extras (genes in the data not in the
  reference).
- Spatial preview — rasterized `scattergl` of every cell's
  `x_centroid` / `y_centroid`.

When to use it: load once at session start. Subsequent reloads of the
same path are `qs_read` cache hits (fast).

### Panel-vs-Data

Cross-tabulates every subpanel against the loaded data. Tells you which
biology is well-detected in *this* sample.

**Sidebar**
- **Min detection % (loaded data)** — slider 0-100. Genes below this
  threshold don't count toward `n_passing` / `frac_passing`.
- **Top-N genes per subpanel (heatmap)** — slider 3-30.
- **Hide xenium5k_in_audit row** — checkbox. The full 4K-gene reference
  row dwarfs everything else; hide it for clarity.

**Main panel**
- **Coverage `DT`** — one row per subpanel, plus rows for `custom-100`
  and `xenium5k_in_audit`. Columns:
  `n_genes`, `n_present`, `n_passing`, `frac_present`, `frac_passing`,
  `mean_detection_pct`, `median_detection_pct`, `mean_expr`,
  `median_expr`. Sortable, filterable.
- **Fraction-passing bar** — horizontal `plotly` bar across subpanels.
- **Top-N heatmap** — `geom_tile` of the top-N genes per subpanel
  ranked by detection %, gene labels overlaid, viridis-C fill.
  White cells = no gene at that rank.

The detection-cutoff and top-N values are stored in
`app_state$compare_min_det` / `compare_topn` so subsequent tabs can
read them.

When to use it: triage. After a load, glance at the bar to see which
subpanels light up; pick those for the Cluster pipeline.

### Cluster

The core analytical tab. Wires every control 1:1 to
`run_cluster_pipeline()` (CLAUDE.md §6).

**Sidebar groups**
- **Feature set**
  - Subpanels selectize (multi) — union of selected subpanels' genes
    becomes the variable feature set.
  - Add individual genes (multi) — typed search over the loaded data's
    rownames.
  - "Use all genes" — overrides the subpanel choice and uses every
    gene in the assay. Sensible default for first runs.
- **Cell filters**
  - `nCount` and `nFeature` range sliders, auto-bounded by the data.
  - Optional categorical metadata filter (any character / factor
    column).
- **Normalization** — `LogNormalize` (default), `SCTransform` (auto-
  warn + fallback to LogNormalize above 200 k cells), or `skip`. The
  scaling checkbox carries a tooltip noting that disabling scaling
  sometimes helps when low-count genes produce UMAP "spikes".
- **PCA** — `npcs` slider 5-50. Internally clamped to
  `min(npcs, n_features - 1, n_cells - 1)`.
- **Batch correction** — `None` or `Harmony`. Harmony reveals
  group-by + theta inputs; gracefully skips with a warning if the
  group-by var has < 2 levels or `harmony` isn't installed.
- **Neighbors / UMAP** — `k.param`, `n.neighbors`, `min.dist`,
  `metric` (cosine / euclidean).
- **Clustering** — Louvain (default) or Leiden. Leiden falls back to
  Louvain with a warning when `leidenAlg` isn't installed.
- **Resolution range + step** — defines the resolutions in the sweep.
  Each resolution gets its own `seurat_clusters_res_<r>` column.
- **Random seed** — wired through `set.seed`, `RunPCA(seed.use)`,
  `RunUMAP(seed.use)`, `FindClusters(random.seed)`.
- **Run pipeline** — triggers the run behind a `waiter` overlay.

**Main panel — after a successful run**
- Status alert summarising cells out × features × resolutions.
- Active-resolution slider (auto-snapped to the median resolution).
- **UMAP** — `plotly` scattergl coloured by the active resolution;
  hover shows cell ID + cluster.
- **Cluster size bar** — `plotly` bar of cells per cluster.
- **Spatial scatter** — `ggrastr::geom_point_rast` past 50 k cells,
  otherwise plain `geom_point`.
- **Cluster × sample stacked bar** — only renders when `orig.ident`
  has ≥ 2 distinct values; otherwise shows a centred message.
- **Run log accordion** — every parameter used, including the
  resolution list and the chosen algorithm.

When to use it: the heart of the app. Pipeline runs are typically a few
seconds on subsetted data, tens of seconds on a full sample.

### Subcluster

Drill into a chosen parent cluster, optionally on a different feature
set, without losing the parent context. State is a stack — push by
drilling in, pop with **Back**.

**Sidebar**
- **Parent resolution** — picks which `seurat_clusters_res_<r>` column
  drives the cluster picker.
- **Cluster(s) to keep** — multi-select; one cluster ID or many.
- **Feature set (optional)** — subpanel selectize + "use all genes"
  toggle. If you don't tick "use all", you must pick at least one
  subpanel.
- **Slim pipeline opts** — `npcs`, resolution range + step, seed.
- **Drill in** / **Back** buttons.

**Main panel**
- **Stack** breadcrumb — ordered list of stack levels with cell counts.
  Active level is highlighted.
- **Status alert** — current depth + level summary.
- **UMAP** + **Spatial scatter** of the active level, driven by an
  active-resolution slider.

The Cluster tab's `xen_clustered` is the immutable root — it never
mutates. When you re-run the Cluster tab the Subcluster stack resets
to that new root.

When to use it: when a cluster looks heterogeneous and you want to know
*why* without rerunning the whole pipeline.

### Markers

`presto::wilcoxauc` marker tables, plus dotplot, mean-expression
heatmap, and a per-gene FeaturePlot panel.

**Sidebar**
- **Source** — Cluster root or Subcluster stack top.
- **Grouping** — picks the `seurat_clusters_res_<r>` column.
- **Filters** — `min_pct_in` (slider 0-100), `max_padj` (slider
  0-0.5), top-N per cluster (slider 3-30).
- **FeaturePlot gene** — selectize over the assay's rownames.
- **Compute markers** — runs presto on the assay's data layer.
  Cached per `(run_id, group_col)`; UI tweaks of the filters don't
  recompute.
- **Download markers CSV** — pulls the unfiltered marker frame.

**Main panel**
- **Top-N table** — filtered `DT` of the top-N genes per cluster.
- **Ranked AUC bar** — horizontal `plotly` bar coloured per cluster.
- **DotPlot** — `Seurat::DotPlot` over the top-N genes, x-axis
  rotated 45°.
- **Heatmap** — `geom_tile` of the mean expression of the top-N
  markers per cluster, viridis-C fill.
- **FeaturePlot panel** — pick any gene in the sidebar to colour the
  UMAP scattergl and the spatial scatter by its expression. Both
  panels rasterize past 50 k cells.

When to use it: after clustering, to characterise clusters and pick
genes for downstream interpretation.

### Clustree

Read-only: consumes the resolution sweep produced by the Cluster tab.
Doesn't trigger reclustering — but does expose a "Use this resolution"
button that hands a chosen value back to the Cluster tab's slider.

**Sidebar**
- **Resolutions** — checkbox group; default = all available.
- **Node colour** — radio: cluster ID (default), `sc3_stability`, or a
  gene from any subpanel. The gene field only shows in gene mode and
  is populated reactively from `rownames(xen_clustered)`. If the gene
  isn't in the assay's data slot, the plot renders an error
  annotation rather than crashing.
- **Edges** — colour by `count` or `in_prop`; `prop_filter` slider
  (default 0.1) to drop weak edges.
- **Node size** — cell count (default) or `sc3_stability`.
- **Use this resolution** — sets `app_state$cluster_jump_res` and
  navigates to the Cluster tab, which snaps its active-resolution
  slider to the chosen value.
- **Downloads** — tree PDF, tree PNG, edge-table CSV.

**Main panel**
- Status card — counts of available resolutions.
- **Resolution-tree** plot.
- **Tree on UMAP** (`clustree::clustree_overlay`) — same tree
  superimposed on UMAP coordinates.
- **Stability summary `DT`** — for each resolution: `n_clusters`,
  `ARI_vs_prev`, `NMI_vs_prev` (NA when `aricode` isn't installed),
  `frac_changed_vs_prev`.

When to use it: deciding which resolution to commit to. Stable
sub-trees that don't reorganise across resolutions are good
candidates.

### Export

Download the session.

**Sidebar**
- **Source** — Cluster root or Subcluster stack top.
- **Seurat object** — `.qs2` of the chosen source.
- **Cluster assignments** — CSV: one row per cell with
  `cell_id`, centroids, `orig.ident`, `nCount_*`, `nFeature_*`, and
  every `seurat_clusters_res_*` column.
- **Marker tables** — combined CSV of every cached marker run, with
  `source_run_id` and `group_col` columns parsed out of the cache key.
- **Session report** — self-contained HTML with audit overview, loaded
  dataset + validation, cluster pipeline run-log, the Subcluster stack,
  the most-recent marker table, and `sessionInfo()`. **No pandoc /
  quarto / commonmark dependency** — the report is written via
  `htmltools::save_html` directly.

When to use it: end of session, or whenever you want to checkpoint
state into `cache/` for hand-off to a separate analysis.

---

## Typical workflows

### A. Compare a fresh dataset against the audit (~5 min)
1. **Load Xenium** → pick the bundle directory.
2. Read the panel-validation alert. ≥95% green is the happy path.
3. **Panel-vs-Data** → set Min detection to 5–10%, glance at
   `frac_passing` to see which subpanels light up.
4. **Export** → session report HTML for the lab notebook.

### B. Recluster against an immune-only feature set (~2 min after load)
1. **Load Xenium** as in A.
2. **Cluster** → in the sidebar:
   - Subpanels: `03_immune_T_cell`, `04_immune_B_plasma`,
     `06_immune_myeloid`, `09_cytokines_chemokines`,
     `29_allograft_autoimmunity_antigen_processing`.
   - Untick "Use all genes".
   - Resolution range 0.2 → 1.0, step 0.1.
   - Click **Run pipeline**.
3. Slide through resolutions on the UMAP card to find sensible cluster
   structure.
4. **Clustree** → confirm resolution stability; "Use this resolution"
   to snap the Cluster tab.

### C. Drill into a cluster you don't recognise (~30 sec)
1. With a Cluster run already in `xen_clustered`, go to **Subcluster**.
2. Pick the parent resolution and the cluster id(s).
3. Optionally switch the feature set (e.g. T-cell only when drilling
   into a putative T-cell cluster).
4. **Drill in**. Inspect the new UMAP + spatial; **Back** to undo.

### D. Find markers, export, write up
1. **Markers** → pick source (Cluster root or Stack top), grouping
   column, click **Compute markers**.
2. Tighten `min_pct_in` / `max_padj` until top-N feels biological.
3. Pick a marker gene of interest and explore the FeaturePlot UMAP +
   spatial.
4. **Download markers CSV** for the supplement.
5. **Export** → session-report HTML; the Markers section will include
   the most-recent computation.

---

## Repo layout

```
xenium-panel-app/
├── app.R                      # `pkgload::load_all(".")` + `xenium_panel_app()`
├── R/
│   ├── globals.R              # paths, palettes
│   ├── load_panels.R          # reads /data/panel_audit/*.csv
│   ├── load_xenium.R          # ingest paths -> Seurat; qs2 cache
│   ├── panel_validate.R       # gene-set comparison
│   ├── panel_compare.R        # subpanel ⇄ data cross-tab
│   ├── cluster_pipeline.R     # the pure run_cluster_pipeline()
│   ├── clustree_module.R      # stability summary helpers
│   ├── markers.R              # presto::wilcoxauc wrapper
│   ├── export.R               # cluster CSV / markers CSV / report
│   └── module_*.R             # one Shiny module per tab
├── data/panel_audit/          # 49 subpanel CSVs + summary + ancillary
├── tests/testthat/            # unit tests + helper-test-data.R fixtures
├── cache/                     # qs2 cache (gitignored)
├── DESCRIPTION / NAMESPACE    # package layout
├── renv.lock / renv/          # reproducible R environment
└── CLAUDE.md / README.md      # build spec + this file
```

---

## Tests

```r
devtools::test()
```

Most tests run on a synthetic 800-cell × 500-gene fixture built in
`tests/testthat/helper-test-data.R` from real audit gene names — so
panel-validation logic is exercised. The `shinytest2` end-to-end smoke
test launches the real app via `AppDriver`, fetches the served HTML,
and asserts every navbar tab renders. It is gated on a detectable
Chrome binary (via `chromote::find_chrome`) and on `NOT_CRAN=true`
(set in `.Renviron`).

---

## Running on HiPerGator

```bash
module load R/4.4 hdf5/1.14
cd /path/to/XeniumPanelExplorer
R -e 'shiny::runApp("app.R", port = 4321, host = "127.0.0.1", launch.browser = FALSE)'
```

Tunnel the port to your laptop (VS Code remote, or
`ssh -L 4321:127.0.0.1:4321 hpg`).

Datasets > 1 M cells: avoid `SCTransform`; the pipeline already gates
that option above 200 k cells and falls back to `LogNormalize` with a
warning. For larger workflows the user has separate GPU tooling
(RAPIDS); this app is intentionally CPU-only.

---

## Known environment notes

- The codebase was developed on a host where
  `/usr/local/lib/R/site-library` contained packages built against an
  older R series. Compilation against R 4.6 sometimes failed because a
  transitive package's source build preferred the stale system include
  path. The fix used during development was
  `install.packages(c(...), lib = .libPaths()[1])` from the project
  root — that forces fresh builds into the renv project library so the
  stale system copies are bypassed at compile time. Document this for
  any host that ships `r-cran-*` apt packages.
- `pandoc` is **not** required. The Export tab's session report is
  built with `htmltools::save_html` directly; no rmarkdown / quarto /
  commonmark dependency on the render path.
- `presto`'s Seurat dispatch still passes the defunct `slot=` argument
  to `GetAssayData()` on SeuratObject ≥ 5.0 — this app calls
  `presto::wilcoxauc` directly on the data-layer matrix to avoid that.

---

## Out of scope (per CLAUDE.md §14)

CellChat / NicheNet, 3D spatial visualisation, multi-dataset integration
beyond Harmony on `orig.ident`, `Azimuth` reference mapping, the
Stellaromics / Pyxa platform, authentication / multi-user state.
