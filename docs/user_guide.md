# Islet Explorer — User Guide

Interactive web application for exploring pancreatic islet CODEX multiplexed imaging data from nPOD donors. Tracks insulin loss, immune infiltration, and cellular microenvironment changes across the ND → Aab+ → T1D disease progression.

**Production URL**: `http://<server-ip>:8080/islet-explorer/`

---

## Table of Contents

1. [Overview](#overview)
2. [Outlier Handling](#outlier-handling)
3. [Plot Tab](#plot-tab)
4. [Trajectory Tab](#trajectory-tab)
5. [Viewer Tab](#viewer-tab)
6. [Statistics Tab](#statistics-tab)
7. [Spatial Tab](#spatial-tab)
8. [Single-Cell Drill-Down](#single-cell-drill-down)
9. [AI Assistant](#ai-assistant)
10. [Data Pipeline](#data-pipeline)
11. [Deployment & Administration](#deployment--administration)
12. [Troubleshooting](#troubleshooting)

---

## Overview

### Data Source

The app loads from a single enriched H5AD file (`data/islet_explorer.h5ad`, ~70 MB) containing:

- **5,214 islets** from 15 nPOD donors (ND, Aab+, T1D)
- **31 protein markers** (mean expression per islet)
- **Pseudotime trajectory** (DPT with INS root, scVI-corrected)
- **21 phenotype proportions** (from single-cell phenotyping)
- **Donor demographics** (age, gender, autoantibody status)
- **62 neighborhood metrics** (peri-islet composition, immune infiltration, enrichment z-scores, distances)
- **QuPath groovy data** (targets, markers, composition, LGALS3)

If the H5AD is unavailable, the app falls back to `data/master_results.xlsx` with reduced functionality (no phenotype proportions, no demographics, no neighborhood metrics).

### Disease Groups

| Group | Color | Description |
|-------|-------|-------------|
| **ND** | Steel Blue (#4477AA) | Non-diabetic control donors |
| **Aab+** | Burnt Umber (#CC6633) | Autoantibody-positive (pre-T1D) |
| **T1D** | Forest Green (#228833) | Type 1 diabetic donors |

Colors are configurable via the **Donor Status Colors** palette selector (Paul Tol default, Bright, Okabe-Ito, Diverging). The selector is available on the Plot sidebar, Trajectory tab, and Spatial tab — all stay synced.

---

## Outlier Handling

The app applies a single, consistent outlier rule across every tab so that scatter plots, trend lines, donor-level means, and statistical tests all agree about which points are "extreme."

### The Rule

**For each donor status (ND, Aab+, T1D), an islet is flagged as an outlier when its robust |z| — the Iglewicz–Hoaglin modified z-score, `0.6745·(value − median)/MAD` — exceeds the value chosen on the per-tab "Outlier z-threshold" slider (default 3, range 0.5–10, step 0.5).** Median and MAD (Median Absolute Deviation) are used instead of mean and SD so the rule is robust to right-skew and detects extremes on **both** the high and the low side of the distribution. When MAD = 0 (zero-inflated groups), the rule falls back to classical mean/SD per group. Computing the threshold per-status (rather than across all islets) prevents the Aab+ and T1D heterogeneity from being mistaken for one giant ND tail. The slider is synced across Plot, Trajectory, and Spatial; Statistics has its own independent slider. Two cards use variants on the same idea:

- **Spatial Card A (Peri-Islet Phenotype Enrichment)**: uses **Tukey 1.5×IQR** per donor status, computed on non-zero values only. Many phenotypes are absent from most islets, so a z-score rule would over-flag the zero mode; the IQR variant tolerates that.
- **Spatial Card C (Phenotype Proximity to Islet)**: applies the |z| > slider rule to `log(1 + distance)` rather than raw distance, because distance distributions have a long right tail.

### The Control

A checkbox **"Exclude outliers"** plus an **"Outlier z-threshold"** slider lives on the Plot sidebar and is mirrored above the Trajectory scatter and inside the Spatial Neighborhood Analysis bar. Toggling the checkbox or moving the slider on any one of them updates the other two, so the Plot, Trajectory, and Spatial tabs always reflect the same setting. **Default: checkbox OFF, slider 3.0.** Hover the label for the full tooltip.

The **Statistics tab has its own dedicated outlier checkbox + slider** (Section 1: Configure Analysis). It uses the same z-score rule but operates independently of the Plot / Trajectory / Spatial synced controls — you can analyse Statistics with outliers excluded at a stricter threshold while keeping the Plot scatter at default, or vice versa. A small caption beneath the Statistics checkbox reports how many islets the most recent Run excluded.

### What you see when the filter is ON

- Excluded points still render on the Plot left scatter, the Plot right distribution panel, the Trajectory scatter, and the Spatial Card C box plot as **faint grey × markers** ("ghost points") so you can see *where* in the distribution the exclusions are landing.
- Plot, Trajectory, Spatial Card A, and Spatial Card C show a **count badge** in or beneath the title: e.g. `"… — 4 outliers excluded (|z| > 3 per donor status)"`. The badge always reflects the live slider value, and the count includes outliers on both ends of the distribution.
- Statistics tab's Methods Reference section names the rule and reports the exact count for the current Run.
- **All four analysis tabs** have an optional **"Show excluded-outlier table"** checkbox: Plot (in the sidebar), Trajectory (above the scatter), Spatial Card A (Peri-Islet Enrichment) and Spatial Card C (Phenotype Proximity) — each next to its own phenotype/metric selector. Each table lists the flagged islets (Case ID, donor status, value, z-score or fence distance) and is shown directly below the corresponding plot.
- Plot points and trajectory ghost-marker positions are stably jittered with a fixed seed, so toggling controls doesn't shuffle the cloud.

### What changes when the filter is OFF

- No points are flagged or hidden; ghost markers disappear; counts go to zero.
- A small banner appears above the Trajectory scatter noting that the Plot/Trajectory/Spatial filter is off; the Statistics tab shows its own yellow "Filter OFF" caption beneath its dedicated checkbox.
- Trend lines, donor-level means, and all statistical tests run on the full data including extreme values.

### Why this matters

With 5,214 islets from only 15 donors, a single extreme islet within one donor can dominate the group mean and pull mixed-effects estimates. Keeping the rule visible (count badge + ghost × markers) means you can audit at a glance whether a striking effect is being driven by one or two unusual islets.

---

## Plot Tab

The Plot tab is the primary exploration interface. It has a **sidebar** (left) and **main panel** (right).

### Sidebar Controls

- **Mode**: Switch between `Scatter` (feature vs diameter) and `Composition` (cell type proportions)
- **Feature/Composition selector**: Choose what to plot
  - *Scatter mode*: Select from markers (31 proteins) or targets (hormone intensities)
  - *Composition mode*: Three option groups (the previous "Immune Metrics" group of 5 hand-curated ratios was removed in Apr 2026 — use the Spatial tab's per-phenotype enrichment / proportion / distance views instead):
    - **Hormone Positivity (QuPath threshold)** — Ins_any, Glu_any, Stt_any
    - **Phenotype Proportions (scVI)** — 21 phenotype proportions from single-cell data
    - **Peri-Islet Proportions** — Proportions of each phenotype in the 20 μm expansion zone around the islet
- **Metric** (Composition / Microenvironment): Picks the y-axis unit. Labels spell out the denominator explicitly so "density" is never ambiguous — see [Density terminology](#density-terminology) below.
- **Region**: Filter by tissue region (core only, peri-islet, or all)
- **Donor Status**: Select disease groups to include (ND, Aab+, T1D)
- **AAb Filter**: Filter by autoantibody positivity (when available)
- **Age Range**: Slider to restrict donor age (H5AD only)
- **Sex**: Checkbox filter (H5AD only)
- **Min cells/islet**: Filter out islets with fewer than N total cells (core + peri). Default: 1 (no filter). Increase to reduce noise from poorly-measured small islets.
- **Diameter Range**: Restrict islet diameter in μm

### Main Panel

- **Scatter (left)**: Interactive Plotly scatter, diameter vs selected feature, with binned summary lines per donor status. Excluded outliers (see [Outlier Handling](#outlier-handling)) appear as faint grey × ghost markers.
- **Distribution (right)**: Violin or box plot of the same feature, grouped by donor status. Has its own **independent** "Exclude zero values" checkbox — toggling it no longer affects the left scatter. Outlier ghost × markers render here too.
- **Outlier Table**: Lists excluded outliers for the selected feature. Hover the "Show excluded-outlier table" checkbox label in the sidebar for the rule.
- **Click Interaction**: Click any data point to open the segmentation viewer (see [Single-Cell Drill-Down](#single-cell-drill-down))

### Tips

- The two "Exclude zero values" checkboxes (one per panel) are independent — set them differently when you want the violin/box to include zeros but the binned scatter to exclude them, or vice-versa.
- Point positions on both panels use seeded jitter, so the cloud doesn't shuffle when you toggle other controls.
- Use Cell Populations mode with "Peri-Islet Proportions" to visualize the cellular microenvironment around each islet. Use "Core+Peri — Phenotype Proportions" to pool both zones into a single per-islet value (cells from core + peri divided by the matching pooled cell count or pooled annotation area).
- The "Immune Metrics" group provides pre-computed ratios that highlight immune infiltration differences between disease stages.

### Density terminology

The word *density* is reused in several places with different denominators. Each axis label now spells the denominator out, but for reference:

- **Cells per islet (core / peri / core+peri)** — a raw count tallied inside the islet's core, the 20 µm peri-islet zone, or both pooled together. Not divided by area. Surfaced by the Metric = "Cells per islet" choice.
- **Cells per mm² (core / peri / core+peri area)** — phenotype cells divided by the QuPath-measured annotation area in mm². Core uses `core_region_um2` (the islet boundary annotation), peri uses `peri_region_um2` (the `_exp20um` ring annotation), and core+peri uses their sum. Surfaced by Metric = "Cells per mm² …" in Cell Populations mode for all three region variants.
- **{Nerve|Capillary|Lymphatic} area (% of region area)** — the only Microenvironment-mode metric. The *fraction* of the selected region (core, peri-islet, or core+peri) occupied by structure annotation, expressed as a percent: `100 × area_um2 / region_um2`. The y-axis names whichever class and region is selected (e.g. `"Nerve area (% of core area)"`). Counts of these structures are NOT exposed — Nerve / Capillary / Lymphatic are area-valued, and a single nerve may be segmented into several pieces, so object counts depend on segmentation morphology, not biology.
- **Cells per µm of distance (Spatial Card C)** — the signed-distance KDE on Card C is a kernel density rescaled by N, so the y-axis is cells per µm of distance from the islet edge — *not* a per-area density. Useful for comparing distance distributions across disease stages; do not read it as cells/mm².
- Log-scale features (e.g., cell counts) use pseudo-log transformation so zero values remain visible at y=0.

---

## Trajectory Tab

Visualizes the disease progression trajectory computed via diffusion pseudotime (DPT) with INS as the root anchor.

### Layout

The Trajectory tab uses a full-width layout (no sidebar). All controls are inline.

### Components

1. **UMAP Scatter** — Interactive Plotly UMAP colored by:
   - Donor status (default)
   - Selected feature (viridis inferno colormap, scaled to data min/max)
   - Pseudotime (DPT)

2. **Pseudotime Strip (Dominant Status × Feature Intensity)** — A single horizontal strip below the UMAP scatter, divided into 25 pseudotime bins. Each bin encodes two variables at once:
   - **Hue** = the *dominant* donor status in that bin, computed from per-status normalized density (count-in-bin divided by total islets for that status). Normalization corrects for cohort imbalance so dominance reflects enrichment relative to each status's pool, not raw count.
   - **Intensity** (saturation toward white) = the bin-mean of the selected feature, robust 5–95th percentile clipped. Floor at 50% saturation so even the dimmest bin remains a clear tint, never white.
   Hover for raw counts per status, normalized density per status, and the feature mean.

3. **Multi-Feature Heatmap** — Z-scored expression of user-selected markers along pseudotime bins. Default markers: INS, GCG, SST, CD3e, CD8a, CD68, CD45, HLADR. Blue-white-red diverging colormap, z-scores clamped to [-2.5, 2.5].

**Y-axis scale on the scatter** — Marker values are **scVI-normalized expression** (bounded roughly `0–1`), aggregated as the **mean across the cells in each islet**. They are *not* raw fluorescence intensity and *not* z-scores. Values almost never reach zero because (1) the scVI output for an individual cell rarely sits at exactly 0, and (2) averaging the cells in an islet smooths out any remaining zeros. The **multi-feature heatmap** below the scatter, by contrast, IS z-scored per marker across pseudotime bins (intentional — it standardizes markers for cross-marker comparison).

4. **Click Interaction** — Click any UMAP point to:
   - Highlight the selected islet
   - Open the segmentation viewer panel
   - Jump to the OME-TIFF Viewer tab for that donor

### Controls

- **Feature selector**: Choose which protein marker to display on the y-axis (drawn from the 31-marker panel only — no aggregate/derived metrics)
- **Trend lines**: None, Overall (single LOESS), or By Donor Status (ND/Aab+/T1D separate LOESS curves)
- **Color points by**: Donor Status (default), Donor ID, or Leiden cluster at any of 4 resolutions (0.3 / 0.5 / 0.8 / 1.0)
- **Pseudotime mode** (Apr 2026): Toggle between two pseudotime axes:
  - *Core only* (default) — DPT computed from raw scVI core means. Strongest INS signal (r ≈ −0.69).
  - *Core + Peri* — DPT computed from a 36-d input combining core scVI (10) + peri scVI (10) + 16 structural cell densities (Neural, Blood Vessel, Endothelial, Lymphatic in core and peri, weighted at α=0.15), with light Harmony (θ=0.1) on imageid. Strengthens immune/vascular marker correlations (CD3e, CD8a, CD31, CD34, PDPN) at a small cost to INS-driven span.
- **Point size by**: Cell Count (default), Islet Diameter, or Uniform
  - *Cell Count* sizes points by `sqrt(total_cells)` — small islets (3-10 cells) appear as tiny dots, large islets (100+ cells) as large dots. This honestly represents measurement quality.
- **Point transparency**: 0 (fully opaque) to 0.95 (highly transparent). Default 0.4. Moving the slider right increases transparency.
- **Point size**: Base point size slider (0.5-5.0)
- **Min cells/islet**: Filter out islets with fewer than N total cells. Default: 1 (no filter). Increase to focus on well-measured islets.

### Cell-Count-Weighted Trends

Trend lines are weighted by `log1p(cell count)` so that well-measured islets (many cells) drive the biological signal while noisy small islets (few cells) contribute proportionally less. This is important because 60% of islets have ≤10 cells — their aggregated expression values are inherently noisier.

Hover over any point to see its cell count in the tooltip (e.g., "Cells: 3"). Small islets at unexpected positions (e.g., ND islets at high pseudotime) typically have very few cells.

### Key Biological Insights

- **INS vs pseudotime**: Strong negative correlation (r = −0.694 in Core mode, −0.51 in Core+Peri mode) — insulin decreases along the trajectory
- **Disease ordering**: ND (early) → Aab+ (middle) → T1D (late) along pseudotime
- **GCG increase**: Glucagon expression increases as insulin decreases, reflecting alpha cell persistence
- **Aab+ heterogeneity**: Aab+ donors genuinely span a heterogeneous range of islet types (not all T1D-like); their pseudotime IQR is wider than ND or T1D — this is biology, not noise.

---

## Viewer Tab

Embedded OME-TIFF viewer (Avivator) for full-resolution multiplexed imaging data.

### Usage

1. Select a donor case from the dropdown, or click an islet in Plot/Trajectory to auto-navigate
2. The viewer loads the corresponding OME-TIFF with configurable channel overlays
3. Channel controls allow toggling individual markers and adjusting intensity ranges

### Requirements

- OME-TIFF files must be accessible at the configured URL
- Avivator runs client-side in the browser (WebGL required)

---

## Statistics Tab

Rigorous statistical analysis of the currently selected feature from the Plot sidebar. The Statistics tab **shares the Plot sidebar** — changing the feature, filters, or disease groups in the sidebar updates both Plot and Statistics simultaneously.

An **"About this tab"** banner at the top summarises what the tab does, the donor cohort (N=15), and the primary test (a donor-aware mixed-effects model with donor as a random intercept). Each section heading carries a small **(?)** icon — hover it for a one-line description of that section's contents.

The tab uses a **5-section narrative layout** that walks you through a logical analytical flow, with numbered headings and brief explanations at each step.

### Section 1: Configure Analysis

- Overview banner showing N islets, global p-value, effect size (η²)
- **Run Statistics** button triggers computation; **Download CSV** exports results
- Test type selector with inline explanations:
  - **Parametric** (ANOVA / t-test): assumes roughly normal data within groups
  - **Non-parametric** (Kruskal-Wallis / Wilcoxon): no normality assumption; based on ranks
- Additional controls: significance level (α), min cells/islet, bin width, diameter range. **Outlier removal** is controlled by a dedicated "Exclude outliers" checkbox + "Outlier z-threshold" slider local to the Statistics tab (z above the slider value per donor status, **off by default**, slider default 3.0). The controls are independent of the Plot / Trajectory / Spatial synced toggle, so you can analyse Statistics with outliers excluded while keeping them visible elsewhere. A small caption under the checkbox reports the most recent exclusion count (or a yellow "Filter OFF" notice). See [Outlier Handling](#outlier-handling) for the full rule.
- **Exclude zero-valued islets** (checkbox, default OFF): filters `value == 0` rows before every Statistics-tab test (mixed-effects, donor-level ANOVA / Kruskal, pairwise Cohen's d, per-bin tests, demographics regressions, covariate-adjusted model). Useful for zero-inflated features — rare phenotype proportions, peri-islet enrichment z-scores, or any metric where many islets register exactly 0. The overview banner shows an orange "Zero-valued excluded (N dropped)" badge when active; the forest-plot title and the Methods Reference paragraph cite the same count. **AUC is read from the Plot tab's binned summary and is not affected by this filter** — toggling re-runs the hypothesis-testing pipeline but leaves the AUC card untouched.

### Section 2: Primary Results

Three equal-height cards showing the core statistical results:

- **Hypothesis Testing** (left): Global test + pairwise comparisons table with Cohen's d effect sizes, 95% CIs, and BH-corrected p-values. The primary "Overall" row uses a mixed-effects model with donor as random intercept. The "Donor-level means" row averages islets within each donor (N=15) as a conservative cross-check.
- **Effect Size Forest Plot** (center): Visual display of pairwise Cohen's d with confidence intervals. Points colored by significance.
- **Area Under Curve** (right): Trapezoidal AUC by disease group with percentage change relative to ND baseline.

### Section 3: Size-Dependent Patterns

Stratified analysis testing whether group differences depend on islet size (effect modification):

- **Stratified Tests by Islet Diameter**: Heatmap of BH-corrected q-values (ANOVA and Kendall τ) computed within each diameter bin. Identifies which size ranges drive or lack group differences. P-values are corrected across bins to control false discovery rate.
- **Trend Analysis (Kendall τ)**: Measures correlation between disease stage (ND=0, Aab+=1, T1D=2) and the feature value within each size bin. τ > 0 means the feature increases with disease progression; τ < 0 means it decreases. Points colored by significance (red = significant, grey = NS).
- **Bin coverage rule**: The heatmap and the trend plot use the same rule — a bin is tested only if at least 2 donor-status groups each have ≥ 2 donors with finite values after filtering. Bins that fail this rule are silently omitted from both plots so the coverage stays in sync.

### Section 4: Confounders & Deeper Analysis

Full-width **Demographics** card (H5AD only, hidden for Excel fallback):

- **Donor Summary Table**: N islets (primary), N donors, age median/range, % male/female per disease group
- **Age vs Feature** (left): Islet-level scatter plot (all ~5,214 islets) with overall linear regression line and Pearson correlation. Subtitle notes that islets are correlated within donors.
- **Sex vs Feature** (right): Box plot of feature value by donor status, faceted by sex (Male/Female). Sex-stratified pairwise tests are **BH-corrected across all 6 contrasts** (3 pairs × 2 sexes); the subtitle reports the minimum adjusted q per sex (low-power sex × group cells, n=2-3 donors, are flagged but still corrected jointly).
- **Autoantibody Profile** (Aab+ only, hidden when no AAb data): Per-donor table showing which AAbs are positive, total AAb count, N islets, and mean feature value. Box plot of feature by number of positive autoantibodies (1, 2, 3+) with Kruskal-Wallis test.
- **Covariate-Adjusted Model**: Donor-level linear model adjusting for age and sex, testing whether donor status retains significance.

### Section 5: Methods Reference

Dynamic text (subdued grey background) describing all statistical tests, corrections, and assumptions. Adjusts based on test type and parameters selected in Section 1.

### Using the Statistics Tab

1. Select a feature in the Plot sidebar (e.g., `prop_CD8a Tcell` or `peri_prop_Macrophage`)
2. Switch to the Statistics tab
3. Adjust test parameters in Section 1 (Configure Analysis)
4. Click **Run Statistics** — all sections populate
5. Review the narrative flow: global test → effect sizes → size-dependent patterns → confounders

---

## Spatial Tab

Tissue-wide spatial visualization and peri-islet microenvironment analysis. Combines single-cell tissue scatter plots with Leiden clustering and neighborhood enrichment metrics.

### 3-Panel Layout

1. **Controls Sidebar** (left panel)
   - Donor selector (15 nPOD donors with disease status labels)
   - Color by: Phenotype (21 cell types) or Leiden cluster
   - Leiden resolution dropdown (0.3, 0.5, 0.8, 1.0) — visible when Leiden coloring selected
   - Region filter: All cells, Core + Peri only, or Core only
   - **Color background cells**: Toggle to show tissue background cells in phenotype colors (dimmed) instead of grey
   - **Phenotype filter**: Show/hide individual phenotypes with checkboxes. Use **All** / **None** links to quickly select or clear all types. Dynamically populated from the selected donor's cell types.
   - Donor status checkboxes
   - Phenotype and Donor Status palette selectors
   - Download CSV button

2. **Tissue Scatter** (center, 800px height)
   - Full tissue-wide scatter plot of ~177K cells per donor
   - Background: tissue cells in light grey by default, or colored by phenotype when "Color background cells" is checked
   - Foreground: core/peri cells colored by phenotype or Leiden cluster
   - Spatial coordinate convention: `coord_fixed()` + `scale_y_reverse()` (microscopy standard)
   - **Zoom**: Drag to select a region, then the view zooms in. Double-click or click **Reset Zoom** to restore the full view.
   - Uses ggplot2 `renderPlot` (not plotly) for performance at >100K points

3. **Leiden Panel** (right)
   - UMAP scatter of 5,214 islets colored by selected Leiden resolution
   - Stacked bar chart of mean phenotype composition per cluster
   - UMAP shows disease-stage separation (uses raw marker PCA visualization coordinates)
   - Hidden with message when Leiden data unavailable in H5AD

### Data Sources

- **Tissue CSVs**: `data/donors/{imageid}.csv` (15 files, ~78 MB total). Columns: X/Y centroids (μm), phenotype, cell_region (core/peri/tissue), islet_name.
- **Leiden clustering**: Stored in `islet_explorer.h5ad` .obs as `leiden_0.3`, `leiden_0.5`, `leiden_0.8`, `leiden_1.0` + UMAP coords `leiden_umap_1`, `leiden_umap_2`.
- **Neighborhood metrics**: ~100 columns merged into H5AD .obs (peri-islet composition, per-phenotype enrichment z-scores and distances, plus aggregate immune metrics).

### Neighborhood Analysis Cards

Below the tissue scatter and Leiden panel, three interactive analysis sections visualize the per-phenotype peri-islet neighborhood metrics. **Apr 2026 update:** all three cards now expose every phenotype uniformly — previous hand-curated 7-immune-type lists replaced with dynamic dropdowns populated from the data. Cards only appear when neighborhood data is available in the H5AD (not in Excel fallback mode).

**Global Controls** (toolbar above the cards):
- **Min cells/islet**: Filter out small islets with few cells (default: 1, i.e., no filter). Higher values (e.g., 50) reduce noise from poorly-measured islets and decrease NA rates in distance metrics.
- **Point size** and **Opacity**: Adjust scatter plot visualization across all cards.
- **Islet count**: Shows how many islets pass the current filters.

**Card A: Peri-Islet Phenotype Enrichment**
- **Left** (`Phenotype (peri enrichment z)` dropdown): Bar chart of donor-level mean/median enrichment z-score across disease stages for the selected phenotype. Shows Kruskal-Wallis p-value. Donor points overlay with summary bars (mean ± SEM or median ± IQR). Outlier exclusion follows the global toggle (Tukey 1.5×IQR per donor status, zero-aware). A **"Show excluded-outlier table"** checkbox next to the phenotype selector reveals the list of excluded islets directly below the bar chart.
- **Right** (`Phenotype (peri vs core)` dropdown): Scatter of peri-zone proportion vs core proportion for the selected phenotype, with dashed y=x diagonal. Sqrt-scale toggle and per-status trend-line option. Points above the line have higher core proportion than peri.
- *Interpretation*: Pick any phenotype (Macrophage, CD8a Tcell, Beta cell, Endothelial, etc.) to see its peri-zone enrichment trajectory across ND→Aab+→T1D.

**Card B: Phenotype Composition & Enrichment**
- **Left**: Grouped bar chart spanning ALL 20 phenotypes (drops "Unknown") × 3 disease stages. Toggle median/mean summary; `Clip extreme z > 5` is a **display saturation** that caps the visible axis so the chart isn't dominated by a few hugely enriched phenotypes — it is NOT outlier exclusion (see [Outlier Handling](#outlier-handling) for the actual rule). Error bars show IQR (median) or SEM (mean).
- **Right**: Heatmap — phenotypes (columns) × disease stages (rows). Numeric values annotated on each cell.
- **Region toggle**: Switch between:
  - *Peri-islet (enrichment z)*: Poisson z-scores comparing peri-islet vs tissue-wide. Diverging colorscale (blue = depleted, red = enriched). Default.
  - *Core (proportion)*: Raw cell type proportions in islet core. Sequential colorscale (white → red).
  - *Peri-islet (proportion)*: Raw cell type proportions in peri-islet zone. Sequential colorscale.
- *Interpretation*: z > 0 means that phenotype is enriched near islets relative to the tissue-wide proportion. Useful for spotting structural changes (Endothelial, Lymphatic, Neural) alongside immune infiltration.

**Card C: Phenotype Proximity to Islet**
- **Left** (`Distance to nearest` dropdown): Box plot of minimum distance (μm) from islet core centroid to nearest peri-zone cells of the selected phenotype. Includes "Immune (all)" aggregate plus every individual phenotype with `min_dist_*` data. Shows non-NA counts per group. Outlier exclusion follows the global toggle (z > 3 per donor status applied to `log(1 + distance)`); excluded points still render as faint × markers. A **"Show excluded-outlier table"** checkbox next to the metric selector lists the excluded islets (Case ID, distance, z-score) below the box plot. See [Outlier Handling](#outlier-handling).
- **Right** (`Phenotype` dropdown): KDE of signed distance from islet boundary for individual cells of the selected phenotype. Negative = inside islet (core), positive = outside (peri-islet). Dashed line at zero marks the boundary.
- *Interpretation*: NA values in the box plot indicate islets with zero cells of that phenotype in the peri-islet zone (biological, not a data error — common for rare phenotypes like CD4 Tcell, B cell). Increase the min cells/islet filter to reduce NAs.

### Data Coverage

- All 5,214 islets have peri-islet data (100% coverage)
- Donor 6533 has 191 islets (fully integrated after Parent annotation fix)
- Biological validation: T1D immune_frac_peri > Aab+ > ND (immune infiltration increases with disease progression)

---

## Single-Cell Drill-Down

Click any islet data point in the Plot or Trajectory tab to inspect individual cells.

### How It Works

1. Click a point in the scatter/UMAP plot
2. A segmentation panel appears below the plot with the title showing the islet name and donor info (case ID, disease status, age, sex)

### Panel Layout

- **Left sidebar**: Layer toggles (Single Cells, Peri Boundary, Structures), Color-by dropdown, Phenotype Palette, Show peri-islet cells checkbox. When no cell data is available, shows a boundary legend instead.
- **Center**: Segmentation map with GeoJSON boundaries and optional single-cell overlay
- **Right sidebar**: Cell composition bar chart and cell count table (visible when Single Cells is enabled)

### Single Cells View

- **Color by**: Phenotype (21 categorical colors) or any of 31 protein markers (viridis inferno continuous scale)
- **Show peri-islet**: Toggle checkbox to include/exclude cells from the 20 μm expansion zone
  - Core cells: filled circles
  - Peri-islet cells: open circles (shape 1)
- **Summary panel**: Horizontal bar chart showing phenotype composition of the islet
- **Count table**: Core vs peri cell counts by phenotype

### Data Source

Per-islet cell CSVs in `data/cells/` (5,214 files, ~203 MB total). Each file contains:
- X/Y centroid coordinates (μm, converted to pixel space for GeoJSON overlay)
- Phenotype label (1 of 21 types)
- Region (core or peri-islet)
- Cell/Nucleus area
- 31 protein marker expression values

### Availability

- All 5,214 islets have cell data
- If cell data is unavailable for a clicked islet, the panel shows "Boundaries" mode with a message

---

## AI Assistant

An embedded AI chat panel powered by the University of Florida Navigator AI Toolkit. Ask questions about your data, plots, statistics, or troubleshooting steps.

### How to Use

1. Type a question in the text area at the bottom of the AI panel (right side of the app)
2. Select a model:
   - **Navigator Fast (gpt-oss-20b)** — Quick responses for simple questions (default)
   - **Navigator Large (gpt-oss-120b)** — More capable reasoning model for complex analysis questions. Shows "Thinking..." while reasoning, then displays the answer.
3. Click **Send** or press Enter
4. Use **New Conversation** to clear chat history and start fresh

### Model Selection Tips

- Use **Fast** for quick lookups: "What markers are in the heatmap?", "How many islets are shown?"
- Use **Large** for analytical questions: "Why does immune infiltration increase in T1D?", "How should I interpret this AUC result?"
- The Large model takes longer to respond because it reasons through the problem before answering

### Token Budget

The app tracks cumulative token usage. If a budget limit is configured (via `OPENAI_TOKEN_BUDGET`), the chat will notify you when the budget is reached. Start a new browser session to reset.

### Requirements

The AI assistant requires:
- `KEY` and `BASE` environment variables configured in `.Renviron`
- The `httr2` R package installed
- Network access to the UF Navigator API endpoint

If credentials are missing, the chat panel displays a configuration message instead of responding.

---

## Data Pipeline

### Rebuilding the App Data

The canonical pipeline rebuilds `data/islet_explorer.h5ad` from source:

```bash
# Activate the Python environment
conda activate scvi-env

# Step 1: Reaggregate islets + trajectory + Leiden (~15 min)
python scripts/reaggregate_islets.py
# → islet_analysis/islets_core_fixed.h5ad (5,214 islets)
# → data/adata_ins_root.h5ad (+ pseudotime + UMAP)
# → islet_analysis/islets_core_clustered.h5ad (+ Leiden at 4 resolutions)

# Step 2: Compute neighborhood metrics (from 2.6M-cell single-cell H5AD)
python scripts/compute_neighborhood_metrics.py
# → data/neighborhood_metrics.csv (5,214 rows × ~100 cols, per-phenotype enrich_z_* and min_dist_*)

# Step 3: Extract per-islet cell CSVs (for drill-down viewer)
python scripts/extract_per_islet_cells.py
# → data/cells/*.csv (5,214 files, ~203 MB)

# Step 4: Extract per-donor tissue CSVs (for Spatial tab scatter, independent of islet filtering)
python scripts/extract_per_donor_tissue.py
# → data/donors/*.csv (15 files, ~78 MB)

# Step 5: Build enriched H5AD (trajectory + groovy + neighborhood + Leiden + donor metadata)
python scripts/build_h5ad_for_app.py
# → data/islet_explorer.h5ad (~70 MB)
```

### Full Lineage

```
CODEX_scvi_BioCov_phenotyped_newDuctal.h5ad  (2.6M cells, canonical single-cell)
  ↓ scripts/reaggregate_islets.py (min_cells=0, require_paired=True)
islets_core_fixed.h5ad  (5,214 islets, proteins + scVI embeddings + trajectory + Leiden)
  ↓ scripts/build_h5ad_for_app.py (+ groovy + neighborhood + donor metadata)
data/islet_explorer.h5ad  (complete app data, incl. Leiden clustering)
```

Branch pipelines:
- `scripts/compute_neighborhood_metrics.py` — peri-islet composition, immune metrics, enrichment z-scores, distances
- `scripts/extract_per_islet_cells.py` — individual cell CSVs for drill-down
- `scripts/extract_per_donor_tissue.py` — per-donor tissue CSVs for Spatial tab scatter
- Leiden clustering from `islet_analysis/islets_core_clustered.h5ad` (4 resolutions) — merged by `build_h5ad_for_app.py`

### Upstream Data Sources

| Source | Location | Description |
|--------|----------|-------------|
| Single-cell H5AD | `single_cell_analysis/CODEX_scvi_BioCov_phenotyped_newDuctal.h5ad` | 2.6M cells, scVI-embedded (note: scVI training was not batch_key-aware; downstream Harmony in reaggregate_islets.py corrects donor variance for clustering and combined-mode pseudotime), phenotyped |
| Groovy TSV exports | `panc_CODEX/results/groovy_exports/` | QuPath islet measurements (15 donors × 4 types) |
| GeoJSON boundaries | `data/json/*.geojson` | Islet segmentation polygons |
| Spatial lookup | `data/islet_spatial_lookup.csv` | Islet centroid coordinates |

See `data/DATA_PROVENANCE.md` for the complete lineage documentation including validation results.

---

## Deployment & Administration

### Architecture

```
nginx (port 8080)
  ↓ reverse-proxy /islet-explorer/
shiny-server (port 3838)
  ↓ symlink /srv/shiny-server/islet-explorer → app/shiny_app/
Shiny app (R worker process)
```

### Updating the App

Code changes take effect when shiny-server spawns a new R worker. No restart is needed for code-only changes.

If a stale R worker persists with old code:
```bash
# Find and kill stale Shiny R workers
ps aux | grep R | grep shiny
kill <PID>
```

Note: `sudo systemctl restart shiny-server` is not available — use the worker kill approach above.

### Local Development

```bash
cd app/shiny_app
Rscript -e 'shiny::runApp(".", port = 7777)'
# Then open http://localhost:7777 in browser
```

Port 7777 is for local testing only and is not accessible to end users.

### R Dependencies

```
shiny, shinyjs, plotly, ggplot2, dplyr, tidyr, readxl, sf, jsonlite,
RColorBrewer, scales, anndata, reticulate, lmerTest, lme4, emmeans
```

### Python Dependencies (pipeline only)

```bash
conda activate scvi-env
# scanpy, anndata, scvi-tools, sklearn, scib-metrics, scipy, pandas, numpy
```

---

## Troubleshooting

### App fails to load

1. Check that `data/islet_explorer.h5ad` exists and is ~70 MB
2. Check that R packages are installed: `library(anndata); library(reticulate)`
3. Check that Python is accessible via reticulate: `reticulate::py_available()`
4. Fall back to Excel: ensure `data/master_results.xlsx` exists (reduced functionality)

### No phenotype or demographic filters

The app is running from Excel fallback. Rebuild the H5AD:
```bash
conda activate scvi-env
python scripts/build_h5ad_for_app.py
```

### No neighborhood metrics in Plot selector

1. Verify `data/neighborhood_metrics.csv` exists (5,214 rows)
2. Rebuild the H5AD: `python scripts/build_h5ad_for_app.py`
3. Kill stale R workers and refresh the browser

### Single-cell drill-down not available

1. Verify `data/cells/` directory contains ~5,023 CSV files
2. If missing, regenerate: `python scripts/extract_per_islet_cells.py`

### Segmentation viewer shows "No GeoJSON found"

1. Check that `data/json/` or `data/gson/` has `.geojson` / `.geojson.gz` files
2. Case ID zero-padding: GeoJSON files use 4-digit padded IDs (`0112.geojson`), data uses unpadded (`112`). The app handles this automatically via `sprintf("%04d", ...)` fallback.

### Statistics tab shows no results

The Statistics tab uses data from the Plot sidebar. Ensure:
1. A feature is selected in the Plot sidebar
2. At least 2 disease groups are checked
3. The diameter range includes some islets

### Spatial tab tissue scatter is empty

1. Verify `data/donors/` directory contains 15 CSV files
2. If missing, regenerate: `python scripts/extract_per_donor_tissue.py`
3. If Leiden panel says "not available", rebuild the H5AD with Leiden data: `python scripts/build_h5ad_for_app.py` (requires `islet_analysis/islets_core_clustered.h5ad`)

### AI assistant shows error or no response

1. **"LLM key not found"**: Set `KEY=your-api-key` and `BASE=https://api.ai.it.ufl.edu` in `~/.Renviron` or `app/shiny_app/.Renviron`
2. **"key not allowed to access model"**: The selected model isn't available. Verify available models with `GET /v1/models`. Current models: `gpt-oss-20b` and `gpt-oss-120b`.
3. **Authentication error (401)**: Confirm the API key is valid and has no extra whitespace
4. **Timeout or empty response**: The Large model (120b) needs time for reasoning. Try the Fast model (20b) for quicker responses, or increase the timeout.
5. **"httr2 package required"**: Install with `install.packages('httr2')` and restart the app
6. **Debug mode**: Set `DEBUG_CREDENTIALS=1` in `.Renviron` and check the R console output for credential loading details

### Leiden UMAP shows a blob

The Leiden UMAP visualization uses raw marker PCA coordinates (same as the Trajectory tab). If it shows a blob, the clustered H5AD may have outdated scVI-based UMAP coordinates. Fix:
1. Copy visualization UMAP from `data/adata_ins_root.h5ad` into `islet_analysis/islets_core_clustered.h5ad`
2. Rebuild: `python scripts/build_h5ad_for_app.py`
3. Kill stale R workers and refresh
