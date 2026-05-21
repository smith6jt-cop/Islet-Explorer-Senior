# Islet Explorer — Competitive Landscape & Feature Analysis

*February 2026*

---

## Part 1: Islet Explorer Feature Inventory

### Data Foundation
- **1,015 pancreatic islets** from nPOD donors across 3 disease stages (ND, Aab+, T1D)
- **54 CODEX protein markers** aggregated to islet level via scVI batch-corrected pipeline
- **Dual data source**: H5AD primary (47 MB, full features) with Excel fallback (graceful feature degradation)
- **Canonical lineage**: single-cell h5ad → fixed islet aggregation → trajectory rebuild → app-ready h5ad

### Tab 1: Plot (Interactive Scatter + Distribution)

| Feature | Details |
|---------|---------|
| 3 analysis modes | Composition (24 features: 3 hormone fractions + 21 cell-type proportions), Markers (% positive or counts), Targets (density or counts) |
| Islet size scatter | X = diameter (um), Y = metric, binned with user-controlled width (1–75 um), summary lines + error bars per donor group |
| Distribution plot | Violin or box plot by donor status, with optional individual points |
| Trend overlay | LOESS smoothing (overall or per-group) |
| Outlier detection | >3 SD flagged, displayed at 1.2x size, optional summary table |
| Scale options | Pseudo-log Y (zeros preserved at 0), log2 X, exclude-zeros toggle |
| Click-to-segmentation | Click any point → embedded QuPath GeoJSON segmentation map with islet/peri-islet/nerve/capillary/lymphatic polygons |
| Filtering | Donor status, autoantibody subtype (GADA/IA2A/ZnT8A/IAA/mIAA), diameter range, age slider, gender checkboxes |
| Color modes | By donor status (fixed palette) or by individual donor ID |
| Summary statistics | Mean±SE, Mean±SD, or Median+IQR |
| Export | Annotated CSV with full filter state metadata |

### Tab 2: Trajectory (Pseudotime + UMAP)

| Feature | Details |
|---------|---------|
| Pseudotime scatter | X = DPT pseudotime (INS-rooted, 0→1 = healthy→diseased), Y = any of 31 protein features |
| Donor status heatmap | 25-bin color bar showing ND→Aab+→T1D progression along pseudotime |
| UMAP pair | Donor status coloring + selected feature coloring (viridis inferno) |
| Multi-feature heatmap | Z-scored expression across pseudotime bins for user-selected markers (8 defaults: INS, GCG, SST, CD3e, CD8a, CD68, CD45, HLADR), blue-white-red diverging, clamped [-2.5, 2.5] |
| Feature categories | Grouped selector: Hormone (3), Immune (8), Vascular (4), Neural (3), Other |
| Point sizing | Uniform or scaled by islet diameter |
| Click-to-segmentation | Same embedded panel as Plot tab (shared root-level output) |
| Outlier handling | >3 SD removed from visualization, shown in expandable table |

### Tab 3: Viewer (OME-TIFF Multi-Channel Imaging)

| Feature | Details |
|---------|---------|
| Avivator integration | Embedded WebGL multi-channel fluorescence viewer via iframe |
| Channel config | Pre-configured 4 primary channels (INS/GCG/SST/DAPI) with full channel list |
| Cross-tab linking | Trajectory tab can force-select an image in the Viewer |
| Environment detection | Auto-adapts URL construction for reverse proxy vs. local dev |

### Tab 4: Statistics (7-Card Analytical Dashboard)

| Card | Content |
|------|---------|
| Overview banner | N islets, global p-value badge, effect size eta-squared, inline controls (test type, alpha, outlier removal, bin width, diameter range) |
| Hypothesis testing | ANOVA or Kruskal-Wallis (global) + pairwise t-test/Wilcoxon with BH correction, Cohen's d + 95% CI, forest plot |
| Per-bin heatmap | Plotly heatmap of -log10(p) across diameter bins x test type, with significance stars |
| Trend analysis | Kendall tau per bin, color-coded significance, overall direction statement |
| Demographics | Age scatter with per-group regression, gender-stratified tests, covariate-adjusted model (conditional on H5AD data) |
| AUC analysis | Trapezoidal AUC by donor group with % change interpretation |
| Methods | Dynamic text describing all tests, corrections, and assumptions used |
| Export | Full annotated CSV with all statistical results |

### Tab 5: AI Assistant ("PANC FLOYD")

| Feature | Details |
|---------|---------|
| Streaming chat | SSE token-by-token display, 15-turn history |
| Dual model | Fast vs. Large with auto-selection by message complexity |
| Token budget | Configurable limit with usage tracking |
| Error handling | 7+ error types with friendly messages (quota, rate limit, auth, DNS, timeout) |
| UF Navigator API | Supports university-hosted LLM endpoint |

### Cross-Cutting Features

- **Shared sidebar** between Plot and Statistics tabs (JS layout switching)
- **Root-level segmentation** shared by Plot + Trajectory (single output binding)
- **GeoJSON caching** across sessions within same R process
- **Synthesized union data** when GeoJSON exports are incomplete (core + band → union)
- **Zero-padded case ID fallback** for GeoJSON file lookup
- **Modular architecture** — 5 independent tab modules + 6 utility/helper files, thin `app.R` entrypoint

---

## Part 2: Competitive Landscape (25 Tools Surveyed)

### Tier 1: Pancreas/Islet-Specific Data Portals

| Tool | Data Type | Islet Quant | Disease Stats | Trajectory | Web |
|------|-----------|-------------|---------------|------------|-----|
| **nPOD Data Portal** | IHC WSI, metadata | No | No | No | Yes (gated) |
| **Nanotomy.org** | EM ultrastructure | No | No | No | Yes (open) |
| **Pancreatlas** | IMC, confocal, WSI images | No | No | No | Yes (open) |
| **HPAP PANC-DB** | scRNA-seq, CyTOF, IMC | No | No | No | Yes (mostly open) |
| **PanKBase** | Multi-modal (planned) | Planned | Planned | Planned | Under development |
| **Islet Regulome Browser** | Epigenomics, GWAS | No | No | No | Yes (open) |
| **Islet Gene View** | RNA-seq (188 donors) | No | T2D only | No | Yes (registration) |
| **TIGER** | RNA-seq eQTL (500+ samples) | No | T2D only | No | Yes (open) |
| **Islet eQTL Explorer** | eQTL variants | No | No | No | Yes (open) |

**Closest competitor: Pancreatlas** — open-access, multi-channel pancreas imaging, React web app with PathViewer. But it is an **image browsing** tool (no quantitative islet-level analysis, no disease-stage statistics, no trajectory).

**HPAP PANC-DB** is the most feature-rich portal but focuses on **scRNA-seq** (CellxGene viewer) — spatial CODEX data is download-only, not interactively explorable.

### Tier 2: General Spatial Biology Platforms

| Tool | Scope | CODEX Support | Analysis | Stats |
|------|-------|---------------|----------|-------|
| **Vitessce** | Framework (any spatial data) | Yes (viewing) | Linked views | No |
| **HuBMAP Portal** | 30+ organs, healthy atlas | Yes (CODEX) | Vitessce viewer | No |
| **CZ CELLxGENE** | scRNA-seq (93M+ cells) | No | DE, UMAP | Minimal |
| **scProAtlas** | 17.5M cells, 15 tissues | Yes (8 tech types) | Neighborhood analysis | No |
| **HTAN Portal** | Cancer multi-modal | CyCIF (not CODEX) | Minerva/CellxGene | No |

### Tier 3: Image Viewers (No Quantitative Analysis)

| Tool | Technology | Multi-Channel | Segmentation Overlay | Statistics |
|------|------------|---------------|---------------------|------------|
| **Avivator** | WebGL, OME-TIFF | Yes | No | No |
| **Minerva Story** | Narrative waypoints | Yes | No | No |
| **TissueViewer** | Server-side mixing | Yes | No | No |
| **cytoviewer** | R/Shiny, Bioconductor | Yes | Yes (masks) | No |
| **Mistic** | Python/Bokeh, local | Yes (tSNE layout) | No | No |

### Tier 4: Commercial/Desktop Tools

| Tool | Access | Key Gap |
|------|--------|---------|
| **Akoya MAV** | Commercial (ImageJ plugin) | Desktop-only, no web, no disease-stage stats |
| **NanoString Spatial Atlas** | Minerva viewer for pancreas | Spatial transcriptomics (not protein), healthy only, 4 donors |

### Tier 5: Raw Datasets (No Interactive Tool)

| Dataset | Data | Gap |
|---------|------|-----|
| **Damond 2019 IMC** | 1,581 islets, 35 markers, T1D | Download-only, no web viewer, IMC not CODEX |
| **Walker T2D CODEX (Zenodo)** | CODEX pancreas, T2D | Download + Pancreatlas viewing only |

---

## Part 3: Competitive Position Assessment

### What Makes Islet Explorer Unique

No other tool in the surveyed landscape combines all of the following in one interactive web application:

1. **Islet-level CODEX protein quantification** — Most tools either show raw images (viewers) or work at single-cell level (CellxGene). Islet Explorer aggregates to the biologically meaningful islet unit.

2. **T1D disease-stage progression analytics** — ND → Aab+ → T1D with proper statistical testing (ANOVA/KW, pairwise with BH correction, effect sizes). Other islet tools focus on T2D (TIGER, IGW) or healthy tissue (HuBMAP, NanoString Atlas).

3. **Pseudotime trajectory from protein imaging data** — scVI batch-corrected DPT with PAGA-UMAP. This is typically done with transcriptomics; doing it with CODEX protein data at islet level is novel.

4. **QuPath segmentation overlay** — Click any data point → see the actual tissue morphology with islet boundaries, vasculature, and nerves. Bridges the gap between quantitative analysis and spatial context.

5. **Integrated statistical dashboard** — 7-card layout with hypothesis testing, per-bin analysis, trend detection, demographics, AUC, and methods documentation. No other pancreas tool provides this level of statistical rigor in a web GUI.

### Feature Gap Analysis: What Others Have That We Don't

| Feature | Who Has It | Relevance |
|---------|-----------|-----------|
| **Single-cell resolution** | CellxGene, scProAtlas, cytoviewer | Our data is islet-aggregated; cell-level drill-down not available |
| **Multi-organ comparison** | HuBMAP, scProAtlas, HPA | We are pancreas-only (by design) |
| **Genomics/eQTL integration** | TIGER, Islet Regulome, IGW | We have protein data, not genomic |
| **Large donor cohorts** | TIGER (500+), IGW (188), nPOD (560+) | We have 1,015 islets but fewer unique donors |
| **Public open access** | Most tools above | We are currently internal (<server-ip>:8080) |
| **Spatial neighborhood analysis** | scProAtlas, scimap | We show segmentation but don't quantify cell-cell spatial relationships |
| **Multi-dataset harmonization** | PanKBase (planned), HCA | We serve one CODEX dataset |
| **Gigapixel whole-slide browsing** | Pancreatlas, nPOD (Aperio) | Our Avivator viewer shows individual OME-TIFFs, not whole slides |

### Overall Assessment

Islet Explorer occupies a unique niche: it is the only web-based tool that provides **quantitative, statistically rigorous, islet-level analysis of CODEX multiplexed protein imaging data across T1D disease stages**, with integrated trajectory analysis and spatial segmentation context.

The closest tools in the landscape are:
- **Pancreatlas** — similar in serving pancreas imaging data interactively, but image-centric (no quantitative analysis)
- **HPAP PANC-DB** — similar in multi-modal pancreas data with disease context, but scRNA-seq-centric (CODEX data is download-only)
- **Damond 2019 IMC dataset** — similar in islet-level T1D multiplexed imaging, but a raw dataset (no interactive tool)

**In terms of analytical depth for T1D islet biology, Islet Explorer has no direct competitor.**

### Potential Enhancement Opportunities (from landscape analysis)

1. **Public access** — Every major competitor is open-access; internal-only limits impact
2. **Spatial neighborhood analysis** — scProAtlas and scimap show cell-cell spatial relationships are increasingly expected
3. **Single-cell drill-down** — CellxGene-style exploration of individual cells within selected islets
4. **Genomics cross-reference** — Link protein findings to eQTL/GWAS (TIGER, Islet Regulome) for multi-modal interpretation
5. **Whole-slide context** — Pancreatlas-style whole-slide browsing alongside islet-level analysis
6. **Multi-dataset support** — If additional CODEX cohorts become available, harmonized cross-dataset queries
