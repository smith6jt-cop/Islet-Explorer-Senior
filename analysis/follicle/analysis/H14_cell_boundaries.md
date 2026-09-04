# H14 — Cell-Marker-Derived Follicle & PALS Boundaries

> **Superseded, 2026-08-11.** This method is no longer how regions are defined.
> Two markers cannot match a classifier trained on 25–33 channels: comparing
> z-scored CD20 and CD3e fields labels anything below CD20's mean as PALS
> whether or not CD3e is elevated, which gave HDL098 a "PALS" whose CD3e median
> (40) was *lower* than its own red pulp (71). Regions now come from the QuPath
> classifier geometry — see `analysis/regions/`,
> `scripts/build_region_layer_qupath.py`, and `docs/HANDOFF.md`.
>
> The notebook still runs and is still worth reading: its side-by-side
> comparison is the evidence for why this route was dropped, and the
> pixel-classifier panel now draws from the region layer rather than from the
> superseded `analysis/geojson/` export. Everything below describes the
> cell-derived arm as built.

## Motivation

The pixel-classifier-based region boundaries used in H1–H5 are trained on
Phenocycler image characteristics and may not generalize to CODEX samples.
H14 derives Follicle and PALS boundaries directly from single-cell marker
expression (CD20 for B-cells, CD3e for T-cells), producing platform-invariant
region definitions. Vessel metrics are then computed within these cell-derived
boundaries, enabling unbiased cross-platform genotype comparisons.

## Data Sources

| Input | Description |
|-------|-------------|
| `Measurements/AnnotationsFinal.csv` | 499K annotation rows, 13 images (SmallVessel, LargeVessel, region annotations) |
| `Measurements/Cells.csv` | 22.6M single cells, 71 columns; CD20 and CD3e mean intensity used here |
| `Groups.xlsx` | Donor metadata, rs3184504 genotypes, clinical covariates |
| `analysis/regions/*.geojson` | the region layer, in µm — the pixel-classifier side of the IoU comparison (was `analysis/geojson/`, in pixels, rescaled by a local pixel-size table that omitted the three ForSH2B3 CODEX sections) |

## Samples

11 donors after exclusions (HDL018, HDL021 excluded for quality; HDL172 excluded — no genotype):

| Genotype | n | Samples |
|----------|---|---------|
| C/C | 3 | HDL011, HDL053, HDL055 |
| C/T | 5 | HDL043, HDL052, HDL070, HDL079, HDL086 |
| T/T | 3 | 1901HBMP004 (CODEX), HDL063, HDL073 (CODEX) |

---

## Analytical Pipeline

### Step 1 — Rectangular Crop Standardization

**Cell 3** | Crops large Phenocycler images to a 5000 x 3000 µm window centred
on the median annotation centroid. CODEX images that already fit within these
dimensions are used in full. Manual crop offset is applied for HDL011 (split
tissue).

- **Crop threshold:** x-extent > 5000 µm or y-extent > 3000 µm triggers cropping
- **Centering:** median of all annotation centroids (region-agnostic, not
  vessel-maximising)
- **Clamping:** window edges constrained to stay within data bounding box

**Purpose:** Ensures comparable tissue areas across platforms (CODEX sections are
physically smaller than Phenocycler whole-slide scans).

### Step 2 — Per-Image Percentile Thresholding

**Cell 4** | Two-pass streaming through the 11.8 GB `Cells.csv`:

**Pass 1 — Intensity histograms:** For each image, accumulates CD20 and CD3e
intensity histograms (2000 bins, range 0–65000) from all cells within the crop
window. Computes percentile thresholds from the cumulative histogram:
- CD20: p90 (top 10%) — selects strongly-expressing B-cells
- CD3e: p95 (top 5%) — T-cells are more diffusely distributed, stricter cutoff
  prevents boundary inflation

**Purpose:** Percentile thresholds are platform-invariant (they adapt to each
image's intensity distribution), unlike Otsu which is sensitive to bimodality
differences between CODEX and Phenocycler.

### Step 3 — Marker-Positive Cell Extraction

**Cell 5** | Pass 2 through `Cells.csv`: extracts (x, y) coordinates for cells
exceeding the per-image thresholds. Also stores raw cell data arrays
`(x, y, cd20, cd3e)` per image for downstream sensitivity analysis.

**Outputs per image:**
- `cd20_coords`: coordinates of CD20-high cells (follicle markers)
- `cd3e_coords`: coordinates of CD3e-high cells (PALS markers)
- `cell_xy`: all cell coordinates (for KDE extent)
- `raw_cell_data`: full intensity arrays (for threshold sweep)

### Step 4 — KDE Smoothing & Polygon Extraction

**Cell 6** | For each image, converts marker-positive cell coordinates into
smoothed density surfaces, then extracts region polygons via contouring.

**Adaptive bandwidth:** Per image, per marker:
```
BW = 4 × median(5th-nearest-neighbor distance)
```
Clamped to [25, 100] µm. Images with tightly clustered cells get smaller
bandwidths (sharper boundaries); dispersed images get larger bandwidths
(smoother boundaries).

**KDE procedure:**
1. Bin marker-positive cell coordinates onto a 10 µm grid
2. Pad grid by 2×BW in each direction (eliminates straight-line edge artifacts)
3. Gaussian-smooth with sigma = BW / grid_spacing
4. Extract contours at a fraction of max density:
   - **Follicle (CD20):** 10% of max density
   - **PALS (CD3e):** 25% of max density (stricter — only densest T-cell cores)
5. Convert contours to Shapely polygons; discard polygons < 8400 µm² (GMM-derived
   minimum from H6)

**PALS trimming:** Raw CD3e polygons are differenced against the Follicle polygon
union, so PALS and Follicle regions are non-overlapping.

**WhitePulp:** Defined as the union of all Follicle + PALS polygons per image.

**Note on polygon count:** The number of polygons per image reflects how
fragmented or consolidated the marker-positive territory is within the crop, not
a count of individual histological follicles. Multiple true follicles may merge
into one polygon if their B-cell zones are close relative to the bandwidth, or a
single follicle may appear as one polygon.

### Step 5 — Vessel Assignment & Metric Computation

**Cell 7** | Assigns SmallVessel and LargeVessel annotations from
`AnnotationsFinal.csv` to cell-derived regions using spatial point-in-polygon
queries (Shapely STRtree).

**Spatial assignment:** Each vessel's centroid is tested against the Follicle,
PALS, and WhitePulp polygon sets. Vessels outside all regions are labelled
"Unassigned."

**LargeVessel exclusion:** LargeVessels (arterioles with visible lumen and
smooth muscle) share CD31/CD34 markers with SmallVessels but are not
microvessels. Their area is subtracted from the region denominator:
```
Effective_Area = Region_Polygon_Area − LargeVessel_Area_Inside_Region
```

**Three vessel metrics computed per image per region:**

| Metric | Formula | Biological meaning |
|--------|---------|-------------------|
| Vessel Area Fraction | SV_area / Effective_Area | Total vascular coverage (dimensionless) |
| Vessel Count Density | SV_count / Effective_Area_mm² | Vessel frequency per unit area |
| Mean Vessel Size | SV_area / SV_count | Individual vessel caliber (µm²) |

Area Fraction = Count Density × Mean Vessel Size. This decomposition
distinguishes whether genotype differences arise from more/fewer vessels versus
larger/smaller vessels.

**Output:** `H14_vessel_area_fractions.csv` — 33 rows (11 images × 3 regions).

### Step 6 — Spatial Visualizations

**Cell 8 — Region boundary maps:** Per-image subplots sorted by genotype
(C/C → C/T → T/T) showing:
- Gray cell dots (subsampled for speed)
- Filled Follicle (blue) and PALS (orange) polygons
- SmallVessel circles (white outlines)
- LargeVessel circles (red, excluded from fractions)
- Colored subplot borders indicating genotype group

**Output:** `H14_vessel_mask_regions.png`

**Cell 9 — Crop context maps:** Full-tissue views showing where the 5000×3000 µm
crop window sits on each image's annotation field. Confirms crops capture
representative tissue.

**Output:** `H14_crop_context.png`

### Step 7 — Vessel Area Verification

**Cells 10–12** | Two-stage verification that vessel area fractions are computed
correctly.

**Visual verification (Cell 11):** Two figures, both sorted by genotype:
- **Full crop views:** SmallVessels drawn as filled circles colored by region
  assignment (blue = Follicle, orange = PALS, gray = unassigned). LargeVessels
  in red. Flagged large SmallVessels (> 1000 µm²) marked as yellow diamonds.
  Per-image text overlay shows area fraction, vessel count, and LargeVessel
  exclusions.
- **Zoomed follicle insets:** Zoomed to the largest follicle polygon per image
  with individual vessel circles at true scale.

**Outputs:** `H14_vessel_area_verification_full.png`,
`H14_vessel_area_verification_zoom.png`, `H14_flagged_large_smallvessels.csv`

**Numerical cross-check (Cell 12):** Independently recomputes vessel area
fractions from scratch per image and compares against the results table.
Maximum absolute discrepancy should be ~0 (machine epsilon).

**Output:** `H14_vessel_fraction_crosscheck.csv`

### Step 8 — Pixel-Classifier Boundary Comparison

**Cell 13** | Loads QuPath GeoJSON exports of pixel-classifier-derived Follicle
and PALS annotations and computes Intersection-over-Union (IoU) against the
cell-derived polygons. Provides a quantitative measure of how well the two
boundary methods agree, and a visual overlay for 4 representative images.

**Outputs:** `H14_iou_comparison.csv`, `H14_boundary_comparison.png`

### Step 9 — Marker Enrichment Validation

**Cell 14** | Validates that the cell-derived boundaries correctly capture the
expected cell types:
- CD20-high cells should be predominantly inside Follicle polygons
- CD3e-high cells should be concentrated inside PALS polygons

Reports the fraction of marker-positive cells falling inside each region.
Expected: CD20-in-Follicle > 50%; CD3e-in-PALS > 30% (lower because T-cells are
more diffuse and PALS polygons are conservatively trimmed).

**Output:** `H14_marker_enrichment.csv`

### Step 10 — Donor Demographics

**Cells 15–16** | Loads clinical metadata from `Groups.xlsx` for the 11 donors
and produces:

- **Individual donor table:** Sample, genotype, platform, age, gender, ethnicity,
  C-peptide, HbA1c, BMI
- **Summary by genotype:** n, mean ± SE for age and BMI, female percentage,
  platform breakdown

**Figure (3 panels):**
1. Age by genotype — bar chart with mean ± SE and individual data points
2. BMI by genotype — bar chart with mean ± SE and individual data points
3. Gender & ethnicity — side-by-side stacked bars per genotype

**Outputs:** `H14_donor_characteristics.csv`, `H14_donor_demographics.png`

### Step 11 — Region Morphometry by Genotype

**Cells 17–18** | Compares the cell-derived region structure itself (not vessels
within regions) across genotypes. Metrics computed from KDE polygons per image:

| Metric | Description |
|--------|-------------|
| Follicle_area_mm² | Total follicle polygon area |
| Follicle_n_polys | Number of distinct follicle polygons |
| Follicle_mean_poly_mm² | Mean area per follicle polygon |
| PALS_area_mm² | Total PALS polygon area |
| PALS_n_polys | Number of distinct PALS polygons |
| WP_area_mm² | Total white pulp area (Follicle + PALS) |
| Follicle_frac_of_WP | Follicle area / WhitePulp area |
| Follicle_PALS_ratio | Follicle area / PALS area |

Each metric tested with Kruskal-Wallis, pairwise Mann-Whitney (with rank-biserial
effect size), and Spearman dosage correlation (C/C=0, C/T=1, T/T=2).

**Outputs:** `H14_region_morphometry.csv`, `H14_region_morphometry_stats.csv`,
`H14_region_morphometry_genotype.png`

### Step 12 — Genotype Analysis of Vessel Metrics

**Cell 19** | Core genotype comparison: tests all three vessel metrics (area
fraction, count density, mean vessel size) across all three regions (Follicle,
PALS, WhitePulp) — 9 metric × region combinations.

**Statistical framework** (per metric × region):
- **Omnibus:** Kruskal-Wallis H-test (3 genotype groups)
- **Pairwise:** Mann-Whitney U with rank-biserial effect size (3 pairs)
- **Gene dosage:** Spearman correlation with ordinal genotype (C/C=0, C/T=1, T/T=2)
- **Unit of analysis:** per-image aggregates (n=11), not individual annotations

**Figure:** 3 × 3 grid of boxplots (rows = metrics, columns = regions) with
individual data points marked by platform (circle = Phenocycler, square = CODEX).

**Outputs:** `H14_genotype_stats.csv`, `H14_genotype_vessel_metrics.png`

### Step 13 — Platform Consistency Check

**Cell 20** | Mann-Whitney U tests comparing CODEX (n=2) versus Phenocycler (n=9)
for all three vessel metrics across all three regions. Tests whether the
cell-derived boundaries successfully eliminate the platform confound that affected
pixel-classifier-based analyses.

> **Corrected 2026-08-05.** This step previously cited "H1 found borderline
> significant platform differences for Follicle density, p=0.055". That value is
> not supported by H1 and appears nowhere in the notebook or its outputs. The
> actual H1 platform term is **p = 0.808** (LMM, n = 2169 follicle annotations
> across 14 samples; `analysis/tables/H1_statistical_tests.csv`). The nearest
> real figure is H8's p = 0.019 for `Follicle_norm_density`, which is a different
> metric in a different notebook. Note also that H14's n=2 CODEX donors are the
> single-cell cohort, where both CODEX donors are T/T — see
> `analysis/PLATFORM_CONFOUND.md` for why that makes platform adjustment
> non-identifiable here, unlike in H1's n=14 annotation cohort.

**Output:** `H14_platform_check.csv`

### Step 14 — Sensitivity Analysis

**Cell 21** | Sweeps the CD20 percentile threshold (80, 85, 90, 95) and
recomputes follicle boundaries, adaptive bandwidth, and vessel metrics at each
setting. Tests whether the results are robust to the choice of threshold.

For each percentile × image: recomputes CD20-positive cells, adaptive bandwidth,
KDE, follicle polygons, and vessel assignment from scratch.

**Tracked metrics across thresholds:**
- Adaptive bandwidth (should decrease with stricter percentile)
- Follicle polygon count and area
- Vessel area fraction, count density, and mean vessel size

**Outputs:** `H14_percentile_sensitivity.csv`, `H14_sensitivity_analysis.png`

---

## Key Design Decisions

| Decision | Rationale |
|----------|-----------|
| Percentile (not Otsu) thresholding | Otsu depends on histogram bimodality which varies by platform; percentile adapts to each image |
| Separate CD20/CD3e percentiles (p90/p95) | T-cells are more diffusely distributed; stricter selection prevents PALS boundary inflation |
| Adaptive bandwidth from k-NN | Fixed bandwidth over-smooths sparse images and under-smooths dense ones |
| Different contour fractions (10%/25%) | CD3e density surfaces are flatter; higher threshold captures only the PALS core |
| PALS differenced from Follicle | Prevents double-counting in the transitional marginal zone |
| Vessel area fraction (not count density alone) | Dimensionless, robust to vessel segmentation granularity; decomposed into count density × mean size |
| LargeVessel area excluded from denominator | Arterioles (CD31+/CD34+) are not microvessels; their area should not inflate the effective region area |
| 5000 × 3000 µm crop | Matches CODEX tissue extent; prevents Phenocycler whole-slide images from dominating area comparisons |
| Minimum polygon area 8400 µm² | From H6 GMM-based follicle size filter; removes noise polygons |

## Output Summary

### Tables (`analysis/tables/`)

| File | Contents |
|------|----------|
| `H14_vessel_area_fractions.csv` | Per-image per-region vessel metrics (33 rows) |
| `H14_vessel_fraction_crosscheck.csv` | Independent recomputation verification |
| `H14_iou_comparison.csv` | Pixel-classifier vs cell-derived boundary IoU |
| `H14_marker_enrichment.csv` | Marker-positive cell fractions inside regions |
| `H14_donor_characteristics.csv` | Clinical metadata for the 11 donors |
| `H14_region_morphometry.csv` | KDE polygon morphometry per image |
| `H14_region_morphometry_stats.csv` | Genotype tests on region morphometry |
| `H14_genotype_stats.csv` | Genotype tests on vessel metrics (3 metrics × 3 regions) |
| `H14_platform_check.csv` | CODEX vs Phenocycler Mann-Whitney on all metrics |
| `H14_percentile_sensitivity.csv` | Threshold sweep results |
| `H14_flagged_large_smallvessels.csv` | SmallVessels > 1000 µm² (potential misclassifications) |

### Figures (`analysis/figures/`)

| File | Contents |
|------|----------|
| `H14_vessel_mask_regions.png` | Cell-derived boundaries + vessel outlines, sorted by genotype |
| `H14_crop_context.png` | Full-tissue views with crop window overlay |
| `H14_vessel_area_verification_full.png` | Full crop vessel area verification, by genotype |
| `H14_vessel_area_verification_zoom.png` | Zoomed follicle insets, by genotype |
| `H14_boundary_comparison.png` | Pixel-classifier vs cell-derived boundary overlay |
| `H14_donor_demographics.png` | Age, BMI, gender/ethnicity by genotype |
| `H14_region_morphometry_genotype.png` | Region morphometry boxplots by genotype |
| `H14_genotype_vessel_metrics.png` | 3×3 vessel metric boxplots by genotype |
| `H14_sensitivity_analysis.png` | CD20 percentile sweep metrics |
