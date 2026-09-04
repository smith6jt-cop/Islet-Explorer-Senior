# Hypothesis notebooks (H1–H20)

Every notebook in `analysis/` follows the same template — load via
`data_utils.load_data()` or `load_all_data()`, aggregate per-image,
run the three-test statistical framework
(Kruskal-Wallis → Mann-Whitney with rank-biserial r → Spearman dosage
trend), and write figures + tables below
`$FOLLICLE_OUTPUT_ROOT/analysis/{figures,tables}/`. They are intended to be
re-executable top-to-bottom after `scripts/preflight.py` passes.

## Index

| # | Notebook | Question | Primary input |
|---|---|---|---|
| H1 | `H1_vessel_density.ipynb` | Small-vessel density by region × rs3184504 (raw + RedPulp-normalized) | `AnnotationsFinal.csv` |
| H2 | `H2_vessel_morphology.ipynb` | Vessel morphology — area, circularity, solidity, elongation — by genotype | `AnnotationsFinal.csv` |
| H3 | `H3_follicle_structure.ipynb` | Follicle count, size, fraction by genotype | `AnnotationsFinal.csv` |
| H4 | `H4_tissue_proportions.ipynb` | Tissue region proportions, white-pulp metrics by genotype | `AnnotationsFinal.csv` |
| H5 | `H5_follicle_vascularization.ipynb` | Per-follicle vessel assignment via cKDTree; scaling relationships | `AnnotationsFinal.csv` |
| H6 | `H6_follicle_reconstruction.ipynb` | Follicle reconstruction: GMM area filter, spatial clustering, WP composite nodules, GeoJSON masks | `AnnotationsFinal.csv` + `analysis/geojson/` |
| H7 | `H7_snp_screen.ipynb` | Multi-SNP screen: 39 SNPs × 5 metrics, BH FDR, platform sensitivity | `Groups.xlsx` |
| H8 | `H8_clinical_covariates.ipynb` | Partial Spearman correlations controlling for age, BMI, HbA1c, C-pep | `Groups.xlsx` |
| H9 | `H9_polygenic_risk.ipynb` | Unweighted PRS from 30+ risk-annotated SNPs vs imaging phenotypes | `Groups.xlsx` |
| H10 | `H10_multivariate.ipynb` | Multivariate per-donor feature matrix; PCA, platform diagnostics | `AnnotationsFinal.csv` |
| H11 | `H11_lda_deconfounded.ipynb` | Platform-deconfounded LDA (5 approaches), permutation testing, bootstrap CIs | `AnnotationsFinal.csv`, `Cells.csv` |
| H12 | `H12_spatial_bootstrap.ipynb` | Spatial-bootstrap CIs for genotype effects | `AnnotationsFinal.csv` |
| H13 | `H13_cropped_replication.ipynb` | Replication on cropped sub-regions to control for tissue-size confounding | `AnnotationsFinal.csv` |
| H14 | `H14_cell_boundaries.ipynb` | Cell-marker-derived Follicle/PALS boundaries (CD20/CD3e KDE) — 14-step pipeline | `Cells.csv` |
| H15 | `H15_single_cell_analysis.ipynb` | Single-cell clustering: arcsinh + z-score → PCA → Harmony → Leiden → cell-type gating | `Cells.csv` |
| H16 | `H16_cell_phenotyping.ipynb` | Cell-type phenotyping: signature scores, regional composition, neighborhood diversity | `Cells.csv` |
| H17 | `H17_perivascular_phenotyping.ipynb` | Perivascular cell-type enrichment by genotype | `Cells.csv` |
| H18 | `H18_HE_follicles.ipynb` | H&E cohort (n=27) follicle morphometry by rs3184504 | `27_donor_spleen_measurements1.csv` |
| H19 | `H19_meta_analysis.ipynb` | Cross-cohort meta-analysis pooling fluorescent (n=14) + H&E (n=27) | `H10_feature_matrix.csv` + `HE_per_donor.csv` |
| H20 | `H20_bcell_density_feasibility.ipynb` | Feasibility of the B-cell density endpoint: hierarchical negative-binomial count rate with a follicle-area offset | `Cells.csv` + region layer |

## Notebook details

### H1 — Vessel density
Counts SmallVessels per region per image, normalizes by region area, and
also expresses each region's density relative to that image's RedPulp
density (to control for tissue quality). Primary headline plot:
**Follicle SmallVessel density × rs3184504**.

### H2 — Vessel morphology
Per-image medians of vessel `Area µm^2`, `Circularity`, `Solidity`, and
`Elongation = Max diameter / Min diameter`. Looks for genotype-linked
shifts in shape distributions, not just counts.

### H3 — Follicle structure
Count of `Classification == "Follicle"` per image; per-image mean area;
follicle area as a fraction of total annotated tissue. Establishes the
follicle-size phenotype that H18 and H19 then replicate in the H&E cohort.

### H4 — Tissue proportions
Per-image fraction of each MAIN_REGION (Follicle, PALS, RedPulp,
Trabeculae); white-pulp aggregate `WP_fraction = (Follicle + PALS) /
total`; `Follicle_PALS_ratio` to detect B-cell-zone vs T-cell-zone
imbalance.

### H5 — Follicle vascularization
Uses `data_utils.assign_vessels_to_follicles` — cKDTree nearest-neighbor
assignment of follicle-parented SmallVessels to follicle centroids. Then
fits the scaling relationship vessel count ∝ follicle area^α and tests
whether α differs by genotype.

### H6 — Follicle reconstruction
For each image: fits a 2-component GMM to follicle areas, drops the
small-component (sub-follicle noise), spatially clusters survivors via
DBSCAN, and builds white-pulp composite nodules. Cross-checks against
the QuPath GeoJSON exports in `analysis/geojson/`.

### H7 — Multi-SNP screen
Reads the full `rs*` panel from `Groups.xlsx`, builds per-SNP dosage
matrices via `data_utils.load_snp_panel`, and runs Spearman across
5 imaging phenotypes for each of ~39 SNPs. Applies Benjamini-Hochberg
FDR. Sensitivity arm: re-runs the screen on Phenocycler-only samples
to test whether rs3184504 effect is platform-specific.

### H8 — Clinical covariates
Loads per-donor Age, Gender, Ethnicity, C-pep, HbA1c, BMI from
`Groups.xlsx` via `data_utils.load_clinical`. Runs partial Spearman
(rank-residualize covariates, then Pearson on residuals) for genotype
effects on each imaging phenotype.

### H9 — Polygenic risk
Builds an unweighted PRS as `sum(dosage) / (2 × N_genotyped)`, dropping
donors with fewer than 10 genotyped SNPs. Tests PRS vs each imaging
phenotype with Spearman.

### H10 — Multivariate feature matrix
Builds the per-donor wide matrix (`H10_feature_matrix.csv`, ~23
columns) via `data_utils.build_feature_matrix`. Runs PCA, plots
genotype trajectories, computes
`data_utils.platform_diagnostic` to flag platform-confounded features.

### H11 — Platform-deconfounded LDA
5 deconfounding approaches: (1) feature residualization via
`data_utils.residualize_platform`, (2) PC-only LDA, (3) hand-curated
clean-feature LDA, (4) shrinkage LDA, (5) single-cell-derived LDA.
Reports permutation-test p-values and bootstrap 95 % CIs.

### H12 — Spatial bootstrap
Resamples annotations within image grids to construct spatial-bootstrap
CIs for the rs3184504 effect, controlling for spatial autocorrelation
that the per-image aggregate doesn't capture.

### H13 — Cropped replication
Crops every image to a fixed tissue-size window and re-runs H1-H4 to
confirm that the genotype effects aren't driven by larger tissue
sections (which Phenocycler tends to have).

### H14 — Cell-marker boundaries (method **abandoned**; kept as the record)
Builds Follicle / PALS regions from cell-marker KDEs (CD20 → Follicle,
CD3e → PALS) rather than QuPath pixel classification. The 14-step
pipeline is documented in `analysis/H14_cell_boundaries.md`.

**This is not how regions are defined any more.** Two markers cannot match a
classifier trained on 25–33 channels: comparing z-scored marker fields labels
anything below CD20's mean as PALS whether or not CD3e is elevated, which gave
HDL098 a "PALS" whose CD3e median (40) was *lower* than its own red pulp (71).
The region layer (`analysis/regions/`) comes from the QuPath classifier
geometry instead. H14 still runs, and its side-by-side comparison — its
pixel-classifier panel now reads the region layer — is the evidence for why the
cell-derived route was dropped.

### H15 — Single-cell analysis
22-marker arcsinh + z-score normalization → 50-PC PCA → Harmony batch
correction (per-Image) → Leiden clustering → manual cell-type gating
on marker centroids. Runs genotype / region / distance / neighborhood
analyses with rs3184504 dosage trends.

### H16 — Cell-type phenotyping (fluorescent cohort)
Overclusters centroids into fine cell types, computes per-type
signature scores, regional cell-type composition by genotype, and
neighborhood diversity (Shannon entropy of nearest-neighbor types).

### H17 — Perivascular phenotyping
Distance-stratified cell-type composition within fixed bands around
SmallVessels. Boxplots + heatmap + stacked bar + spatial maps of
enriched perivascular cell types by genotype.

### H18 — H&E follicle morphometry
27-donor H&E cohort, balanced 9 / 9 / 9 on rs3184504. Tests five
metrics — follicle count, mean area, area fraction, circularity,
nearest-neighbor distance — and reports the original null scan
(closest p = 0.064 for area-fraction C/T vs T/T) alongside a
power-aware section (bootstrap 95 % CIs + MDE r ≈ 0.75 at n = 9 / 9).
Centroid-based spatial features (Clark-Evans R, area Gini, centrality)
appended via `scripts/append_H18_spatial_features.py`. The
power-analysis appendix lives in `scripts/append_H18_power_analysis.py`.
Reads `Spleen_rs3184504_Genotypes.xlsx` for H&E donor genotypes.

### H19 — Cross-cohort meta-analysis
Pools the fluorescent cohort (`H10_feature_matrix.csv`, n = 14) and
the H&E cohort (`HE_per_donor.csv`, n = 27) on four overlapping
follicle metrics. DerSimonian-Laird random-effects meta-analysis on
rank-biserial r per pairwise comparison + Stouffer-combined Spearman
dosage. Built programmatically by `scripts/build_H19_notebook.py`.
Headline: non-additive (heterozygote-dip) pattern in mean follicle
area — C/C vs C/T pooled r = +0.41 (p = 0.069); C/T vs T/T pooled
r = −0.43 (p = 0.064). Four donors overlap both cohorts
(HDL055, HDL063, HDL073, HDL086); sample-level consistency reported
separately from the pooled estimate.

## Statistical framework

All implemented in `analysis/data_utils.py`:

- `run_kruskal(data, value_col)` — Kruskal-Wallis H across genotype.
- `run_pairwise(data, value_col)` — Mann-Whitney U + rank-biserial r
  for all three pairs.
- `run_dosage_trend(data, value_col)` — Spearman ρ vs ordinal
  {C/C: 0, C/T: 1, T/T: 2}.
- `full_stats_table(data, value_col, label)` — glues the three into
  a single 5-row DataFrame per metric.

All three operate on per-image aggregates (one row per sample), not on
individual annotations. That rule is enforced by the helper functions
themselves — they call `.groupby("Sample", observed=True)` before
testing.

## Reproducing

Re-running drops outputs into the same `analysis/figures/` /
`analysis/tables/` paths and overwrites the previous run. Use
`scripts/run_notebooks.sh`, which runs them in five dependency tiers and forces
the inline matplotlib backend — this shell exports `MPLBACKEND=Agg`, under which
`plt.show()` silently emits nothing inside a Jupyter kernel and no figure is
embedded in the `.ipynb`.

Two things break idempotence, both deliberate:

- **The region layer is a persisted artifact.** `analysis/regions/` must be
  built before any notebook runs (`export_all_regions_geojson.groovy` →
  `radial_follicle_profiles.py` → `build_region_layer_qupath.py`), and it is
  gitignored, so a fresh clone has none. Notebooks fail with a `FileNotFoundError`
  naming the builder rather than silently proceeding.
- **Order matters across tiers.** H10 builds the feature matrix H8/H9/H19 read;
  H15/H16 write the cell-type tables H20 reads. Running a later tier first makes
  it consume the previous round's numbers, which is invisible in the output.

H20 additionally caches to `analysis/tables/H20_bcell_*.csv` and
`analysis/.h20_cache/`, and skips recomputation when they exist — delete them
when upstream data changes.
