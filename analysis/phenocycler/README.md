# PhenoCycler single-cell analysis

## Current cell typing and islet trajectories

Open [06 · Distribution-fitted cell typing](notebooks/06_distribution_fitted_phenotyping.ipynb) for current pancreas and lymph-node cell identities, signed protein markers and full per-donor UMAPs. [07 · Per-donor islet trajectories](notebooks/07_islet_trajectories.ipynb) fits pseudotime from selected features, compares measured and excess-zero-adjusted inflammation, and runs Seurat scSLIDE on protein profiles within islets. If you open these notebooks from a separate checkout, set `PHENOCYCLER_DYNAMIC_RELEASE` to the local fitted release.

The [distribution-fitted PhenoCycler notebook](notebooks/06_distribution_fitted_phenotyping.ipynb) opens the November 16 private release. Its per-donor UMAPs include every eligible cell with a complete donor-specific protein profile. The notebook opens donor 115 and links to the full maps for other sections. The cohort UMAP uses a donor-balanced display sample. Set `PHENOCYCLER_DYNAMIC_RELEASE` to the private directory containing `scanpy/release.json`. Unlabeled distributions contribute to final assignments. Immune-lineage recovery declined in the matched-slide check. The November 13 notebooks remain available for comparison.

Open [Cell type phenotyping and annotation](notebooks/02_phenotyping.ipynb) for the cohort UMAP, marker heatmaps, dotplots, marker-pair scatterplots and endocrine quality checks. Change `DONOR`, `DISEASE` or `LEVEL` to inspect a section or hierarchy branch. `WEIGHTING = "donors"` averages measured donor means. Use `"cells"` to pool the measured cells. The default shows the deepest available subtype for every cell. CD68 or CD163 positivity establishes the macrophage parent. Galectin-3 (M2Gal3) separates LGALS3⁺ from LGALS3⁻ macrophages. CD4 contributes to T-cell subtyping. The main expression plots show proteins measured in every donor in a tissue. Supplementary plots group other proteins by the exact set of donors that measured them. The marker sets are fixed across disease groups and retained for donor views. Pancreas has 52 markers shared across 25 donors and 15 subgroup markers. All 53 pLN markers are shared across 19 donors. Each group has separate figures for cell identities and assignment QC or review outcomes. The availability table lists donor IDs for every marker group. Dotplots use Scanpy’s per-marker (`standard_scale="var"`) color scale. Dot area gives the measured fraction ≥0.10. The separate cell-type and QC heatmaps share one range for each marker within the same tissue, donor selection and donor-coverage panel. The lowest and highest measured means across the two figures map to 0 and 1. Single-marker panels show calibrated means on a numeric axis, with the same limits for cell types and QC. If an entire row falls in the highest color band, a quantitative profile shows the original calibrated means. Between-donor heatmaps use the range for the displayed donors of the selected cell type. Hatching identifies a constant measured marker and solid gray tiles identify unavailable measurements. The returned tables retain calibrated means for comparisons across panels.

| Notebook | Analysis |
| --- | --- |
| [01 · Quality control](notebooks/01_preprocessing.ipynb) | Marker availability, missing measurements and mapped-cell coverage |
| [02 · Phenotyping](notebooks/02_phenotyping.ipynb) | Cohort, group and donor signatures, hormone conflicts matched to named cells, and weak named-call support |
| [03 · Spatial analysis](notebooks/03_spatial_analysis.ipynb) | Native section coordinates, cell neighborhoods and independently fitted donor islet trajectories with measured core and peri-islet immune cells per mm² on a shared scale |
| [04 · Group comparisons](notebooks/04_group_comparisons.ipynb) | Donor subtype percentages, disease-group medians and individual donor points, with separate figures for low-frequency subtypes and dominant exocrine cells. Also includes disease-group signatures and a donor heatmap for a selected phenotype |
| [05 · Tissue comparisons](notebooks/05_tissue_comparisons.ipynb) | Paired pancreas and pLN sections, independent embeddings, full-panel protein profiles and a 19-donor audit of unmatched pLN candidates |

Notebook 04 groups pancreatic subtypes by the largest median across ND, AAb+ and T1D. Low-frequency groups have their own figures and percentage axes. Higher individual donor values appear on a labeled upper axis. A zero disease median can accompany a positive donor frequency. Per-donor plots in notebooks 04 and 05 use the selected donor’s observed percentages. Every percentage uses the full admitted cell count for that donor and tissue.

The private Senior installation includes the [pancreas subtype heatmap](rendered/pancreas_subtype_marker_heatmap.png) and [assignment QC heatmap](rendered/pancreas_assignment_qc_shared_heatmap.png).

| Pancreas frequency figures | Disease-group medians and individual donors |
| --- | --- |
| Cell identities | [Exocrine](rendered/pancreas_subtype_proportions_by_disease_dominant.png) · [1% to <10%](rendered/pancreas_subtype_proportions_by_disease.png) · [0.1% to <1%](rendered/pancreas_subtype_proportions_by_disease_below_1pct.png) · [0.01% to <0.1%](rendered/pancreas_subtype_proportions_by_disease_below_0p1pct.png) · [0.001% to <0.01%](rendered/pancreas_subtype_proportions_by_disease_below_0p01pct.png) |
| Assignment QC | [Immune, class unspecified](rendered/pancreas_assignment_qc_proportions_by_disease.png) · [Competing endocrine hormones](rendered/pancreas_assignment_qc_proportions_by_disease_below_0p001pct.png) · [Zero medians with positive donors, unresolved lineage](rendered/pancreas_assignment_qc_proportions_by_disease_zero_median_common.png) · [Zero medians with positive donors, endocrine below floor](rendered/pancreas_assignment_qc_proportions_by_disease_zero_median_below_0p01pct.png) |

## Measured data

Point `PHENOCYCLER_DATA` to the private release folder containing `typing_release.json`. A local `analysis/phenocycler/scanpy/` or `scanpy/` beside the notebooks is also discovered. The required release ID is `io60-phenocycler-20261113-no-epcam`. Every notebook checks the release ID and exclusion list, verifies the actual saved UMAP marker list, and compares the admitted population with the donor manifest. The pancreas hierarchy and its UMAP exclude EpCAM. The dedicated pLN candidate hierarchy also omits EpCAM. Its UMAP was refit without EpCAM using the revised candidate-stratified sample. Its marker calibration stays tissue specific.

| Input under the release folder | Contents |
| --- | --- |
| `pancreas_manifest.csv` and `pancreas/*.h5ad` | Full admitted sections with native cell IDs, routes, refinement states, raw aligned protein and calibrated measurements |
| `pancreas_umap.h5ad` and `pancreas_random.h5ad` | Stratified embedding sample and independent donor-balanced random sample |
| `pancreas_islets.h5ad` | Native islets, measured core inflammation and independently fitted donor graph coordinates and ranks |
| `donor_counts.parquet` and `donor_markers.parquet` | Full-population denominators, availability, comparability, measured cells and positive cells |
| `pancreas_refinement_qc.parquet` | Admitted donor, first-pass and final identities, refinement state/pass, and measured concordance of full-graph neighbors |
| `endocrine_conflict_qc.parquet` | Native object IDs, competing hormones, exclusive pixel counts, whole-interior median brightness and recorded decisions |
| `pln_umap.h5ad` and `pln_no_match_*.parquet` | Dedicated pLN candidate embedding, per-donor gate failures, thresholds and same-cell cluster evidence |
| `pln_lineage_review_cells.parquet` and `pln_unmatched_review_*.parquet` | Native cell IDs, supported broad immune review labels, unresolved reasons, full-section counts and 53-protein summaries |
| `pln_scimap_candidate_*.parquet` and `pln_compartment_*.parquet` | Continuous-score hierarchy, full-coverage proportions and protein summaries, and measurement-compartment sensitivity |
| `islet_inflammation_ranks.parquet`, `islet_inflammation_by_donor.csv` and `islet_inflammation_method.json` | Measured core/peri areas, immune-cell counts, absolute density, cohort and within-donor ranks, coverage and source hashes |

Notebook 05 defaults to the dedicated pLN rule candidates with separate broad review labels for previously unmatched cells whose own immune marker profile agrees with their historical cluster. `PLN_METHOD = "strict"` shows the original rule outcomes. `PLN_METHOD = "scimap"` shows the continuous-score sensitivity. Inspect the original marker evidence before interpreting these alternatives. The review applies the macrophage anchor and Gal3 subtype rule to supported pLN macrophage and myeloid candidates. The dedicated panel lacks CD38, so its five CD38-dependent B-cell and plasma-cell subtypes remain unassessable.

The endocrine quality checks in notebook 02 read three hormone columns and relevant annotations from each section. They compare cells with residual conflicts to named cells in the **same donors and same cytoplasmic measurement compartment**. A separate table reports native-interior area and brightness in conflict candidates. Named-cell summaries distinguish a subthreshold measured hormone, the direct-Single route, and an unavailable refinement bound. Raw protein intensities, calibrated scaled scores and native-interior comparisons retain their own units throughout the analysis.

See the individual trajectories in notebook 03 and the private [one-page-per-donor atlas](rendered/donor_islet_trajectories.pdf). Each donor has its own islet-composition graph, minimum-inflammation root, and measured density curves. The color and density scales are shared across the cohort.

Run the notebooks from this checkout, its notebook directory, or an upstream-style `notebooks/` directory. Keep the accompanying plotting helpers in the notebook directory. Expected dependencies include Scanpy, AnnData, pandas, PyArrow, h5py, NumPy, SciPy and Matplotlib.

The workflow follows the [original PhenoCycler analysis notebooks](https://github.com/smith6jt-cop/Phenocycler_Analysis/tree/ee170f13ebb892bc395bb289625df038ad0fadb1/notebooks).
