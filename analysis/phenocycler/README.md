# PhenoCycler single-cell analysis

Open [Cell type phenotyping and annotation](notebooks/02_phenotyping.ipynb) for the cohort UMAP, marker heatmaps, dotplots, marker-pair scatterplots and endocrine quality checks. Change `DONOR`, `DISEASE` or `LEVEL` to inspect a section or hierarchy branch. `WEIGHTING = "donors"` averages measured donor means. Use `"cells"` to pool the measured cells. The donor comparison labels show how many cells support each assignment.

| Notebook | Analysis |
| --- | --- |
| [01 · Quality control](notebooks/01_preprocessing.ipynb) | Marker availability, missing measurements and mapped-cell coverage |
| [02 · Phenotyping](notebooks/02_phenotyping.ipynb) | Cohort, group and donor signatures, hormone conflicts matched to named cells, and weak named-call support |
| [03 · Spatial analysis](notebooks/03_spatial_analysis.ipynb) | Native section coordinates, cell neighborhoods and per-donor islet trajectories against measured inflammation |
| [04 · Group comparisons](notebooks/04_group_comparisons.ipynb) | Full-section proportions, disease-group signature plots and a donor heatmap for each selected phenotype |
| [05 · Tissue comparisons](notebooks/05_tissue_comparisons.ipynb) | Paired pancreas and pLN sections, independent embeddings and tissue-specific marker panels |

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
| `pln_umap.h5ad` | Dedicated pLN display embedding |

The endocrine quality checks in notebook 02 read three hormone columns and relevant annotations from each section. They compare cells with residual conflicts to named cells in the **same donors and same cytoplasmic measurement compartment**. A separate table reports native-interior area and brightness in conflict candidates. Named-cell summaries distinguish a subthreshold measured hormone, the direct-Single route, and an unavailable refinement bound. Raw protein intensities, calibrated scaled scores and native-interior comparisons retain their own units throughout the analysis.

Run the notebooks from this checkout, its notebook directory, or an upstream-style `notebooks/` directory. Keep `phenotype_signatures.py` alongside the notebooks. Expected dependencies include Scanpy, AnnData, pandas, PyArrow, h5py, NumPy, SciPy and Matplotlib.

The workflow follows the [original PhenoCycler analysis notebooks](https://github.com/smith6jt-cop/Phenocycler_Analysis/tree/ee170f13ebb892bc395bb289625df038ad0fadb1/notebooks).
