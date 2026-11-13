# PhenoCycler single-cell analysis

Start with [Cell type phenotyping and annotation](notebooks/02_phenotyping.ipynb). It shows the shared protein UMAP, a marker heatmap, and disease-group and donor dotplots and scatterplots. Change `DONOR`, `DISEASE`, and `MARKERS` to explore a different section.

| Notebook | Analysis |
| --- | --- |
| [01 · Quality control](notebooks/01_preprocessing.ipynb) | Protein availability and missing measurements across sections |
| [02 · Phenotyping](notebooks/02_phenotyping.ipynb) | Cell identities and marker expression across the pancreas cohort |
| [03 · Spatial analysis](notebooks/03_spatial_analysis.ipynb) | Section maps, cell neighborhoods, per-donor islet pseudotime and measured inflammation |
| [04 · Group comparisons](notebooks/04_group_comparisons.ipynb) | Cell-type proportions and protein signatures by donor and disease group |
| [05 · Tissue comparisons](notebooks/05_tissue_comparisons.ipynb) | Pancreas and dedicated pancreatic lymph-node sections |

## Running the notebooks

Open them in a Jupyter kernel with Scanpy, AnnData, NumPy, pandas, PyArrow, SciPy, Matplotlib, h5py and IPython. On the analysis workstation choose **PhenoCycler (Scanpy)**. The local `rendered/` directory contains copies with the plots already displayed. A notebook opened from this directory, the source notebook directory, or the Senior repository root finds the same `analysis/phenocycler/scanpy/` data. On another installation, place the processed bundle there or set `PHENOCYCLER_DATA` to its directory.

The November 13, 2026 release contains 31,305,205 admitted pancreas cells from 25 donors. Per-donor files under `scanpy/pancreas/` contain the full sections. `pancreas_umap.h5ad` contains the shared display sample, `pancreas_random.h5ad` contains 4,000 randomly sampled cells per donor, and `pancreas_islets.h5ad` contains the native islets and their donor-specific trajectories. The donor counts and marker tables summarize the full admitted populations. The dedicated pLN embedding lives in `pln_umap.h5ad`. The local `scanpy/` directory contains the measured data and stays outside Git.

The analysis follows the [original PhenoCycler analysis notebooks](https://github.com/smith6jt-cop/Phenocycler_Analysis/tree/ee170f13ebb892bc395bb289625df038ad0fadb1/notebooks).
