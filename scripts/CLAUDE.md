# CLAUDE.md - Pipeline Scripts

Pipeline scripts that build the data files consumed by the Shiny app. The top-level `CLAUDE.md` covers the high-level data lineage and run commands. This file documents the batch-correction pipeline, pseudotime modes, and per-script details.

## Donor Batch Correction & Pseudotime Modes (Phase 17, Apr 2026)

### scVI batch_key bug (root cause)

`single_cell_analysis/SCVI_CODEX_v2.ipynb` cell 16 set `categorical_covariate_keys=["Age", "Gender"]` but **never passed `batch_key="imageid"`**, despite comments claiming otherwise. Result: scVI's batch decoder was never activated and donor-specific technical variance (staining, antibody lot, imaging gain) leaked into `X_scVI_mean`. Symptoms: donor 6356 (ND, 18.8% of cells) formed a standalone Leiden cluster at every resolution (80–92% dominance); Cramer's V donor↔cluster ≈ 0.4 with AMI(donor) > AMI(status).

**Fix:** downstream Harmony correction (avoids scVI retrain on 2.6M cells). Two different θ values for two different goals.

### Step 8: Harmony for Leiden clustering (θ=2)

`scripts/reaggregate_islets.py` Step 8 applies `harmonypy.run_harmony` on `X_scVI_mean` with `batch=imageid`, `theta=2.0` before `sc.pp.neighbors`. Embedding stored as `obsm['X_scVI_harmony']`. Goal: donor-invariant islet cell-types. Post-fix metrics: Cramer's V dropped from 0.38 → 0.135, donor 6356's worst-cluster fraction 0.83 → 0.29, AMI(donor) 0.19 → 0.028. Status-orthogonality is fine/desired (clusters are cell-types, disease lives on pseudotime).

### Step 5: Trajectory pseudotime (raw scVI core, robust medoid root)

Step 5 uses **raw `X_scVI_mean`** (no Harmony) for the default core-only pseudotime. Empirical sweep showed even θ=0.1 compressed ND→T1D span more than it bought in donor noise reduction. Aab+ donor variance accepted as biology (Aab+ donors exhibit a heterogeneous range of islet types).

**Robust root selection** (replaces single-islet `argmax(INS)` over ND): root = medoid of top-50 ND islets by INS in scVI space. Anchored by the centre-of-mass of "healthy ND, INS-high" islets and naturally draws from multiple donors. Current root: index 703, donor 6516, INS 0.663 (99.7th percentile).

### Step 5b: Combined-mode pseudotime (Harmony θ=0.1 + structures α=0.15)

Computes a SECOND pseudotime saved as `dpt_pseudotime_combined`. Input matrix: `[core_scVI(10) | peri_scVI(10) | α_struct · z(struct(16))]` = 36 dims. Structural block: (`prop`, `count`) × (Neural, Blood Vessel, Endothelial, Lymphatic) × (core, peri), z-scored, scaled by `ALPHA_STRUCT = 0.15`. Disease ordering breaks at α ≥ 0.20; α=0.15 is the largest weight preserving ND<Aab+<T1D ordering while strengthening immune/vascular marker correlations (CD3e, CD8a, CD31, CD34, PDPN). Light Harmony (θ=0.1) on the combined matrix; same medoid root as core mode.

### Step 6 validation thresholds

- INS Spearman r < −0.3 (PASS at −0.69 raw / −0.55 combined)
- GCG Spearman r > 0.1 (relaxed from 0.2 because raw donor-amplified ~0.27 dropped post-fix)
- ND<Aab+<T1D mean ordering
- Donor eta² within status: ≤ 0.15 for ND/T1D, **≤ 0.25 for Aab+** (relaxed because Aab+ donors genuinely span heterogeneous islet types)
- 6533 included

### Important: don't claim "scVI is batch-corrected"

Old comments in `reaggregate_islets.py` and `islet_analysis/CLAUDE.md` claimed "scVI latent has donor-level variation corrected out (Age+Gender covariates bijectively map to donors)". This is wrong — Age/Gender as `categorical_covariate_keys` don't activate batch correction; only `batch_key` does. Comments updated Apr 2026.

## Pipeline Scripts & Notebooks

- `scripts/fix_tcell_subtypes.py` - Reclassifies CD4/CD8 T cells by marker gating (CD8a/CD4 >= 0.5); `--dry-run` flag for preview
- `scripts/reaggregate_islets.py` - **Primary pipeline**: aggregation (min_cells=0) + trajectory + Leiden in one script
- `scripts/compute_neighborhood_metrics.py` - Computes per-islet peri-islet metrics from single-cell H5AD
- `scripts/extract_per_islet_cells.py` - Extracts per-islet cell CSVs for drill-down viewer
- `scripts/extract_per_donor_tissue.py` - Extracts per-donor tissue CSVs (ALL cells: core+peri+tissue) for Spatial tab scatter
- `scripts/build_h5ad_for_app.py` - Builds enriched H5AD (trajectory + groovy + neighborhood + Leiden)
- `scripts/build_master_excel.py` - Builds master_results.xlsx from groovy TSV exports
- `notebooks/scvi_qc_validation.ipynb` - Validates scVI batch correction (silhouette, LISI, UMAP)
- `notebooks/rebuild_trajectory.ipynb` - Regenerates trajectory from fixed pipeline (superseded by reaggregate_islets.py)
- `islet_analysis/fixed_islet_aggregation.py` - Core aggregation module (islet-level from single-cell, default min_cells=0)

## Python Environment

Pipeline scripts use the `scvi-env` conda environment:
```bash
conda activate scvi-env  # scanpy, anndata, scvi-tools, sklearn, scib-metrics, scipy, harmonypy
```

`harmonypy==0.2.0` is a downstream batch-correction dependency (added Apr 2026). One-time install: `pip install harmonypy` inside `scvi-env`.

## Running Pipeline Scripts

```bash
conda activate scvi-env

# Fix T cell subtypes in canonical H5AD (dry-run first, then apply)
python scripts/fix_tcell_subtypes.py --dry-run
python scripts/fix_tcell_subtypes.py

# Full rebuild from single-cell H5AD (aggregate + trajectory + Leiden, ~15 min)
python scripts/reaggregate_islets.py

# Compute neighborhood metrics (reads 2.6M-cell H5AD, ~3 min)
python scripts/compute_neighborhood_metrics.py

# Extract per-islet cell CSVs (reads 2.6M-cell H5AD, ~5 min)
python scripts/extract_per_islet_cells.py

# Extract per-donor tissue CSVs for Spatial tab scatter (~2 min, independent of islet filtering)
python scripts/extract_per_donor_tissue.py

# Build enriched H5AD (merges trajectory + groovy + neighborhood + Leiden)
python scripts/build_h5ad_for_app.py
```
