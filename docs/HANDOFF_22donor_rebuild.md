# Handoff — 15→22 donor rebuild (2026-07-03)

Committed handoff so this can continue on another machine (e.g. HiPerGator) from git alone.
Local Claude auto-memory does **not** sync — this file is the source of truth.

## ⚠️ Data does not travel via git
`data/` is gitignored (~86 GB) and so are the raw inputs under `~/IO60panc2nd/`
(`Cellmeasurements.csv` 138 GB, `CellmeasurementsBatch2.csv` 60 GB, qptiffs, GeoJSONs).
**None of this is in the repo.** To run the pipeline elsewhere you must transfer the needed
`data/` subdirs and raw inputs separately. What git carries is only the **code + this plan**.

## ⚠️ Working rule (do not repeat this session's mistake)
**Never assume marker/cell-type frequencies or cross-donor uniformity. Donor biology varies.**
- Donor **6476** is the lone **pancreatitis** donor — its extreme CD3e (T-cell infiltration) is
  REAL, not an artifact/outlier to "correct" toward the cohort. A higher CD3e threshold there may
  be correct.
- **CD3e (T cells) is the sparser marker; CD163 (macrophages) is more abundant** (6476 excepted).
- Positivity is **RESTORE-threshold-defined**, never a raw-MFI floor (RESTORE's internal `>50` is
  an implementation constant, not a positive call).
- Justify method changes on **methodological** grounds (stability, representativeness), and validate
  a per-donor threshold/call against **that donor's own data** (its histograms, its negative-control
  population, spatial coherence) — not vs other donors.

## Donors (22 = 15 batch-1 + 7 batch-2)
New: **6450, 6505 = Aab+**; 6523, 6534, 6566, 6591, 6623 = T1D. Split: **ND 8 / Aab+ 3 / T1D 11**.
**Panel mismatch:** batch-1 has CD38, batch-2 has b-Catenin1 (physical panel change) → use the
**58 common markers**; handled in `build_singlecell_anndata.py` by marker-column intersection.
Metadata: `~/IO60panc2nd/donor_metadata_panc.xlsx`.

## Pipeline progress
| Step | Status |
|---|---|
| GeoJSON export (22) | ✅ `scripts/senior/export_new_geojsons.sh`; 6623 needed `~/IO60panc2nd/scripts/export_helpers/export_cells_geojson_safe.groovy` (drops 1 malformed polygon) → `data/redsea_scratch/geojson/` |
| `build_cells_parquet` (22) | ✅ `data/cells/donor_id=*` (both CSVs) |
| REDSEA (22) | ✅ `data/cells_redsea/donor_id=*` — 15 old REUSED (proved deterministic: re-run of 6539 matched June to 2 of 30.2M values), only 7 new computed |
| **RESTORE (22)** | 🔄 in-progress at handoff (`--robust`). **Re-run if incomplete:** `conda run --no-capture-output -n scvi-env python scripts/senior/restore_normalize.py --cells-dir data/cells_redsea --out-dir data/restore_redsea --gated-dir data/restore_gated_redsea --thresh-csv data/restore_thresholds_redsea.csv` |
| `assign_broad_lineage.py` (22) | ⬜ → `data/phenotype/broad/` |
| `build_phenotyped_h5ad.py` (NEW) | ⬜ GAP-#1 merge (broad lineage → `singlecell_protein_phenotyped.h5ad` `obs['phenotype']`, joined by `object_id`) |
| `train_scvi_senior.py` | ⬜ `batch_key="imageid"` (committed earlier), on REDSEA-corrected MFI |
| aggregate → markers/targets → composition/tissue → neighborhood → `compute_trajectory_senior.py` → extract cells → `build_viewer_assets.py` (auto donor list) → `build_app_h5ad.py` | ⬜ |

## RESTORE fix in this session — `scripts/senior/restore_normalize.py`
The flat 15k **pre-subsample** biased sparse-reference pairs. Now **subsample-AFTER-idx_select**:
full per-image data feeds `idx_select`; `patch_restore_cluster_cap` caps the post-idx_select cloud to
`--subsample` (15k) right before clustering. **Critical gotcha fixed:** RESTORE.py:221 runs a
figure-only `scipy.gaussian_kde` **O(n²)** over the full cloud even with `save_figs=False` → a
13-hour stall; `patch_restore_no_kde` neutralizes it (holoviews already stubbed). Result ~1 min/donor.
**SSC (SparseSubspaceClusteringOMP) is O(n²); 15k is its practical ceiling** (crashes ≳60k). NaN:
drop only all-marker-NaN REDSEA edge cells (raw QuPath has 0 marker-NaN), warn on partial-NaN.
Validated: thresholds subsample-stable (6476 CD3e 8-seed CV 6%), unchanged for well-behaved donors.

## Uncommitted-until-now, now in this commit
`restore_normalize.py`, `build_singlecell_anndata.py` (58-common), `build_viewer_assets.py`
(auto donors), `compute_trajectory_senior.py` (data-driven within-status η² gate incl. Aab+ now that
it's 3 donors), NEW `build_phenotyped_h5ad.py` + `archive_15donor.sh` (dry-run default; the 15-donor
archive was NOT yet run) + `export_new_geojsons.sh`, and doc updates. The trajectory batch-correction
work (`train_scvi_senior.py` batch_key etc.) was committed earlier in `61f2006`.

## Side note (this machine only)
`/data` (18 TB HDD) stalled under the OME-TIFF conversion writing ~115 GB files + default
`vm.dirty_ratio=20` on 503 GB RAM → 85 GB dirty froze `sync` (and the Firefox snap update). Not a
disk failure. Mitigate: `sudo sysctl vm.dirty_bytes=4294967296 vm.dirty_background_bytes=1073741824`.
Repo `data/` is on `/home` (nvme), unaffected.
