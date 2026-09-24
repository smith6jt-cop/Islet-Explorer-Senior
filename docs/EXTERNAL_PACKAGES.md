# External and auxiliary packages — the index

**Status:** written 2026-09-02; paths corrected and the typing rows re-stated 2026-09-18. This is
an index, not a second copy of the pins: each pin lives in
exactly one place (a `pyproject.toml` extra, a lane README, or an adapter's import-time assertion) and
this table points at it. The lane that changes a pin updates this row. Companion documents:
`WORKFLOW_2026-09-02_MAP.md`, `PLAN_2026-09-02_AUXILIARY_PACKAGES.md`.

## The convention

- **pip-installable** → pinned **by commit** in a `core/phenocycler/pyproject.toml` optional-dependency
  group (`[dct]`), because a version string can stay `0.1.0` across incompatible commits.
- **not pip-installable** (no packaging metadata, or the code lives in notebooks) → a gitignored
  checkout under `/external/<name>` at a pinned commit, **asserted at import** by the lane's adapter so a
  drifted checkout stops the lane instead of changing a number (`[uniform]`, `[utag]`, KINTSUGI).
- **R packages** → the commit recorded in the lane README or its `setup_envs.sh`.
- Nothing provisions `/external/` on a fresh clone today (Track E.1 of the plan adds it to `setup.sh`).

## In use

| package | role (lane) | pin site → assertion | recorded finding |
|---|---|---|---|
| UniFORM — fork `smith6jt-cop/UniFORM` (unmodified) | cross-donor alignment, celltype Step 2 | `core/phenocycler/scripts/celltype/README.md` § Step 2; `pyproject.toml [uniform]` (deps only); `external/UniFORM @ 7c33396d12e5a19dc0b36e46a94c3293de0cd814` → `uniform_adapter.py:83` | 44/56 markers under the 0.10 bar; 12 shape-limited; CD57 refused as a reagent change |
| UTAG — fork `smith6jt-cop/utag` | tissue domains, celltype Step 4; dormant Xenium benchmark | `pyproject.toml [utag]`; `external/utag @ 056cb87ee90aa850450b806b075033bc8c50703a` → `utag_adapter.py:41` | islet domain stable across resolution; cohort pass ≈ ½ day |
| pSAM — fork `smith6jt-cop/pSAM` (unmodified) | the cosine-scoring idea behind celltype Step 3 (cell-level port; not installed) | `scripts/celltype/score.py:1-3` | encoding deviation measured: signed negatives + shared unit weight resolve 37/37 references |
| KINTSUGI — own repo (STAR Protocols 6, 103976, 2025; Zenodo 10.5281/zenodo.14984518) | the SOM lane's source (analysis end only) | `scripts/kintsugi/README.md`; `external/KINTSUGI @ 716d1ebd87db916ab04910c6374e04d2d9ea508e` (v1.2.3) → `som_adapter.py:55` | pyFlowSOM is not a declared KINTSUGI dependency; the documented conda env failed to solve here |
| pyFlowSOM 0.1.16, stackview | SOM and its frame rendering, kintsugi Step 2 | `external/.venv-som` (numpy 1.26.4, per `external/KINTSUGI/constraints.txt`) | requires numpy < 2, hence the separate venv |
| KODAMA 3.3 — fork `smith6jt-cop/KODAMA` (unmodified fork of `tkcaccia/KODAMAlegacy`) + Rnanoflann | embedding in `scripts/xenium_celltype`; a method in the dormant benchmark | R 4.5.1 user library; HPG pin `d2dbe30e` in `scripts/xenium_celltype/hpg/setup_envs.sh` (quarantined 2026-09-23 under data/ARCHIVE_2026-09-15/quarantine/); **unpinned** in `analysis/xenium/hpc/setup_envs.sh` | `landmarks` is the principled subsample; kNN capped at 500; M-parallelism does not scale (unresolved) |
| VALIS (valis-wsi 1.2.0) — fork unmodified | seed transform, restain lane | `~/.venvs/valis`; `scripts/restain/README.md` | two internal defects worked around (`MicroRigidRegistrar` off by default; `register_micro` reference); seed only — ICP does the rest |
| scvi-tools 1.4.0.post1 / scANVI, scanpy 1.11.5, squidpy 1.6.5, spatialdata 0.5.0, anndata 0.12.4 | the Xenium typer and spatial statistics | conda `scvi-env`; `analysis/xenium/README.md` | `--dl-workers 0` mandatory; rapids_singlecell absent and not needed |
| Nimbus-Inference — fork `smith6jt-cop/Nimbus-Inference` @ `241281ef` (unmodified) | per-marker positivity, broadtype arm E | `external/.venv-nimbus` (py 3.11, numpy 1.26, torch 2.5.1 cu121); `scripts/broadtype/arms/nimbus_rules.py` | public weights (JLrumberger/Nimbus-Inference); magnification declared 20 |
| KRONOS — fork `smith6jt-cop/KRONOS` @ `26f1bc02` (unmodified) | KRONOS v1 cell embeddings, broadtype arm F `kronos-clusters` | `external/.venv-kronos` (torch 2.6.0 cu124, xformers 0.0.29.post3) | **MahmoodLab/KRONOS access still pending review (403) on 2026-09-06**; the arm records `Unavailable(huggingface_access_pending)`. Superseded by the KRONOS2 row: the maintainer's approval is for KRONOS2 (2026-09-06) |
| KRONOS2 — Hugging Face `MahmoodLab/KRONOS2` snapshot `e2a54732` (marker-aware DINOv2 ViT-B/16, 768-d CLS, 256-px reference patches, 268-marker vocabulary; `kronos2_vitb16_teacher.pth` sha256 `afa261ef…`, `marker_metadata.csv` sha256 `1e679144…`), used **through CORAL** — fork `smith6jt-cop/CORAL` @ `2aa29e6d` (2026-08-09, unmodified) | cell embeddings, broadtype arm `kronos2-coral` | assets `external/kronos_assets/KRONOS2/` (gitignored; CC-BY-NC-ND 4.0, do not redistribute); checkout `external/CORAL/`; env `external/.venv-kronos` (transformers 4.56.0 / timm 1.0.19 per the KRONOS2 pins) | maintainer instruction 2026-09-06: use KRONOS2, and use it through the CORAL GitHub |
| TACIT v1.0.0 @ `206ae19c` (installed R package) | signature thresholding, broadtype arm C | R 4.5.1 user library; `scripts/broadtype/arms/tacit_run.R` | installed version uses Louvain and full RunUMAP; real signature `TACIT(data_expression, r, p, Signature)` |
| CELESTA `plevritis-lab/CELESTA` @ `e5e92ffc` + spdep 1.4.2 | spatial EM typing, broadtype arm D | R 4.5.1 user library; `scripts/broadtype/arms/celesta_run.R` | equal-proportions GMM prior inside the method, flagged on its row; default thresholds, no tuning |
| REDSEA — in-house `phenocycler/redsea.py`; reference `literature/external/redseapy` | stage 3 of the pipeline | `config.ini` (α = 1, subtract-only, `use_gpu = false`) | one recorded divergence from upstream (isolated-cell rim) |
| DuckDB, pyarrow / parquet | ingest streaming; every intermediate store; the app's data layer | `pyproject.toml` (unpinned), `apps/islet/R/00_globals.R` | — |
| opencv-headless 5.0, geopandas 1.1.4, shapely 2.1.2, scikit-image 0.26, zarr 3.3, spatialdata 0.8 | integration and restain geometry | `pyproject.toml [integration]` (unpinned); `.venv` | SimpleITK and spatialdata-io are declared but absent from `.venv` |
| TissUUmaps — fork `smith6jt-cop/TissUUmaps` (**+15 commits ahead**) 3.2.1.14; pyvips / libvips 8.15.1 | optional Viewer backend; `.tmap` builder | `deploy/build_tissuumaps_image.sh` (`TISSUUMAPS_REF`); `docs/tissuumaps_evaluation.md` | verified on synthetic slides only; never use the fork's `/slide` route |
| Avivator / viv (vendored static build); bioformats2raw 0.12.1, raw2ometiff 0.10.0 | default Viewer; qptiff → OME-TIFF | `apps/islet/www/avivator`; `scripts/tools/` | Avivator needs HTTP Range serving |
| QuPath 0.7.0 + InstanSeg `fluorescence_nuclei_and_cells-0.1.1` | every upstream measurement and annotation | `core/phenocycler/cohort-spec.json` | annotations are audit context / registration geometry only |
| Globus Connect Personal 3.2.8 | HiPerGator ↔ `/data/xenium` transfers | `~/globusconnectpersonal-3.2.8`; runs under `.venv`'s python | the only working HPG route |
| Shiny, DBI/duckdb 1.5.2, mirai 2.7.1, sf, reticulate 1.44, testthat; edgeR 4.6.3 (hidden tab); rdeck 0.5.2.9000 (off); vitessceR (optional) | `apps/islet` | R 4.5.1 user library; no lockfile; `deploy/environment.yml` omits edgeR | — |
| core/panels = XeniumPanelExplorer @ `076d1192` (Seurat 5.5 app) | the panel-taxonomy authority | subtree; snapshot stamped `git rev-parse HEAD:core/panels/data` (`0348d149…`) | renv deactivated (`renv.lock.r46-backup`) |
| Harmony (python), scVI (app side) | the trajectory build only | `scripts/senior/compute_trajectory_senior.py`, `train_scvi_senior.py` | input-starved; the deployed pseudotime predates the fix |

> **Correction 2026-09-18 — six "In use" rows are no longer in use.** Nimbus-Inference, KRONOS,
> KRONOS2/CORAL, TACIT, CELESTA and the KINTSUGI/pyFlowSOM SOM entered the repository as methods of
> the **broad-typing producer comparison**, whose producer choice was withdrawn 2026-09-08 and which
> is now on the root `CLAUDE.md` never-revive list (see the row added below). Their pins, checkouts,
> venvs and recorded findings stay as written — they are the record of what was measured and with
> what — but no lane in the current workflow calls them. The `arm` identifiers in the role column
> are the on-disk module names (`scripts/broadtype/arms/<name>.py`) and are left alone; in prose,
> say *method* or *the comparison*.

## Evaluated and closed — record kept, not revived

| package | verdict | what survives |
|---|---|---|
| DeepCellTypes — `vanvalenlab/deepcell-types @ 59d74a97…`, `deepcell-auth @ 523a4dd5…`, checkpoints `2026-06-15` and `2026-06-23-ptft` (md5 verified) | lane FAIL 2026-08-24, 0 of 5 held-out donors; excellent epithelial caller; scoped to this model at this checkpoint (`core/phenocycler/docs/DEEPCELL_TYPES_REVIEW_2026-08-24.md`) | `pyproject.toml [dct]`; `gates.py`, `tiff_window.py`, `window_mask.py`, the `label_space.py` idiom; the pre-registration discipline |
| CuPy REDSEA GPU backend (`cupy-cuda13x` installed; the extra declares `cuda12x`) | two showstoppers (`core/phenocycler/docs/REDSEA_GPU_REVIEW_2026-08-27.md`); defect 1 reproduced live 2026-09-01 | off (`use_gpu = false`) |
| MaxFuse (fork **+29, diverged**) | refused as a dependency; algorithm reimplemented in `integration/crossmodal.py`, itself dead as shipped | — |
| SimpleITK 2.5.6, itk-elastix 0.19.2 | present in `~/.venvs/valis` / conda `image2image` only; deliberately not run (the point-set fit beats intensity refinement 1.5–2.9×) | reversible judgement |
| rapids_singlecell | absent locally, not needed (CPU branches everywhere) | HPG only |
| RESTORE (own copy `smith6jt-cop/RESTORE`; `literature/external/RESTORE`) | adopted 06-19, retired 07-23, replaced 07-30, reinstated 08-10, lane retired 08-14 | history only |
| ordered-residual / simultaneous typing, Otsu floors, GMM control purification, the gold-standard lane, Gate 2–5 | measured rejections (`core/phenocycler/docs/METHOD_HISTORY.md`) | — |
| the broad-typing producer comparison (`scripts/broadtype/`, arms A–I incl. Nimbus, TACIT, CELESTA, KRONOS2/CORAL, the SOM and deepcell-types re-runs) | **producer choice withdrawn 2026-09-08** on twelve verified defects (`docs/BROADTYPE_PROTOCOL_2026-09.md` §4.9); on the never-revive list since. `data/development/broadtype/release/` is on disk and is **not** a verified typing | the QC figures P01–P05 / X01–X03, `verdict.py`'s recomputation from records, `gates.require_evidence`, the tiering rule |
| scimap (fork; conda env), Astir (`astir-env`), Vitessce, Streamlit, Lamian (`compute_lamian_pseudotime` never called) | superseded or archived | the scimap env and its 25-phenotype vocabulary linger in the app |

## Dormant here

BANKSY (`Banksy_py` fork), MISTy (`mistyR` fork; absent locally), the benchmark's KODAMA arm, the
`bench_py` conda env, `analysis/xenium/hpc/`, `scripts/xenium_celltype/hpg/`,
`scripts/slurm/run_full_pipeline.sh` (does not pass `--pipelined`), `scripts/senior/*`, the `states` stage.

## Own repositories

| repo | relevance | state |
|---|---|---|
| SpleenFollicleCounterQP | Phase 3b absorbed 2026-09-03 at `analysis/follicle`: H&E and DAPI U-Net lanes (root CI), notebooks H1–H20, H16/H20 retained as historical evidence | source clean on integration commit `1ed38d2` over upstream `6ae30ba`; former 15 dirty entries preserved at `f6cbcd4`. Four required files are absent from the live external `Measurements/` path and fail strict preflight |
| LymphVizkit | Phase 5 domain layer absorbed at `apps/R/ieCore` + `apps/lymph`; no app fork copied | verified bundle source `820bb38`; exact YAML and synthetic-fixture hashes in `apps/R/ieCore/IMPORT_PROVENANCE.md` |
| Islet-Explorer (Junior) | the publication build; deliberately separate; planned 2026-08-16 harvests accounted for | publication `main` at `f2a5c54`; portable validation, redaction, production-mode, and accessibility work from local `e589961` was adapted in D.4a; its Junior-specific OCI release pipeline waits for D.3/D.7 |
| Phenocycler_Analysis | subtree source | published `main` `fc06c0f`; its two unique content commits (`33e4da8`, `ffba0a8`) are ported and tested here; the other two commits are merges |
| Xenium_Analysis | Phase 3 source for `analysis/xenium` | clean published `main` `8062e554`; D.4a accounts for every unchanged, adapted, moved, and deliberately omitted blob |
| XeniumPanelExplorer | vendored source for `core/panels` | published `main` `076d1192`; 282/285 blobs are identical, two nested workflows are omitted in favor of the root path/taxonomy gate, and the README is the sole adaptation; the R suite remains a listed CI gap and passes locally with two explicit environment skips |
| isletscope, spatialbarrier, ADIIDAS, Panc_pLN_Analysis, HPCsizer, KINTSUGI-docs, OMICS, Skills_Registry | prior art or tooling | unreferenced by the repo, plans and memory |

About forty further forks under `smith6jt-cop` are unmodified mirrors with no project reference
(KRONOS is named once in a plan as an evaluation candidate); Celloscope (+4) and samapi (+8) carry
local commits nobody has described.

The exact D.4a audit, URLs, tips, and dispositions are in
`SOURCE_RECONCILIATION_2026-09-04.md`. All six archive-list repositories remain unarchived until
the A.4b identity release makes D.3 and the replacement proof possible.
