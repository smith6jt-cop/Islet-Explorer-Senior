# Workflow map — what runs today, where, and against what

> **Superseded note (2026-09-16):** the content-addressed run-framework sections below (`data/runs/<id>/`,
> run manifests, `adopt`/`resume`/`recover`) describe the workflow as it stood on 2026-09-02 and are no
> longer current — that framework was removed on 2026-09-16. See
> `core/phenocycler/docs/INTEGRATION.md` for the current state. The body of this file is left
> unchanged as a record of what was true when written.

**Status:** written 2026-09-02 against HEAD `1a23105` (branch `feat/dct-step0-pinning`, 56 commits
ahead of `main`, working tree clean). Every status cites a file, a command output or a commit. This
map links to the document that owns each step; it does not replace those documents. The planned
workflow that positions the auxiliary packages is `PLAN_2026-09-02_AUXILIARY_PACKAGES.md`; the
package index is `EXTERNAL_PACKAGES.md`.

Status vocabulary: **runs** (executed, output on disk) · **runs with caveat** · **blocked** (a guard
refuses it) · **dormant** (present, unrun here) · **retired** / **closed** (kept as record only).

## 1. Upstream PhenoCycler workflow — `core/phenocycler/phenocycler/`

Authority: `core/phenocycler/docs/WORKFLOW.md`, `DATA_README.md`, root `CLAUDE.md`.

| stage | reads → writes | command / env | status |
|---|---|---|---|
| 0 cohort manifest | `cohort-spec.json` (20 donors; names `Cellmeasurements_6539parentFixed.csv`) → `data/qupath_manifest.json` | `python -m phenocycler.pipeline manifest create`; repo `.venv` | runs; all 60 recorded source paths resolve |
| 1 `ingest` → 2 `geometry` → 3 `redsea` → 4 `expression` | measurement CSVs, qptiffs, GeoJSONs → `data/runs/<id>/{cells,geometry_qc,redsea,expression}` | `python -m phenocycler.pipeline run --config config.ini --pipelined --through expression`; `.venv`; REDSEA on CPU (`config.ini` `use_gpu = false`; `docs/REDSEA_GPU_REVIEW_2026-08-27.md`) | complete in run **`30f6689e4d27d295f890`** (76 GB; `data/runs/LATEST`; 31,454,041 cells; built 2026-08-28 from `0f5ea04`). The run **carries donor 6539's repaired parent assignment** (core 7,630 / peri 23,455 cells over 569 rings). `expression` is regenerable since `a30a517` |
| `controls` / `calibrate` / `type` / QuPath export | — | — | **blocked `NO_VERIFIED_TYPING`**: `marker_registry.py:323-327` requires the development status, so the gate is latched; lifting it needs a new reviewed registry and roster (root `CLAUDE.md`) |
| `states` (Ki67/PCNA, pair-independent) | expression → `cell_states` | `--through states`; not latched | never run on this cohort; needs a run under current code |
| `status`, `adopt`, `resume`, `recover` | — | `.venv`, from `core/phenocycler/` | `status` computes `run_id 67aec2f7a43de7600662` and reports the four stages **MISSING**: commit `91f5b87` moved the registry fingerprint, so the run-management CLI cannot see the on-disk run. `adopt plan --from-run 30f6689e…` would adopt **zero** stages (ingest/redsea/expression code hashes changed; no ingest transition in `adoption_compatibility.json`). The lanes below are unaffected: they pin the run by content hash (`scripts/celltype/panel.py:45-53`, asserted by tests) |

Retired, not fallbacks (`docs/METHOD_HISTORY.md`; commit `11d80c3`): RESTORE pair webs,
ordered-residual / simultaneous / hierarchical / probabilistic typing, marker-debug and oracle lanes,
Otsu floors, Gate 2–5, the gold-standard annotation lane.

Cohort facts every downstream step inherits are owned elsewhere: the batch alias (memory
`gate5-batch-confound`; `data/donor_acquisition_batch.csv`), donor 6482's 25 µm rings and the 6539
repair (`scripts/kintsugi/README.md` § Step 0), the exclusions (`phenocycler/cohort.py`).

## 2. Development lanes — `core/phenocycler/scripts/*`

All five read run `30f6689e…`, write only the repo-root `data/development/<lane>/`, write no
`StageManifest`, and cannot re-key `run_id` (only `phenocycler/` and `pyproject.toml` are hashed).
**No lane test suite is gated in CI** (`pyproject.toml` `testpaths = ["tests/current"]`). Collected in
`.venv`: celltype 53, kintsugi 22 (28 with `external/.venv-som`), dct 122, xenium_celltype 16 (+4 in
`hpg/`), restain 78.

### 2.1 `celltype/` — the cosine / rule-table lane (`scripts/celltype/README.md`)

| step | what | on disk | status |
|---|---|---|---|
| 1 `run_step1.py`, `figures.py` | per (donor, marker) tail co-occurrence vs derived controls, Moran's I, acquisition-group flag; fail-closed verdict | 20 donors, 2026-08-28 | runs (~10 min); all 1,120 cells `usable` |
| 2 `run_step2.py`, `verify_step2.py` | **UniFORM** background-mode rigid shift | 20 donors, 2026-08-28 | runs; 44/56 under the bar, 12 shape-limited, CD57 refused as a reagent change; `verify_step2` prints only, writes no artifact |
| 3 `run_step3.py` | cosine vs the repo-root `phenotype_workflow.csv` (37 references; pSAM-style scoring at cell level) behind `require_min` | **donor 6450 only; older than the Step 2 regeneration** | runs with caveat; the proliferating leaf is a PCNA background readout (README § 2026-08-28 regeneration) — a maintainer declaration, not a tuning |
| 4 `run_step4.py` | **UTAG** tissue domains, max_dist 20 µm | **6450 only; stale** | runs with caveat; cohort ≈ ½ day |
| `run_islet.py` | `Islet_N` endocrine prior in the core; every reference competes in the ring | 20 donors, 2026-08-28 (1,351,845 islet-associated cells; 6539 carries 18,519 peri) | runs |

### 2.2 `kintsugi/` — the SOM lane, an independent second method against the same rule table (`scripts/kintsugi/README.md`)

| step | what | on disk | status |
|---|---|---|---|
| 0 `run_step0.py` | column-selective extraction of the QuPath signed islet / tissue distances, validated against `cell_region` | 20 donors | runs; 0 unmatched `object_id` |
| 1 `features.py` | feature block on UniFORM-aligned values; arms A/B/C; scaler proven to be feature weighting | library (arm A only ever run) | — |
| 2 `run_step2.py` | **pyFlowSOM** SOM → Leiden meta-clusters | donor 6551, arm A, nodes 9–100 | runs with caveat; needs `external/.venv-som` (numpy < 2) |
| 3 `run_step3.py` | annotate by hierarchy behind the background-mode gate and a matched null | 6551, nodes 9–49 | runs with caveat |
| 4 `run_step4.py` | concordance with 2.1 on `(donor_id, object_id)` | 6551 | 27.6 % label / 62.5 % branch agreement at 25 nodes; the largest disagreement is endocrine resolution |
| 5 | donor → islet → cell hierarchical bootstrap | — | **planned, not built** |

### 2.3 `dct/` — DeepCellTypes pilot (`scripts/dct/README.md`; `docs/DEEPCELL_TYPES_REVIEW_2026-08-24.md`)

**Closed 2026-08-24 at the pre-registered bar: 0 of 5 held-out donors pass** (40 runs, 129,475 cells;
epithelial anchors 0.94–0.97 AUROC, 0 of 13 immune anchors). Scoped to this model at this checkpoint,
never to the panel (standing note of 2026-08-25). Two symbols are live cross-lane dependencies of the
closed lane: `gates.require_evidence` (imported by `celltype/usability.py:38`, never called) and
`tiff_window.read_window` (used by `restain/tiles.py`). The review document's "git revert
`scripts/dct/`" instruction would break both.

### 2.4 `xenium_celltype/` — KODAMA + null-gated group annotation (no README; `docs/HANDOFF_2026-08-26_XENIUM_CELLTYPE.md`)

`xenio.py` loader → PCA → **KODAMA 3.3** (R 4.5.1, `run_kodama.R`) → KMeans → `nullgate.py` (a
lineage is awarded to a *group* only where its identity core beats a size- and detection-matched
random null, as a contrast; negative contrasts refused) → composition with CIs. Run on **donor 6539
only**, a 50,000-cell subsample, k = 20 and 30. `build_phenotyped_base.py --donors …` builds the
established Xenium typer's input (run for 6539 only). `hpg/` is a complete SLURM shard lane that has
never executed (`ssh hpg.rc.ufl.edu` → `Permission denied`). The handoff records: do not report these
compositions; KODAMA's M-parallelism does not scale; `--max-cells` double-subsamples.

### 2.5 `restain/` — post-Xenium IF ↔ Xenium, same section (`scripts/restain/README.md`)

| step | what | on disk | status |
|---|---|---|---|
| `export_dapi.py` → `valis_register.py` (**VALIS 1.2.0**, `~/.venvs/valis`) → `nuclei.py` → `run_refine.py` (`register.icp` affine on the two nucleus fields) → `run_align.py` gate | nucleus-by-nucleus registration | 6436, 6442 | runs: **0.546 / 0.630 µm** median per cell, PASS on 243/243 and 365/365 windows, spatial hold-out within 0.07 µm |
| `restain_registration_inputs.csv`, `registration_inputs.py` | registration-only qptiff ↔ Xenium provenance; qptiff XML supplies calibration and literal DAPI page; no antibody interpretation | 20 recorded / 18 eligible | all 20 local qptiffs reconcile to `xenium_paths.csv`; 6457/6579 are filtered centrally; 6436/6442 pass required-file preflight; 155 synthetic re-stain tests pass |
| `channel_check.py` | islet enrichment per Opal channel in the IF frame, no registration | 5 recorded sections | runs: 6436 Opal 520 (GCG) 66× vs Opal 690 (INS) 2.2×; 6442 401× / 17.1×; Opal 570 is FAILED on the three nPOD sections but not on 115 / 6591 |
| `run_hormone.py --transform refined` (`9e88272`) | AUROC of IF intensity vs the within-cell RNA surrogate contrast, bar 0.75 | `hormone_panc_refined.json` | 6436 INS 0.514 / GCG 0.522; 6442 INS 0.657 / GCG 0.728 — **bar not met**; extraction is still an 8 µm disc; `main` returns 0, no exit codes |
| `protein_rna.py`, `protein_rna_figures.py` (`2dade20`, `1a23105`) | **label-free** rank association (Spearman ρ of IF intensity vs transcript count) against the 18-gene pancreas endocrine identity core from `core/panels`; the AUROC against an invented positive class was removed | both donors | runs; "defining a positive is the maintainer's call" — nothing here labels a cell. 6436 is T1D with focal insulin, so any section-wide statistic on its INS channel says little |
| `run_explorer_export.py` (`9e88272`) | Xenium Explorer bundle (pyramidal OME-TIFF, alignment CSV, GeoJSON); deciles only, no INS+/GCG+ calls; a gated mixture branch exists but does not run | both donors | runs |
| `run_register.py` (tissue-outline ECC), `validate.py` | superseded transform; islet-vs-RNA enrichment | — | retired / not re-run; earlier numbers measured at 31–42 µm are void |

Only 2 of 20 donors have all three inputs (Xenium `morphology_focus`, QuPath Tissue/Islet GeoJSON,
a `restain_sections.csv` row). Recorded rows: 5; GeoJSON exports: 2; IF qptiffs local: 20.

### 2.6 `restain/pheno_annotations.py` — PhenoCycler QuPath annotations ↔ Xenium, **serial** section

A different section pair from 2.5 and a different transform: this one rides the S4 registration in
§4, not the nucleus-matched one above. It shares the Explorer 4.1 GeoJSON contract and nothing else
about what the result may claim.

| step | what | on disk | status |
|---|---|---|---|
| `pheno_residual_search.py` | joint rotation × translation correspondence landscape, translation free at every angle; `peak / p95` against the donor's own landscape | `pheno_residual_search_panc.json`, `figures/pheno_residual_search_*.png` | runs (~2 s/donor, no image I/O). Both donors now peak at **0 deg** within one 50 µm bin, at 6.0× and 6.9× p95. This is what found 6436's 800 µm error; the rotation-only sweep it replaces could not, because rotating about the warped moving centroid fixes that centroid and so cannot express a translation |
| `run_pheno_annotation_export.py --donors 6436 6442` | six QuPath classes onto the Xenium morphology grid, in **level-0 pixels** at 0.2125 µm/px | `explorer_bundle/donor_*/annotations/*_pheno_*.geojson` | runs. **85,026** polygons / 13 files (6436) and **118,645** / 16 files (6442), each validated against the Explorer 4.1 contract on reopen and published under a rollback transaction |
| `PA.xenium_tissue_outline` (new 2026-09-10) | the region annotations are selected against: `rasterize.tissue_mask` on Xenium cell centroids → `cv2.findContours` | `pheno_annotations_qc.json → tissue_outline` | runs. **77.9** and **102.7 mm²**, which is **72.3%** of the morphology rectangle on *both* donors — the exporter previously tested against that rectangle, so more than a quarter of the area an annotation could pass in was blank background. Removed a further 24,867 (6436) and 24,434 (6442) parts that had been landing on it |
| selection rule | keep only polygons lying **entirely** inside the outline; drop straddlers whole and count them | per-class in the QC file | clipping is refused by design — it would invent an anatomical boundary and misstate the annotation's area |
| `pheno_annotation_figures.py`, `pheno_transform_diff_figure.py`, `pheno_capture_box_figure.py`, `pheno_registration_inputs_figure.py` | verification sheets drawn from the **published** files, so a figure cannot agree with the export by construction | `figures/pheno_*` | all regenerated 2026-09-10 against the corrected transforms. Three had been rendering the superseded −9.19 deg transform, and `pheno_capture_box_figure.py` defaulted to drawing a −7 deg candidate that a settled question had left behind |

The bundles remain under the hold in `explorer_bundle/DO_NOT_SEND.md`. Nothing here has been
visually confirmed in Xenium Explorer yet, which is the check that decides it.

## 3. Xenium lane — `analysis/xenium/` and `core/panels/`

Authority: `analysis/xenium/README.md`, `analysis/xenium/CLAUDE.md` (decisions #1–#17). Environment:
conda `scvi-env` (scvi-tools 1.4.0.post1, torch 2.9.0+cu128; verified 2026-09-01), never `xenium_analysis`.

| step | status |
|---|---|
| `core/phenocycler/scripts/xenium_celltype/build_phenotyped_base.py --donors …` → `/data/xenium/processed/<serial>/<roi>/<serial>_phenotyped.h5ad` | runs; 6436 (`0059907`: 641,487 retained) and 6442 (`0059947`: 906,211 retained) rebuilt 2026-09-03; the older 6539 development run also exists |
| `scripts/finalize_panel.py --paths …` — stamps `identity_panel` (478 genes) from the `core/panels` snapshot | runs; snapshot fresh (`SOURCE.txt` tree `0348d149…` = `HEAD:core/panels/data`) |
| `scripts/hierarchical_phenotype.py` — two-tier identity-core → **scANVI** (`--dl-workers 0`) | runs with explicit nested `--paths` and complete-set preflight. Slim sidecars written for 6436 (416,962 anchors / 224,525 propagated; 1,695 low-signal; 25.47% confidence <0.50) and 6442 (755,067 / 151,144; 1,435; 22.08%). IDs align exactly and required fields have no nulls; this is development evidence, not a declared integration release |
| `scripts/insulitis_analysis.py` | nested `--paths` seam fixed; resolves each immune sidecar beside its base file and fails preflight while those immune inputs are absent |
| `scripts/fetch_subpanels.py` | runs; regenerates the gitignored snapshot from `core/panels` |
| `scripts/benchmark/` (BANKSY / UTAG / KODAMA / MISTy), `hpc/`, the thymus lane | dormant by design; `bench_py` and the HPG R library do not exist here |
| `core/panels` Shiny app (Seurat 5.5) | loads under R 4.5.1; renv deactivated (`renv.lock.r46-backup`); 50 testthat blocks unrun |

## 4. Integration — `core/phenocycler/phenocycler/integration/` (S0–S11)

Authority: `core/phenocycler/docs/INTEGRATION.md` (its preamble and paths are corrected in this change).

| stage | status |
|---|---|
| `manifest` S0 (`--report`, `--verify-paths`) | runs. Rebuilt 2026-09-03: 57 rows, **20 paired** pancreas donors, 35 Xenium-only sections and 2 `donor_unknown`; all 20 current-run PhenoCycler qptiff/GeoJSON pairs resolve exactly through `qupath_manifest.json` plus the tracked two-row override file. The missing 6591 pLN source and two unknown Xenium runs remain explicit warnings |
| `export_pheno` S1a | geometry half runs from immutable run `30f6689e4d27d295f890`: **20 `roi=panc` partitions / 25,150,537 analysis-eligible cells**, schema v1 geometry profile, 57 identity-safe features, zero 6457/6579. INS/GCG/SST are withheld; `islet_num` is carried only as `qupath_islet_num`; identity export still requires a validated release |
| `import_xenium` S1b | geometry half runs from per-section bundles. Verified on 6436 (695,113 cells) and 6442 (958,744): zero duplicate cell IDs, blank identity and nullable evidence. Embedded label-shaped columns are ignored; a configured Xenium release must provide an exact cell-ID call set |
| `register`, `qc`, `figures` | geometry stages, no longer identity-blocked. **Rebuilt 2026-09-10**: both donors now `PASS`. 6436's transform was wrong by 12.5 deg and 2.6 mm — its correct coarse hypothesis ranked 4th of 144 while `register_image` refined only the top 3 — and correcting it takes islet RMSE 176.9 → **89.1 µm** and Xenium islets matched within 100 µm from 12 of 136 (chance ≈ 5) to **69**. 6442 is the regression control and barely moves (+5.11 → +5.09 deg; 107.9 → 92.1 µm from the masked ECC alone). Changes: adaptive coarse plateau (`reg_coarse_keep_frac`) instead of a fixed rank cut; ECC masked to the Xenium footprint (`refine_ecc(moving_mask=)`); `reg_nonrigid = false`, on measurement — the field failed the exporter's round trip on 15.0% / 6.1% of real annotation vertices and 6436's earlier field was *folded* (min Jacobian det −1.657 over 0.85%). **The 0.80 tissue-Dice gate is gone**: it is unreachable by extent (ceiling 0.565 / 0.703) and, worse, it moves *against* correctness — 6436's Dice falls 0.549 → 0.530 at the right answer. The gate is now `footprint_ncc_z` (DAPI correlation inside the warped footprint against 120 other placements of that same footprint on the same section: **6.9** and **6.7**, vs 2.3 for 6436's old transform) with the joint rotation+translation islet landscape as corroboration (both peak at 35 µm, ~6× the landscape p95). New diagnostic: `scripts/restain/pheno_residual_search.py` |
| `vocab structures match cellmatch grid donor crossmodal sameslide postxen` | identity stages; blocked unless the required PhenoCycler and/or Xenium releases pass `phenospace.py`. Existing cached structure files are ignored. `hormone_seed_mask` rejects a geometry table, so QuPath context and raw hormone intensities cannot become identity by accident |
| `boundaries`, `rasterize`, `contract`, `tissues`, `sections` | libraries, live (`boundaries.py` and `rasterize.py` extended 2026-09-01 for the restain lane) |
| `crossmodal` S9 | identity-gated and structurally input-starved as shipped (`import_xenium` writes only `feat__score_*`, which `resolve_anchors` excludes as direct gene anchors); no inferred links were produced |

## 5. Application lane — `apps/islet/`, `scripts/senior/`, `deploy/`

Authority: `apps/islet/CLAUDE.md`, `scripts/CLAUDE.md`, `deploy/DEPLOY.md`.

- **Shared R layer (D.2 plus D.4a complete 2026-09-04)**: `apps/R/ieCore` is the one implementation of
  domain loading, typography, safe joins, and statistical helpers. The islet app reads
  `apps/islet/config/islet_default.yml`; `apps/lymph` contains the exact LymphVizkit YAML plus a
  588 KB synthetic contract fixture, not another app fork. Portable hardening from Junior
  `e589961` now validates browser-controlled choices and paths, redacts failures, disables debug/AI
  credential access in production, sanitizes Shiny errors, and adds keyboard/accessibility fixes.
  The 188-assertion suite checks both domains, every tab's UI/server wiring, single helper ownership,
  and these security boundaries.
- **Runtime**: R 4.5.1 (`/opt/R/current` via `/usr/local/bin/Rscript`) with the `~/R/…/4.5` user library —
  the only library holding duckdb, mirai, rdeck, edgeR; `/usr/bin/Rscript` is R 4.6.1 without them.
  Reticulate python: `~/.local/share/islet-explorer-py`. Both pins live only in the gitignored
  `apps/islet/.Renviron`.
- **Data**: `data/app_data/` (1.9 GB) is the **legacy 15-donor Junior/CODEX cohort** (5,214 islets, 21
  phenotypes, built 2026-04-25, no `run_id`). Five of those donors are in the 20-donor run cohort and
  none is ND. The `scripts/senior/*` build chain is input-starved (no `data/singlecell_protein*.h5ad`,
  no `data/parquet/`); Spatial Card C's per-cell distance view is dead.
- **Viewer**: 15 of 22 OME-TIFFs converted (`/data/islet_ome_tiff`); 1 of 22 `.tmap` projects (6539);
  TissUUmaps verified on synthetic slides only.
- **Serving**: shiny-server and nginx are running but serve **no** Islet Explorer. `/srv/shiny-server/islet-explorer`
  and the nginx `local_images` alias point into `~/panc_CODEX`, which moved: Junior's publication
  `main` is `f2a5c54`; the clean local checkout is on review-hardening commit `e589961`. The
  `islet-explorer-merged` entry was never created.
  `deploy/` still carries `<project>` placeholders.
- **CI**: `islet-r.yml` (ieCore install + package/app parse + 188 testthat assertions), `phenocycler-ci.yml` (ruff; current/research/integration suites),
  `phenocycler-release.yml` (tag-triggered; no tests).

## 6. Consolidation status (`docs/CONSOLIDATION.md`, 2026-09-02 section)

Phase 0 done · 1 void · 2 done (no `phenocycler-v*` tag yet; rollback assets retained) · Phase 3 done
(`analysis/xenium` landed 2026-08-28; `analysis/follicle` absorbed and tested 2026-09-03) · 4 rewritten in writing, not executed ·
5 complete 2026-09-04 · 6 source reconciliation complete, archival gated. Phenocycler
`33e4da8`/`ffba0a8` and portable Junior `e589961` work are adapted; Panels/Xenium blob audits and all
six source tips are recorded in `SOURCE_RECONCILIATION_2026-09-04.md`. D.3, repository archival,
the first package tag, fresh-clone release proof, and any production action still wait on A.4b and
the maintainer declarations.

## 7. Environments that actually run things

| environment | runs | note |
|---|---|---|
| repo `.venv` (py 3.12.3, numpy 2.5.1, torch 2.13.0+cu130, CUDA ok) | the pipeline; integration; celltype; kintsugi 0/3/4; dct; xenium_celltype (python); restain except VALIS; both pytest suites | **Globus Connect Personal's daemon runs under this interpreter** — do not recreate it. SimpleITK 2.5.6 and spatialdata-io 0.7.1 were added without changing numpy/torch on 2026-09-03; `pip check` and the formerly skipped SimpleITK test pass |
| `~/.venvs/valis` (valis-wsi 1.2.0, SimpleITK 2.5.6) | `valis_register.py` | its pins are incompatible with `.venv` |
| `external/.venv-som` (numpy 1.26.4, pyFlowSOM 0.1.16) | kintsugi Step 2 | exists because `external/KINTSUGI/envs/env-linux.yml` failed to solve (2026-08-27) |
| conda `scvi-env` (py 3.13.9) | the Xenium typer | the only lane environment recorded in a tracked file |
| R 4.5.1 + `~/R/…/4.5` | `apps/islet`, `core/panels`, KODAMA, Seurat | recorded only in a gitignored file |
| conda `phenocycler-analysis`, `phenocycler_integration` | — | **do not exist on this workstation**; HPG-only recipes (`environment*.yml`) |
| conda `astir-env`, `CODApivot`, `Kintegrate`, `samapi`, `vesselvio`, `XeniumPC_Analysis`, `image2image`, `KINTSUGI` (dangling editable install), `scimap` | — | no operational reference in the repo |
| HiPerGator | `analysis/xenium/hpc/`, `xenium_celltype/hpg/`, `scripts/slurm/` | ssh denied from this workstation; Globus is the only route |

## 8. Read-only checks behind this map (2026-09-01 through 2026-09-04)

```bash
cd core/phenocycler
../../.venv/bin/python -m phenocycler.pipeline status --config config.ini        # run_id 67aec2f7…, 4× MISSING, 3× BLOCKED
../../.venv/bin/python -m phenocycler.integration.pipeline --status
../../.venv/bin/python -m phenocycler.integration.manifest --verify-paths        # 20/20 current-run Pheno qptiff/GeoJSON pairs exact; source catalog gaps remain explicit
../../.venv/bin/python -m pytest -q                    # 322 current tests
../../.venv/bin/python -m pytest -q tests/research    # 18 research tests
../../.venv/bin/python -m pytest -q tests/integration # 303 integration tests
../../.venv/bin/python -m pytest -q scripts/restain/tests # 155 re-stain tests
git -C ../../external/UniFORM rev-parse HEAD    # 7c33396d12e5a19dc0b36e46a94c3293de0cd814
git -C ../../external/utag rev-parse HEAD       # 056cb87ee90aa850450b806b075033bc8c50703a
git -C ../../external/KINTSUGI rev-parse HEAD   # 716d1ebd87db916ab04910c6374e04d2d9ea508e
```

## 9. 2026-09-05 addendum — what changed since this map was written

- New lanes: `core/phenocycler/scripts/broadtype/` (the broad-typing bake-off; protocol
  `docs/BROADTYPE_PROTOCOL_2026-09.md`), `scripts/xenium_typing/` (X1/X2/X3 Xenium producer comparison),
  `scripts/report/` (the two-donor figure set under `data/development/report/`), `scripts/common/` (Track E.3
  executed: gates/tiff_window/window_mask/label_space moved out of the closed dct lane).
- §4: `export_pheno` now exports the complete cell set and attaches a declared release; `import_xenium`
  writes per-gene type anchors and `feat__if_*`; verdict enforcement in match/cellmatch/grid;
  `config.ini` names both release directories. The stale `structures/phenocycler/*/islets.parquet` were removed.
- §3: all 20 Xenium bases exist and are finalized (identity_panel 478, zero state-only overlap); the
  integration manifest resolves `xenium_h5ad` for all 20 paired donors.
- §7: `external/.venv-nimbus`, `external/.venv-kronos` (KRONOS access pending), R CELESTA + spdep installed.
- Xenium Analyzer outputs (`analysis/`, `analysis_summary.html`, `metrics_summary.csv`, `experiment.xenium`) for
  the 16 bundles under `/orange/maigan/Xenium_T1D_PancPLN` landed 2026-09-06 (Globus task `793c7c34…`, 480 files,
  5.8 GB, into the matching bundle dirs under `/data/xenium/HuPrime5kPanTissue/`). The first attempt
  (`b84a08bc…`) had been aimed at endpoint `6ea86da5…` "D6-36_FBMS-J506-W1", which is a *different* workstation, and
  failed every retry with `500 Path not allowed`; this machine is `f1b69b9e…` "LInux Processing" and its Globus
  Connect Personal must be started (`~/globusconnectpersonal-3.2.8/globusconnectpersonal -start`) before a
  transfer. The MM046/MM047 runs (6450, 6505, 6534, 6623) were recorded here as "not found under either
  HiPerGator root" — they were not missing, they were under a **different project root**:
  `/orange/brusko/xenium/HuPrime5kPanTissue/Y9QRJ6/MM046_…` and `…/MM047_…`, not `/orange/maigan/…`
  (maintainer, 2026-09-10). Their DAPI landed 2026-09-10 (Globus task `ef4698c5…`, 4 files, 3.2 GB).
  Two things differ about those two batches and both matter. They use the **newer Xenium channel
  naming** — `morphology_focus/ch0000_dapi.ome.tif`, not `morphology_focus_0000.ome.tif`;
  `rasterize.read_xenium_morphology` already falls back to `*.ome.tif` sorted, so `ch0000_dapi` resolves
  as channel 0 with no renaming, and renaming it to match the old convention would be a lie about which
  file it is. And only the DAPI channel was fetched, so `tifffile` warns that three companion channels
  are missing and zeroes them — the same warning the other 18 donors already produce, and harmless
  because only channel 0 is read. `experiment.xenium` is still absent for these four; nothing in the
  registration path reads it.
