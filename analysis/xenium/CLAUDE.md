# `analysis/xenium` — repository guide

> **Absorbed 2026-08-28** from `smith6jt-cop/Xenium_Analysis@8062e554` under consolidation
> Phase 3. Operational instructions (env, paths, what was dropped and why) are in
> [`README.md`](README.md); this file is the durable architectural guide.
>
> **The numbered decisions #1–#17 below are carried across intact, with their numbering
> preserved.** Other documents in this repo and in `core/phenocycler/` cite them by number,
> so renumbering would break those references. Where a decision no longer describes this
> machine, a dated standing note is added *beside* it — the decision itself is not silently
> rewritten. The Tier A/B/C/D framing at the end does the pruning, not renumbering.

## Standing notes, 2026-08-28 (read with the decisions below)

- **#4 / #7 / #9 — `rapids_singlecell`.** Absent here and **not needed**: `try_import_rsc()`
  returns `None` with a warning and every call site has a complete CPU branch. No code change
  is required. #7's `LD_PRELOAD` recipe is HiPerGator-only and now lives in `hpc/`.
- **#8 — B200 Tensor Cores.** `setup_torch_tf32()` runs unchanged on the RTX 4000 Ada and TF32
  is enabled; the B200-specific *timing* claim is not applicable to this GPU.
- **#14 — the panel audit.** Now at `analysis/xenium/panel_audit.csv` (tracked). The old
  default pointed at `data/processed/panel_audit.csv`, a gitignored path that never held it.
- **#16 / #17 — the vendored snapshot.** No longer fetched over HTTP. It is **generated from
  `core/panels`**, which *is* XeniumPanelExplorer imported as a subtree, verified
  byte-identical to the pin `@076d1192` on 2026-08-28. The snapshot lands at
  `data/development/xenium_panels/subpanels/` and carries a `SOURCE.txt` stamped with
  `git rev-parse HEAD:core/panels/data` (tree `0348d149…`); `tests/` fails if it goes stale.
  The substance of the pin is unchanged.
- **Environment.** `scvi-env`, not `xenium_analysis` — the latter is broken here
  (`_spropack` ImportError). See README.
- **Two roots.** `ROOT` (code) and `DATA_ROOT` (h5ads, default `/data/xenium`) are separate
  and must stay separate; see README for why collapsing them puts 60 GB inside the repo.
- **Cohort.** The pipeline below is written for samples `0041323` / `0041326`. The monorepo
  cohort is the **20 nPOD donors** in `core/phenocycler/cohort-spec.json`, whose Xenium
  bundles are at `/data/xenium/HuPrime5kPanTissue/`. The *method* carries across; the sample
  ids in the tables below do not.

## Pipeline (current, GPU-accelerated post-2026-04-25)

GPU annotations: **[GPU rsc]** uses `rapids_singlecell` (cuML kNN /
PCA / scale / spatial_autocorr); **[GPU cuML]** uses `cuml.cluster` /
`cuml.neighbors` directly; **[GPU torch]** uses scVI on B200 with TF32
Tensor Cores enabled in cell 1; **[CPU]** is intentionally CPU.

```
zarr (data/raw/{sample}.zarr)
  │
  ├─► notebooks/01_preprocessing_v2.ipynb  ← OPTIONAL (stage 02 can
  │       ingest zarr directly), but rerun for clean fresh outputs.
  │       Cell 1:    imports + B200 TF32 setup [GPU torch].
  │       Cell 18:   counts layer → HVG (seurat_v3, metadata only) [CPU]
  │                  → normalize_total + log1p + lognorm layer [CPU, fast]
  │                  → anndata_to_GPU + rsc.pp.scale (on panel_for_embedding) [GPU rsc]
  │                  → rsc.pp.pca(n_comps=50) [GPU rsc] (50 computed for the
  │                    cell-19 scree DIAGNOSTIC; the PC count used downstream
  │                    is fixed at 30 — the cell-19 "elbow" never sets it)
  │                  → anndata_to_CPU(convert_all=True)
  │                  → rsc.pp.neighbors(n_neighbors=10, n_pcs=30, cosine) [GPU rsc]
  │                  → sc.tl.umap(min_dist=0.05, spread=2.5) [CPU, intentional]
  │                  → sc.tl.leiden(resolution=1.5, flavor='igraph') [CPU]
  │       → data/processed/{sample}/{sample}_preprocessed.h5ad
  │
  └─► notebooks/02_phenotyping.ipynb  ← PRIMARY
        Cell 1:    imports + B200 TF32 setup [GPU torch].
        Cell 3:    ingest zarr (or load preprocessed if present),
                   drop control probes, QC filter (min_counts=50,
                   min_genes=20, max_counts=98%-ile, min_cells=100),
                   counts layer preserved.
        Cell 5:    normalize_total + log1p + lognorm layer [CPU, fast].
                   Marker discovery only if leiden_1.5 already exists.
        Cell 10:   panel-aware pancreas marker dict (NO INS/GCG/SST/PPY).
                   LOW_CONFIDENCE_TYPES = {Alpha, Delta, Epsilon, Acinar}.
        Cell 13:   sc.tl.score_genes for each cell type [CPU, fast].
                   Plot deferred to cell 17 (so it lands on scVI UMAP).
        Cell 15:   argmax assignment → ANTI-ACINAR OVERRIDE
                   (per-cell mean of curated specific anti-acinar markers
                   in lognorm; reassign Acinar→alt if margin > 0)
                   → Indeterminate flagging for low-coverage types
                   (median - 1·MAD on celltype_confidence).
        Cell 17:   scVI training [GPU torch, TF32]
                   (n_latent=30, gene_likelihood='nb',
                    continuous_covariate_keys=['cell_area'], NO total_counts)
                   → rsc.pp.neighbors(use_rep='X_scvi', cosine) [GPU rsc]
                   → sc.tl.umap(min_dist=0.05, spread=2.5, random_state=0) [CPU]
                   → sc.tl.leiden(resolution=1.0, flavor='igraph') [CPU]
                   → umap_qc_overlay.png (Spearman ρ < 0.3 gate)
                   → scvi_latent_diagnostic.png (scree of PCA on X_scvi +
                     per-dim usage + library-size leak gate; diagnoses the
                     representation phenotyping ACTUALLY uses. Replaces the
                     old gene-space pca_variance_ratio.png, which measured the
                     z-scored noise floor of the ~5k panel, not the embedding.
                     Re-render: scripts/render_scvi_latent_diagnostic.py)
                   → rsc.pp.neighbors(use_rep='X_pca') [GPU rsc] +
                     PCA-fallback UMAP (X_umap_pca, leiden_pca, both CPU).
                   → celltype_scores_umap.png + predicted_celltypes.png
                     rendered HERE on scVI UMAP.
        Cell 20:   per-cluster consensus annotation
                   (manual_celltype = leiden_scvi → majority predicted_celltype).
        Cell 22:   final celltypes overview, save to {sample}_phenotyped.h5ad.
        │
        ├─► scripts/finalize_panel.py --tissue <t> (stamps identity_panel +
        │     subpanel_<L> = identity_core + is_state_only), THEN
        ├─► scripts/hierarchical_phenotype.py --output inplace --dl-workers 0
        │     (IDENTITY-CORE rewire 2026-06-08 session 2) Merges labels INTO
        │     {sample}_phenotyped.h5ad: per-lineage ABSOLUTE mean lognorm over
        │     identity_core (specificity-filtered) → clean anchors (abs-argmax +
        │     margin + raw-count floor) + 'Unknown' → scANVI on the filtered
        │     feature union PROPAGATES → expression-masked predict() =
        │     celltype_broad + broad_confidence (0–1) + broad_n_identity_genes;
        │     per-lineage detection-gated argmax (identity_core subtype_hint;
        │     MARKER_DICT thin-fallback) → celltype (Indeterminate-free); state
        │     programs → obs['state_*_score']. celltype_flat + X_scvi_broad kept.
        │     See decisions #16/#17. (No --reuse-scvi.)
        │
        └─► notebooks/03_spatial_analysis.ipynb
              Cell 1:    imports including cuKMeans, cuNearestNeighbors.
              Cell 2:    load + slim adata (drop X_pca, X_scvi, etc.).
              Cell 4:    sq.gr.spatial_neighbors(n_neighs=10, delaunay) [CPU, fast].
              Cell 6:    sq.gr.nhood_enrichment(n_perms=100, n_jobs=1) [CPU,
                         documented numba/joblib readonly-array workaround].
              Cell 8:    sq.gr.co_occurrence on 150k subsample [CPU].
              Cell 10:   rsc.gr.spatial_autocorr (Moran's I) [GPU rsc] —
                         attempts full 1.2M cells; 200k subsample retained
                         as fallback if OOM.
              Cell 12:   spatial niches: cuML NearestNeighbors (k=50, brute) [GPU cuML]
                         + cuML KMeans (n_clusters=12, scalable-k-means++) [GPU cuML]
                         + 2-pass numpy majority-vote smoothing [CPU, fast].
              Cell 14:   sq.gr.ligrec on 100k subsample, n_perms=50 [CPU].
              Cell 16:   sq.gr.ripley on 100k subsample [CPU].
              Cell 18:   defensive try/except CSV exports + h5ad save.
```

Stages 04, 05, 06 (group comparisons / cross-tissue / Phenocycler) are
present in `notebooks/` but not part of the per-sample default flow.

## Critical architectural decisions — do not regress

1. **Never pass `total_counts` (or any raw library-size proxy) as a scVI
   *continuous covariate*.** The NB likelihood already models library size via
   its size factor; passing it again double-regresses depth — that was the
   original library-size leak PaCMAP was a workaround for. *This* is the
   invariant. Which *other* nuisance covariates a given scVI model uses is a
   per-model design choice, not a fixed value: the flat full-panel model
   historically used `cell_area`; the per-sample hierarchical **broad** model
   may use within-sample technical covariates (e.g. `segmentation_method`,
   `nucleus_multi`, `z_level_cat`) or none. Choose covariates that vary within
   the sample and aren't collinear with the depth the size factor already
   captures. `cell_area` is cell-type-linked (large acinar vs small
   lymphocytes), so include it only where regressing it out won't erase the
   lineage signal you're trying to embed.
2. **Don't select embedding genes by HVG *variance ranking*.**
   `var['highly_variable']` (seurat_v3 top-2000) is QC/reference metadata only
   and must never subset the PCA/scVI input — ranking ~5k curated probes by
   variance just surfaces the panel's z-scored noise floor (and an old
   "HVG-after-PCA" path was dead code). *This* is the invariant. It does **not**
   forbid purposeful gene-set design: the curated `panel_for_embedding` (full
   panel minus the audited tissue-restricted + sex-chromosome genes, #13/#14) is
   the default, and a deliberately narrowed, lineage-discriminating panel is a
   legitimate evolution — e.g. the broad model's `embedding_panel` (fine-marker
   union + breadth-filtered discriminators, frozen by
   `scripts/select_embedding_features.py`), selected by between-cell-type
   discrimination + detection breadth, NOT by raw variance. The bright line is
   the *selection criterion* (curated / discriminative = OK; HVG-variance
   ranking = not OK), not the gene count. The `highly_variable_genes.png` figure
   carries this disclaimer in-figure.
3. **Squidpy graph stats are subsampled.** Full 1.2M-cell ops with
   permutation tests or all-pairs Dijkstra are infeasible.
   `centrality_scores` skipped entirely. (Moran's I via
   `rsc.gr.spatial_autocorr` runs full-1.2M on GPU now — see #7.)
4. **scanpy CPU UMAP, not rapids/cuML.** Rapids UMAP produced extreme
   outliers on this data. **Other ops (scale, PCA, neighbors, KMeans,
   Moran's I, score_genes) ARE on GPU via `rapids_singlecell` — only
   UMAP and leiden stay CPU.**
5. **`leiden` flavor=`igraph`, `n_iterations=-1`** — fast enough on
   CPU; we don't switch to cugraph leiden because `dask-cuda 24.12`
   has a hard incompat with the installed `dask 2026.1`.
6. **Pin random_state=0 everywhere** — PCA, neighbors, UMAP, leiden,
   scVI training. Reproducibility is non-negotiable.
7. **Driver-script env fix is required for GPU paths.**
   `scripts/env_rapids.sh` (sourced by `scripts/run_pipeline.sh`) exports
   `LD_PRELOAD` (conda `libstdc++.so.6` +
   `/apps/compilers/cuda/12.8.1/lib64/libcudart.so.12`) and adds
   `$CUDA_HOME/nvvm/lib64` to `LD_LIBRARY_PATH`. Without these,
   `cudf` 24.12 calls `pynvjitlink.patch_numba_linker()` and raises
   `RuntimeError` because `numba_cuda` 0.0.17.1 already provides those
   patches — every `import rapids_singlecell` fails and the notebooks
   silently fall back to CPU. The LD_PRELOAD pins runtime libcudart at
   12.8 (matching driver 12.8), so the MVC patch is skipped entirely.
8. **B200 Tensor Cores must be enabled for scVI.** Cell 1 of stages
   01/02/03 sets `torch.set_float32_matmul_precision('high')` and
   `torch.backends.cuda.matmul.allow_tf32 = True`. Without this scVI
   training takes ~30-50% longer on B200 (compute capability 10.0).
9. **`rapids_singlecell == 0.14.1`** — required for scanpy 1.12 API
   compat (the older 0.11.1 calls `_handle_mask_var` with the wrong
   signature for scanpy 1.12+).
10. **Upper QC filter is DENSITY-based (counts/cell_area), not absolute
    counts.** In FFPE Xenium, `cell_area` varies 5-10x and Spearman
    (total_counts, cell_area) ≈ 0.78 in this dataset, so an absolute-counts
    upper cap predominantly removes large *healthy* cells (acinar /
    ductal / large epithelial) — 88% of cells dropped by `counts > q98`
    have completely normal counts/area density. We use `counts/area > q98`
    to drop genuine outlier-density cells (likely doublets) instead.
    Lower bounds remain absolute (`min_counts=50`, `min_genes=20`)
    because those are sanity floors for noise/segmentation artifacts,
    not size-correlated. Implemented in stage 01 cell 14 and stage 02
    cell 3 (standalone branch).
11. **Lineage phenotyping (post-density-filter doublet refinement).**
    Implemented in stage 02 cell 22 (and as a post-hoc script
    `scripts/lineage_phenotype.py`). For each cell, count how many
    distinct markers from each canary lineage panel are detected at
    `lognorm ≥ log1p(3)` ("depth", sample-comparable). When ≥ 2
    mutually-exclusive groups are co-detected, the cell is resolved to
    the dominant lineage if `top_depth − second_depth ≥ 2 AND top_depth ≥ 2`;
    tagged `lineage_phenotyped` and `celltype_lineage` is set to the
    dominant lineage. ≥ 3 mutually-exclusive groups co-detected → tagged
    `doublet_suspected` (segmentation merge), kept as-is.
    `celltype` is NEVER overwritten by lineage refinement —
    `celltype_lineage` is the refined column for downstream consumers.
    Acinar canary panel was expanded from `[AMY1A]` to `[AMY1A, CUZD1]`
    so depth comparisons are fair across all lineages. The exact panels
    live in `scripts/verify_upper_tail_doublets.py::CANARY_PANELS`.

12. **Don't let a thin marker panel bias the argmax — avoid 1-gene "panels".**
    Under detection asymmetry, a marker/score argmax tilts toward the
    least-comprehensive panel: a historical 1-gene Epsilon (`['GHRL']`) beat
    multi-gene panels in the lower-detection 0041326 and mislabeled ~150K cells
    as Epsilon. Enforcement now lives in the code, not a fixed number in this
    doc: `hierarchical_phenotype.py` `build_available_markers` /
    `build_broad_markers` **fail loud** when a fine type has fewer than 2
    on-panel markers, or a broad-lineage union fewer than 3 (verify the live
    thresholds there). Rare types with thin coverage are pinned by **rule-based
    multi-criteria definitions** (`apply_rare_celltype_rules`;
    `scripts/refine_rare_celltypes.py` on the legacy path) — Gamma
    (AQP3 + ETV1 + pan-endocrine, NOT Beta/Delta) and Epsilon
    (GHRL + pan-endocrine, NOT Beta/Alpha/Delta/Gamma) — not by a single-gene
    argmax.

13. **Donor sex is a confound — handle it by gene exclusion + sex-aware
    cross-sample analysis, not reflexively as a covariate.** 0041323 is female,
    0041326 is male; Y-chromosome genes (KDM5D, DDX3Y et al.) are at ~10% /
    ~7.7% in 0041326 vs ~0% in 0041323. Invariants: (a) the full sex-chromosome
    gene set is excluded from cell-type scoring AND from PCA/scVI input via the
    `panel_for_embedding` var-mask (`scripts/finalize_panel.py`); (b) any
    analysis that **pools** male + female cells must control for sex. Adding
    `donor_sex` as a scVI `categorical_covariate` only does anything when one
    model sees both donors (a pooled/integrated run). In a **per-sample** run
    `donor_sex` is constant and the covariate is a silent no-op (the code skips
    constant covariates), so its absence there is correct, not a regression.

14. **Panel audit lives at `data/processed/panel_audit.csv`.** 106
    candidate exclusions across 12 categories (sperm/germline,
    photoreceptor, hepatocyte, cardiac/skeletal muscle, sex chromosome,
    bone, etc.) — genes biologically implausible in pancreas/LN/vessel
    tissue. Two boolean var-masks: `panel_for_scoring` (drops 100
    non-sex tissue-restricted genes) and `panel_for_embedding`
    (additionally drops 6 sex-chromosome genes). Set by
    `scripts/finalize_panel.py`.

15. **Indeterminate is never a final label — but this applies to the LEGACY
    flat-notebook path only.** The canonical hierarchical flow (#16) emits an
    Indeterminate-free `obs['celltype']` directly and does NOT run
    `resolve_indeterminate.py` or `refine_rare_celltypes.py`; the rest of this
    entry is relevant only when running the old score-genes stage-02 notebook
    without the hierarchical merge.

    In that legacy path, the strict margin / panel-coverage / spatial-coherence
    gates flag too many genuine cells as `Indeterminate*` (especially in
    low-detection 0041326). `scripts/resolve_indeterminate.py` then propagates
    labels back via a three-tier fallback: (1) suspected-from-label (the
    `Indeterminate_X*` suffix already encodes the argmax winner — the cell only
    failed a secondary gate); (2) `predicted_celltype` fallback for bare
    `Indeterminate`; (3) scVI-latent kNN (k=15, majority) as last resort.
    Provenance is preserved in `celltype_pre_resolution`, `celltype_resolution`
    ∈ {`confident`, `suspected_from_label`, `predicted_fallback`,
    `knn_propagation`}, and `celltype_neighbor_agreement`. Cluster-majority on
    `leiden_scvi` was tried and rejected (Acinar-dominant clusters hijack
    rare-type cells embedded in them).

16. **Cell typing is TWO-TIER (broad lineage → fine subtype), driven by the
    genes each cell expresses — CANONICAL since 2026-05-27 for all 4 samples.**
    - **Principle (stable).** One flat scVI-over-the-whole-panel → one leiden →
      one global argmax never separated broad lineages before subtyping (blobby
      UMAP + the flat-argmax mislabels #11/#12/#15 it spawned). Instead: (1) a
      broad tier assigns each cell to the lineage whose markers it actually
      expresses; (2) a fine tier subtypes within each lineage. Categorize by
      expressed markers — no proportion priors, no dominant/default fallback
      (Tier-C policy + memory
      `categorize-by-expressed-genes-no-proportion-priors`). Implementation is
      `scripts/hierarchical_phenotype.py`; **treat the code as the source of
      truth for the mechanism — this entry states intent, not a frozen recipe**
      (`docs-drift-verify-against-code`).
    - **Current mechanism (evolves — verify against code; IDENTITY-CORE rewire
      2026-06-08 session 2).** The FIRST scANVI redesign (SEED_MARKERS-seeded scANVI
      whose `predict()` posterior was the label) produced a **CATCH-ALL** (~28% of
      cells labeled lineages they don't express — Immune_B 1%-concordant). It is
      RETIRED. The broad tier is now driven STRICTLY by the upstream ontology-curated
      IDENTITY_CORE (`scripts/tissue_markers/_identity_core`; `var['subpanel_<L>']` =
      the clean per-lineage cores), NOT MARKER_DICT / SEED_MARKERS / the fat subpanels.
      Five hard-won lessons (verify the live params in `hierarchical_phenotype.py`):
      (1+2) per-lineage score = **ABSOLUTE mean lognorm over identity_core, NEVER
      z-scored** — a fat-panel mean is not an identity score, and z-scoring forces a
      uniform argmax that vacuums marker-negatives into rare classes (the catch-all).
      A runtime **specificity filter** (`specificity_filter_cores`, `ENRICH_MIN`) drops
      true-but-non-discriminative core genes (plasma `XBP1` = pan-ER at 73% detection)
      via abundance-independent enrichment vs a rough full-core argmax — the cores are
      cross-lineage-name-clean upstream but still carry these. (3) detection/anchoring
      use a **RAW transcript-count floor on the `counts` layer** (`IC_COUNT_FLOOR`=1; a
      lognorm gate drifts with cell size and low-RNA thymocytes sit at count 1). (4)
      scANVI is used as designed: a **CLEAN ANCHOR** = the abs-score argmax with a
      margin AND ≥`ANCHOR_MIN_GENES` (=2) specific genes detected; ambiguous/low-signal
      → `'Unknown'`; starved lineages rescued; scANVI trains on the **specificity-
      filtered** identity feature union and PROPAGATES to the Unknowns; the label is the
      EXPRESSION-MASKED `predict()` argmax (a cell can only be a lineage it expresses);
      lowsignal cells take the abs-score argmax (NOT scANVI's free pick). (5) NO
      MARKER_DICT/SEED_MARKERS in the broad call. `broad_confidence` = max posterior
      (0–1, reported); `broad_n_identity_genes` = the concordance instrument. Per
      lineage the FINE label is the detection-gated argmax over SPECIFIC subtype scores
      (`DETECTION_SUBTYPE_LINEAGES`), now sourced from the identity_core `subtype_hint`
      panels where coverage ≥ `IC_FINE_MIN` (tissue `IC_HINT_FINE`) and MARKER_DICT only
      as the thin-fine fallback; raw-count gates (`RAW_DETECT_FLOOR`=1; dominant
      Acinar→Ductal canary `DOMINANT_CANARY_FLOOR`=2). State (Tier-3) programs are
      scored into `obs['state_<group>_score']` and NEVER touch identity. Validated on
      50–60k subsamples: weighted marker-concordance 99.7% (pancreas) / 97.6% (thymus),
      catch-all collapsed, Indeterminate-free. (Retired: flat scVI, `embedding_panel`,
      cluster-mean subpanel scoring, z-scored rich-subpanel seeding, SEED_MARKERS
      seeding, kNN denoise, `_identity_trim` from the mask path.)
    - **Hard contracts (firm — downstream depends on these).** `obs['celltype']`
      is final and Indeterminate-free; the hierarchical pass SUBSUMES
      `fix_acinar_argmax.py`, the rare pass of `refine_rare_celltypes.py`, and
      `resolve_indeterminate.py` (#15) — do NOT run those after a hierarchical
      merge. `obs['celltype_lineage']` STAYS FINE (insulitis + stage-07 `isin`
      it against fine endocrine names — never put broad labels there).
      `--output inplace` is additive/non-lossy (flat `celltype_flat` + flat
      `obsm['X_scvi']`/`X_umap` + `obs['leiden_scvi']` preserved) but the writer
      OVERWRITES with **no auto-backup** — back up `{s}_phenotyped.h5ad` first.
      Masks come from the vendored **identity_core** (`data/reference/subpanels/identity_core/`,
      pinned XeniumPanelExplorer @ **076d1192** — PR #7's per-panel ontology-curated
      cores + the `panel_roles.csv`/`panel_gene_roles.csv` manifests) + the tissue
      module's `BROAD_TO_FINE`/`IC_HINT_FINE` (and `MARKER_DICT` only for the thin-fine
      fallback); `finalize_panel.py` stamps `var['subpanel_<L>']` (= identity_core,
      NO trim — the cores are pre-clean), `celltype_panel`, **`identity_panel`**
      (identity_core union; the broad call further specificity-filters it at runtime),
      `is_state_only`, `identity_promoted`. Back up `{s}_phenotyped.h5ad` first
      (`data/processed/_archive_2026-06-08_preredesign2/`).
    - **Validation = structural only** (Tier-C policy — no proportion gates;
      distributions/transitions REPORTED, never gated). The library-size **leak
      gate was REMOVED 2026-06-03** (user-directed): `validate()` no longer computes
      Spearman(broad-UMAP, total_counts). The *design* invariant that prevents the
      leak — never pass `total_counts` as a scVI covariate (#1) — is UNCHANGED; only
      the post-hoc gate is gone.
      `validate()` in `hierarchical_phenotype.py`,
      `scripts/validate_hierarchical_pilot.py`).
    - **Operational gotchas.** Run with `--dl-workers 0` (now the default; any
      `--dl-workers > 0` deadlocks the scANVI semi-supervised finetune at epoch 1
      even with `persistent_workers=False`, GPU idle). The scANVI finetune's
      `n_samples_per_label` is computed **per sample** at runtime
      (`SCANVI_N_PER_LABEL_FRAC=0.10` × the rarest broad `_seed` lineage count,
      floored at 50) — NOT a fixed constant (user-directed 2026-06-03, replacing
      the arbitrary 100). There is **no `--reuse-scvi` flag** — every run retrains
      scANVI from scratch (the feature set is `identity_panel`; a cached latent on a
      different gene set is invalid anyway). `finalize_panel.py` must run BEFORE
      `hierarchical_phenotype.py` to stamp `identity_panel` (else `ensure_celltype_panel`
      falls back to `celltype_panel`). `--output inplace` overwrites with no backup
      → back up first.

17. **3-tier panel framework: identity → mask, subclass → fine, state → score
    (NEVER a lineage axis). CANONICAL since 2026-06-08.** Upstream XeniumPanelExplorer
    (@ **076d1192**) classifies every subpanel as identity / subclass / state (Morris
    2019; Heumos 2023), machine-readable in vendored `panel_roles.csv` +
    `panel_gene_roles.csv`, AND ships a per-panel ontology-curated **identity_core**
    (`data/reference/subpanels/identity_core/`, cols gene,subtype_hint,evidence,note +
    `identity_audit_<t>.csv` det_max). The consumer config is DERIVED from it:
    `scripts/tissue_markers/_panel_roles.py` (allowlist/embedding_extra/state_panels/
    promotions) + `scripts/tissue_markers/_identity_core.py` (the core loader);
    `finalize_panel.py` calls `assert_consistent` (fail-loud if the vendored manifest ≠
    the module's `SUBPANEL_ALLOWLIST` — drift guard).
    - **identity (`broad_mask`)** → `var['subpanel_<L>']` = the clean **identity_core**
      per lineage (the over-inclusive fat subpanels are NO LONGER used; `_identity_trim`
      is RETIRED — the upstream cores are gene-by-gene CL-curated + specificity-audited
      cross-lineage-clean). `identity_panel` = the identity_core union (the scANVI
      feature seed). The broad CALL then applies a runtime data-driven SPECIFICITY
      filter (#16, `specificity_filter_cores`) because the cores still carry true-but-
      non-discriminative genes (plasma `XBP1` = pan-ER); the FINE tier reads the core
      `subtype_hint` (`IC_HINT_FINE`).
    - **state (`module_score`)** → `scripts/state_panels.py` scores each program into
      `obs['state_<group>_score']` (IFN, cytokine, ECM, proliferation, apoptosis,
      pathways, pancreas `38_insulin_biology_T1D`). The GUARD: state-only genes never
      enter `identity_panel` (excluded by construction — they are in no identity panel;
      `var['is_state_only']` audits this) and never seed/argmax a lineage. So an
      IFN-high / stressed / cycling cell keeps its identity and its activation is read
      off the orthogonal `state_*` score.
    - **promotion (Rec 3)** — genes curated `identity` inside a non-mask panel are
      promoted into their mask: thymus mTEC `AIRE`/`FEZF2` (in 11_selection /
      19_AIRE_TRA) → TEC identity (`var['identity_promoted']`).
    - Refresh the vendored snapshot + manifests with `python scripts/fetch_subpanels.py`
      (tier-driven: reads `panel_roles.csv`, fetches identity/state/fine_tier for both
      tissues). Legacy `scripts/program_panels.py` (annotation-only program figure) is
      superseded by `state_panels.py`.

## Tier-C policies are revisable — ask on conflict

The decisions above are not a single tier. **A** = invariants
(correctness / reproducibility / environment); **B** = context-bound lessons;
**C** = revisable policy / philosophy; **D** = catch-ups to the code. "Do not
regress" applies in full force only to Tier A.

The **Tier-C** stances — proportion-free validation and
categorize-by-expressed-genes (memory `no-proportion-based-quality-judgments`,
`categorize-by-expressed-genes-no-proportion-priors`; the validation half of
#16; #15's Indeterminate stance) — encode the user's current *scientific
stance*, not a correctness invariant. **Protocol when a request appears to
conflict with a Tier-C policy:** do not silently comply, and do not silently
refuse on the policy's authority — surface the conflict, ask the user to clarify
intent, AND ask whether the policy itself should be updated. The user is the
authority here; the policy is a default to confirm, not a veto.

## Repo layout (as absorbed; see README for what was dropped)

```
notebooks/
├── 00_ingest_xenium.ipynb       # Xenium output dir → spatialdata zarr (one-time)
├── 01_preprocessing_v2.ipynb    # OPTIONAL — produces _preprocessed.h5ad
├── 02_phenotyping.ipynb         # ingests zarr, scVI, raw labels (FINALIZED by `phenotype`)
├── 03_spatial_analysis.ipynb    # squidpy spatial stats → _spatial_analysis.h5ad
├── 04_group_comparisons.ipynb   # group-level DE (opt-in)
├── 05_tissue_comparisons.ipynb  # cross-tissue integration (opt-in)
├── 06_xenium_phenocycler_integration.ipynb   # (opt-in)
└── 07_immune_islet_proximity.ipynb            # pancreas immune scVI (the `final` stage)

scripts/  — 15 canonical scripts (dev-phase one-shots/probes/renderers pruned 2026-06-09)
├── run_pipeline.sh               # SINGLE driver, raw→final; STAGES="01 02 phenotype 03 final"
├── env_rapids.sh                 # CUDA/rapids LD_PRELOAD env (sourced by the driver)
├── finalize_panel.py             # stamp identity_core masks      ┐ the `phenotype` stage
├── hierarchical_phenotype.py     # two-tier phenotyping; calls    │ (folds in the former
│                                 #   state_panels internally      │  run_identity_core_redesign.sh)
├── state_panels.py               # state-program scores (dep)     │
├── validate_redesign.py          # structural validation (gate)   │
├── render_stage02_phenotyping.py # stage-02 figures (+program_panels.py dep) ┘
├── insulitis_analysis.py         # pancreas `final` (consumes nb07 immune h5ad)
├── thymic_zonation_analysis.py   # thymus `final`
├── fetch_subpanels.py / audit_panel_exclusions.py   # repro: vendor cores / gen panel_audit.csv
├── export_xenium_explorer_groups.py / slim_phenotyped.py / ingest_thymus.py   # utilities
├── tissue_markers/               # per-tissue config: pancreas.py, thymus.py, _identity_core.py,
│                                 #   _panel_roles.py, _shared.py, verify_panel_coverage.py
├── benchmark/                    # spatial-method benchmark sub-project (separate envs)
├── slurm/                        # SLURM array jobs (HPC equivalent of the notebook stages)
└── R/xenium_analysis.R           # benchmark R helper

utils/                            # analysis_utils.py, xenium_explorer_export.py, ligrec_plot.py

data/
├── raw/{sample}.zarr             # Xenium spatialdata zarr stores (the canonical ingest source)
└── processed/{sample}/           # _phenotyped.h5ad, _spatial_analysis.h5ad, _immune_phenotyped.h5ad, CSVs

figures/{sample}/{01_preprocessing,02_phenotyping,03_spatial_analysis,07_*}/   # rendered figures
logs/                             # per-stage .err/.out, pipeline_*.log
notebooks/executed/{sample}/      # nbconvert-executed notebook copies (debug)
```

## Common operations

> **Scope note, 2026-08-28.** This section describes the **HiPerGator** workflow as absorbed,
> and it is still correct there. On this workstation it does not run: `run_pipeline.sh` and
> `env_rapids.sh` moved to `hpc/`, the `xenium_analysis` conda env is broken (see README), and
> `ssh hpg.rc.ufl.edu` currently returns `Permission denied (publickey,keyboard-interactive)`.
> The two stages that *do* run here are `finalize_panel.py` then `hierarchical_phenotype.py`,
> driven directly under `conda run -n scvi-env` — see README (2026-09-02: `--paths` exists on
> `finalize_panel.py` only; the typer's own flag is pending, README § The typer).

`hpc/run_pipeline.sh` is the SINGLE raw→final driver. STAGES tokens:
`01` `02` (notebooks) → `phenotype` (finalize_panel → hierarchical_phenotype
`--output inplace` → validate_redesign → render figures; **auto-runs after `02`**,
suppress with `XENIUM_SKIP_MERGE=1`) → `03` (spatial) → `final` (pancreas: nb07 →
insulitis_analysis.py; thymus: thymic_zonation_analysis.py). Opt-in: `explorer`,
`05`. Default STAGES = `01 02 phenotype 03`. The `phenotype` stage backs up each
`_phenotyped.h5ad` to `data/processed/_archive_<ts>_premerge/` before the in-place
merge (the writer has no auto-backup).

### Run the whole pipeline (raw → final) from scratch
```bash
cd /blue/maigan/smith6jt/Xenium_Analysis
STAGES="01 02 phenotype 03 final" TISSUE=pancreas nohup bash scripts/run_pipeline.sh \
    > logs/pipeline_$(date +%Y%m%d_%H%M%S).log 2>&1 &
STAGES="01 02 phenotype 03 final" TISSUE=thymus   bash scripts/run_pipeline.sh
```
B200 GPU + 64 CPU: stage 01 ≈ 30–40 min/sample, 02 ≈ 60–75, hierarchical ≈ 28
(pancreas), 03 ≈ 15–20, nb07 ≈ 20, insulitis/zonation ≈ 5–15.

### Regenerate only the downstream finals (phenotyping already merged)
```bash
STAGES="03 final" TISSUE=pancreas bash scripts/run_pipeline.sh
STAGES="03 final" TISSUE=thymus   bash scripts/run_pipeline.sh
```

### Re-run phenotyping after a notebook / marker-dict edit
```bash
STAGES="02 phenotype" TISSUE=pancreas bash scripts/run_pipeline.sh   # 02 auto-couples phenotype
```

### Run a single stage on a single sample
```bash
bash scripts/run_pipeline.sh 03 0041323
```

### Run the spatial-method benchmark (separate sub-project)
Compares BANKSY/UTAG (Python), KODAMA/MISTy (R) + the stage-03 niche baseline on
`0041323`+`0041435`. Uses **isolated envs** (`bench_py` + R pkgs in `R_libs`) —
the `xenium_analysis` env is never touched. Build once, then run:
```bash
bash scripts/benchmark/setup_envs.sh                    # one-time env build
STEPS="prep methods eval" SAMPLES="0041323 0041435" VARIANTS="crop50k" \
    bash scripts/benchmark/run_benchmark.sh             # or sbatch slurm/bench.sh
```
KODAMA is the slow method (~45 min/25k cells); cap sweep runs with
`BENCH_TIMEOUT`. Full docs: `scripts/benchmark/README.md`. Outputs:
`data/processed/benchmark/{s}/{s}_benchmark_scorecard_{v}.csv` + `_report_{v}.md`
+ `figures/{s}/benchmark/*.png`.

### Re-render figures from saved h5ad without re-running anything
```bash
source /blue/maigan/smith6jt/miniforge3/etc/profile.d/conda.sh
conda activate xenium_analysis
python3 /tmp/regen_phenotyping_figs.py 0041323  # if the script still exists
```

### Validate a phenotyped h5ad
```python
import scanpy as sc
from scipy.stats import spearmanr
a = sc.read_h5ad("data/processed/0041326/0041326_phenotyped.h5ad")
um = a.obsm.get("X_umap_scvi", a.obsm["X_umap"])
print(spearmanr(um[:, 0], a.obs["total_counts"]))  # |ρ| should be < 0.3
print(a.obs["leiden_scvi"].nunique())              # ≈ 15-20 at resolution 1.0
print(a.obs["celltype"].value_counts())            # cell-type breakdown
```

### Verifying GPU paths actually fired
```bash
# After a run finishes:
grep -E "TF32 enabled|rsc.pp|cuML KMeans|cuNearestNeighbors|GPU:" \
    notebooks/executed/0041323/02_0041323.ipynb | head -10
```
If you see only `sc.pp.*` and no `rsc.*` / `cuML` / `TF32` lines, the
LD_PRELOAD env fix in the driver script wasn't picked up and rapids
imports silently failed (notebook fell back to CPU). Re-source the
driver, then verify:
```bash
LD_PRELOAD=$CONDA_PREFIX/lib/libstdc++.so.6:/apps/compilers/cuda/12.8.1/lib64/libcudart.so.12 \
python -c "import rapids_singlecell as rsc; print(rsc.__version__)"
```

## Background on the original failure mode

The pre-pivot pipeline (April 22-24) had every UMAP dominated by a
library-size gradient — `total_counts` varied monotonically along the
axes. Root causes that were fixed:

1. `MIN_COUNTS=5 / MAX_COUNTS=1500` (cell 14 of old preprocessing) —
   too loose at the bottom (single-transcript spike cells survived) and
   truncating the upper biological tail.
2. HVG selection AFTER PCA — dead code; PCA used all genes.
3. scVI without `cell_area` covariate — segmentation-area effects bled
   into the latent.
4. `svd_solver='arpack'` on 1.2 M cells (initial fix attempt) — Lanczos
   inner loop is single-threaded; one run wedged for 50+ min before
   producing anything. Switched to `svd_solver='randomized'` with
   pinned `random_state=0`.
5. Rapids/cuML UMAP produced extreme outliers — moved to scanpy CPU.
6. Default `n_jobs=1` on Moran's I and centrality — single-threaded
   permutation tests on 1.2 M cells were infeasible. Centrality
   skipped, Moran's I subsampled to 200 k.

Quantitative validation: Spearman correlation of each UMAP axis vs
`total_counts` is now {0.067, 0.185} for 0041323 and {0.125, 0.140}
for 0041326 — both well under the 0.3 threshold for "library size is
NOT driving the embedding."

## When in doubt

- `docs/HANDOFF_2026-07-13.md` (historical) and the repo-root `docs/HANDOFF_*.md` / `docs/WORKFLOW_2026-09-02_MAP.md` → "what just happened, how to resume"
- This file (CLAUDE.md) → "how the repo is structured"
- README.md → user-facing docs
- DATA_README.md → expected data formats, h5ad contracts, headline CSVs
