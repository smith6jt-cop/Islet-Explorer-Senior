"""Stage 2.4 — two-tier hierarchical cell typing (broad lineage first, then subtype).

Replaces the FLAT scheme (one scVI latent over the full ~4.9k panel → one leiden →
one argmax over all fine types) with a hierarchy driven by the vendored cell-type
*subpanels* (data/reference/subpanels/, pinned XeniumPanelExplorer @ 79c5ab9):

  TIER 1 (coarse / broad lineage) — IDENTITY-CORE + clean-anchor scANVI
    (rewire 2026-06-08 session 2; fixes the first redesign's CATCH-ALL — ~28% of
    cells assigned lineages they don't express). Driven STRICTLY by the upstream
    ontology-curated IDENTITY_CORE (var['subpanel_<L>'] = the clean per-lineage cores
    stamped by finalize_panel from XeniumPanelExplorer @076d1192) — NOT MARKER_DICT,
    SEED_MARKERS, or the fat subpanels. Four lessons drive it:
      1+2  per-lineage score = ABSOLUTE mean lognorm over the on-panel identity_core
           (NEVER z-scored — a fat-panel mean is not an identity score, and z-scoring
           it forces a ~uniform argmax that vacuums marker-negatives into rare classes).
      3    anchors/detection use a RAW transcript-count floor on the `counts` layer
           (a lognorm gate drifts with cell size — 81.5% of nonzero entries are single
           transcripts).
      4    a CLEAN ANCHOR detects >= ANCHOR_MIN_GENES of EXACTLY ONE lineage's
           identity_core; ambiguous / low-signal cells -> 'Unknown'. Those anchors
           train a semi-supervised scANVI on the identity_core FEATURE union
           var['identity_panel'] (state/pathway-only genes excluded — the GUARD),
           latent in obsm['X_scvi_broad']; predict(soft=True) PROPAGATES to the
           Unknowns. The label is the EXPRESSION-MASKED predict() argmax — a cell may
           only be a lineage it expresses (>= 1 identity_core gene at the count floor);
           among expressed lineages the posterior arbitrates ('beyond the dictionary').
    broad_confidence = max posterior (0-1); broad_n_identity_genes = identity_core
    genes of the called lineage detected (the concordance instrument). State (Tier-3)
    programs are scored into obs['state_<group>_score'] and NEVER enter the identity
    decision. (Retired: flat scVI, the embedding_panel, cluster-mean subpanel scoring,
    z-scored rich-subpanel seeding, SEED_MARKERS seeding, and kNN denoise.)

  TIER 2 (fine / within-lineage subtype)
    For each broad lineage with >1 fine type and >= N_MIN_SUBCLUSTER cells: a
    SEPARATE sub-PCA + sub-leiden on that lineage's subpanel genes, then assign
    each sub-cluster to the fine type (tight MARKER_DICT panel) with the highest
    cluster-mean fine score. Singleton lineages take their one fine type; tiny
    lineages fall back to per-cell score-argmax. Within Exocrine, an anti-acinar
    (Acinar→Ductal) refinement is applied. After Endocrine subtyping, the
    tissue module's apply_rare_celltype_rules() adds Gamma/Epsilon (the
    multi-criteria logic that fixed the original Epsilon bug).

This subsumes scripts/fix_acinar_argmax.py and the rare-rules pass of
scripts/refine_rare_celltypes.py: it produces a complete, Indeterminate-free,
FINAL celltype, so those two are NOT run after it in the hierarchical chain.

Output (default --output slim): a compact sidecar
data/processed/{s}/{s}_hierarchical.h5ad carrying ONLY the new obs columns +
obsm['X_scvi_broad','X_umap_broad'] (+ the old flat 'X_umap' for side-by-side) —
fast to write/read and non-destructive, for review before the expensive 60 GB
in-place merge. --output inplace merges into {s}_phenotyped.h5ad (preserving the
old label as obs['celltype_flat'] and the flat X_scvi/leiden_scvi/X_umap).

Contract (CLAUDE.md / DATA_README): canonical obs['celltype'] stays a FINE label;
obs['celltype_lineage'] stays a FINE column (insulitis + stage-07 isin it against
fine endocrine names — the binding compatibility constraint). New columns:
celltype_broad, leiden_broad, leiden_sub_<L>, celltype_fine, celltype_hier,
broad_confidence, subtype_confidence; new obsm: X_scvi_broad, X_umap_broad.

Run (monorepo, this workstation):
    conda run -n scvi-env python analysis/xenium/scripts/hierarchical_phenotype.py \
        --tissue pancreas --output slim \
        --subpanel-dir data/development/xenium_panels/subpanels \
        --paths /data/xenium/processed/<serial>/panc/<serial>_phenotyped.h5ad

The env is `scvi-env`, NOT `xenium_analysis`.  Measured 2026-08-28: `xenium_analysis` is
broken here -- `import scanpy` and `import scvi` both raise `ImportError: cannot import name
'_spropack' from 'scipy.sparse.linalg._propack'` from a mixed scipy/numpy install.  `scvi-env`
carries scvi-tools 1.4.0.post1, scanpy 1.11.5, torch 2.9.0+cu128 with CUDA available, and
every API this script calls.

`rapids_singlecell` is ABSENT here and is not needed: `try_import_rsc()` returns None with a
warning and every call site has a complete CPU branch.  `scripts/env_rapids.sh` and its
LD_PRELOAD recipe are HiPerGator-only and now live in `hpc/`.
"""
from __future__ import annotations

import argparse
import importlib
import sys
import time
import os
from pathlib import Path

import numpy as np
import pandas as pd
import scanpy as sc

# Two roots, deliberately separate.  ROOT is where the CODE and the tracked tables live;
# DATA_ROOT is where the multi-GB h5ads live and is NEVER inside the repo.  Collapsing them
# -- as the single HiPerGator ROOT did -- would put data/processed/ under analysis/xenium/,
# and the repo-root .gitignore pattern "/data/" is ANCHORED, so it would not be ignored.
ROOT = Path(__file__).resolve().parents[1]                     # analysis/xenium
REPO = Path(__file__).resolve().parents[3]                     # monorepo root
DATA_ROOT = Path(os.environ.get("XENIUM_DATA_ROOT", "/data/xenium"))
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from scripts.workflow_paths import resolve_phenotyping_jobs

DEFAULT_SAMPLES = ("0041323", "0041326")
DEFAULT_SUBPANEL_DIR = REPO / "data/development/xenium_panels/subpanels"

# ---- fixed parameters (reproducibility: CLAUDE.md #6) ----
RANDOM_STATE = 0
# Broad embedding = scANVI (semi-supervised), seeded by the per-cell marker broad
# call. 2026-06-01 bake-off: scANVI beat vanilla scVI AND scaled-PCA on HELD-OUT
# cells (pancreas .907, thymus .710); vanilla scVI lost (single-donor / no batch →
# scVI's regularization only costs resolution; the marker labels are the lever).
SCVI_N_LAYERS = 2
SCVI_DROPOUT = 0.1
SCVI_MAX_EPOCHS = 400       # unsupervised pretrain (early-stopped)
SCVI_PATIENCE = 40
SCANVI_N_LATENT = 10        # low-dim latent → crisp separation of discrete lineages (bake-off)
SCANVI_MAX_EPOCHS = 60      # semi-supervised finetune (pretrain does the heavy lifting; longer
                            # finetune diverged to NaN on the 1.1M pancreas — keep it short)
# Per-class subsampling per epoch (scvi-tools best practice) so the dominant
# Acinar/DP classes don't swamp the semi-supervised classifier vs rare lineages.
# Determined PER SAMPLE at runtime (user-directed 2026-06-03, replacing the
# arbitrary fixed 100): n_samples_per_label = FRAC * the rarest broad-lineage
# (_seed) class count in that sample, floored. See train_scanvi.
SCANVI_N_PER_LABEL_FRAC = 0.10   # cap per broad lineage = 10% of the rarest lineage
SCANVI_N_PER_LABEL_FLOOR = 50    # never below this (degenerate-class guard)
SCVI_BATCH = 4096           # proven stable on the 1.1M pancreas (May-31); 8192 → NaN divergence
# scVI dataloader workers. 2026-06-02: DL>0 DEADLOCKS the scANVI semi-supervised
# finetune at epoch 1 (worker socket-wait) EVEN with persistent_workers=False - the
# SCVI pretrain survives it, SCANVI.train does not. DL=0 (single-process loader) is
# deadlock-free AND faster here (GPU densify via load_sparse_tensor: pretrain 16.8 vs
# 36 min on the 1.1M pancreas). Default 0; --dl-workers overrides.
DL_WORKERS = 0
# Load sparse count tensors to GPU and densify ON the GPU (avoids the CPU-densify
# bottleneck that idles the B200). GPU-only; disable with --no-sparse-load.
LOAD_SPARSE = True
# NO scVI/scANVI covariates: these are per-sample (donor_sex constant), the NB size
# factor already models depth (#1), and the technical categoricals
# (segmentation_method/nucleus_multi/z_level_cat) DEGRADED separation in the bake-off.
SCVI_CAT_COVARIATES = []
SCVI_CONT_COVARIATES = []
# embedding / clustering
NEIGHBORS_K = 10
UMAP_MIN_DIST = 0.05
UMAP_SPREAD = 2.5
METRIC = "cosine"
UMAP_FIG_SUBSAMPLE = 200_000   # diagnostic broad UMAP is figure-only → subsample (not full 1.2M CPU UMAP)
# Broad diagnostic UMAP params — DELIBERATELY looser than the fine/flat UMAP above.
# The fine params (min_dist=0.05, spread=2.5, k=10) are tuned for the 1.2M-cell fine
# UMAP; reused on the ~8-lineage broad latent they over-compress it into a blob and
# detach micro-clusters into rings. These looser values spread the lineages. The
# broad UMAP is diagnostic-only (downstream drops X_scvi_broad), so this is cosmetic.
BROAD_NEIGHBORS_K = 30
BROAD_UMAP_MIN_DIST = 0.5   # picked from a 4-combo sweep (0041323): 0.5/k=30 gives the
BROAD_UMAP_SPREAD = 1.0     # cleanest lineage territories; k=50 over-fragments Exocrine.
SUB_LEIDEN_RES = 0.5     # within-lineage subtyping resolution
# rare-rules thresholds (mirror scripts/refine_rare_celltypes.py) — LOGNORM, consumed
# by the tissue apply_rare_celltype_rules() multi-criteria pass (tuned; left on
# lognorm — those ultra-rare rules are not the catch-all and not part of the identity
# anchoring lesson-3 scope).
LOG_DETECT = float(np.log1p(3.0))
LOG_HIGH = float(np.log1p(5.0))
LOG_PRESENT = float(np.log1p(1.0))

# IDENTITY-CORE broad call (rewire 2026-06-08 session 2). The broad lineage is driven
# by the upstream ontology-curated identity_core (var['subpanel_<L>'] stamped by
# finalize_panel), NOT MARKER_DICT / SEED_MARKERS / the fat subpanels, and the score
# is the ABSOLUTE mean lognorm (NEVER z-scored — lessons 1+2: a fat-panel mean is not
# an identity score; z-scoring forces a uniform argmax that vacuums marker-negatives
# into rare classes = the catch-all). Anchors/detection use a RAW transcript-count
# floor on the `counts` layer (lesson 3: lognorm gates drift with cell size — 81.5% of
# nonzero gene-in-cell entries are single transcripts):
#   IC_COUNT_FLOOR    raw count >= this to "detect" an identity_core gene. =1: a
#                     non-zero raw count of a SPECIFIC (specificity-filtered) gene is
#                     biologically meaningful in FFPE Xenium (the lesson's own note for
#                     CD3/CD4/CD8; 81.5% of nonzero entries are single transcripts —
#                     and low-RNA thymocytes sit almost entirely at count 1, so a >=2
#                     floor made 45% of thymus 'lowsignal'). ANCHOR_MIN_GENES=2 (>=2
#                     distinct specific genes) supplies the noise protection instead.
#   ANCHOR_MIN_GENES  a CLEAN anchor must detect >= this many of its lineage's
#                     identity_core genes at the count floor (decision #12: no
#                     single-gene call) — the absolute-expression floor on top of the
#                     score argmax.
#   ANCHOR_MARGIN_FRAC a CLEAN anchor's winning lineage absolute score must beat the
#                     runner-up by >= this fraction of the winner (clean separation);
#                     ambiguous / low-signal cells -> 'Unknown' (scANVI places them by
#                     the learned latent — lesson 4: clean anchors + Unknown +
#                     propagate, NOT every-cell seeds). The score is the ABSOLUTE mean
#                     lognorm (lesson 2) — robust to a few broadly-expressed core genes
#                     (XBP1 etc.) which the mean dilutes; raw-detection-count anchoring
#                     is NOT robust to them (the catch-all returned when it was used).
#   SEED_FLOOR        rescue target — any lineage below this many clean anchors is
#                     topped up with its highest-absolute-score Unknown cells (no
#                     lineage silently lost); also the scANVI n_samples_per_label floor.
#   IC_FINE_MIN       a fine type uses its identity_core subtype_hint panel when >=
#                     this many genes are on-panel, else the MARKER_DICT thin-fine
#                     fallback (the plan's sanctioned escape hatch).
#   RAW_DETECT_FLOOR  raw-count floor for the Tier-2 fine detect-gated argmax + cycling
#                     (=1; the original lognorm gate was log1p(1) ≈ count 1, now on the
#                     size-stable raw counts layer — lesson 3).
#   DOMINANT_CANARY_FLOOR  raw-count floor for the dominant-branch Acinar→Ductal canary
#                     reassignment (=2; higher than the fine gate to avoid flipping
#                     acinar to ductal on ambient single-transcript KRT19 overspill).
#   ENRICH_MIN / ENRICH_MIN_KEEP — data-driven SPECIFICITY filter on the identity_core
#     (lesson 2: "absolute mean over SPECIFIC genes"). The upstream cores are
#     cross-lineage-name-clean but still carry true-but-non-discriminative genes
#     (plasma XBP1 / ER RPN1/PDIA4/HYOU1 are detected in 30-73% of ALL pancreas cells
#     and conflated ER-high acinar with ER-high plasma in the scANVI latent). A flat
#     det_max cap can't separate them from abundant-lineage markers (CUZD1 79%); the
#     fix is an abundance-INDEPENDENT enrichment filter: rough-label every cell by its
#     full-core absolute-score argmax, then keep gene g in lineage L iff its mean
#     lognorm in rough-L is >= ENRICH_MIN x its mean elsewhere (always >= ENRICH_MIN_KEEP
#     top-enriched genes per lineage). The filtered cores drive scoring, detection,
#     anchoring AND the scANVI feature union.
IC_COUNT_FLOOR = 1
ANCHOR_MIN_GENES = 2
ANCHOR_MARGIN_FRAC = 0.25
SEED_FLOOR = 50
IC_FINE_MIN = 3
RAW_DETECT_FLOOR = 1
DOMINANT_CANARY_FLOOR = 2
ENRICH_MIN = 2.5
ENRICH_MIN_KEEP = 4


# =============================================================================
# GPU / backend helpers
# =============================================================================
def setup_torch_tf32():
    """Enable B200 TF32 Tensor Cores for scVI (CLAUDE.md #8). Returns use_gpu."""
    try:
        import torch
        if torch.cuda.is_available():
            torch.set_float32_matmul_precision("high")
            torch.backends.cuda.matmul.allow_tf32 = True
            torch.backends.cudnn.allow_tf32 = True
            print(f"GPU: {torch.cuda.get_device_name(0)} (TF32 enabled)")
            return True
    except Exception as exc:  # pragma: no cover
        print(f"[torch] TF32 setup skipped: {exc!r}")
    print("GPU: not available — scVI/embedding on CPU")
    return False


def try_import_rsc():
    try:
        import rapids_singlecell as rsc  # noqa: F401
        print(f"rapids_singlecell {rsc.__version__} (GPU paths active)")
        return rsc
    except Exception as exc:  # pragma: no cover
        print(f"[WARN] rapids_singlecell unavailable ({exc!r}); CPU fallback")
        return None


def build_neighbors(adata, use_rep, key_added, rsc, n_neighbors=None):
    t0 = time.time()
    k = NEIGHBORS_K if n_neighbors is None else n_neighbors
    if rsc is not None:
        rsc.pp.neighbors(adata, use_rep=use_rep, n_neighbors=k,
                         metric=METRIC, random_state=RANDOM_STATE,
                         key_added=key_added)
        backend = "rsc(GPU)"
    else:
        sc.pp.neighbors(adata, use_rep=use_rep, n_neighbors=k,
                        metric=METRIC, random_state=RANDOM_STATE,
                        key_added=key_added)
        backend = "sc(CPU)"
    print(f"      neighbors[{key_added}] {backend}: {time.time()-t0:.1f}s", flush=True)


# =============================================================================
# Panel / marker helpers
# =============================================================================
def ensure_celltype_panel(adata, mod, subpanel_dir):
    """Guarantee var['celltype_panel'] + var['subpanel_<L>'] + var['identity_panel'].

    Prefer the masks stamped by finalize_panel.py; otherwise recompute on the fly
    from the vendored IDENTITY_CORE (so the script is self-contained — same source as
    finalize_panel after the 2026-06-08 session-2 rewire). identity_panel is the scANVI
    FEATURE set; a GUARD sanity check confirms no pure-state-only gene leaked into it.
    """
    if "panel_for_embedding" not in adata.var.columns:
        raise SystemExit("  panel_for_embedding absent — run finalize_panel.py first")
    have_union = "celltype_panel" in adata.var.columns
    lineages = sorted(mod.BROAD_TO_FINE)
    have_all_lineage = all(f"subpanel_{l}" in adata.var.columns for l in lineages)
    embed = adata.var["panel_for_embedding"].to_numpy()
    if have_union and have_all_lineage:
        print(f"  celltype_panel present: "
              f"{int(adata.var['celltype_panel'].sum()):,}/{adata.n_vars:,} genes")
    else:
        print("  celltype_panel/subpanel_* masks missing — recomputing from "
              f"vendored identity_core ({subpanel_dir}/identity_core)")
        from scripts.tissue_markers import _identity_core
        tissue = getattr(mod, "TISSUE", "")
        broad = _identity_core.broad_genes(tissue, mod, subpanel_dir)
        extra = _identity_core.embedding_extra_genes(tissue, mod, subpanel_dir)
        broad_to_genes = {**{L: set(g) for L, g in broad.items()},
                          "_embedding_extra": set(extra)}
        live = set(adata.var_names)
        union = set()
        for lineage, genes in sorted(broad_to_genes.items()):
            present = genes & live
            union |= present
            if lineage.startswith("_"):
                continue
            adata.var[f"subpanel_{lineage}"] = adata.var_names.isin(present)
        adata.var["celltype_panel"] = adata.var_names.isin(union) & embed
        print(f"  celltype_panel (recomputed from identity_core): "
              f"{int(adata.var['celltype_panel'].sum()):,}/{adata.n_vars:,} genes")

    # identity_panel — the scANVI feature set. Prefer finalize_panel's stamp; else
    # fall back to celltype_panel (+ promoted AIRE/FEZF2 for thymus).
    if "identity_panel" not in adata.var.columns:
        from scripts.tissue_markers import _panel_roles
        promoted = set(_panel_roles.promoted_identity_genes(getattr(mod, "TISSUE", "")))
        base = adata.var["celltype_panel"].to_numpy()
        adata.var["identity_panel"] = base | (adata.var_names.isin(promoted) & embed)
        for g in promoted:
            for L in (v for k, v in _panel_roles.promoted_identity_genes(
                    getattr(mod, "TISSUE", "")).items() if k == g):
                col = f"subpanel_{L}"
                if g in adata.var_names and col in adata.var.columns:
                    adata.var.loc[g, col] = True
        print(f"  identity_panel recomputed (fallback): "
              f"{int(adata.var['identity_panel'].sum()):,}/{adata.n_vars:,} genes")
    else:
        print(f"  identity_panel present: "
              f"{int(adata.var['identity_panel'].sum()):,}/{adata.n_vars:,} genes")
    # GUARD sanity: no state-only gene in the identity feature set.
    if "is_state_only" in adata.var.columns:
        overlap = int((adata.var["identity_panel"] & adata.var["is_state_only"]).sum())
        if overlap:
            print(f"  [GUARD][WARN] {overlap} state-only genes present in identity_panel")


def build_available_markers(adata, mod, min_fine=2):
    """MARKER_DICT filtered to genes present AND in panel_for_scoring (cell 10).

    W5 build-time gate (fail-loud): raise if any fine type has 0<survivors<min_fine
    on THIS sample's panel — the silent 1-gene-Acinar class the old version absorbed.
    Off-panel genes are reported (info, not fatal): panels legitimately differ across
    samples, but a thin survivor set must stop the run before the expensive compute.
    """
    if "panel_for_scoring" in adata.var.columns:
        keep = set(adata.var_names[adata.var["panel_for_scoring"].astype(bool)])
    else:
        keep = set(adata.var_names)
    live = set(adata.var_names)
    offpanel = {ct: [g for g in genes if g not in live]
                for ct, genes in mod.MARKER_DICT.items()}
    offpanel = {ct: g for ct, g in offpanel.items() if g}
    if offpanel:
        print(f"  [W5] MARKER_DICT genes off this sample's panel (dropped): {offpanel}")
    out, short = {}, {}
    for ct, genes in mod.MARKER_DICT.items():
        g = [x for x in genes if x in keep]
        if g:
            out[ct] = g
        if 0 < len(g) < min_fine:
            short[ct] = g
    if short:
        raise SystemExit(
            f"[W5] fine types with <{min_fine} on-panel markers after "
            f"panel_for_scoring filter: {short}. Add panel-present markers in "
            f"scripts/tissue_markers, or lower min_fine.")
    return out


def overlay_ic_fine_panels(adata, available, mod, subpanel_dir):
    """Fine score panels: PREFER the identity_core genes whose `subtype_hint` maps
    (tissue module IC_HINT_FINE) when >= IC_FINE_MIN are on-panel; else keep the
    MARKER_DICT panel (the plan's sanctioned thin-fine fallback). Returns the updated
    {fine_type: [genes]} dict and logs each fine type's source.

    The broad tier is identity_core-only; this is the ONLY place MARKER_DICT survives
    (fine subtyping), and only where identity_core's subtype_hint coverage is thin."""
    from scripts.tissue_markers import _identity_core
    live = set(adata.var_names)
    ic = _identity_core.hint_fine_panels(getattr(mod, "TISSUE", ""), mod,
                                         subpanel_dir, live_genes=live)
    out = dict(available)
    src = {}
    n_ic = 0
    for ft in available:
        icg = ic.get(ft, [])
        if len(icg) >= IC_FINE_MIN:
            out[ft] = list(icg)
            src[ft] = f"identity_core({len(icg)})"
            n_ic += 1
        else:
            src[ft] = (f"MARKER_DICT({len(available[ft])})"
                       + (f" [ic={len(icg)}<{IC_FINE_MIN}]" if icg else ""))
    print(f"  fine score panels ({n_ic}/{len(available)} from identity_core, rest "
          f"MARKER_DICT fallback):")
    for ft in sorted(out):
        print(f"      {ft:<16} {src[ft]}")
    return out


def per_cell_mean(adata, genes, layer="lognorm"):
    present = [g for g in genes if g in adata.var_names]
    if not present:
        return None
    X = adata[:, present].layers[layer]
    if hasattr(X, "toarray"):
        return np.asarray(X.mean(axis=1)).ravel()
    return np.asarray(X.mean(axis=1)).ravel()


def detect_count(adata, genes, thr, layer="lognorm"):
    """Per-cell count of `genes` detected at lognorm >= thr (repo convention,
    cf. lineage_phenotype). Returns 1D int array, or None if no gene present."""
    present = [g for g in genes if g in adata.var_names]
    if not present:
        return None
    X = adata[:, present].layers[layer]
    if hasattr(X, "toarray"):
        X = X.toarray()
    return (np.asarray(X) >= thr).sum(axis=1)


# =============================================================================
# Tier 1 — broad lineage embedding + assignment
# =============================================================================
def train_scanvi(adata, gene_mask_col, seed_label_col, use_gpu, latent_key):
    """Semi-supervised scANVI on `gene_mask_col` genes, SEEDED by `seed_label_col`
    (per-cell broad lineage; 'Unknown' for no-marker-signal cells → scANVI places
    them by expression). Stores the latent in obsm[latent_key]; returns
    ``(P, classes)`` — the FULL soft-posterior matrix (n_obs × n_classes) and its
    column order — so the caller argmaxes over KNOWN classes to DRIVE the broad
    label (redesign 2026-06-08: predict() is the label, not just an embedding).

    scANVI beat vanilla scVI AND scaled-PCA on HELD-OUT cells in the 2026-06-01
    bake-off (the marker labels are the lever; vanilla scVI lost on these single-
    donor / no-batch data). No covariates (#1)."""
    import scvi
    scvi.settings.seed = RANDOM_STATE
    if DL_WORKERS and DL_WORKERS > 0:
        scvi.settings.dl_num_workers = int(DL_WORKERS)
        # persistent_workers=True DEADLOCKS at epoch 1 on this stack (CLAUDE.md);
        # False lets DL>0 parallelize the sparse->dense densify without the hang.
        scvi.settings.dl_persistent_workers = False
        print(f"  scVI dataloader: num_workers={DL_WORKERS}, persistent_workers=False")
    mask = adata.var[gene_mask_col].astype(bool).to_numpy()
    n_genes = int(mask.sum())
    print(f"  scANVI on {n_genes:,} {gene_mask_col} genes "
          f"(n_latent={SCANVI_N_LATENT}, nb, seed={seed_label_col})", flush=True)
    sub = adata[:, mask].copy()
    sub.obs["_seed"] = adata.obs[seed_label_col].astype(str).to_numpy()

    def _usable(keys, obs):
        out = []
        for k in keys:
            if k not in obs.columns:
                continue
            if obs[k].nunique(dropna=False) > 1:
                out.append(k)
            else:
                print(f"      [scvi] skipping constant covariate '{k}'", flush=True)
        return out or None
    cat = _usable(SCVI_CAT_COVARIATES, sub.obs)
    n_lab = int((sub.obs["_seed"] != "Unknown").sum())
    print(f"      seed: {n_lab:,}/{sub.n_obs:,} labeled, "
          f"{sub.n_obs - n_lab:,} Unknown; covariates={cat}", flush=True)
    # per-sample n_samples_per_label = FRAC * rarest broad (_seed) class count, floored.
    seed_counts = sub.obs["_seed"].value_counts()
    seed_counts = seed_counts[seed_counts.index != "Unknown"]
    n_per_label = max(SCANVI_N_PER_LABEL_FLOOR,
                      int(round(SCANVI_N_PER_LABEL_FRAC * seed_counts.min())))
    print(f"      n_samples_per_label = {n_per_label:,} "
          f"(= {SCANVI_N_PER_LABEL_FRAC:.0%} of rarest seed class "
          f"'{seed_counts.idxmin()}'={int(seed_counts.min()):,})", flush=True)
    base_kw = dict(accelerator=("gpu" if use_gpu else "cpu"), devices=1,
                   batch_size=SCVI_BATCH)
    # 1) unsupervised pretrain — load_sparse_tensor densifies on GPU (SCVI.train only)
    scvi.model.SCVI.setup_anndata(sub, layer="counts", categorical_covariate_keys=cat)
    vae = scvi.model.SCVI(sub, n_layers=SCVI_N_LAYERS, n_latent=SCANVI_N_LATENT,
                          gene_likelihood="nb", dropout_rate=SCVI_DROPOUT)
    t0 = time.time()
    vae.train(max_epochs=SCVI_MAX_EPOCHS, early_stopping=True,
              early_stopping_patience=SCVI_PATIENCE,
              load_sparse_tensor=bool(LOAD_SPARSE and use_gpu), **base_kw)
    print(f"    scVI pretrain {(time.time()-t0)/60:.1f} min", flush=True)
    # 2) semi-supervised finetune seeded by the marker call. SCANVI.train does NOT
    #    accept load_sparse_tensor; the dataloader-worker densify (dl_num_workers)
    #    keeps it fed instead.
    lvae = scvi.model.SCANVI.from_scvi_model(
        vae, unlabeled_category="Unknown", labels_key="_seed", adata=sub)
    t0 = time.time()
    lvae.train(max_epochs=SCANVI_MAX_EPOCHS, n_samples_per_label=n_per_label,
               early_stopping=True, early_stopping_patience=15, **base_kw)
    print(f"    scANVI finetune {(time.time()-t0)/60:.1f} min", flush=True)
    adata.obsm[latent_key] = np.asarray(lvae.get_latent_representation(), np.float32)
    proba = lvae.predict(soft=True)          # DataFrame: rows=cells, cols=seeded classes
    classes = [str(c) for c in proba.columns]
    P = np.asarray(proba.to_numpy(), np.float32)
    del sub, vae, lvae
    return P, classes


def score_fine(adata, available):
    """Compute tight-panel {ct}_score for each available fine type (on lognorm).
    Shared by the broad tier (broad score = mean of member fine scores) and the
    fine tier (within-lineage argmax)."""
    for ct, genes in available.items():
        sc.tl.score_genes(adata, genes, score_name=f"{ct}_score",
                          layer="lognorm", use_raw=False)


def ic_broad_from_masks(adata, mod):
    """{broad lineage: [on-panel genes]} from the identity_core var['subpanel_<L>']
    masks (stamped by finalize_panel = the clean ontology-curated cores). Excludes
    the cross-lineage Cycling program."""
    out = {}
    for L in sorted(mod.BROAD_TO_FINE):
        if L == "Cycling":
            continue
        c = f"subpanel_{L}"
        if c in adata.var.columns:
            g = adata.var_names[adata.var[c].astype(bool)].tolist()
            if g:
                out[L] = g
    return out


def ic_absolute_scores(adata, ic_broad, lineages):
    """Per-cell, per-lineage ABSOLUTE mean lognorm over the lineage's on-panel
    identity_core genes. NO scaling, NO z-score (lessons 1+2: a fat-panel mean is not
    an identity score, and z-scoring it is the catch-all's silent killer). Returns
    float32 (n_obs, len(lineages))."""
    out = np.zeros((adata.n_obs, len(lineages)), dtype=np.float32)
    for j, L in enumerate(lineages):
        genes = [g for g in ic_broad[L] if g in adata.var_names]
        if not genes:
            continue
        X = adata[:, genes].layers["lognorm"]
        out[:, j] = np.asarray(X.mean(axis=1)).ravel()
    return out


def ic_detect_counts(adata, ic_broad, lineages, floor):
    """Per-cell, per-lineage # identity_core genes detected at RAW counts >= floor
    (the `counts` layer — lesson 3: a lognorm gate drifts with cell size). Returns
    int32 (n_obs, len(lineages))."""
    out = np.zeros((adata.n_obs, len(lineages)), dtype=np.int32)
    for j, L in enumerate(lineages):
        genes = [g for g in ic_broad[L] if g in adata.var_names]
        if not genes:
            continue
        X = adata[:, genes].layers["counts"]
        X = X.toarray() if hasattr(X, "toarray") else np.asarray(X)
        out[:, j] = (np.asarray(X) >= floor).sum(axis=1)
    return out


def specificity_filter_cores(adata, ic_broad, lineages, enrich_min, min_keep):
    """Drop true-but-NON-DISCRIMINATIVE genes from each lineage's identity_core
    (lesson 2: absolute mean over SPECIFIC genes). Abundance-INDEPENDENT enrichment:
    rough-label every cell by its FULL-core absolute-score argmax (robust — the mean
    dilutes broad genes), then keep gene g in lineage L iff
    mean_lognorm(g | rough==L) >= enrich_min * mean_lognorm(g | rough!=L). Always keep
    >= min_keep top-enriched genes (a lineage is never emptied); lineages with < 30
    rough cells keep their full core (too few to estimate enrichment). Returns
    (filtered {L:[genes]}, dropped {L:[genes]})."""
    # rough labels from the full-core absolute scores
    abs_full = ic_absolute_scores(adata, ic_broad, lineages)
    rough = np.array(lineages, dtype=object)[abs_full.argmax(axis=1)]
    filtered, dropped = {}, {}
    for L in lineages:
        genes = ic_broad[L]
        inL = rough == L
        if inL.sum() < 30 or len(genes) <= min_keep:
            filtered[L] = list(genes)
            dropped[L] = []
            continue
        X = adata[:, genes].layers["lognorm"]
        X = X.toarray() if hasattr(X, "toarray") else np.asarray(X)
        mu_in = X[inL].mean(axis=0)
        mu_out = X[~inL].mean(axis=0) + 1e-6
        enr = mu_in / mu_out
        keep_mask = enr >= enrich_min
        if int(keep_mask.sum()) < min_keep:
            keep_mask = np.zeros(len(genes), bool)
            keep_mask[np.argsort(-enr)[:min_keep]] = True
        filtered[L] = [g for g, k in zip(genes, keep_mask) if k]
        dropped[L] = [g for g, k in zip(genes, keep_mask) if not k]
    return filtered, dropped


def seed_clean_anchors(abs_S, det_C, lineages, min_genes, margin_frac, seed_floor):
    """Clean-anchor seeding (lesson 4: clean anchors + 'Unknown' + propagate, NOT
    every-cell seeds). A cell is a CLEAN ANCHOR for lineage L iff:
      * L is the argmax of the ABSOLUTE mean-lognorm identity_core score (lesson 2 —
        robust to a few broadly-expressed core genes, which the mean dilutes; raw
        detection-count anchoring is NOT robust to them — it re-spawned the catch-all);
      * the winner beats the runner-up by >= margin_frac x winner (clean separation);
      * det_C[L] >= min_genes — the cell actually expresses >= min_genes of L's
        identity_core genes at the raw-count floor (decision #12: no single-gene call).
    Cells that fail (ambiguous margin / low signal) -> 'Unknown', placed by the scANVI
    latent. Any lineage below `seed_floor` clean anchors is RESCUED — topped up with
    its highest-absolute-score Unknown cells expressing >= min_genes of its genes — so
    no lineage is silently lost. Returns (seed[str array], counts{L:n}, rescued{L:n})."""
    n, nL = abs_S.shape
    Larr = np.array(lineages, dtype=object)
    rows = np.arange(n)
    order = np.argsort(-abs_S, axis=1)
    top = order[:, 0]
    sec = order[:, 1] if nL > 1 else order[:, 0]
    st = abs_S[rows, top]
    ss = abs_S[rows, sec]
    dt = det_C[rows, top]
    clean = (dt >= min_genes) & (st > 0) & ((st - ss) >= margin_frac * np.maximum(st, 1e-6))
    seed = np.where(clean, Larr[top], "Unknown").astype(object)
    counts = {L: int((seed == L).sum()) for L in lineages}
    rescued = {}
    for j, L in enumerate(lineages):
        if counts[L] >= seed_floor:
            continue
        cand = np.where((det_C[:, j] >= min_genes) & (seed == "Unknown"))[0]
        if cand.size == 0:
            continue
        take = cand[np.argsort(-abs_S[cand, j])[: max(0, seed_floor - counts[L])]]
        seed[take] = L
        rescued[L] = int(take.size)
    return seed, {L: int((seed == L).sum()) for L in lineages}, rescued


def tier1_broad(adata, mod, use_gpu, rsc, make_fig):
    print("\n--- TIER 1: broad lineage (identity_core absolute + clean-anchor scANVI) ---",
          flush=True)
    import anndata as ad
    n = adata.n_obs
    rows = np.arange(n)

    # --- 1. identity_core gene sets per lineage (clean ontology cores) -----------
    ic_full = ic_broad_from_masks(adata, mod)
    lineages = sorted(ic_full)
    nL = len(lineages)
    print("  identity_core lineages: "
          + ", ".join(f"{L}({len(ic_full[L])})" for L in lineages))

    # --- 2. SPECIFICITY filter (lesson 2: absolute mean over SPECIFIC genes) ------
    # Drop true-but-non-discriminative genes (plasma XBP1 / ER RPN1/PDIA4/HYOU1) that
    # the upstream cross-lineage-name-clean cores still carry — they conflate ER-high
    # acinar with ER-high plasma in the scANVI latent (the catch-all). Abundance-
    # independent enrichment vs a rough full-core argmax.
    t0 = time.time()
    ic_broad, ic_dropped = specificity_filter_cores(
        adata, ic_full, lineages, ENRICH_MIN, ENRICH_MIN_KEEP)
    print(f"  specificity filter (enrich>={ENRICH_MIN}, keep>={ENRICH_MIN_KEEP}): "
          f"{time.time()-t0:.1f}s")
    for L in lineages:
        nd = len(ic_dropped[L])
        ex = (f"  -dropped {nd}: {ic_dropped[L][:6]}" + ("..." if nd > 6 else "")) if nd else ""
        print(f"      {L:<12} {len(ic_full[L]):>4} -> {len(ic_broad[L]):>3} specific{ex}")

    # --- 3. ABSOLUTE mean-lognorm score (NO z-score) + RAW-count detection -------
    t0 = time.time()
    abs_S = ic_absolute_scores(adata, ic_broad, lineages)
    det_C = ic_detect_counts(adata, ic_broad, lineages, IC_COUNT_FLOOR)
    print(f"      absolute scores + raw-count detection (specific cores): "
          f"{time.time()-t0:.1f}s", flush=True)

    # --- 4. clean anchors (absolute-score argmax + margin + det floor) + rescue --
    seed, seed_counts, rescued = seed_clean_anchors(
        abs_S, det_C, lineages, ANCHOR_MIN_GENES, ANCHOR_MARGIN_FRAC, SEED_FLOOR)
    n_anchor = int((seed != "Unknown").sum())
    print(f"  clean anchors: {n_anchor:,}/{n:,} ({100*n_anchor/max(n,1):.1f}%); "
          f"{n - n_anchor:,} Unknown -> scANVI propagates")
    for L in lineages:
        extra = f"  (+{rescued[L]} rescued)" if rescued.get(L) else ""
        print(f"      seed {L:<12} {seed_counts[L]:>10,}{extra}")
    adata.obs["_scanvi_seed"] = pd.Categorical(seed,
                                               categories=lineages + ["Unknown"])

    # --- 5. semi-supervised scANVI on the SPECIFIC identity feature union --------
    # The feature set is the union of the FILTERED cores (+ promoted identity genes,
    # e.g. thymus AIRE/FEZF2) — NOT the full identity_panel, so the latent is built on
    # discriminative genes and stops conflating ER-high lineages.
    feat = set().union(*ic_broad.values())
    if "identity_promoted" in adata.var.columns:
        feat |= set(adata.var_names[adata.var["identity_promoted"].astype(bool)])
    feat &= set(adata.var_names)
    adata.var["_identity_specific"] = adata.var_names.isin(feat)
    print(f"  scANVI feature set (specific identity union): {len(feat):,} genes")
    P, classes = train_scanvi(adata, "_identity_specific", "_scanvi_seed",
                              use_gpu, latent_key="X_scvi_broad")
    adata.obs.drop(columns=["_scanvi_seed"], inplace=True)
    adata.var.drop(columns=["_identity_specific"], inplace=True)
    cpos = {c: i for i, c in enumerate(classes)}
    Pord = np.zeros((n, nL), dtype=np.float32)        # posterior aligned to `lineages`
    for j, L in enumerate(lineages):
        if L in cpos:
            Pord[:, j] = P[:, cpos[L]]

    # --- 6. EXPRESSION-MASKED predict(): a cell may only be labeled a lineage it
    # EXPRESSES (>= 1 SPECIFIC identity_core gene at the raw-count floor). Among the
    # expressed lineages the scANVI posterior arbitrates ('beyond the dictionary'); a
    # zero-expression call is structurally impossible. Cells expressing NO specific
    # gene of ANY lineage at the floor (none_elig / lowsignal) take the ABSOLUTE-SCORE
    # argmax (expression-grounded, count-1-level signal) — NOT the scANVI free argmax,
    # which defaults ambiguous low-count cells into one class (a residual mini-catch-all).
    elig = det_C >= 1
    none_elig = ~elig.any(axis=1)
    Pmask = np.where(elig, Pord, -1.0).astype(np.float32)
    lab = Pmask.argmax(axis=1).astype(np.int64)
    lab[none_elig] = abs_S[none_elig].argmax(axis=1)
    Larr = np.array(lineages, dtype=object)
    adata.obs["celltype_broad"] = pd.Categorical(Larr[lab], categories=lineages)
    adata.obs["broad_confidence"] = Pord[rows, lab].astype(np.float32)   # posterior of the call
    adata.obs["broad_confidence_prob"] = adata.obs["broad_confidence"]
    ps = np.sort(Pord, axis=1)
    adata.obs["broad_marker_margin"] = (ps[:, -1] - ps[:, -2]).astype(np.float32)
    adata.obs["celltype_confidence"] = adata.obs["broad_confidence"]
    adata.obs["broad_lowsignal"] = none_elig
    # concordance instrument: # identity_core genes of the CALLED lineage a cell
    # detects (>= IC_COUNT_FLOOR) — the redesign target is "cells labeled L express L".
    adata.obs["broad_n_identity_genes"] = det_C[rows, lab].astype(np.int16)
    overridden = int((Pord.argmax(axis=1) != lab).sum())
    called_det = det_C[rows, lab]
    conc1 = float((called_det >= 1).mean())     # catch-all detector (expresses >=1 specific)
    conc2 = float((called_det >= ANCHOR_MIN_GENES).mean())
    print(f"  expression mask overrode scANVI free-argmax for {overridden:,}/{n:,} "
          f"({100*overridden/max(n,1):.1f}%); {int(none_elig.sum()):,} lowsignal")
    print(f"  marker-concordance (cells labeled L detect >=1 of L's SPECIFIC "
          f"identity_core @ counts>={IC_COUNT_FLOOR}): {100*conc1:.1f}%  "
          f"(>=2: {100*conc2:.1f}%)")
    # per-lineage concordance (>=1 specific) — the redesign's catch-all gate
    for j, L in enumerate(lineages):
        m = lab == j
        if int(m.sum()):
            print(f"      {L:<12} n={int(m.sum()):>9,}  "
                  f">=1 specific: {100*(called_det[m] >= 1).mean():.0f}%")
    print("  celltype_broad distribution:")
    print("    " + adata.obs["celltype_broad"].value_counts().to_string()
          .replace("\n", "\n    "))

    # --- 7. Diagnostic broad UMAP on a subsample (figure only) ------------------
    # NOT downstream-consumed; full-1.2M CPU UMAP idled the GPU ~60 min. The
    # subsample gives the same visual in minutes; stored in uns
    # (index + coords), not obsm (which must span all cells). Still scanpy CPU UMAP (#4).
    if make_fig:
        n_sub = min(UMAP_FIG_SUBSAMPLE, n)
        sidx = np.sort(np.random.default_rng(RANDOM_STATE).choice(n, n_sub, replace=False))
        lat = np.asarray(adata.obsm["X_scvi_broad"][sidx], np.float32)
        usub = ad.AnnData(X=lat)
        usub.obsm["X_scvi_broad"] = lat
        build_neighbors(usub, "X_scvi_broad", key_added="broad", rsc=rsc,
                        n_neighbors=BROAD_NEIGHBORS_K)
        t0 = time.time()
        sc.tl.umap(usub, neighbors_key="broad", min_dist=BROAD_UMAP_MIN_DIST,
                   spread=BROAD_UMAP_SPREAD, random_state=RANDOM_STATE)
        adata.uns["umap_broad"] = {"index": sidx.astype(np.int64),
                                   "coords": np.asarray(usub.obsm["X_umap"], np.float32)}
        print(f"      UMAP[broad] CPU on {n_sub:,}-subsample: "
              f"{time.time()-t0:.1f}s", flush=True)


# =============================================================================
# Tier 2 — within-lineage subtyping
# =============================================================================
def sub_embed_cluster(adata_L, mask_col, rsc, use_gpu):
    """Sub-PCA + sub-leiden on a lineage's subpanel genes. Returns leiden labels
    (str array aligned to adata_L) or None if too few genes."""
    mask = adata_L.var[mask_col].astype(bool).to_numpy()
    if int(mask.sum()) < 3:
        return None
    sub = adata_L[:, mask].copy()
    sub.X = sub.layers["lognorm"].copy()
    ncomp = int(min(30, mask.sum() - 1, sub.n_obs - 1))
    if ncomp < 2:
        return None
    if rsc is not None:
        rsc.get.anndata_to_GPU(sub)
        rsc.pp.scale(sub, max_value=10)
        rsc.pp.pca(sub, n_comps=ncomp, random_state=RANDOM_STATE)
        rsc.get.anndata_to_CPU(sub)
        rsc.pp.neighbors(sub, use_rep="X_pca", n_neighbors=NEIGHBORS_K,
                         metric=METRIC, random_state=RANDOM_STATE)
    else:
        sc.pp.scale(sub, max_value=10)
        sc.pp.pca(sub, n_comps=ncomp, random_state=RANDOM_STATE)
        sc.pp.neighbors(sub, use_rep="X_pca", n_neighbors=NEIGHBORS_K,
                        metric=METRIC, random_state=RANDOM_STATE)
    sc.tl.leiden(sub, resolution=SUB_LEIDEN_RES, key_added="leiden_sub",
                 flavor="igraph", n_iterations=-1, directed=False,
                 random_state=RANDOM_STATE)
    return sub.obs["leiden_sub"].astype(str).to_numpy()


def tier2_fine(adata, mod, available, rsc, use_gpu):
    print("\n--- TIER 2: within-lineage subtyping ---", flush=True)
    # Fine {ct}_score already computed by score_fine() before Tier 1.
    broad = adata.obs["celltype_broad"].astype(str).to_numpy()
    n = adata.n_obs
    celltype_fine = np.array(["Unassigned"] * n, dtype=object)
    predicted = np.array(["Unassigned"] * n, dtype=object)
    subtype_conf = np.full(n, np.nan, dtype=np.float32)
    # Per-lineage sub-cluster ids, NaN/"" outside the lineage.
    sub_cols = {}

    # Dominant type (AMY1A-overspill-prone, e.g. pancreas Acinar): inside its
    # lineage the dominant score_genes score is inflated everywhere by ambient
    # overspill, so cells are labeled by PER-CELL DETECTION of the OTHER fine
    # types' canaries (>=1 canary gene at lognorm >= LOG_DETECT) — default is the
    # dominant type. Per-cell (not cluster-mean) because scattered alt cells in a
    # mostly-dominant section never form pure sub-clusters. Skips sub-embedding.
    dominant_type = getattr(mod, "ANTI_OVERRIDE_DOMINANT_TYPE", None)
    dominant_lineage = next((L for L, fs in mod.BROAD_TO_FINE.items()
                             if dominant_type in fs), None) if dominant_type else None
    canaries = getattr(mod, "LINEAGE_CANARY_PANELS", {})
    # Lineages whose fine call is a detection-gated score-argmax over SPECIFIC
    # subtypes (generic residual dropped from the argmax) — e.g. pancreas
    # Endocrine. Empty for tissues without the attr (thymus unaffected).
    detect_cfg = getattr(mod, "DETECTION_SUBTYPE_LINEAGES", {})

    for lin in sorted(mod.BROAD_TO_FINE):
        fine_types = sorted(mod.BROAD_TO_FINE[lin] & set(available))
        mask = broad == lin
        nL = int(mask.sum())
        if nL == 0:
            continue
        # Singleton / empty-fine lineage: one label, no sub-embedding.
        if len(fine_types) <= 1:
            lab = fine_types[0] if fine_types else lin
            celltype_fine[mask] = lab
            predicted[mask] = lab
            print(f"  {lin:<10} n={nL:>9,}  -> {lab} (singleton)")
            continue

        sub_cols_lin = [f"{ct}_score" for ct in fine_types]
        Ssub = adata.obs.loc[mask, sub_cols_lin].to_numpy()
        # Per-cell fine argmax (predicted_celltype + subtype margin).
        pc = np.array(fine_types)[Ssub.argmax(axis=1)]
        predicted[mask] = pc
        t2 = np.sort(Ssub, axis=1)[:, -2:]
        subtype_conf[mask] = (t2[:, 1] - t2[:, 0]).astype(np.float32)

        if lin == dominant_lineage:
            # per-cell detection: default to dominant_type; a cell is reassigned
            # to a non-dominant fine type when >=1 of that type's canary genes is
            # detected at RAW counts >= RAW_DETECT_FLOOR (lesson 3). Robust to
            # dominant-marker overspill and scattered alt cells; no sub-embedding.
            adata_L = adata[mask]
            lab = np.array([dominant_type] * nL, dtype=object)
            for ft in fine_types:
                if ft == dominant_type:
                    continue
                cnt = detect_count(adata_L, canaries.get(ft, available.get(ft, [])),
                                   DOMINANT_CANARY_FLOOR, layer="counts")
                if cnt is not None:
                    lab[cnt >= 1] = ft
            celltype_fine[mask] = lab
            brk = ", ".join(f"{ft}={int((lab == ft).sum()):,}" for ft in fine_types)
            print(f"  {lin:<10} n={nL:>9,}  detection-based (raw counts): {brk}")
            continue

        if lin in detect_cfg:
            # Detection-GATED score-argmax over SPECIFIC subtypes only. The
            # generic pan-endocrine fine type is intentionally DROPPED from the
            # argmax: its broad score panel is high in every islet cell and
            # out-competed Beta/Delta at cluster-mean argmax, erasing β/δ in the
            # 2026-05-26 pilot (decision #12). A cell is assigned the highest
            # {subtype}_score; cells with NO specific-marker signal at the RAW-count
            # gate (lesson 3) stay the generic `residual`. Mirrors the Exocrine/
            # thymocyte detection logic; no sub-embedding needed.
            cfg = detect_cfg[lin]
            residual = cfg.get("residual", lin)
            sts = [st for st in cfg["subtypes"]
                   if f"{st}_score" in adata.obs.columns]
            adata_L = adata[mask]
            gate = detect_count(adata_L, cfg["gate_markers"],
                                RAW_DETECT_FLOOR, layer="counts")
            has_signal = (gate >= 1) if gate is not None else np.zeros(nL, bool)
            Sspec = adata.obs.loc[mask, [f"{st}_score" for st in sts]].to_numpy()
            pick = np.array(sts)[Sspec.argmax(axis=1)]
            lab = np.where(has_signal, pick, residual).astype(object)
            celltype_fine[mask] = lab
            predicted[mask] = lab
            brk = ", ".join(f"{st}={int((lab == st).sum()):,}" for st in sts)
            brk += f", {residual}(resid)={int((lab == residual).sum()):,}"
            print(f"  {lin:<10} n={nL:>9,}  detection-gated argmax: {brk}")
            del adata_L
            continue

        if nL < mod.N_MIN_SUBCLUSTER:
            # Tiny lineage: per-cell score-argmax (leiden unstable).
            celltype_fine[mask] = pc
            print(f"  {lin:<10} n={nL:>9,}  fine={fine_types}  "
                  f"(tiny → per-cell argmax)")
            continue

        # Sub-embed + sub-cluster, then cluster-mean argmax over fine types.
        adata_L = adata[mask].copy()
        t0 = time.time()
        leiden_sub = sub_embed_cluster(adata_L, f"subpanel_{lin}", rsc, use_gpu)
        if leiden_sub is None:
            celltype_fine[mask] = pc
            print(f"  {lin:<10} n={nL:>9,}  (sub-embed skipped → per-cell argmax)")
            del adata_L
            continue
        tmp = adata_L.obs[sub_cols_lin].copy()
        tmp["_sub"] = leiden_sub
        cmean = tmp.groupby("_sub", observed=True)[sub_cols_lin].mean()
        sub_to_fine = (cmean.idxmax(axis=1).str.replace("_score", "", regex=False)
                       .to_dict())
        lab = pd.Series(leiden_sub).map(sub_to_fine).to_numpy()
        celltype_fine[mask] = lab
        # store sub-cluster ids
        full = np.array([""] * n, dtype=object)
        full[mask] = leiden_sub
        sub_cols[f"leiden_sub_{lin}"] = full
        nsub = len(set(leiden_sub))
        print(f"  {lin:<10} n={nL:>9,}  fine={fine_types}  "
              f"{nsub} sub-clusters in {time.time()-t0:.1f}s", flush=True)
        del adata_L

    # NOTE: within-Exocrine Acinar↔Ductal disambiguation is handled at the
    # sub-cluster level above (detection-based, robust to AMY1A overspill) — no
    # per-cell median flip (that over-flipped ~25% in testing).

    adata.obs["celltype_fine"] = pd.Categorical(celltype_fine)
    # predicted_celltype = the Tier-2 per-cell call (broad-constrained), NOT the
    # raw all-types score-argmax. This is the base the thymus
    # apply_thymocyte_subtyping refines (undetected thymocytes default to the
    # dominant DP, not a noisy score-argmax that over-calls DN/Treg).
    adata.obs["predicted_celltype"] = pd.Categorical(celltype_fine)
    adata.obs["subtype_confidence"] = subtype_conf
    for k, v in sub_cols.items():
        adata.obs[k] = pd.Categorical(v)

    # celltype is the canonical fine label; rare rules will refine it next.
    adata.obs["celltype"] = pd.Categorical(celltype_fine)


def apply_rare_rules(adata, mod):
    """Tissue rare-type rules on celltype (Gamma/Epsilon for pancreas)."""
    print("\n--- rare-celltype rules (within Endocrine) ---", flush=True)
    try:
        counts = mod.apply_rare_celltype_rules(adata, LOG_DETECT, LOG_HIGH,
                                               log_present=LOG_PRESENT)
    except TypeError:
        counts = mod.apply_rare_celltype_rules(adata, LOG_DETECT, LOG_HIGH)
    for ct, c in counts.items():
        print(f"    {ct}: {c:,} cells")
    return counts


# =============================================================================
# Validation + output
# =============================================================================
def validate(adata, mod=None):
    print("\n--- validation ---")
    print("  celltype_broad distribution:")
    print("    " + adata.obs["celltype_broad"].value_counts().to_string()
          .replace("\n", "\n    "))
    vc = adata.obs["celltype"].astype(str).value_counts()
    frac = (vc / vc.sum() * 100).round(2)
    print("  final celltype (%):")
    for k in vc.index[:20]:
        print(f"    {k:<24} {vc[k]:>9,}  {frac[k]:>6.2f}%")
    # scANVI broad_confidence (max posterior) — REPORT ONLY, never gated (Tier-C).
    if "broad_confidence" in adata.obs.columns:
        bcv = adata.obs["broad_confidence"].to_numpy()
        qs = np.quantile(bcv, [0.0, 0.05, 0.25, 0.5, 0.75, 0.95, 1.0])
        print("  broad_confidence (scANVI max posterior) min/5/25/50/75/95/max:")
        print("    " + " ".join(f"{q:.3f}" for q in qs)
              + f"   (<0.50: {int((bcv < 0.5).sum()):,} = {100*(bcv < 0.5).mean():.2f}%)")
        per = adata.obs.groupby("celltype_broad", observed=True)["broad_confidence"].agg(
            ["mean", "size"])
        for L, r in per.iterrows():
            print(f"      {str(L):<12} mean_post={r['mean']:.3f}  n={int(r['size']):,}")
    # Tier-3 state programs — report distribution (never an identity axis).
    state_cols = sorted(c for c in adata.obs.columns
                        if c.startswith("state_") and c.endswith("_score"))
    if state_cols:
        print(f"  state programs scored ({len(state_cols)}):")
        for c in state_cols:
            v = adata.obs[c].to_numpy()
            print(f"      {c:<34} mean={v.mean():+.3f} sd={v.std():.3f} "
                  f"p95={np.quantile(v, 0.95):+.3f}")
    # compatibility check: celltype_lineage must stay FINE (not collapsed to the
    # broad vocabulary, which downstream fine-name seeding relies on).
    if "celltype_lineage" in adata.obs.columns and "celltype_broad" in adata.obs.columns:
        cl = set(adata.obs["celltype_lineage"].astype(str).unique())
        cb = set(adata.obs["celltype_broad"].astype(str).unique())
        if cl and cl <= cb:
            print("  [WARN] celltype_lineage looks collapsed to broad labels — "
                  "downstream fine-name seeding may break!")


def render_figure(adata, sample, figdir, mod):
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except Exception as exc:
        print(f"  [fig] skipped ({exc!r})")
        return

    import anndata as ad
    ub = adata.uns.get("umap_broad")
    if ub is None:
        print("  [fig] no uns['umap_broad'] — skipped")
        return
    # Plot on the diagnostic UMAP SUBSAMPLE (uns index + coords); build a minimal
    # AnnData so we never densify the full 1.2M cells just to draw.
    idx = np.asarray(ub["index"])
    cols = [c for c in ("celltype", "celltype_broad", "broad_confidence",
                        "broad_confidence_prob", "broad_marker_margin")
            if c in adata.obs.columns]
    padata = ad.AnnData(obs=adata.obs.iloc[idx][cols].copy())
    padata.obsm["X_umapfig"] = np.asarray(ub["coords"], np.float32)

    # Use the ACTIVE tissue's palettes, and guarantee every present category has
    # a color (scanpy KeyErrors if a palette dict misses a category).
    def _safe_pal(p, cats):
        if not isinstance(p, dict):
            return None
        return {str(c): p.get(str(c), "#999999") for c in cats}
    pal = _safe_pal(getattr(mod, "PALETTE", None),
                    padata.obs["celltype"].astype(str).unique())
    bpal = _safe_pal(getattr(mod, "BROAD_PALETTE", None),
                     padata.obs["celltype_broad"].astype(str).unique())
    size = max(2.0, 6e5 / max(padata.n_obs, 1))
    figdir = Path(figdir)
    figdir.mkdir(parents=True, exist_ok=True)
    fig, axes = plt.subplots(2, 2, figsize=(20, 18))
    sc.pl.embedding(padata, basis="umapfig", color="celltype_broad", ax=axes[0, 0],
                    show=False, size=size, frameon=False, legend_loc="right margin",
                    palette=bpal, title=f"{sample} — celltype_broad (Tier 1, scANVI)")
    sc.pl.embedding(padata, basis="umapfig", color="celltype", ax=axes[0, 1],
                    show=False, size=size, frameon=False, legend_loc="right margin",
                    palette=pal, title="celltype (Tier 2, fine)")
    # scANVI posterior confidence [0,1] — interpretable (was a right-skewed margin)
    if "broad_confidence_prob" in padata.obs:
        sc.pl.embedding(padata, basis="umapfig", color="broad_confidence_prob",
                        ax=axes[1, 0], show=False, size=size, frameon=False,
                        legend_loc="right margin", vmin=0.0, vmax=1.0,
                        title="broad confidence (scANVI posterior, 0–1)")
    else:
        axes[1, 0].axis("off")
    # scANVI posterior MARGIN (top1−top2), percentile-clipped so the bulk is visible
    mcol = "broad_marker_margin" if "broad_marker_margin" in padata.obs else "broad_confidence"
    bc = padata.obs[mcol].to_numpy()
    vlo, vhi = np.percentile(bc, [2, 98])
    sc.pl.embedding(padata, basis="umapfig", color=mcol, ax=axes[1, 1],
                    show=False, size=size, frameon=False, legend_loc="right margin",
                    vmin=float(vlo), vmax=float(max(vhi, vlo + 1e-6)),
                    title="broad posterior margin (top1−top2, 2–98%ile clip)")
    for ax in axes.ravel():
        for coll in ax.collections:
            coll.set_rasterized(True)
    plt.tight_layout()
    out = figdir / "hierarchical_umap.png"
    plt.savefig(out, dpi=200, bbox_inches="tight")
    plt.close(fig)
    print(f"  [fig] wrote {out}")


# obs columns carried into the slim sidecar
_SLIM_OBS = ["celltype", "celltype_fine", "celltype_lineage", "celltype_broad",
             "celltype_hier", "celltype_cycling", "cycling", "broad_lowsignal",
             "predicted_celltype", "broad_confidence", "broad_confidence_prob",
             "broad_marker_margin", "broad_n_identity_genes", "subtype_confidence",
             "celltype_confidence", "celltype_flat", "donor_sex", "total_counts"]
# (state_<group>_score columns are carried automatically by the `_score` glob in
#  write_output's slim branch.)


def write_output(adata, sample, mode, input_path=None):
    import anndata as ad
    p_pheno = (Path(input_path) if input_path is not None else
               DATA_ROOT / f"processed/{sample}/{sample}_phenotyped.h5ad")
    if mode == "inplace":
        print(f"\n  writing IN-PLACE to {p_pheno.name} ...", flush=True)
        t = time.time()
        adata.write_h5ad(p_pheno)
        print(f"    wrote in {(time.time()-t)/60:.1f} min")
        return p_pheno
    # slim sidecar
    keep = [c for c in _SLIM_OBS if c in adata.obs.columns]
    keep += [c for c in adata.obs.columns if c.startswith("leiden_sub_")
             or c.startswith("broadscore_") or c.endswith("_score")]
    obs = adata.obs[keep].copy()
    obsm = {}
    for k in ("spatial", "X_scvi_broad", "X_umap_broad", "X_umap"):
        if k in adata.obsm:
            obsm[k] = np.asarray(adata.obsm[k])
    slim = ad.AnnData(obs=obs)
    slim.obs_names = adata.obs_names
    for k, v in obsm.items():
        slim.obsm[k] = v
    out = p_pheno.with_name(f"{sample}_hierarchical.h5ad")
    print(f"\n  writing SLIM sidecar {out.name} ...", flush=True)
    t = time.time()
    slim.write_h5ad(out)
    print(f"    wrote {slim.shape} in {time.time()-t:.1f}s")
    return out


# =============================================================================
# Per-sample driver
# =============================================================================
def process_sample(sample, mod, subpanel_dir, use_gpu, rsc, mode, make_fig,
                   input_path=None):
    print(f"\n{'='*72}\n=== {sample}\n{'='*72}", flush=True)
    p = (Path(input_path) if input_path is not None else
         DATA_ROOT / f"processed/{sample}/{sample}_phenotyped.h5ad")
    if not p.exists():
        raise FileNotFoundError(f"{sample}: required phenotyped input not found: {p}")
    t0 = time.time()
    adata = sc.read_h5ad(p)
    print(f"  loaded {adata.shape} in {(time.time()-t0):.1f}s", flush=True)
    if "lognorm" not in adata.layers:
        raise SystemExit(f"  {sample}: layers['lognorm'] missing — required")
    # Drop any stale full-cell broad UMAP obsm (now stored as a uns subsample).
    adata.obsm.pop("X_umap_broad", None)

    # preserve a prior flat label for provenance, if the base still carries one
    if "celltype" in adata.obs.columns:
        adata.obs["celltype_flat"] = adata.obs["celltype"].astype(str)

    ensure_celltype_panel(adata, mod, subpanel_dir)
    # Fine score panels: MARKER_DICT (fail-loud integrity) THEN overlay identity_core
    # subtype_hint panels where coverage is sufficient (the broad tier is
    # identity_core-only; fine keeps MARKER_DICT only as the sanctioned thin-fallback).
    available = build_available_markers(adata, mod)
    available = overlay_ic_fine_panels(adata, available, mod, subpanel_dir)
    print(f"  fine types scored: {sorted(available)}")

    score_fine(adata, available)                       # tight {ct}_score for the fine tier
    # Tier-3 STATE programs → obs['state_<group>_score'] (orthogonal; QC/report only,
    # NEVER an identity axis — decision #3 GUARD). Non-fatal if the snapshot lacks them.
    try:
        from scripts import state_panels
        state_panels.score_state_panels(adata, getattr(mod, "TISSUE", "pancreas"),
                                        reference_dir=subpanel_dir)
    except Exception as exc:
        print(f"  [state] scoring skipped (non-fatal): {exc!r}")
    tier1_broad(adata, mod, use_gpu, rsc, make_fig)    # broad = identity_core only
    tier2_fine(adata, mod, available, rsc, use_gpu)
    apply_rare_rules(adata, mod)

    # finalize hierarchy columns (after rare rules, so fine includes Gamma/Epsilon)
    adata.obs["celltype_fine"] = adata.obs["celltype"].astype(str)
    adata.obs["celltype_lineage"] = adata.obs["celltype"].astype(str)  # FINE (contract)
    adata.obs["celltype_hier"] = pd.Categorical(
        adata.obs["celltype_broad"].astype(str) + "/" + adata.obs["celltype"].astype(str))

    # Cycling is an orthogonal cross-lineage PROGRAM, not a lineage/celltype: flag
    # cells co-expressing >=2 proliferation markers, and expose a composite
    # `celltype_cycling` (e.g. "Ductal_cycling") WITHOUT polluting the canonical
    # celltype / celltype_lineage (which downstream isin against fine names).
    cyc_markers = getattr(mod, "CYCLING_MARKERS", ["MKI67", "TOP2A", "PCNA"])
    cyc = detect_count(adata, cyc_markers, RAW_DETECT_FLOOR, layer="counts")
    is_cyc = (cyc >= 2) if cyc is not None else np.zeros(adata.n_obs, bool)
    adata.obs["cycling"] = is_cyc
    _ct = adata.obs["celltype"].astype(str).to_numpy()
    adata.obs["celltype_cycling"] = pd.Categorical(
        np.where(is_cyc, np.char.add(_ct.astype(str), "_cycling"), _ct))
    print(f"  cycling (>=2 of {cyc_markers} @ raw counts>={RAW_DETECT_FLOOR}): "
          f"{int(is_cyc.sum()):,} cells ({100*int(is_cyc.sum())/adata.n_obs:.2f}%)")

    validate(adata, mod)
    # Save results FIRST — a figure failure must never block the write.
    write_output(adata, sample, mode, input_path=p)
    if make_fig:
        try:
            render_figure(adata, sample, DATA_ROOT / f"figures/{sample}/02_phenotyping", mod)
        except Exception as exc:
            print(f"  [fig] render failed (non-fatal): {exc!r}")
    print(f"\n  {sample} done in {(time.time()-t0)/60:.1f} min", flush=True)


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--tissue", choices=("pancreas", "thymus"), default="pancreas")
    ap.add_argument("--samples", nargs="+", default=list(DEFAULT_SAMPLES))
    ap.add_argument(
        "--paths", nargs="+", default=None,
        help="explicit *_phenotyped.h5ad paths (overrides --samples; supports "
             "processed/<serial>/<roi>/ nesting)",
    )
    ap.add_argument("--subpanel-dir", default=str(DEFAULT_SUBPANEL_DIR))
    ap.add_argument("--output", choices=("slim", "inplace"), default="slim",
                    help="slim sidecar (default, non-destructive) or in-place "
                         "merge into _phenotyped.h5ad")
    ap.add_argument("--no-fig", action="store_true")
    ap.add_argument("--cpu", action="store_true", help="force CPU (no GPU/scVI-GPU)")
    ap.add_argument("--dl-workers", type=int, default=DL_WORKERS,
                    help=f"scVI dataloader worker processes (default {DL_WORKERS}; "
                         f"0 = single-threaded). Higher feeds the GPU faster.")
    ap.add_argument("--no-sparse-load", action="store_true",
                    help="disable load_sparse_tensor (densify on CPU instead of GPU)")
    args = ap.parse_args(argv)
    globals()["DL_WORKERS"] = args.dl_workers
    globals()["LOAD_SPARSE"] = not args.no_sparse_load

    mod = importlib.import_module(f"scripts.tissue_markers.{args.tissue}")
    if not hasattr(mod, "SUBPANEL_ALLOWLIST"):
        raise SystemExit(f"tissue '{args.tissue}' has no SUBPANEL_ALLOWLIST yet")
    try:
        jobs = resolve_phenotyping_jobs(
            samples=args.samples, paths=args.paths, data_root=DATA_ROOT
        )
    except (FileNotFoundError, ValueError) as exc:
        ap.error(str(exc))
    print(f"hierarchical_phenotype: tissue={args.tissue} "
          f"inputs={[str(path) for _, path in jobs]} "
          f"output={args.output}")
    # Input preflight deliberately precedes CUDA/scVI initialization. A typo in
    # one path must never turn a multi-sample invocation into a partial run.
    use_gpu = (not args.cpu) and setup_torch_tf32()
    rsc = None if args.cpu else try_import_rsc()

    t0 = time.time()
    for sample, path in jobs:
        process_sample(sample, mod, Path(args.subpanel_dir), use_gpu, rsc,
                       args.output, not args.no_fig, input_path=path)
    print(f"\n{'='*72}\n=== DONE in {(time.time()-t0)/60:.1f} min\n{'='*72}")


if __name__ == "__main__":
    main()
