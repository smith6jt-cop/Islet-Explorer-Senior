"""Rule-based gamma-delta (γδ) T-cell identification for the thymus Xenium runs.

WHY RULE-BASED (not score_genes argmax)
----------------------------------------
γδ T cells are a *rare* lineage and the discriminating probes are *sparse*
on this panel (per panel_audit_thymus.csv, detection across the two donors):
    TRGC2  0.89% / 1.58%      (γ constant — custom add-on, NOT on base 5K)
    TRDC   0.55% / 0.69%      (δ constant — base 5K panel)
    TRAC   31.4% / 40.7%      (α constant — add-on; αβ contrast)
    TRBC2  27.4% / 36.5%      (β constant — add-on; αβ contrast)
    TRAV24 0.33% / 0.52%      (iNKT invariant Vα — add-on; confound exclusion)
    TRAV1-2 0.40% / 0.88%     (MAIT invariant Vα — add-on; confound exclusion)
Score_genes argmax over a 1-2 gene "panel" against detection asymmetry is
exactly the failure mode CLAUDE.md decisions #12/#15 forbid. Instead we use an
explicit multi-criteria co-detection rule, identical in spirit to the Gamma /
Epsilon rules in scripts/refine_rare_celltypes.py::apply_rare_celltype_rules.

THE CALL (two tiers + optional cluster corroboration)
-----------------------------------------------------
Markers (only ones present on this panel are used; CD3D is NOT on the panel):
    GD_CHAINS       = TRGC2, TRDC            (γ / δ constant)
    CD3_COMPLEX     = CD3E, CD3G, CD247      (bona-fide T-cell evidence)
    AB_CHAINS       = TRAC, TRBC2            (αβ TCR — must be LOW for γδ)
    INNATE_INVARIANT= TRAV24, TRAV1-2        (iNKT / MAIT — must be ABSENT)

For each cell we count "depth" = number of distinct markers in a group detected
at lognorm ≥ a presence threshold (sample-comparable, like decision #11):
    gd_depth, cd3_depth, ab_depth ; innate = any INNATE_INVARIANT detected

Because the TCR chains are ultra-sparse, the standard log1p(3) "depth" gate
(decision #11) would gut sensitivity — a true γδ cell usually shows only ONE
transcript of ONE chain. We therefore relax the presence gate for the sparse
TCR locus to DETECT_TCR = log1p(1) (any real count), and keep the abundant
CD3 context at DETECT_CD3 = log1p(2). Both are module constants — tune freely.

    Tier 1  CONFIDENT  : T-lineage & gd_depth≥1 & cd3_depth≥1 & ab_depth==0
                         & not innate                       → "T_gd"
             (gd_depth==2, both chains, flagged confident_both_chains)
    Tier 2  CANDIDATE  : T-lineage & gd_depth≥1 & cd3_depth≥1 & not innate
                         & gd_depth ≥ ab_depth (γδ evidence ≥ αβ evidence)
                                                             → "T_gd_candidate"
    Pass C  CLUSTER    : (optional, needs a cluster column) a cd3+ / ab-low cell
             that DROPPED its own γδ transcript but sits in a cluster whose
             γδ-positive fraction is ≥ CLUSTER_ENRICH_FOLD× the global rate
             is upgraded                                     → "T_gd_cluster"

celltype is NEVER overwritten (decision #11). The γδ call lands in
`celltype_lineage` (created from `celltype` if absent) plus audit columns.

Usage
-----
    python3 scripts/refine_gd_tcells.py                       # dry-run, all thymus h5ads
    python3 scripts/refine_gd_tcells.py PATH.h5ad [...]       # explicit files
    python3 scripts/refine_gd_tcells.py --cluster-key leiden_scvi --write
"""
from __future__ import annotations

import argparse
import glob
import time
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd

# ---------------------------------------------------------------- parameters
DETECT_TCR = float(np.log1p(1.0))   # ≈ 0.693 — sparse TCR-locus presence gate
DETECT_CD3 = float(np.log1p(2.0))   # ≈ 1.099 — CD3 context (more abundant)
CLUSTER_ENRICH_FOLD = 3.0           # cluster γδ+ fraction vs global to corroborate
CLUSTER_MIN_GD = 3                  # min γδ+ cells in a cluster to trust enrichment

GD_CHAINS = ("TRGC2", "TRDC")
CD3_COMPLEX = ("CD3E", "CD3G", "CD247")
AB_CHAINS = ("TRAC", "TRBC2")
INNATE_INVARIANT = ("TRAV24", "TRAV1-2")

# Existing-label prefixes that count a cell as already T-lineage. If none of
# these obs columns exist, the T-lineage gate falls back to cd3_depth ≥ 1.
T_LINEAGE_COLS = ("celltype_lineage", "celltype", "immune_subtype")
T_LINEAGE_PREFIXES = ("T_", "T cell", "Tcell", "T-cell", "Thymocyte", "DP_", "DN_")

DEFAULT_GLOB = "data/processed/*THYMUS*/*phenotyped.h5ad"
FALLBACK_GLOB = "data/processed/*THYMUS*/*preprocessed.h5ad"


# Names stage 01 / other pipelines have used for the log-normalised layer.
LOGNORM_LAYERS = ("lognorm", "log_normalized", "lognormalized", "log1p")


def lognorm_source(adata):
    """Return (matrix, name) for the log-normalised layer, or (.X, None)."""
    for name in LOGNORM_LAYERS:
        if name in adata.layers:
            return adata.layers[name], name
    return adata.X, None


# ----------------------------------------------------------------- helpers
def get_lognorm_col(adata, gene):
    """(n_cells,) lognorm vector for a gene, or None if the gene is absent."""
    if gene not in adata.var_names:
        return None
    i = adata.var_names.get_loc(gene)
    L = lognorm_source(adata)[0]
    col = L[:, i]
    if hasattr(col, "toarray"):
        col = col.toarray()
    return np.asarray(col).ravel()


def group_depth(adata, genes, thresh):
    """Per-cell count of distinct `genes` detected at lognorm ≥ thresh, plus
    the list of genes actually found on the panel."""
    n = adata.n_obs
    depth = np.zeros(n, dtype=np.int16)
    present = []
    for g in genes:
        v = get_lognorm_col(adata, g)
        if v is None:
            continue
        present.append(g)
        depth += (v >= thresh).astype(np.int16)
    return depth, present


def t_lineage_mask(adata, cd3_depth):
    """True where the cell is already labelled T-lineage by an upstream column,
    OR (fallback) expresses the CD3 complex."""
    n = adata.n_obs
    mask = np.zeros(n, dtype=bool)
    used = None
    for col in T_LINEAGE_COLS:
        if col in adata.obs.columns:
            vals = adata.obs[col].astype(str).values
            hit = np.zeros(n, dtype=bool)
            for pref in T_LINEAGE_PREFIXES:
                hit |= pd.Series(vals).str.startswith(pref).values
            mask |= hit
            used = col if used is None else used
    # Always allow CD3+ cells through (rescues mislabelled / unphenotyped data)
    mask |= cd3_depth >= 1
    return mask, used


# -------------------------------------------------------------- the rule
def apply_gd_tcell_rules(adata, cluster_key=None):
    print(f"  γδ rule  (DETECT_TCR={DETECT_TCR:.3f}, DETECT_CD3={DETECT_CD3:.3f})")

    gd_depth, gd_present = group_depth(adata, GD_CHAINS, DETECT_TCR)
    cd3_depth, cd3_present = group_depth(adata, CD3_COMPLEX, DETECT_CD3)
    ab_depth, ab_present = group_depth(adata, AB_CHAINS, DETECT_TCR)
    inv_depth, inv_present = group_depth(adata, INNATE_INVARIANT, DETECT_TCR)
    innate = inv_depth >= 1

    if not gd_present:
        print("    no γδ-chain genes (TRGC2/TRDC) on panel — nothing to do")
        return 0
    print(f"    panel: gd={gd_present} cd3={cd3_present} "
          f"ab={ab_present} invariant={inv_present}")

    adata.obs["gd_chain_depth"] = gd_depth
    adata.obs["cd3_complex_depth"] = cd3_depth
    adata.obs["ab_chain_depth"] = ab_depth
    adata.obs["gd_innate_invariant"] = innate

    t_mask, used_col = t_lineage_mask(adata, cd3_depth)
    print(f"    T-lineage gate via "
          f"{used_col or '(no label col)'} + CD3 fallback → "
          f"{int(t_mask.sum()):,} candidate T cells")

    base = t_mask & (gd_depth >= 1) & (cd3_depth >= 1) & (~innate)
    confident = base & (ab_depth == 0)
    candidate = base & (~confident) & (gd_depth >= ab_depth)
    both_chains = confident & (gd_depth == 2)

    call = np.array(["none"] * adata.n_obs, dtype=object)
    resolution = np.array([""] * adata.n_obs, dtype=object)
    call[candidate] = "T_gd_candidate"
    resolution[candidate] = "permissive"
    call[confident] = "T_gd"
    resolution[confident] = "confident"
    resolution[both_chains] = "confident_both_chains"

    # ---- optional Pass C: cluster corroboration (recover dropped chains) ----
    if cluster_key and cluster_key in adata.obs.columns:
        clusters = adata.obs[cluster_key].astype(str).values
        gd_pos = (gd_depth >= 1)
        global_rate = gd_pos.mean()
        df = pd.DataFrame({"c": clusters, "gd": gd_pos})
        by = df.groupby("c")["gd"].agg(["mean", "sum"])
        enriched = set(by.index[(by["mean"] >= CLUSTER_ENRICH_FOLD * global_rate)
                                & (by["sum"] >= CLUSTER_MIN_GD)])
        # cd3+ / ab-low T cells in an enriched cluster, not already called,
        # that simply lost their own γδ transcript
        in_enriched = np.isin(clusters, list(enriched))
        upgrade = (t_mask & in_enriched & (call == "none")
                   & (cd3_depth >= 1) & (ab_depth == 0) & (~innate))
        call[upgrade] = "T_gd_cluster"
        resolution[upgrade] = "cluster_corroborated"
        print(f"    cluster pass: global γδ+ rate {global_rate:.4f}; "
              f"{len(enriched)} enriched clusters; "
              f"{int(upgrade.sum()):,} cells upgraded")
    elif cluster_key:
        print(f"    cluster pass skipped — '{cluster_key}' not in obs")

    adata.obs["gd_call"] = pd.Categorical(call)
    adata.obs["gd_resolution"] = pd.Categorical(resolution)

    # Write the lineage label — celltype is NEVER overwritten (decision #11)
    if "celltype_lineage" not in adata.obs.columns:
        src = "celltype" if "celltype" in adata.obs.columns else None
        adata.obs["celltype_lineage"] = (
            adata.obs[src].astype(str).values if src
            else np.array(["unknown"] * adata.n_obs, dtype=object))
    cl = adata.obs["celltype_lineage"].astype(str).values.copy()
    is_gd = np.isin(call, ["T_gd", "T_gd_cluster"])  # candidates NOT relabelled
    cl[is_gd] = "T_gd"
    adata.obs["celltype_lineage"] = pd.Categorical(cl)

    n_conf = int(confident.sum())
    n_both = int(both_chains.sum())
    n_cand = int(candidate.sum())
    n_clu = int((call == "T_gd_cluster").sum())
    print(f"    → T_gd confident: {n_conf:,} ({n_both:,} both-chain) | "
          f"candidate: {n_cand:,} | cluster: {n_clu:,}")
    return n_conf + n_clu


# -------------------------------------------------------------- driver
def process_file(path, cluster_key, write):
    print(f"\n{'='*72}\n=== {path}\n{'='*72}")
    t = time.time()
    a = ad.read_h5ad(path)
    print(f"  loaded {a.shape} in {time.time()-t:.1f}s")
    src_name = lognorm_source(a)[1]
    if src_name:
        print(f"  expression source: layers['{src_name}']")
    else:
        print("  WARN: no known log-normalised layer "
              f"({LOGNORM_LAYERS}) — using .X (ensure it is log-normalised)")
    apply_gd_tcell_rules(a, cluster_key=cluster_key)
    print("\n  gd_call tally:")
    print(a.obs["gd_call"].value_counts().to_string())
    if not write:
        print("\n  --dry-run (default): not writing. Pass --write to persist.")
        return
    print(f"\n  writing {Path(path).name}...")
    t = time.time()
    a.write(path, compression="gzip")
    print(f"    saved in {time.time()-t:.1f}s")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("paths", nargs="*", help="h5ad files (default: thymus glob)")
    ap.add_argument("--cluster-key", default=None,
                    help="obs column for Pass-C cluster corroboration "
                         "(e.g. leiden_scvi)")
    ap.add_argument("--write", action="store_true",
                    help="persist results (default is dry-run)")
    args = ap.parse_args()

    paths = args.paths
    if not paths:
        paths = sorted(glob.glob(DEFAULT_GLOB)) or sorted(glob.glob(FALLBACK_GLOB))
    if not paths:
        raise SystemExit("no h5ad files found — pass paths explicitly")
    for p in paths:
        process_file(p, args.cluster_key, args.write)


if __name__ == "__main__":
    main()
