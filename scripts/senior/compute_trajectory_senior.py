#!/usr/bin/env python3
"""
Senior islet trajectory (diffusion pseudotime) + islet Leiden.

BATCH-CORRECTION + ROBUSTNESS FIX (2026-07). An audit of the previously-deployed
trajectory found the pseudotime was a donor/slide (batch) axis, not disease
progression: batch eta^2 = 0.988 vs disease eta^2 = 0.103; within-status donor
eta^2 ND 0.99 / T1D 0.97; the ND<Aab+<T1D ordering FAILED (Aab+>T1D); the root was
a saturated 3-cell islet. Root causes and the fixes applied here:

  1. NO batch correction on the trajectory graph — neighbors were built on the RAW
     X_scVI_mean. FIX: light Harmony(batch=imageid) theta-SWEEP on X_scVI_mean, pick
     the smallest theta whose within-status donor eta^2 (ND,T1D) <= 0.15 while INS
     stays anti-correlated and ND<T1D holds. (Pairs with train_scvi_senior.py, which
     now sets batch_key="imageid" so the latent is corrected at the source.)
  2. Brittle root = single argmax(INS) over ND -> a 3-cell islet. FIX: medoid of the
     top-K ND islets by INS in the trajectory latent, restricted to islets with
     total_cells >= ROOT_MIN_CELLS (Junior recipe, reaggregate_islets.py).
  3. No minimum-cell QC. FIX: root cell floor + a `low_cell_qc` obs flag + a
     sensitivity check on the >=5-cell subset (dataset size unchanged so core<->peri
     pairing and the app's per-islet keying stay intact).
  4. Validation gates were print-only and ran AFTER write. FIX: hard gates that RAISE
     (write *.rejected.h5ad + exit 1) BEFORE overwriting the deployed file; metrics
     persisted to uns['trajectory_validation'] for the app to surface.
  5. "Core + Peri" combined pseudotime now computed (Step-5b port) so the app's
     selector is real, not a silent no-op.

COHORT CAVEAT (surfaced, not "fixed"): Aab+ = a single donor (6521, 399 islets).
Its pseudotime reflects one individual, so Aab+ is EXCLUDED from the hard ordering /
eta^2 gates (reported informationally only).

Islet Leiden (Spatial tab) unchanged: Harmony(theta=2, batch=imageid) on X_scVI_mean
-> separate graph -> leiden at 4 resolutions. Run in `scimap` (or scvi-env) env.
"""
import argparse
import sys
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import scanpy as sc
from scipy.stats import spearmanr, mannwhitneyu

SD = Path(__file__).resolve().parents[2] / "data"
F = SD / "islets_core_senior.h5ad"
PERI = SD / "islets_peri_senior.h5ad"

# --- tunables -----------------------------------------------------------------
THETA_SWEEP = [0.0, 0.1, 0.5, 1.0]   # smallest passing theta wins (less span compression)
N_NEIGHBORS = 30
K_ROOT = 50
ROOT_MIN_CELLS = 20
LOW_CELL_QC = 10                     # obs flag threshold (noisy small-islet marker)
SENSITIVITY_MIN_CELLS = 5           # subset for the pollution sensitivity check
ETA2_GATE = 0.15                    # within-status donor eta^2 ceiling (ND, T1D)
ALPHA_STRUCT = 0.15                 # combined-mode structural block weight
# Broad-lineage structural phenotypes for combined mode; reconciled to whatever is
# actually present in core.uns['phenotype_names'] at runtime (taxonomy may drift).
STRUCT_PHENOTYPES = ["Endothelial", "Fibroblast", "Muscle", "Blood vessel",
                     "Blood Vessel", "Neural", "Lymphatic"]


def eta2(y, groups):
    s = pd.Series(np.asarray(y, dtype=float)); grand = s.mean(); sst = ((s - grand) ** 2).sum()
    ssb = sum(len(v) * (v.mean() - grand) ** 2 for _, v in s.groupby(pd.Series(np.asarray(groups))))
    return float(ssb / sst) if sst > 0 else float("nan")


def within_status_batch_eta2(pt, status, donor, statuses=("ND", "T1D"), min_islets=3):
    """Per-status donor eta^2 using only donors with >= min_islets islets."""
    out = {}
    for stt in statuses:
        m = np.asarray(status) == stt
        pt_s, d_s = np.asarray(pt)[m], np.asarray(donor)[m]
        keep = np.isin(d_s, [d for d in np.unique(d_s) if (d_s == d).sum() >= min_islets])
        if keep.sum() < min_islets or len(np.unique(d_s[keep])) < 2:
            out[stt] = float("nan")
        else:
            out[stt] = eta2(pt_s[keep], d_s[keep])
    return out


def _harmony_rep(a, theta):
    """Return (rep_name) after (optionally) Harmony-correcting X_scVI_mean on imageid.
    theta==0 -> use raw X_scVI_mean. harmonypy Z_corr is (d, N) -> transpose to (N, d)."""
    if theta <= 0:
        return "X_scVI_mean"
    import harmonypy as hm
    a.obs["_iid"] = a.obs["imageid"].astype(str).astype("category")
    ho = hm.run_harmony(np.asarray(a.obsm["X_scVI_mean"], dtype=np.float64),
                        a.obs, "_iid", theta=theta, max_iter_harmony=20, verbose=False)
    Z = np.asarray(ho.Z_corr)
    if Z.shape[0] != a.n_obs:       # harmonypy returns (d, N); we want (N, d)
        Z = Z.T
    a.obsm["X_scVI_harmony_pt"] = np.ascontiguousarray(Z.astype(np.float32))
    del a.obs["_iid"]
    return "X_scVI_harmony_pt"


def _robust_root(a, rep_name, ins_idx):
    """Medoid of the top-K ND islets by INS (in `rep_name`), restricted to islets
    with total_cells >= ROOT_MIN_CELLS. Falls back to unfiltered ND if too few."""
    ds = a.obs["donor_status"].astype(str).to_numpy()
    tc = a.obs["total_cells"].to_numpy()
    nd = np.where((ds == "ND") & (tc >= ROOT_MIN_CELLS))[0]
    if nd.size < 10:
        nd = np.where(ds == "ND")[0]
        if nd.size == 0:
            raise ValueError("No ND islets found for root selection")
    ins = np.asarray(a.X[nd, ins_idx]).ravel()
    top = nd[np.argsort(-ins)[:min(K_ROOT, nd.size)]]
    rep = np.asarray(a.obsm[rep_name])
    centroid = rep[top].mean(axis=0)
    root = int(top[np.argmin(np.linalg.norm(rep[top] - centroid, axis=1))])
    pool = a.obs["imageid"].iloc[top].astype(str).value_counts().to_dict()
    return root, pool


def _clear_traj(a):
    for k in ("X_diffmap",):
        a.obsm.pop(k, None)
    for k in ("diffmap_evals", "neighbors", "paga", "iroot"):
        a.uns.pop(k, None)


def compute_core_dpt(a, theta, ins_idx):
    """Build neighbors on the (Harmony) latent, pick a robust root, run diffmap+DPT.
    Mutates `a`; returns the arrays + root metadata needed to set them as default."""
    _clear_traj(a)
    rep_name = _harmony_rep(a, theta)
    sc.pp.neighbors(a, n_neighbors=N_NEIGHBORS, use_rep=rep_name, metric="cosine")
    root, pool = _robust_root(a, rep_name, ins_idx)
    a.uns["iroot"] = root
    sc.tl.diffmap(a, n_comps=10)
    sc.tl.dpt(a)
    return dict(
        pt=a.obs["dpt_pseudotime"].to_numpy().copy(),
        dc1=a.obsm["X_diffmap"][:, 0].copy(),
        dc2=a.obsm["X_diffmap"][:, 1].copy(),
        diffmap=a.obsm["X_diffmap"].copy(),
        neighbors_uns=a.uns.get("neighbors", None),
        rep_name=rep_name, root=root, pool=pool, theta=theta,
    )


def evaluate_gates(a, pt, ins_idx, root, exclude_aab=True):
    """Return a metrics dict + per-gate booleans + hard_pass. Aab+ excluded from the
    hard ordering / eta^2 gates (single donor -> structurally unreliable)."""
    ds = a.obs["donor_status"].astype(str).to_numpy()
    don = a.obs["imageid"].astype(str).to_numpy()
    tc = a.obs["total_cells"].to_numpy()
    ins = np.asarray(a.X[:, ins_idx]).ravel()
    v = np.isfinite(pt)
    r_ins = float(spearmanr(pt[v], ins[v])[0])
    r_gcg = float("nan")
    if "GCG" in a.var_names:
        gcg = np.asarray(a.X[:, list(a.var_names).index("GCG")]).ravel()
        r_gcg = float(spearmanr(pt[v], gcg[v])[0])
    nd_pt, t1d_pt = pt[ds == "ND"], pt[ds == "T1D"]
    med_nd = float(np.median(nd_pt)) if nd_pt.size else float("nan")
    med_t1d = float(np.median(t1d_pt)) if t1d_pt.size else float("nan")
    try:
        mwu_p = float(mannwhitneyu(t1d_pt, nd_pt, alternative="greater").pvalue)
    except ValueError:
        mwu_p = float("nan")
    wse = within_status_batch_eta2(pt, ds, don, ("ND", "T1D"))
    means = {s: float(pt[ds == s].mean()) for s in ["ND", "Aab+", "T1D"] if (ds == s).any()}

    g_ins = r_ins < -0.3
    g_gcg = (not np.isfinite(r_gcg)) or (r_gcg > 0.1)
    g_order = (med_t1d > med_nd) and (np.isfinite(mwu_p) and mwu_p < 0.05)
    g_eta = all((not np.isfinite(wse[s])) or (wse[s] <= ETA2_GATE) for s in ("ND", "T1D"))
    g_root = (ds[root] == "ND") and (tc[root] >= ROOT_MIN_CELLS)
    hard = bool(g_ins and g_gcg and g_order and g_eta and g_root)

    metrics = dict(
        spearman_INS=r_ins, spearman_GCG=r_gcg,
        median_pt_ND=med_nd, median_pt_T1D=med_t1d, mwu_T1D_gt_ND_p=mwu_p,
        within_status_eta2_ND=wse["ND"], within_status_eta2_T1D=wse["T1D"],
        eta2_disease=eta2(pt, ds), eta2_batch=eta2(pt, don),
        mean_pt_by_status=means, root_index=int(root),
        root_donor=str(don[root]), root_cells=int(tc[root]),
        gate_INS=bool(g_ins), gate_GCG=bool(g_gcg), gate_ordering=bool(g_order),
        gate_within_status_eta2=bool(g_eta), gate_root=bool(g_root), hard_pass=hard,
        aab_excluded_from_hard_gates=bool(exclude_aab),
    )
    return metrics


def _print_gates(tag, m):
    print(f"\n[{tag}] validation")
    print(f"  spearman(pt,INS) = {m['spearman_INS']:+.3f}  (<-0.3) -> {'PASS' if m['gate_INS'] else 'FAIL'}")
    print(f"  spearman(pt,GCG) = {m['spearman_GCG']:+.3f}  (>0.1)  -> {'PASS' if m['gate_GCG'] else 'FAIL'}")
    print(f"  median pt  T1D {m['median_pt_T1D']:.3f} > ND {m['median_pt_ND']:.3f} "
          f"(MWU p={m['mwu_T1D_gt_ND_p']:.1e}) -> {'PASS' if m['gate_ordering'] else 'FAIL'}")
    print(f"  within-status donor eta^2  ND={m['within_status_eta2_ND']:.3f} "
          f"T1D={m['within_status_eta2_T1D']:.3f}  (<= {ETA2_GATE}) -> {'PASS' if m['gate_within_status_eta2'] else 'FAIL'}")
    print(f"  root: donor={m['root_donor']} cells={m['root_cells']} -> {'PASS' if m['gate_root'] else 'FAIL'}")
    print(f"  [info] eta2 disease={m['eta2_disease']:.3f} batch={m['eta2_batch']:.3f} "
          f"means={ {k: round(v,3) for k,v in m['mean_pt_by_status'].items()} }")


def choose_theta(a, ins_idx):
    """Sweep THETA_SWEEP on copies; pick the smallest theta passing the hard gates.
    If none pass, pick the one minimising max within-status eta^2 subject to
    INS<-0.3 & ordering; else the largest theta. Returns (theta, metrics_by_theta)."""
    results = {}
    for theta in THETA_SWEEP:
        try:
            r = compute_core_dpt(a.copy(), theta, ins_idx)
            m = evaluate_gates(a, r["pt"], ins_idx, r["root"])
            results[theta] = m
            _print_gates(f"theta={theta}", m)
        except Exception as e:           # a bad Harmony run shouldn't abort the sweep
            print(f"[theta={theta}] ERROR: {type(e).__name__}: {e}")
    if not results:
        raise RuntimeError("theta sweep produced no results")
    passing = [t for t in THETA_SWEEP if results.get(t, {}).get("hard_pass")]
    if passing:
        chosen = min(passing)
    else:
        ok = [t for t in THETA_SWEEP if t in results
              and results[t]["gate_INS"] and results[t]["gate_ordering"]]
        if ok:
            chosen = min(ok, key=lambda t: max(
                results[t]["within_status_eta2_ND"] if np.isfinite(results[t]["within_status_eta2_ND"]) else 1.0,
                results[t]["within_status_eta2_T1D"] if np.isfinite(results[t]["within_status_eta2_T1D"]) else 1.0))
        else:
            chosen = max(results)
    print(f"\n[theta] chosen = {chosen} "
          f"({'hard gates pass' if results[chosen]['hard_pass'] else 'best available — will be gated at write'})")
    return chosen, results


def sensitivity_check(a, pt, ins_idx):
    """Compare INS rho + within-status eta^2 on the >=SENSITIVITY_MIN_CELLS subset vs
    all islets. A large improvement means tiny islets pollute the manifold."""
    tc = a.obs["total_cells"].to_numpy()
    sub = tc >= SENSITIVITY_MIN_CELLS
    ds = a.obs["donor_status"].astype(str).to_numpy()
    don = a.obs["imageid"].astype(str).to_numpy()
    ins = np.asarray(a.X[:, ins_idx]).ravel()
    full_eta = within_status_batch_eta2(pt, ds, don)
    sub_eta = within_status_batch_eta2(pt[sub], ds[sub], don[sub])
    r_full = float(spearmanr(pt, ins)[0]); r_sub = float(spearmanr(pt[sub], ins[sub])[0])
    drops = [full_eta[s] - sub_eta[s] for s in ("ND", "T1D")
             if np.isfinite(full_eta[s]) and np.isfinite(sub_eta[s])]
    worst_drop = max(drops) if drops else 0.0
    print(f"\n[sensitivity] >= {SENSITIVITY_MIN_CELLS} cells (n={int(sub.sum())}/{a.n_obs}): "
          f"INS rho {r_full:+.3f} -> {r_sub:+.3f}; within-status eta^2 "
          f"drop(max)={worst_drop:.3f}")
    if worst_drop > 0.05:
        print("  NOTE: tiny islets materially inflate batch eta^2 — consider the graph-drop "
              "fallback (build graph on >=5-cell islets, transfer pt by NN in the latent).")
    return dict(subset_n=int(sub.sum()), INS_rho_full=r_full, INS_rho_subset=r_sub,
                eta2_worst_drop=float(worst_drop))


def _struct_block(adata, names_present, idx_map, label):
    props = np.asarray(adata.obsm["phenotype_proportions"])[:, idx_map]
    counts = np.asarray(adata.obsm["phenotype_counts"])[:, idx_map]
    z = lambda x: (x - x.mean(0)) / (x.std(0) + 1e-8)
    print(f"  [combined] {label} structural phenotypes: {names_present}")
    return np.hstack([z(props), z(counts)])


def compute_combined_dpt(a, root_idx, ins_idx):
    """Step-5b: pseudotime on [core scVI(10) | peri scVI(10) | alpha*z(struct)].
    Returns comb_pt/dc1/dc2 or None if inputs are unavailable."""
    if not PERI.exists():
        print(f"[combined] {PERI} missing — skipping combined mode"); return None
    peri = ad.read_h5ad(PERI)
    if "X_scVI_mean" not in a.obsm or "X_scVI_mean" not in peri.obsm:
        print("[combined] X_scVI_mean absent in core/peri — skipping"); return None
    try:
        peri_idx = pd.Series(range(peri.n_obs), index=peri.obs_names).loc[a.obs_names].to_numpy()
    except KeyError:
        print("[combined] cannot align peri to core by obs_names (require_paired) — skipping")
        return None
    core_scvi = np.asarray(a.obsm["X_scVI_mean"])
    peri_scvi = np.asarray(peri.obsm["X_scVI_mean"])[peri_idx]

    parts = [core_scvi, peri_scvi]
    names = list(a.uns.get("phenotype_names", []))
    present = [p for p in STRUCT_PHENOTYPES if p in names]
    # de-dup case variants (e.g. Blood vessel / Blood Vessel) keeping first hit
    seen, present_u = set(), []
    for p in present:
        key = p.lower()
        if key not in seen:
            seen.add(key); present_u.append(p)
    if present_u:
        idx_map = [names.index(p) for p in present_u]
        parts.append(_struct_block(a, present_u, idx_map, "core") * ALPHA_STRUCT)
        if "phenotype_proportions" in peri.obsm and set(present_u).issubset(peri.uns.get("phenotype_names", [])):
            pnames = list(peri.uns["phenotype_names"])
            pidx = [pnames.index(p) for p in present_u]
            parts.append(_struct_block(peri[peri_idx], present_u, pidx, "peri") * ALPHA_STRUCT)
    else:
        print("[combined] no structural phenotypes matched taxonomy — core+peri only")
    X_comb = np.hstack(parts)
    print(f"[combined] matrix {X_comb.shape} (alpha_struct={ALPHA_STRUCT})")

    import harmonypy as hm
    a.obs["_iid_c"] = a.obs["imageid"].astype(str).astype("category")
    ho = hm.run_harmony(X_comb.astype(np.float64), a.obs, "_iid_c",
                        theta=0.1, max_iter_harmony=20, verbose=False)
    Z = np.asarray(ho.Z_corr)
    if Z.shape[0] != a.n_obs:
        Z = Z.T
    a.obsm["X_scVI_combined_harmony_pt"] = np.ascontiguousarray(Z.astype(np.float32))
    del a.obs["_iid_c"]

    _clear_traj(a)
    sc.pp.neighbors(a, n_neighbors=N_NEIGHBORS, use_rep="X_scVI_combined_harmony_pt", metric="cosine")
    a.uns["iroot"] = root_idx
    sc.tl.diffmap(a, n_comps=10)
    sc.tl.dpt(a)
    return dict(pt=a.obs["dpt_pseudotime"].to_numpy().copy(),
                dc1=a.obsm["X_diffmap"][:, 0].copy(),
                dc2=a.obsm["X_diffmap"][:, 1].copy())


def compute_leiden(a):
    try:
        import harmonypy as hm
        meta = a.obs.copy(); meta["_imageid_str"] = pd.Categorical(a.obs["imageid"].astype(str))
        ho = hm.run_harmony(np.asarray(a.obsm["X_scVI_mean"], dtype=np.float64),
                            meta, "_imageid_str", theta=2.0, max_iter_harmony=20)
        Z = np.asarray(ho.Z_corr)
        if Z.shape[0] != a.n_obs:
            Z = Z.T
        a.obsm["X_scVI_harmony"] = np.ascontiguousarray(Z)
        sc.pp.neighbors(a, use_rep="X_scVI_harmony", n_neighbors=15, metric="cosine",
                        key_added="leiden_nn")
        for res in [0.3, 0.5, 0.8, 1.0]:
            sc.tl.leiden(a, resolution=res, key_added=f"leiden_{res}", neighbors_key="leiden_nn")
        print("[leiden] Harmony(theta=2) on X_scVI_mean -> leiden_{0.3,0.5,0.8,1.0}")
        return True
    except Exception as e:
        print("[leiden skipped]:", type(e).__name__, e)
        return False


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--force", action="store_true",
                    help="write even if hard validation gates fail (debug only)")
    args = ap.parse_args()

    a = ad.read_h5ad(F)
    assert "X_scVI_mean" in a.obsm, (
        "X_scVI_mean missing from islets_core_senior.h5ad — re-run aggregate_islets.py "
        "with the scVI latent attached (train_scvi_senior.py with batch_key must run first)."
    )
    ins_idx = list(a.var_names).index("INS")

    # QC flag (dataset size unchanged — preserves core<->peri pairing + app keying).
    a.obs["low_cell_qc"] = (a.obs["total_cells"].to_numpy() < LOW_CELL_QC)

    # --- Core pseudotime: theta-sweep, then materialise the chosen theta on `a` ---
    print("=" * 60 + "\nSTEP 5: Core pseudotime (Harmony theta-sweep)\n" + "=" * 60)
    theta, sweep = choose_theta(a, ins_idx)
    core = compute_core_dpt(a, theta, ins_idx)          # populate real `a` cleanly
    a.obs["dpt_pseudotime"] = core["pt"]
    a.obs["DC1"], a.obs["DC2"] = core["dc1"], core["dc2"]
    core_metrics = evaluate_gates(a, core["pt"], ins_idx, core["root"])
    _print_gates("core (chosen theta)", core_metrics)
    print(f"[root] medoid of top-{K_ROOT} ND-by-INS (>= {ROOT_MIN_CELLS} cells) in {core['rep_name']}: "
          f"idx={core['root']} donor pool={core['pool']}")
    sens = sensitivity_check(a, core["pt"], ins_idx)

    # --- Combined (core+peri+structures) pseudotime — Step 5b ---
    print("\n" + "=" * 60 + "\nSTEP 5b: Combined (Core+Peri) pseudotime\n" + "=" * 60)
    # stash core defaults (combined recompute overwrites neighbors/diffmap on `a`)
    core_diffmap, core_neigh = core["diffmap"], core["neighbors_uns"]
    comb = compute_combined_dpt(a, core["root"], ins_idx)
    combined_metrics = None
    if comb is not None:
        combined_metrics = evaluate_gates(a, comb["pt"], ins_idx, core["root"])
        _print_gates("combined", combined_metrics)
        # restore core as the default axis
        a.obs["dpt_pseudotime"] = core["pt"]; a.obs["DC1"], a.obs["DC2"] = core["dc1"], core["dc2"]
        a.obsm["X_diffmap"] = core_diffmap
        if core_neigh is not None:
            a.uns["neighbors"] = core_neigh
        if combined_metrics["hard_pass"] or args.force:
            a.obs["dpt_pseudotime_combined"] = comb["pt"]
            a.obs["DC1_combined"], a.obs["DC2_combined"] = comb["dc1"], comb["dc2"]
            why = "gates pass" if combined_metrics["hard_pass"] else "--force override (gates FAIL)"
            print(f"[combined] {why} -> writing dpt_pseudotime_combined")
        else:
            print("[combined] gates FAIL -> NOT writing combined column (app falls back to core)")
    else:
        a.obsm["X_diffmap"] = core_diffmap
        if core_neigh is not None:
            a.uns["neighbors"] = core_neigh

    # --- Islet Leiden (Spatial tab) ---
    leiden_ok = compute_leiden(a)

    # --- Display UMAP from raw-marker PCA (clearer than a donor-separated scVI blob;
    #     DPT uses X_diffmap, so this is display-only). Junior recipe. ---
    print("\n[umap] visualization UMAP from raw-marker PCA")
    sc.pp.pca(a, n_comps=15)
    sc.pp.neighbors(a, n_neighbors=30, use_rep="X_pca", metric="euclidean", key_added="neighbors_viz")
    sc.tl.umap(a, neighbors_key="neighbors_viz", min_dist=0.5, spread=2.0)

    # --- Persist validation for the app + final hard gate before write ---
    a.uns["trajectory_validation"] = {
        "core": core_metrics,
        "combined": combined_metrics if combined_metrics is not None else {},
        "chosen_theta": float(theta),
        "theta_sweep": {str(t): sweep[t] for t in sweep},
        "sensitivity": sens,
        "aab_note": "Aab+ is a single donor (6521); excluded from hard ordering/eta^2 gates.",
    }

    if not core_metrics["hard_pass"] and not args.force:
        rej = F.with_suffix(".rejected.h5ad")
        a.write_h5ad(rej)
        print(f"\nHARD GATES FAILED — NOT overwriting {F.name}. Wrote {rej.name} for inspection.")
        print("Re-run with --force to override (debug), or fix the batch correction / inputs.")
        sys.exit(1)

    a.write_h5ad(F)
    print(f"\n[done] {F}  (theta={theta}, leiden={'yes' if leiden_ok else 'no'}, "
          f"combined={'yes' if (combined_metrics and (combined_metrics['hard_pass'] or args.force)) else 'no'})")


if __name__ == "__main__":
    main()
