#!/usr/bin/env python3
"""Step 4 of the spatial-method benchmark: evaluate every completed method on a
(sample, variant) and emit a scorecard CSV + markdown report + figures.

Runs in `xenium_analysis`. REPORT-ONLY (Tier-C): it tabulates objective metrics
and lets the USER rank — no proportion-plausibility verdict, no "winner" cell.
Mirrors the report-not-gate stance of scripts/validate_hierarchical_pilot.py.

Two tracks:
  * DOMAIN  (banksy/utag/kodama/baseline): coherence + separation + agreement
    + stability + biological recovery.
  * INTERACTION (misty): R2/gain/contributions/importances + a gain-vs-Moran bridge.

Usage:
  python scripts/benchmark/evaluate.py --sample 0041435 --variant crop25k
  python scripts/benchmark/evaluate.py --sample 0041435 --variant crop50k --sweep-variants crop25k crop50k crop100k
"""
from __future__ import annotations

import argparse
import json
import sys
import os
from pathlib import Path

import numpy as np
import pandas as pd

# Two roots, deliberately separate.  ROOT is where the CODE and the tracked tables live;
# DATA_ROOT is where the multi-GB h5ads live and is NEVER inside the repo.  Collapsing them
# -- as the single HiPerGator ROOT did -- would put data/processed/ under analysis/xenium/,
# and the repo-root .gitignore pattern "/data/" is ANCHORED, so it would not be ignored.
ROOT = Path(__file__).resolve().parents[2]                     # analysis/xenium
REPO = Path(__file__).resolve().parents[4]                     # monorepo root
DATA_ROOT = Path(os.environ.get("XENIUM_DATA_ROOT", "/data/xenium"))
for p in (str(ROOT), str(ROOT / "scripts" / "benchmark")):
    if p not in sys.path:
        sys.path.insert(0, p)

import bench_common as bc  # noqa: E402
from eval_lib import (  # noqa: E402
    agreement,
    bio_recovery,
    domain_metrics,
    figures,
    interaction_metrics,
)

DOMAIN_QUALITY_COLS = ["n_domains", "size_gini", "chaos", "pas", "moran_domains",
                       "silhouette_scvi", "davies_bouldin_scvi", "calinski_harabasz_scvi"]


def _load_input(sample, variant):
    import scanpy as sc
    p = bc.input_paths(sample, variant)["h5ad"]
    return sc.read_h5ad(p)


def _aligned_labels(sample, variant, ref_ids):
    """{method: labels aligned to ref_ids} for every domain method whose labels
    file exists (seed 0)."""
    out = {}
    ref = pd.Index(ref_ids.astype(str))
    for m in bc.DOMAIN_METHODS:
        f = bc.labels_path(sample, variant, m)
        if not f.exists():
            continue
        df = pd.read_csv(f, dtype={"cell_id": str})
        s = df.set_index("cell_id")["domain"].reindex(ref)
        if s.isna().any():
            print(f"  [warn] {m}: {int(s.isna().sum())} cells unmatched — skipping")
            continue
        out[m] = s.astype(str).to_numpy()
    return out


def _seeded_labels(sample, variant, ref_ids):
    """{method: {seed: labels}} from results/stability/*_seed*.csv (optional)."""
    ref = pd.Index(ref_ids.astype(str))
    sdir = bc.results_dir(sample) / "stability"
    out: dict = {}
    if not sdir.exists():
        return out
    for f in sdir.glob(f"*_{variant}_seed*_labels.csv"):
        name = f.stem
        method = name.split(f"_{variant}_seed")[0]
        seed = int(name.split("seed")[-1].split("_")[0])
        df = pd.read_csv(f, dtype={"cell_id": str}).set_index("cell_id")["domain"].reindex(ref)
        if df.isna().any():
            continue
        out.setdefault(method, {})[seed] = df.astype(str).to_numpy()
    # include seed-0 main labels so stability has >=2 seeds
    for m in bc.DOMAIN_METHODS:
        f = bc.labels_path(sample, variant, m)
        if f.exists():
            df = pd.read_csv(f, dtype={"cell_id": str}).set_index("cell_id")["domain"].reindex(ref)
            if not df.isna().any():
                out.setdefault(m, {})[0] = df.astype(str).to_numpy()
    return out


def _metrics_json(sample, variant, method):
    f = bc.metrics_path(sample, variant, method)
    if f.exists():
        try:
            return json.loads(f.read_text())
        except Exception:  # noqa: BLE001
            return {}
    return {}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--sample", required=True)
    ap.add_argument("--variant", default=bc.HEADLINE_VARIANT)
    ap.add_argument("--tissue")
    ap.add_argument("--sweep-variants", nargs="*", default=None,
                    help="variants to include in the scalability curves")
    args = ap.parse_args()
    tissue = args.tissue or bc.SAMPLE_TISSUE.get(args.sample, "?")
    sample, variant = args.sample, args.variant
    figd = bc.fig_dir(sample)
    figd.mkdir(parents=True, exist_ok=True)

    print(f"== evaluate {sample} ({tissue}) :: {variant} ==")
    adata = _load_input(sample, variant)
    ids = (adata.obs["cell_id"].astype(str).to_numpy()
           if "cell_id" in adata.obs.columns else adata.obs_names.astype(str).to_numpy())
    xy = np.asarray(adata.obsm["spatial"], dtype=np.float64)
    X_scvi = np.asarray(adata.obsm["X_scvi"]) if "X_scvi" in adata.obsm else None
    knn = domain_metrics.spatial_knn_idx(xy, bc.EVAL_SPATIAL_K)

    label_dict = _aligned_labels(sample, variant, ids)
    print(f"  domain methods present: {list(label_dict)}")

    # ---- domain-track metrics ----
    rows = {}
    for m, lab in label_dict.items():
        print(f"  metrics: {m}")
        rows[m] = domain_metrics.compute_domain_metrics(lab, xy, X_scvi, knn, seed=bc.SEED)
    metrics_df = pd.DataFrame(rows).T if rows else pd.DataFrame()

    # ---- agreement ----
    mats = agreement.pairwise_matrices(label_dict) if len(label_dict) >= 2 else {}
    consensus = (agreement.mean_offdiag_agreement(mats["ami"]) if mats else pd.Series(dtype=float))
    for nm, mat in mats.items():
        mat.to_csv(bc.results_dir(sample) / f"agreement_{nm}_{variant}.csv")

    # ---- stability (optional) ----
    stab = agreement.stability(_seeded_labels(sample, variant, ids))

    # ---- biological recovery ----
    bio = {}
    ct_fine = adata.obs["celltype"].astype(str).to_numpy() if "celltype" in adata.obs else None
    ct_broad = (adata.obs["celltype_broad"].astype(str).to_numpy()
                if "celltype_broad" in adata.obs else ct_fine)
    if ct_fine is not None:
        for m, lab in label_dict.items():
            _, _, pear = bio_recovery.enrichment(lab, ct_fine)
            pear.to_csv(bc.results_dir(sample) / f"enrichment_pearson_{m}_{variant}.csv")
            figures.enrichment_heatmap(pear, figd / f"domain_celltype_enrichment_{m}.png", method=m)
            bio[m] = bio_recovery.label_nmi(lab, ct_broad)

    # ---- interaction track (MISTy) ----
    misty = {}
    mcsvs = interaction_metrics.load_misty_csvs(bc.results_dir(sample), variant)
    if mcsvs:
        print("  interaction: MISTy results found")
        improv = interaction_metrics.summarize_improvements(mcsvs)
        contrib = interaction_metrics.summarize_contributions(mcsvs)
        fp = interaction_metrics.fingerprint(improv, contrib)
        misty = {"improv": improv, "contrib": contrib, "fingerprint": fp}
        # bridge: per-marker gain.R2 vs per-marker Moran's I (on lognorm)
        bridge = {"available": False}
        if improv.get("available") and "lognorm" in adata.layers:
            markers = [g for g in improv["_gain"].index if g in set(map(str, adata.var_names))]
            if markers:
                sub = adata[:, markers]
                Xm = sub.layers["lognorm"]
                Xm = Xm.toarray() if hasattr(Xm, "toarray") else np.asarray(Xm)
                moran = dict(zip(markers, interaction_metrics.moran_per_feature(Xm, knn)))
                bridge = interaction_metrics.bridge_gain_vs_moran(improv["_gain"], moran)
        misty["bridge"] = bridge
        # figures
        figures.misty_improvement(improv.get("_gain"), figd / "misty_improvement.png")
        figures.misty_contributions(contrib.get("by_view", {}), figd / "misty_view_contributions.png")
        for vw in ("juxta", "para"):
            mm = interaction_metrics.importances_matrix(mcsvs, view_filter=vw)
            figures.misty_importance_heatmap(mm, figd / f"misty_importance_{vw}.png", view=vw)
        if bridge.get("available"):
            figures.misty_bridge_scatter(bridge["_table"], bridge["spearman_rho"],
                                         figd / "misty_vs_moran_bridge.png")

    # ---- computational + scalability ----
    comp_rows = {}
    for m in bc.ALL_METHODS:
        j = _metrics_json(sample, variant, m)
        if j:
            comp_rows[m] = {k: j.get(k) for k in ("n_cells", "wall_s", "peak_rss_mb",
                                                  "status", "exit_code")}
    comp_df = pd.DataFrame(comp_rows).T if comp_rows else pd.DataFrame()

    sweep_long = []
    for v in (args.sweep_variants or [variant]):
        for m in bc.ALL_METHODS:
            j = _metrics_json(sample, v, m)
            if j and j.get("status") == "ok":
                sweep_long.append({"method": m, "variant": v, "n_cells": j.get("n_cells"),
                                   "wall_s": j.get("wall_s"), "peak_rss_mb": j.get("peak_rss_mb")})
    sweep_df = pd.DataFrame(sweep_long)

    # ---- assemble master scorecard ----
    score = _assemble_scorecard(label_dict, metrics_df, consensus, stab, bio,
                                comp_df, misty)
    out_csv = bc.bench_dir(sample) / f"{sample}_benchmark_scorecard_{variant}.csv"
    score.to_csv(out_csv)
    print(f"  scorecard -> {out_csv}")

    # ---- figures ----
    if label_dict:
        figures.spatial_domain_maps(xy, label_dict, figd / "spatial_domains_all_methods.png",
                                    title=f"{sample} ({tissue}) {variant}: spatial domains")
        figures.domain_quality_bars(metrics_df, ["chaos", "pas", "moran_domains",
                                                 "silhouette_scvi"], figd / "domain_quality_bars.png")
    if mats:
        figures.agreement_heatmap(mats["ami"], figd / "crossmethod_AMI_heatmap.png",
                                  title=f"{sample} {variant}: cross-method AMI")
        figures.agreement_heatmap(mats["ari"], figd / "crossmethod_ARI_heatmap.png",
                                  title=f"{sample} {variant}: cross-method ARI")
    figures.stability_bar(stab, figd / "stability.png")
    if not sweep_df.empty:
        figures.scalability_curves(sweep_df, figd / "scalability_runtime_vs_N.png")

    # ---- report ----
    _write_report(sample, tissue, variant, score, mats, misty, sweep_df, figd, adata)
    print("== done ==")


def _assemble_scorecard(label_dict, metrics_df, consensus, stab, bio, comp_df, misty):
    methods = list(dict.fromkeys(list(label_dict) + (["misty"] if misty else [])))
    rec = {}
    for m in methods:
        r = {"kind": "interaction" if m == "misty" else "domain"}
        for c in ("n_cells", "wall_s", "peak_rss_mb", "status"):
            r[c] = comp_df.loc[m, c] if (not comp_df.empty and m in comp_df.index) else np.nan
        if m in metrics_df.index:
            for c in DOMAIN_QUALITY_COLS:
                r[c] = metrics_df.loc[m, c]
        r["consensus_ami"] = float(consensus[m]) if m in getattr(consensus, "index", []) else np.nan
        r["stability_mean_ari"] = stab.get(m, {}).get("mean_ari", np.nan)
        r["bio_nmi_broad"] = bio.get(m, {}).get("nmi", np.nan)
        r["bio_ami_broad"] = bio.get(m, {}).get("ami", np.nan)
        if m == "misty":
            fp = misty.get("fingerprint", {})
            r["misty_median_gain_R2"] = fp.get("median_gain_R2", np.nan)
            r["misty_n_gain_gt_0.05"] = fp.get("n_gain_gt_0.05", np.nan)
            r["misty_mean_para_fraction"] = fp.get("mean_para_fraction", np.nan)
            r["misty_bridge_rho"] = misty.get("bridge", {}).get("spearman_rho", np.nan)
        rec[m] = r
    return pd.DataFrame(rec).T


def _write_report(sample, tissue, variant, score, mats, misty, sweep_df, figd, adata):
    L = [f"# Spatial-method benchmark — {sample} ({tissue}), variant `{variant}`", "",
         "REPORT-ONLY: objective metrics are tabulated; the **user** ranks the "
         "methods. No proportion-plausibility judgement is made (Tier-C policy). "
         "Each domain-quality metric has a known bias — read them together, and "
         "always with `n_domains`.", "",
         "## Master scorecard", "",
         score.to_markdown(floatfmt=".3f"), "",
         "### Domain-quality metric bias notes", ""]
    for k, v in domain_metrics.METRIC_BIAS.items():
        L.append(f"- **{k}**: {v}")
    L += ["", "## Cross-method agreement (AMI — k-robust; ARI also written)", ""]
    if mats:
        L.append(mats["ami"].to_markdown(floatfmt=".3f"))
    else:
        L.append("_fewer than 2 domain methods present_")
    if misty:
        fp = misty.get("fingerprint", {})
        br = misty.get("bridge", {})
        im = misty.get("improv", {})
        L += ["", "## Interaction track (MISTy)", "",
              f"- targets modelled: {im.get('n_targets', 'NA')}",
              f"- median gain.R2: {fp.get('median_gain_R2', float('nan')):.3f}",
              f"- targets with gain.R2 > 0.05: {fp.get('n_gain_gt_0.05', 'NA')}",
              f"- mean para-view fraction: {fp.get('mean_para_fraction', float('nan')):.3f}",
              f"- bridge Spearman(gain.R2, Moran's I): "
              f"{br.get('spearman_rho', float('nan')):.3f} "
              f"(p={br.get('p_value', float('nan')):.3g}, n={br.get('n', 0)})"]
    if not sweep_df.empty:
        L += ["", "## Scalability (status=ok runs)", "", sweep_df.to_markdown(index=False)]
    L += ["", "## Figures", ""]
    for p in sorted(figd.glob("*.png")):
        L.append(f"- `{p.relative_to(ROOT)}`")
    out = bc.bench_dir(sample) / f"{sample}_benchmark_report_{variant}.md"
    out.write_text("\n".join(L) + "\n")
    print(f"  report -> {out}")


if __name__ == "__main__":
    main()
