"""Interaction-track metrics for MISTy (NOT a domain method — evaluated on its
own terms). Reads the tidy CSVs written by run_misty.R and summarises:
  * per-target intra.R2 / multi.R2 / gain.R2 (variance explained, multiview gain)
  * per-view contribution fractions (intra / juxta / para)
  * predictor->target importances (kept for the heatmap)
  * 3-number fingerprint (median gain.R2, n(gain>.05), mean para fraction)

BRIDGE to a standard spatial statistic: Spearman( per-marker gain.R2 ,
per-marker Moran's I ). Markers MISTy finds predictable from spatial views
should be the spatially-autocorrelated markers — a unifying sanity link between
the interaction track and the domain track's spatial graph. (MISTy here models
markers, so a celltype-level nhood_enrichment bridge would need a separate
abundance-mode run; noted as a future extension.)
"""
from __future__ import annotations

import glob
import re
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import sparse
from scipy.stats import spearmanr


def load_misty_csvs(results_dir: Path, variant: str) -> dict:
    out = {}
    for f in glob.glob(str(results_dir / f"misty_{variant}_*.csv")):
        name = re.sub(rf"^misty_{re.escape(variant)}_", "", Path(f).stem)
        try:
            out[name] = pd.read_csv(f)
        except Exception as e:  # noqa: BLE001
            print(f"    [warn] could not read {f}: {e}")
    return out


def _find(dfs: dict, *needles):
    """First dataframe whose key contains all needles (case-insensitive)."""
    for k, v in dfs.items():
        kl = k.lower()
        if all(n in kl for n in needles):
            return v
    return None


def _measure_col(df):
    for c in df.columns:
        if c.lower() in ("measure", "metric"):
            return c
    return None


def summarize_improvements(dfs: dict) -> dict:
    """median/IQR of multi.R2 & gain.R2 + per-target wide table."""
    imp = _find(dfs, "improvements")
    # prefer the raw per-target 'improvements' over the '.stats' aggregate
    for k, v in dfs.items():
        if k.lower() == "improvements":
            imp = v
            break
    if imp is None:
        return {"available": False}
    mcol = _measure_col(imp)
    tcol = next((c for c in imp.columns if c.lower() in ("target", "predictor_target")), None)
    vcol = next((c for c in imp.columns if c.lower() in ("value", "mean")), None)
    if not (mcol and tcol and vcol):
        return {"available": False, "note": f"unexpected cols {list(imp.columns)}"}
    wide = imp.pivot_table(index=tcol, columns=mcol, values=vcol, aggfunc="mean")
    wide.columns = [str(c) for c in wide.columns]

    def col(name):
        # exact measure match first (avoids 'gain' also matching 'p.gain'),
        # then a contains-fallback on the full measure token.
        for c in wide.columns:
            if c.lower() == name:
                return wide[c]
        for c in wide.columns:
            if name in c.lower():
                return wide[c]
        return pd.Series(dtype=float)

    gain = col("gain.r2").dropna()
    multi = col("multi.r2").dropna()
    res = {
        "available": True,
        "n_targets": int(wide.shape[0]),
        "median_gain_R2": float(gain.median()) if len(gain) else float("nan"),
        "iqr_gain_R2": float(gain.quantile(.75) - gain.quantile(.25)) if len(gain) else float("nan"),
        "median_multi_R2": float(multi.median()) if len(multi) else float("nan"),
        "n_gain_gt_0.01": int((gain > 0.01).sum()),
        "n_gain_gt_0.05": int((gain > 0.05).sum()),
        "n_gain_gt_0.10": int((gain > 0.10).sum()),
        "_wide": wide, "_gain": gain,
    }
    return res


def summarize_contributions(dfs: dict) -> dict:
    con = None
    for k, v in dfs.items():
        if k.lower() == "contributions":
            con = v
            break
    con = con if con is not None else _find(dfs, "contributions")
    if con is None:
        return {"available": False}
    vcol = next((c for c in con.columns if c.lower() in ("value", "mean", "fraction")), None)
    vview = next((c for c in con.columns if c.lower() == "view"), None)
    if not (vcol and vview):
        return {"available": False}
    # contributions can be signed/raw; report mean per view normalised to fractions
    g = con.groupby(vview)[vcol].mean()
    g = g.clip(lower=0)
    frac = (g / g.sum()) if g.sum() else g
    return {"available": True, "by_view": frac.to_dict()}


def importances_matrix(dfs: dict, view_filter: str | None = None) -> pd.DataFrame | None:
    agg = _find(dfs, "importances", "aggregated")
    if agg is None:
        agg = _find(dfs, "importances")
    if agg is None:
        return None
    cols = {c.lower(): c for c in agg.columns}
    pcol = cols.get("predictor")
    tcol = cols.get("target")
    icol = cols.get("importance") or cols.get("value") or cols.get("mean")
    vcol = cols.get("view")
    if not (pcol and tcol and icol):
        return None
    d = agg
    if view_filter and vcol:
        d = d[d[vcol].astype(str).str.contains(view_filter, case=False)]
    return d.pivot_table(index=pcol, columns=tcol, values=icol, aggfunc="mean")


def moran_per_feature(X: np.ndarray, knn_idx: np.ndarray) -> np.ndarray:
    """Moran's I for every column of X on the kNN graph."""
    X = np.asarray(X, dtype=float)
    n, p = X.shape
    rows = np.repeat(np.arange(n), knn_idx.shape[1])
    cols = knn_idx.ravel()
    W = sparse.csr_matrix((np.ones(len(rows)), (rows, cols)), shape=(n, n))
    S0 = W.sum()
    Z = X - X.mean(axis=0, keepdims=True)
    num = (Z * (W @ Z)).sum(axis=0)
    den = (Z * Z).sum(axis=0)
    with np.errstate(divide="ignore", invalid="ignore"):
        I = (n / S0) * (num / den)
    return np.asarray(I)


def bridge_gain_vs_moran(gain: pd.Series, moran: dict) -> dict:
    """Spearman( per-marker gain.R2 , per-marker Moran's I )."""
    if gain is None or len(gain) == 0:
        return {"available": False}
    common = [g for g in gain.index if g in moran]
    if len(common) < 3:
        return {"available": False, "n": len(common)}
    g = np.array([gain[c] for c in common])
    m = np.array([moran[c] for c in common])
    rho, p = spearmanr(g, m)
    return {"available": True, "n": len(common), "spearman_rho": float(rho),
            "p_value": float(p),
            "_table": pd.DataFrame({"marker": common, "gain_R2": g, "moran_I": m})}


def fingerprint(improv: dict, contrib: dict) -> dict:
    para = float("nan")
    if contrib.get("available"):
        for v, f in contrib["by_view"].items():
            if "para" in str(v).lower():
                para = float(f)
    return {"median_gain_R2": improv.get("median_gain_R2", float("nan")),
            "n_gain_gt_0.05": improv.get("n_gain_gt_0.05", 0),
            "mean_para_fraction": para}
