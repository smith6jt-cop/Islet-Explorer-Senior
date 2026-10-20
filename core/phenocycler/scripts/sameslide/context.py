"""Cross-fitted protein-guided RNA support; reported separately from independent RNA evidence."""
from __future__ import annotations

import numpy as np
import pandas as pd
import scipy.sparse as sp
from scipy.stats import poisson
from evidence import detected as independent_detection, ALPHA

#: A class must have at least this many cells before it gets its own prior. Below it, the
#: pooled prior is used and the cell is recorded as having had no usable context -- a small
#: class would otherwise estimate pi from a handful of cells and state it with confidence.
MIN_CLASS_CELLS = 200

#: Purity below this and the protein call is not a usable context: the Xenium cell straddles
#: several PhenoCycler cells that disagree, so its "class" is an artefact of the overlap.
#: These cells fall back to the pooled prior rather than being dropped -- they still have
#: counts, and a pooled prior is exactly the uniform assumption, which is the honest default.
MIN_CALL_PURITY = 0.8
MIN_OVERLAP = 0.5

#: Floor on the expressed-state mean, so a class where a gene is background-only still has a
#: defined alternative hypothesis rather than a degenerate one.
MU_FLOOR = 0.05

POOLED = "__pooled__"


class ContextError(ValueError):
    pass


def context_labels(calls: pd.DataFrame, *, min_purity: float = MIN_CALL_PURITY,
                   min_overlap: float = MIN_OVERLAP,
                   min_cells: int = MIN_CLASS_CELLS) -> tuple[pd.Series, pd.Series]:
    """Protein call per cell, blanked to POOLED where it is not usable as context.

    Two ways a call fails to be context: the transfer was impure (the Xenium cell overlaps
    PhenoCycler cells that disagree), or the class is too small to estimate a prior from.
    Both fall back to POOLED rather than being dropped.
    """
    for threshold in (min_purity, min_overlap):
        if not 0 <= threshold <= 1:
            raise ContextError("purity and overlap thresholds must be between zero and one")
    lab = calls["call"].fillna(POOLED).astype(str)
    reason = pd.Series("", index=calls.index, dtype=object)
    invalid = lab.str.lower().isin({"", "nan", "none", "ambiguous", "unknown", "unavailable", POOLED})
    invalid |= lab.str.startswith("Unavailable")
    if "status" in calls:
        invalid |= calls["status"].fillna("unavailable").ne("awarded")
    reason.loc[invalid] = "unresolved_protein_call"
    for column, threshold, message in (("call_purity", min_purity, "low_call_purity"),
                                        ("overlap_fraction", min_overlap, "low_overlap")):
        values = pd.to_numeric(calls.get(column, pd.Series(np.nan, index=calls.index)), errors="coerce")
        bad = ~np.isfinite(values) | (values < threshold)
        reason.loc[bad & reason.eq("")] = message
    lab = lab.where(reason.eq(""), POOLED)
    counts = lab.value_counts()
    small = set(counts[counts < min_cells].index) - {POOLED}
    reason.loc[lab.isin(small)] = "small_protein_class"
    return lab.where(~lab.isin(small), POOLED), reason


def class_labels(calls: pd.DataFrame, *, min_purity: float = MIN_CALL_PURITY,
                 min_overlap: float = MIN_OVERLAP,
                 min_cells: int = MIN_CLASS_CELLS) -> pd.Series:
    return context_labels(calls, min_purity=min_purity, min_overlap=min_overlap,
                          min_cells=min_cells)[0]


def background_detection_rate(lam: np.ndarray) -> float:
    """Mean P(at least one background count) over a set of cells."""
    return float(np.mean(1.0 - np.exp(-np.asarray(lam, dtype=float))))


def class_tables(counts: sp.csr_matrix, lam: np.ndarray, labels: pd.Series,
                 *, min_cells: int = MIN_CLASS_CELLS
                 ) -> tuple[pd.DataFrame, pd.DataFrame, dict]:
    """Per (class, gene): the expressing prior pi, and the expressed-state mean mu.

    ``pi`` is the background-CORRECTED detection rate:

        r_obs = pi * P(detect | expressed) + (1 - pi) * r_bg

    with ``P(detect | expressed)`` taken as 1, which makes pi a LOWER bound -- a gene that is
    expressed but often undetected reads as less prevalent than it is. Conservative in the
    right direction: it understates the prior rather than inventing one.

    ``mu`` is the background-corrected mean count, floored so the alternative hypothesis is
    always defined.
    """
    X = sp.csr_matrix(counts)
    if X.shape[0] != len(lam) or X.shape[0] != len(labels):
        raise ContextError(f"counts has {X.shape[0]} rows, lambda {len(lam)}, "
                           f"labels {len(labels)}")
    lam = np.asarray(lam, dtype=float)
    det = independent_detection(X, lam)
    pis, mus, meta = {}, {}, {}
    if not len(lam) or not np.isfinite(lam).all() or (lam < 0).any():
        raise ContextError("background rates require nonempty finite nonnegative training data")
    groups = {POOLED: np.arange(X.shape[0])}
    for cls in labels.unique():
        rows = np.flatnonzero(labels.to_numpy() == cls)
        if cls != POOLED and len(rows) >= min_cells:
            groups[cls] = rows
    for cls, rows in groups.items():
        r_obs = np.asarray(det[rows].mean(axis=0)).ravel()
        # The prior uses the same prespecified detection event as all other RNA stages.
        # Raw count>0 prevalence is retained elsewhere as an audit statistic only.
        critical = poisson.isf(ALPHA, lam[rows]) + 1
        r_bg = float(np.mean(poisson.sf(critical - 1, lam[rows])))
        pi = np.clip((r_obs - r_bg) / max(1.0 - r_bg, 1e-9), 0.0, 1.0)
        mean_obs = np.asarray(X[rows].mean(axis=0)).ravel()
        mu = np.maximum(mean_obs - lam[rows].mean(), MU_FLOOR)
        pis[cls], mus[cls] = pi, mu
        meta[cls] = {"n_cells": int(len(rows)), "r_bg": r_bg,
                     "pi_median": float(np.median(pi)),
                     "genes_pi_zero": int((pi <= 0).sum())}
    return pd.DataFrame(pis).T, pd.DataFrame(mus).T, meta


def detected(counts: sp.csr_matrix, lam: np.ndarray, labels: pd.Series,
             pi: pd.DataFrame, mu: pd.DataFrame, *, odds: float = 1.0) -> sp.csr_matrix:
    """Boolean CSR: is each observed count better explained by expression than by background?

    Evaluated only where ``k > 0`` -- a zero count is not a detection under any rule, and
    restricting to the nonzeros keeps this to the ~215M entries that exist rather than the
    3.6 billion that do not.
    """
    X = sp.csr_matrix(counts).tocoo()
    lam = np.asarray(lam, dtype=float)
    cls_of_row = labels.to_numpy()
    order = {c: i for i, c in enumerate(pi.index)}
    ci = np.array([order[c] for c in cls_of_row], dtype=np.int64)

    k = X.data.astype(np.float64)
    p_i, g_i = ci[X.row], X.col
    prior = pi.to_numpy()[p_i, g_i]
    mu_e = mu.to_numpy()[p_i, g_i]
    lam_c = lam[X.row]

    # Poisson on both sides: background at lambda_c, expression at mu. Both are counts of
    # the same kind, so the likelihood ratio is the only thing that differs between them.
    if odds <= 0:
        raise ContextError("odds must be positive")
    with np.errstate(divide="ignore"):
        ll_expressed = poisson.logpmf(k, mu_e) + np.log(prior)
        ll_background = poisson.logpmf(k, lam_c) + np.log1p(-prior)
    keep = (k > 0) & (ll_expressed > np.log(odds) + ll_background)
    return sp.csr_matrix(
        (np.ones(int(keep.sum()), bool), (X.row[keep], X.col[keep])), shape=X.shape)


def crossfit(counts, lam, calls, xy, *, min_cells=MIN_CLASS_CELLS,
             min_purity=MIN_CALL_PURITY, min_overlap=MIN_OVERLAP):
    """Two spatial folds; every row uses only the other fold's RNA to fit its prior."""
    X = sp.csr_matrix(counts)
    xy = np.asarray(xy, dtype=float)
    if X.shape[0] != len(calls) or len(lam) != len(calls) or xy.shape != (len(calls), 2):
        raise ContextError("counts, background, calls and coordinates must share cell order")
    if len(calls) < 2 or not np.isfinite(xy).all():
        raise ContextError("cross-fitting requires at least two cells with finite coordinates")
    axis = 0 if np.ptp(xy[:, 0]) > 0 else 1
    split = float(np.median(xy[:, axis]))
    folds = (xy[:, axis] >= split).astype(int)
    if len(np.unique(folds)) < 2:
        # The median can equal the minimum when coordinates have many ties.
        split = float((xy[:, axis].min() + xy[:, axis].max()) / 2)
        folds = (xy[:, axis] >= split).astype(int)
    if len(np.unique(folds)) < 2:
        raise ContextError("no spatial split exists for these coordinates")
    labels, reasons = context_labels(calls.reset_index(drop=True), min_cells=0,
                                    min_purity=min_purity, min_overlap=min_overlap)
    resolved = labels.copy()
    matrices, row_indices, training = [], [], {}
    lam = np.asarray(lam)
    for fold in (0, 1):
        test = np.flatnonzero(folds == fold)
        train = np.flatnonzero(folds != fold)
        pi, mu, meta = class_tables(X[train], lam[train], labels.iloc[train].reset_index(drop=True),
                                   min_cells=min_cells)
        lab = labels.iloc[test].reset_index(drop=True)
        missing = ~lab.isin(pi.index)
        resolved.iloc[test[missing.to_numpy()]] = POOLED
        reasons.iloc[test[missing.to_numpy()]] = "class_absent_or_small_in_training_fold"
        lab = lab.where(~missing, POOLED)
        matrices.append(detected(X[test], lam[test], lab, pi, mu).tocoo())
        row_indices.append(test)
        training[str(fold)] = meta
    rows = np.concatenate([idx[m.row] for idx, m in zip(row_indices, matrices)])
    cols = np.concatenate([m.col for m in matrices])
    support = sp.csr_matrix((np.ones(len(rows), bool), (rows, cols)), shape=X.shape)
    per_cell = pd.DataFrame({"context_class": resolved, "context_fallback": reasons,
                             "evaluation_fold": folds})
    return support, per_cell, {"split_axis": int(axis), "split_um": split,
                              "min_class_cells": min_cells, "min_call_purity": min_purity,
                              "min_overlap": min_overlap, "training_classes_by_test_fold": training,
                              "method": "protein_guided_crossfit", "posterior_calibrated": False}


def evaluate(base, calls, out, *, min_cells=MIN_CLASS_CELLS,
             min_purity=MIN_CALL_PURITY, min_overlap=MIN_OVERLAP):
    """Persist support and axes separately; never replace/subtract the raw count matrix."""
    import json
    from pathlib import Path
    import anndata as ad
    from analysis_input import load_base

    base, calls, out = Path(base), Path(calls), Path(out)
    a = ad.read_h5ad(base)
    prep = a.uns.get("prep", {})
    load_base(base, donor=prep.get("donor_id"), roi=prep.get("roi"), serial=prep.get("serial"))
    c = pd.read_parquet(calls)
    if "cell_id" not in c or c.cell_id.astype(str).duplicated().any():
        raise ContextError("protein calls require unique cell_id")
    c["cell_id"] = c.cell_id.astype(str)
    c = c.set_index("cell_id").reindex(a.obs_names)
    fitting = a.var.identity_eligible.to_numpy(bool)
    fitted_support, cells, meta = crossfit(a.layers["counts"][:, fitting], a.obs.gdna_lambda.to_numpy(), c,
                                   a.obsm["spatial"], min_cells=min_cells,
                                   min_purity=min_purity, min_overlap=min_overlap)
    coo = fitted_support.tocoo()
    support = sp.csr_matrix((coo.data, (coo.row, np.flatnonzero(fitting)[coo.col])), shape=a.shape)
    meta['excluded_control_genes'] = a.var_names[~fitting].tolist()
    meta['excluded_control_support_interpretation'] = 'not evaluated by context; see independent fixed-threshold evidence'
    cells.insert(0, "cell_id", a.obs_names.astype(str).to_numpy())
    meta.update({"base": str(base.resolve()), "calls": str(calls.resolve()),
                 "n_cells": a.n_obs, "n_genes": a.n_vars, "n_supported_entries": int(support.nnz),
                 "donor_id": prep.get("donor_id"), "roi": prep.get("roi"), "serial": prep.get("serial"),
                 "fallback_counts": cells.context_fallback.value_counts().to_dict()})
    out.mkdir(parents=True, exist_ok=True)
    sp.save_npz(out / "support.npz", support)
    cells[["cell_id"]].to_csv(out / "cells.csv", index=False)
    pd.DataFrame({"gene": a.var_names.astype(str)}).to_csv(out / "genes.csv", index=False)
    cells.to_parquet(out / "context_cells.parquet", index=False)
    (out / "context.json").write_text(json.dumps(meta, indent=2))
    return meta


def main(argv=None):
    import argparse
    parser = argparse.ArgumentParser(description="Cross-fitted protein-guided RNA support on a common QC base")
    parser.add_argument("--base", required=True)
    parser.add_argument("--calls", required=True)
    parser.add_argument("--out", required=True)
    parser.add_argument("--min-cells", type=int, default=MIN_CLASS_CELLS)
    parser.add_argument("--min-purity", type=float, default=MIN_CALL_PURITY)
    parser.add_argument("--min-overlap", type=float, default=MIN_OVERLAP)
    a = parser.parse_args(argv)
    evaluate(a.base, a.calls, a.out, min_cells=a.min_cells, min_purity=a.min_purity,
             min_overlap=a.min_overlap)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

