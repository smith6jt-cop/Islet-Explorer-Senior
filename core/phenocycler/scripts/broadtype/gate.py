"""The label-free gate.  Every statistic here is a measurement against a null drawn from the
same cells; none consults an expected frequency, a QuPath class, or another arm's answer.

Blocking (per arm, per donor):
  * held-out marker enrichment  -- labels with >= 2 declared markers, one ablation per marker
  * spatial coherence           -- every declared label
  * stability                   -- verdict invariance across replicates (stochastic arms)
Descriptive (reported beside, never blocking):
  * cross-lineage tail contamination
  * QuPath islet-core enrichment of Endocrine (a CHECK that owes nothing to the markers)
  * cross-arm agreement (ARI / AMI)
Unjudgeable is never a pass (``gates.require_evidence`` idiom).
"""
from __future__ import annotations

import pathlib
import sys

import numpy as np
import pandas as pd

HERE = pathlib.Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import pins  # noqa: E402,F401
from criteria import CRITERIA  # noqa: E402
from records import AMBIGUOUS, UNAVAILABLE  # noqa: E402
import gates  # noqa: E402  (scripts/common/gates.py)
import usability  # noqa: E402  (scripts/celltype/usability.py)


# --------------------------------------------------------------------------- tails
def tail_mask(values: np.ndarray, q: float = CRITERIA["tail_q"], min_k: int = CRITERIA["tail_min_k"]) -> np.ndarray:
    """Step 1's top-``q`` tail (k = max(min_k, q*n) over finite values); NaN never in the tail."""
    v = np.asarray(values, dtype=float)
    finite = np.isfinite(v)
    n = int(finite.sum())
    out = np.zeros(v.shape[0], dtype=bool)
    if n == 0:
        return out
    k = min(n, max(min_k, int(q * n)))
    vv = np.where(finite, v, -np.inf)
    idx = np.argpartition(-vv, k - 1)[:k]
    out[idx] = True
    return out


def hash_subsample(ids: np.ndarray, n: int, salt: str = "broadtype") -> np.ndarray:
    """Deterministic, content-keyed subsample (indices), independent of row order."""
    import hashlib
    if len(ids) <= n:
        return np.arange(len(ids))
    keys = np.fromiter((int(hashlib.sha256(f"{i}|{salt}".encode()).hexdigest()[:16], 16) for i in ids),
                       dtype=np.uint64, count=len(ids))
    return np.sort(np.argsort(keys)[:n])


# --------------------------------------------------------------------------- held-out enrichment
def _group_matched_draw(rng, groups: np.ndarray, target_counts: dict, eligible: np.ndarray) -> np.ndarray:
    """One random cell set with the same per-group size as the tail, among eligible cells."""
    pick = []
    for g, cnt in target_counts.items():
        pool = np.flatnonzero((groups == g) & eligible)
        if cnt > len(pool):
            cnt = len(pool)
        if cnt:
            pick.append(rng.choice(pool, cnt, replace=False))
    return np.concatenate(pick) if pick else np.zeros(0, dtype=int)


def heldout_enrichment(labels: np.ndarray, label: str, tail: np.ndarray, groups: np.ndarray | None,
                       *, draws: int = CRITERIA["null_draws"], q: float = CRITERIA["null_quantile"],
                       seed: int = 0) -> dict:
    """h = P(label | tail cell) / P(label | judged cell) against a size- and group-matched null.

    Judged cells are those the arm did not mark Unavailable.  The withheld marker's tail is
    computed on the ORIGINAL values (never on the ablated input); the labels are the ablated
    arm's calls.  Returns the statistic, the null p95, and whether it was judgeable.
    """
    labels = np.asarray(labels).astype(str)
    judged = labels != UNAVAILABLE
    if groups is None:
        groups = np.zeros(len(labels), dtype=int)
    groups = np.asarray(groups)
    is_L = labels == label
    n_tail = int((tail & judged).sum())
    n_L = int((is_L & judged).sum())
    reason = gates.require_evidence(n_tail, CRITERIA["min_tail_in_pool"], f"{label}: tail cells in pool") \
        or gates.require_evidence(n_L, CRITERIA["min_calls"], f"{label}: calls")
    base = is_L[judged].mean() if judged.any() else np.nan
    if reason or not np.isfinite(base) or base <= 0:
        return {"label": label, "h": np.nan, "null_p": np.nan, "n_tail": n_tail, "n_calls": n_L,
                "judgeable": False, "clear": False, "reason": reason or "no calls for label"}
    h = is_L[tail & judged].mean() / base
    rng = np.random.default_rng(seed)
    counts = pd.Series(groups[tail & judged]).value_counts().to_dict()
    null = np.empty(draws)
    for i in range(draws):
        pick = _group_matched_draw(rng, groups, counts, judged)
        null[i] = is_L[pick].mean() / base if len(pick) else np.nan
    p = float(np.nanquantile(null, q))
    return {"label": label, "h": float(h), "null_p": p, "null_mean": float(np.nanmean(null)),
            "n_tail": n_tail, "n_calls": n_L, "judgeable": True, "clear": bool(h > p), "reason": ""}


# --------------------------------------------------------------------------- spatial coherence
def spatial_coherence(labels: np.ndarray, label: str, xy: np.ndarray, groups: np.ndarray | None, *,
                      k: int = CRITERIA["morans_k"], draws: int = CRITERIA["null_draws"],
                      q: float = CRITERIA["null_quantile"], seed: int = 0, max_cells: int = 200_000,
                      ids: np.ndarray | None = None) -> dict:
    """Moran's I of the one-hot label against within-group label permutations."""
    labels = np.asarray(labels).astype(str)
    judged = labels != UNAVAILABLE
    idx = np.flatnonzero(judged)
    if len(idx) > max_cells:
        sub = hash_subsample(ids[idx] if ids is not None else idx.astype(str), max_cells, salt=f"moran|{label}")
        idx = idx[sub]
    lab = labels[idx]
    n_L = int((lab == label).sum())
    reason = gates.require_evidence(n_L, CRITERIA["min_calls"], f"{label}: calls")
    if reason:
        return {"label": label, "I": np.nan, "null_p": np.nan, "n_calls": n_L, "judgeable": False, "clear": False, "reason": reason}
    pts = np.asarray(xy, dtype=float)[idx]
    nbr = usability.neighbour_graph(pts, k=k)
    one = (lab == label).astype(float)
    I_obs = float(usability.morans_i(one, nbr))
    g = (np.asarray(groups)[idx] if groups is not None else np.zeros(len(idx), dtype=int))
    rng = np.random.default_rng(seed)
    null = np.empty(draws)
    for d in range(draws):
        perm = one.copy()
        for gv in np.unique(g):
            sel = np.flatnonzero(g == gv)
            perm[sel] = one[rng.permutation(sel)]
        null[d] = usability.morans_i(perm, nbr)
    p = float(np.nanquantile(null, q))
    return {"label": label, "I": I_obs, "null_p": p, "null_mean": float(np.nanmean(null)), "n_calls": n_L,
            "n_cells_used": int(len(idx)), "judgeable": True, "clear": bool(I_obs > p), "reason": ""}


# --------------------------------------------------------------------------- descriptive statistics
def contamination(labels: np.ndarray, label: str, other_tails: dict, *, draws: int = CRITERIA["null_draws"],
                  q: float = CRITERIA["null_quantile"], seed: int = 0) -> dict:
    """Share of the label's calls that sit in ANOTHER label's marker tail, against a size-matched null."""
    labels = np.asarray(labels).astype(str)
    judged = labels != UNAVAILABLE
    is_L = (labels == label) & judged
    n_L = int(is_L.sum())
    if n_L < CRITERIA["min_calls"] or not other_tails:
        return {"label": label, "share": np.nan, "null_p": np.nan, "n_calls": n_L, "judgeable": False}
    any_other = np.zeros(len(labels), dtype=bool)
    for m, t in other_tails.items():
        any_other |= np.asarray(t, dtype=bool)
    share = float(any_other[is_L].mean())
    rng = np.random.default_rng(seed)
    pool = np.flatnonzero(judged)
    null = np.array([any_other[rng.choice(pool, n_L, replace=False)].mean() for _ in range(draws)])
    return {"label": label, "share": share, "null_p": float(np.quantile(null, q)),
            "null_mean": float(null.mean()), "n_calls": n_L, "judgeable": True}


def islet_core_audit(labels: np.ndarray, cell_region: np.ndarray, label: str = "Endocrine") -> dict:
    """AUDIT ONLY.  Enrichment of the label inside QuPath islet cores vs the donor base rate."""
    labels = np.asarray(labels).astype(str)
    judged = labels != UNAVAILABLE
    core = np.asarray(cell_region).astype(str) == "core"
    base = (labels[judged] == label).mean() if judged.any() else np.nan
    inside = (labels[judged & core] == label).mean() if (judged & core).any() else np.nan
    return {"label": label, "base_rate": float(base), "rate_in_core": float(inside),
            "enrichment": float(inside / base) if base and np.isfinite(base) and base > 0 else np.nan,
            "n_core": int((judged & core).sum()), "audit_only": True}


def agreement(a: np.ndarray, b: np.ndarray) -> dict:
    from sklearn.metrics import adjusted_mutual_info_score, adjusted_rand_score
    a = np.asarray(a).astype(str); b = np.asarray(b).astype(str)
    both = (a != UNAVAILABLE) & (b != UNAVAILABLE)
    if both.sum() < CRITERIA["min_calls"]:
        return {"ari": np.nan, "ami": np.nan, "n": int(both.sum())}
    return {"ari": float(adjusted_rand_score(a[both], b[both])),
            "ami": float(adjusted_mutual_info_score(a[both], b[both])), "n": int(both.sum())}


# --------------------------------------------------------------------------- cluster annotation
def standardise_within_donor(X: np.ndarray) -> np.ndarray:
    mu = np.nanmean(X, axis=0)
    sd = np.nanstd(X, axis=0)
    sd = np.where(sd > 0, sd, 1.0)
    Z = (X - mu) / sd
    return np.where(np.isfinite(Z), Z, 0.0)


def null_gated_cluster_annotation(Z: np.ndarray, clusters: np.ndarray, declared: dict, *,
                                  draws: int = CRITERIA["null_draws"], q: float = CRITERIA["null_quantile"],
                                  margin_frac: float = 0.25, seed: int = 0,
                                  min_cells: int = CRITERIA["min_calls"]) -> pd.DataFrame:
    """Award a broad label to a cluster only where its declared-marker contrast beats a size-matched
    null AND ranks first by a margin.  ``declared`` maps label -> column indices into ``Z``.

    contrast_L(c) = mean over M_L of ( mean Z[c, m] - mean Z[~c, m] ).  Null: ``draws`` random
    cell sets of the cluster's size.  Multiple labels above their null without a clear winner ->
    Ambiguous(multi_lineage); none -> Unavailable(no_lineage_beats_null); too small ->
    Unavailable(cluster_below_min_cells).  The false-award rate on synthetic noise is measured by
    the tests and recorded beside every table that uses this.
    """
    rng = np.random.default_rng(seed)
    clusters = np.asarray(clusters)
    n = Z.shape[0]
    col_mean = Z.mean(axis=0)
    rows = []
    labels = list(declared)
    for c in np.unique(clusters):
        inc = clusters == c
        nc = int(inc.sum())
        if nc < min_cells:
            rows.append({"cluster": c, "n": nc, "label": UNAVAILABLE, "reason": "cluster_below_min_cells"})
            continue
        in_mean = Z[inc].mean(axis=0)
        out_mean = (col_mean * n - in_mean * nc) / max(n - nc, 1)
        contrast = {L: float(np.mean(in_mean[cols] - out_mean[cols])) for L, cols in declared.items()}
        null = {L: np.empty(draws) for L in labels}
        for d in range(draws):
            pick = rng.choice(n, nc, replace=False)
            pm = Z[pick].mean(axis=0)
            om = (col_mean * n - pm * nc) / max(n - nc, 1)
            for L, cols in declared.items():
                null[L][d] = np.mean(pm[cols] - om[cols])
        p95 = {L: float(np.quantile(null[L], q)) for L in labels}
        above = [L for L in labels if contrast[L] > max(p95[L], 0.0)]
        ranked = sorted(labels, key=lambda L: -contrast[L])
        top, second = ranked[0], (ranked[1] if len(ranked) > 1 else None)
        rec = {"cluster": c, "n": nc, **{f"contrast_{L}": contrast[L] for L in labels},
               **{f"p95_{L}": p95[L] for L in labels}, "n_above_null": len(above)}
        if not above:
            rec.update(label=UNAVAILABLE, reason="no_lineage_beats_null")
        elif top in above and (second is None or contrast[top] - contrast[second] >= margin_frac * abs(contrast[top])):
            rec.update(label=top, reason="")
        else:
            rec.update(label=AMBIGUOUS, reason="multi_lineage")
        rows.append(rec)
    return pd.DataFrame(rows)


# --------------------------------------------------------------------------- the judgement
def judge_label(label: str, primary_labels: np.ndarray, ablations: dict, replicates: dict, *,
                values: pd.DataFrame, xy: np.ndarray, groups: np.ndarray | None, ids: np.ndarray,
                declared_markers: list, held_out_gated: bool, other_tails: dict,
                cell_region: np.ndarray | None, seed: int = 0) -> dict:
    """Blocking + descriptive statistics for one label on one donor.

    ``ablations`` maps withheld marker -> the ablated arm's labels (aligned to ``ids``);
    ``replicates`` maps replicate name -> labels.  ``values`` holds the ORIGINAL aligned marker
    values (columns = markers) for the same cells.
    """
    out = {"label": label, "declared_markers": list(declared_markers), "held_out_gated": held_out_gated}
    # 1. held-out enrichment
    ho = []
    if held_out_gated:
        for m in declared_markers:
            if m not in ablations:
                ho.append({"marker": m, "judgeable": False, "clear": False, "reason": "ablation_not_run"})
                continue
            if m not in values.columns:
                ho.append({"marker": m, "judgeable": False, "clear": False, "reason": "marker_not_acquired"})
                continue
            tail = tail_mask(values[m].to_numpy())
            r = heldout_enrichment(np.asarray(ablations[m]), label, tail, groups, seed=seed)
            ho.append({"marker": m, **r})
        judged = [r for r in ho if r["judgeable"]]
        cleared = [r for r in judged if r["clear"]]
        need = int(np.ceil(CRITERIA["heldout_frac"] * len(judged))) if judged else None
        out["heldout"] = {"per_marker": ho, "n_judged": len(judged), "n_clear": len(cleared),
                          "n_required": need,
                          "verified": bool(judged and len(judged) >= 2 and len(cleared) >= need),
                          "reason": "" if (judged and len(judged) >= 2) else "fewer than 2 judgeable held-out markers"}
    else:
        out["heldout"] = {"per_marker": [], "verified": None,
                          "reason": "not held-out gated (mutually exclusive subtypes or a single declared marker)"}
    # 2. spatial coherence on the primary calls
    out["spatial"] = spatial_coherence(primary_labels, label, xy, groups, ids=ids, seed=seed)
    # 3. stability: the primary verdict components must hold on every replicate
    stab = []
    for name, lab in replicates.items():
        s = spatial_coherence(np.asarray(lab), label, xy, groups, ids=ids, seed=seed)
        stab.append({"replicate": name, "spatial_clear": s["clear"], "judgeable": s["judgeable"], "n_calls": s["n_calls"]})
    out["stability"] = {"per_replicate": stab,
                        "invariant": all(r["spatial_clear"] == out["spatial"]["clear"] for r in stab) if stab else None}
    # descriptive
    out["contamination"] = contamination(primary_labels, label, other_tails, seed=seed)
    if cell_region is not None and label == "Endocrine":
        out["islet_core_audit"] = islet_core_audit(primary_labels, cell_region, label)
    # verdict for the label
    blocking = [out["spatial"]["clear"]]
    if held_out_gated:
        blocking.append(bool(out["heldout"]["verified"]))
    if stab:
        blocking.append(bool(out["stability"]["invariant"]))
    out["verified"] = bool(blocking) and all(blocking)
    return out
