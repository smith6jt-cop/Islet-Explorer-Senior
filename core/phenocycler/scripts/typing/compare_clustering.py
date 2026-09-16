#!/usr/bin/env python3
"""Choose the lane's clustering method on evidence the methods cannot tune.

The maintainer asked for the better of Leiden, Latent Dirichlet Allocation and KODAMA.  Each is a
separate module (``cluster_leiden.py``, ``cluster_lda.py``, ``cluster_kodama.py``) exposing
``fit_predict(X, Z, *, seed=0, **params)``; this file runs them on the same cells and scores them.
Backends are imported defensively -- a missing or broken one is reported and skipped, because a
comparison that dies on the module that is not written yet is useless.

**Nothing here is a criterion a clustering method could optimise against.**  That is the whole
point, and it is why the obvious scores are absent: there is no accuracy against a label set, no
silhouette or Calinski-Harabasz (a method that minimises within-cluster distance wins those by
construction), and no comparison to an expected cell-type frequency -- the lane's rules forbid the
last one outright.  What is scored instead:

``islet agreement``
    QuPath's islet polygons are hand-drawn geometry.  No marker, no model and no clustering saw
    them.  Islet-core cells are ~5.5% of the test window, so a partition that concentrates them
    into few clusters that are mostly islet has found real structure.  Reported three ways --
    which clusters they went to, how many clusters hold 90% of them against how many the same
    cells scattered at random would need, and the purity of that set -- because each alone is
    gameable: purity by over-splitting, concentration by one giant cluster.

``rare-population recovery``
    41 cells in the test window carry high SST.  Delta cells are the rarest thing the panel can
    name, and the study is about a minority population (immune cells near islets), so burying a
    0.35% population in a big cluster is the failure mode that matters most here.  The reference
    set is the widest gap between SST values above the median -- a property of the distribution,
    not a cutoff chosen by anyone, and independent of how many cells are in the table.  On the
    test window it falls between 0.499 and 0.436, which is where Nimbus's own decision boundary
    sits, so the rule and the model pick out the same cells.  It scores a partition; it is never a
    feature and no method sees it.

``held-out marker agreement``
    The sharpest check.  call.py's naming statistic reads the top ``TOP_N`` = 3 of a lineage's
    markers, so a lineage that declares more than three has markers to spare.  Here each such
    lineage is cut at random into a naming subset of 3 and a held-out remainder, the cluster goes
    to whichever lineage scores highest on its subset, and the remainder is then asked whether it
    is also elevated in that cluster.  Markers the naming never saw agreeing with the name is
    evidence the cluster is a cell type rather than an artefact of the features.  Scored against
    every marker set of the same size from outside the lineage where that set can be enumerated
    (30 subsets for Vascular's one held-out marker, 3,276 for Exocrine's three; only Myeloid's
    five overflow and are sampled), and averaged over many random splits so no single lucky cut
    decides it.  Read the ordering and not the level: the level still moves with the split count,
    while a structured partition beats a random one of the same k at every count tried.  Only the
    mean of the held-out markers is used, never their top few: taking a maximum over the remainder
    would smuggle the naming statistic's own selection back in.  This is call.py's statistic and
    not its Benjamini-Hochberg policy -- the policy answers "is this cluster safe to name", which
    is the ambiguous row's question, and on this window it named 4 of 48 clusters and left nothing
    to score.

``spatial coherence``
    Cells of one type sit together.  Share of each cell's 12 nearest neighbours carrying the same
    label, against labels shuffled *within* 500 um blocks -- local, so the null destroys cell-scale
    structure while preserving where the tissue's compartments are.  Reported as the ratio to the
    null rather than raw, because the null preserves each method's own k and cluster sizes and the
    raw share falls with k for any method.

``ambiguous fraction``
    call.py's own naming, run unchanged on each partition: how many cells sit in a cluster no
    lineage can claim.  Lower is better **only if the other criteria hold** -- a method that names
    everything confidently and wrongly scores best on this one alone, which is why it is one row of
    seven and never the decision.

``stability``
    Adjusted Rand index between three seeds.  A partition that does not survive its own random
    seed is not a finding, whatever it scores elsewhere.

No criterion decides on its own: the recommendation is the best mean rank across all seven, every
per-criterion value and null is printed, and any criterion where the recommended method is last is
called out so the trade-off is visible rather than hidden in an aggregate.

    ../../.venv/bin/python scripts/typing/compare_clustering.py --donors 6442
"""

from __future__ import annotations

import argparse
import glob
import importlib
import inspect
import json
import sys
import time
import traceback
from itertools import combinations
from math import comb
from pathlib import Path

import numpy as np
import pandas as pd
import pyarrow.parquet as pq
from sklearn.metrics import adjusted_rand_score
from sklearn.neighbors import NearestNeighbors

sys.path.insert(0, str(Path(__file__).parent))
import lineages as L  # noqa: E402
from call import TOP_N, _statistic, name_clusters  # noqa: E402

REPO = Path(__file__).resolve().parents[4]
DEV = REPO / "data/development/typing"
#: The validated Nimbus window: one 2048 px tile of 6442 centred on islets, 13,648 cells of which
#: 11,659 are analysis-eligible and 638 sit inside a QuPath islet polygon.  Kept as a source in its
#: own right because it is the only positivity that has been checked against the polygons
#: (GCG AUC 0.996, INS 0.991), and because a whole slide takes ~2.5 h of GPU to produce.
WINDOW = Path("/data/nimbus_typing_test")

METHODS = ("leiden", "lda", "kodama")
#: Three seeds is the minimum that gives more than one pairwise ARI and so distinguishes "unstable"
#: from "one unlucky run".
SEEDS = (0, 1, 2)

K_NEIGHBOURS = 12          # validate.py's, so the two files' coherence numbers are comparable
BLOCK_UM = 500.0           # validate.py's
COVER = 0.90               # "how few clusters hold 90% of the islet cells"
#: Marker whose rare high tail is used for the rare-population check.  SST because delta cells are
#: the least abundant thing the panel can name -- 41 of 11,659 cells in the test window -- and a
#: method that dissolves them is exactly what this criterion is for.
RARE_MARKER = "SST"

#: call.py's naming is a family of ``k x lineages`` tests corrected with Benjamini-Hochberg, and
#: it refuses to run on fewer than ``m_tests / alpha`` null draws.  The harness sizes the draws to
#: the partition it is scoring, because a method that returns many clusters would otherwise be
#: skipped for a reason that has nothing to do with the quality of its clusters.  Past this cap it
#: is reported as unaffordable instead: at 1500 clusters call.py needs 270,000 draws per naming,
#: and this harness names each partition once plus once per held-out split.
NAMING_DRAW_CAP = 120_000

#: The held-out markers' null is the set of all same-size marker subsets from outside the lineage,
#: and for most lineages that set is small enough to enumerate: 30 subsets for Vascular's one
#: held-out marker, 406 for T_NK's two, 3,276 for Exocrine's three.  Enumerated, the null is exact
#: and has no draw count to argue about.  Only Myeloid's five held-out markers overflow this cap
#: (65,780 subsets), and only that one is sampled.
EXACT_NULL_MAX = 20_000
#: Draws for the one lineage whose subset space is too large to enumerate.  With the null drawing
#: from its own generator, the count barely matters: on a 40-cluster partition of the test window
#: the cell-weighted excess reads +0.0672 (250 draws), +0.0670 (500), +0.0684 (1000), +0.0683
#: (2000), +0.0680 (4000), +0.0677 (8000) -- a spread of 0.0014 over a 32-fold range, and 2,000
#: costs 0.03 s.  (Before the generators were separated this looked like a drift from +0.089 to
#: +0.028; that was the null sampling consuming the stream the marker splits were drawn from.)
HELDOUT_NULL_DRAWS = 2_000
#: Resamples of the marker splits, used only to put a spread on the held-out excess.  200 costs
#: 0.4 s at 2,000 scored cluster-splits and settles the 5th/95th percentiles to ~0.002.
BOOT_DRAWS = 200

#: (key, description, +1 if larger is better).  Every one of these is computed from something the
#: methods cannot see or influence.
CRITERIA = (
    ("islet_purity",        "purity of the clusters holding 90% of islet-core cells", +1),
    ("islet_concentration", "islet cells packed into fewer clusters than chance, x",   +1),
    ("rare_f1",             f"best cluster F1 for the {RARE_MARKER}-high cells",       +1),
    ("heldout_excess",      "held-out marker elevation above its null (SD units)",     +1),
    ("spatial_enrichment",  "same-label neighbour share / block-shuffled null",        +1),
    ("ambiguous_frac",      "cells in clusters call.py cannot name",                   -1),
    ("stability_ari",       f"mean ARI across seeds {list(SEEDS)}",                    +1),
)


# --------------------------------------------------------------------------- backends

def load_backend(name: str):
    """``cluster_<name>``, or a recorded reason it is unusable.  Never raises."""
    try:
        mod = importlib.import_module(f"cluster_{name}")
    except Exception as exc:                       # not written yet, or its dependency is missing
        return None, f"import failed: {type(exc).__name__}: {exc}"
    fn = getattr(mod, "fit_predict", None)
    if not callable(fn):
        return None, "module has no callable fit_predict"
    return mod, None


def check_labels(raw, n: int) -> np.ndarray:
    """The interface, enforced.  A backend that breaks it is skipped, not silently repaired."""
    a = np.asarray(raw)
    if a.ndim != 1 or a.shape[0] != n:
        raise ValueError(f"expected {n} labels, got shape {a.shape}")
    if not np.issubdtype(a.dtype, np.integer):
        raise ValueError(f"labels must be integer, got {a.dtype}")
    u = np.unique(a)
    if u.size == 0 or u[0] != 0 or not np.array_equal(u, np.arange(u.size)):
        raise ValueError(f"labels must be 0..k-1 with no gaps; got {u[:8]}... ({u.size} distinct)")
    return a.astype(np.int32)


# --------------------------------------------------------------------------- the cells

def _dedupe(pos: pd.DataFrame) -> tuple[pd.DataFrame, int]:
    """One row per cell.  A cell whose polygon spans two Nimbus tiles is rasterized into both, so
    the whole-slide table can carry it twice; the tile holding more of it is the one to keep."""
    if not pos["object_id"].duplicated().any():
        return pos, 0
    kept = pos.sort_values("area", ascending=False).drop_duplicates("object_id", keep="first")
    return kept.sort_index(), int(len(pos) - len(kept))


def _window_table(donor: str) -> pd.DataFrame:
    masks = sorted(glob.glob(str(WINDOW / f"mask_{donor}_*.npz")))
    if not masks:
        raise FileNotFoundError(f"no validated window for donor {donor} under {WINDOW}")
    tab = pd.read_parquet(WINDOW / "nimbus_out/cell_table.parquet")
    oids = np.load(masks[0], allow_pickle=True)["oids"]
    lab = tab["label"].to_numpy(np.int64)
    if lab.min() < 1 or lab.max() > len(oids):
        raise ValueError(f"labels {lab.min()}..{lab.max()} do not index {len(oids)} oids")
    tab["object_id"] = np.asarray(oids, dtype=object)[lab - 1].astype(str)
    return tab


def _positivity_path(donor: str) -> Path:
    return DEV / f"positivity/donor_id={donor}/data.parquet"


def choose_source(donor: str, source: str) -> str:
    """``auto`` takes whichever table covers more cells, and says so.

    The whole-slide run overwrites the same path a smoke test wrote, so "does the file exist" is
    not the question -- "does it have more cells in it than the validated window" is.
    """
    if source != "auto":
        return source
    p = _positivity_path(donor)
    n_pos = pq.ParquetFile(p).metadata.num_rows if p.exists() else 0
    n_win = 0
    if glob.glob(str(WINDOW / f"mask_{donor}_*.npz")):
        n_win = pq.ParquetFile(WINDOW / "nimbus_out/cell_table.parquet").metadata.num_rows
    if n_pos == 0 and n_win == 0:
        raise FileNotFoundError(f"donor {donor}: neither {p} nor a window under {WINDOW}")
    pick = "positivity" if n_pos >= n_win else "window"
    print(f"[{donor}] source=auto -> {pick}  (positivity {n_pos:,} rows, window {n_win:,} rows)")
    return pick


def prepare(donor: str, source: str, max_cells: int, seed: int) -> dict:
    """Join positivity to the substrate and build the two matrices the interface promises."""
    src = choose_source(donor, source)
    pos = _window_table(donor) if src == "window" else pd.read_parquet(_positivity_path(donor))
    pos["object_id"] = pos["object_id"].astype(str)
    pos, n_dup = _dedupe(pos)

    sub = pd.read_parquet(DEV / f"substrate/donor_id={donor}/data.parquet",
                          columns=["object_id", "eligible", "cell_region", "islet_num",
                                   "X_centroid", "Y_centroid"])
    sub["object_id"] = sub["object_id"].astype(str)
    markers = [c for c in pos.columns if c in set(L.identity_markers())]
    df = sub.merge(pos[["object_id", "area", *markers]], on="object_id", how="inner",
                   validate="one_to_one")
    # a positivity row whose object_id is in no substrate row means the join key is wrong, which
    # is the failure this lane is most careful about -- say it rather than quietly losing cells
    if len(df) < len(pos):
        print(f"[{donor}] {len(pos) - len(df):,} of {len(pos):,} positivity rows matched no "
              f"substrate object_id")

    X = df[markers].to_numpy(np.float64)
    eligible = df["eligible"].to_numpy(bool)
    if not eligible.any():
        raise ValueError(f"donor {donor}: no eligible cell has a positivity row")

    # A channel that is absent on this donor, or that Nimbus returned flat, carries no information
    # and would make Z all-NaN for every cell.  Drop the channel, never the cells.
    Xe = X[eligible]
    keep = np.isfinite(Xe).any(axis=0)
    keep[keep] = np.nanstd(Xe[:, keep], axis=0) > 0
    if not keep.all():
        print(f"[{donor}] dropped {(~keep).sum()} channel(s) with no usable variation: "
              f"{[m for m, k in zip(markers, keep) if not k]}")
        markers = [m for m, k in zip(markers, keep) if k]
        X = X[:, keep]

    # A cell with no centroid cannot be scored spatially and has no position to score; it is
    # dropped with the others that carry an incomplete measurement, not filled and not silently
    # passed to a nearest-neighbour search that would raise inside sklearn.
    placed = np.isfinite(df[["X_centroid", "Y_centroid"]].to_numpy(float)).all(axis=1)
    if (eligible & ~placed).any():
        print(f"[{donor}] {(eligible & ~placed).sum():,} eligible cells have no centroid")
    have = eligible & placed & np.isfinite(X).all(axis=1)
    df, X = df.loc[have].reset_index(drop=True), X[have]

    if max_cells and len(df) > max_cells:
        # Subsampling changes what there is to find, so it is opt-in and announced, never silent.
        idx = np.random.default_rng(seed).choice(len(df), max_cells, replace=False)
        idx.sort()
        df, X = df.iloc[idx].reset_index(drop=True), X[idx]

    # Same standardisation as call.py: channels Nimbus scores broadly positive (TP63 in 98.6% of
    # cells, FOXP3 in 90.4%) contribute their variation and not their level.
    mu, sd = X.mean(axis=0), X.std(axis=0)
    Z = (X - mu) / sd

    core = (df["cell_region"] == "core").to_numpy()
    rare, rare_meta = rare_reference(X[:, markers.index(RARE_MARKER)]) \
        if RARE_MARKER in markers else (np.zeros(len(df), bool), {"n": 0})
    meta = {"source": src, "n_positivity_rows": int(len(pos)), "n_duplicate_rows": n_dup,
            "n_joined": int(len(eligible)),
            "n_positivity_unmatched": int(len(pos) - len(eligible)),
            "n_eligible": int(eligible.sum()),
            "n_typed": int(len(df)), "n_dropped_incomplete": int(eligible.sum() - len(df)),
            "n_markers": len(markers), "markers": markers,
            "n_core": int(core.sum()), "core_share": round(float(core.mean()), 5),
            "rare_marker": RARE_MARKER, "rare": rare_meta}
    print(f"[{donor}] {len(df):,} cells x {len(markers)} markers | islet-core {core.sum():,} "
          f"({core.mean()*100:.2f}%) | {RARE_MARKER}-high {rare_meta['n']}")
    if 0 < rare_meta["n"] < 10:
        print(f"[{donor}] WARNING: the {RARE_MARKER} gap left only {rare_meta['n']} cells; the "
              f"rare-population criterion is noise at that count -- read its null beside it")
    return {"donor": donor, "df": df, "X": X, "Z": Z, "markers": markers,
            "core": core, "rare": rare, "meta": meta}


def rare_reference(v: np.ndarray) -> tuple[np.ndarray, dict]:
    """The high mode of one marker, cut at the widest gap between values above its median.

    Not a threshold on the confidences: it is read off this donor's own distribution after the
    fact, and it never enters a feature, a call or a parameter -- it only scores a partition that
    was built without it.

    The search runs over every value above the median rather than a fixed number of top cells, so
    it does not depend on how many cells are in the table.  On the test window it cuts at
    0.499 -> 0.436, a gap of 0.063 where the gaps at the very top of the same tail are 0.006 --
    the trough between two modes, not sparsity at the extreme.  That is 41 cells, the same 41 a
    top-200-rank search returns and the same place Nimbus's own decision boundary falls; at 30%
    and 50% subsamples the rule holds the share (0.52%, 0.36%, 0.35%).  A denser table deepens
    the contrast rather than blurring it, because the tail fills in and the trough does not.
    """
    s = np.sort(v)[::-1]
    s = s[s >= np.median(v)]
    if len(s) < 2:
        return np.zeros(len(v), bool), {"n": 0, "error": "no variation above the median"}
    gaps = s[:-1] - s[1:]
    i = int(np.argmax(gaps))
    cut = float(s[i])
    sel = v >= cut
    return sel, {"n": int(sel.sum()), "rank_of_gap": i + 1, "cut": cut,
                 "next_value": float(s[i + 1]), "gap": float(gaps[i]),
                 "n_searched": int(len(s))}


# --------------------------------------------------------------------------- the criteria

def _cover(in_set: np.ndarray, sizes: np.ndarray, frac: float) -> tuple[int, float]:
    """Fewest clusters holding ``frac`` of a marked population, and their pooled purity."""
    total = in_set.sum()
    if total == 0:
        return 0, float("nan")
    order = np.argsort(-in_set)
    n = int(np.searchsorted(np.cumsum(in_set[order]), frac * total) + 1)
    take = order[:n]
    return n, float(in_set[take].sum() / sizes[take].sum())


def marked_recovery(labels: np.ndarray, marked: np.ndarray, k: int, draws: int,
                    seed: int) -> dict:
    """How a partition holds a population it never saw, against the same population scattered.

    The null draws the same number of cells uniformly, so it carries this method's own k and
    cluster sizes; it is what "no relationship at all" looks like for *this* partition, which is
    the only fair reference when the methods return different numbers of clusters.
    """
    sizes = np.bincount(labels, minlength=k).astype(np.float64)
    in_set = np.bincount(labels[marked], minlength=k).astype(np.float64)
    total = float(marked.sum())
    n90, purity = _cover(in_set, sizes, COVER)

    prec = np.divide(in_set, sizes, out=np.zeros(k), where=sizes > 0)
    rec = in_set / total if total else np.zeros(k)
    f1 = np.divide(2 * prec * rec, prec + rec, out=np.zeros(k), where=(prec + rec) > 0)
    b = int(np.argmax(f1))

    rng = np.random.default_rng(seed)
    n_null, p_null, f_null = np.empty(draws), np.empty(draws), np.empty(draws)
    for d in range(draws):
        pick = rng.choice(len(labels), int(total), replace=False)
        c = np.bincount(labels[pick], minlength=k).astype(np.float64)
        n_null[d], p_null[d] = _cover(c, sizes, COVER)
        pr = np.divide(c, sizes, out=np.zeros(k), where=sizes > 0)
        re = c / total if total else np.zeros(k)
        f_null[d] = np.max(np.divide(2 * pr * re, pr + re, out=np.zeros(k), where=(pr + re) > 0))
    # how much tighter than chance: the null carries this method's own k and cluster sizes, so
    # the ratio is comparable across methods that return very different numbers of clusters --
    # a raw count of clusters would hand the criterion to whichever method returns fewest.
    conc = float(np.percentile(n_null, 5) / n90) if n90 else float("nan")
    return {"n_marked": int(total), "clusters_90": n90, "concentration": round(conc, 3),
            "purity_90": round(purity, 4), "best_cluster": b,
            "best_precision": round(float(prec[b]), 4),
            "best_recall": round(float(rec[b]), 4), "best_f1": round(float(f1[b]), 4),
            "best_cluster_size": int(sizes[b]),
            "null_clusters_90_p05": float(np.percentile(n_null, 5)),
            "null_purity_90_p95": round(float(np.percentile(p_null, 95)), 4),
            "null_best_f1_p95": round(float(np.percentile(f_null, 95)), 4),
            "null_draws": draws}


def composition(labels: np.ndarray, marked: np.ndarray, k: int, top: int = 6) -> list[dict]:
    """Where the marked cells actually went: the largest clusters holding them, with purity."""
    sizes = np.bincount(labels, minlength=k)
    in_set = np.bincount(labels[marked], minlength=k)
    order = np.argsort(-in_set)[:top]
    return [{"cluster": int(c), "n_marked": int(in_set[c]), "cluster_size": int(sizes[c]),
             "purity": round(float(in_set[c] / sizes[c]), 4) if sizes[c] else None,
             "share_of_marked": round(float(in_set[c] / max(marked.sum(), 1)), 4)}
            for c in order if in_set[c] > 0]


def cluster_means(Z: np.ndarray, labels: np.ndarray, k: int) -> np.ndarray:
    """Per-cluster mean of every marker.

    One weighted bincount per marker rather than k boolean masks: at 1.4M cells the mask form is k
    full passes over the matrix, and k is 20-40 here.  Z carries no NaN -- cells with an incomplete
    measurement were excluded in :func:`prepare`, never filled -- so a plain mean is the nanmean.
    """
    n = np.bincount(labels, minlength=k).astype(np.float64)
    M = np.empty((k, Z.shape[1]))
    for j in range(Z.shape[1]):
        s = np.bincount(labels, weights=Z[:, j], minlength=k)
        M[:, j] = np.divide(s, n, out=np.full(k, np.nan), where=n > 0)
    return M


class NamingTooExpensive(RuntimeError):
    """This partition needs more null draws than call.py's naming can afford."""


def name_from_means(M: np.ndarray, markers: list[str], floor: int, margin: float, alpha: float,
                    seed: int):
    """call.py's naming, unchanged, on one row per cluster.

    ``name_clusters`` takes the mean of each cluster as its first act, so handing it the means it
    would compute is the same computation on a ``k x markers`` matrix rather than an
    ``n x markers`` one -- 0.6 s at 1.4M cells, and identical output.

    Bound by keyword and checked against the live signature on purpose: call.py is under active
    change (it grew ``alpha`` and Benjamini-Hochberg during this comparison), and a positional call
    turns that into a cryptic TypeError inside a backend's result instead of a legible failure
    here.
    """
    lineages = [lin for lin, ms in L.LINEAGE_MARKERS.items() if any(m in markers for m in ms)]
    draws = max(floor, int(np.ceil(len(M) * len(lineages) / alpha)))
    if draws > NAMING_DRAW_CAP:
        raise NamingTooExpensive(
            f"{len(M)} clusters x {len(lineages)} lineages at alpha={alpha} needs {draws:,} null "
            f"draws for Benjamini-Hochberg, over the {NAMING_DRAW_CAP:,} this harness will spend")
    have = {"draws": draws, "margin": margin, "alpha": alpha, "seed": seed}
    params = inspect.signature(name_clusters).parameters
    missing = [n for n, prm in params.items()
               if n not in ("Z", "labels", "markers", "k") and n not in have
               and prm.default is inspect.Parameter.empty]
    if missing:
        raise RuntimeError(f"call.name_clusters now requires {missing}, which this harness has no "
                           f"value for; update compare_clustering.py deliberately")
    named, why, *_ = name_clusters(M, np.arange(len(M)), markers, len(M),
                                   **{n: v for n, v in have.items() if n in params})
    return named, why, draws


def heldout_agreement(M: np.ndarray, markers: list[str], sizes: np.ndarray, splits: int,
                      seed: int) -> dict:
    """Name each cluster from 3 of a lineage's markers; ask whether the rest agree.

    Only lineages declaring more than TOP_N markers can be cut this way -- on this panel that is
    Exocrine (6), Vascular (4), T_NK (5) and Myeloid (8), leaving 11 of the 34 channels held out.
    A lineage with three or fewer markers still competes for the name; it just contributes no
    held-out evidence, which is honest rather than fixable.

    This uses call.py's naming *statistic* and not its naming *policy*: the cluster goes to
    whichever lineage scores highest on its naming subset, with no null and no FDR.  The question
    here is whether a cluster is a cell type, and call.py's Benjamini-Hochberg bar answers a
    different one -- measured on this window it named 4 of 48 leiden clusters, which left 7
    cluster-splits to score and made the criterion noise.  The ambiguous-fraction criterion is
    where call.py's policy is judged.
    """
    idx = {m: i for i, m in enumerate(markers)}
    present = {lin: [m for m in ms if m in idx] for lin, ms in L.LINEAGE_MARKERS.items()}
    present = {lin: ms for lin, ms in present.items() if ms}
    if not any(len(ms) > TOP_N for ms in present.values()):
        return {"error": f"no lineage declares more than TOP_N={TOP_N} markers on this panel"}
    lins = list(present)

    # two generators, not one: if the null sampling drew from the same stream, changing how many
    # null draws it takes would silently change which markers the splits hold out, and the two
    # would be impossible to tell apart in the result
    rng = np.random.default_rng(seed)
    rng_null = np.random.default_rng(seed + 1_000_003)
    rows, nulls, per_split = [], {}, []
    for s in range(splits):
        naming, held = {}, {}
        for lin, ms in present.items():
            pick = (list(rng.choice(np.asarray(ms, dtype=object), TOP_N, replace=False))
                    if len(ms) > TOP_N else list(ms))
            naming[lin], held[lin] = pick, [m for m in ms if m not in pick]
        S = np.vstack([_statistic(M, [idx[m] for m in naming[lin]]) for lin in lins]).T
        winner = np.asarray(lins, dtype=object)[np.nanargmax(np.nan_to_num(S, nan=-np.inf), 1)]

        cleared, scored = 0, 0
        for c, lin in enumerate(winner):
            cols = [idx[m] for m in held[lin]]
            if not cols:
                continue
            obs = float(np.nanmean(M[c, cols]))
            # the null is random markers of the same count from *outside* this lineage, so it asks
            # "would any set this size look this elevated here", not "is this cluster bright".
            # It depends only on the cluster and the lineage -- every split holds out the same
            # number of markers for a given lineage -- so it is computed once per pair.
            if (c, lin) not in nulls:
                pool = np.asarray([i for m, i in idx.items() if m not in present[lin]])
                exact = comb(len(pool), len(cols)) <= EXACT_NULL_MAX
                take = (np.asarray(list(combinations(range(len(pool)), len(cols)))) if exact else
                        np.argsort(rng_null.random((HELDOUT_NULL_DRAWS, len(pool))),
                                   axis=1)[:, :len(cols)])
                nulls[(c, lin)] = (float(np.nanpercentile(np.nanmean(M[c, pool[take]], axis=1),
                                                         95)), len(take), exact)
            p95, n_null, exact = nulls[(c, lin)]
            rows.append({"split": s, "cluster": int(c), "lineage": lin, "n_cells": int(sizes[c]),
                         "n_naming": len(naming[lin]), "n_heldout": len(cols),
                         "observed": round(obs, 4), "null_p95": round(p95, 4),
                         "null_subsets": n_null, "null_exact": exact,
                         "excess": round(obs - p95, 4), "cleared": bool(obs > p95)})
            cleared += obs > p95
            scored += 1
        per_split.append({"split": s, "n_scored": scored, "n_cleared": int(cleared)})

    if not rows:
        return {"error": "no cluster went to a lineage that had markers left over"}
    w = np.array([r["n_cells"] for r in rows], float)
    ex = np.array([r["excess"] for r in rows], float)
    ob = np.array([r["observed"] for r in rows], float)
    nu = np.array([r["null_p95"] for r in rows], float)
    cl = np.array([r["cleared"] for r in rows], float)

    # The excess is an average over random marker splits, and which three markers name a lineage
    # genuinely changes it.  Resampling the splits gives the spread that average carries, so a
    # reader can see whether two methods differ by more than the criterion's own noise -- without
    # it a 0.003 gap reads as a rank.
    by_split: dict[int, list] = {}
    for r in rows:
        by_split.setdefault(r["split"], []).append(r)
    ids, rs = list(by_split), np.random.default_rng(seed + 7)
    boot = np.empty(BOOT_DRAWS)
    for b in range(BOOT_DRAWS):
        sel = [r for i in rs.choice(ids, len(ids), replace=True) for r in by_split[i]]
        boot[b] = np.average([r["excess"] for r in sel],
                             weights=[r["n_cells"] for r in sel])
    lo, hi = (float(x) for x in np.percentile(boot, [5, 95]))

    return {"splits": splits, "n_scored": len(rows),
            "excess_cellweighted": round(float(np.average(ex, weights=w)), 4),
            "excess_p05": round(lo, 4), "excess_p95": round(hi, 4),
            "excess_mean": round(float(ex.mean()), 4),
            "observed_cellweighted": round(float(np.average(ob, weights=w)), 4),
            "null_p95_cellweighted": round(float(np.average(nu, weights=w)), 4),
            "share_cleared": round(float(cl.mean()), 4),
            "share_cleared_cellweighted": round(float(np.average(cl, weights=w)), 4),
            "null_subsets_min": int(min(r["null_subsets"] for r in rows)),
            "null_subsets_max": int(max(r["null_subsets"] for r in rows)),
            "exact_null_share": round(float(np.mean([r["null_exact"] for r in rows])), 3),
            "per_split": per_split,
            "worst": sorted(rows, key=lambda r: r["excess"])[:5],
            "best": sorted(rows, key=lambda r: -r["excess"])[:5]}


def spatial_coherence(xy: np.ndarray, labels: np.ndarray, draws: int, seed: int,
                      ind: np.ndarray | None = None) -> tuple[dict, np.ndarray]:
    """validate.py's check, cell-weighted over the whole partition.

    Reported as a ratio to the null because the null shuffles *within* the same partition, so it
    already carries this method's k and cluster sizes; the raw share falls with k for any method
    and would rank a coarse partition above a correct one.
    """
    if ind is None:
        nn = NearestNeighbors(n_neighbors=K_NEIGHBOURS + 1, n_jobs=-1).fit(xy)
        ind = nn.kneighbors(xy, return_distance=False)[:, 1:]
    obs = float((labels[ind] == labels[:, None]).mean())

    block = (xy // BLOCK_UM).astype(np.int64)
    key = block[:, 0] * 1_000_003 + block[:, 1]
    order = np.argsort(key, kind="stable")
    groups = [g for g in np.split(order, np.flatnonzero(np.diff(key[order])) + 1) if g.size > 1]
    rng = np.random.default_rng(seed)
    null = np.empty(draws)
    for d in range(draws):
        perm = labels.copy()
        for g in groups:
            perm[g] = labels[rng.permutation(g)]
        null[d] = (perm[ind] == perm[:, None]).mean()
    m = float(null.mean())
    return ({"observed": round(obs, 4), "null_mean": round(m, 4),
             "null_p95": round(float(np.percentile(null, 95)), 4),
             "enrichment": round(obs / m, 3) if m > 0 else None,
             "n_blocks": len(groups), "block_um": BLOCK_UM, "k_neighbours": K_NEIGHBOURS,
             "null_draws": draws}, ind)


def stability(runs: dict[int, np.ndarray]) -> dict:
    seeds = sorted(runs)
    pairs = [(a, b) for i, a in enumerate(seeds) for b in seeds[i + 1:]]
    ari = {f"{a}v{b}": round(float(adjusted_rand_score(runs[a], runs[b])), 4) for a, b in pairs}
    v = list(ari.values())
    return {"pairs": ari, "mean_ari": round(float(np.mean(v)), 4) if v else None,
            "min_ari": round(float(np.min(v)), 4) if v else None,
            "k_per_seed": {str(s): int(runs[s].max() + 1) for s in seeds}}


# --------------------------------------------------------------------------- one method

def evaluate(name: str, mod, prep: dict, args, ind: np.ndarray | None) -> tuple[dict, np.ndarray]:
    X, Z, df = prep["X"], prep["Z"], prep["df"]
    xy = df[["X_centroid", "Y_centroid"]].to_numpy(float)
    # which implementation was actually scored: these modules are being written while this runs,
    # and a result is worth nothing without the file it came from
    out: dict = {"method": name, "module": getattr(mod, "__file__", None),
                 "params": dict(getattr(mod, "PARAMS", {})),
                 "declared_name": getattr(mod, "NAME", None), "criteria_unavailable": {}}

    runs, timing = {}, {}
    for s in args.seeds:
        t0 = time.time()
        runs[s] = check_labels(mod.fit_predict(X, Z, seed=s), len(X))
        timing[str(s)] = round(time.time() - t0, 1)
    labels = runs[args.seeds[0]]
    k = int(labels.max() + 1)
    sizes = np.bincount(labels, minlength=k)
    out.update({"seconds_per_seed": timing, "k": k,
                "largest_cluster_share": round(float(sizes.max() / len(labels)), 4),
                "n_clusters_under_10_cells": int((sizes < 10).sum()),
                "median_cluster_size": int(np.median(sizes))})

    out["islet"] = marked_recovery(labels, prep["core"], k, args.perm_draws, args.seed)
    out["islet_composition"] = composition(labels, prep["core"], k)
    out["rare"] = marked_recovery(labels, prep["rare"], k, args.perm_draws, args.seed)
    out["rare_composition"] = composition(labels, prep["rare"], k)

    M = cluster_means(Z, labels, k)
    out["heldout"] = heldout_agreement(M, prep["markers"], sizes, args.splits, args.seed)
    if "error" in out["heldout"]:
        out["criteria_unavailable"]["heldout_excess"] = out["heldout"]["error"]

    # Only the ambiguous row goes through call.py's naming policy.  A partition that policy cannot
    # afford to name loses that one row and keeps the other six, rather than dropping out of the
    # comparison: not being nameable at a reviewable cost is a finding about the method, not a
    # reason to hide it.
    named = None
    try:
        named, why, n_draws = name_from_means(M, prep["markers"], args.null_draws, args.margin,
                                              args.fdr, args.seed)
    except NamingTooExpensive as exc:
        out["naming"] = {"error": str(exc)}
        out["criteria_unavailable"]["ambiguous_frac"] = str(exc)
    if named is not None:
        amb = np.array([n == L.AMBIGUOUS for n in named])
        reasons: dict[str, int] = {}
        for c in np.flatnonzero(amb):
            reasons[why[c]] = reasons.get(why[c], 0) + int(sizes[c])
        out["naming"] = {"ambiguous_clusters": int(amb.sum()),
                         "ambiguous_cells": int(sizes[amb].sum()),
                         "ambiguous_frac": round(float(sizes[amb].sum() / len(labels)), 4),
                         "naming_draws": n_draws, "reasons_by_cells": reasons,
                         "lineage_cells": {lin: int(sizes[[n == lin for n in named]].sum())
                                           for lin in dict.fromkeys(named)}}

    coh, ind = spatial_coherence(xy, labels, args.coherence_draws, args.seed, ind)
    out["spatial"] = coh
    if named is not None:
        code = pd.factorize(pd.Series(named)[labels].to_numpy())[0].astype(np.int32)
        out["spatial_lineage"], _ = spatial_coherence(xy, code, args.coherence_draws,
                                                      args.seed, ind)
    out["stability"] = stability(runs)

    out["criteria"] = {"islet_purity": out["islet"]["purity_90"],
                       "islet_concentration": out["islet"]["concentration"],
                       "rare_f1": out["rare"]["best_f1"],
                       "heldout_excess": out["heldout"].get("excess_cellweighted"),
                       "spatial_enrichment": out["spatial"]["enrichment"],
                       "ambiguous_frac": out["naming"].get("ambiguous_frac"),
                       "stability_ari": out["stability"]["mean_ari"]}
    return out, ind


# --------------------------------------------------------------------------- reporting

def _jsonable(o):
    """NaN and inf are not JSON.  A metric that could not be computed is written as null."""
    if isinstance(o, dict):
        return {k: _jsonable(v) for k, v in o.items()}
    if isinstance(o, (list, tuple)):
        return [_jsonable(v) for v in o]
    if isinstance(o, (np.integer,)):
        return int(o)
    if isinstance(o, (float, np.floating)):
        return float(o) if np.isfinite(o) else None
    if isinstance(o, (np.bool_,)):
        return bool(o)
    if isinstance(o, pd.DataFrame):
        return o.to_dict("records")
    return o


def rank(results: dict[str, dict]) -> tuple[dict, dict, list[str]]:
    """Rank per criterion (1 = best), and say which criteria could not be computed at all.

    A criterion no method produced a value for is dropped from the mean rank rather than ranked
    arbitrarily -- on the first run of this harness call.py named nothing on any partition, which
    silently handed two of the seven rows to whichever method happened to be listed first.  A
    criterion missing for *some* methods still counts, and the ones missing it rank last on it.
    """
    ok = {m: r for m, r in results.items() if "criteria" in r}
    ranks: dict[str, dict[str, float | None]] = {}
    dropped: list[str] = []
    for key, _desc, direction in CRITERIA:
        vals = {m: r["criteria"].get(key) for m, r in ok.items()}
        vals = {m: (None if v is None or not np.isfinite(v) else float(v) * -direction)
                for m, v in vals.items()}
        have = {m: v for m, v in vals.items() if v is not None}
        if not have:
            ranks[key] = {m: None for m in ok}
            dropped.append(key)
            continue
        order = sorted(have, key=lambda m: have[m])
        r: dict[str, float | None] = {}
        for m in order:
            ties = [x for x in order if have[x] == have[m]]      # ties share the better rank
            r[m] = float(min(order.index(x) for x in ties) + 1)
        for m in ok:
            r.setdefault(m, float(len(ok)))
        ranks[key] = r
    counted = [k for k, _d, _s in CRITERIA if k not in dropped]
    mean = {m: round(float(np.mean([ranks[k][m] for k in counted])), 3) for m in ok} \
        if counted else {m: float("nan") for m in ok}
    return ranks, mean, dropped


def _across_donors(results: list[dict], key: str) -> float | None:
    """A criterion averaged over the donors that produced it, or None if none did."""
    v = [r["criteria"][key] for r in results
         if r["criteria"].get(key) is not None and np.isfinite(r["criteria"][key])]
    return float(np.mean(v)) if v else None


def _where(rows: list[dict], n_marked: int) -> str:
    return "  ".join(f"c{r['cluster']}: {r['n_marked']}/{r['cluster_size']} cells "
                     f"({r['purity']*100:.0f}% pure, {r['share_of_marked']*100:.0f}% of them)"
                     for r in rows[:3]) or f"nowhere ({n_marked} cells)"


def report(donor_results: dict, out_path: Path, args) -> dict:
    per_method: dict[str, dict[str, dict]] = {}
    for donor, res in donor_results.items():
        for m, r in res["methods"].items():
            per_method.setdefault(m, {})[donor] = r

    pooled: dict[str, dict] = {}
    for m, byd in per_method.items():
        good = [r for r in byd.values() if "criteria" in r]
        if not good:
            pooled[m] = {"error": "; ".join(str(r.get("error")) for r in byd.values())}
            continue
        pooled[m] = {"criteria": {k: _across_donors(good, k) for k, _d, _s in CRITERIA},
                     "k": [g["k"] for g in good],
                     "seconds": [sum(float(v) for v in g["seconds_per_seed"].values())
                                 for g in good]}
    ranks, mean_rank, dropped = rank(pooled)

    lines = []
    lines.append("\n" + "=" * 100)
    lines.append("CLUSTERING METHOD COMPARISON   "
                 f"donors={list(donor_results)}  seeds={args.seeds}")
    lines.append("=" * 100)

    for donor, res in donor_results.items():
        meta = res["meta"]
        lines.append(f"\n[{donor}] {meta['n_typed']:,} cells x {meta['n_markers']} markers "
                     f"from {meta['source']}   "
                     f"({meta['n_dropped_incomplete']:,} eligible cells had an incomplete "
                     f"measurement and are excluded, never filled)")
        lines.append(f"       islet-core cells {meta['n_core']:,} ({meta['core_share']*100:.2f}% "
                     f"- the null a cluster has to beat on purity)")
        r = meta["rare"]
        lines.append(f"       {meta['rare_marker']}-high cells {r['n']} "
                     f"(gap at rank {r.get('rank_of_gap')}: {r.get('cut'):.3f} -> "
                     f"{r.get('next_value'):.3f})")
        for m, x in res["methods"].items():
            if "criteria" not in x:
                lines.append(f"  {m:<10} SKIPPED: {x.get('error')}")
                continue
            isl, rr, hd, sp, st = (x["islet"], x["rare"], x["heldout"], x["spatial"],
                                   x["stability"])
            lines.append(f"  {m:<10} k={x['k']:<4} largest={x['largest_cluster_share']*100:4.1f}% "
                         f"<10 cells={x['n_clusters_under_10_cells']:<3} "
                         f"{sum(float(v) for v in x['seconds_per_seed'].values()):.0f}s/"
                         f"{len(x['seconds_per_seed'])} seeds")
            lines.append(f"             islet   {isl['clusters_90']} clusters hold 90% "
                         f"(null p05 {isl['null_clusters_90_p05']:.0f} = "
                         f"{isl['concentration']}x tighter than chance), purity "
                         f"{isl['purity_90']:.3f} (null p95 {isl['null_purity_90_p95']:.3f}), "
                         f"best F1 {isl['best_f1']:.3f} (null p95 "
                         f"{isl['null_best_f1_p95']:.3f})")
            lines.append(f"                     where they went -- "
                         f"{_where(x['islet_composition'], isl['n_marked'])}")
            lines.append(f"             rare    best cluster {rr['best_cluster']} "
                         f"n={rr['best_cluster_size']}, purity {rr['best_precision']:.3f}, "
                         f"completeness {rr['best_recall']:.3f}, F1 {rr['best_f1']:.3f} "
                         f"(null p95 {rr['null_best_f1_p95']:.3f}); "
                         f"{rr['clusters_90']} clusters hold 90%")
            lines.append(f"                     where they went -- "
                         f"{_where(x['rare_composition'], rr['n_marked'])}")
            if "error" in hd:
                lines.append(f"             heldout {hd['error']}")
            else:
                lines.append(f"             heldout {hd['n_scored']} cluster-splits scored, "
                             f"{hd['share_cleared']*100:.0f}% clear their null "
                             f"({hd['share_cleared_cellweighted']*100:.0f}% cell-weighted); "
                             f"excess {hd['excess_cellweighted']:+.3f} SD "
                             f"[{hd['excess_p05']:+.3f}, {hd['excess_p95']:+.3f}] "
                             f"(obs {hd['observed_cellweighted']:+.3f} vs null p95 "
                             f"{hd['null_p95_cellweighted']:+.3f}; "
                             f"{hd['exact_null_share']*100:.0f}% of nulls enumerated exactly)")
            lines.append(f"             spatial {sp['observed']:.3f} vs block-shuffled null "
                         f"{sp['null_mean']:.3f} (p95 {sp['null_p95']:.3f}, "
                         f"{sp['n_blocks']} blocks) = {sp['enrichment']}x; lineage-level "
                         f"{x.get('spatial_lineage', {}).get('enrichment')}x"
                         + ("  [few blocks: on a window this small the null is close to a "
                            "global shuffle]" if sp["n_blocks"] < 25 else ""))
            nm = x["naming"]
            lines.append(f"             naming  {nm['error']}" if "error" in nm else
                         f"             naming  {nm['ambiguous_cells']:,} cells "
                         f"({nm['ambiguous_frac']*100:.1f}%) in {nm['ambiguous_clusters']} "
                         f"clusters no lineage claims, over {nm['naming_draws']:,} null draws")
            lines.append(f"             stable  mean ARI {st['mean_ari']} "
                         f"(min {st['min_ari']}, pairs {st['pairs']})")

    methods = [m for m in pooled if "criteria" in pooled[m]]
    if not methods:
        # still write the file: the reasons each backend was skipped are the result of this run
        lines.append("\nNo method produced a usable partition; nothing to recommend.")
        print("\n".join(lines))
        payload = {"generated": time.strftime("%Y-%m-%dT%H:%M:%S"), "seeds": list(args.seeds),
                   "donors": donor_results, "pooled": pooled, "recommendation": None,
                   "reason": "no backend produced a usable partition", "text": "\n".join(lines)}
        out_path.parent.mkdir(parents=True, exist_ok=True)
        out_path.write_text(json.dumps(_jsonable(payload), indent=2, default=str))
        print(f"\nwrote {out_path}")
        return payload

    w = max(len(k) for k, _d, _s in CRITERIA) + 2
    lines.append("\n" + "-" * 100)
    lines.append("PER-CRITERION  (value, then rank in brackets; averaged over donors)")
    lines.append("-" * 100)
    lines.append(f"{'criterion':<{w}}" + "".join(f"{m:>22}" for m in methods) + "   better")
    for key, desc, direction in CRITERIA:
        row = f"{key:<{w}}"
        for m in methods:
            v, rk = pooled[m]["criteria"].get(key), ranks[key][m]
            row += (f"{v:>16.3f} [{int(rk)}]" if v is not None and rk is not None
                    else f"{'not computed':>22}")
        lines.append(row + ("   higher" if direction > 0 else "   lower"))
        lines.append(f"{'':<{w}}{desc}")
    lines.append(f"\n{'mean rank':<{w}}" + "".join(f"{mean_rank[m]:>22.2f}" for m in methods))
    for key in [k for k, _d, _s in CRITERIA if len(methods) > 1 and k not in dropped
                and len({ranks[k][x] for x in methods}) == 1]:
        v = pooled[methods[0]]["criteria"].get(key)
        lines.append(f"  {key}: every method ties at {v}; this row separates nothing here")
    for key in dropped:
        why = sorted({r["criteria_unavailable"].get(key, "") for d in donor_results.values()
                      for r in d["methods"].values() if "criteria_unavailable" in r} - {""})
        lines.append(f"  {key}: no method produced a value, so it is out of the mean rank"
                     + (f" -- {why[0]}" if why else ""))

    best = min(methods, key=lambda m: mean_rank[m])
    if len(methods) == 1:
        # a rank over one method is not a comparison, and printing "leads on everything" would
        # read as though it had beaten something
        lines.append("\n" + "=" * 100)
        lines.append(f"NOT A COMPARISON: {best} is the only method that produced a usable "
                     f"partition. The rows above are its measurements; the other methods are "
                     f"listed with the reason they were skipped. Re-run when they exist.")
        lines.append("=" * 100)
        print("\n".join(lines))
        payload = {"generated": time.strftime("%Y-%m-%dT%H:%M:%S"), "seeds": list(args.seeds),
                   "donors": donor_results, "pooled": pooled, "ranks": ranks,
                   "recommendation": None, "reason": f"only {best} was usable",
                   "text": "\n".join(lines)}
        out_path.parent.mkdir(parents=True, exist_ok=True)
        out_path.write_text(json.dumps(_jsonable(payload), indent=2, default=str))
        print(f"\nwrote {out_path}")
        return payload

    # only where the criterion actually separates the methods: a row every method ties on says
    # nothing about any of them, and reading it as a lead or a trail is worse than not reading it
    split = [k for k, _d, _s in CRITERIA if k not in dropped
             and len({ranks[k][x] for x in methods}) > 1]
    leads = [k for k in split if ranks[k][best] == 1]
    trails = [k for k in split if ranks[k][best] == max(ranks[k][x] for x in methods)]
    lines.append("\n" + "=" * 100)
    lines.append(f"RECOMMENDATION: {best}   (mean rank {mean_rank[best]:.2f} over "
                 f"{len(CRITERIA) - len(dropped)} of {len(CRITERIA)} criteria and "
                 f"{len(methods)} usable methods)")
    lines.append(f"  leads on : {', '.join(leads) if leads else 'no criterion'}")
    lines.append(f"  trails on: {', '.join(trails) if trails else 'no criterion'}")
    tie = [m for m in methods if abs(mean_rank[m] - mean_rank[best]) < 0.5 and m != best]
    if tie:
        lines.append(f"  within half a rank of it: {', '.join(tie)} -- the mean rank does not "
                     f"separate these; read the per-criterion rows and decide on the ones that "
                     f"matter for the question being asked.")
    # the only criterion that carries its own spread; say when a rank on it is inside that spread
    span = {m: [(r["heldout"]["excess_p05"], r["heldout"]["excess_p95"])
                for d in donor_results.values() for n, r in d["methods"].items()
                if n == m and "excess_p05" in r.get("heldout", {})] for m in methods}
    span = {m: (min(a for a, _ in v), max(b for _, b in v)) for m, v in span.items() if v}
    order = sorted(span, key=lambda m: -pooled[m]["criteria"]["heldout_excess"])
    if len(order) > 1 and span[order[0]][0] < span[order[1]][1]:
        lines.append(f"  heldout_excess does not separate {order[0]} and {order[1]} here: "
                     f"resampling the marker splits gives {order[0]} "
                     f"[{span[order[0]][0]:+.3f}, {span[order[0]][1]:+.3f}] and {order[1]} "
                     f"[{span[order[1]][0]:+.3f}, {span[order[1]][1]:+.3f}], which overlap. "
                     f"Its rank still counts once in the mean; read it as a tie.")
    lines.append("  Read ambiguous_frac last, and only once the rest hold: a method that names "
                 "every cluster confidently and wrongly leads that row on its own.")
    lines.append("  None of the seven is decisive by itself and this file does not weight them; "
                 "the recommendation is an unweighted mean rank and the rows above are the "
                 "evidence.")
    lines.append("=" * 100)
    print("\n".join(lines))

    payload = {"generated": time.strftime("%Y-%m-%dT%H:%M:%S"),
               "seeds": list(args.seeds), "criteria": [{"key": k, "description": d,
                                                        "higher_is_better": s > 0}
                                                       for k, d, s in CRITERIA],
               "donors": donor_results, "pooled": pooled, "ranks": ranks,
               "mean_rank": mean_rank, "recommendation": best,
               "leads_on": leads, "trails_on": trails, "text": "\n".join(lines)}
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text(json.dumps(_jsonable(payload), indent=2, default=str))
    print(f"\nwrote {out_path}")
    return payload


# --------------------------------------------------------------------------- cli

def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--donors", nargs="+", required=True)
    ap.add_argument("--methods", nargs="+", default=list(METHODS))
    ap.add_argument("--out", type=Path, default=DEV / "clustering_comparison.json")
    ap.add_argument("--source", choices=("auto", "window", "positivity"), default="auto")
    # every seed re-runs every backend, so this sets the cost: on the test window Leiden is ~10 s
    # a seed and KODAMA is ~4 min, while the harness's own scoring is ~4 s per method whatever n
    ap.add_argument("--seeds", nargs="+", type=int, default=list(SEEDS),
                    help="two or more, or there is no stability to measure")
    ap.add_argument("--seed", type=int, default=0, help="seed for the nulls, not for the methods")
    # 80 random marker splits.  Which three markers name a lineage genuinely changes the answer,
    # so the excess settles slowly: on a 40-cluster partition of the test window it reads +0.189
    # (5 splits), +0.083 (10), +0.068 (20), +0.108 (40), +0.102 (80), +0.090 (160), +0.084 (320).
    # The level is still moving at 80 and this is honest about that -- what is stable is the
    # ordering the criterion is used for: that same partition beats a random one of the same k at
    # every split count tried (+0.083 vs -0.015 at 10, +0.102 vs -0.012 at 80).  80 costs 0.05 s.
    ap.add_argument("--splits", type=int, default=80, help="random naming/held-out marker splits")
    # call.py's own CLI defaults, so the ambiguous fraction reported here is the one call.py
    # would produce on that partition.  --null-draws is a floor: the naming needs
    # clusters x lineages / fdr draws, and is given that many when it is more.
    ap.add_argument("--null-draws", type=int, default=10_000, help="floor; call.py's default")
    ap.add_argument("--margin", type=float, default=0.15, help="call.py's default")
    ap.add_argument("--fdr", type=float, default=0.05, help="call.py's default")
    # both nulls are converged well below their defaults on the test window: the islet
    # concentration reads 3.244 at 20 permutation draws and 3.125 at 50, 100 and 200, and the
    # coherence enrichment reads 1.962 / 1.959 / 1.948 / 1.952 at 5 / 10 / 20 / 50 block shuffles.
    ap.add_argument("--perm-draws", type=int, default=50, help="islet and rare permutation nulls")
    ap.add_argument("--coherence-draws", type=int, default=20, help="block-shuffle draws")
    ap.add_argument("--max-cells", type=int, default=0, help="subsample; 0 = every cell")
    a = ap.parse_args()

    backends = {}
    for m in a.methods:
        mod, err = load_backend(m)
        backends[m] = (mod, err)
        print(f"backend {m:<10} " + (f"OK  {mod.__file__}" if mod else f"SKIPPED - {err}"))

    donor_results = {}
    for d in a.donors:
        try:
            prep = prepare(d, a.source, a.max_cells, a.seed)
        except Exception as exc:                   # a donor with no positivity yet must not take
            print(f"[{d}] SKIPPED - {type(exc).__name__}: {exc}")   # the donors that do have it
            continue
        res: dict = {"meta": prep["meta"], "methods": {}}
        ind = None
        for m, (mod, err) in backends.items():
            if mod is None:
                res["methods"][m] = {"method": m, "error": err}
                continue
            try:
                t0 = time.time()
                r, ind = evaluate(m, mod, prep, a, ind)
                r["total_seconds"] = round(time.time() - t0, 1)
                res["methods"][m] = r
                print(f"[{d}] {m}: k={r['k']} in {r['total_seconds']:.0f}s")
            except Exception as exc:                       # one bad backend must not end the run
                tb = traceback.format_exc()
                # some libraries raise with no message at all, and "NotImplementedError:" on its
                # own tells a reader nothing; carry the line it came from
                msg = (f"{type(exc).__name__}: {exc}" if str(exc) else
                       f"{type(exc).__name__} (no message) at {tb.splitlines()[-2].strip()}")
                res["methods"][m] = {"method": m, "error": msg, "traceback": tb[-2000:]}
                print(f"[{d}] {m}: FAILED - {msg}")
        donor_results[d] = res
    if not donor_results:
        raise SystemExit("no donor had both a substrate and a positivity table; "
                         "nothing to compare")
    report(donor_results, a.out, a)


if __name__ == "__main__":
    main()
