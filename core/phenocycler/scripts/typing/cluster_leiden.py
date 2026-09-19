#!/usr/bin/env python3
"""Leiden communities on a k-nearest-neighbour graph of the Nimbus confidences.

RETIRED -- NOT A FALLBACK.  Recorded here 2026-09-18 because the docstring below reads as a live
specification and is not one.  This module belongs to the broad-typing producer comparison
WITHDRAWN 2026-09-08 (``cluster_leiden.py`` / ``cluster_lda.py`` / ``cluster_kodama.py`` scored by
``compare_clustering.py``).  It did not produce the calls on disk, and the current typing is the
subpanel classification in ``subpanel.py``; see ``README.md``, "Retired, and not to be revived as
fallbacks".  What follows is kept as a record of what was measured and why.

One of the three clustering methods ``compare_clustering.py`` scores against each other.  Leiden
is the one that is not told how many groups to look for: it reads local density on a kNN graph,
so an exocrine community of several thousand cells and a delta-cell community of fifty can both
survive the same partition, and ``resolution`` sets granularity rather than a hard count.  That is
the property this study needs -- its question is about a minority population near islets, so a
method that cannot hold a rare community is failing at the thing it is for.

Every number quoted below was measured on the validated window of donor 6442: 11,659
analysis-eligible cells, 34 identity markers, 638 of them inside a QuPath islet polygon and 41
carrying high SST.  The polygons are hand-drawn geometry that no marker, no model and no
clustering ever sees, so they are the check and not an input.  The SST-high cells are the rarest
thing this panel can name -- 0.35% of the window -- and are the population every parameter here
was watched against.

Nothing in this file compares a confidence to a cutoff.  The confidences are coordinates.

At whole-slide size, measured.  Nothing here subsamples -- every eligible cell is a vertex.  On
this machine (64 cores, other jobs running, so these are pessimistic): 1.2 s for the 11,659-cell
window, 4.4 s at 100k, 18.5 s at 400k, 90 s at 1.4M (37 s kNN + 53 s Leiden, 12.8 GiB peak RSS)
and 154 s at 1.9M (16.7 GiB).  Roughly linear -- 44 us/cell at 100k against 81 us/cell at 1.9M --
and the memory is the graph: 39M connectivity entries at 1.4M cells, which scanpy turns into a
Python list of edge tuples before igraph sees it.  Add ~25 s once per process for numba to compile
pynndescent.  Those figures come from the test window resampled to size with 0.10 SD of jitter,
because the whole-slide positivity table is still being produced; the timing and the memory carry
over, the CLUSTER COUNT does not -- the jitter manufactures ~11,500 tight balls of near-duplicates
and Leiden dutifully finds them.  For k, the real subsamples are the guide: 33, 39, 44, 48 clusters
at 1,457, 2,914, 5,829 and 11,659 real cells, about five more per doubling, which extrapolates to
80-90 clusters on a 1.4M-cell donor.  Leiden's resolution is not scale-free, so that is a
projection and not a promise; check k on the first whole-slide run.
"""

from __future__ import annotations

import numpy as np

NAME = "leiden"

#: Defaults, each chosen from a property of the data that is not a composition.  ``fit_predict``
#: reads this dict for its own defaults, so a caller that passes nothing gets what was measured.
PARAMS = {
    #: kNN degree.  Recovery of the 41 SST-high cells is not what separates these: every
    #: neighbourhood from 10 to 50 puts all 41 in one community.  What varies is whether that
    #: community is theirs -- purity 0.759 at 10, 15 and 20, 0.788 at 30, 0.891 at 50 -- and how
    #: the rest of the partition behaves.  10 is out: at this resolution it already emits an
    #: 8-cell cluster (smallest cluster 8 against 45 at 20), and call.py names a cluster from its
    #: mean profile over 34 markers, which 8 cells cannot support.  50 is out for a reason the
    #: coarser resolutions expose.  (An earlier note here said a 50-neighbour ball at resolution
    #: 2.0 recovers only 37 of the 41 delta cells.  Independent re-measurement did NOT reproduce
    #: that: at the shipped n_jobs it is 41/41 at all three seeds, worst case 40/41 at n_jobs=1.
    #: The partition is thread-dependent, so the original figure was probably taken at an
    #: unrecorded thread count.  The exclusion of 50 stands on stability alone -- see below --
    #: not on that number.)  Among what is left, stability decides:
    #: mean ARI across seeds 0.638 at 20 against 0.573 (15), 0.624 (30), and the cell-weighted
    #: bootstrap Jaccard is 0.533 +- 0.028 at 20 against 0.509 +- 0.024 at 15.
    "n_neighbors": 20,
    #: Granularity.  Two criteria that pull opposite ways, neither of them a composition.
    #: Held-out marker agreement -- cluster on 28 markers, predict the 6 withheld, scored on cells
    #: withheld as well -- rises with resolution and always will: 0.359, 0.376, 0.387, 0.395,
    #: 0.414 at 1.0, 2.0, 3.0, 4.0, 6.0, because more clusters means more free means.  What bounds
    #: it is reproducibility, measured as Hennig's clusterwise bootstrap Jaccard.  THE LEVELS
    #: FIRST CLAIMED HERE DID NOT SURVIVE INDEPENDENT RE-MEASUREMENT and the rule built on them
    #: is withdrawn.  Two implementations, 10 draws each, give 0.450 +- 0.022 / 0.449 +- 0.011 /
    #: 0.438 +- 0.013 / 0.412 +- 0.014 (counting multiplicities) and 0.454 / 0.452 / 0.441 /
    #: 0.413 (fpc::clusterboot's unique-point convention) at 1.5, 2.0, 3.0, 4.0 -- 0.09-0.14
    #: below the figures first recorded.  The ORDERING reproduces (monotone decreasing) but NO
    #: resolution reaches Hennig's 0.5, so "3.0 is the last resolution a full SD clear of the
    #: level" is false, and 1.5-2.0 are marginally MORE stable than 3.0 rather than less.
    #:
    #: What this does and does not settle.  It does not make 3.0 wrong: held-out marker agreement
    #: rises monotonically with resolution and stability falls monotonically, so any value in
    #: 1.5-3.0 is defensible and none is picked out by the data.  It does mean the choice is a
    #: judgement inside a shallow band, not a landmark -- so record the value with the run.  In
    #: this lane the final granularity is set by the merge stage (scripts/typing/metacluster.py)
    #: rather than by this resolution, which is why a weakly-determined fine resolution is
    #: tolerable here and would not be in a one-stage design.
    "resolution": 3.0,
    #: The modern scanpy path.  ``flavor="igraph"`` calls igraph's own C Leiden with the modularity
    #: objective; 2 refinement iterations is what scanpy pairs with it.  Confirmed rather than
    #: assumed: running to convergence (-1) buys 0.8% of modularity (0.6819 -> 0.6874) for 1.8x
    #: the Leiden time, leaves k at 48 and the SST-high recovery at 41/41 either way.  2 is taken
    #: because the extra iterations change nothing that is being measured.
    "n_iterations": 2,
    #: Threads for the kNN build.  scanpy's own default is ``settings.n_jobs = 1`` and it offers
    #: no way to change it except a global, so the transformer is built here instead.  At 11,659
    #: cells this buys nothing -- 1.4 s on one thread, 1.4 s on 64 -- and it is put in for the
    #: whole slide, where the kNN build is the run.  -1 is also the only value pynndescent does
    #: NOT push into ``numba.set_num_threads``, so it is the only one that leaves the process's
    #: thread count alone for whatever the comparison runs next.  The partition does depend on it:
    #: labels are bit-identical run to run at a fixed seed and thread count, but ARI between the
    #: 1-thread build and 8, 32 or 64 threads is 0.62-0.67, so the thread count belongs in the run
    #: record next to the seed.  What the threads do not move is the biology -- the SST-high cells
    #: come back 41/41 in a community of 75.9% purity at every thread count tried.
    "n_jobs": -1,
}


def fit_predict(X: np.ndarray, Z: np.ndarray, *, seed: int = 0, **params) -> np.ndarray:
    """X: raw confidences, cells x markers, non-negative, in [0,1], no NaN.
       Z: the same matrix standardised per marker (zero mean, unit variance per column).
       returns: int32 array, one cluster label per cell, labels 0..k-1, no gaps.

    Uses **Z**.  A kNN graph is euclidean distance, so whichever channel spans the widest range
    decides who is whose neighbour, and on raw X that is the abundant ones: Pan-Cytokeratin runs
    0.02 to 0.78 across this window (SD 0.256) while SST runs 0.02 to 0.16 for 99% of cells
    (SD 0.047) and only reaches 0.8 on the delta cells themselves.  Built on X, all 41 SST-high
    cells do end up together -- inside a 409-cell mixed endocrine community, 10.0% of it, which is
    not a delta-cell cluster but an islet one.  Standardising puts SST on the same footing as the
    rest (a delta cell then sits 14.4 SD out on it) and the same 41 cells hold a community of
    their own: 54 cells, 75.9% purity, and all 54 of them inside a QuPath islet polygon.  It is
    also what the lane already does for a separate reason -- Nimbus scores TP63 above its decision
    boundary in 98.6% of cells and FOXP3 in 90.4%, and standardising makes those channels
    contribute their variation rather than their level.

    No PCA.  34 dimensions is already narrow, and the unprojected graph is the only
    representation whose rare community survives the whole resolution range: purity 0.759, 0.774,
    0.759, 0.759 at resolutions 1.5, 2.0, 3.0, 4.0.  Every projection tried collapses somewhere --
    10 components give 0.076, 0.098, 0.889, 0.263 over those same resolutions, and even 30
    components drop to 0.387 at 1.5.  The reason is visible before any clustering: 10 components
    keep 71.5% of the total variance but only 11.1% of SST's, the worst-retained channel of the
    34, which is what a variance ranking does to a marker carried by 0.35% of the cells.  It takes
    30 of the 34 components to keep every channel, at which point the projection is a rotation --
    nothing gained, and a rare community to lose.
    """
    import anndata as ad
    import scanpy as sc
    from pynndescent import PyNNDescentTransformer

    p = dict(PARAMS)
    unknown = sorted(set(params) - set(p))
    if unknown:                       # a typo silently ignored is a parameter silently unset
        raise TypeError(f"cluster_leiden.fit_predict got unknown parameter(s) {unknown}; "
                        f"known: {sorted(p)}")
    p.update(params)

    F = np.ascontiguousarray(Z, dtype=np.float32)
    if F.ndim != 2:
        raise ValueError(f"Z must be cells x markers, got shape {F.shape}")
    n = F.shape[0]
    if np.asarray(X).shape[0] != n:
        raise ValueError(f"X has {np.asarray(X).shape[0]} rows, Z has {n}")
    if not np.isfinite(F).all():
        raise ValueError("Z carries NaN or inf.  A missing measurement is unavailable, not zero: "
                         "drop the cell upstream and record the reason, do not impute it here")
    nn = int(p["n_neighbors"])
    if not 2 <= nn < n:
        raise ValueError(f"n_neighbors must be in [2, {n}), got {nn}")

    # scanpy passes `settings.n_jobs` (1 by default) into its own transformer and offers no other
    # way to change it, so the transformer is built here instead.  n_trees and n_iters repeat
    # scanpy's UMAP-derived formulas verbatim: at n_jobs=1 this reproduces the stock path exactly
    # (ARI 1.000 on the test window), so what is gained is the thread count and nothing else.
    knn = PyNNDescentTransformer(
        n_neighbors=nn, metric="euclidean", random_state=seed, n_jobs=p["n_jobs"],
        n_trees=min(64, 5 + round(n ** 0.5 / 20.0)), n_iters=max(5, round(np.log2(n))))

    a = ad.AnnData(F)
    sc.pp.neighbors(a, use_rep="X", transformer=knn, random_state=seed)
    sc.tl.leiden(a, resolution=float(p["resolution"]), key_added="lc", random_state=seed,
                 flavor="igraph", n_iterations=int(p["n_iterations"]), directed=False)
    raw = a.obs["lc"].to_numpy().astype(np.int64)

    # Largest community first, so cluster 0 means the same kind of thing between runs and an
    # empty id can never leave a gap in the labels the caller checks for.
    counts = np.bincount(raw)
    order = np.argsort(-counts, kind="stable")
    order = order[counts[order] > 0]
    remap = np.full(counts.size, -1, np.int32)
    remap[order] = np.arange(order.size, dtype=np.int32)
    out = remap[raw]
    if (out < 0).any():
        raise RuntimeError("relabelling lost a cell")
    return out
