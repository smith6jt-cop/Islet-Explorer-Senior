#!/usr/bin/env python3
"""Cluster the cells on the graph KODAMA learns, never on a projection of it.

RETIRED -- NOT A FALLBACK.  Recorded here 2026-09-18 because the docstring below reads as a live
specification and is not one.  This module belongs to the broad-typing producer comparison
WITHDRAWN 2026-09-08 (``cluster_leiden.py`` / ``cluster_lda.py`` / ``cluster_kodama.py`` scored by
``compare_clustering.py``).  It did not produce the calls on disk, and the current typing is the
subpanel classification in ``subpanel.py``; see ``README.md``, "Retired, and not to be revived as
fallbacks".  What follows is kept as a record of what was measured and why.

KODAMA is a representation learner, not a clusterer.  ``KODAMA.matrix`` runs M Monte-Carlo
iterations -- each draws ``landmarks`` landmark cells, trains a cross-validated PLS-DA on them
and projects every other cell onto it, giving one label per cell -- and then rewrites its own
k-nearest-neighbour arrays as

    d_KODAMA(i, j) = (1 + d_euclidean(i, j)) / p_ij ** 2

where ``p_ij`` is the fraction of the M iterations that put i and j in the same class, and
re-sorts each row by it.  Two cells the iterations keep separating are pushed apart however close
they started; two they keep together are pulled in.  That rewritten array *is* the KODAMA
representation: ``KODAMA.visualization`` does nothing else with the fit, it hands those same two
arrays to UMAP.

**So this module clusters the graph, and builds no embedding at all.**  Measured on the 6442
test window (11,659 eligible cells, 34 markers, M=100, landmarks=500, and the SAME Leiden --
degree 15, resolution 2.0 -- on all three), seeds 0 and 1, with both KODAMA's seed and Leiden's
moving together:

    representation          clusters     seed-to-seed ARI   islet cells spread over
    KODAMA kNN graph        31, 33              0.761        6 clusters
    KODAMA UMAP, 10 dims    90, 88              0.681       10 clusters
    KODAMA UMAP,  2 dims   102, 102             0.556       10 clusters

Laying the graph out first shatters the section into three times as many communities and makes
the answer less reproducible; ten dimensions repair only part of what two destroy.  The
higher-dimensional representation the task asked about does exist -- ``config.umap.default`` takes
``n_components`` -- but the graph beneath it is better still, and skipping the layout also saves
the 30 s it costs.  Nothing is thrown away by stopping at the graph.

What the metric buys, on the same window and the same Leiden, over seeds 0/1/2: mean pairwise ARI
0.699 on the KODAMA graph (33-35 clusters) against 0.593 on the Euclidean one (37-40 clusters).
Both reach purity 1.00 against the islet geometry, so what KODAMA adds here is reproducibility,
not resolution of the islet -- for about ten times the wall clock.

**The clusterer is Leiden at the degree and resolution ``call.py`` already runs.**  KODAMA must be
given a clusterer, and giving it a different one would measure the clusterer rather than the
representation.  With ``graph_k = 15`` and ``resolution = 2.0`` fixed on both sides, the only
thing that differs between this module and the Leiden module is the metric the graph is built in,
which is exactly the thing KODAMA contributes.

Parameters, and the non-composition evidence each was chosen from -- all on the test window, three
seeds each, and all measured as *stability of the representation*: the mean per-cell Jaccard
overlap of a cell's 15 nearest KODAMA neighbours between two independent seeds.

* ``M = 100``.  Jaccard 0.860 (M=25) -> 0.941 (50) -> 0.972 (100) -> 0.981 (200).  ``p_ij`` is a
  mean of M Bernoulli draws, so it converges like 1/sqrt(M) and the curve says so: doubling 50 to
  100 buys 0.031, doubling 100 to 200 buys 0.009 for twice the wall clock.  100 is the knee.
* ``landmarks = 500``.  Jaccard 0.947 (200) -> 0.972 (500) -> 0.973 (1000) -> 0.962 (2000): flat
  from 500 up, and 500 is the cheapest point on that flat (160 s, against 237 s at 1000 and 204 s
  at 2000, M=100 throughout).  It is also exactly where KODAMA's own neighbour count saturates --
  ``neighbors = min(landmarks, 0.75n - 1, 500)`` -- so above 500 the landmarks buy a finer
  per-iteration partition but not a wider graph, and below it the graph itself narrows: at
  landmarks=200 the run returns only a 200-neighbour graph to choose the 15 from.
* ``ncomp = 50`` is not a choice here.  KODAMA's ``quality_control`` caps the PLS-DA component
  count at ``ncol(data)``, so on a 34-marker panel every run uses all 34 and the cap never binds.
  It is passed as 50 only so the argument keeps its meaning if the panel ever widens.
* ``graph_k = 15`` and ``resolution = 2.0`` are ``call.py``'s own defaults, taken to keep the
  comparison controlled (above).  15 is not a fragile point, and the trade it sits on is visible:
  mean seed-to-seed cluster ARI 0.668 / 0.736 / 0.791 / 0.834 at degree 10 / 15 / 30 / 50, rising
  because more edges smooth the graph, against granularity falling the same way (37 / 32 / 26 / 24
  clusters).  There is no interior optimum to find, which is why this one is taken from the
  comparison rather than tuned.
* ``n_cores = 5``.  KODAMA's M loop dispatches chunks of ``ceiling(M / (4 * n.cores))`` to
  ``mclapply``, so it cannot fill more than that many workers: 5 cores is the most M=100 can
  occupy (it needs M >= 4 * cores^2).  Asking for more is measurably slower, not faster --
  recorded in ``scripts/xenium_celltype/hpg/kodama_shard.R``, which is where this gate comes from,
  and re-enforced in the R half.  The 64-core box cannot be filled by one KODAMA run.

Cost is driven by input WIDTH, and the panel is narrow.  ``run_kodama.R`` measured 5,101 genes at
1583.2 s against 44.7 s for 50 principal components at the same n / landmarks / M -- 35x for width
alone, which is 316 s per Monte-Carlo iteration on 5,000 cells.  Confirmed here at 34 markers: one
iteration costs 2.91 s on 11,659 cells (single core, from M=4 at 51.8 s and M=12 at 75.1 s), two
orders of magnitude less on twice the cells.  Width is not what stops this method at cohort scale;
the per-cell R loops are.

Uses ``Z``, the per-marker standardised matrix, not ``X``.  Every step inside KODAMA that touches
the data is Euclidean -- the landmark draw is ``kmeans(data, landmarks)``, the neighbour search is
``Rnanoflann::nn(..., "euclidean")``, and PLS-DA works on the raw spread -- so on ``X`` a marker's
influence would be its own standard deviation, which on these confidences ranges from 0.047 to
0.270, a 5.8x span.  Standardising is the same decision ``call.py`` makes for the same reason: a
channel that is broadly positive (Nimbus scores TP63 above its boundary in 98.6% of cells) should
contribute its variation and not its level.  Run both ways on the test window, ``Z`` is the more
reproducible representation -- seed-to-seed neighbour agreement 0.989 against 0.939 -- and the
cheaper fit, 160 s against 206 s.  Islet purity is 1.00 either way, so the choice rests on
stability and cost, not on what the clusters turned out to contain.  The confidences are never
thresholded anywhere in this module.

What a whole slide will cost, measured rather than assumed.  Nothing in this module is size
limited -- it streams one parquet each way and never materialises an n x n object -- but KODAMA's
own fixed cost is a 501-neighbour Euclidean search over every cell, and a kd-tree in 34 dimensions
is well past the point where it helps.  ``Rnanoflann::nn(Z, Z, 501)`` on 16 cores, measured:

    cells      11,659    50,000   150,000   400,000
    seconds       1.2      17.0     225.3    1971.2
    exponent            1.82      2.36      2.21     (over the preceding step)

That is n^2.2, not n log n, and it is the term that decides the answer.  Extrapolating the last
step to 1.4 M cells gives ~8.6 h on 16 cores, and the module runs it at ``n_cores = 5`` (the M
loop's gate), which measured 4x rather than 10.6x parallel speed-up: call it 22-35 h.  The
Monte-Carlo loop is nearly linear by comparison -- 2.91 s per iteration at 11,659 cells, 14.2 s at
50,000, n^1.09 -- so ~535 s per iteration at 1.4 M, about 3 h for M=100 at 5-way occupancy, and it
is the ONLY one of the three that shards across nodes (``hpg/kodama_shard.R`` is exact).  The
dissimilarity rewrite is an R-level loop over n cells x 500 neighbours x M draws -- 7e10
comparisons at 1.4 M -- but it is linear and adds under an hour.  Total: **a day to a day and a
half per donor, two to three days for the pair**, against minutes for the same Leiden on the
Euclidean graph.  Memory is not the limit: the two neighbour arrays are 1.4e6 x 501 x 8 B x 2 =
10.5 GiB (2.99 GiB measured at 400,000 cells, RSS 3.3 GB), transiently doubled by KODAMA's own
``[, -1]`` copies, plus 2 x 1.1 GiB for ``res`` and ``res_constrain`` at M=100 -- budget 40 GB
against the box's 503.

Two things about these parameters do NOT transfer unchanged to a whole slide, and they are
properties of the method, not of the settings.  KODAMA never adds an edge: it re-weights and
re-sorts the Euclidean 500-nearest-neighbour list and nothing else, so a cell whose real relatives
lie outside its 500 Euclidean neighbours can never be reconnected to them.  On this window 500
candidates are 4.3% of the section; on 1.4 M cells they are 0.036%, so the metric has far less
room to correct anything.  And ``landmarks = 500`` is 4.3% of this window but 0.036% of a slide,
which means a rare population -- a few hundred delta cells in 1.4 M -- would appear in essentially
no landmark set and be projected rather than modelled.  Both argue for landmarks scaling with n,
which makes the cost below worse, not better.  Neither can be fixed from outside the package.

The one lever inside the package does not work.  Dropping ``landmarks`` narrows ``neighbors`` in
proportion and so cuts both the search and the rewrite, but it raises the per-iteration cost:
measured at 50,000 cells, landmarks=100 is 306 s + 18.4 s per iteration against landmarks=500's
578 s + 14.2 s, which crosses over well before M=100 and ends up slower there (2,146 s against
1,998 s).  So the honest position is that this method runs unchanged at cohort scale and is not
cheap enough to run casually; it is worth its day when the question is whether a structure
survives the metric, and not worth it as the routine path.

The R half is ``cluster_kodama.R`` beside this file; they exchange parquet in a scratch directory,
which is the contract ``scripts/broadtype/arms/celesta.py`` and ``tacit.py`` already use for R
methods.  Set ``TYPING_KEEP_SCRATCH=1`` to keep it.
"""

from __future__ import annotations

import json
import os
import shutil
import subprocess
import tempfile
from pathlib import Path

import numpy as np

NAME = "kodama"

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[3]
R_SCRIPT = HERE / "cluster_kodama.R"
SCRATCH_ROOT = REPO / "data/development/typing/scratch/kodama"

#: R 4.5.1 with the user library that holds KODAMA.  ``/usr/bin/Rscript`` on this machine is
#: 4.6.1 and cannot see it, so the path is pinned rather than resolved from PATH.
RSCRIPT = "/usr/local/bin/Rscript"

#: The install this module was written against.  An unmodified fork of tkcaccia/KODAMAlegacy.
KODAMA_REMOTE = "smith6jt-cop/KODAMA"
KODAMA_VERSION = "3.3"
KODAMA_SHA = "d2dbe30ee66509b82c084616d6961c5b292cd059"

PARAMS = {
    # Monte-Carlo iterations. p_ij is a mean of M Bernoulli draws; seed-to-seed neighbour
    # agreement is 0.860 / 0.941 / 0.972 / 0.981 at M = 25 / 50 / 100 / 200, so 100 is the knee.
    "M": 100,
    # Landmark cells per iteration. Neighbour agreement 0.947 / 0.972 / 0.973 / 0.962 at
    # 200 / 500 / 1000 / 2000, and 500 is where KODAMA's own neighbour count hits its 500 cap.
    "landmarks": 500,
    # PLS-DA components. Capped to ncol(data) by KODAMA itself, so on a 34-marker panel this is
    # always 34; carried at 50 so the argument survives a wider panel.
    "ncomp": 50,
    # kNN graph degree handed to Leiden -- call.py's --neighbours, so the comparison isolates
    # the metric. Seed-to-seed cluster ARI 0.668 / 0.736 / 0.791 at degree 10 / 15 / 30, against
    # 37 / 32 / 26 clusters: stability and granularity trade monotonically, with no interior best.
    "graph_k": 15,
    # Leiden granularity -- call.py's --resolution, for the same reason.
    "resolution": 2.0,
    # KODAMA's M loop cannot occupy more than ceiling(M / (4 * n_cores)) workers, so M=100 fills
    # exactly 5. The R half refuses to start if this pair leaves cores idle.
    "n_cores": 5,
}


def _run_r(scratch: Path, seed: int, p: dict) -> dict:
    """Fit KODAMA in R and read back its kNN graph; a failure raises with R's own message."""
    cmd = [RSCRIPT, str(R_SCRIPT), "--scratch", str(scratch),
           "--M", str(int(p["M"])), "--landmarks", str(int(p["landmarks"])),
           "--ncomp", str(int(p["ncomp"])), "--graph-k", str(int(p["graph_k"])),
           "--n-cores", str(int(p["n_cores"])), "--seed", str(int(seed))]
    # R links a threaded BLAS and KODAMA's inner products are small, so the default thread pool
    # spends its time contending: 47.3 s -> 15.9 s for the same fit once it is pinned to one.
    env = {**os.environ, "OPENBLAS_NUM_THREADS": "1", "OMP_NUM_THREADS": "1",
           "MKL_NUM_THREADS": "1"}
    r = subprocess.run(cmd, capture_output=True, text=True, env=env)
    if r.returncode != 0:
        raise RuntimeError(f"{R_SCRIPT.name} failed (rc={r.returncode}):\n{r.stderr[-4000:]}")
    meta = json.loads((scratch / "kodama.json").read_text())
    if meta["kodama_version"] != KODAMA_VERSION or meta["kodama_sha"] != KODAMA_SHA:
        raise RuntimeError(f"R ran KODAMA {meta['kodama_version']} {meta['kodama_sha']}, "
                           f"this module is pinned to {KODAMA_VERSION} {KODAMA_SHA} "
                           f"({KODAMA_REMOTE})")
    return meta


def fit_predict(X: np.ndarray, Z: np.ndarray, *, seed: int = 0, **params) -> np.ndarray:
    """X: raw confidences, cells x markers, non-negative, in [0, 1], no NaN.
       Z: the same matrix standardised per marker (zero mean, unit variance per column).
       returns: int32 array, one cluster label per cell, labels 0..k-1, no gaps.

    Uses ``Z``.  Every step inside KODAMA is Euclidean -- the landmark draw is ``kmeans``, the
    neighbour search is ``Rnanoflann::nn(..., "euclidean")``, PLS-DA reads raw spread -- so on
    ``X`` each marker's weight would be its own standard deviation, 0.047 to 0.270 on these
    confidences.  ``X`` is accepted, and only used to check that the two agree in shape.
    """
    import anndata as ad
    import pyarrow as pa
    import pyarrow.parquet as pq
    import scanpy as sc
    import scipy.sparse as sp

    p = dict(PARAMS)
    unknown = set(params) - set(p)
    if unknown:
        raise TypeError(f"unknown parameter(s) for {NAME}: {sorted(unknown)}")
    p.update(params)

    Z = np.ascontiguousarray(np.asarray(Z, dtype=np.float64))
    if Z.ndim != 2:
        raise ValueError(f"Z must be cells x markers, got shape {Z.shape}")
    if np.asarray(X).shape != Z.shape:
        raise ValueError(f"X {np.asarray(X).shape} and Z {Z.shape} disagree in shape")
    if not np.isfinite(Z).all():
        # KODAMA stops on NA itself; stopping here says which column, and a missing measurement
        # must never arrive as an imputed value.
        bad = np.flatnonzero(~np.isfinite(Z).all(axis=0))
        raise ValueError(f"Z has non-finite values in column(s) {bad.tolist()}; a missing "
                         f"measurement is unavailable, not a number")
    n = Z.shape[0]

    SCRATCH_ROOT.mkdir(parents=True, exist_ok=True)
    scratch = Path(tempfile.mkdtemp(prefix=f"{NAME}_{os.getpid()}_", dir=SCRATCH_ROOT))
    try:
        pq.write_table(pa.table({f"m{i:02d}": Z[:, i] for i in range(Z.shape[1])}),
                       scratch / "Z.parquet", compression="zstd")
        meta = _run_r(scratch, seed, p)
        k = int(p["graph_k"])
        # Column by column into numpy, not through Python objects: at 1.4 M cells the graph is
        # 42 M values and to_pydict() would box every one of them.
        knn = pq.read_table(scratch / "knn.parquet")
        idx = np.column_stack([knn.column(f"i{j + 1:02d}").to_numpy(zero_copy_only=False)
                               for j in range(k)]).astype(np.int64) - 1
        dst = np.column_stack([knn.column(f"d{j + 1:02d}").to_numpy(zero_copy_only=False)
                               for j in range(k)]).astype(np.float64)
        del knn
    finally:
        if not os.environ.get("TYPING_KEEP_SCRATCH"):
            shutil.rmtree(scratch, ignore_errors=True)
        else:
            print(f"{NAME}: scratch kept at {scratch}")

    if idx.shape[0] != n:
        raise RuntimeError(f"R returned {idx.shape[0]} rows for {n} cells")

    # scanpy's documented contract for a precomputed graph: n_neighbors - 1 non-zero entries per
    # row, self excluded. KODAMA already dropped self and re-sorted each row by its own distance.
    D = sp.csr_matrix((dst.ravel(), (np.repeat(np.arange(n), k), idx.ravel())), shape=(n, n))
    a = ad.AnnData(np.zeros((n, 1), dtype=np.float32))
    sc.pp.neighbors(a, n_neighbors=k + 1, distances=D)
    sc.tl.leiden(a, resolution=float(p["resolution"]), key_added="lc", random_state=seed,
                 flavor="igraph", n_iterations=2, directed=False)
    labels = a.obs["lc"].to_numpy().astype(np.int32)

    _, labels = np.unique(labels, return_inverse=True)   # 0..k-1 with no gaps, whatever igraph did
    labels = labels.astype(np.int32)
    print(f"[{NAME}] {n:,} cells -> {labels.max() + 1} clusters | KODAMA {meta['matrix_seconds']:.1f} s "
          f"(M={meta['M']}, landmarks={meta['landmarks']}, mean CV accuracy "
          f"{meta['mean_cv_accuracy']:.4f})")
    return labels
