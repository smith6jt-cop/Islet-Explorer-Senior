#!/usr/bin/env python3
"""Group cells by Latent Dirichlet Allocation over their Nimbus marker confidences.

RETIRED -- NOT A FALLBACK.  Recorded here 2026-09-18 because the docstring below reads as a live
specification and is not one.  This module belongs to the broad-typing producer comparison
WITHDRAWN 2026-09-08 (``cluster_leiden.py`` / ``cluster_lda.py`` / ``cluster_kodama.py`` scored by
``compare_clustering.py``).  It did not produce the calls on disk, and the current typing is the
subpanel classification in ``subpanel.py``; see ``README.md``, "Retired, and not to be revived as
fallbacks".  What follows is kept as a record of what was measured and why.

One of the clustering methods offered to `call.py`; the caller names whatever comes back, so
nothing here consults the lineage table and no operating point is chosen on a confidence.

**What LDA assumes, and what that costs on this data.**  LDA is a mixture model over
non-negative counts: a cell is a *composition* over latent topics, a topic is a distribution over
markers, and the per-cell counts are drawn from that mixture.  Two consequences matter here and
both were measured rather than assumed.

*It normalises every cell.*  A cell whose markers all sit at 0.05 and a cell whose markers all sit
at 0.70 have the same composition and are the same document to LDA.  Nimbus confidences are
absolute per-marker probabilities -- that is the whole reason this lane uses them instead of
intensities -- so LDA discards the one property that made them worth having.  An encoding that
puts absence back in the vocabulary (each marker as two entries, ``conf`` and ``1 - conf``) was
tried for exactly this reason.  At k=20 it scores 0.3439 mean / 0.2923 median on held-out marker
agreement against 0.3485 / 0.3025 for the encoding below -- the same, for twice the vocabulary
and 1.5x the runtime -- and in its unscaled form it spends two of its live topics on pure
*absence* programs (~CD79a/~FOXP3/~TP63, correlation -0.62 with a cell's total confidence mass)
and leaves 5.5 of 20 topics holding no cell at all.  Not worth it.

*It is driven by whichever markers carry the mass.*  On the test window the five highest-mean
channels -- TP63 0.73, FOXP3 0.67, CD79a 0.62, E-cadherin 0.48, EpCAM 0.45 -- hold 42% of all
confidence mass, and TP63's standard deviation is 0.080, so on raw confidences TP63 is a stopword:
present in every document at the same rate.  Dividing each marker by its own mean fixes that (see
``_counts``): at a matched document length and k=20 it is worth 0.3274 against 0.3148 for raw
confidences, and at its own best document length 0.3485.

**Is LDA the right model for this?  Measured, and the honest answer is: it works, and it does not
beat the simpler alternative.**  Against a partition of the same size, on the same criterion
(fit on 33 markers, ask what fraction of the 34th marker's variance the partition explains --
chance 0.0017), LDA scores 0.2609 at k=8 where MiniBatchKMeans on the standardised matrix scores
0.2775, and 0.3080 at k=12 where k-means scores 0.3092.  LDA is at or just below k-means at every
matched k tried, and it costs 9.6 s against k-means' 0.10 s on the same 11,659 x 34 matrix, a
factor of 96.  The thing LDA is supposed to add in exchange -- an honest mixture description of a
doublet or an overlapped cell -- does not materialise here either: at the settings below the
per-cell posterior is nearly one-hot (median top topic 0.964, median gap to the second 0.959,
only 2.8% of cells within 0.05 of a tie), because 34 counts
against a document-topic prior of total mass 1 concentrate it.  On this data LDA is a working
clustering method, not a better one, and the mixture semantics that motivate it are inert.

**Do the topics correspond to cell types, or to marker programs?  Both, and the split is
measurable.**  Top markers of each topic by lift over the corpus marginal, at the shipped
setting: one topic is Endocrine (GCG, INS, SST) and is 98.5% islet-core against QuPath polygons
the model never saw; others carry Caveolin/CD34/CD31 (vascular), Iba1/CD206/CD68 (macrophage),
MPO/CD68/Iba1 (granulocyte-myeloid), CD66/MPO (granulocyte) and CD3e/CD8 (T, mixed with
SMA/Podoplanin/Vimentin -- at this many topics a pure T topic is not affordable).  The remaining
two topics hold 75% of the cells and are the exocrine bulk split by *which non-epithelial channel
is elevated alongside the keratins* -- CD11b/M2Gal3 in one, Pan-Cytokeratin/Ker8_18 in the other.
That is a co-expression program, not an exocrine subtype, and it is where three quarters of this
window ends up.  It is not a defect of LDA: the panel carries no acinar-specific and no
ductal-specific marker, so no method has an exocrine split to find, and the topics spend
themselves on programs instead.  The *general brightness* topic this encoding is supposed to risk
did not appear -- no topic correlates above 0.20 with a cell's total confidence mass.

**The hard label is the argmax topic, and nothing more.**  `call.py` owns the only ambiguity
mechanism in this lane (`Ambiguous`, decided per *cluster* against a null of random marker sets);
a per-cell "mixed" outcome emitted here would arrive at `name_clusters` as one more cluster to
name, and on the measurement above there would be almost nothing in it.

**The partition is only moderately reproducible, and the instability has an address.**  Mean
pairwise adjusted Rand over 8 seeds is 0.588 +- 0.143 across all 11,659 cells but 0.716 +- 0.090
across the 2,881 cells that fall outside the two exocrine-bulk topics: the seeds disagree about
where to cut the exocrine mass, which is the part of the section that has no split to find, and
agree much better everywhere else.

**The islet claim first made here was too strong and is corrected.**  It said every one of 8 seeds
returns exactly one cluster above 90% islet-core at purity 0.983 +- 0.007.  Independent
re-measurement over 16 seeds says otherwise: k = 8 in 16/16, but the islet mass is isolated above
90% purity in 13 of 16 seeds (11 as a single cluster, 2 split across two, union purity >= 0.98,
recall 0.82-0.85) and NOT isolated in 3 of 16, where measured purity is 0.729-0.745.  Over seeds
0-7 that is purity 0.892 +- 0.119, not 0.983 +- 0.007.  Recall (0.836 +- 0.005) is the part that
holds.

The failure is diagnosable, not noise, and it is this model's characteristic risk on this data:
in the three failing seeds the endocrine topic FUSES WITH A T-CELL PROGRAM.  Its top markers by
lift become GCG 8.4, INS 7.7, CD3e 4.0, CD8 3.8, SST 2.1 -- against GCG 11.1, INS 10.1, SST 2.3
and no immune marker at seed 0 -- and it absorbs ~180 non-islet cells.  A topic model finds marker
PROGRAMS, and a program is not obliged to be a cell type.

**Reproducing a run needs the thread count, not just the seed.**  ``n_jobs`` changes the partition
at a fixed seed: ARI against n_jobs=8 at seed 0 is 0.965 (1), 0.992 (2), 0.931 (4), 0.996 (16).
sklearn's ``_e_step`` slices documents with ``gen_even_slices(n, n_jobs)`` and each slice draws its
doc-topic initialisation from a copy of the same RandomState, so the per-document init depends on
how the rows were sliced.  Record n_jobs beside the seed, and note that raising it to shorten a
whole-slide run changes the clustering as well as the clock.

**Whole-slide cost, measured rather than extrapolated.**  ``fit_predict`` at these settings on
1,900,000 rows x 34 markers took **29.8 minutes** of wall clock at ``n_jobs=8``, peak RSS
2.1 GiB, and populated all 8 topics -- on a machine simultaneously running the Nimbus GPU job and
two other clustering methods, so it is an upper bound rather than a best case.  That is 1.14x the
straight-line extrapolation from the 9.6 s test-set fit, i.e. the cost is linear in cells with no
surprise term: online variational Bayes visits each mini-batch once per pass and the per-pass work
is O(cells x markers x topics).  Nothing needs subsampling at this size; the whole donor fits in
memory twice over (a 1.9M x 34 float64 matrix is 0.48 GiB, and ``_counts`` makes one copy of it).
Halving the wall clock is a matter of raising ``n_jobs`` or lowering ``max_iter`` to 20, which the
convergence numbers below say costs nothing.  If a future donor ever does exceed memory,
``LatentDirichletAllocation.partial_fit`` takes mini-batches from disk with the same update, at
the price of a second code path; it is not needed for 1.4-1.9M cells.
"""

from __future__ import annotations

import numpy as np

NAME = "lda"

PARAMS: dict[str, object] = {
    #: Topics.  Held-out perplexity on a 20% cell holdout, averaged over 8 seeds: 31.129 +- 0.070
    #: at 5, 31.084 +- 0.051 at 6, **31.057 +- 0.082 at 8**, 31.098 +- 0.211 at 10,
    #: 31.253 +- 0.131 at 12, 31.635 +- 0.175 at 16, and on a single seed it goes on rising --
    #: 31.9 at 20, 32.4 at 25, 33.7 at 40, 34.7 at 50.  8 is the minimum; 6 and 10 are inside one
    #: standard deviation of it and 12 is outside.  Nothing about the number of cell types on the
    #: panel entered this -- and 8 is fewer than the panel declares, which is the finding rather
    #: than a target to be met.
    "n_components": 8,
    #: Mean pseudo-count budget per cell.  34 is the value at which the marker matrix needs no
    #: rescaling at all after ``_counts`` divides by the column means (mean row sum = one count
    #: per marker), and it is where held-out marker agreement peaks: median eta^2 across the 34
    #: leave-one-marker-out fits is 0.2574 at 8.5, 0.2746 at 17, **0.3025 at 34**, 0.2632 at 51,
    #: 0.2494 at 70, 0.2567 at 140 and 0.2271 at 350.  Inflating the counts sharpens the
    #: posterior but sharpens it onto the abundant markers.
    "doc_length": 34.0,
    #: Online (mini-batch) variational Bayes rather than batch, so one code path serves 11.6k
    #: cells and 1.9M.  It costs nothing here: held-out perplexity averaged over 3 seeds is
    #: 31.111 for full batch against 30.837 / 31.003 / 31.146 for online at batch sizes
    #: 1024 / 4096 / 16384 -- the whole spread sits inside the 0.082 seed-to-seed standard
    #: deviation measured at this K.
    "learning_method": "online",
    "batch_size": 4096,
    #: Passes over the data.  Held-out perplexity is 31.138 after 5, 31.037 after 10, 30.995
    #: after 20, 30.989 after 30, 30.986 after 50 and 30.990 after 80 -- flat from 20 on, which
    #: is also where the partition stops moving (ARI against the 5-iteration fit settles at
    #: 0.86).  30 clears the plateau, and the cost is linear in it, which is what decides it at
    #: whole-slide size.
    "max_iter": 30,
    #: sklearn parallelises the E-step over cells.  Left explicit because the machine is shared
    #: with the Nimbus GPU job's loaders; raise it for a whole-slide run.
    "n_jobs": 8,
}


def _counts(X: np.ndarray, doc_length: float) -> np.ndarray:
    """Confidences -> the fractional counts LDA is defined over.

    Each marker is divided by its own mean across the cells being clustered, so every channel
    contributes the same total mass to the corpus and none of them is a stopword.  This is the
    non-negative counterpart of the standardisation `call.py` applies for the graph methods --
    same purpose, stated in the file as "contribute their variation and not their level" -- done
    as a scaling rather than a shift because a shift would put negative values into a model that
    has no meaning for them.

    The mean is taken over the matrix handed in, so the statistic never crosses a donor.
    """
    mu = X.mean(axis=0)
    dead = np.flatnonzero(mu <= 0)
    if dead.size:
        raise ValueError(
            f"columns {dead.tolist()} have zero mean confidence; a marker that is absent for "
            f"every cell must be dropped upstream, not divided by"
        )
    return (X / mu) * (doc_length / X.shape[1])


def _fit(X: np.ndarray, seed: int, p: dict):
    from sklearn.decomposition import LatentDirichletAllocation

    model = LatentDirichletAllocation(
        n_components=int(p["n_components"]),
        # doc_topic_prior and topic_word_prior are left at sklearn's 1/n_components.  With
        # doc_length=34 the observed counts outweigh the document-topic prior's total mass of 1
        # by 34:1, so the budget above is the knob that governs posterior sharpness and the
        # prior is not a second one to tune.
        learning_method=str(p["learning_method"]),
        batch_size=int(p["batch_size"]),
        max_iter=int(p["max_iter"]),
        n_jobs=int(p["n_jobs"]),
        evaluate_every=-1,
        random_state=seed,
    )
    return model, model.fit_transform(X)


def _theta(X, *, seed: int = 0, **params) -> np.ndarray:
    """The per-cell topic posterior, for diagnostics.  Not part of the clustering interface."""
    p = {**PARAMS, **params}
    _, theta = _fit(_counts(np.asarray(X, dtype=np.float64), float(p["doc_length"])), seed, p)
    return theta


def fit_predict(X, Z, *, seed=0, **params) -> np.ndarray:
    """X: raw confidences, cells x markers, non-negative, in [0,1], no NaN.
       Z: the same matrix standardised per marker (zero mean, unit variance per column).
       returns: int32 array, one cluster label per cell, labels 0..k-1, no gaps.

    Uses **X**.  LDA's likelihood is multinomial, so its features have to be non-negative counts;
    ``Z`` is centred and roughly half of it is negative, which the model has no meaning for.  The
    equalisation ``Z`` provides is done here instead by dividing each marker by its own mean,
    which keeps every value non-negative (see ``_counts``).

    Returns fewer than ``n_components`` labels when a topic wins no cell's argmax -- it is a
    mixture component, not a partition cell, and there is no guarantee it takes any cell.  Labels
    are compacted so the caller's ``range(k)`` loops stay correct.
    """
    p = {**PARAMS, **params}
    X = np.asarray(X, dtype=np.float64)
    if X.ndim != 2 or X.shape[0] == 0 or X.shape[1] == 0:
        raise ValueError(f"X must be a non-empty 2-D cells x markers matrix, got {X.shape}")
    if not np.isfinite(X).all():
        raise ValueError(
            "X carries a non-finite value; a missing measurement is unavailable, not zero, and "
            "must be excluded by the caller before it reaches the clustering"
        )
    if X.min() < 0.0 or X.max() > 1.0:
        raise ValueError(
            f"X must be Nimbus confidences in [0, 1], got [{X.min():.4g}, {X.max():.4g}]"
        )

    _, theta = _fit(_counts(X, float(p["doc_length"])), int(seed), p)
    # argmax, then renumber: the caller needs 0..k-1 with no gaps, and unpopulated topics leave
    # gaps.  No per-cell "mixed" outcome is invented here -- `call.py` owns Ambiguous.
    _, labels = np.unique(theta.argmax(axis=1), return_inverse=True)
    return labels.astype(np.int32)
