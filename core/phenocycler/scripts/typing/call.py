#!/usr/bin/env python3
"""Name every cell, from Nimbus per-marker positivity.

Two steps kept apart on purpose: **find structure, then name it.**  Cells are clustered on their
marker confidences with no reference to the lineage table, so the groups are whatever the data
contains.  Only then is each cluster named, by asking which lineage's markers stand out in it --
a comparison against a null, never a threshold.

Why confidences and not intensities.  Four intensity representations were tried on these two
donors and each failed for a reason worth recording, because they are the reasons this problem
has been hard:

* averaging a lineage's markers assumes a cell expresses all of them.  It expresses a subset -- a
  macrophage is CD68+/Iba1+ and negative for the other six myeloid channels -- so the score fell
  with the number of markers a lineage declares and only one- and two-marker lineages could win;
* standardising a marker across clusters deletes whatever is abundant: most of a pancreas section
  is exocrine, so most clusters are exocrine, so the real exocrine clusters score zero;
* within-donor ranking hides abundant cell types by construction.  A typical acinar cell sits at
  rank 0.5 on Pan-Cytokeratin *because* half the cells are acinar, while a T cell at 1% abundance
  sits at 0.99 on CD3e.  Composition came back with T cells inflated and acinar tissue missing;
* and underneath all of it, these distributions are unimodal in log space.  There is no threshold
  in the intensity to find, for most of the panel.

Nimbus confidences have none of those defects: they are per-marker probabilities on one [0, 1]
scale, absolute rather than relative to the population, and the model has already accounted for
brightness and for a neighbour's bleed.  They are used as features and never thresholded, so no
operating point is chosen anywhere in this file.
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
sys.path.insert(0, str(Path(__file__).parent))
import compartments as C  # noqa: E402
import lineages as L  # noqa: E402
import metacluster as MC  # noqa: E402

REPO = Path(__file__).resolve().parents[4]
DEV = REPO / "data/development/typing"

#: The most markers of one lineage that are asked to agree.  Evidence is read at every k from 1
#: to this, and the lineage is scored at whichever k speaks loudest -- see :func:`name_clusters`.
TOP_N = 3


def _stat_k(Z: np.ndarray, cols: list[int], k: int) -> np.ndarray:
    """Mean of a lineage's ``k`` strongest markers, per row of ``Z``."""
    sub = Z[:, cols]
    part = -np.sort(-np.nan_to_num(sub, nan=-np.inf), axis=1)[:, :k]
    return np.nanmean(np.where(np.isfinite(part), part, np.nan), axis=1)


def cluster(Z: np.ndarray, method: str, resolution: float, k: int, seed: int,
            neighbours: int) -> np.ndarray:
    """Group cells by their marker profile, with no reference to the lineage table.

    ``leiden`` (default) builds a k-nearest-neighbour graph and finds communities in it.  It is
    the right tool here and k-means is not: k-means minimises Euclidean distance to k centroids,
    which assumes clusters are round and of similar size.  Cell-type clusters are neither -- an
    acinar cluster and a delta-cell cluster differ in size by three orders of magnitude -- so
    k-means tends to split the abundant type and swallow the rare one, and it forces k to be
    guessed in advance.  Leiden reads local density instead, so a rare population survives as its
    own community and the resolution parameter sets granularity rather than a hard count.

    ``kmeans`` is kept so the two can be compared on evidence neither can tune: purity against the
    islet geometry, and how many cells land in Ambiguous because their cluster is a mixture.
    """
    if method == "kmeans":
        from sklearn.cluster import MiniBatchKMeans
        return MiniBatchKMeans(n_clusters=k, random_state=seed, n_init=10,
                               batch_size=10_000, max_iter=300).fit_predict(Z).astype(np.int32)
    if method != "leiden":
        raise ValueError(f"unknown clustering method {method!r}")
    import anndata as ad
    import scanpy as sc
    a = ad.AnnData(np.ascontiguousarray(Z, dtype=np.float32))
    # the confidences are already one comparable scale, so no PCA and no further scaling
    sc.pp.neighbors(a, n_neighbors=neighbours, use_rep="X", random_state=seed)
    sc.tl.leiden(a, resolution=resolution, key_added="lc", random_state=seed,
                 flavor="igraph", n_iterations=2, directed=False)
    return a.obs["lc"].to_numpy().astype(np.int32)


def load(donor: str) -> tuple[pd.DataFrame, list[str]]:
    pos = pd.read_parquet(DEV / f"positivity/donor_id={donor}/data.parquet")
    sub = pd.read_parquet(DEV / f"substrate/donor_id={donor}/data.parquet",
                          columns=["object_id", "eligible", "unavailable_reason",
                                   "has_nucleus", "has_cytoplasm",
                                   "cell_region", "islet_num", "X_centroid", "Y_centroid"])
    markers = [c for c in pos.columns if c in set(L.identity_markers())]
    df = sub.merge(pos[["object_id", "area", *markers]], on="object_id", how="left")
    return df, markers


def name_clusters(Z: np.ndarray, Xraw: np.ndarray, labels: np.ndarray, markers: list[str], k: int,
                  draws: int, margin: float, alpha: float, seed: int):
    """Cluster -> lineage.  Every lineage is read at its best k, against a null matched to both.

    WHY BEST-k AND NOT A FIXED ONE.  A cell type expresses a subset of what its lineage can
    express, and the size of that subset differs by cell type, not just by lineage.  Fixing the
    statistic at the mean of a lineage's top three markers penalises exactly the cell types whose
    identity rests on one: measured on donor 6442, a cluster with Podoplanin +4.87 -- lymphatic
    endothelium, which is genuinely Podoplanin-high and CD31-low -- was named *Neural*, because
    averaging Podoplanin with Vascular's two quiet markers dropped it below B3TUBB's lone +4.17.
    The same defect put 99,410 cells into Neural on a single marker at +1.90.  So the evidence is
    computed at k = 1, 2, 3 and the lineage is scored at whichever k is most surprising.

    Reading at several k means more chances to be surprised, so each (lineage, k) is standardised
    against its OWN null before they are compared, and the family-wise bar below is computed on
    the same best-of-k statistic under the null -- so the extra looks are paid for, not ignored.

    THE NULL randomises which CELLS form a cluster, not which markers form a lineage.  Marker
    permutation fails three ways here, each of which only shows up on a positive control: a random
    four-of-thirty-four set captures three of a lineage's four markers ~0.26% of the time, which
    floors the p-value; excluding a lineage's own markers removes the floor but depresses the null
    whenever that lineage sits above average (9 of 30 pure-noise clusters named); and a marker-set
    null is inflated by the very signal it is meant to detect.  Under random assignment a cluster
    of n cells has a mean vector about the global mean with covariance Sigma (N - n) / (n (N - 1)),
    which is sampled directly here, keeping the marker-to-marker correlation.

    MULTIPLICITY is a max-statistic (Westfall-Young) bar over the whole n_clusters x n_lineages
    family, so the ~13,500 tests that made a per-test p95 meaningless at 1500 clusters become ~280
    at a reviewable count, and the bar is one number rather than a per-test threshold.
    """
    M = np.full((k, len(markers)), np.nan)
    for c in range(k):
        sel = labels == c
        if sel.any():
            M[c] = np.nanmean(Z[sel], axis=0)

    idx = {m: i for i, m in enumerate(markers)}
    names, cols, carry_cols = [], [], []
    for lin, ms in L.LINEAGE_MARKERS.items():
        use = [idx[m] for m in ms if m in idx]
        if use:
            names.append(lin); cols.append(use)
            carry_cols.append([idx[m] for m in ms
                               if m in idx and m not in L.SUPPORTING_ONLY])

    rng = np.random.default_rng(seed)
    mu = np.nanmean(Z, axis=0)
    Sigma = np.cov(np.nan_to_num(Z - mu, nan=0.0), rowvar=False)
    chol = np.linalg.cholesky(Sigma + 1e-9 * np.eye(len(markers)))
    N = len(Z)
    sizes = np.array([max(int((labels == c).sum()), 1) for c in range(k)])

    S = np.full((k, len(cols)), np.nan)          # best-of-k evidence, in null SDs
    best_k = np.zeros((k, len(cols)), dtype=int)
    # How far each lineage's best marker rises above the donor's own baseline for that marker.
    base = np.nanmean(Xraw, axis=0)
    positive = Xraw >= MIN_CONF
    Rlift = np.full((k, len(cols)), np.nan)
    Rpos = np.full((k, len(cols)), np.nan)
    for c in range(k):
        sel = labels == c
        if sel.any():
            m_c = np.nanmean(Xraw[sel], axis=0)
            f_c = np.nanmean(positive[sel], axis=0)
            for j, col in enumerate(cols):
                # the same marker must satisfy both, so pick the best marker by lift and read
                # its positive fraction rather than mixing two different markers' evidence.
                # The marker that ADMITS the label must be one the lineage can be named on --
                # see lineages.SUPPORTING_ONLY.  A supporting marker still raises the z
                # statistic; it just cannot be the one that carries the call.
                carry = carry_cols[j]
                if not carry:
                    continue
                bestm = carry[int(np.nanargmax(m_c[carry] - base[carry]))]
                Rlift[c, j] = m_c[bestm] - base[bestm]
                Rpos[c, j] = f_c[bestm]
    fam_max = np.full(draws, -np.inf)
    for c in range(k):
        n_c = sizes[c]
        scale = np.sqrt(max((N - n_c) / (n_c * max(N - 1, 1)), 1e-12))
        sim = mu + (rng.standard_normal((draws, len(markers))) @ chol.T) * scale
        z_draw = np.full((draws, len(cols)), -np.inf)
        for j, col in enumerate(cols):
            for kk in range(1, min(TOP_N, len(col)) + 1):
                null = _stat_k(sim, col, kk)
                sd = null.std()
                if not np.isfinite(sd) or sd <= 0:
                    continue
                z_obs = (_stat_k(M[c][None, :], col, kk)[0] - null.mean()) / sd
                if not np.isfinite(S[c, j]) or z_obs > S[c, j]:
                    S[c, j], best_k[c, j] = z_obs, kk
                np.maximum(z_draw[:, j], (null - null.mean()) / sd, out=z_draw[:, j])
        np.maximum(fam_max, z_draw.max(axis=1), out=fam_max)
    bar = float(np.nanpercentile(fam_max, 100 * (1 - alpha)))
    # Both conditions: unusual for this donor AND actually raised above that donor's baseline.
    # z alone named 139,040 hormone-negative cells Endocrine on the donor that lost its beta cells;
    # an absolute floor alone handed every cluster to whichever lineage carried TP63.
    clears = (S > bar) & (Rlift >= MIN_RAW_LIFT) & (Rpos >= MIN_POSITIVE_FRACTION)

    named, why, best, runner = [], [], np.full(k, np.nan), np.full(k, np.nan)
    for c in range(k):
        cand = np.where(clears[c] & np.isfinite(S[c]))[0]
        if cand.size == 0:
            named.append(L.AMBIGUOUS)
            why.append("below_bar_or_not_positive")
            continue
        cand = cand[np.argsort(-S[c, cand])]
        best[c] = S[c, cand[0]]
        runner[c] = S[c, cand[1]] if cand.size > 1 else -np.inf
        if cand.size > 1 and (best[c] - runner[c]) < margin * abs(best[c]):
            named.append(L.AMBIGUOUS); why.append(f"tied:{names[cand[0]]}|{names[cand[1]]}")
        else:
            named.append(names[cand[0]]); why.append("")
    return named, why, S, np.full_like(S, bar), best, runner, names, M


#: The model's own decision boundary, used in the two places where the question is "is this
#: actually positive" rather than "is this unusual".  It is the one boundary in this lane and it
#: is Nimbus's, not one fitted here.
#:
#: WHY A FLOOR IS NECESSARY AT ALL.  Every other statistic in this file is a within-donor
#: comparison, and a within-donor comparison has no floor: standardising divides by the marker's
#: spread across that slide, so on a donor that has LOST a cell type the marker is absent
#: everywhere, its spread collapses, and a trivial absolute elevation becomes a large z.  Measured
#: on donor 6436 (T1D, beta cells destroyed): 139,040 cells were named Endocrine while sitting
#: outside any islet, with median GCG confidence 0.026 and INS 0.338 -- five standard deviations
#: above that donor's mean and biologically absent at the same time.
MIN_CONF = 0.5
SUBTYPE_MIN_CONF = MIN_CONF

#: How far a lineage's best marker must rise ABOVE THE DONOR'S OWN BASELINE, in raw confidence.
#:
#: A floor on the absolute level fails in both directions and was tried first: it hands a free
#: pass to any lineage carrying a broadly-positive marker -- Nimbus scores TP63 positive in 98.6%
#: of all cells, so Exocrine cleared an absolute 0.5 in every cluster on the slide and won a
#: cluster whose INS was +5.98 and whose Endocrine score was seventeen times higher -- while
#: blocking genuinely specific markers whose cluster MEAN sits just under the line (INS and GCG
#: exceed 0.5 in only 55-59% of true endocrine cells).
#:
#: Elevation separates all three cases in one criterion, in the model's own units: the false T1D
#: endocrine cluster lifts GCG 0.026 against a donor mean near 0.02 (~0.006); a real endocrine
#: cluster lifts INS ~0.55 against the same baseline (~0.53); TP63 lifts ~0.90 against a baseline
#: of ~0.90 (~0.00).  Two orders of magnitude between the real case and either failure.
MIN_RAW_LIFT = 0.10

#: And what fraction of the cluster's cells must actually be positive for that marker.
#:
#: Lift alone is not enough either, and the two failures are complementary rather than
#: alternatives.  On donor 6436 the false endocrine clusters lift INS to 0.338 -- genuinely raised
#: above a T1D donor's near-zero baseline, so lift passes -- but 0.338 is below the model's own
#: boundary, so the cells are not insulin-positive, and only 28-33% of them clear it.  A real
#: endocrine cluster on 6442 clears it in 55-59%.  Measured across the three cases this has to
#: separate:
#:
#:     lift    positive-fraction
#:     0.53    0.55-0.59   real endocrine (6442)      -> named
#:     0.32    0.28-0.33   false endocrine (6436)     -> refused, too few cells positive
#:     0.00    0.986       TP63, always-on            -> refused, never raised
#:
#: Lift rules out a permanently high marker; the positive fraction rules out a weakly raised one.
MIN_POSITIVE_FRACTION = 0.5

#: How discriminating a marker must be before a POSITIVE call on it can name a single cell.
#:
#: A flat boundary treats every marker alike, and on this panel that is useless for per-cell work:
#: measured on donor 6442, TP63 is positive in 95.5% of all cells, FOXP3 in 89.2%, CD79a in 82.1%,
#: so 98.8% of the unnamed cells are "positive" for Exocrine and 96.0% for T_NK and almost nothing
#: is positive for exactly one lineage.  Weighting by specificity -- one minus the marker's own
#: donor-wide positive fraction -- restores the distinction the reviewer makes by eye, because the
#: panel separates cleanly on this donor with a visible gap:
#:
#:     <=10% positive   INS, GCG, CD8, CD3e, CD68, SMA, CD31, CD20, MPO, Podoplanin, CD206,
#:                      CD209, CD11c, CD34, Caveolin, SST, CD56, Iba1, CD66, B3TUBB
#:     15-35% positive  Vimentin, Pan_Cytokeratin, Keratin_5, Collagen_IV, CD11b, M2Gal3, CD163
#:     >50% positive    CD4, E_cadherin, Ker8_18, EpCAM, CD79a, FOXP3, TP63  <- settle nothing
#:
#: The gap between CD163 (0.65) and CD4 (0.48) is where this sits.  Specificity is recomputed per
#: donor, never assumed, and the markers falling either side are recorded with the run.
MIN_SPECIFICITY = 0.60


def resolve_unnamed(df: pd.DataFrame, Xraw: np.ndarray, markers: list[str],
                    nuclear_only: np.ndarray, margin: float) -> tuple[np.ndarray, np.ndarray]:
    """Give a cell its own say when its cluster could not be named.

    Naming is done per cluster, so a cluster that two lineages claim equally leaves every one of
    its cells unnamed -- including cells that are individually unambiguous.  On donor 6442 that is
    199,951 cells, 99,410 of them in a single Myeloid-versus-Stromal tie, and a reviewer looking at
    the image can see that many are plainly positive for one myeloid marker or for Vimentin and
    negative for everything else that could compete.

    This resolves those cells one at a time, on their own markers, weighted by how discriminating
    each marker actually is on this donor (see :data:`MIN_SPECIFICITY`).  A marker positive in 96%
    of cells cannot name anything and is ignored here; a marker positive in 3% names a cell on its
    own.  Cells whose only positive markers are the always-on ones stay unnamed, and say so.

    Cells with a nucleus but no cytoplasm are read on their nuclear markers alone, because a
    surface or cytoplasmic marker measured on a mask that is only nucleus is measured in a
    compartment that cell does not have.
    """
    idx = {m: i for i, m in enumerate(markers)}
    spec = 1.0 - np.nanmean(Xraw >= MIN_CONF, axis=0)          # per marker, this donor
    nuclear = {m for m in C.NUCLEUS}
    lin_names, lin_cols = [], []
    for lin, ms in L.LINEAGE_MARKERS.items():
        use = [idx[m] for m in ms if m in idx and spec[idx[m]] >= MIN_SPECIFICITY]
        if use:
            lin_names.append(lin); lin_cols.append(use)
    if not lin_cols:
        return df["lineage"].to_numpy().astype(object), df["reason"].to_numpy().astype(object)

    positive = Xraw >= MIN_CONF
    # a marker a cell cannot legitimately be scored on contributes nothing
    allowed = np.ones(len(markers), bool)
    nuc_mask = np.array([markers[i] in nuclear for i in range(len(markers))])
    ev = np.where(positive, spec[None, :], 0.0)
    ev_nuc = np.where(positive & nuc_mask[None, :], spec[None, :], 0.0)

    lineage = df["lineage"].to_numpy().astype(object).copy()
    reason = df["reason"].to_numpy().astype(object).copy()
    unnamed = (lineage == L.AMBIGUOUS)
    if not unnamed.any():
        return lineage, reason
    rows = np.where(unnamed)[0]
    src = np.where(nuclear_only[rows][:, None], ev_nuc[rows], ev[rows])
    scores = np.vstack([src[:, c].max(axis=1) for c in lin_cols]).T
    order = np.argsort(-scores, axis=1)
    top = scores[np.arange(len(scores)), order[:, 0]]
    two = scores[np.arange(len(scores)), order[:, 1]] if scores.shape[1] > 1 else np.zeros(len(scores))
    ok = (top >= MIN_SPECIFICITY) & ((top - two) >= margin * top)
    named = np.array(lin_names, dtype=object)[order[:, 0]]
    lineage[rows[ok]] = named[ok]
    reason[rows[ok]] = "resolved_per_cell"
    reason[rows[~ok & nuclear_only[rows]]] = "nucleus_only_no_discriminating_marker"
    return lineage, reason


def subtype(df: pd.DataFrame, X: np.ndarray, markers: list[str], sep: float) -> np.ndarray:
    """Resolve subtypes within a parent lineage, per cell, on the RAW confidences.

    ON RAW AND NOT STANDARDISED VALUES.  This is the opposite of what the clustering needs and the
    reason is worth stating, because using the standardised matrix here silently destroyed a whole
    cell type.  Standardising divides each marker by its spread across the entire section: INS and
    GCG are near zero in ~98.5% of a pancreas, so their spread is tiny and an islet cell scores
    z = +8 to +13 on them, while SST is more broadly expressed, its spread is larger, and its z
    ceiling is +4.4.  Argmax over z therefore ranks markers by how RARE they are rather than by
    how positive the cell is: measured on donor 6442, 4,357 of the 4,414 SST-high endocrine cells
    had INS or GCG "higher" in z, and Delta came back as 41 cells in 1.7 million.  The Nimbus
    confidences are already one common scale -- each is the probability that this cell is positive
    for that marker -- which is exactly the comparison a subtype call needs.  On raw values the
    same donor gives Delta 1,858 cells, 3.8% of the endocrine cells inside islet cores.

    A cell whose best two subtypes are closer than ``sep`` keeps the parent lineage rather than
    being forced to choose, and so does one whose winner is below the model's own boundary.  Both
    are honest outcomes: an islet cell that is genuinely beta-and-alpha-ambiguous is Endocrine.
    """
    idx = {m: i for i, m in enumerate(markers)}
    out = df["lineage"].to_numpy().astype(object).copy()
    for parent, subs in L.SUBTYPE_MARKERS.items():
        sel = (df["lineage"] == parent).to_numpy()
        if not sel.any():
            continue
        usable = {s: [idx[m] for m in ms if m in idx] for s, ms in subs.items()}
        usable = {s: c for s, c in usable.items() if c}
        if len(usable) < 2:
            continue
        names = list(usable)
        sc = np.vstack([np.nanmax(X[np.ix_(sel, usable[s])], axis=1) for s in names]).T
        order = np.argsort(-np.nan_to_num(sc, nan=-np.inf), axis=1)
        top = sc[np.arange(len(sc)), order[:, 0]]
        two = sc[np.arange(len(sc)), order[:, 1]]
        ok = (np.isfinite(top) & (top >= SUBTYPE_MIN_CONF)
              & (top - np.where(np.isfinite(two), two, -np.inf) >= sep))
        lab = np.array(names, dtype=object)[order[:, 0]]
        out[np.where(sel)[0][ok]] = lab[ok]
    return out


def modifiers(df: pd.DataFrame, X: np.ndarray, markers: list[str]) -> np.ndarray:
    """Refinements OF a subtype, scored inside the cells that already carry it.

A refinement is not an alternative: a FOXP3-high CD4 T cell is still a CD4 T cell.  Scoring
    such a thing as a peer in the same argmax gave 11,463 "Tregs" on donor 6442 -- a third of
    every T cell -- because Nimbus scores FOXP3 above its own boundary in 90.4% of cells, so it
    won the argmax nearly everywhere.

    The flag marks the upper quartile of the parent population on the modifier's markers, so its
    SIZE IS 25% OF THE PARENT BY CONSTRUCTION and means nothing on its own; what it carries is
    which cells rank highest, against their own kind rather than against the whole section.  It is
    named for the measurement (``FOXP3_high``) and not for the cell type someone might infer from
    it, because this panel cannot support a Treg call and a percentile wearing a cell-type name
    would read as one.
    """
    idx = {m: i for i, m in enumerate(markers)}
    out = np.full(len(df), "", dtype=object)
    ct = df["cell_type"].to_numpy()
    for parent, (name, ms) in L.SUBTYPE_MODIFIERS.items():
        cols = [idx[m] for m in ms if m in idx]
        sel = ct == parent
        if not cols or sel.sum() < 100:
            continue
        v = np.nanmax(X[np.ix_(np.where(sel)[0], cols)], axis=1)
        finite = np.isfinite(v)
        if finite.sum() < 100:
            continue
        # upper quartile of the parent population, against that population and nothing else
        cut = np.nanpercentile(v[finite], 75)
        hit = np.where(sel)[0][finite & (v > cut)]
        out[hit] = name
    return out


def run(donor: str, method: str, resolution: float, neighbours: int,
        k: int, seed: int, draws: int, margin: float, alpha: float, sep: float,
        merge_lo: int, merge_hi: int) -> Path:
    df, markers = load(donor)
    X = df[markers].to_numpy(np.float64)
    have = df["eligible"].to_numpy(bool) & np.isfinite(X).all(axis=1)

    # Standardise each marker across cells so channels that are broadly positive -- Nimbus scores
    # TP63 above 0.5 in 98.6% of cells, FOXP3 in 90.4% -- contribute their variation and not
    # their level.  No marker is dropped for this; the panel decides that, not the statistic.
    mu, sd = np.nanmean(X[have], axis=0), np.nanstd(X[have], axis=0)
    Z = (X - mu) / np.where(sd > 0, sd, np.nan)

    fine = np.full(len(df), -1, np.int32)
    fine[have] = cluster(Z[have], method, resolution, k, seed, neighbours)
    n_fine = int(fine[have].max()) + 1

    # TWO STAGES.  Cluster fine so a rare population survives -- the delta cells are 44 of 11,659
    # in the validated window, and at a few dozen clusters a population that size is absorbed into
    # a larger one -- then merge to a count a person can read, because the reviewer looking at
    # these calls on the images is the only real check on whether they are right.  The merge runs
    # on the profile matrix, one row per cluster, so a 44-cell cluster weighs as much as a
    # 40,000-cell one; that equal weighting is what carries the rare type through instead of
    # dissolving it.  Naming then happens at the merged granularity, where the test family is a
    # few hundred rather than the ~13,500 that made the bar meaningless.
    sweep = None
    if merge_hi and n_fine > merge_lo:
        P = MC.profiles(Z[have], fine[have], n_fine)
        sweep = MC.sweep(P, merge_lo, min(merge_hi, n_fine - 1), draws=15, seed=seed)
        n_final = MC.choose(sweep, merge_lo, min(merge_hi, n_fine - 1))
        mapping = MC.merge(P, n_final)
        labels = np.where(fine >= 0, mapping[np.clip(fine, 0, None)], -1).astype(np.int32)
        n_clusters = n_final
        print(f"[{donor}] {n_fine} fine clusters -> {n_final} merged "
              f"(stability {float(sweep.set_index('n_final').loc[n_final, 'stability_ari']):.3f})",
              flush=True)
    else:
        # A method whose granularity is set by its own parameter (LDA's topic count, say) may
        # produce fewer clusters than the reviewable band.  Merging is then meaningless and the
        # method is simply doing both stages at once -- recorded, not worked around.
        labels, n_clusters = fine, n_fine
        print(f"[{donor}] {n_fine} clusters; no merge (at or below the reviewable band)",
              flush=True)

    named, why, S, qval, best, runner, lin_names, M = name_clusters(
        Z[have], X[have], labels[have], markers, n_clusters, draws, margin, alpha, seed)

    lineage = np.full(len(df), L.UNAVAILABLE, dtype=object)
    reason = df["unavailable_reason"].to_numpy().astype(object).copy()
    lineage[have] = np.array(named, dtype=object)[labels[have]]
    reason[have] = np.array(why, dtype=object)[labels[have]]
    missing = df["eligible"].to_numpy(bool) & ~have
    reason[missing] = "no_positivity_row"

    df["lineage"] = lineage
    df["reason"] = reason

    # Cells with a nucleus but no cytoplasm may only be read on nuclear markers.
    nuc_only = (df.get("has_cytoplasm", pd.Series(True, index=df.index)).fillna(True)
                .astype(bool).to_numpy() == False)
    X_typed = X.copy()
    if nuc_only.any():
        non_nuclear = [i for i, m in enumerate(markers) if m not in set(C.NUCLEUS)]
        X_typed[np.ix_(nuc_only, non_nuclear)] = np.nan
        print(f"[{donor}] {int(nuc_only.sum()):,} cells have a nucleus but no cytoplasm; "
              f"read on nuclear markers only", flush=True)

    df["lineage"], df["reason"] = resolve_unnamed(df, X_typed, markers, nuc_only, margin)
    n_res = int((df["reason"] == "resolved_per_cell").sum())
    print(f"[{donor}] {n_res:,} cells resolved individually after their cluster could not be named",
          flush=True)
    df["cluster"] = labels
    df["fine_cluster"] = fine
    df["cell_type"] = subtype(df, X_typed, markers, sep)
    df["modifier"] = modifiers(df, X, markers)

    dest = DEV / f"calls/donor_id={donor}"; dest.mkdir(parents=True, exist_ok=True)
    keep = ["object_id", "cluster", "fine_cluster", "lineage", "cell_type", "modifier",
            "reason", "eligible",
            "cell_region", "islet_num", "X_centroid", "Y_centroid"]
    df[keep].to_parquet(dest / "data.parquet", index=False)
    pd.DataFrame({"cluster": np.arange(n_clusters), "lineage": named, "reason": why,
                  "best": best, "runner_up": runner,
                  "n_cells": np.bincount(labels[have], minlength=n_clusters),
                  **{f"S__{n}": S[:, i] for i, n in enumerate(lin_names)},
                  **{f"bar__{n}": qval[:, i] for i, n in enumerate(lin_names)},
                  **{f"conf__{m}": M[:, i] for i, m in enumerate(markers)}}
                 ).to_csv(dest / "clusters.csv", index=False)
    (dest / "params.json").write_text(json.dumps(
        {"donor_id": donor, "method": method, "resolution": resolution,
         "n_neighbours": neighbours, "n_clusters": n_clusters,
         "k_requested": k if method == "kmeans" else None,
         "n_fine_clusters": n_fine, "merge_band": [merge_lo, merge_hi],
         "merge_sweep": (sweep.to_dict("records") if sweep is not None else None),
         "seed": seed, "null_draws": draws, "margin": margin, "fdr_alpha": alpha,
         "n_tests": int(n_clusters * len(lin_names)),
         "subtype_separation": sep, "markers": markers, "lineages": lin_names,
         "n_objects": int(len(df)), "n_typed": int(have.sum())}, indent=2))
    return dest / "data.parquet"


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--donors", nargs="+", required=True)
    ap.add_argument("--method", choices=("leiden", "kmeans"), default="leiden")
    ap.add_argument("--resolution", type=float, default=2.0, help="leiden granularity")
    ap.add_argument("--neighbours", type=int, default=15, help="leiden kNN graph degree")
    ap.add_argument("--k", type=int, default=1500, help="kmeans only")
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--null-draws", type=int, default=2_000,
                    help="draws for the max-statistic bar")
    ap.add_argument("--fdr", type=float, default=0.05,
                    help="family-wise error rate for the max-statistic bar")
    ap.add_argument("--margin", type=float, default=0.15,
                    help="top-to-runner-up separation required, as a fraction of the top score")
    ap.add_argument("--subtype-sep", type=float, default=0.05,
                    help="separation required between the best two subtypes, in CONFIDENCE units "
                         "(0-1), not standard deviations")
    ap.add_argument("--merge-lo", type=int, default=30,
                    help="smallest reviewable merged cluster count")
    ap.add_argument("--merge-hi", type=int, default=40,
                    help="largest; 0 disables the merge")
    a = ap.parse_args()
    for d in a.donors:
        p = run(d, a.method, a.resolution, a.neighbours, a.k, a.seed,
                a.null_draws, a.margin, a.fdr, a.subtype_sep, a.merge_lo, a.merge_hi)
        t = pd.read_parquet(p)
        print(f"\n[{d}] {p}")
        for col in ("lineage", "cell_type", "modifier"):
            print(f"  -- {col} --")
            for v, n in t[col].value_counts().items():
                print(f"   {v:<16} {n:>9,}  {n/len(t)*100:>5.1f}%")


if __name__ == "__main__":
    main()
