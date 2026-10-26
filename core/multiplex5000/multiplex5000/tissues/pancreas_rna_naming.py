"""Pancreas Xenium RNA naming decisions used by the release worker and Stage 06.

These functions are moved verbatim from the current section worker. The worker
imports them here so there is one authority for the naming and admission rules.
"""

from __future__ import annotations

import numpy as np

MIN_NAMING_STANDOUT = 1
SHARED_SECRETORY_NAMING_GENES = frozenset({"XBP1", "FKBP11", "DERL3"})
NAMING_POLICIES = ("native", "shared_er_exact_v1")
REDIRECT_REASON = "cd3_negative_nk_positive_redirect"


def naming_policy_panels(panels: dict, policy: str) -> dict:
    """Keep shared secretory genes measured as support, never as B identity.

    The native policy is retained for provenance/reproduction. The selectable
    repair changes only B naming-score/admission roles; embedding features and
    full-core support diagnostics are unchanged.
    """
    if policy not in NAMING_POLICIES:
        raise ValueError(f"unknown naming policy: {policy}")
    if policy == "native":
        return panels
    out = dict(panels)
    out["cores"] = dict(panels["cores"])
    genes, naming = panels["cores"]["Immune_B_Plasma"]
    out["cores"]["Immune_B_Plasma"] = (genes, [g for g in naming if g not in SHARED_SECRETORY_NAMING_GENES])
    return out

def exact_detection_q(dets: np.ndarray, sizes: np.ndarray) -> np.ndarray:
    """Section-wide one-sided cluster-versus-rest detection test, BH corrected.

    Every group × every measured naming-basis feature belongs to the declared
    family. Clusters were selected from expression, so this is an exploratory
    contrast after selection, not a guaranteed lineage false-discovery rate.
    Statistical significance remains separate from program coherence.
    """
    from scipy.stats import hypergeom
    sizes = np.asarray(sizes, dtype=np.int64)
    if sizes.shape != (dets.shape[0],) or (sizes <= 0).any():
        raise ValueError("Exact detection contrast requires positive exact group sizes")
    if not np.isfinite(dets).all() or (dets < 0).any() or (dets > 1).any():
        raise ValueError("Detection fractions must be finite and lie between zero and one")
    detected = np.rint(dets * sizes[:, None]).astype(np.int64)
    # Fresh native statistics use float32. Valid ratios in large parents can
    # differ from their exact integer by >.01 purely through float32 rounding.
    # The tolerance follows the input precision and stays strictly below half a
    # count, so nearest-integer recovery remains unambiguous.
    precision = np.finfo(dets.dtype).eps
    tolerance = np.minimum(.49, np.maximum(.01, precision * sizes[:, None]))
    if (np.abs(detected - dets * sizes[:, None]) > tolerance).any():
        raise ValueError("Stored fractions do not recover exact integer detection counts")
    total = detected.sum(axis=0)
    p = hypergeom.sf(detected - 1, sizes.sum(), total[None, :], sizes[:, None])
    p[sizes.sum() - sizes == 0] = 1
    flat = p.ravel()
    order = np.argsort(flat, kind="stable")
    ranked = flat[order]
    values = np.minimum.accumulate((ranked * len(flat) / np.arange(1, len(flat) + 1))[::-1])[::-1]
    q = np.empty_like(values)
    q[order] = np.minimum(values, 1)
    return q.reshape(dets.shape)

def decide(z, dets, d_other, thr, core_names, score_cols, core_naming_cols, t_cols, nk_hint_cols,
           *, exact_detection_fdr=None, cluster_sizes=None):
    """Name every cluster of one partition by canonical enrichment against its siblings.

    A gene STANDS OUT in cluster k iff its cluster-mean lognorm is >= `thr` SD above the mean
    over that partition's own clusters AND its detection fraction in k exceeds its detection
    fraction in the rest of the partition.

    A core's SCORE is the mean z over `score_cols[core]` -- the genes of that core that may name
    its lineage.  Scoring the whole core instead biases the argmax by core size; the measurement
    behind that is in the module docstring.

    A cluster takes the argmax core's label only when all three hold:
      * the winning score is above the partition mean at all (z > 0) -- otherwise nothing is
        enriched here and the argmax is picking the least-negative of eleven negatives;
      * at least one gene that MAY NAME the winning lineage stands out (core/panels' own
        `evidence` column, values canonical or both);
      * that naming gene is detected above the cluster's OWN background: `bg` is the MEDIAN
        of the cluster's detection fractions over the measured set -- the 485 identity-core
        genes carried in this section, not the whole 5K panel -- and the gene must clear it
        strictly, `dets > bg`.

    WHERE THIS GUARD IS DEGENERATE.  `bg` is zero whenever at most half of those 485 genes are
    detected in any cell of the cluster (measured: every zero-background cluster in the release
    has 242 or fewer of the 485 detected, median 177).  The test then reduces to `dets > 0`, so
    the naming gene need be seen in ONE cell -- for those clusters the guard admits exactly what
    a background is supposed to exclude, and the claim it once carried here, that a near-empty
    cluster cannot be named by a gene seen in eight of its cells, was false for them.
    Recomputed 2026-09-18 from `_work/*.gene_stats.parquet` and `_work/*.subclusters.parquet`
    of `/data/xenium_typing/panc_2026-09-17/` (21 sections; the reconstruction reproduces
    `n_naming_above_background_of_winning_core` exactly for all 4,914 labelled subclusters):
      * 544 of 4,914 labelled subclusters (11.1 %) have `bg == 0`, covering 53,333 called cells
        (0.35 % of the 15,318,868 cells in labelled subclusters).  Per section the share of
        labelled subclusters with a zero background runs 3.5 % (6593, 10 of 289) to 29.9 %
        (6523, 61 of 204); 6539 is 22.1 % and 6380 20.3 %.
      * they are small: median 61 cells against 1,924 for the labelled subclusters whose
        background is above zero.
      * 72 of the 544 (10,354 cells) were admitted on a SINGLE naming gene above that zero
        background, and for 37 of those (4,706 cells, largest 563) the one gene was detected in
        exactly one cell -- e.g. 6523_c19_s07, 563 cells, Neural on TH in 1 cell; 6523_c06_s01,
        410 cells, Granulocyte on CXCR2 in 1 cell; HDL112_c10_s08, 360 cells, Granulocyte on
        PRG2 in 1 cell.
      * across all 544, 185 carry an admitted naming gene detected in exactly one cell, and 86
        (11,549 cells) have no admitted naming gene detected in more than two.
    A non-degenerate floor is a METHOD DECISION and is deliberately not taken in this function;
    the guard as written is a no-op above zero for those 544.
    A cluster that fails is uncalled, with the reason recorded. It is NEVER reassigned to the
    runner-up and never to a larger core.

    NEGATIVE EVIDENCE IS EVIDENCE: a cluster whose argmax is Immune_T but in which neither CD3E
    nor CD3G stands out, while NK genes do and an NK naming gene clears the bar, is NK.
    """
    n_groups = z.shape[0]
    so = (z >= thr) & (dets > d_other)
    bg = np.median(dets, axis=1)
    above_bg = dets > bg[:, None]
    native_above_bg = above_bg.copy()
    detection_q = None
    if exact_detection_fdr is not None:
        if not 0 < exact_detection_fdr < 1 or cluster_sizes is None:
            raise ValueError("Exact admission requires declared FDR and exact cluster sizes")
        detection_q = exact_detection_q(dets, cluster_sizes)
        above_bg &= (bg[:, None] > 0) | (detection_q <= exact_detection_fdr)
    scores = np.column_stack([
        z[:, score_cols[c]].mean(axis=1) if score_cols[c].size else np.full(n_groups, -np.inf)
        for c in core_names
    ])
    naming_so = np.column_stack([
        so[:, core_naming_cols[c]].sum(axis=1) if core_naming_cols[c].size
        else np.zeros(n_groups, dtype=int) for c in core_names
    ])
    naming_ok = np.column_stack([
        (so[:, core_naming_cols[c]] & above_bg[:, core_naming_cols[c]]).sum(axis=1)
        if core_naming_cols[c].size else np.zeros(n_groups, dtype=int) for c in core_names
    ])
    naming_before_exact = np.column_stack([
        (so[:, core_naming_cols[c]] & native_above_bg[:, core_naming_cols[c]]).sum(axis=1)
        if core_naming_cols[c].size else np.zeros(n_groups, dtype=int) for c in core_names
    ])
    win = scores.argmax(axis=1)
    srt = np.sort(scores, axis=1)
    margin = srt[:, -1] - srt[:, -2]
    nk_j = core_names.index("Immune_NK_ILC")
    calls, reasons, redirected = [], [], []
    for k in range(n_groups):
        w = core_names[win[k]]
        moved = False
        if w == "Immune_T":
            t_out = bool(so[k, t_cols].any()) if t_cols.size else False
            nk_out = bool(so[k, nk_hint_cols].any()) if nk_hint_cols.size else False
            if (not t_out) and nk_out and naming_before_exact[k, nk_j] >= MIN_NAMING_STANDOUT:
                w, moved = "Immune_NK_ILC", True
        j = core_names.index(w)
        if scores[k, j] <= 0:
            calls.append("Uncalled_no_enrichment")
            reasons.append(f"winning_core_{w}_score_{scores[k, j]:.3f}_not_above_partition_mean")
        elif naming_ok[k, j] >= MIN_NAMING_STANDOUT:
            calls.append(w)
            reasons.append(REDIRECT_REASON if moved else "")
        elif naming_so[k, j] >= MIN_NAMING_STANDOUT:
            calls.append("Uncalled_no_naming_evidence")
            suffix = ("at_or_below_cluster_own_background" if detection_q is None
                      else "at_or_below_cluster_own_background_or_failed_zero_background_exact_detection")
            reasons.append(f"winning_core_{w}_naming_standouts_{suffix}")
        else:
            calls.append("Uncalled_no_naming_evidence")
            reasons.append(f"winning_core_{w}_has_zero_standout_naming_genes")
        redirected.append(moved)
    return {"standout": so, "scores": scores, "naming_so": naming_so, "naming_ok": naming_ok,
            "win": win, "margin": margin, "calls": calls, "reasons": reasons,
            "redirected": redirected, "exact_detection_q": detection_q}

def redirect_counts(dec: dict, sizes) -> dict:
    """The negative-evidence redirect, counted over ONE partition -- the one it is given.

    Two different events, and they are not interchangeable:

    * APPLIED -- the argmax was Immune_T, CD3E/CD3G did not stand out, NK genes did, so the
      winner moved.  ``dec["redirected"]``.
    * SHIPPED -- that moved winner then cleared the enrichment gate and the admission rule, so
      the cluster went out labelled Immune_NK_ILC with :data:`REDIRECT_REASON` as its reason.
      A moved winner whose NK score is not above the partition mean is uncalled instead, and
      counting only ``redirected`` would report it as a label that never existed.

    ``sizes`` is that partition's own cell count per group, so the cell number belongs to the
    same partition as the cluster numbers beside it.  Pass the parent partition and you get the
    parent's counts; pass the final subcluster partition and you get the counts of the labels
    that ship.  The caller must say which in the field name.
    """
    shipped = [j for j, why in enumerate(dec["reasons"]) if why == REDIRECT_REASON]
    sizes = np.asarray(sizes)
    return {
        "clusters_applied": int(sum(bool(x) for x in dec["redirected"])),
        "clusters_shipped": len(shipped),
        "cells_shipped": int(sizes[shipped].sum()) if shipped else 0,
    }
