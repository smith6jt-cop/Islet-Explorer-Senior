"""Does the registration put the IF's islets where Xenium's endocrine RNA is?

The check is independent by construction, and that is the whole point of doing it this way.
The registration aligned **tissue outlines only** -- it never saw an islet and never saw a
gene.  The IF islets come from a pixel classifier on protein.  The Xenium score comes from the
18-gene ``01_pancreas_endocrine`` identity core.  So agreement between them is two instruments
finding the same structures, not one being told where the other's are.

It is also the smallest real cross-modality measurement this project can make: it needs no
cell typing on either side, no verified marker pair, and no hormone threshold.

Measured 2026-08-28, cells within 25 um of a kept IF islet centroid:

    donor 6442 (ND) : endocrine core 2.0340 inside vs 0.1176 elsewhere = 17.3x
    donor 6436 (T1D): endocrine core 0.3577 inside vs 0.1265 elsewhere =  2.8x
    random-placement null on both: ~0.117 (p95 ~0.16)

**Do not read 6436's weaker enrichment as beta-cell loss without more.**  Split by the core's
own ``subtype_hint``, 6436 gives alpha 2.83x and beta 2.89x -- moving together.  Beta-cell
loss depresses beta *relative to* alpha; a uniform scale-down of every endocrine gene set is
what a mask landing partly off-target does, and 6436 is also the donor with the worse
registration (Dice 0.890 vs 0.951) and a 4% tissue-extent mismatch.  6442 by contrast gives
beta 25.4x against alpha 6.3x, which is the pattern an intact islet should show.  n=1 per
group; this is a registration diagnostic, not a disease result, and nothing here may be
ordered by disease state.

``HHEX`` (delta) reads *below* background in both donors and is not interpretable: it is a
one-gene panel, which is exactly the thin-panel trap the Xenium lane's decision #12 exists to
prevent.
"""

from __future__ import annotations

import csv
from pathlib import Path

import numpy as np

#: Cells this far from a kept islet centroid count as "inside".  ~ one islet radius for the
#: smallest object kept (50 um diameter), so it is the neighbourhood the landmark actually
#: asserts rather than a tuned window.
NEAR_UM = 25.0

#: Draws of the placement null.  Each draw re-places the same number of points uniformly over
#: the section's bounding box; points off tissue simply match no cells, so the null is
#: effectively "the same number of arbitrary on-tissue neighbourhoods".
NULL_DRAWS = 20


def endocrine_core(subpanel_dir: Path) -> tuple[list[str], dict[str, list[str]]]:
    """The 18-gene pancreas endocrine identity core, and its ``subtype_hint`` split."""
    root = Path(subpanel_dir)
    candidates = (
        root / "tissues" / "pancreas" / "identity_core" / "01_pancreas_endocrine.csv",
        root / "identity_core" / "pancreas" / "01_pancreas_endocrine.csv",
        root / "pancreas" / "identity_core" / "01_pancreas_endocrine.csv",
    )
    p = next((candidate for candidate in candidates if candidate.is_file()), candidates[0])
    rows = list(csv.DictReader(p.read_text().splitlines()))
    by_hint: dict[str, list[str]] = {}
    for r in rows:
        by_hint.setdefault(r.get("subtype_hint") or "pan", []).append(r["gene"])
    return [r["gene"] for r in rows], by_hint


def gene_set_score(counts, genes: list[str], gene_index: dict[str, int]) -> np.ndarray:
    """Per-cell mean over a gene set, each gene scaled to unit SD.

    The SD scaling is the same correction the Xenium lane's ``nullgate.standardise_genes``
    applies, and for the same measured reason: without it one bright gene carries the set, and
    the score stops being about the set.  Only the scale is applied, never the centring, so a
    sparse matrix stays sparse.
    """
    import scipy.sparse as sp
    cols = [gene_index[g] for g in genes if g in gene_index]
    if not cols:
        raise ValueError("no gene of this set is on the panel")
    X = sp.csr_matrix(counts)
    mean = np.asarray(X.mean(axis=0)).ravel()
    sd = np.sqrt(np.maximum(np.asarray(X.multiply(X).mean(axis=0)).ravel() - mean ** 2, 0.0))
    inv = np.where(sd > 0, 1.0 / np.where(sd > 0, sd, 1.0), 0.0)
    return np.asarray(X[:, cols].multiply(inv[cols]).mean(axis=1)).ravel()


def near_mask(xy_um: np.ndarray, points_um: np.ndarray, *, radius_um: float = NEAR_UM):
    from scipy.spatial import cKDTree
    tree = cKDTree(np.asarray(xy_um, dtype=float))
    m = np.zeros(len(xy_um), dtype=bool)
    for idx in tree.query_ball_point(np.asarray(points_um, dtype=float), radius_um):
        m[idx] = True
    return m, tree


def enrichment(score: np.ndarray, xy_um: np.ndarray, points_um: np.ndarray, *,
               radius_um: float = NEAR_UM, draws: int = NULL_DRAWS, seed: int = 0) -> dict:
    """Score inside the landmark neighbourhoods vs elsewhere, against a placement null."""
    inside, tree = near_mask(xy_um, points_um, radius_um=radius_um)
    if not inside.any():
        raise ValueError("no cell lies within the radius of any landmark")
    rng = np.random.default_rng(seed)
    lo, hi = np.asarray(xy_um).min(0), np.asarray(xy_um).max(0)
    nulls = []
    for _ in range(draws):
        pts = rng.uniform(lo, hi, size=np.shape(points_um))
        m, _ = near_mask(xy_um, pts, radius_um=radius_um)
        if m.any():
            nulls.append(float(score[m].mean()))
    a, b = float(score[inside].mean()), float(score[~inside].mean())
    return {"n_cells_inside": int(inside.sum()), "mean_inside": a, "mean_outside": b,
            "ratio": a / b if b else float("nan"),
            "null_mean": float(np.mean(nulls)) if nulls else float("nan"),
            "null_p95": float(np.percentile(nulls, 95)) if nulls else float("nan"),
            "null_draws": len(nulls),
            "clears_null": bool(nulls and a > np.percentile(nulls, 95))}
