"""Protein INS and GCG on the Xenium cells that were measured for RNA.

This is the one measurement in the cohort where protein and RNA are on ONE CELL, and it is
the only check that exists on Xenium's endocrine tier -- which rests entirely on surrogates,
because INS, GCG, SST and PPY are all off the 5K panel. ``PDX1``, one of its five beta
surrogates, is detected in 1.18% of cells.

**It derives no threshold.** The question is a ranking:

    does protein INS, measured on the same cell, rank the cells whose RNA looks beta above
    the cells whose RNA looks alpha -- and does GCG rank them the other way?

That is an AUROC, so no hormone cutoff is needed on the protein side and none is invented.
And the RNA side needs no cutoff either: among cells already inside an islet, the label is
whichever surrogate set reads higher *in that cell*, which is a contrast, not a gate.

**Restricting to islet residents is what makes it a test of SUBTYPE.** Run over the whole
section, both hormones would separate islet from non-islet and both AUROCs would look
excellent while saying nothing about alpha versus beta. Inside islets that free signal is
gone: a registration that merely lands on islets scores 0.5 here.

**The channel map is not in the image.** The qptiff carries channel NAMES
(``DAPI / Opal 570 / Opal 690 / Opal 520 / Sample AF``) and no ``<Biomarker>`` tag, so which
Opal is which antibody exists only in ``data/integration/restain_sections.csv``, hand-recorded.
``Opal 570`` is recorded ``FAILED`` and this module refuses to read a channel so marked rather
than leaving it available to a caller who does not know.

**Autofluorescence: reported both ways, and the default follows the instruction.** ``Sample
AF`` is recorded ``AF`` and the maintainer's instruction is to disregard it, so the headline
numbers are raw. But measured on islet-resident cells it is not inert, and that is recorded
rather than argued:

    donor 6442: INS~GCG rho  +0.085 raw  ->  -0.006 after regressing AF out
    donor 6436: INS~GCG rho  -0.387 raw  ->  -0.305 after regressing AF out
    Sample AF ~ Opal 690 (INS): rho 0.50 (6442), 0.44 (6436)
    Opal 570 (FAILED) ~ Opal 690: rho 0.71 (6442), 0.29 (6436)

Two hormones on different cells should read NEGATIVELY correlated. In 6442 the raw value is
positive and only the AF regression removes that, which is what shared background does. And
``Opal 570`` is recorded FAILED -- so whatever it holds is background, and Opal 690 sharing a
0.71 rank correlation with it in that section is a caution about the INS channel there, not a
result.

**``Sample AF`` is CLIPPED at 255 in 18-29% of islet-resident cells**, so a linear AF
correction is not even well defined for those cells: the covariate is censored exactly where
it is largest. Both facts are emitted by :func:`channel_diagnostics` beside every AUROC, so a
reader never has to ask.
"""

from __future__ import annotations

import re
from pathlib import Path

import numpy as np
import pandas as pd

#: Radius around a warped cell centroid over which a channel is averaged. The cohort's median
#: Xenium cell area is ~50 um^2, so ~4 um is a cell radius: big enough to average over the
#: cell rather than one pixel, small enough not to reach a neighbour. Polygon-exact extraction
#: is available via `integration/boundaries.py` and would matter for an intensity threshold;
#: for a RANK statistic the window is second-order.
SAMPLE_RADIUS_UM = 4.0

#: Channel roles that are never read, whatever a caller asks for.
UNUSABLE = {"FAILED", "AF"}


class HormoneError(RuntimeError):
    """A channel cannot be read as the marker it is claimed to be."""


def load_section_map(csv_path: Path) -> pd.DataFrame:
    df = pd.read_csv(csv_path, comment="#", dtype=str)
    required = ["image_file", "donor_id", "roi", "xenium_serial", "pixel_size_um"]
    for col in required:
        if df[col].isna().any() or (df[col].astype(str).str.strip() == "").any():
            raise HormoneError(
                f"{csv_path} has a blank {col!r}. This file is hand-recorded because the "
                "facts are in no image; a blank is refused, never guessed.")
    return df


def channel_index(qptiff: Path, marker: str, row: pd.Series) -> int:
    """Index of the channel that stained ``marker``, matched by the qptiff's own names.

    The index is DERIVED from the file's channel names rather than hardcoded, so a section
    acquired with a different channel order does not silently return someone else's stain.
    """
    role_col = {v: k for k, v in row.items()
                if isinstance(k, str) and k.startswith("channel_")}
    if marker not in role_col:
        raise HormoneError(f"{marker!r} is not recorded for {row['image_file']}; "
                           f"recorded: {sorted(v for v in role_col if v)}")
    col = role_col[marker]
    if str(row[col]).strip().upper() in UNUSABLE:
        raise HormoneError(f"{col} is recorded {row[col]!r} for {row['image_file']} "
                           "and must not be read")
    return _page_for_column(qptiff, col)


def _page_for_column(qptiff: Path, col: str) -> int:
    """Page index whose qptiff ``<Name>`` matches the registry column ``col``."""
    want = col.replace("channel_", "")                 # e.g. 'opal690'
    names = channel_names(qptiff)
    norm = [re.sub(r"[^a-z0-9]", "", n.lower()) for n in names]
    for i, n in enumerate(norm):
        if n == want:
            return i
    raise HormoneError(f"{qptiff.name} has channels {names}; none matches {want!r}")


def channel_index_for_role(qptiff: Path, row: pd.Series, role: str) -> int | None:
    """Page index of the channel recorded with role ``role`` (e.g. ``"AF"``), or ``None``.

    Separate from :func:`channel_index` because ``UNUSABLE`` means "never read as a MARKER",
    not "never read at all": the AF channel is a covariate and the FAILED one is a diagnostic,
    and both must still be found the way every other channel is -- by the qptiff's own
    ``<Name>`` XML, never by position.

    This replaces a zip of the file's page order against a hardcoded list of registry columns.
    That was correct only for the five sections recorded so far, and the first section
    acquired in a different channel order would have regressed the wrong covariate out with no
    error and no test to catch it.
    """
    hits = [k for k, v in row.items()
            if isinstance(k, str) and k.startswith("channel_")
            and str(v).strip().upper() == role.strip().upper()]
    if not hits:
        return None
    if len(hits) > 1:
        raise HormoneError(
            f"{row['image_file']} records role {role!r} on {len(hits)} channels ({sorted(hits)}); "
            "the registry is hand-recorded and an ambiguous role is refused, never picked.")
    return _page_for_column(qptiff, hits[0])


def channel_names(qptiff: Path) -> list[str]:
    """Channel names in page order, from the qptiff's own XML."""
    import tifffile
    with tifffile.TiffFile(str(qptiff)) as tif:
        out = []
        for page in tif.series[0].pages:
            m = re.search(r"<Name>([^<]*)</Name>", page.description or "")
            out.append(m.group(1) if m else "")
    return out


def sample_cells(image: np.ndarray, um_per_px: float, xy_um: np.ndarray, *,
                 radius_um: float = SAMPLE_RADIUS_UM) -> np.ndarray:
    """Mean channel value in a disk around each point, in the IMAGE's own frame.

    Points outside the image return NaN rather than 0 -- a cell that fell off the edge is
    unavailable, not dark, and averaging a 0 into a rank would move it to the bottom.
    """
    h, w = image.shape
    xy_px = np.asarray(xy_um, dtype=float) / um_per_px
    r = max(1, int(round(radius_um / um_per_px)))
    dy, dx = np.mgrid[-r:r + 1, -r:r + 1]
    disk = (dx ** 2 + dy ** 2) <= r ** 2
    oy, ox = np.nonzero(disk)
    oy, ox = oy - r, ox - r

    col = np.rint(xy_px[:, 0]).astype(np.int64)
    row = np.rint(xy_px[:, 1]).astype(np.int64)
    out = np.full(len(xy_px), np.nan, dtype=np.float64)
    inside = (col >= r) & (col < w - r) & (row >= r) & (row < h - r)
    if not inside.any():
        return out
    rr = row[inside][:, None] + oy[None, :]
    cc = col[inside][:, None] + ox[None, :]
    out[inside] = image[rr, cc].mean(axis=1)
    return out


def regress_out(values: np.ndarray, covariate: np.ndarray) -> np.ndarray:
    """``values`` with the linear part explained by ``covariate`` removed.

    Used for the autofluorescence arm only, and never as the default -- see the module note.
    Rows where either is non-finite return NaN rather than 0, because a cell whose covariate
    is unavailable has an unavailable residual, not a zero one.
    """
    v = np.asarray(values, dtype=float)
    c = np.asarray(covariate, dtype=float)
    ok = np.isfinite(v) & np.isfinite(c)
    out = np.full(v.shape, np.nan)
    if ok.sum() < 3:
        return out
    A = np.vstack([np.ones(ok.sum()), c[ok]]).T
    coef, *_ = np.linalg.lstsq(A, v[ok], rcond=None)
    out[ok] = v[ok] - A @ coef
    return out


def channel_diagnostics(per_channel: dict[str, np.ndarray]) -> dict:
    """What every hormone AUROC must be read beside.

    Rank correlations between all channels (two hormones on different cells should be
    NEGATIVE), and the saturated fraction of each -- a channel clipped at its ceiling is
    censored, and a covariate that is censored exactly where it is largest cannot be
    regressed out honestly.
    """
    import itertools
    names = list(per_channel)
    df = pd.DataFrame({k: np.asarray(v, dtype=float) for k, v in per_channel.items()}).dropna()
    rho = {f"{a}~{b}": float(df[a].corr(df[b], method="spearman"))
           for a, b in itertools.combinations(names, 2)} if len(df) > 2 else {}
    sat = {k: float(np.nanmean(np.asarray(v, dtype=float) >= 255.0))
           for k, v in per_channel.items()}
    return {"n_cells": int(len(df)), "spearman": rho, "saturated_fraction_at_255": sat}


def agreement_by_scale(xy_um: np.ndarray, rna: np.ndarray, protein: np.ndarray, *,
                       radii_um=(0, 10, 20, 30, 50, 80, 120)) -> list[dict]:
    """Spearman(RNA contrast, protein contrast) as the neighbourhood is widened.

    This is what says how to READ a cell-level AUROC, and it distinguishes two very different
    situations that a single number cannot:

      * a fixed registration OFFSET -- agreement climbs and then PLATEAUS once the radius
        exceeds the offset, because averaging has absorbed the misplacement;
      * a noisy or structure-level correspondence -- agreement keeps climbing past the size of
        the structures themselves, which means what agrees is islet-to-islet, not cell-to-cell.

    Measured on islet-resident cells: 6442 climbs monotonically 0.074 -> 0.287 and is still
    climbing at 120 um, which is wider than an islet; 6436 is flat at 0.03-0.05 at every scale.
    Neither is the plateau an offset would give, so neither donor's weak per-cell AUROC is
    explained by a constant misalignment.
    """
    from scipy.spatial import cKDTree
    xy = np.asarray(xy_um, dtype=float)
    rna = np.asarray(rna, dtype=float)
    protein = np.asarray(protein, dtype=float)
    tree = cKDTree(xy)
    out = []
    for r in radii_um:
        if r <= 0:
            a, b, k = rna, protein, 1.0
        else:
            idx = tree.query_ball_point(xy, float(r))
            k = float(np.mean([len(i) for i in idx]))
            a = np.array([rna[i].mean() for i in idx])
            b = np.array([protein[i].mean() for i in idx])
        ok = np.isfinite(a) & np.isfinite(b)
        rho = (float(pd.Series(a[ok]).corr(pd.Series(b[ok]), method="spearman"))
               if ok.sum() > 2 else float("nan"))
        out.append({"radius_um": float(r), "mean_neighbours": k, "spearman": rho})
    return out


def auroc(score: np.ndarray, positive: np.ndarray) -> float:
    """Rank-based AUROC. No threshold anywhere, which is the point.

    Ties are handled by average ranks, so a channel that is flat over most cells scores 0.5
    rather than being rewarded for its tie structure.
    """
    score = np.asarray(score, dtype=float)
    positive = np.asarray(positive, dtype=bool)
    ok = np.isfinite(score)
    score, positive = score[ok], positive[ok]
    n1, n0 = int(positive.sum()), int((~positive).sum())
    if n1 == 0 or n0 == 0:
        return float("nan")
    ranks = pd.Series(score).rank(method="average").to_numpy()
    return float((ranks[positive].sum() - n1 * (n1 + 1) / 2.0) / (n1 * n0))


def auroc_with_null(score: np.ndarray, positive: np.ndarray, *, draws: int = 200,
                    seed: int = 0) -> dict:
    """AUROC beside the null that destroys the PAIRING but keeps both distributions.

    Permuting the RNA labels among the same paired cells leaves the protein values and the
    label frequencies exactly as they are, so anything the null cannot reach is the
    correspondence itself rather than a property of either measurement.
    """
    real = auroc(score, positive)
    rng = np.random.default_rng(seed)
    pos = np.asarray(positive, dtype=bool).copy()
    nulls = []
    for _ in range(draws):
        rng.shuffle(pos)
        nulls.append(auroc(score, pos))
    nulls = np.asarray([v for v in nulls if np.isfinite(v)])
    return {"auroc": real, "n_positive": int(np.asarray(positive).sum()),
            "n_negative": int((~np.asarray(positive, dtype=bool)).sum()),
            "null_mean": float(nulls.mean()) if len(nulls) else float("nan"),
            "null_p2_5": float(np.percentile(nulls, 2.5)) if len(nulls) else float("nan"),
            "null_p97_5": float(np.percentile(nulls, 97.5)) if len(nulls) else float("nan"),
            "null_draws": int(len(nulls)),
            "clears_null": bool(len(nulls) and real > np.percentile(nulls, 97.5))}
