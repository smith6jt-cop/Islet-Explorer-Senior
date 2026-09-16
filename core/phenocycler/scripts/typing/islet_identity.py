#!/usr/bin/env python3
"""Phase 1: is the UNDECIDED islet-core population a mixture, and do CD99 and shape see it?

WHAT THIS ASKS, AND WHY THE FIRST VERSION ASKED THE WRONG THING.

The first version compared two labelled groups -- hormone-called endocrine against everything
else in the islet -- by median and Cliff's delta, and concluded "no separation". Three things
were wrong with that framing and each of them mattered:

  * AVERAGES ANSWER THE WRONG QUESTION. Whether two groups differ on average says nothing about
    whether a given cell can be assigned. A mixture hiding inside one group is invisible to a
    median.
  * THE SETTLED CELLS DO NOT NEED DECIDING. A cell with unambiguous hormone signal is already
    an alpha, beta or delta cell. Including 185,000 of them in the comparison buries the few
    thousand whose identity is actually open, and makes the aggregate a statement about the
    settled majority.
  * THE STRUCTURAL CLASSES WERE NEVER THIS QUESTION'S BUSINESS. Blood vessel, nerve and
    lymphatic are decided by centroid-in-a-pixel-classifier-object plus that object's own
    markers. They are removed here rather than compared.

So: the REFERENCE is cells with unambiguous hormone signal. The TARGET is islet-core cells that
are neither structural nor unambiguously hormone-positive -- the population whose identity is
genuinely open. The question is whether that target is ONE population or TWO, and whether CD99
and shape are what splits it.

HOW IT IS TESTED -- dispersion and mixture, not location:

  variance ratio      spread of the target over spread of the reference, per feature. A target
                      that is a mixture of endocrine-like and non-endocrine-like cells is WIDER
                      than either component. >1 is the signature; ~1 says one population.
  bimodality          Hartigan-style dip on the per-cell Mahalanobis distance from the reference
                      centroid, plus a 1-vs-2-component Gaussian mixture compared by BIC. Two
                      components winning by BIC is the direct evidence of a mixture.
  per-cell distance   Mahalanobis in the CD99+shape space, covariance estimated on the reference.
                      Reported as a distribution, so a bimodal target is visible as two humps
                      rather than averaged into one number.

SHAPE IS NUCLEUS-FIRST. Cell-boundary features were measured across six donors and every class
-- immune, exocrine, alpha, beta, delta -- sat at cell circularity 0.514-0.523, a 1% spread.
Real immune and acinar cells are not the same shape; that is the signature of boundaries cut
against neighbours in packed tissue, so the cell outline is a tessellation artefact and carries
little morphology. Nuclei are segmented from DAPI directly and are not tessellated, so they lead.
Cell features are still reported, marked, so the claim stays checkable.

SCOPE. Only donors whose islet annotations are settled. The batch-2 donors are being re-exported
against a new tissue annotation that re-derives the islet set, so their `cell_region` answers a
question about islets that no longer exist. Refused by name, not skipped silently.

    python core/phenocycler/scripts/typing/islet_identity.py [--donors 115]
"""
from __future__ import annotations

import argparse
import pathlib
import sys

import numpy as np
import pandas as pd
import pyarrow.parquet as pq

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))
import subpanel_io as S

#: Donors whose islet identity cannot be asked yet. Empty since 2026-09-15: the batch-2 donors
#: were re-exported, rebuilt, re-typed against islet membership derived from the outlines, and
#: their missed-islet sweep returns zero at every setting. The gate the maintainer set -- newest
#: typing AND current islet geojson, per donor -- is met for all twenty.
ISLET_IDENTITY_IN_FLUX: set = set()
STRUCTURAL = {"Blood vessel", "Nerve", "Lymphatic"}
HORMONES = ["INS", "GCG", "SST"]
MARKER = "CD99"

#: Nucleus shape leads -- see the header. Cell features follow and are labelled.
NUCLEUS_SHAPE = ["nucleus_circularity", "nucleus_solidity", "nucleus_area_um2",
                 "nucleus_elongation", "nucleus_length2_over_area"]
CELL_SHAPE = ["cell_circularity", "cell_solidity", "cell_elongation"]

#: A cell is UNAMBIGUOUSLY endocrine when its top hormone clears this on the scaled [0,1] values
#: and the runner-up is well below it. Deliberately strict: this set defines what endocrine
#: LOOKS like, so admitting a doubtful cell corrupts the reference the whole test rests on.
REFERENCE_MIN = 0.35
REFERENCE_MARGIN = 2.0

#: The undecided set excludes anything with even moderate hormone signal -- those are plausibly
#: endocrine and belong to neither side cleanly.
UNDECIDED_MAX = 0.10


def derive_shape(f: pd.DataFrame) -> pd.DataFrame:
    for name, (num, den) in (("cell_elongation", ("cell_max_diam_um", "cell_min_diam_um")),
                             ("nucleus_elongation", ("nucleus_max_diam_um", "nucleus_min_diam_um"))):
        f[name] = f[num] / f[den].where(f[den] > 0)
    for name, (ln, ar) in (("nucleus_length2_over_area", ("nucleus_length_um", "nucleus_area_um2")),):
        f[name] = (f[ln] ** 2) / f[ar].where(f[ar] > 0)
    return f


def dip_statistic(x: np.ndarray) -> float:
    """Hartigan's dip, computed as the max gap between the ECDF and its greatest convex minorant.

    Small and self-contained rather than a dependency. The VALUE matters less than the
    comparison the caller makes: the same statistic on the reference gives the null for a
    population that is not a mixture.
    """
    x = np.sort(x[np.isfinite(x)])
    n = len(x)
    if n < 20:
        return float("nan")
    ecdf = np.arange(1, n + 1) / n
    # convex minorant via a monotone lower hull on (x, ecdf)
    hull = [0]
    for i in range(1, n):
        while len(hull) >= 2:
            a, b = hull[-2], hull[-1]
            if (ecdf[b] - ecdf[a]) * (x[i] - x[a]) >= (ecdf[i] - ecdf[a]) * (x[b] - x[a]):
                hull.pop()
            else:
                break
        hull.append(i)
    gcm = np.interp(x, x[hull], ecdf[hull])
    return float(np.max(np.abs(ecdf - gcm)))


def mixture_bic(x: np.ndarray, seed: int = 0) -> tuple[float, float, float]:
    """(bic_1, bic_2, weight of the smaller component) for a 1- vs 2-component Gaussian fit."""
    from sklearn.mixture import GaussianMixture
    x = x[np.isfinite(x)].reshape(-1, 1)
    if len(x) < 200:
        return float("nan"), float("nan"), float("nan")
    if len(x) > 50000:
        x = np.random.default_rng(seed).choice(x, 50000, replace=False).reshape(-1, 1)
    g1 = GaussianMixture(1, random_state=seed).fit(x)
    g2 = GaussianMixture(2, random_state=seed, n_init=3).fit(x)
    return float(g1.bic(x)), float(g2.bic(x)), float(min(g2.weights_))


def mahalanobis(target: np.ndarray, reference: np.ndarray) -> np.ndarray:
    """Per-cell distance from the reference centroid, covariance from the REFERENCE.

    Reference covariance, not pooled: the question is how far each undecided cell sits from
    what endocrine looks like, and pooling would let the target's own spread define the metric
    it is being judged by.
    """
    mu = np.nanmean(reference, axis=0)
    cov = np.cov(reference[~np.isnan(reference).any(axis=1)].T)
    cov += np.eye(cov.shape[0]) * 1e-9
    inv = np.linalg.pinv(cov)
    d = target - mu
    ok = ~np.isnan(d).any(axis=1)
    out = np.full(len(target), np.nan)
    out[ok] = np.sqrt(np.einsum("ij,jk,ik->i", d[ok], inv, d[ok]))
    return out


def load(donor: str):
    paths = {
        "cells": S.OUT / f"cells/donor_id={donor}/data.parquet",
        "shape": S.OUT / f"qupath_fundamentals/donor_id={donor}/data.parquet",
        "marker": S.ALIGNED / f"donor_id={donor}/data.parquet",
    }
    for p in paths.values():
        if not p.exists():
            return None
    have = set(pq.ParquetFile(paths["marker"]).schema_arrow.names)
    need = [MARKER] + [h for h in HORMONES if h in have]
    if MARKER not in have or len(need) < 3:
        return None

    cells = pq.read_table(
        paths["cells"],
        columns=["object_id", "cell_class", "cell_region_20um", "islet_num"]).to_pandas()
    marker = pq.read_table(paths["marker"], columns=["object_id"] + need).to_pandas()
    shape = derive_shape(pq.read_table(paths["shape"]).to_pandas())

    df = cells.merge(marker, on="object_id", how="inner", validate="one_to_one")
    df = df.merge(shape, on="object_id", how="left", validate="one_to_one")
    df = df[(df["cell_region_20um"] == "core") & (~df["cell_class"].isin(STRUCTURAL))
            & (df["cell_class"] != "Excluded")]
    if df.empty:
        return None

    # Islet size, counted over the cells this analysis actually admits -- core, non-structural,
    # not Excluded -- so it is the size of the thing being looked at, not QuPath's detection
    # count. The endocrine shape difference is a property of real islets; a two-cell classifier
    # fragment has no architecture to impose one, so size has to be available to stratify on.
    df["islet_cells"] = df.groupby("islet_num")["object_id"].transform("size")

    scale = S.load_scale()
    horm = [h for h in HORMONES if h in df.columns]
    scaled = pd.DataFrame(S.apply_scale(df[horm].to_numpy(), horm, scale),
                          columns=horm, index=df.index)
    top = scaled.max(axis=1)
    second = scaled.apply(lambda r: np.sort(r.to_numpy())[-2], axis=1)
    df["hormone_top"] = top
    df["reference"] = (top >= REFERENCE_MIN) & (top >= REFERENCE_MARGIN * second.clip(lower=1e-6))
    df["undecided"] = top < UNDECIDED_MAX
    return df


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--donors", nargs="*", default=None)
    args = ap.parse_args()

    available = sorted(p.name.split("=")[1] for p in (S.OUT / "cells").glob("donor_id=*"))
    wanted = args.donors or available
    refused = sorted(set(wanted) & ISLET_IDENTITY_IN_FLUX)
    if refused:
        print(f"REFUSED, islet identity in flux (re-exported against a new tissue "
              f"annotation): {' '.join(refused)}\n")
    donors = [d for d in wanted if d not in ISLET_IDENTITY_IN_FLUX]

    features = [MARKER] + NUCLEUS_SHAPE + CELL_SHAPE
    rows = []
    for donor in donors:
        df = load(donor)
        if df is None:
            print(f"  {donor}: inputs incomplete"); continue
        ref = df[df["reference"]]
        und = df[df["undecided"]]
        if len(ref) < 500 or len(und) < 200:
            print(f"  {donor}: reference {len(ref):,} / undecided {len(und):,} -- too few")
            continue

        print(f"\n{donor}  reference (unambiguous hormone) {len(ref):,} | "
              f"undecided {len(und):,} | other islet-core {len(df) - len(ref) - len(und):,}")
        print(f"   {'feature':<26}{'ref sd':>10}{'und sd':>10}{'var ratio':>11}"
              f"{'ref med':>10}{'und med':>10}")
        for f in features:
            if f not in df.columns:
                continue
            r = ref[f].to_numpy(dtype=float); u = und[f].to_numpy(dtype=float)
            r = r[np.isfinite(r)]; u = u[np.isfinite(u)]
            if len(r) < 100 or len(u) < 100:
                continue
            vr = float(np.var(u) / np.var(r)) if np.var(r) > 0 else np.nan
            tag = "" if f in (MARKER,) + tuple(NUCLEUS_SHAPE) else "  [cell]"
            print(f"   {f:<26}{np.std(r):>10.3f}{np.std(u):>10.3f}{vr:>11.2f}"
                  f"{np.median(r):>10.3f}{np.median(u):>10.3f}{tag}")
            rows.append({"donor_id": donor, "feature": f, "var_ratio": vr,
                         "n_ref": len(r), "n_und": len(u),
                         "median_ref": float(np.median(r)), "median_und": float(np.median(u))})

        space = [MARKER] + [f for f in NUCLEUS_SHAPE if f in df.columns]
        R = ref[space].to_numpy(dtype=float); U = und[space].to_numpy(dtype=float)
        d_und = mahalanobis(U, R); d_ref = mahalanobis(R, R)
        b1, b2, w = mixture_bic(d_und)
        rb1, rb2, rw = mixture_bic(d_ref)
        print(f"   Mahalanobis in ({', '.join(space)}):")
        print(f"     undecided : median {np.nanmedian(d_und):6.2f}  dip {dip_statistic(d_und):.4f}"
              f"  BIC 1-comp {b1:,.0f} vs 2-comp {b2:,.0f}"
              f"  -> {'TWO components' if b2 < b1 else 'one'}"
              + (f", smaller weight {w:.3f}" if b2 < b1 else ""))
        print(f"     reference : median {np.nanmedian(d_ref):6.2f}  dip {dip_statistic(d_ref):.4f}"
              f"  BIC 1-comp {rb1:,.0f} vs 2-comp {rb2:,.0f}"
              f"  -> {'TWO components' if rb2 < rb1 else 'one'}   (the null)")

    if rows:
        out = S.OUT / "islet_identity_separation.parquet"
        pd.DataFrame(rows).to_parquet(out, index=False)
        agg = (pd.DataFrame(rows).groupby("feature")["var_ratio"]
                 .median().sort_values(ascending=False))
        print("\n" + "=" * 74)
        print("VARIANCE RATIO across donors (undecided spread / reference spread).")
        print(">1 means the undecided population is WIDER than the settled endocrine one,")
        print("which is what a mixture of two cell types looks like on that feature.\n")
        print(agg.to_string(float_format=lambda v: f"{v:6.2f}"))
        print(f"\n-> {out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
