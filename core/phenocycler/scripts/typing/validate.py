#!/usr/bin/env python3
"""Checks on a finished call set, none of which the workflow can tune itself against.

Every check here compares the calls to something the typing never saw: the QuPath islet polygons
(geometry, drawn by hand), and each cell's spatial neighbours.  None of them is a threshold the
workflow optimises, and none compares a donor to an expected frequency -- the numbers are reported
per donor, against that donor's own data, and are for the reviewer to weigh.

What each check is for:

``islet composition``
    Islet cores should be overwhelmingly Endocrine and nearly free of Exocrine.  The annotation is
    geometric, so this is independent of every marker the call used.
``endocrine placement``
    Where the Endocrine calls physically sit.  A typing that scatters them through acinar tissue
    is wrong however good its marker statistics look.
``spatial coherence``
    Cells of one type sit together.  Measured as the fraction of each cell's 12 nearest neighbours
    sharing its label, against a null that shuffles labels within 500 um blocks -- local, so it
    destroys cell-scale structure while preserving where the tissue's compartments are.

    READ THE RATIO, NOT ``clears``.  At a million cells the null is narrow enough that every
    label clears its p95, so the boolean separates nothing; the ratio of observed to null is what
    carries the information.  Two things bound how it should be read.  A label is not credible
    below the ratio the ``Ambiguous`` bucket achieves -- on donor 6442 that is 1.43, and Neural
    came in at 1.26, i.e. less spatially organised than the cells the workflow declined to name,
    for a lineage whose defining structure is nerve bundles.  And the test has NO POWER for an
    abundant label: within a 500 um block most cells genuinely are exocrine, so the null already
    contains mostly exocrine and the ratio is compressed toward 1 by construction (Exocrine 1.22,
    Myeloid 1.13).  Those say nothing either way.
``anchor share``
    How much of each lineage rests on a cell's own marker positivity rather than on the profile of
    its cluster.  A lineage carried entirely by inference is one to look at first in QuPath.
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.neighbors import NearestNeighbors

sys.path.insert(0, str(Path(__file__).parent))
import lineages as L  # noqa: E402

REPO = Path(__file__).resolve().parents[4]
DEV = REPO / "data/development/typing"
K_NEIGHBOURS = 12
BLOCK_UM = 500.0


def coherence(df: pd.DataFrame, labels: np.ndarray, seed: int, draws: int) -> pd.DataFrame:
    xy = df[["X_centroid", "Y_centroid"]].to_numpy(float)
    ok = np.isfinite(xy).all(1) & ~pd.isna(labels)
    xy, lab = xy[ok], np.asarray(labels)[ok]
    nn = NearestNeighbors(n_neighbors=K_NEIGHBOURS + 1).fit(xy)
    _, ind = nn.kneighbors(xy)
    ind = ind[:, 1:]
    same = (lab[ind] == lab[:, None]).mean(1)

    block = (xy // BLOCK_UM).astype(np.int64)
    key = block[:, 0] * 1_000_003 + block[:, 1]
    rng = np.random.default_rng(seed)
    order = np.argsort(key, kind="stable")
    bounds = np.flatnonzero(np.diff(key[order])) + 1
    groups = np.split(order, bounds)

    rows = []
    null = np.empty((draws, len(np.unique(lab))))
    uniq = np.unique(lab)
    for d in range(draws):
        perm = lab.copy()
        for g in groups:
            if g.size > 1:
                perm[g] = lab[rng.permutation(g)]
        s = (perm[ind] == perm[:, None]).mean(1)
        null[d] = [s[perm == u].mean() if (perm == u).any() else np.nan for u in uniq]
    for i, u in enumerate(uniq):
        obs = same[lab == u].mean()
        p95 = float(np.nanpercentile(null[:, i], 95))
        rows.append({"label": u, "n": int((lab == u).sum()), "observed": round(float(obs), 4),
                     "null_p95": round(p95, 4),
                     "ratio": round(float(obs / p95), 2) if p95 > 0 else np.nan,
                     "clears": bool(obs > p95)})
    return pd.DataFrame(rows).sort_values("ratio", ascending=False)


def islet_composition(donor: str) -> pd.DataFrame:
    """Endocrine composition per islet, which is what survives the segmentation.

    Beta-versus-alpha is NOT reliably resolvable per cell on this segmentation and the export says
    so: 31.2% of donor 6442's islet-core cells carry both INS and GCG above the model's own
    boundary, and a cell cannot be both.  The cause is geometric, not a model failure -- Nimbus
    tracks the data (INS confidence correlates +0.805 with measured INS intensity) but islet cells
    are packed tightly and the polygons are cytoplasm-poor (median cytoplasm 19-20 um2, N:C 0.67
    against 0.41-0.45 for the Xenium segmentation on the adjacent section), so a beta cell's mask
    contains its alpha neighbour's glucagon.  An argmax between two contaminated signals is close
    to a coin flip, which is why alpha edges out beta on an ND donor even though INS is four times
    brighter in raw intensity (median 2,418 against 560 in islet cores).

    Averaging over an islet helps but does NOT cancel, and the residue has a direction that
    matters.  Measured on donor 6442 over 293 islets with at least ten endocrine cells:
    corr(log islet size, double-positive fraction) = +0.52, and corr(log islet size, INS-dominant
    fraction) = -0.39.  Bigger islets are packed more densely, bleed more, and read as less beta.

    **This is a trap for the disease comparison, so it is stated here rather than discovered
    later.**  T1D islets are smaller and atrophic.  Compared naively against ND islets, the size
    artefact alone would make them look MORE beta-rich -- the opposite of the biology, with a
    convincing effect size.  Any cross-donor comparison of islet composition must condition on
    islet size (match the size distributions, or model size as a covariate), and a comparison that
    does not is measuring segmentation.

    Reported per islet with its cell count, never pooled into one number for a donor.
    """
    calls = pd.read_parquet(DEV / f"calls/donor_id={donor}/data.parquet")
    pos = pd.read_parquet(DEV / f"positivity/donor_id={donor}/data.parquet",
                          columns=["object_id", "INS", "GCG", "SST"])
    d = calls.merge(pos, on="object_id", how="left")
    # islet_num carries an empty STRING for cells outside any islet, not NaN, so notna() alone
    # groups every non-islet endocrine cell into one giant pseudo-islet (2,533 cells on 6442).
    isl_id = d["islet_num"].astype("string").str.strip()
    d = d[(d["lineage"] == "Endocrine") & isl_id.notna() & (isl_id != "") & (isl_id != "nan")]
    if d.empty:
        return pd.DataFrame()
    h = d[["INS", "GCG", "SST"]].to_numpy(float)
    dom = np.array(["INS", "GCG", "SST"])[np.nanargmax(np.nan_to_num(h, nan=-1), axis=1)]
    d = d.assign(dominant=dom,
                 double_positive=(d["INS"] > 0.5) & (d["GCG"] > 0.5))
    g = d.groupby("islet_num")
    out = pd.DataFrame({
        "n_endocrine": g.size(),
        "frac_INS_dominant": g["dominant"].apply(lambda v: float((v == "INS").mean())).round(3),
        "frac_GCG_dominant": g["dominant"].apply(lambda v: float((v == "GCG").mean())).round(3),
        "frac_SST_dominant": g["dominant"].apply(lambda v: float((v == "SST").mean())).round(3),
        "frac_double_positive": g["double_positive"].mean().round(3),
    }).reset_index()
    return out[out["n_endocrine"] >= 10].sort_values("n_endocrine", ascending=False)


def report(donor: str, seed: int, draws: int) -> dict:
    df = pd.read_parquet(DEV / f"calls/donor_id={donor}/data.parquet")
    n = len(df)
    out: dict = {"donor_id": donor, "n_objects": n}

    print(f"\n{'='*82}\nDONOR {donor}   {n:,} objects")
    print("\n-- composition (all objects) --")
    for v, c in df["lineage"].value_counts().items():
        print(f"   {v:<14}{c:>9,}{c/n*100:>7.1f}%")
    out["lineage_counts"] = df["lineage"].value_counts().to_dict()

    core = df[df["cell_region"] == "core"]
    if len(core):
        print(f"\n-- islet cores ({len(core):,} cells; geometry the typing never saw) --")
        comp = (core["lineage"].value_counts() / len(core) * 100).round(1)
        for v, c in comp.items():
            print(f"   {v:<14}{c:>7.1f}%")
        out["islet_core_pct"] = comp.to_dict()

    for lin in ("Endocrine", "Exocrine"):
        sel = df[df["lineage"] == lin]
        if len(sel):
            pct = (sel["cell_region"].value_counts() / len(sel) * 100).round(1).to_dict()
            print(f"\n-- where {lin} calls sit --  " +
                  "  ".join(f"{k}={v}%" for k, v in pct.items()))
            out[f"{lin}_placement_pct"] = pct

    print(f"\n-- spatial coherence: share of {K_NEIGHBOURS} nearest neighbours sharing the label,")
    print(f"   against labels shuffled within {BLOCK_UM:.0f} um blocks --")
    coh = coherence(df, df["lineage"].to_numpy(), seed, draws)
    print(coh.to_string(index=False))
    out["coherence"] = coh.to_dict("records")

    isl = islet_composition(donor)
    if len(isl):
        print(f"\n-- endocrine composition per islet ({len(isl)} islets with >=10 endocrine "
              f"cells; the per-cell beta/alpha call is NOT reliable, see islet_composition) --")
        print(f"   frac_INS_dominant  median {isl.frac_INS_dominant.median():.2f}  "
              f"IQR {isl.frac_INS_dominant.quantile(.25):.2f}-{isl.frac_INS_dominant.quantile(.75):.2f}")
        print(f"   frac_GCG_dominant  median {isl.frac_GCG_dominant.median():.2f}")
        print(f"   frac_double_positive median {isl.frac_double_positive.median():.2f}"
              f"  <- how much of the islet cannot be resolved per cell")
        isl.to_csv(DEV / f"calls/donor_id={donor}/islet_composition.csv", index=False)
        out["n_islets"] = int(len(isl))
        out["islet_frac_INS_dominant_median"] = float(isl.frac_INS_dominant.median())
        out["islet_frac_double_positive_median"] = float(isl.frac_double_positive.median())

    print("\n-- subtype composition --")
    for v, c in df["cell_type"].value_counts().items():
        print(f"   {v:<16}{c:>9,}{c/n*100:>7.1f}%")
    out["cell_type_counts"] = df["cell_type"].value_counts().to_dict()
    return out


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--donors", nargs="+", required=True)
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--null-draws", type=int, default=50)
    a = ap.parse_args()
    reports = [report(d, a.seed, a.null_draws) for d in a.donors]
    dest = DEV / "validation.json"
    dest.write_text(json.dumps(reports, indent=2, default=str))
    print(f"\nwrote {dest}")


if __name__ == "__main__":
    main()
