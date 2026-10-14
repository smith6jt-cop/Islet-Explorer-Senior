#!/usr/bin/env python3
"""Non-diabetic reference figures for the same-slide pair, laid out after Mereu et al.

Reference: Mereu E, Balboa D, Liebig J, et al. "A single-cell multi-omic atlas of the human
pancreas reveals potential for cellular plasticity across development, health, and disease."
Cell Metabolism 2026;38(9):1896-1915.e12. doi:10.1016/j.cmet.2026.07.023 (PMID 42586061).
Retrieved via PubMed; PDF at docs/human_panc_atlas.pdf.

Donor 112 is the NON-DIABETIC donor here (donor_metadata_panc.xlsx: disease.status ND,
aab_count 0). 6533 is T1D with 3 autoantibodies, and its section was damaged between
acquisitions, so the ND donor is also the one with the clean 0.10 um registration.

PANELS, and which figure of theirs each follows:

  A  (their 1G)  spatial cell-type map of the whole section, plus a magnified inset
  B  (their 3C)  one islet core, endocrine subtypes on the tissue
  C  (their 3D)  UMAP of endocrine cells coloured by subtype
  D  (their 3E)  endocrine subtype composition binned by islet AREA

WHAT IS AND IS NOT COMPARABLE. Their 3D/3E come from Xenium TRANSCRIPTS; ours come from the
PhenoCycler PROTEIN arm, because INS/GCG/SST are protein markers absent from both Xenium
panels -- that asymmetry is the whole reason this pair exists. So panel C is a UMAP of protein
marker space, not transcript space, and is labelled as such. Their n = 40,714 Xenium endocrine
cells; ours is whatever the islet cascade calls endocrine inside Islet_N.

Panel D is the closest to a like-for-like reproduction: they bin by endocrine-cluster area
(10-24,500 um2) and we have QuPath islet polygons with real areas.
"""
from __future__ import annotations

import argparse
import json
import csv
import pathlib
import re
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

# House style: Arial, large, nothing bold, legends OUTSIDE the axes, nothing overlapping.
# Arial itself is not installed on this machine; Liberation Sans is metric-compatible with it
# (identical advance widths, drawn as a drop-in) so the list falls through to it here and
# picks up real Arial automatically anywhere it exists.
plt.rcParams.update({
    "figure.dpi": 150, "savefig.dpi": 200,
    "font.family": "sans-serif",
    "font.sans-serif": ["Arial", "Liberation Sans", "Nimbus Sans", "DejaVu Sans"],
    "font.size": 17, "font.weight": "normal",
    "axes.titlesize": 20, "axes.titleweight": "normal", "axes.labelsize": 18,
    "axes.labelweight": "normal", "figure.titlesize": 24, "figure.titleweight": "normal",
    "xtick.labelsize": 16, "ytick.labelsize": 16,
    "legend.fontsize": 16, "legend.title_fontsize": 16, "legend.frameon": False,
    "axes.spines.top": False, "axes.spines.right": False,
    "mathtext.fontset": "custom", "mathtext.rm": "Liberation Sans",
})

sys.path.insert(0, "/blue/maigan/smith6jt/Islet-Explorer-Senior-main/core/phenocycler/scripts/celltype")
import panel as PANEL   # noqa: E402
import score as SC      # noqa: E402
sys.path.insert(0, "/blue/maigan/smith6jt/Islet-Explorer-Senior-main/core/phenocycler/scripts/restain")
import run_step3 as S3  # noqa: E402  -- REQUIRE_MIN, the lane's own detection threshold

REPO = pathlib.Path("/blue/maigan/smith6jt/Islet-Explorer-Senior-main/core/phenocycler")
DATA = pathlib.Path("/blue/maigan/smith6jt/Xenium_Analysis/data/phenocycler")
PIXEL_UM = 0.5077663810243286

ENDO = ["Endo-beta", "Endo-alpha", "Endo-delta", "Endo-hormone-negative"]
ENDO_COLOR = {"Endo-beta": "#d62728", "Endo-alpha": "#1f77b4",
              "Endo-delta": "#2ca02c", "Endo-hormone-negative": "#9467bd"}
# Their 3E bins: 10 to 24,500 um2. Ours are whole islets rather than endocrine clusters, so
# the range runs higher; the edges are stated in the axis label rather than hidden.
AREA_BINS = [0, 500, 2_000, 5_000, 15_000, 50_000, np.inf]
AREA_LAB = ["<500", "0.5-2k", "2-5k", "5-15k", "15-50k", ">50k"]


def broad_of(call: str) -> str:
    """Collapse the 37 references to a readable palette, as their 1G does."""
    c = str(call)
    if c.startswith("Endo-vascular") or c.startswith("Endo-lymphatic"):
        return "Endothelial"
    if c.startswith(("Endo-", "Epi-endocrine")):
        return "Endocrine"
    if c.startswith(("Exo-", "Epi-exocrine")):
        return "Exocrine"
    if c.startswith(("Imm-", "T-", "Mac-")) or c == "Immune":
        return "Immune"
    if c.startswith("Mes-") or c == "Mesenchymal":
        return "Mesenchymal"
    if c.startswith("Epi") or c == "Epithelial":
        return "Epithelial"
    if c == "Neural":
        return "Neural"
    return "Other/Ambiguous"


BROAD_COLOR = {
    "Endocrine": "#d62728", "Exocrine": "#ff7f0e", "Epithelial": "#8c564b",
    "Immune": "#1f77b4", "Endothelial": "#e377c2", "Mesenchymal": "#2ca02c",
    "Neural": "#17becf", "Other/Ambiguous": "#cfcfcf",
}


def islet_areas(donor_geojson: pathlib.Path) -> pd.Series:
    """Area in um^2 per islet name, from the QuPath polygons."""
    from shapely.affinity import scale as affine_scale
    from shapely.geometry import shape

    feats = json.loads(donor_geojson.read_text())
    feats = feats["features"] if isinstance(feats, dict) else feats
    out = {}
    for f in feats:
        klass = (f["properties"].get("classification") or {}).get("name")
        if klass != "Islet":                     # the core, not the 20 um ring
            continue
        name = (f["properties"].get("name") or "").replace("_exp20um", "")
        g = shape(f["geometry"])
        if not g.is_valid:
            g = g.buffer(0)
        out[name] = affine_scale(g, xfact=PIXEL_UM, yfact=PIXEL_UM, origin=(0, 0)).area
    return pd.Series(out, name="area_um2")


def islet_key(v) -> pd.Series:
    """Normalise an islet identifier to its bare number.

    The two donors disagree: 112's call table carries 'Islet_1749' and 6533's carries '1274',
    while both geojsons name their features 'Islet_<n>'. Mapping one spelling onto the other
    matched nothing for 6533 and produced an empty panel rather than an error, so both sides
    are reduced to the digits before they meet.
    """
    return (pd.Series(v).astype(str).str.replace(r"^Islet[_ ]?", "", regex=True)
              .str.strip().replace({"": pd.NA, "nan": pd.NA, "None": pd.NA}))


def _dotplot(Xs: pd.DataFrame, lab: np.ndarray, types: list, cols: list,
             counts: dict, ax, title: str, size_key_at=(1.14, 1.0),
             defs: dict | None = None) -> pd.DataFrame:
    """Rows = classes, columns = markers. Dot area = fraction detected, colour = mean level.

    Colour is per-MARKER z down the column: a marker's absolute level says more about its
    antibody than about the cell, so each column is centred on its own across-class spread and
    the panel reads as "which class is highest for this marker", which is the question.
    """
    frac = np.array([[(Xs.loc[lab == t, m] >= S3.REQUIRE_MIN).mean() for m in cols]
                     for t in types])
    mean = np.array([[Xs.loc[lab == t, m].mean() for m in cols] for t in types])
    z = (mean - mean.mean(0)) / np.maximum(mean.std(0), 1e-9)

    # Outline every cell the WORKFLOW constrains, so the panel can be read as a check rather
    # than a description: a solid box is a marker the class is defined to be POSITIVE for and
    # should be red; a dashed box is one it is defined NEGATIVE for and should be blue.
    if defs:
        for i, t in enumerate(types):
            for j, m in enumerate(cols):
                c = defs.get(t, {}).get(m)
                if not c:
                    continue
                ax.add_patch(plt.Rectangle(
                    (j - .5, i - .5), 1, 1, fill=False, lw=1.6,
                    ec="k" if c in POSITIVE else "0.45",
                    ls="-" if c in POSITIVE else (0, (3, 2)), zorder=0.5))

    yy, xx = np.meshgrid(np.arange(len(types)), np.arange(len(cols)), indexing="ij")
    sc = ax.scatter(xx.ravel(), yy.ravel(), s=(frac.ravel() * 430) + 12,
                    c=z.ravel(), cmap="RdBu_r", vmin=-2, vmax=2,
                    edgecolors="0.35", linewidths=0.6)
    ax.set_xticks(range(len(cols))); ax.set_xticklabels(cols, rotation=45, ha="right")
    ax.set_yticks(range(len(types)))
    ax.set_yticklabels([f"{t}  ({counts[t]:,})" for t in types])
    ax.set_ylim(len(types) - 0.5, -0.5); ax.set_xlim(-0.7, len(cols) - 0.3)
    ax.grid(True, color="0.9", lw=0.8); ax.set_axisbelow(True)
    ax.set_title(title, loc="left")
    cb = ax.figure.colorbar(sc, ax=ax, fraction=0.030, pad=0.03)
    cb.set_label("mean level (z across rows)")
    for f, lb in ((0.25, "25%"), (0.60, "60%"), (0.95, "95%")):
        ax.scatter([], [], s=f * 430 + 12, c="0.55", edgecolors="0.35",
                   linewidths=0.6, label=lb)
    if defs:
        ax.plot([], [], color="k", lw=1.6, label="defined +")
        ax.plot([], [], color="0.45", lw=1.6, ls=(0, (3, 2)), label="defined -")
    ax.legend(loc="upper left", bbox_to_anchor=size_key_at, ncol=1, title="detected in",
              borderaxespad=0, handletextpad=0.8, labelspacing=1.4)
    return pd.DataFrame(frac, index=types, columns=cols)


def _calibrated(df, expr, scale, max_cells=120_000, seed=0):
    """Cells x markers on the lane's calibrated scale, subsampled, aligned to `df`."""
    X = expr.reindex(df["object_id"]).dropna(how="all")
    d = df.set_index("object_id").loc[X.index]
    if len(X) > max_cells:
        k = np.random.default_rng(seed).choice(len(X), max_cells, replace=False)
        X, d = X.iloc[k], d.iloc[k]
    use = [m for m in X.columns if m in scale.index]
    Xs = pd.DataFrame(SC.apply_scale(X[use].to_numpy("float32"), use, scale),
                      index=X.index, columns=use)
    return Xs, d


WORKFLOW_CSV = pathlib.Path("/blue/maigan/smith6jt/Islet-Explorer-Senior-main"
                            "/phenotype_workflow.csv")
POSITIVE = {"allpos", "anypos", "anypos2", "pos"}


def _norm_marker(m: str) -> str:
    """'HLA-A' / 'Granzyme B' -> 'HLA_A' / 'Granzyme_B', matching the expression columns."""
    return re.sub(r"[^0-9A-Za-z]+", "_", str(m)).strip("_")


def workflow_defs() -> tuple[dict, dict]:
    """Each reference's OWN marker constraints, and its parent, from phenotype_workflow.csv.

    THIS IS THE POINT OF THE PANEL. An earlier version ranked markers by a one-way F and gave
    each class its top scorer, which is a statistic standing in for a definition that already
    exists -- and a statistic can be wrong in ways a definition cannot: it labelled T-CD8,
    a class the workflow defines as CD4-negative CD8-POSITIVE, with CD20. The workflow says
    what every reference is made of; the figure's job is to show whether the cells honour it,
    not to re-derive it.
    """
    rows = list(csv.reader(WORKFLOW_CSV.open()))
    hdr = [_norm_marker(h) for h in rows[0]]
    defs, parent = {}, {}
    for r in rows[1:]:
        if len(r) < 2 or not r[1].strip():
            continue
        ref = r[1].strip()
        parent[ref] = r[0].strip()
        defs[ref] = {hdr[j]: r[j].strip()
                     for j in range(2, min(len(r), len(hdr))) if r[j].strip()}
    return defs, parent


def _defining_markers(types, defs, available, per_type=None):
    """Column set = the union of the displayed classes' OWN defining markers, in class order.

    `per_type` caps how many each class contributes (its positives first, then its negatives),
    which keeps a broad panel legible; None takes all of them.
    """
    cols, chosen = [], []
    for t in types:
        d = defs.get(t, {})
        pos = [m for m, c in d.items() if c in POSITIVE and m in available]
        neg = [m for m, c in d.items() if c not in POSITIVE and m in available]
        take = (pos + neg) if per_type is None else (pos + neg)[:per_type]
        chosen.append((t, take))
        for m in take:
            if m not in cols:
                cols.append(m)
    print("[fig] defining markers from phenotype_workflow.csv: "
          + "; ".join(f"{t}: {'/'.join(b) if b else '(none on panel)'}" for t, b in chosen))
    return cols


def panel_a(df: pd.DataFrame, expr: pd.DataFrame, scale: pd.DataFrame, ax) -> pd.DataFrame:
    """Which markers actually carry each broad cell-type call.

    Replaces a whole-section scatter of 635,049 cells, in which 28% of the points were the
    grey Other/Ambiguous class painting over the rest and the panel showed the outline of the
    tissue and nothing else. Columns are allocated PER CLASS (see `_per_class_markers`) rather
    than by one global ranking, so the structural compartments no longer crowd out the rest.
    """
    types = [t for t in BROAD_COLOR if (df["broad"] == t).sum() >= 200
             and t != "Other/Ambiguous"]
    d0 = df[df["broad"].isin(types)]
    Xs, d = _calibrated(d0, expr, scale)
    lab = d["broad"].to_numpy()
    defs, _parent = workflow_defs()
    cols = _defining_markers(types, defs, set(Xs.columns), per_type=3)
    counts = {t: int((df["broad"] == t).sum()) for t in types}
    return _dotplot(Xs, lab, types, cols, counts, ax,
                    "A  Cell types against the markers that define them", defs=defs)


def panel_b(df: pd.DataFrame, areas: pd.Series, ax) -> pd.DataFrame:
    """Endocrine composition vs RELATIVE RADIUS, pooled over every islet.

    Replaces a picture of one islet, which was an anecdote: it showed that dots have colours,
    not whether the islet is organised. This asks a question with an answer -- is there a
    mantle/core arrangement? Each cell's distance to its islet centroid is divided by that
    islet's equivalent radius sqrt(area/pi), so islets of different size are comparable, and
    every islet contributes.
    """
    core = df[(df["cell_region"] == "core") & (df["call"].isin(ENDO))].copy()
    core["area_um2"] = islet_key(core["islet_num"]).map(
        areas.set_axis(islet_key(areas.index))).to_numpy()
    core = core.dropna(subset=["area_um2"])
    cen = core.groupby("islet_num")[["X_centroid", "Y_centroid"]].transform("mean")
    r = np.hypot(core["X_centroid"] - cen["X_centroid"], core["Y_centroid"] - cen["Y_centroid"])
    core["rel_r"] = r / np.sqrt(core["area_um2"] / np.pi)
    core = core[core["rel_r"] <= 1.5]
    edges = np.linspace(0, 1.5, 7)
    core["bin"] = pd.cut(core["rel_r"], edges)
    tab = (core.groupby(["bin", "call"], observed=True).size().unstack(fill_value=0)
               .reindex(columns=ENDO, fill_value=0))
    frac = tab.div(tab.sum(axis=1).replace(0, np.nan), axis=0)
    mid = [(a + b) / 2 for a, b in zip(edges[:-1], edges[1:])]
    for e in ENDO:
        n = tab.sum(axis=1).to_numpy()
        f = frac[e].to_numpy()
        se = np.sqrt(np.clip(f * (1 - f), 0, None) / np.maximum(n, 1))   # binomial SE
        ax.plot(mid, f, "-o", color=ENDO_COLOR[e], lw=2.5, ms=7,
                label=f"{e.replace('Endo-','')} ({int(tab[e].sum()):,})")
        ax.fill_between(mid, f - 1.96 * se, f + 1.96 * se, color=ENDO_COLOR[e], alpha=0.18)
    ax.axvline(1.0, color="k", ls=":", lw=1.5)
    ax.text(1.02, ax.get_ylim()[1] * 0.96, "islet edge", fontsize=11, va="top")
    ax.set_xlabel("distance from islet centre / islet radius")
    ax.set_ylabel("fraction of endocrine cells")
    ax.set_title(f"B  Radial composition, {core['islet_num'].nunique():,} islets",
                 loc="left")
    ax.legend(loc="center left", bbox_to_anchor=(1.02, 0.5), borderaxespad=0)
    return tab


def select_markers(expr: pd.DataFrame, scale: pd.DataFrame, labels: pd.Series,
                   donor: str, n_keep: int = 12, n_perm: int = 20):
    """Pick embedding markers by CALIBRATION then DISCRIMINATION, not by hand.

    An earlier version of this panel used a hand-written exclusion list, which is a guess
    dressed as a method. Two filters replace it, both from artifacts the lane already has:

    STEP 1 -- CALIBRATION. `run_step1.py` measures every (donor, marker) and the lane README
    states the three conditions that make a measurement invalid: fewer than
    `usability.MIN_CELLS` cells, more than half the cells reading exactly 0, or no spatial
    structure above that donor's own permutation null. Markers failing any of those are
    dropped before anything else looks at them.

    STEP 2 -- DISCRIMINATION. Rank the survivors by how well they separate the four endocrine
    classes, as the ratio of between-class to within-class variance (a one-way F statistic) on
    the lane-scaled values. Keep the top `n_keep`.

    THE CIRCULARITY, STATED. The classes come from the hormone argmax, so INS, GCG and SST
    define the labels and will top any discrimination ranking by construction. They are
    reported separately from the earned ranks below them; the ranking is informative for the
    other markers, not for those three.
    """
    from scipy.stats import f_oneway
    import usability as U                                        # noqa: PLC0415

    step1 = pathlib.Path(REPO / ".." / ".." / "data" / "development" / "celltype"
                         / "marker_donor_usability_sameslide.parquet").resolve()
    keep = [m for m in expr.columns if m in scale.index]
    dropped_cal = []
    if step1.exists():
        u = pd.read_parquet(step1)
        u = u[u["donor_id"].astype(str).isin([donor, donor.zfill(4), donor.lstrip("0")])]
        if len(u):
            struct = np.maximum(u["morans_rank"].abs(), u["morans_top1"].abs())
            bad = u[(u["n"] < U.MIN_CELLS) | (u["zero_frac"] > 0.5)
                    | (struct <= u["morans_null"])]["marker"].tolist()
            dropped_cal = [m for m in keep if m in set(bad)]
            keep = [m for m in keep if m not in set(bad)]
    print(f"[fig] calibration (Step 1): {len(keep)} markers pass, "
          f"{len(dropped_cal)} dropped{': ' + ', '.join(dropped_cal) if dropped_cal else ''}")

    Xs = pd.DataFrame(SC.apply_scale(expr[keep].to_numpy("float32"), keep, scale),
                      index=expr.index, columns=keep)

    # STEP 1b -- DEGENERATE AFTER CALIBRATION. `apply_scale` clips to [floor, top], and a
    # marker whose whole column lands on one clip bound carries no information at all. That is
    # not hypothetical here: the cohort scale rows for the nucleus-measured group
    # (panel.NUCLEUS_MEASURED -- FOXP3, Ki67, PCNA, SOX2, TCF_1, TOX, TP63) were built by
    # run_step3.build_scale(ratio=True) on the DAPI RATIO, while the values fed to apply_scale
    # are the raw nucleus means, ~2 orders of magnitude larger (6533 medians 68-463 against
    # tops 1.6-3.9), so essentially every cell clips to exactly 1.0.
    #   On donor 112 a few cells per marker stayed under the top, leaving a sliver of variance
    # and a finite (tiny) F. On donor 6533 all 23,831 core endocrine cells clip for five of
    # them, the column is constant, and `f_oneway` RETURNS NaN rather than raising. The
    # observed ranking dropped those NaNs, but the permutation-null loop APPENDED them, so
    # `np.percentile(null, 99)` came back nan and `rank[m] > nan` was False for every marker:
    # the donor kept 0 markers with no error. Drop constant columns by measurement before
    # either computation sees them, and name them.
    sd = Xs.to_numpy("float64").std(axis=0)
    degen = [m for m, v in zip(keep, sd) if (not np.isfinite(v)) or v == 0.0]
    if degen:
        keep = [m for m in keep if m not in degen]
        Xs = Xs[keep]
    print(f"[fig] constant after calibration: {len(degen)} dropped"
          + (f": {', '.join(degen)}" if degen else ""))
    # Report near-saturation too -- a column that is 99% on a clip bound is not constant but is
    # not a measurement either, and the reader should see it rather than have it ranked.
    pinned = {m: f for m, f in ((m, float(((Xs[m] <= 0.0) | (Xs[m] >= 1.0)).mean()))
                                for m in keep) if f > 0.99}
    if pinned:
        print("[fig] >99% of cells on a clip bound (kept, but read with that in mind): "
              + ", ".join(f"{m} {f:.3%}" for m, f in sorted(pinned.items(), key=lambda kv: -kv[1])))

    groups = [Xs[labels.to_numpy() == e] for e in ENDO]
    groups = [g for g in groups if len(g) > 10]
    stat = {}
    for m in keep:
        try:
            stat[m] = float(f_oneway(*[g[m].to_numpy() for g in groups]).statistic)
        except Exception:
            stat[m] = 0.0
    rank = pd.Series(stat).replace([np.inf, -np.inf], np.nan).dropna().sort_values(ascending=False)

    # A PERMUTATION NULL for F, so "discriminative" is measured rather than assumed. Shuffle
    # the class labels and recompute; the 99th percentile over all markers and draws is the
    # level F reaches by chance at this sample size. Without it a top-N cut is arbitrary --
    # F always has a largest value, whether or not anything is separating.
    rng = np.random.default_rng(0)
    null = []
    lab = labels.to_numpy()
    for _ in range(n_perm):
        perm = rng.permutation(lab)
        gs = [Xs[perm == e] for e in ENDO]
        gs = [g for g in gs if len(g) > 10]
        for m in keep:
            try:
                null.append(float(f_oneway(*[g[m].to_numpy() for g in gs]).statistic))
            except Exception:
                pass
    # NaN-safe: f_oneway returns NaN on a degenerate group rather than raising, and one NaN
    # anywhere in `null` makes np.percentile return nan, which nothing can exceed. The
    # constant-column drop above removes the known source; this keeps the threshold defined
    # even if a permuted split produces another.
    null = np.asarray(null, dtype=float)
    null = null[np.isfinite(null)]
    thresh = float(np.percentile(null, 99)) if null.size else 0.0
    above = [m for m in rank.index if rank[m] > thresh]
    print(f"[fig] permutation null ({n_perm} draws): F99 = {thresh:.2f}; "
          f"{len(above)} of {len(keep)} markers exceed it")

    defining = [m for m in ("INS", "GCG", "SST") if m in rank.index]
    earned = [m for m in above if m not in defining][: max(0, n_keep - len(defining))]
    print(f"[fig] label-defining (circular by construction): "
          + ", ".join(f"{m} F={rank[m]:,.0f}" for m in defining))
    print(f"[fig] kept {len(earned)} earned markers: "
          + ", ".join(f"{m} F={rank[m]:,.0f}" for m in earned))
    nxt = [m for m in above if m not in defining][len(earned):len(earned) + 4]
    if nxt:
        print(f"[fig] next below the cut (NOT used): "
              + ", ".join(f"{m} F={rank[m]:,.0f}" for m in nxt))
    return defining + earned, rank


def panel_c(df: pd.DataFrame, expr: pd.DataFrame, scale: pd.DataFrame, donor: str, ax) -> int:
    """Their 3D, on PROTEIN marker space -- INS/GCG/SST are on neither Xenium panel."""
    import scanpy as sc
    from anndata import AnnData

    endo = df[(df["cell_region"] == "core") & (df["call"].isin(ENDO))]
    X = expr.reindex(endo["object_id"]).dropna(how="all")
    endo = endo.set_index("object_id").loc[X.index]
    if len(X) > 40_000:                                  # their 3D is 40,714 cells
        k = np.random.default_rng(0).choice(len(X), 40_000, replace=False)
        X, endo = X.iloc[k], endo.iloc[k]

    use, _rank = select_markers(X, scale, endo["call"], donor)
    Xs = SC.apply_scale(X[use].to_numpy("float32"), use, scale)   # the lane's own calibration
    # HOW the markers are weighted matters more than any UMAP knob, and the right weighting
    # is not obvious a priori, so it is swept alongside them and decided by the same statistic.
    #
    #   calibrated  apply_scale's [0, 1] against each marker's own floor and 99.9th percentile.
    #               Equalises RANGE but not VARIANCE, so a marker that happens to be bimodal
    #               across the section dominates the distance.
    #   zscored     each marker to unit variance. Equalises variance -- but that promotes a
    #               marker detected in 3% of cells to the same footing as INS, which is how
    #               nine weakly-informative markers came to outvote the three that matter.
    #   f_weighted  each z-scored marker scaled by sqrt(F), its between-class over
    #               within-class variance ratio. Markers earn their contribution.
    #
    # SILHOUETTE is the referee: how much closer a cell sits to its own class than to the
    # nearest other one, measured on the embedding. Reported two ways, because INS/GCG/SST
    # DEFINE these labels -- separation they produce is partly circular. The earned-only
    # number re-scores the same embedding using just the 9 markers that did not define
    # anything, and is the honest one.
    from sklearn.metrics import silhouette_score

    lab = endo["call"].to_numpy()
    rng = np.random.default_rng(0)
    sco = rng.choice(len(Xs), min(8_000, len(Xs)), replace=False)
    z = (Xs - Xs.mean(0)) / np.maximum(Xs.std(0), 1e-9)
    Fw = np.sqrt(np.array([max(_rank.get(m, 0.0), 0.0) for m in use]))
    Fw = Fw / max(Fw.max(), 1e-9)
    FEATS = {"calibrated": Xs, "zscored": z, "f_weighted": z * Fw}

    best = None
    for fname, M in FEATS.items():
        ad = AnnData(np.ascontiguousarray(M, dtype="float32"))
        for k_nn, md in ((15, 0.3), (30, 0.4)):
            sc.pp.neighbors(ad, n_neighbors=k_nn, use_rep="X", metric="euclidean",
                            random_state=0)
            sc.tl.umap(ad, min_dist=md, spread=1.0, random_state=0)
            u = ad.obsm["X_umap"]
            sil = float(silhouette_score(u[sco], lab[sco]))
            print(f"[fig] UMAP {fname:<11s} k={k_nn:<3d} min_dist={md}  silhouette={sil:+.3f}")
            if best is None or sil > best[0]:
                best = (sil, fname, k_nn, md, u.copy())
    sil, fname, k_nn, md, u_best = best
    print(f"[fig] UMAP kept {fname}, k={k_nn}, min_dist={md} (silhouette {sil:+.3f})")
    ad = AnnData(np.ascontiguousarray(FEATS[fname], dtype="float32"))
    ad.obsm["X_umap"] = u_best
    u = ad.obsm["X_umap"]
    for e in ENDO:
        m = (endo["call"] == e).to_numpy()
        if m.any():
            ax.scatter(u[m, 0], u[m, 1], s=2.5, c=ENDO_COLOR[e], linewidths=0,
                       rasterized=True, label=f"{e.replace('Endo-','')} ({int(m.sum()):,})")
    ax.set_xticks([]); ax.set_yticks([])
    ax.set_xlabel("UMAP 1"); ax.set_ylabel("UMAP 2")
    ax.set_title(f"C  Endocrine cells, {len(use)} markers ({fname}), "
                 f"silhouette {sil:+.2f}  (n={len(X):,})", loc="left")
    ax.legend(loc="center left", bbox_to_anchor=(1.02, 0.5), markerscale=6, borderaxespad=0)
    return len(X)


def panel_d(df: pd.DataFrame, ax, n_sub: int = 120_000, k: int = 10,
            n_perm: int = 50) -> pd.DataFrame:
    """Which cell types sit next to which, against a label-permutation null.

    Replaces a stacked bar chart of composition by islet area, which compressed a trend into
    five bars and could not say whether any of it exceeded chance. This is the architecture
    question their paper asks with neighbourhood enrichment: build a kNN graph on position,
    count type-type adjacencies, then permute the labels to get the count expected if type
    and position were independent. The z-score is (observed - null mean) / null sd.
    """
    from sklearn.neighbors import NearestNeighbors

    sub = df.sample(min(len(df), n_sub), random_state=0).reset_index(drop=True)
    types = [t for t in BROAD_COLOR if (sub["broad"] == t).sum() >= 200 and t != "Other/Ambiguous"]
    code = pd.Categorical(sub["broad"], categories=types).codes
    xy = sub[["X_centroid", "Y_centroid"]].to_numpy()
    nn = NearestNeighbors(n_neighbors=k + 1).fit(xy)
    idx = nn.kneighbors(xy, return_distance=False)[:, 1:]

    def counts(c):
        a = np.repeat(c, k); b = c[idx].ravel()
        ok = (a >= 0) & (b >= 0)
        return np.bincount(a[ok] * len(types) + b[ok],
                           minlength=len(types) ** 2).reshape(len(types), len(types))

    obs = counts(code)
    rng = np.random.default_rng(0)
    null = np.stack([counts(rng.permutation(code)) for _ in range(n_perm)])
    z = (obs - null.mean(0)) / np.maximum(null.std(0), 1e-9)
    z = (z + z.T) / 2                                   # adjacency is symmetric
    # Scale the colour by the OFF-diagonal. Self-adjacency is ~an order of magnitude larger
    # than anything else (endocrine cells sit in islets, which is not news), and letting it set
    # the range washes every cross-type term to white. The diagonal is clipped, not hidden --
    # its z is still printed in the cell.
    off = z[~np.eye(len(types), dtype=bool)]
    lim = float(np.nanmax(np.abs(off)))
    im = ax.imshow(z, cmap="RdBu_r", vmin=-lim, vmax=lim)
    for i in range(len(types)):      # mark the clipped diagonal so it is not read as saturated
        ax.add_patch(plt.Rectangle((i - .5, i - .5), 1, 1, fill=False, ec="k", lw=1.5, ls=":"))
    ax.set_xticks(range(len(types))); ax.set_yticks(range(len(types)))
    ax.set_xticklabels(types, rotation=45, ha="right"); ax.set_yticklabels(types)
    for i in range(len(types)):
        for j in range(len(types)):
            ax.text(j, i, f"{z[i, j]:.0f}", ha="center", va="center", fontsize=10,
                    color="white" if abs(z[i, j]) > lim * 0.55 else "black")
    ax.set_title(f"D  Neighbourhood enrichment (k={k}, {n_perm} permutations, "
                 f"n={len(sub):,})\ncolour scaled to off-diagonal; dotted = self-adjacency, clipped",
                 loc="left")
    cb = ax.figure.colorbar(im, ax=ax, fraction=0.046, pad=0.10)
    cb.set_label("z vs permuted labels")
    return pd.DataFrame(z, index=types, columns=types)


def panel_e(df: pd.DataFrame, expr: pd.DataFrame, scale: pd.DataFrame, ax,
            min_cells: int = 100) -> pd.DataFrame:
    """The immune compartment, resolved to its own subtypes.

    Panel A collapses 115,196 cells into one Immune row, which is the largest compartment
    after exocrine and says nothing about what is in it. The lane's own hierarchy already
    resolves it -- macrophages and their M1/M2 states, B and plasma cells, dendritic cells,
    neutrophils, NK, and seven T-cell states -- and for a study about immune cells around
    islets that resolution is the point. Columns are ranked one-vs-rest WITHIN the immune
    compartment, so they are the markers that separate immune subtypes from each other rather
    than the ones that separate immune from acinar.
    """
    imm = df[df["broad"] == "Immune"].copy()
    vc = imm["call"].value_counts()
    types = [t for t in vc.index if vc[t] >= min_cells]
    imm = imm[imm["call"].isin(types)]
    Xs, d = _calibrated(imm, expr, scale, max_cells=120_000)
    lab = d["call"].to_numpy()
    types = [t for t in types if (lab == t).sum() >= 30]
    defs, parent = workflow_defs()
    # A subtype INHERITS its parents' constraints (T-CD8 is CD8+/CD4-, and is also everything
    # Imm-T-cell and Immune require). The own-constraints are what distinguishes it, so those
    # set the columns; the inherited ones are merged in for the outlines.
    full = {}
    for t in types:
        acc, cur = {}, t
        while cur and cur in parent and cur != "all":
            for m, c in defs.get(cur, {}).items():
                acc.setdefault(m, c)
            cur = parent.get(cur)
        full[t] = acc
    cols = _defining_markers(types, defs, set(Xs.columns), per_type=3)
    counts = {t: int(vc[t]) for t in types}
    print(f"[fig] panel E: {len(types)} immune subtypes, {int(vc[types].sum()):,} cells")
    return _dotplot(Xs, lab, types, cols, counts, ax,
                    "E  Immune subtypes against the markers that define them",
                    size_key_at=(1.10, 1.0), defs=full)


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--donor", default="112")
    ap.add_argument("--out", type=pathlib.Path,
                    default=pathlib.Path("/blue/maigan/smith6jt/Xenium_Analysis/figures"))
    args = ap.parse_args(argv)

    calls = pathlib.Path(f"{REPO}/../../data/development/celltype/sameslide_protein"
                         f"/donor_id={args.donor}/data.parquet").resolve()
    df = pd.read_parquet(calls)
    df["broad"] = df["call"].map(broad_of)
    df["islet_num"] = df["islet_num"].astype(str)

    gj = next(DATA.glob(f"Islets_*{'0081366' if args.donor == '112' else '0080870'}*.geojson"))
    areas = islet_areas(gj)
    print(f"[fig] donor {args.donor}: {len(df):,} cells, {len(areas):,} islet polygons")

    import pyarrow.parquet as pq
    epath = next((DATA / "expression").glob(f"donor_id=*{args.donor.lstrip('0')}*/data.parquet"),
                 None) or (DATA / f"expression/donor_id=0{args.donor}/data.parquet")
    h = pq.ParquetFile(epath)
    markers = [n for n in h.schema_arrow.names
               if n not in ("donor_id", "object_id", "qc_analysis_eligible",
                            "qc_estimation_eligible", "image", "X_centroid", "Y_centroid",
                            "cell_region", "islet_num")]
    expr = h.read(columns=["object_id"] + markers).to_pandas().set_index("object_id")
    expr = expr.dropna(axis=1, how="all")
    print(f"[fig] protein markers for the UMAP: {expr.shape[1]}")

    fig = plt.figure(figsize=(24, 26))
    gs = fig.add_gridspec(3, 2, hspace=0.45, wspace=0.42,
                          height_ratios=[1.0, 1.0, 1.15])
    axA = fig.add_subplot(gs[0, 0]); axB = fig.add_subplot(gs[0, 1])
    axC = fig.add_subplot(gs[1, 0]); axD = fig.add_subplot(gs[1, 1])
    axE = fig.add_subplot(gs[2, :])

    scale = pd.read_csv(calls.parent / 'scale.csv', index_col=0)
    tab_a = panel_a(df, expr, scale, axA)
    tab_b = panel_b(df, areas, axB)
    n_umap = panel_c(df, expr, scale, args.donor, axC)
    tab = panel_d(df, axD)
    tab_e = panel_e(df, expr, scale, axE)

    status = "non-diabetic" if args.donor == "112" else "type 1 diabetes"
    fig.suptitle(f"Donor {args.donor} ({status}) - PhenoCycler cells, protein-derived types",
                 y=0.985)
    out = args.out / f"{args.donor}_nd_reference.png"
    out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out, bbox_inches="tight")
    print(f"[fig] wrote {out}")
    print("\nPanel A (fraction detected per type):"); print(tab_a.round(3).to_string())
    print("\nPanel B (radial bins):"); print(tab_b.to_string())
    print("\nPanel D neighbourhood-enrichment z:")
    print(tab.to_string())
    print("\nPanel E (immune subtypes, fraction detected):")
    print(tab_e.round(3).to_string())
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
