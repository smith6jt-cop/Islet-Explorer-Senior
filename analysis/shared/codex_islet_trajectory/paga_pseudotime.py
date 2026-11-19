"""
PAGA pseudotime + force-directed layout for the Islet Explorer trajectory tab.

Ports Barlow et al., eLife 2024;13:RP100535, Methods -> "PAGA analysis":

    "The PAGA embedding was computed using the default parameters except for the
     following: The neighborhood search was performed using cosine distance and
     15 nearest neighbors; Leiden clustering used a resolution of 1."

and their Figure 2B, "PAGA-force directed layout of islets colored by pseudotime.
Each point represents an islet."

The force-directed layout is the point of this module. A feature-vs-pseudotime
scatter piles ~85-96% of islets onto already-inked pixels because most cell-type
proportions are exactly zero (95.6% for CD8a T cells, 98.2% for B cells) and no
axis transform moves an exact zero. The PAGA layout is a 2-D graph embedding, so
the layout algorithm spreads islets apart by construction: measured on the same
occupancy grid, 8.1% overplotted vs 84-96% for the scatter, densest pixel 3
islets vs 125-129.

The graph inputs are log1p(cells / endocrine cells) per class and region in the
current typing format (owner, 2026-11-18; "Graph inputs" below): the endocrine
classes in the core, exocrine in the peri-islet band, macrophages by Gal3, and
no QC class.

Called from R via reticulate (app/shiny_app/R/paga_layout.R). Kept as a real
module rather than inline R strings so it can be run and tested standalone:

    python app/shiny_app/python/paga_pseudotime.py --min-cells 10
"""

import os
import re
from collections import OrderedDict
import numpy as np
import pandas as pd

# Barlow's stated parameters.
N_NEIGHBORS = 15
METRIC = "cosine"
LEIDEN_RESOLUTION = 1.0
N_DIFFMAP_COMPS = 10

# ---------------------------------------------------------------------------
# Graph inputs: the current typing format (owner, 2026-11-18)
# ---------------------------------------------------------------------------
# "The trajectory/pseudotime needs to be updated to the current typing format.
#  There should not be peri-islet features for the endocrine classes, core
#  exocrine, or 12 macrophage subtypes. Also the qc classes should not be
#  entered into the pseudotime pipeline."
#
#   endocrine classes   core only
#   exocrine classes    peri-islet band only
#   macrophages         the two Gal3 subclasses; the twelve marker-profile labels
#                       of an older root are summed by their Gal3 sign
#   QC classes          never: Conflict, Doublet suspected, Tissue, Unclassified
#                       and the unresolved lineage residuals
#   every other named class   core and peri-islet band
#
# R/barlow_pseudotime.R applies the same rules for the destiny fallback; the
# label sets below are compared with the R source by a regression test.
TRAJECTORY_INPUT_POLICY = "phenotype-passthrough/2026-11-19"

# Counted in the core only (uniform_release.json label_groups.endocrine).
ENDOCRINE_CLASSES = frozenset({
    "Beta cell", "Alpha cell", "Delta cell", "CHGA+ hormone-negative",
})
# Counted in the peri-islet band only: the uniform roots' Exocrine, its Acinar /
# Ductal split in the November revision, and the frozen typing's Ductal.
EXOCRINE_CLASSES = frozenset({"Exocrine", "Acinar cell", "Ductal cell", "Ductal"})
# names(UNCLASSIFIED_RAW_TO_KEY) in R/00_globals.R: the app's Unclassified family.
QC_CLASSES = frozenset({
    "Tissue", "Unclassified", "Doublet suspected",
    "Hematopoietic CD45+", "Hematopoietic CD45-", "T cell",
})
# PHENOTYPE_RETIRED in R/00_globals.R. These frozen-typing labels are not
# comparable trajectory features across the original and added donors.
RETIRED_PHENOTYPES = frozenset({
    "Endocrine", "Acinar", "Endothelial", "Immune", "Unknown",
    "ECAD+", "Structural", "APCs",
})
MACROPHAGE_LABELS = ("Macrophage_LGALS3+", "Macrophage_LGALS3-")
# MACROPHAGE_LABEL_PATTERN in R/00_globals.R: the two current labels and the
# three- or four-marker labels of immutable older releases.
MACROPHAGE_LABEL_PATTERN = (r"^(Macrophage_Gal3[+-]|(CD68[+]CD163[+-]|CD68-CD163[+])?"
                            r"Macrophage_LGALS3[+-](CD4[+-])?)$")
_CONFLICT_PATTERN = re.compile(r"^Conflict($|[\s_])", re.IGNORECASE)

_FEATURES_CSV = None
_FEATURES = None
_FIT_CACHE = OrderedDict()
_FIT_CACHE_LIMIT = 2


def trajectory_class(label):
    """The current-format class of a data label: a macrophage label becomes its Gal3 subclass."""
    label = str(label)
    if re.match(MACROPHAGE_LABEL_PATTERN, label):
        return "Macrophage_LGALS3" + re.search(r"(?:LGALS3|Gal3)([+-])", label).group(1)
    return label


def is_qc_class(label):
    """Conflict, Doublet suspected, Tissue, Unclassified and the unresolved lineage residuals."""
    label = str(label).strip()
    return label in QC_CLASSES or bool(_CONFLICT_PATTERN.match(label))


def trajectory_feature_allowed(region, label):
    """Whether a class counted in a region ("core" or "peri") is a graph input."""
    cls = trajectory_class(label)
    if is_qc_class(cls) or cls in RETIRED_PHENOTYPES:
        return False
    if cls in ENDOCRINE_CLASSES:
        return region == "core"
    if cls in EXOCRINE_CLASSES:
        return region == "peri"
    return True


def composition_features(columns):
    """The stored core_/peri_ ratio columns that feed the graph, in file order.

    Twelve-label macrophage columns are returned as stored; trajectory_inputs()
    sums them by Gal3 sign.
    """
    return [c for c in columns if c.startswith(("core_", "peri_"))
            and trajectory_feature_allowed(*c.split("_", 1))]


def phenotype_inputs(frame):
    """The graph-input names and matrix for a feature table.

    A class with one stored column is used as stored, so its values and the
    stored endocrine denominator are unchanged. A class with several source
    labels (the pooled macrophages) is log1p(summed count / endocrine cells) from
    the n_core_ / n_peri_ counts: the log of a sum is not the sum of the logs.
    """
    members = OrderedDict()
    for column in composition_features(frame.columns):
        region, label = column.split("_", 1)
        cls = trajectory_class(label)
        if cls not in ENDOCRINE_CLASSES and cls not in EXOCRINE_CLASSES:
            region = "combined"
        members.setdefault(f"{region}_{cls}", []).append(column)
    names, values = [], []
    for name, columns in members.items():
        if columns == [name]:
            value = frame[name].to_numpy(dtype=float)
        else:
            if name.startswith("combined_"):
                labels = list(dict.fromkeys(c.split("_", 1)[1] for c in columns))
                columns = [f"{region}_{label}" for label in labels for region in ("core", "peri")]
            counts = ["n_" + c for c in columns]
            missing = [c for c in counts + ["n_endocrine_core"] if c not in frame.columns]
            if missing:
                raise ValueError(f"{name} pools {len(columns)} source labels and needs their "
                                 f"counts; missing {missing}")
            total = frame[counts].to_numpy(dtype=float).sum(axis=1)
            endocrine = np.maximum(frame["n_endocrine_core"].to_numpy(dtype=float), 1.0)
            value = np.log1p(total / endocrine)
        names.append(name)
        values.append(value)
    if not names:
        return names, np.empty((len(frame), 0))
    return names, np.column_stack(values)


OVERLAP_STRUCTURES = ("nerve", "lymphatic", "vessel")


def trajectory_inputs(frame):
    names, counts = phenotype_inputs(frame)
    overlap_names = [f"{name}__{structure}" for name in names for structure in OVERLAP_STRUCTURES]
    raw_columns = ["overlap__" + name for name in overlap_names]
    missing = [name for name in raw_columns if name not in frame]
    if missing:
        raise ValueError(f"Missing phenotype-specific passthrough measurements: {missing}")
    raw = frame[raw_columns].to_numpy(dtype=float)
    if not (np.isfinite(raw) & (raw >= 0)).all():
        raise ValueError("Invalid phenotype-specific passthrough totals")
    denominator = np.maximum(frame["n_endocrine_core"].to_numpy(dtype=float), 1)
    return names + overlap_names, np.column_stack([counts, np.log1p(raw / denominator[:, None])])


def _features_path():
    """Resolve data/barlow_islet_features.csv from either the app dir or repo root."""
    global _FEATURES_CSV
    if _FEATURES_CSV and os.path.exists(_FEATURES_CSV):
        return _FEATURES_CSV
    here = os.path.dirname(os.path.abspath(__file__))
    for cand in (
        os.path.join(os.environ.get("ISLET_DATA_ROOT", ""), "barlow_islet_features.csv"),
        os.path.join(here, "..", "..", "..", "data", "barlow_islet_features.csv"),
        os.path.join(here, "..", "..", "data", "barlow_islet_features.csv"),
        os.path.join("data", "barlow_islet_features.csv"),
    ):
        if os.path.exists(cand):
            _FEATURES_CSV = os.path.abspath(cand)
            return _FEATURES_CSV
    return None


def _range_values(min_cells, max_cells):
    try:
        mc = max(1, int(float(min_cells)))
    except (TypeError, ValueError):
        mc = 1
    try:
        xc = None if max_cells is None else int(float(max_cells))
    except (TypeError, ValueError):
        xc = None
    return mc, xc


def _classes(object_classes):
    """The endocrine-object classes to fit on (Islet / EndocrineCluster); None = all."""
    if object_classes is None:
        return None
    if isinstance(object_classes, str):
        object_classes = [object_classes]
    return tuple(sorted(str(c) for c in object_classes))


def _context(min_cells, max_cells, object_classes=None):
    """Cache the selected cohort and its neighbor graph for DPT and PAGA."""
    global _FEATURES
    mc, xc = _range_values(min_cells, max_cells)
    classes = _classes(object_classes)
    key = (mc, xc, classes)
    if key in _FIT_CACHE:
        _FIT_CACHE.move_to_end(key)
        return _FIT_CACHE[key]

    path = _features_path()
    if path is None:
        return {"error": "data/barlow_islet_features.csv not found "
                         "(run scripts/build_barlow_features.py)"}
    try:
        import scanpy as sc
        import anndata as ad
    except ImportError as e:
        return {"error": f"scanpy/anndata unavailable in this interpreter: {e}"}

    if _FEATURES is None:
        _FEATURES = pd.read_csv(path)
    b = _FEATURES
    n = b["total_cells_core"]
    sel = n >= mc
    if xc is not None:
        sel &= n <= xc
    if classes is not None and "object_class" in b.columns:
        sel &= b["object_class"].astype(str).isin(classes)
    sub = b[sel].reset_index(drop=True)
    if len(sub) < 50:
        return {"error": f"only {len(sub)} islets in cell range "
                         f"[{mc}, {xc if xc is not None else 'inf'}]; "
                         "a diffusion map needs more"}

    try:
        feat, X = trajectory_inputs(sub)
    except ValueError as e:
        return {"error": str(e)}
    # An object with no cell of any input class has no make-up to place and no direction for the cosine
    # metric. Like an object outside the cell range it gets no pseudotime; the number is reported.
    placed = (X != 0).any(axis=1)
    n_without_inputs = int((~placed).sum())
    if n_without_inputs:
        sub, X = sub[placed].reset_index(drop=True), X[placed]
        if len(sub) < 50:
            return {"error": f"only {len(sub)} objects in the cell range have a cell of an input class"}
    # Constant columns carry no information and make the cosine metric unstable.
    keep = X.std(axis=0) > 0
    X = X[:, keep]
    if X.shape[1] < 3:
        return {"error": "fewer than 3 informative composition features"}

    # The graph is built on the DISTINCT composition profiles and mapped back to every object (`inv`).
    # Since 2026-10-17 the endocrine-object classes, not a cells-per-islet range, choose the cohort, and Singles
    # and Clusters share few profiles (all 9,496 objects of the uniform root have 3,773): identical points form
    # disconnected k-NN cliques, the diffusion map's eigensolver does not converge (ARPACK, 3.4 min), and the
    # app fell back to destiny, which overflowed R's C stack. Identical profiles get identical pseudotime anyway.
    U, inv = np.unique(np.round(X, 9), axis=0, return_inverse=True)
    inv = np.asarray(inv).reshape(-1)
    status = sub["donor_status"].astype(str).to_numpy()
    has_nd = np.zeros(len(U), bool)
    has_nd[inv[status == "ND"]] = True
    first = np.full(len(U), -1)
    for i in range(len(inv) - 1, -1, -1):
        first[inv[i]] = i
    obs = pd.DataFrame({"has_nd": has_nd,
                        "case_id": sub["case_id"].astype(str).to_numpy()[first],
                        "n_objects": np.bincount(inv, minlength=len(U))},
                       index=pd.Index([str(i) for i in range(len(U))], name="profile_id"))
    A = ad.AnnData(X=U, obs=obs)
    sc.pp.neighbors(A, n_neighbors=min(N_NEIGHBORS, len(U) - 1), metric=METRIC, use_rep="X")
    ctx = {"sc": sc, "A": A, "sub": sub, "inv": inv, "mc": mc, "xc": xc,
           "features": [c for c, used in zip(feat, keep) if used],
           "n_without_inputs": n_without_inputs,
           "root": None, "pt": None, "layout": None}
    _FIT_CACHE[key] = ctx
    while len(_FIT_CACHE) > _FIT_CACHE_LIMIT:
        _FIT_CACHE.popitem(last=False)
    return ctx


def _ensure_dpt(ctx):
    if ctx.get("pt") is not None:
        return ctx["pt"]
    sc, A, sub = ctx["sc"], ctx["A"], ctx["sub"]
    sc.tl.diffmap(A, n_comps=min(N_DIFFMAP_COMPS, A.n_obs - 2))
    dm = np.asarray(A.obsm["X_diffmap"])[:, 1:4]
    nd = np.where(A.obs["has_nd"].to_numpy())[0]
    if nd.size == 0:
        return {"error": "no ND objects in this selection; no root available"}
    cen = dm[nd].mean(axis=0)
    root = int(nd[np.argmin(np.linalg.norm(dm[nd] - cen, axis=1))])
    A.uns["iroot"] = root
    sc.tl.dpt(A)
    pt = A.obs["dpt_pseudotime"].to_numpy(dtype=float)[ctx["inv"]]
    beta = sub["core_Beta cell"].to_numpy(dtype=float)
    fin = np.isfinite(pt) & np.isfinite(beta)
    flipped = False
    # Keep the existing beta-ratio orientation convention. DPT is a distance
    # from the root, so an eigenvector sign change cannot reverse it; this
    # convention also does not guarantee a monotonic insulin relationship.
    if fin.sum() > 10:
        from scipy.stats import spearmanr
        if spearmanr(pt[fin], beta[fin]).statistic > 0:
            pt = np.nanmax(pt[np.isfinite(pt)]) - pt
            flipped = True
    ctx["root"], ctx["pt"], ctx["flipped"] = root, pt, flipped
    return pt


def compute_pseudotime(min_cells=1, max_cells=None, object_classes=None):
    """Return only the DPT axis; leave clustering and force layout deferred."""
    ctx = _context(min_cells, max_cells, object_classes)
    if "error" in ctx:
        return ctx
    pt = _ensure_dpt(ctx)
    if isinstance(pt, dict):
        return pt
    sub = ctx["sub"]
    return {"islet_id": sub["islet_id"].astype(str).tolist(),
            "pseudotime": pt.astype(float).tolist(),
            "features": list(ctx["features"]),
            "input_policy": TRAJECTORY_INPUT_POLICY,
            "n_without_inputs": int(ctx["n_without_inputs"])}


def compute(min_cells=1, max_cells=None, object_classes=None):
    """Return DPT plus the PAGA-seeded force-directed layout and diagnostics."""
    ctx = _context(min_cells, max_cells, object_classes)
    if "error" in ctx:
        return ctx
    sc, A, sub = ctx["sc"], ctx["A"], ctx["sub"]
    mc, xc = ctx["mc"], ctx["xc"]

    if ctx["layout"] is None:
        # flavor="igraph" is scanpy's supported non-deprecated path and is what the
        # future default will be; n_iterations=2 matches that default.
        sc.tl.leiden(A, resolution=LEIDEN_RESOLUTION, key_added="leiden",
                     flavor="igraph", n_iterations=2, directed=False, random_state=0)
        sc.tl.paga(A, groups="leiden")
        # PAGA connectivities must be computed before they can seed the layout.
        sc.pl.paga(A, plot=False)
        sc.tl.draw_graph(A, init_pos="paga", random_state=0)
        lay_key = next(k for k in A.obsm if k.startswith("X_draw_graph"))
        xy = np.asarray(A.obsm[lay_key])[ctx["inv"]]
        ctx["layout"] = (lay_key, xy)
    else:
        lay_key, xy = ctx["layout"]

    pt = _ensure_dpt(ctx)
    if isinstance(pt, dict):
        return pt
    root = ctx["root"]
    flipped = ctx.get("flipped", False)

    return {
        "islet_id": sub["islet_id"].astype(str).tolist(),
        "case_id": sub["case_id"].astype(str).tolist(),
        "donor_status": sub["donor_status"].astype(str).tolist(),
        "leiden": A.obs["leiden"].astype(str).to_numpy()[ctx["inv"]].tolist(),
        "fa_1": xy[:, 0].astype(float).tolist(),
        "fa_2": xy[:, 1].astype(float).tolist(),
        "pseudotime": pt.astype(float).tolist(),
        "n_islets": int(len(sub)),
        "n_profiles": int(A.n_obs),
        "n_clusters": int(A.obs["leiden"].nunique()),
        "layout_key": str(lay_key),
        "min_cells": mc,
        "max_cells": xc,
        "root_index": root,
        "root_donor": str(A.obs["case_id"].iloc[root]),
        "flipped": bool(flipped),
        "features": list(ctx["features"]),
        "input_policy": TRAJECTORY_INPUT_POLICY,
        "n_without_inputs": int(ctx["n_without_inputs"]),
    }


if __name__ == "__main__":
    import argparse, json
    ap = argparse.ArgumentParser()
    ap.add_argument("--min-cells", type=int, default=10)
    ap.add_argument("--max-cells", type=int, default=None)
    a = ap.parse_args()
    r = compute(a.min_cells, a.max_cells)
    if "error" in r:
        print("ERROR:", r["error"])
    else:
        print(json.dumps({k: v for k, v in r.items()
                          if not isinstance(v, list)}, indent=2))
