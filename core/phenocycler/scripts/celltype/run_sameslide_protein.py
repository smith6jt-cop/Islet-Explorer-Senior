"""Protein-arm cell typing for the two same-slide Pheno-Xeno donors (6533, 0112).

This is Stage 5 of the streamlined Xenium plan, and under the user-directed order it runs
**first** -- before the RNA arm, not after it.  The reason is arithmetic: these bundles carry
~126 estimated false-positive transcripts against a median 433 per cell, so the RNA needs a
context prior and the protein is the only measurement dense enough to supply one.  70-plex on
the same physical section, against ~433 sparse transcripts of which ~29% are gDNA background.

It is the stock :mod:`celltype` lane with the two cohort-wide steps removed, because these
two donors are not part of the cohort that was aligned:

* **No UniFORM.**  ``run_step3`` reads ``data/development/celltype/uniform`` -- the Step-2
  aligned partitions, which exist only for the 22 cohort donors.  These two are not in it.
* **No cohort scale.**  ``run_step3.build_scale`` pools every donor onto one floor and one
  top per marker.  Pooling two donors whose INS whole-cell means differ 10x (11.6 against
  113.0) would put one of them entirely under the other's floor.  A per-donor scale is used
  instead, and the cost is stated rather than hidden: **6533 and 0112 are not on a common
  scale, so their scaled values are not comparable cell-for-cell between donors.**  Within a
  donor they are, which is what a per-donor call needs.

Everything that decides a cell type is reused unchanged -- ``rules.flatten`` /
``reference_matrix`` / ``requirements``, ``score.score_and_gate`` / ``call_from_scores``, and
``run_islet.intruder_references`` / ``subtype_endocrine``.  In particular the islet prior uses
``run_islet.intruder_references``, which excludes the **whole epithelial branch**, ancestors
and siblings, from the set allowed to displace the prior.  A hand-rolled filter that excluded
only endocrine-prefixed names was tried and measured: the generic ``Epithelial`` parent then
took 21,450 of 65,363 islet cells, reproducing exactly the failure that function exists to
prevent.

Input is the same-slide lane's own output, which already follows the stock ``expression``
contract -- REDSEA spillover-corrected whole-cell means for the 64 corrected markers, and the
RAW QuPath *nucleus* mean for the eight in ``panel.NUCLEUS_MEASURED`` (verified column by
column against ``redsea/`` and ``cells/``).

    python scripts/celltype/run_sameslide_protein.py                 # both donors
    python scripts/celltype/run_sameslide_protein.py --donors 6533
"""
from __future__ import annotations

import argparse
import json
import os
import pathlib
import sys
import time

import numpy as np
import pandas as pd
import pyarrow.parquet as pq

HERE = pathlib.Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import panel                                     # noqa: E402
import rules as R                                # noqa: E402
import run_islet as RI                           # noqa: E402
import run_step3 as S3                           # noqa: E402
import score as SC                               # noqa: E402

#: Where the same-slide lane wrote its per-donor partitions.  Not under this repo: the raw
#: PhenoCycler and Xenium inputs for these two donors live in the Xenium_Analysis tree, and
#: the lane wrote its outputs beside them.
DATA_ROOT = pathlib.Path(os.environ.get(
    "SAMESLIDE_PHENO_DATA", "/blue/maigan/smith6jt/Xenium_Analysis/data/phenocycler"))

OUT = panel.OUT / "sameslide_protein"

DONORS = ("6533", "0112")

#: PhenoCycler qptiff pixel size, from the image metadata both donors declare.  The islet
#: GeoJSONs are exported in qptiff PIXEL space while the cell table is already in microns, so
#: one of the two has to be converted and this is the factor.
PIXEL_UM = 0.5077663810243286

#: Islet annotations, per donor, by the route each one actually reached disk by.
#:
#: 6533's QuPath export carried the islet assignment in the ``Parent`` column, so its
#: partitions already have ``cell_region`` in {core, peri, tissue} and no geometry is needed.
#: 0112's TSV carries no islet parent -- every one of its 758,451 cells reads ``tissue`` --
#: but its islet GeoJSON is intact, so membership is recomputed by polygon containment.
#: (6533's own islets GeoJSON was overwritten by a cells export and cannot be used.)
ISLET_GEOJSON = {
    "0112": "Islets_0081366_HDL112_pLNpanc_Scan3.er.qptiff - resolution #1.geojson",
}


def expression_path(donor: str) -> pathlib.Path:
    return DATA_ROOT / "expression" / f"donor_id={donor}" / "data.parquet"


def load_expression(donor: str, markers: list) -> tuple:
    """One donor's selected expression, split into the rule markers it MEASURED and the rest.

    A rule marker can be absent two ways and both mean the same thing: the column is not in
    the file, or the column is there and every value is NaN.  The second is the live case --
    this 70-plex panel carries neither ``Ker8_18`` nor ``IAPP``, and the writer emitted them
    as all-NaN placeholders, so a naive read would put a column of NaN into the cosine.

    Nothing is inferred from an unmeasured marker.  It is removed from the marker space
    entirely, and :func:`audit_unmeasured` then reports exactly which constraints that costs
    -- a marker silently left in place would enter the angle as a zero, which reads as
    "measured and absent" rather than "not measured".
    """
    path = expression_path(donor)
    handle = pq.ParquetFile(path)
    on_disk = {f.name for f in handle.schema_arrow}
    keep = ["object_id", "qc_analysis_eligible", "cell_region", "islet_num",
            "X_centroid", "Y_centroid"]
    keep = [c for c in keep if c in on_disk]
    present = [m for m in markers if m in on_disk]
    # The file, not the dataset: the partition directory carries donor_id as a dictionary
    # key while the file carries it as a string, and merging the two schemas fails.
    df = handle.read(columns=list(dict.fromkeys(present + keep))).to_pandas()
    measured = [m for m in present
                if np.isfinite(df[m].to_numpy("float64")).any()]
    unmeasured = [m for m in markers if m not in measured]
    return df, measured, unmeasured


def audit_unmeasured(donor: str, unmeasured: list, Rm: pd.DataFrame, reqs: dict) -> list:
    """Report what each unmeasured marker costs, and drop the references it kills.

    Three different consequences, and conflating them is how a panel gap turns into a wrong
    call:

    * **A positive OR-group loses one alternative.**  The epithelial gate is
      "any of E-cadherin | Ker8/18 | Pan-CK | EpCAM", so losing Ker8/18 leaves three
      measured alternatives and the gate still means what it says.
    * **A negative constraint disappears.**  Every immune reference carries ``Ker8_18`` at
      -1 in the signed encoding; without it they are less able to push an epithelial cell
      away.  A loss of specificity, reported, not fatal.
    * **A mandatory constraint or a whole OR-group disappears.**  Then the reference cannot
      be evaluated as written, and it is DROPPED rather than scored on what is left -- which
      would silently weaken the rule into something the maintainer never wrote.
    """
    if not unmeasured:
        return list(Rm.columns)
    dead, weakened_pos, weakened_neg = [], {}, {}
    gone = set(unmeasured)
    for ph in Rm.columns:
        all_of, groups = reqs.get(ph, ([], []))
        if set(all_of) & gone or any(set(g) <= gone for g in groups):
            dead.append(ph)
            continue
        lost_pos = [m for g in groups for m in g if m in gone]
        if lost_pos:
            weakened_pos[ph] = lost_pos
        lost_neg = [m for m in gone if m in Rm.index and Rm.loc[m, ph] < 0]
        if lost_neg:
            weakened_neg[ph] = lost_neg
    print(f"  {donor}: {len(unmeasured)} rule marker(s) not measured on this panel: "
          f"{unmeasured}")
    if weakened_pos:
        print(f"    OR-group alternative lost (gate still satisfiable) on "
              f"{len(weakened_pos)} reference(s): {sorted(weakened_pos)}")
    if weakened_neg:
        print(f"    negative constraint lost (less specific) on {len(weakened_neg)} "
              f"reference(s): {sorted(weakened_neg)}")
    if dead:
        print(f"    DROPPED -- a mandatory constraint or a whole OR-group is unmeasurable: "
              f"{dead}")
    return [p for p in Rm.columns if p not in dead]


def donor_scale(df: pd.DataFrame, markers: list) -> pd.DataFrame:
    """Per-donor floor and robust top, on this donor's analysis-eligible cells only.

    The floor is the mode of ``log1p(value)`` over positive values -- the same estimator
    ``run_step3.build_scale`` uses, kept identical so the two lanes mean the same thing by
    "background" -- and the top is the 99.9th percentile rather than the maximum, so one hot
    pixel cannot set the scale.  What differs from ``build_scale`` is only the pool: this
    donor, not the cohort.
    """
    eligible = df[df["qc_analysis_eligible"]] if "qc_analysis_eligible" in df else df
    pool, floors = {}, {}
    for m in markers:
        v = eligible[m].to_numpy("float64")
        v = v[np.isfinite(v)]
        pool[m] = v
        pos = v[v > 0]
        if len(pos) < 1000:
            continue
        h, e = np.histogram(np.log1p(pos), bins=512)
        i = min(int(np.argmax(h[1:])) + 1, len(e) - 2)   # right edge of the modal bin
        floors[m] = float(np.expm1(e[1:][i]))
    return SC.marker_scale(pool, floors)


def _load_islet_polygons(path: pathlib.Path) -> tuple:
    """``(core, expanded)`` shapely geometries in MICRONS, keyed by islet name.

    Same reader contract as ``kintsugi/ring_width.py``: the class is in
    ``properties.classification.name`` and the 20 um ring repeats the core's name with an
    ``_exp20um`` suffix.
    """
    from shapely.geometry import shape
    from shapely.affinity import scale as affine_scale

    features = json.loads(path.read_text())
    features = features["features"] if isinstance(features, dict) else features
    core, expanded = {}, {}
    for f in features:
        name = (f["properties"].get("name") or "").replace("_exp20um", "")
        klass = (f["properties"].get("classification") or {}).get("name")
        g = shape(f["geometry"])
        if not g.is_valid:
            g = g.buffer(0)
        g = affine_scale(g, xfact=PIXEL_UM, yfact=PIXEL_UM, origin=(0, 0))
        (expanded if klass == "IsletExpanded" else core)[name] = g
    return core, expanded


def islet_regions(donor: str, df: pd.DataFrame) -> pd.DataFrame:
    """``cell_region`` in {core, peri, tissue} plus ``islet_num``, by whichever route exists.

    Returns the column untouched when the QuPath parent already carried it (6533).  When it
    did not (0112, whose every cell reads ``tissue``), membership is recomputed from the islet
    GeoJSON by point-in-polygon -- the same prior, arrived at geometrically.  A cell inside a
    core polygon is ``core``; inside the 20 um ring but not the core is ``peri``.
    """
    have_regions = set(df["cell_region"].dropna().unique()) & {"core", "peri"}
    if have_regions:
        return df[["cell_region", "islet_num"]].copy()

    name = ISLET_GEOJSON.get(donor)
    if name is None:
        print(f"  {donor}: no islet annotation by either route -- no prior, every reference "
              f"competes on every cell")
        return df[["cell_region", "islet_num"]].copy()

    from shapely import STRtree, points as shp_points

    core, expanded = _load_islet_polygons(DATA_ROOT / name)
    xy = shp_points(df["X_centroid"].to_numpy("float64"), df["Y_centroid"].to_numpy("float64"))
    region = np.full(len(df), "tissue", dtype=object)
    islet = np.full(len(df), "", dtype=object)
    for label, geoms in (("peri", expanded), ("core", core)):      # core wins over peri
        if not geoms:
            continue
        names = list(geoms)
        tree = STRtree([geoms[n] for n in names])
        hit_pt, hit_poly = tree.query(xy, predicate="within")
        region[hit_pt] = label
        islet[hit_pt] = np.asarray(names, dtype=object)[hit_poly]
    print(f"  {donor}: islet membership from {name} -- "
          f"{int((region == 'core').sum()):,} core, {int((region == 'peri').sum()):,} peri, "
          f"{len(core):,} islet polygons")
    return pd.DataFrame({"cell_region": region, "islet_num": islet}, index=df.index)


def call_donor(donor: str, Rm: pd.DataFrame, reqs: dict, intruders: list,
               require_min: float) -> tuple:
    """One donor's protein call: two regions, two reference sets, one prior.

    Identical in structure to ``run_islet.run_donor`` -- inside ``Islet_N`` only intruders may
    displace the endocrine prior and a cell showing none is subtyped by hormone; in the ring
    and the rest of the tissue there is no prior and every reference competes.  What differs
    is the scale (per donor) and that the whole section is called, not only the islets, since
    for these two donors the protein arm is the prior the RNA arm will consume.
    """
    df, have, unmeasured = load_expression(donor, list(Rm.index))
    alive = audit_unmeasured(donor, unmeasured, Rm, reqs)
    Rm = Rm.loc[have, alive]
    Rm = Rm.loc[:, (Rm != 0).any(axis=0)]
    regions = islet_regions(donor, df)
    df = df.assign(cell_region=regions["cell_region"], islet_num=regions["islet_num"])
    df = df[df["qc_analysis_eligible"]].reset_index(drop=True)
    if df.empty:
        return pd.DataFrame(), pd.DataFrame()

    scale = donor_scale(df, have)
    unscalable = [m for m in have if not np.isfinite(scale.loc[m, "top"])]
    if unscalable:
        raise SystemExit(f"donor {donor}: no usable scale for {unscalable}")
    Xs = SC.apply_scale(df[have].to_numpy("float32"), have, scale)

    # The DAPI ratio is a per-reference decision (run_step3.NUCLEAR_RATIO).  On the current
    # rule table no reference qualifies -- flatten() composes every leaf with pan-epithelial
    # or pan-immune ancestors -- so this is inert here and says so rather than being dropped.
    nuc_all = SC.nuclear_only_references(Rm, panel.NUCLEUS_MEASURED)
    nuc_by_ref = dict(zip(Rm.columns, nuc_all))
    if nuc_all.any():
        dfr, _ = SC.nuclear_ratio(df.copy(), have, "DAPI", panel.NUCLEUS_MEASURED)
        Xs_nuc = SC.apply_scale(dfr[have].to_numpy("float32"), have,
                                donor_scale(dfr, have))
    else:
        Xs_nuc = None

    usable = list(Rm.columns)
    intruder = [p for p in usable if p in intruders]
    print(f"  {donor}: {len(usable)} references evaluable, {len(intruder)} of them may "
          f"displace the endocrine prior")
    core = (df["cell_region"] == "core").to_numpy()

    def _call(rows, refs):
        if not rows.any():
            return None
        nuc_ref = np.array([nuc_by_ref.get(r, False) for r in refs], dtype=bool)
        S, allowed = SC.score_and_gate(Xs[rows], Xs_nuc[rows] if Xs_nuc is not None else None,
                                       Rm.loc[have, refs].to_numpy("float32"),
                                       have, reqs, refs, require_min, nuc_ref)
        o = SC.call_from_scores(S, refs, allowed=allowed)
        o.index = np.flatnonzero(rows)
        return o

    parts = [x for x in (_call(core, intruder), _call(~core, usable)) if x is not None]
    out = pd.concat(parts).sort_index().reset_index(drop=True)

    # Inside Islet_N a cell showing no intruder is endocrine BY THE PRIOR, and is subtyped by
    # hormone rather than left Ambiguous.
    sub = RI.subtype_endocrine(Xs, have, require_min)
    fell_through = core & (out["call"] == "Ambiguous").to_numpy()
    for c in ("call", "runner_up"):
        out[c] = np.where(fell_through, sub[c].to_numpy(dtype=object),
                          out[c].to_numpy(dtype=object))
    for c in ("score", "margin"):
        out[c] = np.where(fell_through, sub[c].to_numpy("float64"),
                          out[c].to_numpy("float64"))
    out["reason"] = np.where(fell_through, "endocrine_by_islet_prior",
                             out.get("reason", pd.Series([""] * len(out))).to_numpy(dtype=object))

    out.insert(0, "object_id", df["object_id"].to_numpy())
    out.insert(0, "donor_id", donor)
    out["cell_region"] = df["cell_region"].to_numpy()
    out["islet_num"] = df["islet_num"].to_numpy()
    out["X_centroid"] = df["X_centroid"].to_numpy()
    out["Y_centroid"] = df["Y_centroid"].to_numpy()
    return out, scale


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--donors", nargs="*", default=list(DONORS))
    ap.add_argument("--require-min", type=float, default=S3.REQUIRE_MIN)
    ap.add_argument("--encoding", default=S3.ENCODING, choices=("zero", "signed"))
    args = ap.parse_args()

    raw = R.load_rules()
    flat = R.flatten(raw)
    Rm, _ = R.reference_matrix(flat, panel.markers(include_variant=True), args.encoding,
                               strict=False)
    Rm = Rm.loc[(Rm != 0).any(axis=1)]
    reqs = R.requirements(raw)
    intruders = RI.intruder_references(raw, Rm.columns)
    excluded = [c for c in Rm.columns if c not in intruders]
    print(f"protein arm — {Rm.shape[1]} references over {Rm.shape[0]} markers, "
          f"encoding={args.encoding}, require_min={args.require_min}")
    print(f"  inside Islet_N, {len(intruders)} references may displace the endocrine prior")
    print(f"  excluded as implied by it ({len(excluded)}): {excluded}")

    OUT.mkdir(parents=True, exist_ok=True)
    for d in args.donors:
        t = time.time()
        calls, scale = call_donor(d, Rm, reqs, intruders, args.require_min)
        if calls.empty:
            print(f"  {d}: no analysis-eligible cells"); continue
        p = OUT / f"donor_id={d}"; p.mkdir(parents=True, exist_ok=True)
        calls.to_parquet(p / "data.parquet", index=False)
        scale.to_csv(p / "scale.csv")
        print(f"  {d}: {len(calls):,} cells called in {time.time()-t:.0f}s -> {p}")
        for region in ("core", "peri", "tissue"):
            g = calls[calls.cell_region == region]
            if g.empty:
                continue
            print(f"\n  --- {region}, {len(g):,} cells ---")
            for c, n in g["call"].value_counts().head(12).items():
                print(f"      {c:<40}{n:>10,}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
