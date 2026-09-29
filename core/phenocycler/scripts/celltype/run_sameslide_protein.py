"""Protein-arm cell typing for the two same-slide Pheno-Xeno donors (6533, 0112).

This is Stage 5 of the streamlined Xenium plan, and under the user-directed order it runs
**first** -- before the RNA arm, not after it.  The reason is arithmetic: these bundles carry
~126 estimated false-positive transcripts against a median 433 per cell, so the RNA needs a
context prior and the protein is the only measurement dense enough to supply one.  70-plex on
the same physical section, against ~433 sparse transcripts of which ~29% are gDNA background.

**The lane's full step order runs.** Step 1 (marker usability, per donor) gates the markers
scored here and rejects 0 of 58 on both donors. Step 2 (UniFORM) was re-run with these two
donors inside the cohort -- ``shifts_22donor.json``, 22 donors, ``uniform_commit
c750a9a7...``, the commit ``uniform_adapter`` pins -- so this module reads the ALIGNED
partitions and takes the COHORT scale from them, exactly as ``run_step3`` and ``run_islet``
do. Step 3 (cosine behind the requirement gate) is reused unchanged.

An earlier version of this module derived a per-donor scale and carried the caveat that the
two donors were therefore not comparable to each other. Step 2 removes it, and the alignment
is visible in the values: median CD99 194.3 against 192.1 and median INS 9.8 against 11.2,
where the uncorrected intensities differed by up to 10x. Per-marker, the shift table records
the mode SD falling from 0.368 to 0.028.

Both steps reach these donors through the lane's own resolution, not through a bypass.
``panel.expression_path`` selects the stable expression root per donor and normalises ``0112``
and ``112`` to one donor; ``run_step3.ALIGNED`` is the lane's ``uniform`` directory. Where a
checkout keeps that data outside the repo, it is wired in by symlink rather than by an
override in the code. The one accommodation that remains in code is that a donor absent from
the acquisition table resolves its markers from its own partition -- see
``panel.markers_on_disk`` -- because acquisition is a fact about the run and the table has no
row for a donor acquired later.

Everything that decides a cell type is reused unchanged -- ``rules.flatten`` /
``reference_matrix`` / ``requirements``, ``score.score_and_gate`` / ``call_from_scores``, and
``run_islet.intruder_references`` / ``subtype_endocrine``.

**Inside ``Islet_N`` the decision is an ordered cascade**, not a competition -- see
:func:`_islet_core_cascade`.  The polygons were drawn by a QuPath pixel classifier on INS,
GCG, SST and CD99, so the cascade reads those four in the order the annotation was built:
hormone argmax, then CD99 for the hormone-negative endocrine cells, then marker profile for
the CD99-negative remainder.  That last step scores the non-endocrine references, which
``run_islet.intruder_references`` supplies -- it excludes the **whole epithelial branch**,
ancestors and siblings, because an endocrine cell satisfies those rules by construction.  A
hand-rolled filter that excluded only endocrine-prefixed names was tried and measured: the
generic ``Epithelial`` parent then took 21,450 of 65,363 islet cells.

The ring and the rest of the tissue have no prior and every reference competes, unchanged.

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

DONORS = ("6533", "112")

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

#: The structures that legitimately traverse an islet, as their own QuPath annotation layers,
#: and the reference each one makes a cell a candidate for. A vessel, a lymphatic and a nerve
#: are objects that were drawn, so membership is read by point-in-polygon exactly as islet
#: membership is -- not inferred from a marker and not inferred from a neighbourhood.
#:
#: The marker still has to agree. Being inside a drawn vessel makes a cell a CANDIDATE for
#: ``Endo-vascular``; it is called that only if it also meets what that reference requires
#: (CD31, plus one of CD34 / CD99 / Caveolin). Centroid and corresponding marker, both.
OBJECT_LAYERS = {"blood_vessel": "Blood_vessel", "lymphatic": "Lymphatic", "nerve": "Nerve"}
STRUCTURE_REFERENCE = {"blood_vessel": "Endo-vascular", "lymphatic": "Endo-lymphatic",
                       "nerve": "Neural"}


#: Step 1's output for these donors. The lane's step order is not optional: Step 1 decides
#: which measurements the later steps may lean on, and running Step 3 without it means
#: scoring on measurements nothing has vetted.
#:
#:     python scripts/celltype/run_step1.py --donors 6533 112 \
#:       --out ../../data/development/celltype/marker_donor_usability_sameslide.parquet
STEP1 = panel.OUT / "marker_donor_usability_sameslide.parquet"


def step1_unusable(donor: str) -> list:
    """Markers Step 1 says are not a valid measurement on this donor.

    ``usability.py`` stopped emitting a label because on the cohort every marker passed and
    the label carried no information -- but the three criteria it was built from are stated
    in the lane README and are applied here to Step 1's own columns, so the gate is the
    lane's and not this module's:

    * fewer than ``usability.MIN_CELLS`` cells,
    * more than half the cells reading exactly 0,
    * no spatial structure above this donor's own permutation null (neither Moran's I on
      ranks nor on the bright tail clears it).

    Measured on 6533 and 0112: **0 of 58 markers unusable on either donor**, the same answer
    Step 1 gave on all 1,120 cohort (donor, marker) cells. The gate excludes nothing here.
    That is a result, not a reason to skip the step.
    """
    if not STEP1.exists():
        raise SystemExit(
            f"Step 1 has not been run for these donors: {STEP1} is absent. Run\n"
            f"  python scripts/celltype/run_step1.py --donors 6533 112 --out {STEP1}\n"
            "first -- Step 3 is not entitled to measurements Step 1 has not vetted.")
    import usability as U

    d = pd.read_parquet(STEP1)
    d = d[d["donor_id"].astype(str).isin(panel._donor_spellings(donor))]
    if d.empty:
        raise SystemExit(f"{STEP1} carries no rows for donor {donor}; re-run Step 1 with "
                         f"--donors including it.")
    structure = np.maximum(d["morans_rank"].abs(), d["morans_top1"].abs())
    bad = ((d["n"] < U.MIN_CELLS) | (d["zero_frac"] > 0.5)
           | (structure <= d["morans_null"]))
    return sorted(d.loc[bad, "marker"].tolist())

def expression_path(donor: str) -> pathlib.Path:
    """The raw expression partition, resolved by ``panel`` -- stable root first, retired
    fallback second, with ``112`` and ``0112`` normalised to one donor there."""
    return panel.expression_path(donor)


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
    keep = ["object_id", "qc_analysis_eligible", "cell_region", "islet_num",
            "X_centroid", "Y_centroid"]
    df, measured, _ = S3._read_donor(donor, markers, keep, ratio=False)
    measured = [m for m in measured if np.isfinite(df[m].to_numpy("float64")).any()]
    unmeasured = [m for m in markers if m not in measured]
    # The object layers are keyed by IMAGE and the aligned partition does not carry it, so it
    # comes from the expression partition -- one value per donor, read once.
    path = expression_path(donor)
    if path.exists():
        df["image"] = str(pq.ParquetFile(path).read(columns=["image"]).to_pandas()
                          ["image"].iloc[0])
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


def _layer_path(image: str, kind: str) -> pathlib.Path | None:
    """The object layer for one image, named as QuPath exported it.

    QuPath writes ``<Class>__<image name>.geojson`` with the image's spaces and ``#``
    flattened to underscores, so the name is derivable from the ``image`` column rather than
    hard-coded per donor. A glob is tried as well, because the export names have already
    changed once in this project.
    """
    stem = image.replace(" ", "_").replace("#", "_")
    exact = DATA_ROOT / f"{OBJECT_LAYERS[kind]}__{stem}.geojson"
    if exact.exists():
        return exact
    hits = sorted(DATA_ROOT.glob(f"{OBJECT_LAYERS[kind]}__*.geojson"))
    hits = [h for h in hits if image.split(".")[0] in h.name]
    return hits[0] if hits else None


def _load_annotation_polygons(path: pathlib.Path) -> list:
    """Every annotation in a layer, as shapely geometries in MICRONS."""
    from shapely.affinity import scale as affine_scale
    from shapely.geometry import shape

    features = json.loads(path.read_text())
    features = features["features"] if isinstance(features, dict) else features
    out = []
    for f in features:
        g = shape(f["geometry"])
        if not g.is_valid:
            g = g.buffer(0)
        out.append(affine_scale(g, xfact=PIXEL_UM, yfact=PIXEL_UM, origin=(0, 0)))
    return out


def structure_membership(donor: str, df: pd.DataFrame) -> np.ndarray:
    """Which drawn structure each cell's centroid falls inside, or ``''``.

    Stamped for EVERY cell, in every region, because it is a fact about the cell and costs
    one spatial query. It is only allowed to DECIDE inside ``Islet_N`` -- see
    :func:`_islet_core_cascade` -- where the question "is this an intruding vessel, lymphatic
    or nerve" is the one the islet prior leaves open. Outside the islet nothing needs the
    shortcut: every reference already competes there, and the membership is reported beside
    the call instead so the two can be compared.

    A cell inside more than one layer takes the more specific one -- nerve, then lymphatic,
    then blood vessel -- because blood vessel is by far the largest layer (111,806 polygons
    against 11,293 and 8,535 on 6533) and would otherwise absorb the overlaps.
    """
    from shapely import STRtree, points as shp_points

    image = str(df["image"].iloc[0]) if "image" in df.columns else ""
    out = np.full(len(df), "", dtype=object)
    if not image:
        return out
    xy = shp_points(df["X_centroid"].to_numpy("float64"), df["Y_centroid"].to_numpy("float64"))
    found = []
    for kind in ("blood_vessel", "lymphatic", "nerve"):       # least to most specific
        path = _layer_path(image, kind)
        if path is None:
            continue
        geoms = _load_annotation_polygons(path)
        if not geoms:
            continue
        hit, _ = STRtree(geoms).query(xy, predicate="within")
        out[hit] = kind
        found.append(f"{kind} {len(geoms):,} polys -> {len(np.unique(hit)):,} cells")
    print(f"  {donor}: object layers -- " + ("; ".join(found) if found else "none found"))
    return out


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

    name = next((ISLET_GEOJSON[d] for d in panel._donor_spellings(donor)
                 if d in ISLET_GEOJSON), None)
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



#: The two reasons that mean "this cell is endocrine", as opposed to a marker-profile call.
ENDOCRINE_REASONS = ("hormone_argmax", "cd99_islet_endocrine")


def _islet_core_cascade(Xs: np.ndarray, markers: list, core: np.ndarray, intruders: list,
                        require_min: float, score_rows, structure: np.ndarray,
                        available: set) -> pd.DataFrame | None:
    """Inside `Islet_N`: hormone, then CD99, then everything else. In that order.

    The islet polygons were drawn by a QuPath pixel classifier on **INS, GCG, SST and CD99**,
    so those four are what the annotation is made of and the cascade reads them in the order
    the annotation was built:

    1. **Hormone.** The argmax over INS/GCG/SST, if it clears the floor, IS the call --
       ``Endo-beta`` / ``Endo-alpha`` / ``Endo-delta``. Nothing else is consulted, because
       inside a hormone-drawn polygon a hormone-bright cell is not an open question.
    2. **CD99.** A cell with no hormone but bright CD99 is still endocrine, hormone-negative
       -- a real category in T1D islets, and the one the pixel classifier's fourth channel
       exists to catch. CD99 earns this: on the per-donor scale it reads positive in 82.7%
       (6533) and 65.0% (0112) of core cells against 4.6% / 3.8% in the surrounding ring,
       an 18-75x enrichment, so it separates islet-endocrine from islet-non-endocrine rather
       than merely being present.
    3. **Structure, by centroid and by marker.** A vessel, a lymphatic and a nerve traverse
       islets, and each was drawn as its own QuPath annotation layer. A CD99-negative cell
       whose centroid sits inside one is a candidate for that structure's reference -- and is
       called it only if it also meets what the reference requires. Centroid AND corresponding
       marker; neither alone.
    4. **Everything else** is typed by marker profile against the non-endocrine reference set.

    This replaces an earlier arrangement in which endocrine references were removed from
    contention inside the core and a cell became endocrine only by meeting NO intruder rule.
    That ordering had no endocrine escape hatch, and it showed: 98.9% of the cells it called
    ``Imm-NK-cell`` inside 6533's cores cleared the hormone floor, 1,969 of them SST-top, at
    a median Granzyme B and CD3e of 0.00. A delta cell that happened to satisfy
    CD56 + CD57 + any one of nine pan-immune markers had nowhere else to go. Here the hormone
    is read first and the question never arises.
    """
    if not core.any():
        return None
    idx = np.flatnonzero(core)
    m = {k: i for i, k in enumerate(markers)}

    sub = RI.subtype_endocrine(Xs[idx], markers, require_min)
    by_hormone = (sub["call"] != "Endo-hormone-negative").to_numpy()
    cd99 = (Xs[idx, m["CD99"]] >= require_min) if "CD99" in m else np.zeros(len(idx), bool)
    by_cd99 = ~by_hormone & cd99
    rest = ~by_hormone & ~cd99

    out = sub.copy()
    # "pending" is deliberate: step 4 claims every row it is handed, so a surviving "pending"
    # after the cascade would be a routing bug, not a cell type. The assert below says so.
    out["reason"] = np.where(by_hormone, "hormone_argmax",
                             np.where(by_cd99, "cd99_islet_endocrine", "pending"))
    out["n_candidates"] = np.nan
    def _write(rows_in_core: np.ndarray, refs: list, reason: str) -> np.ndarray:
        """Score `rows_in_core` against `refs`; return who was actually called."""
        if not rows_in_core.any() or not refs:
            return np.zeros(len(idx), dtype=bool)
        mask = np.zeros(len(core), dtype=bool)
        mask[idx[rows_in_core]] = True
        scored = score_rows(mask, refs)
        scored.index = np.flatnonzero(rows_in_core)
        # A candidacy step only claims the cells it actually called. The last step claims
        # every cell it was given, Ambiguous included -- there is nothing after it.
        keep = ((scored["call"] != "Ambiguous").to_numpy() if reason != "marker_profile"
                else np.ones(len(scored), dtype=bool))
        take = scored.index[keep]
        for c in out.columns:
            if c in scored.columns:
                out.loc[take, c] = scored.loc[take, c].to_numpy()
        got = scored.loc[take, "reason"].to_numpy(dtype=object)
        out.loc[take, "reason"] = np.where(got == "", reason, got)
        done = np.zeros(len(idx), dtype=bool)
        done[take] = True
        return done

    if rest.any():
        # Step 3 -- the drawn structures, one reference each, each behind its own requirement
        # gate. A cell inside a vessel that shows no CD31 is NOT called a vessel; it falls
        # through to step 4 with every other reference still available to it.
        placed = np.zeros(len(idx), dtype=bool)
        for kind, ref in STRUCTURE_REFERENCE.items():
            if ref not in available:
                continue
            sel = rest & ~placed & (structure[idx] == kind)
            placed |= _write(sel, [ref], f"structure_{kind}")
        # Step 4 -- marker profile. `intruders` is run_islet.intruder_references: the whole
        # epithelial branch is out, because an endocrine cell satisfies those rules by
        # construction and this branch has already established the cell is not endocrine.
        _write(rest & ~placed, intruders, "marker_profile")
    stranded = int((out["reason"] == "pending").sum())
    if stranded:
        raise AssertionError(f"{stranded} core cell(s) left the cascade unrouted")
    out.index = idx
    return out

def call_donor(donor: str, Rm: pd.DataFrame, reqs: dict, scale: pd.DataFrame,
               intruders: list, require_min: float) -> tuple:
    """One donor's protein call: two regions, two reference sets, one prior.

    Identical in structure to ``run_islet.run_donor`` -- inside ``Islet_N`` only intruders may
    displace the endocrine prior and a cell showing none is subtyped by hormone; in the ring
    and the rest of the tissue there is no prior and every reference competes.  What differs
    is the scale (per donor) and that the whole section is called, not only the islets, since
    for these two donors the protein arm is the prior the RNA arm will consume.
    """
    df, have, unmeasured = load_expression(donor, list(Rm.index))
    rejected = [m for m in step1_unusable(donor) if m in have]
    print(f"  {donor}: Step 1 vetted {len(have)} rule markers, rejected {len(rejected)}"
          + (f": {rejected}" if rejected else ""))
    if rejected:
        have = [m for m in have if m not in rejected]
        unmeasured = sorted(set(unmeasured) | set(rejected))
    alive = audit_unmeasured(donor, unmeasured, Rm, reqs)
    Rm = Rm.loc[have, alive]
    Rm = Rm.loc[:, (Rm != 0).any(axis=0)]
    regions = islet_regions(donor, df)
    df = df.assign(cell_region=regions["cell_region"], islet_num=regions["islet_num"],
                   structure=structure_membership(donor, df))
    df = df[df["qc_analysis_eligible"]].reset_index(drop=True)
    if df.empty:
        return pd.DataFrame(), pd.DataFrame()

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
                                S3._nuclear_scale(scale, [donor], have))
    else:
        Xs_nuc = None

    usable = list(Rm.columns)
    intruder = [p for p in usable if p in intruders]
    print(f"  {donor}: {len(usable)} references evaluable, {len(intruder)} of them non-endocrine")
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

    core_out = _islet_core_cascade(Xs, have, core, intruder, require_min, _call,
                                   df["structure"].to_numpy(dtype=object), set(usable))
    parts = [x for x in (core_out, _call(~core, usable)) if x is not None]
    out = pd.concat(parts).sort_index().reset_index(drop=True)
    out["islet_endocrine"] = out["reason"].isin(ENDOCRINE_REASONS).to_numpy()

    out.insert(0, "object_id", df["object_id"].to_numpy())
    out.insert(0, "donor_id", donor)
    out["cell_region"] = df["cell_region"].to_numpy()
    out["islet_num"] = df["islet_num"].to_numpy()
    out["structure"] = df["structure"].to_numpy()
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
    print(f"  inside Islet_N: hormone -> CD99 -> marker profile; the last step scores the "
          f"{len(intruders)} non-endocrine references")
    print(f"  not in that set, endocrine + epithelial branch ({len(excluded)}): {excluded}")

    # The cohort floor and top, on the ALIGNED values -- the whole point of having run Step 2.
    #
    # Built over EVERY aligned donor on disk, never over `--donors`. Deriving it from the
    # donors being called would make `--donors 6533` silently produce a per-donor scale --
    # the exact thing Step 2 exists to remove -- and it would do so without changing a single
    # message. The scale is a property of the alignment, not of the request.
    scale_donors = sorted(p.name.removeprefix("donor_id=") for p in S3.ALIGNED.glob("donor_id=*")
                          if (p / "data.parquet").is_file())
    if not scale_donors:
        raise SystemExit(f"no aligned partitions under {S3.ALIGNED}; Step 2 must run first")
    missing = [d for d in args.donors
               if not any(panel._same_donor(d, s) for s in scale_donors)]
    if missing:
        raise SystemExit(f"donor(s) {missing} have no aligned partition under {S3.ALIGNED}")
    scale = S3.build_scale(scale_donors, list(Rm.index), ratio=S3.NUCLEAR_RATIO)
    print(f"  cohort scale over {len(scale_donors)} aligned donor(s) {scale_donors}, "
          f"{int(scale['top'].notna().sum())} markers")

    OUT.mkdir(parents=True, exist_ok=True)
    for d in args.donors:
        t = time.time()
        calls, scale_used = call_donor(d, Rm, reqs, scale, intruders, args.require_min)
        if calls.empty:
            print(f"  {d}: no analysis-eligible cells"); continue
        p = OUT / f"donor_id={d}"; p.mkdir(parents=True, exist_ok=True)
        calls.to_parquet(p / "data.parquet", index=False)
        scale_used.to_csv(p / "scale.csv")
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
