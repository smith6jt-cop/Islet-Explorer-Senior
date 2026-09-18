#!/usr/bin/env python
"""Protein INS and GCG inside each Xenium cell's OWN polygon, against that cell's own RNA.

    cd core/phenocycler
    ../../.venv/bin/python scripts/restain/protein_rna.py --donors 6436 6442

**What changed, and why.** The earlier arm sampled a 4 um disc at a warped centroid and then
restricted to cells inside a QuPath ``Islet`` polygon. The maintainer's instruction of
2026-09-01 removes both: the QuPath islet annotations are **not to be used on Xenium or on the
post-Xenium IF**, and the comparison is to run on the Xenium cell and nucleus segmentations
against the IF insulin and glucagon channels only.

That is a stronger design, not merely a different one. Every object here is Xenium's own --
``cell_boundaries.parquet`` and ``nucleus_boundaries.parquet`` are segmentations of the same
DAPI the registration was fitted on -- so a cell's protein and a cell's RNA are the same cell
by construction, and no annotation, no islet call and no third party's polygon enters the path.
A disc also straddles neighbours at 9 um spacing; a polygon does not.

**Scored against the subpanel, not the whole panel.** ``core/panels``' 18-gene pancreas
endocrine identity core is the authority, and most of the 5,101-gene panel is not expressed in
pancreas at all. An earlier version ranked every gene and reported the top 25; that list was
mostly noise, and reading a cell type off the top of it is what produced a wrong conclusion.

**No cell is labelled here, and no AUROC is computed.** An AUROC needs a confirmed positive
class, and this lane has none: nobody has verified that a cell with one detected transcript is a
cell of that type. An earlier version defined "positive = count > 0", scored an AUROC against
it, and reported the result as though the label were truth -- an invented ground truth, which
this repository's rules forbid. What is reported is Spearman rho between the per-cell IF
intensity and the per-cell transcript count, plus the median IF intensity at each count level.
Those are statements about two measurements agreeing; they assert nothing about what a cell is.

**Defining a positive is the maintainer's call, not this file's.** If a truth set is ever
supplied, a classification statistic can be computed against it here.

**A whole-section AUROC is the wrong instrument for a focal signal, and this file cannot fix
that on its own.** Donor 6436 is T1D: insulin survives in a minority of lobules and is absent
from the rest, which is what the disease is. Averaged over the section its INS channel looks
inert -- and it is not, it is bright and cellular in the islets that still have beta cells. Any
number in this file that is computed over all cells inherits that, so read it beside
``protein_rna_figures.py --window INS``, which shows one islet at a time.
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(HERE.parents[1]))
sys.path.insert(0, str(HERE.parents[1] / "scripts" / "xenium_celltype"))

import hormone as HZ                                             # noqa: E402
import xenio                                                     # noqa: E402
from run_hormone import IMAGES, SECTIONS                         # noqa: E402
from run_register import OUT, XENIUM_PATHS, discover, donor_of    # noqa: E402
from run_register import GEOJSON_DIR                             # noqa: E402

#: Resolution the protein is read at, in microns per pixel. 0.5 um/px puts ~280 pixels inside a
#: 9 um cell and ~110 inside a 6 um nucleus, which is a mean over a cell rather than over a
#: handful of pixels. It is a real pyramid level of these qptiffs (level 1), so no resampling
#: happens and the box-mean origin rule stays exact.
PROTEIN_UM_PER_PX = 0.5

#: The subpanel this lane scores against: ``core/panels``' pancreas endocrine identity core,
#: read through the regenerable snapshot. Not the whole 5,101-gene panel -- most of it is not
#: expressed in pancreas, and ranking all of it produced a list that invited over-reading.
SUBPANEL = "identity_core/pancreas/01_pancreas_endocrine.csv"


class ProteinRnaError(RuntimeError):
    """The per-cell comparison cannot be made as asked -- never approximated."""


def warp_polygons(polys, transform):
    """Xenium-micron polygons -> IF-micron polygons.

    ``transform.apply`` is the forward direction: this lane fits Xenium -> IF, and the cells
    being warped are Xenium's. Getting this backwards does not error, it puts every polygon
    ~2 x (714, 282) um away on 6436 and every protein value becomes a sample of empty slide.
    """
    import shapely

    return [shapely.transform(g, lambda a: transform.apply(a)) for g in polys]


def pyramid_level_for(qptiff: Path, um_per_px: float, base_um_per_px: float) -> tuple[int, int]:
    """The level whose scale is exactly ``um_per_px``, and its integer factor.

    Exact rather than nearest. A level that is merely close would put the label raster and the
    image on different grids, and the resulting error grows across the slide instead of showing
    up as an offset.
    """
    import export_dapi as ED

    import tifffile
    with tifffile.TiffFile(str(qptiff)) as tif:
        series = tif.series[0]
        shapes = [(int(lv.shape[-2]), int(lv.shape[-1]))
                  for lv in getattr(series, "levels", [series])]
    want = int(round(um_per_px / base_um_per_px))
    for lvl in range(len(shapes)):
        f = ED.pyramid_factor(shapes, lvl)
        if f == want:
            return lvl, f
    raise ProteinRnaError(
        f"{qptiff.name} has no pyramid level at {um_per_px} um/px "
        f"(base {base_um_per_px}, factors {[ED.pyramid_factor(shapes, l) for l in range(len(shapes))]})")


def mean_in_labels(image: np.ndarray, labels: np.ndarray, n: int) -> np.ndarray:
    """Mean of ``image`` over each label 1..n, NaN where a label drew no pixel.

    NaN rather than 0: a cell that fell outside the imaged area has an unavailable protein
    value, not a dark one, and a 0 would rank it at the bottom of every statistic below.
    """
    from scipy import ndimage

    counts = np.bincount(labels.ravel(), minlength=n + 1)[1:]
    sums = ndimage.sum_labels(image, labels, index=np.arange(1, n + 1))
    out = np.full(n, np.nan, dtype=np.float64)
    hit = counts > 0
    out[hit] = sums[hit] / counts[hit]
    return out


def gene_rank_association(counts, values: np.ndarray, gene_cols: list[int]) -> dict:
    """Association between two MEASUREMENTS on the same cell. No labels, no truth set.

    **Why this is not an AUROC.** An AUROC scores a ranking against a positive class, and this
    lane has no confirmed positive class: nobody has verified that a cell with one detected
    transcript of a gene is a cell of that type. An earlier version of this file defined
    "positive = count > 0", computed an AUROC against it, and reported the result as if the
    label were ground truth. It was a definition invented here, and the repository's own rule
    is that existing classifications are audit context and never a gate.

    So what is reported instead needs no class at all: Spearman rho between the per-cell IF
    intensity and the per-cell transcript count for one gene, plus the median IF intensity at
    each count level. Both are statements about two measurements agreeing or not. Neither
    asserts that any cell is anything, and the count-stratified medians are there so the shape
    can be read directly rather than through a coefficient.

    Cells with no IF value are dropped, never zero-filled: nothing was measured for them.
    """
    import scipy.sparse as sp
    from scipy.stats import rankdata

    v = np.asarray(values, dtype=np.float64)
    ok = np.isfinite(v)
    if ok.sum() < 1000:
        raise ProteinRnaError(f"only {int(ok.sum())} cells carry a protein value")
    r_if = rankdata(v[ok])
    X = sp.csc_matrix(counts)[ok][:, gene_cols].toarray()

    out = {}
    for j, col in enumerate(gene_cols):
        c = X[:, j].astype(np.float64)
        r_c = rankdata(c)
        rho = float(np.corrcoef(r_if, r_c)[0, 1]) if c.std() > 0 else float("nan")
        strata = {}
        for k in (0, 1, 2, 3):
            m = (c >= k) if k == 3 else (c == k)
            key = "3+" if k == 3 else str(k)
            strata[key] = {"n": int(m.sum()),
                           "median_if": float(np.median(v[ok][m])) if m.any() else None}
        out[col] = {"spearman_rho": rho, "n_cells": int(ok.sum()),
                    "n_cells_detected": int((c > 0).sum()), "median_if_by_count": strata}
    return out





def rank_contrast(a: np.ndarray, b: np.ndarray) -> np.ndarray:
    """``rank(a) - rank(b)`` per cell, on the cells where both were measured.

    **This is the statistic that separates alpha from beta, and neither channel alone can.**
    Insulin and glucagon are both confined to islets, so each on its own ranks *any* islet gene
    near the top -- which is why 6442's GCG channel returns a beta gene list. Subtracting the
    two within one cell removes everything the two channels share, which is islet membership,
    and leaves only what differs between them, which is subtype.

    Percentile ranks rather than raw values because the two channels are different fluorophores
    at different exposures on an 8-bit scanner: their units are not comparable and a difference
    of raw intensities would be dominated by whichever was brighter. A cell missing either side
    is NaN, never zero -- zero is the middle of this scale and would read as "equally both".
    """
    from scipy.stats import rankdata

    a = np.asarray(a, dtype=np.float64)
    b = np.asarray(b, dtype=np.float64)
    both = np.isfinite(a) & np.isfinite(b)
    out = np.full(a.shape, np.nan)
    if both.sum() < 2:
        return out
    n = float(both.sum())
    out[both] = rankdata(a[both]) / n - rankdata(b[both]) / n
    return out


def score_against_genes(values: np.ndarray, dd, by_hint: dict) -> dict:
    """One per-cell IF vector against the SUBPANEL's endocrine genes, grouped by subtype hint.

    Restricted to ``core/panels``' 18-gene pancreas endocrine identity core, which is the panel
    taxonomy's one authority. Ranking all 5,101 genes was a mistake and is gone: most of the
    panel is not expressed in pancreas, so a ranked list of everything is mostly noise dressed
    as a result.

    ``subtype_hint`` is carried through from the subpanel and is a HINT, in its own file's
    words. It is not used to label any cell here, and the per-hint means below are means of
    association statistics -- not a classification score, and not evidence that a cell is a
    beta cell.
    """
    gi = {g: i for i, g in enumerate(dd.genes)}
    wanted = [(g, hint or "pan") for hint, genes in by_hint.items() for g in genes if g in gi]
    cols = [gi[g] for g, _ in wanted]
    assoc = gene_rank_association(dd.counts, values, cols)

    out = {"n_cells_measured": int(np.isfinite(values).sum()),
           "statistic": "spearman_rho(IF intensity, transcript count) per cell",
           "no_ground_truth": ("no cell is labelled here; there is no confirmed positive class "
                               "for cell type in this lane and none is invented"),
           "genes": {}}
    for (g, hint), col in zip(wanted, cols):
        out["genes"][g] = {"subtype_hint": hint, **assoc[col]}
    for hint in ("beta", "alpha"):
        vals = [v["spearman_rho"] for v in out["genes"].values()
                if v["subtype_hint"] == hint and np.isfinite(v["spearman_rho"])]
        out[f"mean_rho_{hint}"] = float(np.mean(vals)) if vals else None
    return out




def polygons_per_cell(bundle: Path, compartment: str, cell_ids) -> list:
    """One geometry per Xenium cell, in Xenium microns, or ``None`` where there is none.

    Neither boundary table is one row per cell, and the two are not irregular in the same way:

    * ``cell_boundaries.parquet`` can carry more than one ring for a cell;
    * ``nucleus_boundaries.parquet`` is per **nucleus** -- a cell may have two, or none.

    Parts belonging to one cell are combined and drawn at a single label, so a two-nucleus cell
    contributes both nuclei to that cell's mean rather than whichever one a naive reindex
    happened to keep. A cell with no nucleus is ``None`` and ends up NaN, which is the honest
    value: nothing was measured for it, and it must not rank as dark.
    """
    import shapely

    from phenocycler.integration.boundaries import (load_xenium_nucleus_polygons,
                                                    load_xenium_polygons)

    if compartment == "cell":
        series = load_xenium_polygons(bundle / "cell_boundaries.parquet", cell_ids=cell_ids)
        owner = np.asarray(series.index, dtype=object)
        geoms = list(series.to_numpy())
    elif compartment == "nucleus":
        series, owner = load_xenium_nucleus_polygons(bundle / "nucleus_boundaries.parquet")
        geoms = list(series.to_numpy())
        owner = np.asarray(owner, dtype=object)
        want = set(map(str, cell_ids))
        keep = np.array([str(o) in want for o in owner], dtype=bool)
        geoms = [g for g, k in zip(geoms, keep) if k]
        owner = owner[keep]
    else:
        raise ProteinRnaError(f"unknown compartment {compartment!r}")

    by: dict[str, list] = {}
    for o, g in zip(owner, geoms):
        if g is None or g.is_empty:
            continue
        by.setdefault(str(o), []).append(g)

    out = []
    for cid in map(str, cell_ids):
        parts = by.get(cid)
        if not parts:
            out.append(None)
        elif len(parts) == 1:
            out.append(parts[0])
        else:
            # A MultiPolygon, deliberately not a union. `rasterize_labels` draws every part of
            # a multipart geometry at the same label, so the union buys nothing -- and
            # `union_all` raises TopologyException on real Xenium rings ("unable to assign free
            # hole to a shell"), which would drop a cell for a reason that does not matter here.
            flat = []
            for g in parts:
                flat.extend(list(g.geoms) if g.geom_type == "MultiPolygon" else [g])
            out.append(shapely.MultiPolygon(flat))
    return out


def measure_one(serial: str, sections: pd.DataFrame, *, roi: str = "panc",
                out_dir: Path = OUT, um_per_px: float = PROTEIN_UM_PER_PX) -> dict:
    import tifffile
    import zarr

    from phenocycler.integration.boundaries import load_xenium_polygons
    from phenocycler.integration.rasterize import GridSpec
    from rasterize_labels import rasterize_labels
    import validate as VAL
    from run_align import _refined_transform
    from run_register import SUBPANELS

    _core, by_hint = VAL.endocrine_core(SUBPANELS)
    donor = donor_of(serial, XENIUM_PATHS)
    row = sections[sections.xenium_serial == serial].iloc[0]
    qptiff = IMAGES / row.image_file
    if not qptiff.exists():
        raise ProteinRnaError(f"{qptiff} is absent")
    base_um = float(row.pixel_size_um)
    transform = _refined_transform(donor, roi, out_dir)

    dd = xenio.load_donor(donor, XENIUM_PATHS, roi=roi, prefer="bundle")
    bundle = xenio._bundle_dir(donor, XENIUM_PATHS, roi)
    print(f"[protein-rna] {donor}: {len(dd.cell_id):,} cells, {len(dd.genes):,} genes "
          f"({dd.source})", flush=True)

    lvl, factor = pyramid_level_for(qptiff, um_per_px, base_um)
    with tifffile.TiffFile(str(qptiff)) as tif:
        series = tif.series[0]
        H = int(series.levels[lvl].shape[-2]); W = int(series.levels[lvl].shape[-1])
    # The grid IS the image: same origin, same pixel size, same extent, so a label index and a
    # pixel index are the same thing and no second alignment can creep in between them.
    grid = GridSpec(origin_x_um=0.5 * um_per_px, origin_y_um=0.5 * um_per_px,
                    pixel_um=um_per_px, height=H, width=W)

    rec = {"donor_id": donor, "xenium_serial": serial, "roi": roi,
           "qptiff": str(qptiff), "xenium_bundle": str(bundle),
           "transform_source": "refined", "protein_um_per_px": um_per_px,
           "pyramid_level": lvl, "pyramid_factor": factor,
           "n_cells": int(len(dd.cell_id)), "n_genes": int(len(dd.genes)),
           "uses_qupath_annotations": False, "subpanel": SUBPANEL,
           "note": ("no QuPath annotation of any class is read here; every object is Xenium's "
                    "own segmentation, per the maintainer's instruction of 2026-09-01"),
           "compartments": {}}

    channels = {}
    for marker in ("INS", "GCG"):
        page = HZ.channel_index(qptiff, marker, row)
        channels[marker] = (page, HZ.channel_names(qptiff)[page])

    for compartment in ("cell", "nucleus"):
        geoms = polygons_per_cell(bundle, compartment, dd.cell_id)
        keep = np.array([g is not None for g in geoms])
        missing = int((~keep).sum())
        idx = np.nonzero(keep)[0]
        warped = warp_polygons([geoms[i] for i in idx], transform)
        lab = rasterize_labels(warped, grid, order="input")
        print(f"[protein-rna] {donor}: {compartment} {len(idx):,} polygons rasterised "
              f"({missing:,} absent)", flush=True)

        comp = {"n_polygons": int(len(idx)), "n_absent": missing, "markers": {}}
        measured: dict[str, np.ndarray] = {}
        with tifffile.TiffFile(str(qptiff)) as tif:
            store = tif.series[0].aszarr(level=lvl)
            try:
                z = zarr.open(store, mode="r")
                for marker, (page, _chname) in channels.items():
                    img = np.asarray(z[page]).astype(np.float32)
                    vals_sub = mean_in_labels(img, lab.labels, len(idx))
                    del img
                    vals = np.full(len(dd.cell_id), np.nan)
                    vals[idx] = vals_sub
                    measured[marker] = vals
            finally:
                store.close()

        scored = dict(measured)
        scored["INS_minus_GCG"] = rank_contrast(measured["INS"], measured["GCG"])

        for name, vals in scored.items():
            # The per-cell vector is kept because everything downstream -- cell groups, figures,
            # any re-scoring -- must use the SAME numbers these AUROCs were computed from.
            vp = out_dir / "protein_rna" / f"{donor}_{serial}_{compartment}_{name}.npy"
            vp.parent.mkdir(parents=True, exist_ok=True)
            np.save(vp, vals)
            entry = score_against_genes(vals, dd, by_hint)
            entry["values_npy"] = str(vp)
            if name in channels:
                entry.update({"qptiff_channel": channels[name][1],
                              "channel_page": int(channels[name][0]),
                              "recorded_role": name,
                              "median": float(np.nanmedian(vals)),
                              "p99": float(np.nanpercentile(vals, 99))})
            else:
                entry["definition"] = ("percentile rank of INS minus percentile rank of GCG, "
                                       "within one cell; removes what the two channels share "
                                       "(islet membership) and keeps what differs (subtype)")
            comp["markers"][name] = entry
            b, a = entry.get("mean_rho_beta"), entry.get("mean_rho_alpha")
            print(f"[protein-rna]   {compartment}/{name}: mean rho beta {b:+.3f} / "
                  f"alpha {a:+.3f}   "
                  + ", ".join(f"{g} {v['spearman_rho']:+.3f}"
                              for g, v in entry["genes"].items()
                              if v["subtype_hint"] in ("beta", "alpha")), flush=True)
        del lab
        rec["compartments"][compartment] = comp

    return rec


def registered_serials(roi: str = "panc") -> dict[str, str]:
    """serial -> donor, for every donor whose registration PASSED its gate.

    This is what the module actually needs: a refined transform that cleared the release gate.
    It used to enumerate from the QuPath GeoJSON directory -- three files -- while stating in three
    places that it reads no annotation, so ten donors with a PASS gate were skipped in silence and
    the lane looked like a registration failure. The gate verdict is the honest predicate, and it
    lives in the donor directory next to the transform it certifies.
    """
    import json as _json
    out: dict[str, str] = {}
    root = OUT / "valis"
    if not root.is_dir():
        return out
    for donor_dir in sorted(root.glob("donor_id=*")):
        donor = donor_dir.name.split("=", 1)[1]
        if not (donor_dir / "refined.json").exists():
            continue
        verdicts = []
        for gate in sorted(donor_dir.glob("alignment_*.json")):
            try:
                j = _json.loads(gate.read_text())
            except (OSError, ValueError):
                continue
            v = j.get("verdict") or (j.get("alignment") or {}).get("verdict")
            if v:
                verdicts.append(str(v))
        if "PASS" not in verdicts:
            continue
        try:
            export = _json.loads((donor_dir / "export.json").read_text())
        except (OSError, ValueError):
            continue
        if export.get("roi") != roi:
            continue
        serial = str(export.get("serial") or "")
        if serial:
            out[serial] = donor
    return out


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--donors", nargs="*", default=None)
    ap.add_argument("--roi", default="panc")
    ap.add_argument("--out", type=Path, default=OUT)
    ap.add_argument("--um-per-px", type=float, default=PROTEIN_UM_PER_PX)
    args = ap.parse_args(argv)

    sections = HZ.load_section_map(SECTIONS)
    (args.out / "protein_rna").mkdir(parents=True, exist_ok=True)
    recs = {}
    for serial, donor in sorted(registered_serials(args.roi).items()):
        if args.donors and donor not in args.donors:
            continue
        rec = measure_one(serial, sections, roi=args.roi, out_dir=args.out,
                          um_per_px=args.um_per_px)
        rec.setdefault("serial", serial)
        rec.setdefault("donor", donor)
        recs[serial] = rec
    if not recs:
        print("no donor matched", file=sys.stderr)
        return 1
    p = args.out / f"protein_rna_{args.roi}.json"
    # MERGE, never replace. `--donors 6539` used to rewrite this file with that one record and
    # drop 6436 and 6442 from it -- the per-donor .npy arrays survived, so nothing was
    # measured twice, but the index silently lost every donor not named on the command line.
    # Keyed by serial so a re-run of a donor replaces its own record and no other.
    if p.exists():
        try:
            prior = json.loads(p.read_text())
        except (OSError, json.JSONDecodeError):
            prior = []
        # Records written before `serial` existed carry None, so they cannot be matched by
        # key. Identify those by what they measured: a prior record whose n_cells equals one
        # just written is the same section, and the fresh record wins.
        measured_now = {r.get("n_cells") for r in recs.values()}
        for rec in prior if isinstance(prior, list) else []:
            key = rec.get("serial")
            if key is None:
                if rec.get("n_cells") in measured_now:
                    continue
                key = f"n_cells={rec.get('n_cells')}"
            recs.setdefault(key, rec)
    merged = sorted(recs.values(), key=lambda r: str(r.get("serial") or ""))
    p.write_text(json.dumps(merged, indent=2, allow_nan=False))
    print(f"wrote {p}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
