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

**Nothing here is told which genes are endocrine.** The headline is, for every gene on the
panel, the AUROC with which the IF channel ranks the cells that express it above the cells that
do not -- a Mann-Whitney statistic, so there is no threshold on the protein side (it is a rank)
and none on the RNA side (detected or not). If the registration and the channel assignment are
both right, the beta genes must come out at the top for insulin and the alpha genes for
glucagon, *without having been named*. That makes this a test the data can fail.

It is also, incidentally, what settles which antibody a channel holds -- ``channel_check.py``
showed that 6436's only live hormone channel is the one recorded GCG, and a ranked gene list is
evidence about that recording which a mean intensity is not.
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

#: How many top-ranked genes to report per channel. Enough to see whether a coherent cell type
#: comes out rather than one lucky gene.
TOP_N = 25


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


def gene_auroc(counts, values: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Per gene: AUROC with which ``values`` ranks expressing cells above non-expressing ones.

    Mann-Whitney U over ranks, computed for all genes at once as one sparse mat-vec rather than
    a loop -- 5,101 genes against ~700k cells is not a loop worth writing.

    No threshold anywhere. The protein side is a rank, so its scale and its background level
    are irrelevant; the RNA side is detected-or-not, which is the natural dichotomy for counts
    this sparse and needs no cutoff to be chosen.
    """
    import scipy.sparse as sp
    from scipy.stats import rankdata

    v = np.asarray(values, dtype=np.float64)
    ok = np.isfinite(v)
    if ok.sum() < 1000:
        raise ProteinRnaError(f"only {int(ok.sum())} cells carry a protein value")
    r = rankdata(v[ok])                                   # 1-based, ties averaged
    X = sp.csc_matrix(counts[ok])
    X.data = np.ones_like(X.data)                         # detected or not
    X = X.tocsc()
    n_pos = np.asarray(X.sum(axis=0)).ravel()
    rank_sum = X.T.dot(r)
    n = int(ok.sum())
    n_neg = n - n_pos
    auroc = np.full(X.shape[1], np.nan)
    good = (n_pos > 0) & (n_neg > 0)
    auroc[good] = ((rank_sum[good] - n_pos[good] * (n_pos[good] + 1) / 2.0)
                   / (n_pos[good] * n_neg[good]))
    return auroc, n_pos


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


def score_against_genes(values: np.ndarray, dd, top_n: int) -> dict:
    """One per-cell vector -> the genes it ranks highest and lowest, and how many cells it had."""
    auroc, n_pos = gene_auroc(dd.counts, values)
    order = np.argsort(np.where(np.isfinite(auroc), auroc, -np.inf))[::-1]
    def rows(ix):
        return [{"gene": dd.genes[i], "auroc": float(auroc[i]),
                 "n_cells_expressing": int(n_pos[i])} for i in ix]
    finite = np.isfinite(auroc)
    return {"n_cells_measured": int(np.isfinite(values).sum()),
            "n_genes_scored": int(finite.sum()),
            "top_genes": rows(order[:top_n]),
            "bottom_genes": rows([i for i in order[::-1] if finite[i]][:top_n])}


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
                out_dir: Path = OUT, um_per_px: float = PROTEIN_UM_PER_PX,
                top_n: int = TOP_N) -> dict:
    import tifffile
    import zarr

    from phenocycler.integration.boundaries import load_xenium_polygons
    from phenocycler.integration.rasterize import GridSpec
    from rasterize_labels import rasterize_labels
    from run_align import _refined_transform

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
           "uses_qupath_annotations": False,
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
            entry = score_against_genes(vals, dd, top_n)
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
            top = entry["top_genes"]
            print(f"[protein-rna]   {compartment}/{name}: "
                  + ", ".join(f"{t['gene']} {t['auroc']:.3f}" for t in top[:6]), flush=True)
            if name == "INS_minus_GCG":
                print("[protein-rna]                      lowest: "
                      + ", ".join(f"{t['gene']} {t['auroc']:.3f}"
                                  for t in entry["bottom_genes"][:6]), flush=True)
        del lab
        rec["compartments"][compartment] = comp

    return rec


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--donors", nargs="*", default=None)
    ap.add_argument("--roi", default="panc")
    ap.add_argument("--geojson-dir", type=Path, default=GEOJSON_DIR,
                    help="only used to enumerate which sections have a refined transform; no "
                         "annotation inside these files is read")
    ap.add_argument("--out", type=Path, default=OUT)
    ap.add_argument("--um-per-px", type=float, default=PROTEIN_UM_PER_PX)
    ap.add_argument("--top-n", type=int, default=TOP_N)
    args = ap.parse_args(argv)

    sections = HZ.load_section_map(SECTIONS)
    (args.out / "protein_rna").mkdir(parents=True, exist_ok=True)
    recs = []
    for serial, _gj in discover(args.geojson_dir).items():
        donor = donor_of(serial, XENIUM_PATHS)
        if args.donors and donor not in args.donors:
            continue
        recs.append(measure_one(serial, sections, roi=args.roi, out_dir=args.out,
                                um_per_px=args.um_per_px, top_n=args.top_n))
    if not recs:
        print("no donor matched", file=sys.stderr)
        return 1
    p = args.out / f"protein_rna_{args.roi}.json"
    p.write_text(json.dumps(recs, indent=2, allow_nan=False))
    print(f"wrote {p}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
