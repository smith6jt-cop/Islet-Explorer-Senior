#!/usr/bin/env python
"""One islet, six ways: the IF as imaged, the protein per Xenium cell, and that cell's own RNA.

    cd core/phenocycler
    ../../.venv/bin/python scripts/restain/protein_rna_figures.py --donors 6436 6442

The point of the sheet is that the middle row is drawn from **Xenium's own segmentation**,
warped onto the IF by the registration, and filled with the number `protein_rna.py` actually
measured for that cell -- the same array its gene statistics were computed from, loaded from
disk rather than re-extracted. So what is on the page is what the numbers say, and the two can
be compared by eye on the same cells.

No QuPath annotation of any class is read. The window is chosen by **Xenium RNA alone** -- the
densest neighbourhood of endocrine-committed cells -- so nothing about where to look comes from
the IF or from a drawn polygon.

The RNA panel is detection, not a threshold: a cell is beta-committed if it detects a beta gene
of the identity core and no alpha one, alpha-committed the other way round, and both/neither are
drawn as themselves rather than forced into one of the two.
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(HERE.parents[1]))
sys.path.insert(0, str(HERE.parents[1] / "scripts" / "xenium_celltype"))

import display as DISP                                           # noqa: E402
import hormone as HZ                                             # noqa: E402
import protein_rna as PR                                         # noqa: E402
import validate as VAL                                           # noqa: E402
import xenio                                                     # noqa: E402
from run_hormone import IMAGES, SECTIONS                         # noqa: E402
from run_register import (GEOJSON_DIR, OUT, SUBPANELS,           # noqa: E402
                          XENIUM_PATHS, discover, donor_of)

#: Half-width of the drawn window, in microns. An islet is 50-250 um across; this frames one
#: with exocrine around it, so a channel that is merely bright everywhere cannot look specific.
HALF_UM = 130.0

BETA_C, ALPHA_C, BOTH_C, NEITHER_C = "#1c7ed6", "#f76707", "#ae3ec9", "#dee2e6"

#: Percentile of the per-cell values, over every cell in the section, that tops the colour
#: scale. 99.5 rather than 100 because a handful of cells reach 200+ counts and would flatten
#: everything else; stated on every panel so the reader knows what saturating means.
VALUE_PERCENTILE = 99.5


def committed(dd, by_hint) -> tuple[np.ndarray, np.ndarray]:
    """Beta- and alpha-committed cells, by detection of the identity core only."""
    gi = {g: i for i, g in enumerate(dd.genes)}
    X = dd.counts > 0

    def anyof(genes):
        cols = [gi[g] for g in genes if g in gi]
        if not cols:
            raise PR.ProteinRnaError("no gene of this set is on the panel")
        return np.asarray(X[:, cols].sum(axis=1)).ravel() > 0

    return anyof(by_hint["beta"]), anyof(by_hint["alpha"])


def figure_one(serial: str, sections, *, roi="panc", out_dir: Path = OUT,
               half_um: float = HALF_UM) -> dict:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import shapely
    import tifffile
    import zarr
    from matplotlib.collections import PolyCollection

    from run_align import _refined_transform

    donor = donor_of(serial, XENIUM_PATHS)
    row = sections[sections.xenium_serial == serial].iloc[0]
    qptiff = IMAGES / row.image_file
    base_um = float(row.pixel_size_um)
    transform = _refined_transform(donor, roi, out_dir)
    dd = xenio.load_donor(donor, XENIUM_PATHS, roi=roi, prefer="bundle")
    _core, by_hint = VAL.endocrine_core(SUBPANELS)
    beta, alpha = committed(dd, by_hint)

    # --- where to look: chosen by RNA alone ---------------------------------------------- #
    from scipy.spatial import cKDTree
    endo = beta | alpha
    tree = cKDTree(dd.xy[endo])
    n_near = tree.query_ball_point(dd.xy[endo], r=half_um * 0.6, return_length=True)
    centre = dd.xy[endo][int(np.argmax(n_near))]

    sel = (np.abs(dd.xy - centre) < half_um).all(axis=1)
    ids = dd.cell_id[sel]
    print(f"[fig] {donor}: window at ({centre[0]:.0f}, {centre[1]:.0f}) um, "
          f"{sel.sum():,} cells, {int((beta & sel).sum())} beta- / "
          f"{int((alpha & sel).sum())} alpha-committed", flush=True)

    lvl, factor = PR.pyramid_level_for(qptiff, PR.PROTEIN_UM_PER_PX, base_um)
    s = PR.PROTEIN_UM_PER_PX
    bundle = xenio._bundle_dir(donor, XENIUM_PATHS, roi)

    polys = {}
    for comp in ("cell", "nucleus"):
        g = PR.polygons_per_cell(bundle, comp, ids)
        polys[comp] = [None if x is None else shapely.transform(x, transform.apply) for x in g]

    def rings(geoms):
        out = []
        for g in geoms:
            if g is None:
                out.append(None); continue
            parts = list(g.geoms) if g.geom_type == "MultiPolygon" else [g]
            big = max(parts, key=lambda p: p.area)
            out.append(np.asarray(shapely.get_coordinates(big.exterior)) / s - 0.5)
        return out

    cell_rings, nuc_rings = rings(polys["cell"]), rings(polys["nucleus"])
    allpts = np.vstack([r for r in cell_rings if r is not None])
    x0, y0 = np.floor(allpts.min(0)).astype(int) - 4
    x1, y1 = np.ceil(allpts.max(0)).astype(int) + 4

    names = HZ.channel_names(qptiff)
    pages = {m: HZ.channel_index(qptiff, m, row) for m in ("INS", "GCG")}
    dapi = names.index("DAPI")

    # The display range comes from the per-CELL values over the WHOLE section, not from the
    # pixels and not from this window. Both halves matter.
    #
    # Not the pixels: 97% of a pancreas section is exocrine, so a p99.5 over pixels lands at 12
    # counts on 6442 while an islet cell reads 40-240. Every endocrine cell then saturates to
    # one colour and the within-islet comparison this sheet exists for is erased -- which is
    # what the first version of this figure did.
    #
    # Not this window: a range taken from what is drawn rescales each panel to its own contents
    # and makes a dead channel look exactly like a live one. vmin stays 0 in every panel, and 0
    # is where the scanner clamps, never a percentile black point.
    vals_all, vmaxes = {}, {}
    for m in ("INS", "GCG"):
        vals_all[m] = np.load(out_dir / "protein_rna" / f"{donor}_{serial}_cell_{m}.npy")
        vmaxes[m] = float(np.nanpercentile(vals_all[m], VALUE_PERCENTILE))
    vals = {m: vals_all[m][sel] for m in ("INS", "GCG")}

    crops = {}
    with tifffile.TiffFile(str(qptiff)) as tif:
        store = tif.series[0].aszarr(level=lvl)
        try:
            z = zarr.open(store, mode="r")
            for key, page in (("DAPI", dapi), ("INS", pages["INS"]), ("GCG", pages["GCG"])):
                crops[key] = np.asarray(z[page, y0:y1, x0:x1]).astype(np.float32)
            vmaxes["DAPI"] = DISP.unsigned_display_range(
                np.asarray(z[dapi, ::23, ::23]).ravel()).vmax
        finally:
            store.close()

    fig, axes = plt.subplots(2, 3, figsize=(16.2, 11.4))
    ext = [x0, x1, y1, y0]

    def img_panel(ax, key, title, outline=None, ec="#f76707"):
        ax.imshow(crops[key], cmap="gray", vmin=0.0, vmax=vmaxes[key],
                  interpolation="nearest", extent=ext)
        if outline is not None:
            ax.add_collection(PolyCollection([r for r in outline if r is not None],
                                             facecolors="none", edgecolors=ec, linewidths=0.6))
        ax.set_xlim(x0, x1); ax.set_ylim(y1, y0); ax.set_xticks([]); ax.set_yticks([])
        ax.set_title(title, fontsize=10)

    img_panel(axes[0, 0], "DAPI", "IF DAPI + Xenium NUCLEUS polygons\n"
                                  "(the registration, at 0.5 µm/px)", nuc_rings, "#f76707")
    img_panel(axes[0, 1], "INS", f"IF {names[pages['INS']]} — recorded INS\n"
                                 f"vmin 0, vmax {vmaxes['INS']:.0f}"
                                 f" (p{VALUE_PERCENTILE} of all cells in the section)")
    img_panel(axes[0, 2], "GCG", f"IF {names[pages['GCG']]} — recorded GCG\n"
                                 f"vmin 0, vmax {vmaxes['GCG']:.0f}"
                                 f" (p{VALUE_PERCENTILE} of all cells in the section)")

    def filled(ax, v, title, vmax):
        keep = [i for i, r in enumerate(cell_rings) if r is not None and np.isfinite(v[i])]
        pc = PolyCollection([cell_rings[i] for i in keep], array=np.asarray(v)[keep],
                            cmap="magma", edgecolors="#00000022", linewidths=0.3)
        pc.set_clim(0.0, vmax)
        ax.add_collection(pc)
        ax.set_facecolor("#f8f9fa")
        ax.set_xlim(x0, x1); ax.set_ylim(y1, y0); ax.set_xticks([]); ax.set_yticks([])
        ax.set_title(title, fontsize=10)
        fig.colorbar(pc, ax=ax, fraction=0.046, pad=0.02)

    filled(axes[1, 0], vals["INS"],
           "Xenium CELL polygons filled with the measured\n"
           f"{names[pages['INS']]} (recorded INS) — the number the stats used", vmaxes["INS"])
    filled(axes[1, 1], vals["GCG"],
           "Xenium CELL polygons filled with the measured\n"
           f"{names[pages['GCG']]} (recorded GCG)", vmaxes["GCG"])

    ax = axes[1, 2]
    cls = np.where(beta[sel] & ~alpha[sel], BETA_C,
                   np.where(alpha[sel] & ~beta[sel], ALPHA_C,
                            np.where(beta[sel] & alpha[sel], BOTH_C, NEITHER_C)))
    keep = [i for i, r in enumerate(cell_rings) if r is not None]
    ax.add_collection(PolyCollection([cell_rings[i] for i in keep],
                                     facecolors=cls[keep], edgecolors="#00000022",
                                     linewidths=0.3))
    ax.set_facecolor("#f8f9fa")
    ax.set_xlim(x0, x1); ax.set_ylim(y1, y0); ax.set_xticks([]); ax.set_yticks([])
    ax.set_title("the SAME cells, by their own Xenium RNA\n"
                 "blue = beta-committed, orange = alpha-committed,\n"
                 "purple = both, grey = neither (detection only)", fontsize=10)

    fig.suptitle(
        f"Donor {donor} · Xenium {serial} · protein and RNA on the same cell, by cell_id\n"
        f"Xenium segmentation warped onto the post-Xenium IF. No QuPath annotation is used; "
        f"the window was chosen by Xenium RNA alone.", fontsize=12)
    fig.tight_layout(rect=(0, 0, 1, 0.94))
    p = out_dir / "protein_rna" / f"protein_rna_{donor}_{serial}.png"
    fig.savefig(p, dpi=145)
    plt.close(fig)
    print(f"[fig] {p}", flush=True)
    return {"donor_id": donor, "serial": serial, "figure": str(p),
            "window_centre_um": [float(centre[0]), float(centre[1])],
            "n_cells_drawn": int(sel.sum()),
            "n_beta_committed": int((beta & sel).sum()),
            "n_alpha_committed": int((alpha & sel).sum())}


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--donors", nargs="*", default=None)
    ap.add_argument("--roi", default="panc")
    ap.add_argument("--geojson-dir", type=Path, default=GEOJSON_DIR)
    ap.add_argument("--out", type=Path, default=OUT)
    args = ap.parse_args(argv)

    sections = HZ.load_section_map(SECTIONS)
    recs = []
    for serial, _gj in discover(args.geojson_dir).items():
        donor = donor_of(serial, XENIUM_PATHS)
        if args.donors and donor not in args.donors:
            continue
        recs.append(figure_one(serial, sections, roi=args.roi, out_dir=args.out))
    if not recs:
        print("no donor matched", file=sys.stderr)
        return 1
    (args.out / "protein_rna" / "figures.json").write_text(json.dumps(recs, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
