#!/usr/bin/env python
"""Nucleus centroids on both sides of the re-stain, in each image's own micron frame.

    cd core/phenocycler
    ../../.venv/bin/python scripts/restain/nuclei.py --donors 6436 6442

The IF re-stained the same physical section Xenium measured, so the two DAPI images show the
SAME nuclei. That makes the registration a point-set problem rather than an intensity one, and
a point-set problem is the easier of the two: it needs no shared intensity model, it converges
by construction once the seed is inside half the nearest-neighbour distance (~5 um here), and
its residual is a per-cell distance rather than a per-window statistic.

Xenium's side needs no detection -- ``nucleus_boundaries.parquet`` is a segmentation of Xenium's
own DAPI, so its polygon centroids ARE the nuclei. This detects the IF side.

**Detection, not segmentation.** A band-pass at the nucleus scale, local maxima at a spacing no
smaller than a nucleus, then a parabolic sub-pixel fit. No cell objects are created and nothing
here is a cell type: these points exist only to be matched.

**The micron frame is taken from ``export.json``, never assumed.** The two pyramids are built
differently -- the Xenium morphology is decimated and the IF qptiff is box-averaged (measured;
see ``export_dapi.py``) -- so ``micron = pixel * um_per_px`` is wrong for both, by different
amounts. ``origin_um`` is the micron position of pixel coordinate 0 and it is not optional: a
constant offset is exactly the error this lane is trying to remove, so introducing one here
would be self-defeating.
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

import geo                                                        # noqa: E402
import tiles as TL                                                # noqa: E402
from run_register import GEOJSON_DIR, OUT, XENIUM_PATHS, discover, donor_of   # noqa: E402

#: A pancreatic nucleus is ~6 um across. The band-pass keeps that scale and rejects both the
#: pixel noise below it and the tissue-level brightness gradient above it.
NUCLEUS_UM = 6.0
#: Peaks closer than this are one nucleus counted twice. Slightly under a nucleus diameter so
#: two touching nuclei can still both be found.
MIN_SEPARATION_UM = 4.0
#: Robust z of the band-pass response, over the whole image, below which a peak is noise.
MIN_RESPONSE_Z = 3.0


class NucleiError(RuntimeError):
    """Nuclei cannot be detected as the geometry the caller claims."""


def band_pass(image: np.ndarray, um_per_px: float, *, nucleus_um: float = NUCLEUS_UM):
    """Difference of Gaussians centred on the nucleus scale, as float32.

    ``low`` blurs away the pixel noise without touching a nucleus; ``high`` blurs a nucleus
    away entirely, so the difference is what is nucleus-sized. Working on the difference rather
    than the raw DAPI is what makes one threshold valid across a section whose staining and
    tissue density vary several-fold.
    """
    from scipy.ndimage import gaussian_filter

    img = np.asarray(image, dtype=np.float32)
    low = max(0.8, 0.25 * nucleus_um / um_per_px)
    high = max(low * 2.0, 0.75 * nucleus_um / um_per_px)
    return gaussian_filter(img, low) - gaussian_filter(img, high)


def _robust_scale(a: np.ndarray, *, stride: int = 17) -> tuple[float, float]:
    """(median, sigma) of ``a`` from a strided subsample, sigma via the MAD.

    Strided rather than random so it is reproducible without a seed, and the MAD rather than
    the standard deviation because the nuclei themselves are the outliers this is meant to
    ignore -- an SD computed over a nucleus-rich field would set the threshold above them.
    """
    s = a[::stride, ::stride].ravel()
    s = s[np.isfinite(s)]
    if s.size == 0:
        raise NucleiError("band-pass response is entirely non-finite")
    med = float(np.median(s))
    sigma = float(1.4826 * np.median(np.abs(s - med)))
    if sigma <= 0.0:
        # More than half the sampled pixels share one value: a background with no variation at
        # all, which a real DAPI plane does not have but a synthetic one does. The SD is then
        # the only scale available; it is inflated by the nuclei themselves, so this is a
        # conservative fallback rather than a better estimator.
        sigma = float(s.std())
    return med, sigma


def local_maxima(a: np.ndarray, separation_px: int, threshold: float) -> tuple[np.ndarray, np.ndarray]:
    """Rows and columns of the strict local maxima of ``a`` at least ``separation_px`` apart.

    A maximum filter and an equality test, not ``skimage.feature.peak_local_max``. The two find
    the same peaks -- ``peak_local_max`` runs exactly this test first -- but it then calls
    ``ensure_spacing`` on every candidate whenever ``min_distance > 1``, which builds repeated
    KD-trees over the full candidate set. On a 180 megapixel section with ~10^6 nuclei that step
    dominates the run by orders of magnitude and enforces a spacing the maximum filter has
    already imposed.

    The border is dropped because the caller fits a parabola through each peak's neighbours.
    """
    from scipy.ndimage import maximum_filter

    mx = maximum_filter(a, size=2 * separation_px + 1, mode="nearest")
    peak = (a == mx) & (a >= threshold)
    peak[0, :] = peak[-1, :] = peak[:, 0] = peak[:, -1] = False
    return np.nonzero(peak)


def _parabolic_shift(vm1: np.ndarray, v0: np.ndarray, vp1: np.ndarray) -> np.ndarray:
    """Sub-sample peak offset in [-1, 1] from three samples; 0 where the fit is degenerate."""
    denom = vm1 - 2.0 * v0 + vp1
    out = np.zeros_like(v0, dtype=np.float64)
    ok = np.isfinite(denom) & (denom != 0.0)
    out[ok] = 0.5 * (vm1[ok] - vp1[ok]) / denom[ok]
    out[np.abs(out) > 1.0] = 0.0
    return out


def detect(image: np.ndarray, um_per_px: float, origin_um: float, *,
           nucleus_um: float = NUCLEUS_UM,
           min_separation_um: float = MIN_SEPARATION_UM,
           min_response_z: float = MIN_RESPONSE_Z,
           valid: np.ndarray | None = None) -> tuple[np.ndarray, np.ndarray]:
    """DAPI plane -> ``(xy_um, response_z)``, centroids in the image's own micron frame.

    ``origin_um`` is where pixel coordinate 0 sits in microns and comes from ``export.json``.
    Sub-pixel refinement is a parabolic fit on the band-pass response, which matters: at 1 um/px
    an integer peak carries +/-0.5 um of quantisation, and while that averages out of a fitted
    transform it does NOT average out of the per-cell distance this lane has to report.
    """
    dog = band_pass(image, um_per_px, nucleus_um=nucleus_um)
    med, sigma = _robust_scale(dog)
    if not sigma > 0:
        raise NucleiError("band-pass response has no spread -- image blank or constant?")
    thresh = med + min_response_z * sigma

    sep = max(1, int(round(min_separation_um / um_per_px)))
    r, c = local_maxima(dog, sep, thresh)
    if len(r) == 0:
        return np.empty((0, 2)), np.empty(0)

    if valid is not None:
        keep = valid[r, c]
        r, c = r[keep], c[keep]
        if len(r) == 0:
            return np.empty((0, 2)), np.empty(0)

    v0 = dog[r, c].astype(np.float64)
    dr = _parabolic_shift(dog[r - 1, c].astype(np.float64), v0, dog[r + 1, c].astype(np.float64))
    dc = _parabolic_shift(dog[r, c - 1].astype(np.float64), v0, dog[r, c + 1].astype(np.float64))

    xy = np.column_stack([(c + dc) * um_per_px + origin_um,
                          (r + dr) * um_per_px + origin_um])
    return xy, (v0 - med) / sigma


def xenium_nucleus_centroids(bundle: Path) -> tuple[np.ndarray, np.ndarray]:
    """``nucleus_boundaries.parquet`` -> ``(xy_um, cell_id)``, one row per NUCLEUS.

    Already in microns and already a segmentation of the DAPI these points are matched against,
    so there is nothing to detect and nothing to scale. A cell with two nuclei contributes two
    points and a cell with none contributes none, which is what makes this the nucleus field
    rather than the cell field.
    """
    import shapely

    from phenocycler.integration.boundaries import _read_boundary_rings

    p = Path(bundle) / "nucleus_boundaries.parquet"
    if not p.exists():
        raise NucleiError(f"{p} is absent -- it ships in the Xenium output bundle")
    cell_id, _labels, geoms = _read_boundary_rings(p)
    cen = shapely.centroid(geoms)
    return np.column_stack([shapely.get_x(cen), shapely.get_y(cen)]), np.asarray(cell_id)


def if_plane(meta: dict, target_um: float | None):
    """The IF DAPI plane to detect on, with its micron frame.

    ``None`` uses the plane ``export_dapi.py`` already wrote (1 um/px). A finer ``target_um``
    goes back to the qptiff for a lower pyramid level, deriving the factor and the origin by the
    SAME rules ``export_dapi`` used -- the point of the option is to find out whether the
    residual is limited by how precisely a 1 um/px peak can be located, and that answer is only
    meaningful if nothing else changes with it.
    """
    import tifffile

    import export_dapi as ED
    import hormone as HZ
    from run_hormone import IMAGES, SECTIONS

    qptiff = Path(meta["qptiff"])
    row = HZ.load_section_map(SECTIONS)
    row = row[row.image_file == qptiff.name].iloc[0]
    pixel_um = float(row.pixel_size_um)
    with tifffile.TiffFile(str(qptiff)) as tif:
        series = tif.series[0]
        levels = list(getattr(series, "levels", [series]))
        shapes = [(lv.shape[-2], lv.shape[-1]) for lv in levels]
        lvl = TL.level_for(shapes, pixel_um, target_um)
        f = ED.pyramid_factor(shapes, lvl)
        img = levels[lvl].pages[HZ.channel_index(qptiff, "NUCLEUS", row)].asarray()
    return img, pixel_um * f, ED.origin_um(pixel_um, f, averaged=True), int(lvl)


def detect_one(donor_dir: Path, *, gj: Path | None = None,
               min_response_z: float = MIN_RESPONSE_Z,
               if_target_um: float | None = None) -> dict:
    import tifffile

    donor_dir = Path(donor_dir)
    meta = json.loads((donor_dir / "export.json").read_text())
    for key in ("if_origin_um", "xenium_origin_um"):
        if key not in meta:
            raise NucleiError(
                f"{donor_dir/'export.json'} predates the pixel-origin fix (no {key}) -- "
                f"re-run export_dapi.py before detecting nuclei")

    if if_target_um is None:
        if_img = tifffile.imread(str(donor_dir / "slides" / "02_if_dapi.ome.tiff"))
        if_um, if_origin, lvl = (float(meta["if_um_per_px"]), float(meta["if_origin_um"]),
                                 int(meta.get("if", {}).get("level", -1)))
    else:
        if_img, if_um, if_origin, lvl = if_plane(meta, if_target_um)
    xy_if, z = detect(if_img, if_um, if_origin, min_response_z=min_response_z)
    del if_img

    rec = {"donor_id": meta["donor_id"], "serial": meta["serial"], "roi": meta["roi"],
           "if": {"um_per_px": if_um, "origin_um": if_origin, "level": lvl,
                  "n_detected": int(len(xy_if)),
                  "min_response_z": float(min_response_z),
                  "nucleus_um": NUCLEUS_UM, "min_separation_um": MIN_SEPARATION_UM}}

    if gj is not None:
        # Keep only what the maintainer drew as tissue. The IF scans the whole slide, which
        # carries a pancreatic lymph node with no counterpart in the __Panc__ Xenium bundle
        # (xenium_paths.csv rows 33-34); those nuclei can only ever be unmatched.
        import shapely
        tissue = geo.tissue_polygon(geo.load_annotations(gj))
        pts = shapely.points(xy_if[:, 0], xy_if[:, 1])
        shapely.prepare(tissue)                  # prepares in place and returns None
        inside = shapely.contains(tissue, pts)
        rec["if"]["n_in_tissue"] = int(inside.sum())
        xy_if, z = xy_if[inside], z[inside]

    np.save(donor_dir / "if_nuclei_um.npy", xy_if)
    np.save(donor_dir / "if_nuclei_z.npy", z.astype(np.float32))

    bundle = Path(meta["xenium_bundle"])
    xy_xen, cid = xenium_nucleus_centroids(bundle)
    np.save(donor_dir / "xenium_nuclei_um.npy", xy_xen)
    np.save(donor_dir / "xenium_nuclei_cell_ids.npy", cid)
    rec["xenium"] = {"n_nuclei": int(len(xy_xen)), "source": str(bundle / "nucleus_boundaries.parquet")}

    n_if = rec["if"].get("n_in_tissue", rec["if"]["n_detected"])
    rec["density_ratio_if_over_xenium"] = float(n_if / max(1, len(xy_xen)))
    print(f"[{rec['donor_id']}] IF nuclei {n_if:,} (of {rec['if']['n_detected']:,} detected) "
          f"vs Xenium nuclei {len(xy_xen):,}  ratio {rec['density_ratio_if_over_xenium']:.3f}",
          flush=True)
    return rec


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--donors", nargs="*", default=None)
    ap.add_argument("--roi", default="panc")
    ap.add_argument("--geojson-dir", type=Path, default=GEOJSON_DIR)
    ap.add_argument("--out", type=Path, default=OUT)
    ap.add_argument("--min-response-z", type=float, default=MIN_RESPONSE_Z)
    ap.add_argument("--if-target-um", type=float, default=None,
                    help="detect on the qptiff pyramid level nearest this instead of the "
                         "exported 1 um/px plane (e.g. 0.5), to test whether the residual is "
                         "limited by peak localisation")
    ap.add_argument("--whole-slide", action="store_true",
                    help="do not restrict IF detections to the drawn Tissue polygon")
    args = ap.parse_args(argv)

    rows = []
    for serial, gj in sorted(discover(args.geojson_dir).items()):
        donor = donor_of(serial, XENIUM_PATHS)
        if args.donors and donor not in set(args.donors):
            continue
        d = args.out / "valis" / f"donor_id={donor}"
        if not (d / "export.json").exists():
            print(f"[{donor}] no export.json -- run export_dapi.py first", file=sys.stderr)
            continue
        rows.append(detect_one(d, gj=None if args.whole_slide else gj,
                               min_response_z=args.min_response_z,
                               if_target_um=args.if_target_um))
    if not rows:
        print("no donor matched", file=sys.stderr)
        return 1
    tag = "" if args.if_target_um is None else f"_{args.if_target_um:g}um"
    p = args.out / f"nuclei_{args.roi}{tag}.json"
    p.write_text(json.dumps(rows, indent=2))
    print(f"-> {p}", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
