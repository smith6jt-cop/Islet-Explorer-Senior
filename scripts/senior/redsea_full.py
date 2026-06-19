#!/usr/bin/env python3
"""
Scalable pixel-level REDSEA (lateral spillover correction) for the 59-plex CODEX /
PhenoCycler pancreas data (23.3M cells, 15 donor images).

Why: global clustering on this panel produces a large fake "Unassigned" fraction driven by
a ~22% "pan-positive bright" population = lateral signal bleeding across shared cell
boundaries. `redsea_approx.py` proved spillover is real (impossible co-positives collapse)
but the centroid-distance table approximation over-erodes true single-positives. REDSEA
(Bai et al. 2021) corrects it at the source using the segmentation boundary, but the
reference `redseapy` does not scale (dense cellNum^2 matrix, whole-image imread, pure-python
pixel loops). This reimplements the SAME math scalably:

  data[i,ch]    = whole-cell SUM of channel ch over cell i's pixels        (np.bincount)
  edge[i,ch]    = SUM over cell i's boundary band (within ~2px of a DIFFERENT label)
  contact[i,j]  = # of 8-connected pixel adjacencies between cells i and j  (sparse COO)
  F             = row-normalize(contact)            (F[i,j] = contact[i,j]/rowsum_i)
  corrected_sum = clip( data + edge - (F @ edge), 0 )      # subtract+reinforce, recipient form
  corrected_mean= corrected_sum / cell_area_px

Mask source: per-cell GeoJSON exported from QuPath (scripts/export_cells_geojson.groovy),
full-res pixel coords, feature id == object_id (UUID). We rasterize to an int32 instance
mask (cv2 cannot safely fill uint32; skimage.draw is pure-numpy). Read the qptiff ONE
CHANNEL AT A TIME (series 0 / level 0). Join corrected MFI back to object_id.

Run in `scvi-env`. Reads (read-only): the GeoJSON, the qptiff, data/cells/ (sanity only).
Writes: data/cells_redsea/donor_id=<id>/data_0.parquet (+ optional persisted mask).

  python scripts/senior/redsea_full.py --donor 6539                 # smoke test (smallest)
  python scripts/senior/redsea_full.py --donor 6539 --keep-mask     # also persist the mask
  python scripts/senior/redsea_full.py --all                        # every donor with a geojson
"""
from __future__ import annotations
import argparse
import glob
import json
import re
import time
from pathlib import Path

import numpy as np
import pandas as pd
import scipy.sparse as sp
import tifffile
from skimage.draw import polygon as draw_polygon
from skimage.segmentation import find_boundaries
from skimage.morphology import binary_dilation, diamond

BASE = Path(__file__).resolve().parents[2] / "data"
IMAGES_DIR = Path("/home/smith6jt/IO60panc2nd/Images")
SCRATCH = BASE / "redsea_scratch"
GEOJSON_DIR = SCRATCH / "geojson"
MASK_DIR = SCRATCH / "masks"
INTER_DIR = SCRATCH / "intermediates"
OUT_DIR = BASE / "cells_redsea"

SANITIZE = lambda s: re.sub(r"[ /\-]", "_", s)  # qptiff channel name -> parquet column name


def log(msg):
    print(msg, flush=True)


# --------------------------------------------------------------------------- discovery
def donor_image(donor: str) -> str:
    """The QuPath image string for this donor, e.g. '6539_Scan1.er.qptiff - resolution #1'."""
    f = sorted(glob.glob(str(BASE / "cells" / f"donor_id={donor}" / "*.parquet")))
    if not f:
        raise SystemExit(f"[err] no data/cells/donor_id={donor}")
    return str(pd.read_parquet(f[0], columns=["image"]).image.iloc[0])


def resolve_paths(donor: str):
    img = donor_image(donor)
    tag = re.sub(r"[^A-Za-z0-9._-]", "_", img)              # matches the Groovy sanitization
    geojson = GEOJSON_DIR / f"cells__{tag}.geojson"
    qptiff = IMAGES_DIR / img.split(" - ")[0].strip()       # '...qptiff - resolution #1' -> '...qptiff'
    if not geojson.exists():
        raise SystemExit(f"[err] missing geojson: {geojson}\n      run export_cells_geojson.groovy for {img}")
    if not qptiff.exists():
        raise SystemExit(f"[err] missing qptiff: {qptiff}")
    return img, geojson, qptiff


# --------------------------------------------------------------------------- qptiff
def qptiff_channels(qptiff: Path, level: int):
    """Return (level_pages, channel_names, (H,W), pixel_size_um) for the chosen pyramid level."""
    tf = tifffile.TiffFile(str(qptiff))
    series = tf.series[0]
    lv = series.levels[level]
    pages = lv.pages
    names = []
    for pg in pages:
        d = pg.description or ""
        m = re.search(r"<Biomarker>([^<]+)</Biomarker>", d) or re.search(r"<Name>([^<]+)</Name>", d)
        names.append(m.group(1).strip() if m else f"ch{len(names)}")
    H, W = lv.shape[-2], lv.shape[-1]
    # pixel size (um/px) from level-0 XResolution (pixels per ResolutionUnit; 3 == cm)
    px_um = None
    try:
        p0 = series.levels[0].pages[0]
        xr = p0.tags["XResolution"].value
        unit = p0.tags["ResolutionUnit"].value  # 3==cm, 2==inch
        ppu = xr[0] / xr[1]
        per_cm = ppu if int(unit) == 3 else ppu / 2.54
        px_um = 1e4 / per_cm * (2 ** level)
    except Exception:
        pass
    return tf, pages, names, (H, W), px_um


def read_channel(pages, ch, dtype=np.float32):
    """One full-res channel as a 2D array (float32 keeps bincount weights light)."""
    return np.asarray(pages[ch].asarray(), dtype=dtype)


# --------------------------------------------------------------------------- rasterize
def rasterize_mask(geojson: Path, shape, downsample: float):
    """Rasterize per-cell polygons into an int32 instance mask. label = feature_index+1.
    Returns (mask, object_ids[list]) where object_ids[label-1] is the cell's UUID."""
    log(f"  [raster] loading {geojson.name} ...")
    feats = json.load(open(geojson))["features"]
    n = len(feats)
    mask = np.zeros(shape, dtype=np.int32)
    oids = [None] * n
    H, W = shape
    s = 1.0 / downsample
    n_holes = 0
    t0 = time.time()
    for i, ft in enumerate(feats):
        oids[i] = ft.get("id")
        label = i + 1
        g = ft["geometry"]
        gtype = g["type"]
        polys = g["coordinates"] if gtype == "MultiPolygon" else [g["coordinates"]]
        for poly in polys:                       # poly = [exterior, hole1, ...]
            ext = np.asarray(poly[0], dtype=np.float64)
            if ext.shape[0] < 3:
                continue
            cc = ext[:, 0] * s                    # x -> col
            rr = ext[:, 1] * s                    # y -> row
            pr, pc = draw_polygon(rr, cc, shape=shape)
            mask[pr, pc] = label
            for hole in poly[1:]:                 # carve interior rings back to background
                h = np.asarray(hole, dtype=np.float64)
                if h.shape[0] < 3:
                    continue
                n_holes += 1
                hr, hc = draw_polygon(h[:, 1] * s, h[:, 0] * s, shape=shape)
                mask[hr, hc] = 0
        if (i + 1) % 100000 == 0:
            log(f"  [raster] {i+1:,}/{n:,} ({(time.time()-t0):.0f}s)")
    log(f"  [raster] {n:,} cells, {n_holes} holes carved, {time.time()-t0:.0f}s")
    return mask, oids


# --------------------------------------------------------------------------- contact matrix
def contact_matrix(mask: np.ndarray, n: int):
    """8-connected cell-cell adjacency counts as a symmetric sparse CSR (NxN)."""
    H, W = mask.shape
    keys = []
    # 4 direction vectors cover all 8 neighbours without double counting
    for dy, dx in [(0, 1), (1, 0), (1, 1), (1, -1)]:
        if dx >= 0:
            a = mask[: H - dy, : W - dx]
            b = mask[dy:, dx:]
        else:
            a = mask[: H - dy, -dx:]
            b = mask[dy:, : W + dx]
        v = (a > 0) & (b > 0) & (a != b)
        if not v.any():
            continue
        aa = a[v].astype(np.int64)
        bb = b[v].astype(np.int64)
        lo = np.minimum(aa, bb)
        hi = np.maximum(aa, bb)
        keys.append(lo * (np.int64(n) + 1) + hi)
    if not keys:
        return sp.csr_matrix((n, n), dtype=np.float64)
    keys = np.concatenate(keys)
    uk, cnt = np.unique(keys, return_counts=True)
    lo = (uk // (n + 1)).astype(np.int64)
    hi = (uk % (n + 1)).astype(np.int64)
    # to 0-based labels (labels are 1..n)
    M = sp.coo_matrix((cnt.astype(np.float64), (lo - 1, hi - 1)), shape=(n, n)).tocsr()
    M = M + M.T                                   # symmetric
    return M


def contact_matrix_gapbridge(mask, n, grow):
    """Like contact_matrix, but first fills thin (<=grow px) background gaps between cells
    with the nearest cell label, so neighbours separated by a 1-2 px gap still count as
    contacting (only ~23% of pixels are foreground here -- cells are NOT gapless, so direct
    adjacency alone leaves ~8% of cells with no neighbour). The grown mask is used ONLY for
    the contact graph; sums/edge use the original mask."""
    if grow <= 0:
        return contact_matrix(mask, n)
    from scipy.ndimage import distance_transform_edt
    bg = mask == 0
    dist, idx = distance_transform_edt(bg, return_distances=True, return_indices=True)
    grown = mask[tuple(idx)]
    grown[bg & (dist > grow + 0.5)] = 0           # only bridge gaps within `grow` px
    M = contact_matrix(grown, n)
    del grown, idx, dist
    return M


# --------------------------------------------------------------------------- compensate
def compensate(data, edge, sizes, M, comp_mode, alpha):
    """REDSEA compensation. corrected_sum = data + comp_mode*edge - alpha*(F @ edge), clamp>=0.
      comp_mode=1 : subtract + reinforce (the redseapy default; reinforces own boundary)
      comp_mode=0 : subtract only (remove incoming neighbour spillover; no reinforce)
    F is row-normalized contact (F[i,j]=contact[i,j]/rowsum_i); isolated cells -> F row 0
    (passthrough). Returns (corrected_mean[N,C], clamp_frac[C], n_isolated)."""
    rowsum = np.asarray(M.sum(1)).ravel()
    inv = np.zeros_like(rowsum)
    nz = rowsum > 0
    inv[nz] = 1.0 / rowsum[nz]
    F = sp.diags(inv) @ M                           # rows sum to 1 (recipient form); F @ edge, NOT F.T
    corrected_sum = data + comp_mode * edge - alpha * (F @ edge)
    np.clip(corrected_sum, 0, None, out=corrected_sum)
    with np.errstate(invalid="ignore", divide="ignore"):
        corrected = corrected_sum / sizes[:, None]
    corrected[sizes == 0, :] = np.nan
    clamp_frac = (corrected_sum == 0).mean(axis=0)
    return corrected, clamp_frac, int((~nz).sum())


def save_intermediates(donor, oids, cols, data, edge, sizes, M):
    INTER_DIR.mkdir(parents=True, exist_ok=True)
    p = INTER_DIR / f"{donor}.npz"
    M = M.tocsr()
    np.savez_compressed(p, oids=np.array(oids, dtype=object), cols=np.array(cols),
                        data=data.astype(np.float32), edge=edge.astype(np.float32),
                        sizes=sizes, m_data=M.data, m_indices=M.indices, m_indptr=M.indptr,
                        m_shape=np.array(M.shape))
    log(f"[{donor}] saved intermediates -> {p}")


def load_intermediates(donor):
    z = np.load(INTER_DIR / f"{donor}.npz", allow_pickle=True)
    M = sp.csr_matrix((z["m_data"], z["m_indices"], z["m_indptr"]), shape=tuple(z["m_shape"]))
    return (list(z["oids"]), list(z["cols"]), z["data"].astype(np.float64),
            z["edge"].astype(np.float64), z["sizes"].astype(np.float64), M)


def write_corrected(donor, oids, cols, corrected, sizes):
    out = pd.DataFrame({"object_id": oids, "donor_id": donor, "cell_area_px": sizes})
    for k, c in enumerate(cols):
        out[c] = corrected[:, k]
    out_dir = OUT_DIR / f"donor_id={donor}"
    out_dir.mkdir(parents=True, exist_ok=True)
    out.to_parquet(out_dir / "data_0.parquet", index=False)
    log(f"[{donor}] wrote {out_dir/'data_0.parquet'} ({len(out):,} rows)")


# --------------------------------------------------------------------------- main per donor
def process_donor(donor: str, args):
    t_all = time.time()
    img, geojson, qptiff = resolve_paths(donor)
    level = int(round(np.log2(args.downsample)))
    if 2 ** level != int(args.downsample):
        raise SystemExit(f"[err] --downsample must be a power of 2 (got {args.downsample})")
    log(f"[{donor}] image={img}")

    tf, pages, ch_names, (H, W), px_um = qptiff_channels(qptiff, level)
    cols = [SANITIZE(c) for c in ch_names]
    log(f"[{donor}] level={level} dims={H}x{W} channels={len(cols)} px_um={px_um}")

    # 1) rasterize instance mask
    mask, oids = rasterize_mask(geojson, (H, W), args.downsample)
    n = len(oids)
    sizes = np.bincount(mask.ravel(), minlength=n + 1)[1:].astype(np.float64)  # cell area in px
    n_empty = int((sizes == 0).sum())
    log(f"[{donor}] rasterized {n:,} cells | {n_empty} with 0 px | total fg px={int(sizes.sum()):,}")

    # 2) boundary band (within ~edge_radius px of a different label), 8-connected
    t0 = time.time()
    inner = find_boundaries(mask, mode="inner", connectivity=2, background=0)
    band = binary_dilation(inner, diamond(args.edge_radius)) & (mask > 0)
    band_idx = np.flatnonzero(band.ravel())
    mask_flat = mask.ravel()
    mask_band = mask_flat[band_idx]
    log(f"[{donor}] boundary band: {band_idx.size:,} px ({100*band_idx.size/mask.size:.1f}% of image), {time.time()-t0:.0f}s")
    del inner, band

    # 3) per-channel whole-cell sums + boundary sums (streamed)
    nch = len(cols)
    data = np.zeros((n, nch), dtype=np.float64)
    edge = np.zeros((n, nch), dtype=np.float64)
    t0 = time.time()
    for k in range(nch):
        ch = read_channel(pages, k).ravel()
        data[:, k] = np.bincount(mask_flat, weights=ch, minlength=n + 1)[1:]
        edge[:, k] = np.bincount(mask_band, weights=ch[band_idx], minlength=n + 1)[1:]
        del ch
        if (k + 1) % 10 == 0 or k == nch - 1:
            log(f"[{donor}] channels {k+1}/{nch} ({time.time()-t0:.0f}s)")
    tf.close()

    # 4) contact matrix (gap-bridged: cells here are not gapless)
    t0 = time.time()
    M = contact_matrix_gapbridge(mask, n, args.gap_bridge)
    log(f"[{donor}] contact matrix (gap_bridge={args.gap_bridge}): {M.nnz:,} contacts, {time.time()-t0:.0f}s")
    if args.save_intermediates:
        save_intermediates(donor, oids, cols, data, edge, sizes, M)

    # 5) compensate
    corrected, clamp_frac, n_isolated = compensate(data, edge, sizes, M, args.comp_mode, args.alpha)
    log(f"[{donor}] compensate(mode={args.comp_mode}, alpha={args.alpha}): {n_isolated} isolated cells")

    # ---- alignment sanity vs the existing parquet (correlation of UNcorrected mean) ----
    uncorr = np.where(sizes[:, None] > 0, data / np.where(sizes[:, None] == 0, 1, sizes[:, None]), np.nan)
    sanity = alignment_check(donor, oids, cols, uncorr, sizes, px_um)

    # 6) write corrected parquet (one row per object_id)
    write_corrected(donor, oids, cols, corrected, sizes)

    if args.keep_mask:
        MASK_DIR.mkdir(parents=True, exist_ok=True)
        tifffile.imwrite(MASK_DIR / f"{donor}.tif", mask, compression="zlib")
        log(f"[{donor}] persisted mask -> {MASK_DIR/f'{donor}.tif'}")
    del mask

    log(f"[{donor}] DONE in {time.time()-t_all:.0f}s | top-clamped: " +
        ", ".join(f"{cols[i]}={clamp_frac[i]*100:.1f}%" for i in np.argsort(clamp_frac)[-5:][::-1]))
    return sanity


def alignment_check(donor, oids, cols, uncorr_mean, sizes, px_um):
    """Correlate my UNcorrected full-res mean against the existing QuPath parquet means.
    A high correlation (~1) confirms mask<->image alignment before we trust the correction."""
    f = sorted(glob.glob(str(BASE / "cells" / f"donor_id={donor}" / "*.parquet")))[0]
    probe = [c for c in ["DAPI", "INS", "Pan_Cytokeratin", "CD3e", "Vimentin", "CD31"] if c in cols]
    ref = pd.read_parquet(f, columns=["object_id", "cell_area"] + probe)
    mine = pd.DataFrame({"object_id": oids, "_area_px": sizes})
    for c in probe:
        mine[c + "_mine"] = uncorr_mean[:, cols.index(c)]
    j = ref.merge(mine, on="object_id", how="inner")
    log(f"[{donor}] alignment vs parquet ({len(j):,}/{len(ref):,} matched on object_id):")
    out = {}
    for c in probe:
        a = j[c].to_numpy(np.float64)
        b = j[c + "_mine"].to_numpy(np.float64)
        ok = np.isfinite(a) & np.isfinite(b)
        r = np.corrcoef(a[ok], b[ok])[0, 1] if ok.sum() > 10 else np.nan
        ratio = np.nanmedian(b[ok] / np.where(a[ok] == 0, np.nan, a[ok]))
        out[c] = (r, ratio)
        log(f"    {c:16s} pearson_r={r:.3f}  median(mine/parquet)={ratio:.3f}")
    if px_um:
        ar = j["_area_px"].to_numpy(np.float64) * (px_um ** 2)
        rr = np.nanmedian(ar / j["cell_area"].to_numpy(np.float64))
        log(f"    cell_area: median(px_area*px_um^2 / parquet_um2)={rr:.3f}  (px_um={px_um:.4f})")
    return out


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--donor", default=None)
    ap.add_argument("--all", action="store_true")
    ap.add_argument("--downsample", type=float, default=1.0, help="power of 2; 1=full res (level 0)")
    ap.add_argument("--edge-radius", type=int, default=0,
                    help="boundary band dilation radius (px). 0 = the 1-px cell rim (any label "
                         "transition). A THIN band is essential: a thick band over-erodes here "
                         "because cells are small. (Original redseapy element_size~2 assumes big MIBI cells.)")
    ap.add_argument("--comp-mode", type=int, default=0, choices=[0, 1],
                    help="0=subtract only (default). 1=subtract+reinforce is the redseapy default "
                         "but is a structural NO-OP in this confluent tissue: directly-touching "
                         "cells share an ambiguous boundary, so reinforce (+edge) cancels the "
                         "subtraction (-F@edge). Reinforce only works for background-separated MIBI cells.")
    ap.add_argument("--alpha", type=float, default=1.0,
                    help="spillover subtraction strength (1.0 = natural full subtraction at the "
                         "1-px boundary: collapses impossible co-positives ~83%% while preserving "
                         "single-lineage cells; see redsea_full_validate.py)")
    ap.add_argument("--gap-bridge", type=int, default=1,
                    help="bridge background gaps <= N px when building the contact graph "
                         "(cells are not gapless here; 0 disables)")
    ap.add_argument("--keep-mask", action="store_true", help="persist the int32 instance mask")
    ap.add_argument("--save-intermediates", action="store_true",
                    help="dump data/edge/sizes/contact to scratch for fast --from-intermediates re-runs")
    ap.add_argument("--from-intermediates", action="store_true",
                    help="skip raster+channels; reload saved intermediates and just (re)compensate")
    a = ap.parse_args(argv)

    if a.all:
        donors = sorted(p.name.split("=")[1] for p in (BASE / "cells").glob("donor_id=*"))
    elif a.donor:
        donors = [a.donor]
    else:
        raise SystemExit("specify --donor <id> or --all")

    for d in donors:
        try:
            if a.from_intermediates:
                oids, cols, data, edge, sizes, M = load_intermediates(d)
                if a.gap_bridge > 0 and (MASK_DIR / f"{d}.tif").exists():
                    mask = tifffile.imread(MASK_DIR / f"{d}.tif")
                    M = contact_matrix_gapbridge(mask, len(oids), a.gap_bridge)
                    del mask
                    log(f"[{d}] recomputed gap-bridged contacts: {M.nnz:,}")
                corrected, clamp_frac, n_iso = compensate(data, edge, sizes, M, a.comp_mode, a.alpha)
                log(f"[{d}] recompensate(mode={a.comp_mode}, alpha={a.alpha}): {n_iso} isolated; "
                    "top-clamped: " + ", ".join(f"{cols[i]}={clamp_frac[i]*100:.1f}%"
                                                 for i in np.argsort(clamp_frac)[-5:][::-1]))
                write_corrected(d, oids, cols, corrected, sizes)
            else:
                process_donor(d, a)
        except SystemExit as e:
            log(f"[{d}] SKIP: {e}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
