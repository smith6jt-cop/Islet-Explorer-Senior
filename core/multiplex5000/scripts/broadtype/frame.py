"""The frozen screen frame: 16 hash-drawn 1024-px cores per donor, one mask-derived cell set.

    cd core/multiplex5000 && ../../.venv/bin/python scripts/broadtype/frame.py --donors 6436 6442

Per donor this streams the cell GeoJSON twice (once to rasterize the 2048-px windows around
each core, once inside ``verify_against_archive``) and writes

    data/development/broadtype/frame/regions.json                the draw, with window coords
    data/development/broadtype/frame/cells/donor_id=<d>.parquet  object_id, region_idx, draw, in_expression
    data/development/broadtype/frame/masks/donor_id=<d>/w<idx>.npy  int32 window labels (global)
    data/development/broadtype/frame/verify_against_archive.json  the bit-identity gate
    data/development/broadtype/frame/tails_in_frame.csv          Step-1 tail cells inside the frame

The islet-anchored draw uses QuPath parentage (``cell_region == 'core'``) as GEOMETRY only; it is
content-blind and cannot see any call.  Labels are global (``object_ids[label - 1]`` is the
UUID), so a cell in the frame joins the run's expression table on ``object_id``.  A donor whose
archive mask disagrees with the fresh raster stops the lane.
"""
from __future__ import annotations

import argparse
import json
import pathlib
import sys
from dataclasses import asdict

import numpy as np
import pandas as pd

HERE = pathlib.Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import pins  # noqa: E402
from criteria import CRITERIA  # noqa: E402
import regions as dct_regions  # noqa: E402  (scripts/dct/regions.py -- the hash rule)
import window_mask  # noqa: E402  (scripts/common)

FRAME = pins.OUT / "frame"


def image_shape(donor: str) -> tuple:
    p = pins.archive_mask(donor)
    if p.exists():
        import tifffile
        with tifffile.TiffFile(p) as tf:
            return tuple(int(v) for v in tf.pages[0].shape)
    from phenocycler.redsea import qptiff_channels
    _tiff, _pages, _names, shape, _px = qptiff_channels(pins.qptiff_path(donor), 0)
    return tuple(int(v) for v in shape)


def windows_for(regs: list) -> list:
    pad = (CRITERIA["window_px"] - CRITERIA["core_px"]) // 2
    return [window_mask.Window(int(r.y0) - pad, int(r.x0) - pad, int(CRITERIA["window_px"])) for r in regs]


def build_donor(donor: str, *, verify: bool = True) -> dict:
    shape = image_shape(donor)
    regs = dct_regions.select_regions(
        donor, pins.expression_parquet(donor), pins.PIXEL_SIZE_UM, shape,
        core_px=CRITERIA["core_px"], n_blind=CRITERIA["n_blind"], n_islet=CRITERIA["n_islet"],
        max_window=CRITERIA["window_px"])
    wins = windows_for(regs)
    masks, oids = window_mask.rasterize_windows(pins.geojson_path(donor), wins)
    oids = np.asarray(oids, dtype=object)
    pad = (CRITERIA["window_px"] - CRITERIA["core_px"]) // 2
    mdir = FRAME / "masks" / f"donor_id={donor}"
    mdir.mkdir(parents=True, exist_ok=True)
    rows = []
    for i, (m, r) in enumerate(zip(masks, regs)):
        np.save(mdir / f"w{i:02d}.npy", m)
        core = m[pad:pad + CRITERIA["core_px"], pad:pad + CRITERIA["core_px"]]
        labels = np.unique(core[core > 0])
        for lab in labels:
            rows.append((str(oids[int(lab) - 1]), i, r.draw))
    cells = pd.DataFrame(rows, columns=["object_id", "region_idx", "draw"])
    dup = cells["object_id"].duplicated(keep="first")
    if dup.any():
        # a cell straddling two adjacent cores is scored once, in the first core that holds it
        cells = cells.loc[~dup].reset_index(drop=True)
    expr_ids = pd.read_parquet(pins.expression_parquet(donor), columns=["object_id"])["object_id"].astype(str)
    cells["in_expression"] = cells["object_id"].isin(set(expr_ids))
    cdir = FRAME / "cells"
    cdir.mkdir(parents=True, exist_ok=True)
    cells.to_parquet(cdir / f"donor_id={donor}.parquet", index=False)

    ver = None
    if verify and pins.archive_mask(donor).exists():
        ok, per = window_mask.verify_against_archive(pins.geojson_path(donor), pins.archive_mask(donor), wins)
        ver = {"donor_id": donor, "archive": str(pins.archive_mask(donor)), "all_identical": bool(ok), "windows": per}
        if not ok:
            raise RuntimeError(f"donor {donor}: fresh raster differs from the archived mask -- the lane stops here")
    return {
        "donor_id": donor,
        "image_shape": list(shape),
        "regions": [{**asdict(r), "window": {"y0": w.y0, "x0": w.x0, "size": w.size}} for r, w in zip(regs, wins)],
        "n_cells_frame": int(len(cells)),
        "n_cells_in_expression": int(cells["in_expression"].sum()),
        "n_geojson_only": int((~cells["in_expression"]).sum()),
        "verify": ver,
    }


def tails_in_frame(donors: list) -> pd.DataFrame:
    """Per (donor, marker): Step-1 tail size k and how many tail cells the frame holds."""
    q, min_k = CRITERIA["tail_q"], CRITERIA["tail_min_k"]
    out = []
    for d in donors:
        frame_ids = set(pd.read_parquet(FRAME / "cells" / f"donor_id={d}.parquet")["object_id"].astype(str))
        tab = pd.read_parquet(pins.uniform_parquet(d))
        markers = [c for c in tab.columns if c not in ("donor_id", "object_id", "qc_analysis_eligible",
                                                      "X_centroid", "Y_centroid", "cell_region", "islet_num")]
        ids = tab["object_id"].astype(str).to_numpy()
        for m in markers:
            v = tab[m].to_numpy(float)
            finite = np.isfinite(v)
            n = int(finite.sum())
            if n == 0:
                out.append((d, m, 0, 0, 0, False)); continue
            k = max(min_k, int(q * n))
            vv = np.where(finite, v, -np.inf)
            idx = np.argpartition(-vv, k - 1)[:k]
            in_frame = int(sum(1 for i in idx if ids[i] in frame_ids))
            out.append((d, m, n, k, in_frame, in_frame >= CRITERIA["min_tail_in_pool"]))
    df = pd.DataFrame(out, columns=["donor_id", "marker", "n_finite", "k_tail", "n_tail_in_frame", "judgeable_in_frame"])
    df.to_csv(FRAME / "tails_in_frame.csv", index=False)
    return df


def main(argv=None) -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--donors", nargs="+", default=list(pins.DONORS))
    ap.add_argument("--no-verify", action="store_true")
    a = ap.parse_args(argv)
    pins.assert_pins()
    FRAME.mkdir(parents=True, exist_ok=True)
    summary = {"method_version": "broadtype-frame-2026-09-05", "criteria": CRITERIA,
               "rule": "dct regions.select_regions hash rule; 8 blind + 8 islet-anchored 1024-px cores; "
                       "one mask-derived cell set per donor for every arm", "donors": {}}
    verify_rows = []
    for d in a.donors:
        res = build_donor(d, verify=not a.no_verify)
        summary["donors"][d] = {k: v for k, v in res.items() if k != "verify"}
        if res["verify"] is not None:
            verify_rows.append(res["verify"])
        print(f"donor {d}: {len(res['regions'])} regions, {res['n_cells_frame']:,} frame cells "
              f"({res['n_cells_in_expression']:,} in expression; {res['n_geojson_only']:,} GeoJSON-only), "
              f"archive identical={res['verify']['all_identical'] if res['verify'] else 'n/a'}", flush=True)
    (FRAME / "regions.json").write_text(json.dumps(summary, indent=1))
    if verify_rows:
        (FRAME / "verify_against_archive.json").write_text(json.dumps(verify_rows, indent=1))
    t = tails_in_frame(a.donors)
    print(f"tails_in_frame: {int(t['judgeable_in_frame'].sum())}/{len(t)} (donor, marker) pairs judgeable in the frame")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
