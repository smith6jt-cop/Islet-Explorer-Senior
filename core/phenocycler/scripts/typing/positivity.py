#!/usr/bin/env python3
"""Whole-slide per-marker positivity from Nimbus, tile by tile.

RETIRED -- NOT A FALLBACK.  Recorded here 2026-09-18 because the docstring below reads as a live
specification and is not one.  This module belongs to the Nimbus/Leiden chain
(``substrate.py`` -> ``positivity.py`` -> ``call.py`` -> ``validate.py`` -> ``export_qupath.py``,
with ``metacluster.py``) that the typing README used to describe.  It did not produce the calls on
disk.  The current typing is the subpanel classification in ``subpanel.py`` / ``build_cells.py``;
see ``README.md``, "Retired, and not to be revived as fallbacks".  What follows is kept as a record
of what was measured and why, not as a lane to fall back to.

Why a model and not a threshold: on this panel the per-cell mean intensity has no separable
positive population for most markers -- the distributions are unimodal in log space, and the
repo's own earlier work measured the same thing.  The information that distinguishes real
staining from a neighbour's bleed is *spatial* (is the signal on this cell's membrane?), and a
per-cell mean has thrown it away.  Nimbus reads the image around each cell and returns a
per-marker confidence, so no operating point has to be chosen at all: the confidences are used
as features downstream, never thresholded.

Validated before this was run, on one islet-centred window of donor 6442 against QuPath islet
annotations the model never saw: GCG AUC 0.996, INS 0.991, and Ker8_18 correctly *depleted*
inside islets at 0.130.

Runs the model in ``external/.venv-nimbus`` through the existing worker contract; this driver
stays in the repo ``.venv`` and exchanges files, because the two environments pin incompatible
torch/numpy stacks.
"""

from __future__ import annotations

import argparse
import json
import os
import pathlib
import shutil
import subprocess
from concurrent.futures import ThreadPoolExecutor
import sys
import time

import numpy as np
import pandas as pd
import tifffile

REPO = pathlib.Path(__file__).resolve().parents[4]
sys.path.insert(0, str(REPO / "core/phenocycler/scripts/common"))
sys.path.insert(0, str(REPO / "core/phenocycler/scripts/typing"))
from tiff_window import read_window  # noqa: E402
from window_mask import Window, rasterize_windows  # noqa: E402
import lineages as L  # noqa: E402

NIMBUS_PY = REPO / "external/.venv-nimbus/bin/python"
WORKER = REPO / "core/phenocycler/scripts/broadtype/arms/nimbus_worker.py"
OUT = REPO / "data/development/typing/positivity"
#: Pixel scratch goes to /data: tile TIFFs are ~8 MB per channel and /home has no room.
SCRATCH = pathlib.Path("/data/typing_nimbus")
TILE = 2048
#: A cell is ~9 px in radius, so only cells within a few px of a tile edge are truncated.  They
#: are flagged rather than dropped, so any statistic can state whether it excluded them.
EDGE_PX = 16

PARAMS = {"dataset_magnification": 20, "model_magnification": 10, "normalization_quantile": 0.999,
          "test_time_aug": True, "batch_size": 4, "input_shape": [1024, 1024], "device": "auto"}
PIN = {"checkpoint": "V1.pt",
       "checkpoint_sha256": "0651a3852fbafc03be415a1ec80ef947cc1913a081cf0bc649301ade887a82a9",
       "nimbus_version": "0.0.5", "torch": "2.5.1+cu121", "numpy": "1.26.4",
       "nimbus_commit": "241281ef85ba34d9857cc03886cc45aedd0f0b34"}


def manifest_image(donor: str) -> dict:
    man = json.loads((REPO / "core/phenocycler/data/qupath_manifest.json").read_text())
    hits = [i for i in man["images"] if str(i.get("donor_id")) == str(donor)]
    if not hits:
        raise KeyError(f"donor {donor} not in the QuPath manifest")
    return hits[0]


def _path(entry) -> pathlib.Path:
    return pathlib.Path(entry["path"] if isinstance(entry, dict) else entry)


def populated_tiles(donor: str, px_um: float) -> list[tuple[int, int]]:
    """Tiles holding at least one eligible cell centroid.  Roughly 40% of the grid is empty."""
    sub = pd.read_parquet(OUT.parent / f"substrate/donor_id={donor}/data.parquet",
                          columns=["eligible", "X_centroid", "Y_centroid"])
    sub = sub[sub["eligible"]]
    x = (sub["X_centroid"].to_numpy() / px_um).astype(np.int64)
    y = (sub["Y_centroid"].to_numpy() / px_um).astype(np.int64)
    return sorted({(int(a), int(b)) for a, b in zip(x // TILE, y // TILE)})


def run_batch(job_dir: pathlib.Path, fovs: list[dict], channels: list[str]) -> pd.DataFrame:
    job = {"fovs": fovs, "channels": channels, "output_dir": str(job_dir / "nimbus_out"),
           **PARAMS, "pin": PIN}
    (job_dir / "job.json").write_text(json.dumps(job, indent=1))
    r = subprocess.run([str(NIMBUS_PY), str(WORKER), str(job_dir / "job.json")],
                       capture_output=True, text=True)
    if r.returncode != 0:
        raise RuntimeError(f"nimbus worker failed:\n{r.stdout[-2000:]}\n{r.stderr[-3000:]}")
    return pd.read_parquet(job_dir / "nimbus_out" / "cell_table.parquet")


def run(donor: str, batch: int, limit: int | None, workers: int) -> pathlib.Path:
    img = manifest_image(donor)
    px_um = float(img["pixel_size_um_x"])
    qptiff, geojson = _path(img["qptiff"]), _path(img["cell_geojson"])
    by_marker = {c["marker"]: int(c["channel_index"]) for c in img["channel_map"] if c.get("marker")}
    markers = [m for m in L.identity_markers() if m in by_marker]

    tiles = populated_tiles(donor, px_um)
    if limit:
        tiles = tiles[:limit]
    windows = [Window(ty * TILE, tx * TILE, TILE) for tx, ty in tiles]
    print(f"[{donor}] {len(tiles)} populated tiles, {len(markers)} channels", flush=True)

    # Scratch is keyed on the PROCESS, not just the donor.  Keyed on the donor alone, a second
    # run of the same donor deletes the first run's FOV tiles while it is still reading them:
    # measured, two concurrent runs produced "FileNotFoundError: M2Gal3.tiff" and a worker crash,
    # and the timings taken from them were meaningless.  Checkpoints stay under OUT and are shared
    # deliberately, so a restart still resumes; only the pixel scratch is per-process.
    scratch = SCRATCH / f"{donor}__pid{os.getpid()}"
    if scratch.exists():
        shutil.rmtree(scratch)
    scratch.mkdir(parents=True)

    t0 = time.time()
    masks, oids = rasterize_windows(geojson, windows)
    oids = np.asarray([str(o) for o in oids], dtype=object)
    print(f"[{donor}] rasterized {len(windows)} tiles in {time.time() - t0:.0f}s "
          f"({len(oids):,} features)", flush=True)

    # Checkpoint each batch.  The first run of this driver lost 35 minutes of GPU when it died on
    # tile 36 of 197 with everything still in memory; a whole-slide pass is hours, so a crash must
    # cost one batch and not the donor.
    parts = OUT / f"donor_id={donor}" / "parts"
    parts.mkdir(parents=True, exist_ok=True)

    # CONCURRENCY.  The batch wall time (529-695 s per 12 tiles) is almost exactly the inference
    # time measured alone (52 s/tile), so the GPU is the whole cost -- and it sat at 5-23%
    # utilisation while doing it.  Nimbus runs eight test-time-augmentation passes per channel,
    # one channel at a time, on a small model, so the card idles between kernel launches: the job
    # is latency-bound, not throughput-bound.  Several workers on one card therefore overlap
    # almost perfectly.  Each thread opens its OWN TiffFile handle -- a shared handle is not safe
    # across concurrent reads -- and the subprocess call releases the GIL, so threads are enough.
    done_parts, todo = [], []
    for start in range(0, len(windows), batch):
        part = parts / f"batch{start:05d}.parquet"
        (done_parts if part.exists() else todo).append((start, part))
    for start, part in done_parts:
        print(f"[{donor}] batch{start:05d} already done, reusing", flush=True)

    n_done = [len(done_parts) * batch]

    def do_batch(item):
        start, part = item
        grp = list(range(start, min(start + batch, len(windows))))
        job_dir = scratch / f"batch{start:05d}"
        fovs = []
        with tifffile.TiffFile(qptiff) as tf_local:
            pages = tf_local.series[0].pages
            for i in grp:
                w = windows[i]
                fov = job_dir / "fovs" / f"t{i:05d}"
                fov.mkdir(parents=True, exist_ok=True)
                for m in markers:
                    arr = read_window(pages[by_marker[m]], w.y0, w.y1, w.x0, w.x1)
                    # Tiles on the slide's right and bottom edges are clipped by the image, while
                    # the mask is always a full tile, and Nimbus refuses the shape mismatch.
                    # Pad with zeros: beyond the slide there is no image and no polygon, so the
                    # mask is already 0 there and label alignment is preserved exactly.
                    if arr.shape != (TILE, TILE):
                        pad = np.zeros((TILE, TILE), dtype=arr.dtype)
                        pad[:arr.shape[0], :arr.shape[1]] = arr
                        arr = pad
                    tifffile.imwrite(fov / f"{m}.tiff", np.ascontiguousarray(arr),
                                     photometric="minisblack", compression="zlib")
                mp = job_dir / "masks" / f"t{i:05d}.tiff"
                mp.parent.mkdir(parents=True, exist_ok=True)
                tifffile.imwrite(mp, np.ascontiguousarray(masks[i].astype(np.int32)),
                                 photometric="minisblack", compression="zlib")
                fovs.append({"name": f"t{i:05d}", "path": str(fov), "mask": str(mp),
                             "window": {"y0": w.y0, "x0": w.x0, "size": TILE},
                             "n_labels": int(np.unique(masks[i][masks[i] > 0]).size)})
        tab = run_batch(job_dir, fovs, markers)
        tab["tile"] = tab["fov"].str.slice(1).astype(int)
        tab.to_parquet(part, index=False)
        shutil.rmtree(job_dir)                          # keep /data bounded
        n_done[0] += len(grp)
        print(f"[{donor}] {n_done[0]}/{len(windows)} tiles  {time.time() - t0:.0f}s elapsed",
              flush=True)
        return tab

    frames = [pd.read_parquet(p_) for _, p_ in done_parts]
    if todo:
        with ThreadPoolExecutor(max_workers=workers) as pool:
            frames.extend(pool.map(do_batch, todo))

    tab = pd.concat(frames, ignore_index=True)
    tab["object_id"] = oids[tab["label"].to_numpy().astype(np.int64) - 1]

    # A cell straddling a tile boundary has its pixels drawn into BOTH tiles, so it comes back as
    # two rows sharing one object_id, each scored on a fragment of the cell.  Left in, that
    # duplicates the cell and fans out every downstream join on object_id.  Keep the tile holding
    # the larger fragment -- the row scored on the more complete view of the cell -- and mark the
    # cell as split, so any statistic can state whether it excluded them.  Cells are ~9 px in
    # radius against a 2048 px tile, so this is a small population, but it is not empty.
    split_ids = set(tab.loc[tab.duplicated("object_id", keep=False), "object_id"])
    tab = (tab.sort_values("area", ascending=False)
              .drop_duplicates("object_id", keep="first")
              .reset_index(drop=True))
    tab["edge_cell"] = tab["object_id"].isin(split_ids)
    print(f"[{donor}] {len(split_ids):,} cells straddled a tile boundary; kept the larger "
          f"fragment of each", flush=True)

    dest = OUT / f"donor_id={donor}"
    dest.mkdir(parents=True, exist_ok=True)
    keep = ["object_id", "tile", "label", "area", "edge_cell", *markers]
    tab[keep].to_parquet(dest / "data.parquet", index=False)
    (dest / "params.json").write_text(json.dumps(
        {"donor_id": donor, "tile_px": TILE, "n_tiles": len(tiles), "markers": markers,
         "pixel_size_um": px_um, "params": PARAMS, "pin": PIN,
         "seconds": round(time.time() - t0, 1), "n_rows": int(len(tab))}, indent=2))
    shutil.rmtree(scratch, ignore_errors=True)      # pixel scratch; checkpoints stay
    print(f"[{donor}] wrote {dest/'data.parquet'}  rows={len(tab):,} "
          f"in {(time.time()-t0)/60:.1f} min", flush=True)
    return dest / "data.parquet"


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--donors", nargs="+", required=True)
    ap.add_argument("--batch", type=int, default=12, help="tiles per Nimbus invocation")
    ap.add_argument("--limit", type=int, default=None, help="first N tiles only (smoke test)")
    ap.add_argument("--workers", type=int, default=4,
                    help="concurrent Nimbus workers on the one GPU; it is latency-bound, so "
                         "several overlap almost perfectly")
    a = ap.parse_args()
    for d in a.donors:
        run(d, a.batch, a.limit, a.workers)


if __name__ == "__main__":
    main()
