#!/usr/bin/env python
"""KRONOS2 run the way CORAL documents it, end to end, for the two representative donors.

    cd core/multiplex5000 && ../../.venv/bin/python scripts/kronos2_doc/run.py --donors 6436 6442

Three steps, each in the environment that owns it:

  1. this driver (``.venv``)          -- cut the frozen frame's sixteen 2048-px windows out of the
                                         qptiff as multi-page TIFFs (DAPI + the 28 declared identity
                                         markers) and write each window's cell-mask labels;
  2. ``scripts/kronos2_doc/worker.py``  (``external/.venv-kronos``) -- CORAL Tutorial 3 verbatim:
                                         ingest, import the cell mask, cell-centred 64-px patches,
                                         ``mean_marker`` baseline and KRONOS2 features;
  3. ``scripts/kronos2_doc/cluster.py`` (``scvi-env``)              -- CORAL Tutorial 4 verbatim:
                                         neighbours(15, use_rep='X') + Leiden(0.5, igraph) + UMAP,
                                         then each cluster's mean-marker profile.

This is NOT a bake-off arm: no ablations, no null gate, no verdict.  It is the documented pipeline
run as documented, so its clusters and their marker profiles can be read directly.  Outputs go to
``data/development/kronos2_documented/``; the arm ``kronos2-coral`` and its records are untouched.
"""
from __future__ import annotations

import argparse
import json
import pathlib
import subprocess
import sys
import time

import numpy as np
import pandas as pd

HERE = pathlib.Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent / "broadtype"))
sys.path.insert(0, str(HERE.parent / "broadtype" / "arms"))
import declarations  # noqa: E402
import pins  # noqa: E402
import kronos2_coral as k2  # noqa: E402  (CANONICAL, NOVEL_ROWS, channel_plan, write_window_tiffs)
import kronos_clusters as v1  # noqa: E402  (frame regions, mask paths, GeoJSON ids, frame cells)

OUT = pins.OUT.parent / "kronos2_documented"
KRONOS_PY = pins.REPO / "external" / ".venv-kronos" / "bin" / "python"
SCVI_PY = pathlib.Path("/home/smith6jt/miniconda3/envs/scvi-env/bin/python")
WORKER = HERE / "worker.py"
CLUSTER = HERE / "cluster.py"
PATCH_PX = 64          # CORAL Tutorial 3: "coral patch --mode cell --patch-size 64"


def prepare(donor: str, jobdir: pathlib.Path) -> dict:
    """Window TIFFs + per-window mask-label lists for the frozen frame."""
    dec = declarations.load()
    markers = sorted(dec["identity_markers"])
    regions = v1.frame_regions(donor)
    oids = v1.geojson_object_ids(donor)
    fc = v1.frame_labels(donor, oids)
    fc = fc.loc[fc["in_expression"]]
    need = sorted(int(r) for r in fc["region_idx"].unique())
    raw, canonical = k2.channel_plan(donor, markers)
    jobdir.mkdir(parents=True, exist_ok=True)
    k2.write_window_tiffs(donor, raw, regions, need, jobdir)
    windows = []
    for r in need:
        cells = fc.loc[fc["region_idx"] == r]
        lp = jobdir / f"w{r:02d}_labels.npy"
        np.save(lp, cells["label"].to_numpy(np.int64))
        windows.append({"idx": r, "tif": str(jobdir / f"w{r:02d}.tif"),
                        "mask": str(v1.mask_path(donor, r)), "labels": str(lp)})
    return {"donor_id": donor, "jobdir": str(jobdir), "windows": windows,
            "channel_names": canonical, "raw_markers": raw, "nuclear_marker": "DAPI",
            "novel_rows": [r for r in k2.NOVEL_ROWS if r["marker_name"] in set(markers)],
            "model_dir": str(k2.KRONOS2_DIR), "mpp": float(pins.PIXEL_SIZE_UM),
            "patch_size": PATCH_PX, "num_workers": 8, "batch_size": 16,
            "out": str(jobdir / "features.npz"), "result": str(jobdir / "features_record.json")}, fc


def run(cmd: list, log: pathlib.Path) -> None:
    print(f"[doc] {' '.join(str(c) for c in cmd[:3])} ... -> {log}", flush=True)
    with open(log, "w") as fh:
        p = subprocess.run([str(c) for c in cmd], stdout=fh, stderr=subprocess.STDOUT)
    if p.returncode != 0:
        raise RuntimeError(f"step failed (rc={p.returncode}); see {log}")


def main(argv=None) -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--donors", nargs="+", default=list(pins.DONORS))
    ap.add_argument("--resolution", type=float, default=0.5)
    ap.add_argument("--skip-features", action="store_true")
    a = ap.parse_args(argv)
    pins.assert_pins()
    OUT.mkdir(parents=True, exist_ok=True)
    logs = OUT / "logs"; logs.mkdir(exist_ok=True)
    summary = {"method": "CORAL Tutorial 3 (per-cell KRONOS2) + Tutorial 4 (Leiden), as documented",
               "written_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
               "scope": "the frozen screen frame: 16 x 2048-px windows per donor, the same cell set every arm was judged on",
               "resolution": a.resolution, "donors": {}}
    for d in a.donors:
        t0 = time.time()
        jobdir = OUT / f"donor_id={d}" / "coral"
        job, fc = prepare(d, jobdir)
        (jobdir / "job.json").write_text(json.dumps(job, indent=1))
        if not a.skip_features or not pathlib.Path(job["out"]).exists():
            run([KRONOS_PY, WORKER, "--job", jobdir / "job.json"], logs / f"features_{d}.log")
        cdir = OUT / f"donor_id={d}" / "clusters"
        run([SCVI_PY, CLUSTER, "--features", job["out"], "--out", cdir, "--resolution", a.resolution],
            logs / f"cluster_{d}.log")
        # join the documented clusters back to the QuPath object_id (the only join key)
        key = fc[["object_id", "label", "region_idx"]].copy()
        key["label"] = key["label"].astype(np.int64); key["region_idx"] = key["region_idx"].astype(np.int32)
        for space in ("kronos2", "mean_marker"):
            cl = pd.read_parquet(cdir / f"clusters__{space}.parquet")
            m = cl.merge(key, on=["label", "region_idx"], how="left", validate="one_to_one")
            if m["object_id"].isna().any():
                raise RuntimeError(f"donor {d} {space}: {int(m['object_id'].isna().sum())} clustered cells lack an object_id")
            m[["object_id", "region_idx", "cluster", "umap_1", "umap_2"]].to_parquet(
                cdir / f"clusters__{space}__by_object_id.parquet", index=False)
        rec = json.loads((jobdir / "features_record.json").read_text())
        crec = json.loads((cdir / "cluster_record.json").read_text())
        summary["donors"][d] = {"features": {k: rec[k] for k in ("n_cells", "embed_dim", "batch_size_used",
                                                                 "patch_size", "patch_um", "seconds_total")},
                                "clusters": crec["spaces"], "kronos2_vs_mean_marker": crec["kronos2_vs_mean_marker"],
                                "seconds": round(time.time() - t0, 1)}
        print(f"[doc] donor {d} done in {(time.time() - t0) / 60:.1f} min", flush=True)
    (OUT / "summary.json").write_text(json.dumps(summary, indent=1, default=str))
    print(json.dumps(summary, indent=1, default=str))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
