#!/usr/bin/env python
"""KRONOS2 per-cell features exactly as CORAL's Tutorial 3 documents them.  Runs in
``external/.venv-kronos`` (CORAL + the kronos2 extra).

The documented cell pipeline, verbatim from ``tutorials/3-Cell-Segmentation-and-Feature-Extraction.ipynb``::

    coral tissue --job-dir <dir>                                   # Tutorial 2, a stated prerequisite
    coral cell   --job-dir <dir> --custom-mask-path <dir>          # Approach B: import an existing mask
    coral patch  --job-dir <dir> --mode cell --patch-size 64
    coral extract --job-dir <dir> --extractor mean_marker --patches cell_<mpp>mpp_64px --subset panel.yaml
    coral extract --job-dir <dir> --extractor KRONOS2     --patches cell_<mpp>mpp_64px --subset panel.yaml \
                  --batch-size 1024 --gpu 0

This module calls the CoralSlide API the tutorial documents as "the *same* steps on a single slide,
in-process", with the tutorial's own arguments.  Two deviations, both recorded in the run record:

* **Patch size.**  ``--patch-size 64``, as documented.  The tutorial notes 64 px is "~24 um here",
  true at the cHL slide's 0.37 um/px; at this cohort's 0.5078 um/px the same 64 px is 32.5 um.  The
  documented number is kept and the micron size recorded.
* **Scope.**  The tutorial ingests one slide.  A PhenoCycler slide here is 23,760 x 24,960 px x 59
  channels and ``convert_to_canonical`` reads it as one array (~70 GB), so the unit is instead the
  frozen screen frame: sixteen 2048-px windows per donor, the same cell set every other arm was
  judged on.  Each window is one CORAL slide store.
* **Batch size.**  The tutorial gives two values: ``--batch-size 1024`` on its CLI line and
  ``batch_size=16`` in the Python section of the same page.  CORAL's own source settles it -- "the
  gold-standard features were produced at batch 16 ... reproduce the gold standard with
  ``batch_size=16``" (``features/kronos2.py``), which also states batch invariance is GPU-dependent
  (drift up to ~1e-4).  16 is therefore what runs here: the reproducible documented value rather than
  the throughput one.  1024 would not fit a 20 GB card at 29 markers x 64 px in any case, because the
  marker-aware ViT tokenises every marker.

Nothing else is changed: patch size 64 px (cell-centred, cell-isolated), the panel subset, the model's
own marker-aware normalisation, and the mean-marker baseline the tutorial pairs with it.
"""
from __future__ import annotations

import argparse
import json
import logging
import pathlib
import sys
import time

import numpy as np
import pandas as pd

logging.getLogger("xformers").setLevel(logging.ERROR)
DOC_BATCH = 16          # CORAL features/kronos2.py gold standard; the CLI line's 1024 is throughput
DOC_BATCH_CLI = 1024


def stored_names(slide, wanted: list) -> list:
    from coral.markers.normalize import normalize_key
    by_key = {}
    for m in slide.markers:
        by_key.setdefault(normalize_key(m), m)
    out, missing = [], []
    for w in wanted:
        m = by_key.get(normalize_key(w))
        (out if m else missing).append(m or w)
    if missing:
        raise RuntimeError(f"{slide.path.name}: markers absent from the store: {missing}")
    return out


def ensure_marker_map(job: dict, jobdir: pathlib.Path):
    """CORAL's cohort marker map: registry hits RESOLVED, everything else NOVEL under its own name."""
    from coral.io.ingest import resolve_cohort_markers
    from coral.markers.marker_map import read_marker_map, write_marker_map
    from coral.markers.normalize import NOVEL_TOKEN, REVIEW_TOKEN

    first = pathlib.Path(job["windows"][0]["tif"])
    res, n_review = resolve_cohort_markers([first], jobdir, channel_names=job["channel_names"],
                                           nuclear_marker=job["nuclear_marker"])
    if n_review:
        df = read_marker_map(jobdir)
        rev = df["status"] == REVIEW_TOKEN
        df.loc[rev, "mapped_canonical_name"] = df.loc[rev, "original_name"]
        df.loc[rev, "status"] = NOVEL_TOKEN
        write_marker_map(df, jobdir)
        res, n_review = resolve_cohort_markers([first], jobdir, channel_names=job["channel_names"],
                                               nuclear_marker=job["nuclear_marker"])
        if n_review:
            raise RuntimeError(f"{n_review} markers still need review in {jobdir/'marker_map.csv'}")
    return res


def build_store(job: dict, w: dict, jobdir: pathlib.Path, res, cfg):
    """coral ingest -> coral tissue -> coral cell --custom-mask-path -> coral patch --mode cell --patch-size 64.

    Tutorial 3 lists Tutorial 2 (tissue) as a prerequisite, and although Approach B's imported mask is
    authoritative and not tissue-restricted, CORAL's novel-marker prepare pass reads the tissue mask
    ("Each needs a tissue mask -- run detect_tissue first"), which this panel needs for the seven
    markers outside KRONOS2's vocabulary.  Otsu is the CLI's default method.
    """
    from coral import CoralSlide
    from coral.io.ingest import convert_to_canonical
    from coral.tissue import OtsuTissueSegmenter

    store = jobdir / f"{pathlib.Path(w['tif']).stem}.zarr"
    if store.exists() and CoralSlide.open(store).status()["ingest"]["status"] == "completed":
        slide = CoralSlide.open(store)
    else:
        slide = convert_to_canonical(w["tif"], jobdir, resolution=res, mpp=float(job["mpp"]),
                                     channel_names=job["channel_names"],
                                     nuclear_marker=job["nuclear_marker"], quiet=True)
    if slide.status().get("tissue[otsu]", {}).get("status") != "completed":
        slide.detect_tissue(OtsuTissueSegmenter(), viz=False)                    # Tutorial 2, default method
    if slide.status()["cells"]["status"] != "completed":
        slide.import_cell_mask(np.load(w["mask"]).astype(np.int32), viz=False)   # ids preserved, authoritative
    slug = cfg.resolved_slug(slide._mpp())
    if slide.status().get(f"patch[{slug}]", {}).get("status") != "completed":
        slide.extract_patches(cfg, viz=False)
    return slide, slug


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description="KRONOS2 the documented way (CORAL Tutorial 3)")
    ap.add_argument("--job", required=True)
    a = ap.parse_args(argv)
    job = json.loads(pathlib.Path(a.job).read_text())

    import torch
    import yaml
    from coral.config.patch import PatchConfig
    from coral.config.subset import Selection
    from coral.features import MeanMarkerExtractor
    from coral.features.kronos2 import Kronos2Extractor
    from coral.features.prepare import prepare_cohort_stats

    t0 = time.time()
    jobdir = pathlib.Path(job["jobdir"]); jobdir.mkdir(parents=True, exist_ok=True)
    res = ensure_marker_map(job, jobdir)
    cfg = PatchConfig(patch_size=int(job["patch_size"]), mode="cell_centered")   # Tutorial 3 step 2
    slides, slugs = [], {}
    for w in job["windows"]:
        s, slug = build_store(job, w, jobdir, res, cfg)
        slides.append(s); slugs[int(w["idx"])] = slug
    t_prep = time.time() - t0

    names = stored_names(slides[0], job["channel_names"])
    (jobdir / "panel.yaml").write_text(yaml.safe_dump({"name": "panel", "channels": {"include": names}}, sort_keys=False))
    sel = Selection(include=names)

    ext = Kronos2Extractor.from_pretrained(job["model_dir"], device="cuda" if torch.cuda.is_available() else "cpu")
    add_csv = jobdir / "additional_markers.csv"
    if job.get("novel_rows"):
        if not add_csv.exists():
            pd.DataFrame(job["novel_rows"]).to_csv(add_csv, index=False)
        prepare_cohort_stats(ext, slides, add_csv, channels=sel)
    novel = ext.novel_markers(names, job["nuclear_marker"])
    if novel:
        raise RuntimeError(f"markers still novel to KRONOS2 after registration: {novel}")

    batch = int(job.get("batch_size", DOC_BATCH))
    print(f"[doc] batch {batch} (CORAL gold standard {DOC_BATCH}; its CLI line's throughput value is {DOC_BATCH_CLI})", flush=True)

    parts = []
    for w, s in zip(job["windows"], slides):
        s.encode_features(MeanMarkerExtractor(), cfg, channels=sel, suffix="markers_panel")     # baseline
        s.encode_features(ext, cfg, channels=sel, suffix="markers_panel", batch_size=batch,
                          num_workers=int(job.get("num_workers", 8)))                            # KRONOS2
        k = s.features(ext, cfg, channels=sel, suffix="markers_panel")
        m = s.features(MeanMarkerExtractor(), cfg, channels=sel, suffix="markers_panel")
        K = np.asarray(k["features"].values, dtype=np.float32)
        M = np.asarray(m["features"].values, dtype=np.float32)
        cid = np.asarray(k["cell_id"].values, dtype=np.int64)
        want = np.load(w["labels"]).astype(np.int64)
        pos = pd.Index(cid).get_indexer(pd.Index(want))
        if (pos < 0).any():
            raise RuntimeError(f"window {w['idx']}: {(pos < 0).sum()} frame cells have no CORAL embedding")
        parts.append({"idx": int(w["idx"]), "K": K[pos], "M": M[pos], "labels": want,
                      "markers": [str(x) for x in np.asarray(m["marker"].values)]})
        print(f"[doc] window {w['idx']:02d}: {len(want):,} cells", flush=True)

    K = np.concatenate([p["K"] for p in parts]); M = np.concatenate([p["M"] for p in parts])
    lab = np.concatenate([p["labels"] for p in parts])
    reg = np.concatenate([np.full(len(p["labels"]), p["idx"], dtype=np.int32) for p in parts])
    np.savez(job["out"], kronos2=K, mean_marker=M, labels=lab, region_idx=reg,
             markers=np.array(parts[0]["markers"], dtype=object), allow_pickle=True)
    rec = {"ok": True, "donor_id": job["donor_id"], "n_cells": int(len(lab)), "embed_dim": int(K.shape[1]),
           "markers": parts[0]["markers"], "batch_size_used": batch,
           "batch_size_gold_standard": DOC_BATCH, "batch_size_cli_line": DOC_BATCH_CLI,
           "patch_size": int(job["patch_size"]), "patch_um": round(int(job["patch_size"]) * float(job["mpp"]), 1),
           "patch_um_in_tutorial": 24, "mode": "cell_centered",
           "mpp": float(job["mpp"]), "slugs": slugs, "panel": names,
           "seconds_prepare": round(t_prep, 1), "seconds_total": round(time.time() - t0, 1),
           "device": str(next(ext._model.parameters()).device),
           "documented_calls": [
               "coral cell --custom-mask-path <dir>  (Approach B: imported mask, ids preserved, authoritative)",
               "coral patch --mode cell --patch-size 64",
               "coral extract --extractor mean_marker --patches cell_<mpp>mpp_64px --subset panel.yaml",
               "coral extract --extractor KRONOS2 --patches cell_<mpp>mpp_64px --subset panel.yaml --batch-size 1024 --gpu 0",
               "(run at batch 16: CORAL's own gold-standard value, features/kronos2.py)"],
           "versions": versions()}
    pathlib.Path(job["result"]).write_text(json.dumps(rec, indent=1, default=str))
    print(json.dumps({k: v for k, v in rec.items() if k not in ("markers", "panel", "batch_probe", "versions")}))
    return 0


def versions() -> dict:
    out = {"python": sys.version.split()[0], "numpy": np.__version__}
    for n in ("torch", "transformers", "timm", "coral", "xformers"):
        try:
            out[n] = str(getattr(__import__(n), "__version__", "installed"))
        except Exception as exc:  # noqa: BLE001
            out[n] = f"missing ({type(exc).__name__})"
    return out


if __name__ == "__main__":
    raise SystemExit(main())
