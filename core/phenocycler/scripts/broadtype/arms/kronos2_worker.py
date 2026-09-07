#!/usr/bin/env python
"""KRONOS2 worker for arm ``kronos2-coral`` -- runs in ``external/.venv-kronos`` (CORAL + KRONOS2 pins).

One job = one (donor, config): the arm hands over per-window TIFFs (written once per donor), label masks,
the marker names in TIFF channel order (registry names AND their CORAL canonical names), the novel-marker
text rows, and the channel selection for this config.  The worker drives CORAL's public API exactly as
its cell tutorial does:

    resolve_cohort_markers -> convert_to_canonical            (one OME-Zarr store per frame window)
    slide.import_cell_mask(mask)                              (labels preserved; = GeoJSON feature index + 1)
    slide.detect_tissue(OtsuTissueSegmenter())                (the prepare pass needs a tissue mask)
    slide.extract_patches(PatchConfig(patch_size, "cell_centered"))
    prepare_cohort_stats(extractor, slides, additional_markers.csv)   (novel markers: (mean, std) from the
                                                               donor's own windows, registered into the model)
    slide.encode_features(Kronos2Extractor, cfg, channels=Selection(include=...), suffix=<config>)
    slide.features(...) -> (cell, 768) + cell_id

and then the arm's own kNN/Leiden (``kronos_worker.cluster_embeddings``) on the CLS embeddings.
Nothing here reads a QuPath class, a legacy label or an expected frequency.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import logging
import pathlib
import sys
import time

import numpy as np
import pandas as pd

HERE = pathlib.Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import kronos_worker  # noqa: E402  (pure clustering helpers; its torch/kronos imports stay lazy)

logging.basicConfig(level=logging.WARNING, format="%(levelname)s %(name)s: %(message)s")
log = logging.getLogger("kronos2_worker")


def _sha256(path: pathlib.Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 22), b""):
            h.update(chunk)
    return h.hexdigest()


def ensure_marker_map(job: dict, slides_dir: pathlib.Path):
    """CORAL's marker map for this donor's windows: registry hits RESOLVED, the rest NOVEL (own name kept)."""
    from coral.io.ingest import resolve_cohort_markers
    from coral.markers.marker_map import read_marker_map, write_marker_map
    from coral.markers.normalize import NOVEL_TOKEN, REVIEW_TOKEN

    first = pathlib.Path(job["windows"][0]["tif"])
    resolution, n_review = resolve_cohort_markers([first], slides_dir, channel_names=job["channel_names"],
                                                  nuclear_marker=job["nuclear_marker"])
    if n_review:
        df = read_marker_map(slides_dir)
        review = df["status"] == REVIEW_TOKEN
        df.loc[review, "mapped_canonical_name"] = df.loc[review, "original_name"]
        df.loc[review, "status"] = NOVEL_TOKEN
        write_marker_map(df, slides_dir)
        resolution, n_review = resolve_cohort_markers([first], slides_dir, channel_names=job["channel_names"],
                                                      nuclear_marker=job["nuclear_marker"])
        if n_review:
            raise RuntimeError(f"{n_review} marker(s) still need review after marking novel: {slides_dir / 'marker_map.csv'}")
    return resolution, {k: [v[0], v[1], bool(v[2])] for k, v in resolution.items()}


def ensure_slide(job: dict, w: dict, slides_dir: pathlib.Path, resolution, cfg):
    """Ingest -> mask -> tissue -> cell patches for one window, each step skipped when already complete."""
    from coral import CoralSlide
    from coral.io.ingest import convert_to_canonical
    from coral.tissue import OtsuTissueSegmenter

    store = slides_dir / f"{pathlib.Path(w['tif']).stem}.zarr"
    if store.exists() and CoralSlide.open(store).status()["ingest"]["status"] == "completed":
        slide = CoralSlide.open(store)
    else:
        slide = convert_to_canonical(w["tif"], slides_dir, resolution=resolution, mpp=float(job["mpp"]),
                                     channel_names=job["channel_names"], nuclear_marker=job["nuclear_marker"], quiet=True)
    st = slide.status()
    if st["cells"]["status"] != "completed":
        mask = np.load(w["mask"]).astype(np.int32)
        slide.import_cell_mask(mask, viz=False)
    if st.get("tissue[otsu]", {}).get("status") != "completed":
        slide.detect_tissue(OtsuTissueSegmenter(), viz=False)
    slug = cfg.resolved_slug(slide._mpp())
    if slide.status().get(f"patch[{slug}]", {}).get("status") != "completed":
        slide.extract_patches(cfg, viz=False)
    return slide, slug


def stored_names(slide, wanted: list) -> list:
    """The slide's own spelling of each wanted marker (CORAL stores cleaned names: 'DC-SIGN' -> 'dc_sign'),
    matched by CORAL's separator- and case-insensitive key; a miss is an error, never a fuzzy guess."""
    from coral.markers.normalize import normalize_key
    by_key = {}
    for m in slide.markers:
        by_key.setdefault(normalize_key(m), m)
    out, missing = [], []
    for w in wanted:
        m = by_key.get(normalize_key(w))
        (out if m else missing).append(m or w)
    if missing:
        raise RuntimeError(f"{slide.path.name}: requested markers absent from the store: {missing}; stored: {slide.markers}")
    return out


def embed_job(job: dict) -> dict:
    import torch
    from coral.config.patch import PatchConfig
    from coral.config.subset import Selection
    from coral.features.kronos2 import Kronos2Extractor
    from coral.features.prepare import prepare_cohort_stats

    t0 = time.time()
    slides_dir = pathlib.Path(job["slides_dir"]); slides_dir.mkdir(parents=True, exist_ok=True)
    resolution, res_rec = ensure_marker_map(job, slides_dir)
    cfg = PatchConfig(patch_size=int(job["patch_size"]), mode="cell_centered")
    slides, slugs = [], {}
    for w in job["windows"]:
        s, slug = ensure_slide(job, w, slides_dir, resolution, cfg)
        slides.append(s); slugs[int(w["idx"])] = slug
    t_ingest = time.time() - t0

    ext = Kronos2Extractor.from_pretrained(job["model_dir"], device="cuda" if torch.cuda.is_available() else "cpu")
    all_names = stored_names(slides[0], list(job["all_channels"]))
    all_sel = Selection(include=all_names)
    add_csv = slides_dir / "additional_markers.csv"
    if job.get("novel_rows"):
        if not add_csv.exists():
            pd.DataFrame(job["novel_rows"]).to_csv(add_csv, index=False)
        prepared = prepare_cohort_stats(ext, slides, add_csv, channels=all_sel)
        novel_after = ext.novel_markers(all_names, job["nuclear_marker"])
        stats = pd.read_csv(add_csv)[["marker_name", "mean", "std"]].to_dict("records") if prepared else []
    else:
        novel_after, stats = ext.novel_markers(all_names, job["nuclear_marker"]), []
    if novel_after:
        raise RuntimeError(f"markers still novel to KRONOS2 after registration: {novel_after}")
    t_prep = time.time() - t0 - t_ingest

    sel = Selection(include=stored_names(slides[0], list(job["select"])))
    embs, labels, regions = [], [], []
    for w, s in zip(job["windows"], slides):
        # num_workers is CORAL's loader parallelism (patch reads overlap the forward pass); it changes wall time only
        s.encode_features(ext, cfg, channels=sel, suffix=job["variant"], batch_size=int(job["batch_size"]),
                          num_workers=int(job.get("num_workers", 0)))
        ds = s.features(ext, cfg, channels=sel, suffix=job["variant"])
        feats = np.asarray(ds["features"].values if "features" in ds else ds.to_array().values[0], dtype=np.float32)
        cid = np.asarray(ds["cell_id"].values, dtype=np.int64)
        want = np.load(w["labels"]).astype(np.int64)
        pos = pd.Index(cid).get_indexer(pd.Index(want))
        if (pos < 0).any():
            raise RuntimeError(f"window {w['idx']}: {(pos < 0).sum()} requested labels have no CORAL embedding")
        embs.append(feats[pos]); labels.append(want); regions.append(np.full(len(want), int(w["idx"]), dtype=np.int32))
    E = np.concatenate(embs); L = np.concatenate(labels); R = np.concatenate(regions)
    np.savez(job["out"], emb=E, labels=L, region_idx=R)
    res = {"ok": True, "n_cells": int(len(L)), "embed_dim": int(E.shape[1]), "seconds_ingest": round(t_ingest, 1),
           "seconds_prepare": round(t_prep, 1), "seconds_total": round(time.time() - t0, 1),
           "model_sha256": _sha256(pathlib.Path(job["model_dir"]) / "kronos2_vitb16_teacher.pth"),
           "marker_metadata_sha256": _sha256(pathlib.Path(job["model_dir"]) / "marker_metadata.csv"),
           "resolution": res_rec, "novel_stats": stats, "slugs": slugs, "device": str(next(ext._model.parameters()).device),
           "versions": versions()}
    if job.get("result"):
        pathlib.Path(job["result"]).write_text(json.dumps(res, indent=1, default=str))
    return res


def versions() -> dict:
    out = {"python": sys.version.split()[0], "numpy": np.__version__}
    for name in ("torch", "xformers", "transformers", "timm", "coral", "huggingface_hub", "pynndescent", "leidenalg", "igraph"):
        try:
            mod = __import__(name)
            out[name] = str(getattr(mod, "__version__", "installed"))
        except Exception as exc:  # noqa: BLE001
            out[name] = f"missing ({type(exc).__name__})"
    try:
        import subprocess
        out["coral_commit"] = subprocess.run(["git", "-C", str(HERE.parents[3] / "external" / "CORAL"), "rev-parse", "HEAD"],
                                            capture_output=True, text=True, check=True).stdout.strip()
    except Exception:  # noqa: BLE001
        pass
    return out


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description="KRONOS2 / CORAL worker (external/.venv-kronos)")
    sub = ap.add_subparsers(dest="cmd", required=True)
    sub.add_parser("versions")
    for name in ("embed", "cluster"):
        p = sub.add_parser(name); p.add_argument("--job", required=True)
    a = ap.parse_args(argv)
    if a.cmd == "versions":
        print(json.dumps(versions())); return 0
    job = json.loads(pathlib.Path(a.job).read_text())
    res = embed_job(job) if a.cmd == "embed" else kronos_worker.cluster_job(job)
    print(json.dumps({k: v for k, v in res.items() if k not in ("resolution", "n_clusters", "novel_stats")}, default=str))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
