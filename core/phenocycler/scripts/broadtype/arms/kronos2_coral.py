"""Arm F2 ``kronos2-coral``: KRONOS2 (marker-aware DINOv2 ViT-B/16) per-cell CLS embeddings obtained
THROUGH CORAL (mahmoodlab/CORAL, the maintainer's fork), -> kNN/Leiden clusters -> null-gated broad-label
annotation.  Development only; frame scope.  Declared 2026-09-06 on the maintainer's instruction that the
Hugging Face approval is for KRONOS2 and that it is to be used through the CORAL GitHub; arm F
``kronos-clusters`` (KRONOS v1, access still pending) stays as recorded.

Method (frozen before the first record; any change is a new arm name)
  1. Assets.  ``external/kronos_assets/KRONOS2/`` = Hugging Face ``MahmoodLab/KRONOS2`` snapshot
     ``e2a54732`` (weights ``kronos2_vitb16_teacher.pth``, ``marker_metadata.csv``, modelling code);
     CORAL checkout ``external/CORAL`` (fork ``smith6jt-cop/CORAL``); both shas recorded per run.
  2. Markers.  The identity input is ``ctx.markers`` (registry type markers only, minus any ablation)
     plus the DAPI plane (CORAL requires the nuclear stain in every selection; DAPI is a structural
     channel, neither identity nor state).  Registry names are mapped to CORAL's canonical registry
     by antigen (``CANONICAL``; e.g. Pan_Cytokeratin -> CYTOKERATIN, CD209 -> DC-SIGN, SMA -> A-SMA).
     Seven markers are NOVEL to the KRONOS2 vocabulary (Ker8_18, Keratin_5, Caveolin, TP63, GCG, INS,
     SST): each gets an ``additional_markers.csv`` row whose compartment/family REUSE existing KRONOS2
     categories (the model's own rule) and whose (mean, std) CORAL computes from this donor's frame
     windows (tissue-masked pixels, its prepare pass) -- recorded per run.  No alias/fuzzy step.
  3. Windows -> CORAL slides.  Per frame region: the 2048-px level-0 qptiff planes of every identity
     marker on the donor + DAPI as one multi-page TIFF (written once per donor), ingested by
     ``convert_to_canonical``; the frame mask (GLOBAL labels = GeoJSON feature index + 1) imported by
     ``import_cell_mask`` (labels preserved); Otsu tissue mask; ``PatchConfig(64 px, cell_centered)``
     -- 64 px at 0.5078 um/px = 32 um cell tiles, CORAL's cell-tutorial size (64 px at 0.37 um/px);
     patches are isolated to the cell (neighbours + background zeroed) by CORAL's CoralDataset.
  4. Embeddings.  ``Kronos2Extractor`` (fp32) on the selected channels, batch 128 (CORAL's cell tutorial
     uses 1024; its batch-invariance note concerns ~1e-4 effects irrelevant to clustering), 768-d CLS
     per cell, cached under ``scratch/kronos2-coral/cache/`` keyed by donor / region / channel set /
     model sha / patch / batch, so seed and subsample replicates reuse them and an ablation never does.
  5. Clusters + annotation: identical to arm F (pynndescent kNN k = 15 + Leiden RBConfiguration at
     resolutions {0.3, 0.5, 1.0} x 5 seeds; ``gate.null_gated_cluster_annotation`` on the MARKER values;
     resolution by seed stability of the annotated labels).

Run from ``core/phenocycler``::

    ../../.venv/bin/python scripts/broadtype/arms/kronos2_coral.py --donor 6436 --scope frame
"""
from __future__ import annotations

import argparse
import hashlib
import json
import os
import pathlib
import shutil
import subprocess
import sys
import time

import numpy as np
import pandas as pd

HERE = pathlib.Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(HERE.parent))
import arm_base  # noqa: E402
import gate  # noqa: E402
import pins  # noqa: E402
import records  # noqa: E402
from criteria import CRITERIA  # noqa: E402
from records import UNAVAILABLE  # noqa: E402
import kronos_clusters as v1  # noqa: E402  (frame/mask/id plumbing, clustering + annotation; nothing KRONOS-v1 is called)
from som_kintsugi import award_rate, shuffle_columns  # noqa: E402  (the false-award helpers, unchanged)

ARM = "kronos2-coral"
STOCHASTIC = True
WORKER = HERE / "kronos2_worker.py"
KRONOS_VENV_PYTHON = v1.KRONOS_VENV_PYTHON
KRONOS2_DIR = pins.REPO / "external" / "kronos_assets" / "KRONOS2"
CORAL_DIR = pins.REPO / "external" / "CORAL"
KRONOS2_HF_REPO = "MahmoodLab/KRONOS2"
KRONOS2_HF_REVISION = "e2a5473206c02d9973885b6445ef783f7e082951"
SCRATCH = pins.OUT / "scratch" / ARM
CACHE = SCRATCH / "cache"
SLIDES = SCRATCH / "coral"          # one CORAL job dir per donor: w<idx>.tif/.zarr, marker_map.csv, additional_markers.csv
NUCLEAR = "DAPI"

REASON_FRAME_ONLY = v1.REASON_FRAME_ONLY
REASON_NO_EMBEDDING = v1.REASON_NO_EMBEDDING
REASON_ASSETS = "kronos2_assets_unavailable"

#: registry marker -> CORAL canonical registry name (same antigen), or the marker's own name when NOVEL
CANONICAL = {
    "CD11c": "CD11C", "CD163": "CD163", "CD20": "CD20", "CD206": "CD206", "CD209": "DC-SIGN", "CD31": "CD31",
    "CD3e": "CD3e", "CD4": "CD4", "CD68": "CD68", "CD79a": "CD79A", "CD8": "CD8", "E_cadherin": "E-CADHERIN",
    "EpCAM": "EPCAM", "FOXP3": "FOXP3", "Iba1": "IBA1", "M2Gal3": "GALECTIN-3", "MPO": "MPO",
    "Pan_Cytokeratin": "CYTOKERATIN", "Podoplanin": "PODOPLANIN", "SMA": "A-SMA", "Vimentin": "VIMENTIN",
    "DAPI": "DAPI",
    # novel to the KRONOS2 vocabulary: own name kept, text row below, (mean, std) from the donor's windows
    "Ker8_18": "Ker8_18", "Keratin_5": "Keratin_5", "Caveolin": "Caveolin", "TP63": "TP63",
    "GCG": "GCG", "INS": "INS", "SST": "SST",
}
#: additional_markers.csv text rows (compartment/family reuse existing KRONOS2 categories; mean/std left
#: blank for CORAL's prepare pass).  Family choices are recorded assumptions, not measurements.
NOVEL_ROWS = [
    {"marker_name": "Ker8_18", "compartment": "cytoplasm", "family": "epithelial_keratin",
     "compartment_desc": "cytoplasmic protein", "family_desc": "Epithelial keratin (cytoskeletal filaments)",
     "marker_full_name": "Cytokeratin 8/18 (KRT8, KRT18)", "mean": "", "std": "", "pretraining": "no"},
    {"marker_name": "Keratin_5", "compartment": "cytoplasm", "family": "epithelial_keratin",
     "compartment_desc": "cytoplasmic protein", "family_desc": "Epithelial keratin (cytoskeletal filaments)",
     "marker_full_name": "Cytokeratin 5 (KRT5)", "mean": "", "std": "", "pretraining": "no"},
    {"marker_name": "Caveolin", "compartment": "membrane", "family": "endothelial_marker",
     "compartment_desc": "membrane protein", "family_desc": "Endothelial cell marker",
     "marker_full_name": "Caveolin-1 (CAV1)", "mean": "", "std": "", "pretraining": "no"},
    {"marker_name": "TP63", "compartment": "nucleus", "family": "transcription_factor",
     "compartment_desc": "nuclear protein", "family_desc": "Nuclear transcription factor",
     "marker_full_name": "Tumor protein p63 (TP63)", "mean": "", "std": "", "pretraining": "no"},
    {"marker_name": "GCG", "compartment": "cytoplasm", "family": "epithelial_marker",
     "compartment_desc": "cytoplasmic protein", "family_desc": "Epithelial marker",
     "marker_full_name": "Glucagon (GCG), pancreatic alpha-cell hormone", "mean": "", "std": "", "pretraining": "no"},
    {"marker_name": "INS", "compartment": "cytoplasm", "family": "epithelial_marker",
     "compartment_desc": "cytoplasmic protein", "family_desc": "Epithelial marker",
     "marker_full_name": "Insulin (INS), pancreatic beta-cell hormone", "mean": "", "std": "", "pretraining": "no"},
    {"marker_name": "SST", "compartment": "cytoplasm", "family": "epithelial_marker",
     "compartment_desc": "cytoplasmic protein", "family_desc": "Epithelial marker",
     "marker_full_name": "Somatostatin (SST), pancreatic delta-cell hormone", "mean": "", "std": "", "pretraining": "no"},
]
NOVEL = {r["marker_name"] for r in NOVEL_ROWS}

DEFAULT_PARAMS: dict = {
    "model": "KRONOS2 marker-aware DINOv2 ViT-B/16 (768-d CLS), through CORAL Kronos2Extractor (fp32)",
    "patch_px": 64,
    "patch_note": "CORAL PatchConfig(patch_size=64, mode='cell_centered'): 32 um at 0.5078 um/px; cell-isolated tiles (CoralDataset)",
    "normalisation": "CORAL dtype scaling (/65535) then KRONOS2 model.preprocess per-marker z-score (marker_metadata.csv; novel markers: (mean, std) from the donor's tissue-masked frame windows via CORAL prepare_cohort_stats)",
    "nuclear_marker": NUCLEAR,
    "mask": "CORAL CoralDataset: raw box * (cell_mask == cell_id), then /65535, then model.preprocess z-score -- background pixels sit at (0 - mean)/std (CORAL's semantics, recorded; arm F masked AFTER normalisation)",
    "false_award_measurement": "annotation-only: gate.null_gated_cluster_annotation on a column-shuffled copy of the standardised marker block with the PRIMARY clusters (awarded / eligible clusters); the arm's clusters are not re-fitted on shuffled images",
    "marker_mapping": "CANONICAL (antigen match to CORAL registry_v1; exact separator-insensitive key inside KRONOS2); novel rows NOVEL_ROWS",
    "embedding": "x_norm_clstoken, 768-d; per-marker/spatial tokens discarded",
    "batch_size": 128,
    "num_workers": 8,          # CORAL loader subprocesses reading patches ahead of the forward pass (wall time only)
    "knn_k": 15, "knn_metric": "euclidean",
    "knn_graph": "pynndescent k+1 neighbours, self removed, symmetrised union, unweighted",
    "leiden": "leidenalg.RBConfigurationVertexPartition, n_iterations=-1 (to convergence)",
    "resolutions": [0.3, 0.5, 1.0],
    "n_seeds": CRITERIA["n_replicates"],
    "seeds": "ctx.seed .. ctx.seed + n_seeds - 1 (Leiden); knn_seed = ctx.seed",
    "resolution_choice": "argmax over resolutions of mean pairwise ARI of the ANNOTATED per-cell labels across seeds; tie -> lowest resolution",
    "annotation": "gate.null_gated_cluster_annotation on gate.standardise_within_donor(ctx.X) (marker values, not the embedding)",
    "null_draws": CRITERIA["null_draws"], "null_quantile": CRITERIA["null_quantile"], "margin_frac": 0.25, "min_cells": CRITERIA["min_calls"],
    "score": "contrast - null p95 of the cluster's awarded label (Ambiguous: ranked-first label); NaN for Unavailable",
    "scope_supported": ["frame"],
    "declared_after": "arm F kronos-clusters (KRONOS v1) never ran: MahmoodLab/KRONOS access pending; the maintainer's approval is for KRONOS2 (2026-09-06), to be used through CORAL",
    "pins": {"kronos2_hf_repo": KRONOS2_HF_REPO, "kronos2_hf_revision": KRONOS2_HF_REVISION,
             "coral": "smith6jt-cop/CORAL", "worker_python": str(KRONOS_VENV_PYTHON)},
}

KEEP_SCRATCH = False
PURGE_CACHE = False


def _sha256(path: pathlib.Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 22), b""):
            h.update(chunk)
    return h.hexdigest()


def ensure_assets() -> dict:
    need = [KRONOS2_DIR / "kronos2_vitb16_teacher.pth", KRONOS2_DIR / "marker_metadata.csv", KRONOS2_DIR / "modeling_kronos2.py",
            CORAL_DIR / "src" / "coral" / "__init__.py"]
    missing = [str(p) for p in need if not p.exists()]
    if missing:
        return {"ok": False, "reason": f"{REASON_ASSETS}:{len(missing)}_files_missing", "missing": missing}
    coral_commit = subprocess.run(["git", "-C", str(CORAL_DIR), "rev-parse", "HEAD"], capture_output=True, text=True).stdout.strip()
    return {"ok": True, "sha256": {"kronos2_vitb16_teacher.pth": _sha256(need[0]), "marker_metadata.csv": _sha256(need[1])},
            "coral_commit": coral_commit, "kronos2_dir": str(KRONOS2_DIR)}


def channel_plan(donor: str, markers: list) -> tuple:
    """(registry names in TIFF order, canonical names in the same order) = every identity marker on the donor + DAPI."""
    full = list(markers)
    unknown = [m for m in full if m not in CANONICAL]
    if unknown:
        raise KeyError(f"{ARM}: registry markers without a CORAL/KRONOS2 mapping: {unknown}")
    raw = [NUCLEAR] + full
    return raw, [CANONICAL[m] for m in raw]


def cache_key(donor: str, region: dict, select: list, assets: dict) -> str:
    w = region["window"]
    payload = json.dumps({"donor": donor, "y0": w["y0"], "x0": w["x0"], "size": w["size"], "select": sorted(select),
                          "model": assets["sha256"]["kronos2_vitb16_teacher.pth"][:16], "coral": assets["coral_commit"][:12],
                          "patch": DEFAULT_PARAMS["patch_px"], "batch": DEFAULT_PARAMS["batch_size"], "nuclear": NUCLEAR}, sort_keys=True)
    return hashlib.sha256(payload.encode()).hexdigest()[:16]


def write_window_tiffs(donor: str, raw_markers: list, regions: list, need: list, sdir: pathlib.Path) -> list:
    """One multi-page uint16 TIFF per needed frame window (channels = raw_markers order), written once."""
    import tifffile
    todo = [r for r in need if not (sdir / f"w{r:02d}.tif").exists() and not (sdir / f"w{r:02d}.zarr").exists()]
    if todo:
        tf, pages, idx, _names, _shape = v1.open_planes(donor, raw_markers)
        try:
            for r in todo:
                planes = v1.read_planes(pages, idx, regions[r]["window"])
                tifffile.imwrite(sdir / f"w{r:02d}.tif", planes, photometric="minisblack", metadata={"axes": "CYX"})
                del planes
        finally:
            tf.close()
    return todo


def embed_frame(ctx: arm_base.ArmContext, assets: dict, exchange: pathlib.Path) -> tuple:
    donor = ctx.donor_id
    regions = v1.frame_regions(donor)
    oids = v1.geojson_object_ids(donor)
    fc = v1.frame_labels(donor, oids)
    fc = fc.loc[fc["in_expression"]]
    need = sorted(set(int(r) for r in np.unique(ctx.region_idx)))
    full_markers = list(ctx.markers) + ([ctx.ablated] if ctx.ablated else [])
    raw, canonical = channel_plan(donor, sorted(full_markers))
    select = [CANONICAL[m] for m in [NUCLEAR] + list(ctx.markers)]
    sdir = SLIDES / f"donor_id={donor}"
    sdir.mkdir(parents=True, exist_ok=True)
    keys = {r: cache_key(donor, regions[r], select, assets) for r in need}
    cdir = CACHE / f"donor_id={donor}"
    cdir.mkdir(parents=True, exist_ok=True)
    todo = [r for r in need if not (cdir / f"w{r:02d}__{keys[r]}.npz").exists()]
    info = {"regions": need, "cache_hits": [r for r in need if r not in todo], "channels_raw": raw, "channels_canonical": canonical,
            "select": select, "novel": sorted(NOVEL & set(full_markers)), "slides_dir": str(sdir)}
    if todo:
        written = write_window_tiffs(donor, raw, regions, todo, sdir)
        info["tiffs_written"] = written
        exchange.mkdir(parents=True, exist_ok=True)
        windows = []
        for r in todo:
            cells = fc.loc[fc["region_idx"] == r]
            np.save(exchange / f"w{r:02d}_labels.npy", cells["label"].to_numpy(np.int64))
            windows.append({"idx": r, "tif": str(sdir / f"w{r:02d}.tif"), "mask": str(v1.mask_path(donor, r)),
                            "labels": str(exchange / f"w{r:02d}_labels.npy")})
        job = {"donor_id": donor, "config_id": ctx.config_id, "slides_dir": str(sdir), "windows": windows,
               "channel_names": canonical, "raw_markers": raw, "nuclear_marker": NUCLEAR, "all_channels": canonical,
               "novel_rows": [r_ for r_ in NOVEL_ROWS if r_["marker_name"] in set(full_markers)],
               "select": select, "variant": ctx.config_id.replace(":", "_"), "model_dir": str(KRONOS2_DIR),
               "mpp": float(pins.PIXEL_SIZE_UM), "patch_size": DEFAULT_PARAMS["patch_px"], "batch_size": DEFAULT_PARAMS["batch_size"],
               "num_workers": DEFAULT_PARAMS["num_workers"],
               "out": str(exchange / "emb.npz"), "result": str(exchange / "embed_result.json")}
        res = run_worker("embed", job, exchange / "embed_job.json")
        info["embedded"] = {k: v for k, v in res.items() if k != "resolution"}
        z = np.load(exchange / "emb.npz")
        for r in todo:
            sel = z["region_idx"] == r
            cells = fc.loc[fc["region_idx"] == r]
            lab = z["labels"][sel]
            if len(lab) != len(cells) or not np.array_equal(np.sort(lab), np.sort(cells["label"].to_numpy(np.int64))):
                raise RuntimeError(f"region {r}: worker returned {len(lab)} embeddings for {len(cells)} cells")
            order = pd.Index(lab).get_indexer(pd.Index(cells["label"].to_numpy(np.int64)))
            np.savez(cdir / f"w{r:02d}__{keys[r]}.npz", emb=z["emb"][sel][order], labels=cells["label"].to_numpy(np.int64),
                     object_ids=cells["object_id"].to_numpy().astype("S36"))
    parts_e, parts_id = [], []
    for r in need:
        zz = np.load(cdir / f"w{r:02d}__{keys[r]}.npz")
        parts_e.append(zz["emb"]); parts_id.append(zz["object_ids"].astype(str))
    E_all = np.concatenate(parts_e); ids_all = np.concatenate(parts_id)
    pos = pd.Index(ids_all).get_indexer(pd.Index(ctx.ids.astype(str)))
    have = pos >= 0
    E = np.full((len(ctx.ids), E_all.shape[1]), np.nan, dtype=np.float32)
    E[have] = E_all[pos[have]]
    info["cache_keys"] = {int(k): v for k, v in keys.items()}
    info["n_with_embedding"] = int(have.sum())
    return E, have, info


def run_worker(cmd: str, job: dict, job_path: pathlib.Path) -> dict:
    job_path.write_text(json.dumps(job, indent=1))
    proc = subprocess.run([str(KRONOS_VENV_PYTHON), str(WORKER), cmd, "--job", str(job_path)], capture_output=True, text=True)
    if proc.returncode != 0:
        raise RuntimeError(f"{ARM} worker {cmd} failed (rc={proc.returncode}):\n{proc.stderr[-4000:]}")
    if job.get("result") and pathlib.Path(job["result"]).exists():
        return json.loads(pathlib.Path(job["result"]).read_text())
    return json.loads(proc.stdout.strip().splitlines()[-1])


def _sidecar(ctx: arm_base.ArmContext, payload: dict) -> pathlib.Path:
    d = records.out_dir(ARM, ctx.donor_id, ctx.scope)
    d.mkdir(parents=True, exist_ok=True)
    p = d / f"kronos2_details__{ctx.config_id.replace(':', '_')}.json"
    p.write_text(json.dumps(payload, indent=1, default=str))
    return p


def run(ctx: arm_base.ArmContext) -> pd.DataFrame:
    t0 = time.time()
    DEFAULT_PARAMS["run"] = {"donor_id": ctx.donor_id, "config_id": ctx.config_id, "seed": int(ctx.seed),
                             "identity_markers": list(ctx.markers), "ablated": ctx.ablated}
    if ctx.scope != "frame":
        DEFAULT_PARAMS["run"]["deviation"] = REASON_FRAME_ONLY
        return v1._all(ctx, UNAVAILABLE, REASON_FRAME_ONLY)
    if PURGE_CACHE and CACHE.exists():
        shutil.rmtree(CACHE)
    assets = ensure_assets()
    DEFAULT_PARAMS["run"]["assets"] = assets
    if not assets["ok"]:
        return v1._all(ctx, UNAVAILABLE, assets["reason"])
    DEFAULT_PARAMS["external_pins"] = {"kronos2_weights_sha256": assets["sha256"]["kronos2_vitb16_teacher.pth"],
                                       "kronos2_marker_metadata_sha256": assets["sha256"]["marker_metadata.csv"],
                                       "coral_commit": assets["coral_commit"]}
    if not KRONOS_VENV_PYTHON.exists():
        raise FileNotFoundError(f"worker interpreter missing: {KRONOS_VENV_PYTHON}")
    exchange = SCRATCH / f"run__{ctx.donor_id}__{ctx.config_id.replace(':', '_')}__{os.getpid()}"
    try:
        E, have, einfo = embed_frame(ctx, assets, exchange)
        DEFAULT_PARAMS["run"].update({k: einfo[k] for k in ("channels_raw", "channels_canonical", "select", "novel", "cache_keys", "n_with_embedding")})
        if "embedded" in einfo:
            DEFAULT_PARAMS["run"]["embedded"] = einfo["embedded"]
        if have.sum() <= DEFAULT_PARAMS["knn_k"] + 1:
            raise RuntimeError(f"only {have.sum()} cells carry an embedding; cannot build the graph")
        seeds = [int(ctx.seed) + i for i in range(int(DEFAULT_PARAMS["n_seeds"]))]
        exchange.mkdir(parents=True, exist_ok=True)
        np.save(exchange / "embeddings.npy", E[have])
        job = {"embeddings": str(exchange / "embeddings.npy"), "out": str(exchange / "fits.npz"),
               "resolutions": DEFAULT_PARAMS["resolutions"], "seeds": seeds, "k": DEFAULT_PARAMS["knn_k"],
               "knn_seed": int(ctx.seed), "metric": DEFAULT_PARAMS["knn_metric"], "n_iterations": -1,
               "result": str(exchange / "cluster_result.json")}
        cinfo = run_worker("cluster", job, exchange / "cluster_job.json")
        fits = {k: v for k, v in np.load(exchange / "fits.npz").items()}
    finally:
        if not KEEP_SCRATCH and exchange.exists():
            shutil.rmtree(exchange, ignore_errors=True)
    Z_all = gate.standardise_within_donor(ctx.X.to_numpy(dtype=float))
    declared = v1.declared_columns(list(ctx.markers))
    lab, reason, score, details = v1.annotate_and_choose(Z_all[have], fits, DEFAULT_PARAMS["resolutions"], seeds, declared,
                                                         primary_seed=int(ctx.seed))
    # annotation-only false-award rate: the gate on a structureless twin of the marker block, same clusters
    prim = fits[details["primary_fit"]]
    ann_only = gate.null_gated_cluster_annotation(shuffle_columns(Z_all[have], int(ctx.seed)), prim, declared, seed=int(ctx.seed))
    DEFAULT_PARAMS["false_award"] = {"annotation_only": award_rate(ann_only), "false_award_rate": award_rate(ann_only)["rate"],
                                     "definition": DEFAULT_PARAMS["false_award_measurement"]}
    calls = v1._all(ctx, UNAVAILABLE, REASON_NO_EMBEDDING)
    idx = np.flatnonzero(have)
    calls.loc[idx, "label"] = lab
    calls.loc[idx, "reason"] = reason
    calls.loc[idx, "score"] = score
    DEFAULT_PARAMS["run"].update({"chosen_resolution": details["chosen_resolution"], "primary_fit": details["primary_fit"],
                                  "stability_mean_pairwise_ari": details["stability_mean_pairwise_ari"],
                                  "n_clusters": details["n_clusters"], "graph": {k: v for k, v in cinfo.items() if k != "n_clusters"},
                                  "declared_columns": declared, "seconds": round(time.time() - t0, 1)})
    _sidecar(ctx, {"arm": ARM, "params": DEFAULT_PARAMS, "details": details, "embed": einfo})
    return calls


def _preparse() -> None:
    global KEEP_SCRATCH, PURGE_CACHE
    ap = argparse.ArgumentParser(add_help=False)
    ap.add_argument("--keep-scratch", action="store_true")
    ap.add_argument("--purge-cache", action="store_true")
    a, rest = ap.parse_known_args()
    KEEP_SCRATCH, PURGE_CACHE = a.keep_scratch, a.purge_cache
    sys.argv = [sys.argv[0]] + rest


if __name__ == "__main__":
    _preparse()
    raise SystemExit(arm_base.cli(ARM, run, env=str(KRONOS_VENV_PYTHON), default_params=DEFAULT_PARAMS))
