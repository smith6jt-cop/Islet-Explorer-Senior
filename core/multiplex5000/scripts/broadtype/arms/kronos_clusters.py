"""Arm F ``kronos-clusters``: KRONOS ViT-S/16 cell embeddings -> kNN/Leiden clusters -> null-gated
broad-label annotation.  Development only; frame scope.

Method (frozen 2026-09-05; any change is a new arm name)
  1. Assets.  ``external/kronos_assets/{kronos_vits16_model.pt, marker_metadata.csv}`` (and
     ``config.json`` when obtainable) from Hugging Face ``MahmoodLab/KRONOS`` at the pinned revision
     ``KRONOS_HF_REVISION``.  Missing files are fetched ONCE per run with the cached token.  The repo
     is gated and the maintainer's request was "awaiting a review from the repo authors" on
     2026-09-05, so a 403 writes ``Unavailable(huggingface_access_pending)`` for every cell and the
     record is still written -- the screen sees a row per cell, never a silently missing arm.
  2. Marker mapping.  ``ctx.markers`` (registry type markers only, already minus any ablation) ->
     ``marker_metadata.csv`` rows: uppercased name, the alias table ``kronos_marker_aliases.json``, a
     punctuation-free comparison, then (second pass) a declared family fallback.  Unmatched markers
     get an unused ODD ``marker_id`` adjacent to the metadata marker with the highest name similarity
     -- the tutorial's own suggestion rule; odd ids are the ones KRONOS reserves for markers it was
     not trained on -- and mean/std computed from this donor's frame window planes on the 0-1 scale.
     Every assignment and the step that produced it is recorded.
  3. Patches.  Per frame region: the 2048-px level-0 qptiff planes of exactly ``ctx.markers`` via
     ``common/tiff_window.read_window``; the frame mask ``frame/masks/donor_id=<d>/w<idx>.npy`` with
     GLOBAL labels (``object_ids[label - 1]`` is the UUID; the id list is streamed from the donor
     GeoJSON and checked against the manifest's ``object_id_sha256`` and feature count, and every
     frame cell's label must be present in its window).  The worker (``kronos_worker.py`` in
     ``external/.venv-kronos``) cuts a 64-px cell-centred patch per cell, normalises per marker,
     multiplies by the cell's own mask (the tutorial idiom) and keeps the 384-d CLS embedding.
     Embeddings are cached under ``scratch/kronos-clusters/cache/`` keyed by donor / region /
     marker list + ids + normalisation + model sha, so seed and subsample replicates re-use them and
     an ablation (a different marker list) never does.
  4. Clusters.  pynndescent kNN (k = 15, euclidean, symmetrised, unweighted) + Leiden
     (RBConfiguration, run to convergence) at resolutions {0.3, 0.5, 1.0} x seeds
     {ctx.seed .. ctx.seed + 4}, in the worker.
  5. Annotation.  Every fit is annotated with ``gate.null_gated_cluster_annotation`` on
     ``Z = gate.standardise_within_donor(ctx.X)`` -- the MARKER values of the same cells, never the
     embedding -- so a cluster is awarded a label only where its declared-marker contrast beats a
     size-matched null and leads by the 25 % margin.  The resolution is chosen by SEED STABILITY of
     the annotated labels (mean pairwise ARI across the five seeds; ties -> the lowest resolution).
     Primary = the chosen resolution at ``ctx.seed``; ``score`` = contrast - null p95 of the cluster's
     awarded (or, for Ambiguous, ranked-first) label; reasons come straight from the annotation
     (``multi_lineage``, ``no_lineage_beats_null``, ``cluster_below_min_cells``).

Nothing here consults a QuPath class, a legacy label, another arm, a state marker or an expected
frequency; the number of clusters is whatever Leiden returns and the resolution is chosen by
stability, never by agreement with anything.

Deviations from the contract, recorded here and in the report:
  * ``--scope donor`` is not run.  The protocol runs image arms on the frame and whole-slide only
    as survivors, and whole-slide KRONOS embedding of ~1.9M cells is hours of GPU; every cell is
    written ``Unavailable(image_arm_frame_only)``.
  * The per-run scratch exchange directory is deleted afterwards (``--keep-scratch`` keeps it);
    the embedding CACHE under ``scratch/kronos-clusters/cache/`` persists by design because it is
    content-keyed (``--purge-cache`` removes it first).
  * A download failure other than the gated 403 is recorded under its own reason
    (``kronos_assets_unavailable:<ExceptionName>``), not as access pending.

Run from ``core/multiplex5000``::

    ../../.venv/bin/python scripts/broadtype/arms/kronos_clusters.py --donor 6436 --scope frame
"""
from __future__ import annotations

import argparse
import hashlib
import json
import os
import pathlib
import re
import shutil
import subprocess
import sys
import time
from difflib import SequenceMatcher

import numpy as np
import pandas as pd

HERE = pathlib.Path(__file__).resolve().parent            # scripts/broadtype/arms
sys.path.insert(0, str(HERE.parent))                        # scripts/broadtype
sys.path.insert(0, str(HERE))
import arm_base  # noqa: E402
import declarations  # noqa: E402
import gate  # noqa: E402
import pins  # noqa: E402
import records  # noqa: E402
from criteria import CRITERIA  # noqa: E402
from records import AMBIGUOUS, UNAVAILABLE  # noqa: E402
import kronos_worker  # noqa: E402  (imported for its pure helpers; torch stays lazy)

ARM = "kronos-clusters"
STOCHASTIC = True   # honours ctx.seed; the screen runs seeds 1..5 as replicates

# --------------------------------------------------------------------------- pins and paths
KRONOS_VENV_PYTHON = pins.REPO / "external" / ".venv-kronos" / "bin" / "python"
KRONOS_ASSETS = pins.REPO / "external" / "kronos_assets"
KRONOS_HF_REPO = "MahmoodLab/KRONOS"
#: The model repo commit listed by ``HfApi().model_info`` on 2026-09-05 (files: config.json,
#: kronos_vits16_model.pt, marker_metadata.csv, pytorch_model.bin, meta.yaml).  Downloads pin it.
KRONOS_HF_REVISION = "e1c7d9aab62c98a8341c3f6777cdf8f1339b07d0"
#: The KRONOS code installed in ``external/.venv-kronos`` (fork smith6jt-cop/KRONOS).
KRONOS_CODE_COMMIT = "26f1bc02"
KRONOS_DIST_VERSION = "0.1.0"
ASSET_FILES = ("kronos_vits16_model.pt", "marker_metadata.csv")
OPTIONAL_ASSET_FILES = ("config.json",)
WORKER = HERE / "kronos_worker.py"
ALIASES = HERE / "kronos_marker_aliases.json"
SCRATCH = pins.OUT / "scratch" / ARM
CACHE = SCRATCH / "cache"
HF_TOKEN_FILE = pathlib.Path.home() / ".cache" / "huggingface" / "token"

REASON_ACCESS_PENDING = "huggingface_access_pending"
REASON_FRAME_ONLY = "image_arm_frame_only"
REASON_NO_EMBEDDING = "no_embedding_for_cell"

#: Recorded with every record.  ``run`` is filled in by :func:`run` (arm_base.cli hands this very
#: dict to ``write`` AFTER ``run`` returns, so the run-specific assignments travel with it).
DEFAULT_PARAMS: dict = {
    "model": "kronos vits16 (ViT-S/16, DINOv2-style, 16 register tokens)",
    "model_type": "vits16",
    "token_overlap": True,
    "patch_px": 64,
    "marker_max_value": 65535.0,
    "normalisation": "x / marker_max_value -> (x - marker_mean) / marker_std per marker (KRONOS marker_metadata; dataset frame-window stats for unmatched)",
    "mask": "patch * (mask == own label), applied after normalisation (tutorial feature_extraction idiom)",
    "cell_centre": "integer centroid of the cell's own mask pixels in the frame window",
    "embedding": "patch_features = CLS token (x_norm_clstoken), 384-d; per-marker features discarded",
    "batch_size": 128,
    "knn_k": 15,
    "knn_metric": "euclidean",
    "knn_graph": "pynndescent k+1 neighbours, self removed, symmetrised union, unweighted",
    "leiden": "leidenalg.RBConfigurationVertexPartition, n_iterations=-1 (to convergence)",
    "resolutions": [0.3, 0.5, 1.0],
    "n_seeds": CRITERIA["n_replicates"],
    "seeds": "ctx.seed .. ctx.seed + n_seeds - 1 (Leiden); knn_seed = ctx.seed",
    "resolution_choice": "argmax over resolutions of mean pairwise ARI of the ANNOTATED per-cell labels across seeds; tie -> lowest resolution",
    "annotation": "gate.null_gated_cluster_annotation on gate.standardise_within_donor(ctx.X) (marker values, not the embedding)",
    "null_draws": CRITERIA["null_draws"],
    "null_quantile": CRITERIA["null_quantile"],
    "margin_frac": 0.25,
    "min_cells": CRITERIA["min_calls"],
    "score": "contrast - null p95 of the cluster's awarded label (Ambiguous: ranked-first label); NaN for Unavailable",
    "scope_supported": ["frame"],
    "pins": {
        "kronos_code": f"smith6jt-cop/KRONOS@{KRONOS_CODE_COMMIT} (dist KRONOS {KRONOS_DIST_VERSION})",
        "kronos_hf_repo": KRONOS_HF_REPO,
        "kronos_hf_revision": KRONOS_HF_REVISION,
        "worker_python": str(KRONOS_VENV_PYTHON),
    },
}

KEEP_SCRATCH = False
PURGE_CACHE = False


# --------------------------------------------------------------------------- assets
def _hf_download(filename: str, dest: pathlib.Path) -> str:
    """One ``hf_hub_download`` at the pinned revision into ``dest``.  Separate so tests can
    monkeypatch it; any exception propagates to :func:`ensure_assets`."""
    from huggingface_hub import hf_hub_download
    token = HF_TOKEN_FILE.read_text().strip() if HF_TOKEN_FILE.exists() else None
    return hf_hub_download(KRONOS_HF_REPO, filename, revision=KRONOS_HF_REVISION, local_dir=str(dest), token=token)


def _sha256(path: pathlib.Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 22), b""):
            h.update(chunk)
    return h.hexdigest()


def _is_gated(exc: BaseException) -> bool:
    text = f"{type(exc).__name__} {exc}"
    return "GatedRepoError" in text or "awaiting a review" in text or "403" in text


def ensure_assets(assets_dir: pathlib.Path = KRONOS_ASSETS, *, fetch: bool = True) -> dict:
    """Make the two required assets present, fetching once if allowed.  Never raises on a
    download failure: it reports ``ok=False`` with the reason the record will carry."""
    assets_dir.mkdir(parents=True, exist_ok=True)
    out = {"dir": str(assets_dir), "ok": False, "reason": "", "fetched": [], "errors": {}, "sha256": {}}
    missing = [f for f in ASSET_FILES if not (assets_dir / f).exists()]
    if fetch and missing:
        for f in list(missing) + [f for f in OPTIONAL_ASSET_FILES if not (assets_dir / f).exists()]:
            try:
                _hf_download(f, assets_dir)
                out["fetched"].append(f)
            except Exception as exc:  # noqa: BLE001  -- recorded, then decides the reason
                out["errors"][f] = f"{type(exc).__name__}: {str(exc).splitlines()[0][:200]}"
                if _is_gated(exc):
                    out["gated"] = True
    missing = [f for f in ASSET_FILES if not (assets_dir / f).exists()]
    if missing:
        errs = out["errors"]
        if out.get("gated") or not errs and not fetch:
            out["reason"] = REASON_ACCESS_PENDING
        elif errs:
            first = next(iter(errs.values())).split(":")[0]
            out["reason"] = f"kronos_assets_unavailable:{first}"
        else:
            out["reason"] = REASON_ACCESS_PENDING
        out["missing"] = missing
        return out
    for f in ASSET_FILES + OPTIONAL_ASSET_FILES:
        p = assets_dir / f
        if p.exists():
            out["sha256"][f] = _sha256(p)
    out["ok"] = True
    return out


# --------------------------------------------------------------------------- marker mapping
def load_marker_metadata(path: pathlib.Path) -> pd.DataFrame:
    """``marker_metadata.csv`` indexed by UPPERCASE ``marker_name``; the four tutorial columns."""
    df = pd.read_csv(path)
    need = ("marker_name", "marker_id", "marker_mean", "marker_std")
    missing = [c for c in need if c not in df.columns]
    if missing:
        raise ValueError(f"{path}: marker_metadata is missing columns {missing}")
    df["marker_name"] = df["marker_name"].astype(str).str.strip().str.upper()
    df["marker_id"] = df["marker_id"].astype(int)
    if df["marker_name"].duplicated().any():
        raise ValueError(f"{path}: duplicate marker_name rows")
    return df.set_index("marker_name")


def load_aliases(path: pathlib.Path = ALIASES) -> dict:
    d = json.loads(path.read_text())
    return {"aliases": {k: [str(v).upper() for v in vs] for k, vs in d.get("aliases", {}).items()},
            "fallbacks": {k: [str(v).upper() for v in vs] for k, vs in d.get("fallbacks", {}).items()}}


def _squash(s: str) -> str:
    return re.sub(r"[^A-Z0-9]", "", str(s).upper())


def _odd_id_near(base: int, taken: set) -> int:
    """The first unused ODD id in 1..512 closest to ``base`` (base+1, base-1, base+3, ...)."""
    for d in range(1, 513):
        for c in (base + d, base - d):
            if 1 <= c <= 512 and c % 2 == 1 and c not in taken:
                return c
    raise RuntimeError("no odd marker_id left in 1..512")


def map_markers(markers: list, meta: pd.DataFrame, aliases: dict) -> pd.DataFrame:
    """One row per marker in ``markers`` order: ``marker, metadata_name, marker_id, mean, std, mode``.

    ``mode`` is one of ``exact`` / ``alias`` / ``squashed`` / ``fallback`` / ``odd_id``.  For
    ``odd_id`` rows ``mean``/``std`` are NaN here; :func:`fill_dataset_stats` supplies them from
    the frame planes, and ``suggestion`` / ``similarity`` record the tutorial-style name match
    that placed the id.  An id is never handed to two channels.
    """
    names = list(meta.index)
    squashed = {}
    for n in names:
        squashed.setdefault(_squash(n), n)
    taken = set()
    rows = {m: None for m in markers}

    def claim(m, name, mode):
        mid = int(meta.loc[name, "marker_id"])
        if mid in taken:
            return False
        taken.add(mid)
        rows[m] = {"marker": m, "metadata_name": name, "marker_id": mid,
                   "mean": float(meta.loc[name, "marker_mean"]), "std": float(meta.loc[name, "marker_std"]),
                   "mode": mode, "suggestion": "", "similarity": np.nan}
        return True

    # pass 1: the marker's own identity (exact, alias, squashed)
    for m in markers:
        cands = [(m.upper(), "exact")] + [(a, "alias") for a in aliases["aliases"].get(m, [])]
        for name, mode in cands:
            if name in meta.index and claim(m, name, mode):
                break
        if rows[m] is None:
            for name, _ in cands:
                hit = squashed.get(_squash(name))
                if hit is not None and claim(m, hit, "squashed"):
                    break
    # pass 2: declared family fallbacks, only for markers still unmatched
    for m in markers:
        if rows[m] is None:
            for name in aliases["fallbacks"].get(m, []):
                if name in meta.index and claim(m, name, "fallback"):
                    break
    # pass 3: unused odd id next to the most name-similar metadata marker (tutorial suggestion rule)
    meta_ids = set(int(v) for v in meta["marker_id"])
    for m in markers:
        if rows[m] is None:
            sims = np.array([SequenceMatcher(None, m.upper(), n).ratio() for n in names])
            j = int(np.argmax(sims))
            mid = _odd_id_near(int(meta.loc[names[j], "marker_id"]), taken | meta_ids)
            taken.add(mid)
            rows[m] = {"marker": m, "metadata_name": None, "marker_id": mid, "mean": np.nan, "std": np.nan,
                       "mode": "odd_id", "suggestion": names[j], "similarity": float(sims[j])}
    df = pd.DataFrame([rows[m] for m in markers])
    if df["marker_id"].duplicated().any():
        raise RuntimeError("marker_id assigned twice")
    if not df["marker_id"].between(1, 512).all():
        raise RuntimeError("marker_id outside 1..512")
    return df


def fill_dataset_stats(mapping: pd.DataFrame, stats: dict) -> pd.DataFrame:
    """Give every ``odd_id`` row the mean/std computed from this donor's frame planes (0-1 scale)."""
    out = mapping.copy()
    for i, r in out.iterrows():
        if r["mode"] == "odd_id":
            s = stats[r["marker"]]
            if not (np.isfinite(s["std"]) and s["std"] > 0):
                raise RuntimeError(f"{r['marker']}: dataset std is {s['std']} -- cannot normalise")
            out.loc[i, "mean"] = s["mean"]
            out.loc[i, "std"] = s["std"]
    return out


# --------------------------------------------------------------------------- frame geometry
def frame_regions(donor: str) -> list:
    reg = json.loads((pins.OUT / "frame" / "regions.json").read_text())
    return reg["donors"][str(donor)]["regions"]


def mask_path(donor: str, idx: int) -> pathlib.Path:
    return pins.OUT / "frame" / "masks" / f"donor_id={donor}" / f"w{idx:02d}.npy"


def geojson_object_ids(donor: str, cache_dir: pathlib.Path = CACHE) -> np.ndarray:
    """The donor GeoJSON's feature ids in stream order (``object_ids[label - 1]`` = UUID), cached
    by the manifest's GeoJSON sha and verified against its feature count and ``object_id_sha256``."""
    from phenocycler.artifacts import hash_identifiers, iter_qupath_geojson_features
    im = pins.manifest_image(donor)
    gj_sha = str(im["cell_geojson"]["sha256"])
    want_n = int(im["cell_geometry"]["feature_count"])
    want_hash = str(im["cell_geometry"]["object_id_sha256"])
    p = cache_dir / f"donor_id={donor}" / f"geojson_object_ids__{gj_sha[:12]}.npy"
    if p.exists():
        ids = np.load(p)
    else:
        ids = np.array([str(ft.get("id")) for ft in iter_qupath_geojson_features(pins.geojson_path(donor))], dtype="S36")
        p.parent.mkdir(parents=True, exist_ok=True)
        np.save(p, ids)
    if len(ids) != want_n:
        raise RuntimeError(f"donor {donor}: {len(ids)} GeoJSON features, manifest says {want_n}")
    got_hash = hash_identifiers(x.decode() for x in ids)
    if got_hash != want_hash:
        raise RuntimeError(f"donor {donor}: object_id hash {got_hash} != manifest {want_hash}")
    return ids


def frame_labels(donor: str, oids: np.ndarray) -> pd.DataFrame:
    """Frame cells with their GLOBAL mask label, and the check that every label is in its window."""
    fc = arm_base.frame_cells(donor).copy()
    index = pd.Index(oids.astype(str))
    pos = index.get_indexer(fc["object_id"].astype(str))
    if (pos < 0).any():
        raise RuntimeError(f"donor {donor}: {(pos < 0).sum()} frame cells are not in the GeoJSON id list")
    fc["label"] = (pos + 1).astype(np.int64)
    for r, grp in fc.groupby("region_idx"):
        m = np.load(mask_path(donor, int(r)))
        present = np.isin(grp["label"].to_numpy(), np.unique(m))
        if not present.all():
            raise RuntimeError(f"donor {donor} region {r}: {(~present).sum()} frame cells have no pixels in the window mask "
                               "-- the label invariant is broken; the lane stops here")
    return fc


def open_planes(donor: str, markers: list):
    """(tiff, pages, page index per marker, channel names) for exactly ``markers``, level 0."""
    from phenocycler.redsea import qptiff_channels
    im = pins.manifest_image(donor)
    cm = {c["marker"]: c for c in im["channel_map"]}
    missing = [m for m in markers if m not in cm]
    if missing:
        raise KeyError(f"donor {donor}: markers without a qptiff channel: {missing}")
    tf, pages, names, shape, _px = qptiff_channels(pins.qptiff_path(donor), 0)
    idx = []
    for m in markers:
        i = int(cm[m]["channel_index"])
        if names[i] != cm[m]["channel_name"]:
            tf.close()
            raise RuntimeError(f"donor {donor}: page {i} is {names[i]!r}, manifest says {cm[m]['channel_name']!r}")
        idx.append(i)
    return tf, pages, idx, [names[i] for i in idx], shape


def read_planes(pages, idx: list, window: dict) -> np.ndarray:
    """``(M, size, size)`` uint16 level-0 planes of one frame window (``tiff_window.read_window``)."""
    from tiff_window import read_window
    y0, x0, size = int(window["y0"]), int(window["x0"]), int(window["size"])
    out = np.zeros((len(idx), size, size), dtype=np.uint16)
    for k, i in enumerate(idx):
        a = read_window(pages[i], y0, y0 + size, x0, x0 + size)
        if a.dtype != np.uint16:
            raise TypeError(f"plane dtype {a.dtype}; marker_max_value {DEFAULT_PARAMS['marker_max_value']} assumes uint16")
        out[k, :a.shape[0], :a.shape[1]] = a
    return out


def dataset_stats(pages, idx: list, windows: list, markers: list, max_value: float) -> dict:
    """Mean/std of every listed marker over all pixels of the frame windows, on the 0-1 scale --
    the tutorial's instruction for markers KRONOS never saw ('intensities in range 0-1')."""
    acc = {m: [0.0, 0.0, 0] for m in markers}
    for w in windows:
        planes = read_planes(pages, idx, w).astype(np.float64) / max_value
        for k, m in enumerate(markers):
            x = planes[k]
            acc[m][0] += float(x.sum()); acc[m][1] += float((x * x).sum()); acc[m][2] += int(x.size)
    out = {}
    for m, (s, ss, n) in acc.items():
        mean = s / n
        var = max(ss / n - mean * mean, 0.0)
        out[m] = {"mean": mean, "std": float(np.sqrt(var)), "n_pixels": n, "source": "frame_windows_all_pixels_0-1_scale"}
    return out


# --------------------------------------------------------------------------- embeddings
def cache_key(donor: str, region: dict, mapping: pd.DataFrame, model_sha: str) -> str:
    spec = {"donor": str(donor), "window": region["window"], "digest": region.get("digest"),
            "markers": mapping[["marker", "marker_id", "mean", "std"]].round(9).to_dict("records"),
            "patch_px": DEFAULT_PARAMS["patch_px"], "token_overlap": DEFAULT_PARAMS["token_overlap"],
            "model_type": DEFAULT_PARAMS["model_type"], "marker_max_value": DEFAULT_PARAMS["marker_max_value"],
            "mask": DEFAULT_PARAMS["mask"], "model_sha256": model_sha}
    return hashlib.sha256(json.dumps(spec, sort_keys=True).encode()).hexdigest()[:16]


def worker_versions() -> dict:
    r = subprocess.run([str(KRONOS_VENV_PYTHON), str(WORKER), "versions"], capture_output=True, text=True, check=True)
    return json.loads(r.stdout.strip().splitlines()[-1])


def run_worker(cmd: str, job: dict, job_path: pathlib.Path) -> dict:
    job_path.write_text(json.dumps(job, indent=1))
    subprocess.run([str(KRONOS_VENV_PYTHON), str(WORKER), cmd, "--job", str(job_path)], check=True, cwd=str(pins.PKG))
    return json.loads(pathlib.Path(job["result"]).read_text())


def embed_frame(ctx: arm_base.ArmContext, mapping: pd.DataFrame, assets: dict, exchange: pathlib.Path) -> tuple:
    """Embeddings for every ctx cell: ``(E (n, dim) float32 with NaN rows where absent, have, info)``."""
    donor = ctx.donor_id
    regions = frame_regions(donor)
    oids = geojson_object_ids(donor)
    fc = frame_labels(donor, oids)
    fc = fc.loc[fc["in_expression"]]                    # only cells with a measurement row can be ctx cells
    need = sorted(set(int(r) for r in np.unique(ctx.region_idx)))
    model_sha = assets["sha256"]["kronos_vits16_model.pt"]
    max_value = float(DEFAULT_PARAMS["marker_max_value"])
    tf, pages, idx, names, _shape = open_planes(donor, list(ctx.markers))
    info = {"regions": need, "cache_hits": [], "embedded": [], "channel_names": names, "dataset_stats": {}}
    try:
        odd = mapping.loc[mapping["mode"] == "odd_id", "marker"].tolist()
        if odd:
            # one pass over every frame window for the unmatched markers only (16 reads per marker)
            sub_idx = [idx[list(ctx.markers).index(m)] for m in odd]
            stats = dataset_stats(pages, sub_idx, [r["window"] for r in regions], odd, max_value)
            mapping = fill_dataset_stats(mapping, stats)
            info["dataset_stats"] = stats
        keys = {r: cache_key(donor, regions[r], mapping, model_sha) for r in need}
        cdir = CACHE / f"donor_id={donor}"
        cdir.mkdir(parents=True, exist_ok=True)
        todo = [r for r in need if not (cdir / f"w{r:02d}__{keys[r]}.npz").exists()]
        info["cache_hits"] = [r for r in need if r not in todo]
        if todo:
            exchange.mkdir(parents=True, exist_ok=True)
            windows = []
            for r in todo:
                cells = fc.loc[fc["region_idx"] == r]
                labels = cells["label"].to_numpy(np.int64)
                planes = read_planes(pages, idx, regions[r]["window"])
                np.save(exchange / f"w{r:02d}_planes.npy", planes)
                np.save(exchange / f"w{r:02d}_labels.npy", labels)
                windows.append({"idx": r, "planes": str(exchange / f"w{r:02d}_planes.npy"), "mask": str(mask_path(donor, r)),
                                "labels": str(exchange / f"w{r:02d}_labels.npy"), "out": str(exchange / f"w{r:02d}_emb.npy")})
                del planes
            job = {"model_path": str(KRONOS_ASSETS / "kronos_vits16_model.pt"), "model_type": DEFAULT_PARAMS["model_type"],
                   "token_overlap": DEFAULT_PARAMS["token_overlap"], "patch_px": DEFAULT_PARAMS["patch_px"],
                   "batch_size": DEFAULT_PARAMS["batch_size"], "marker_max_value": max_value,
                   "markers": mapping[["marker", "marker_id", "mean", "std"]].to_dict("records"),
                   "windows": windows, "result": str(exchange / "embed_result.json")}
            res = run_worker("embed", job, exchange / "embed_job.json")
            if res.get("random_init"):
                raise RuntimeError("worker reports random_init -- refusing to record")
            info["embedded"] = res
            for r in todo:
                cells = fc.loc[fc["region_idx"] == r]
                emb = np.load(exchange / f"w{r:02d}_emb.npy")
                if emb.shape[0] != len(cells):
                    raise RuntimeError(f"region {r}: {emb.shape[0]} embeddings for {len(cells)} cells")
                np.savez(cdir / f"w{r:02d}__{keys[r]}.npz", emb=emb, labels=cells["label"].to_numpy(np.int64),
                         object_ids=cells["object_id"].to_numpy().astype("S36"))
                (exchange / f"w{r:02d}_planes.npy").unlink(missing_ok=True)
    finally:
        tf.close()
    parts_e, parts_id = [], []
    for r in need:
        z = np.load(cdir / f"w{r:02d}__{keys[r]}.npz")
        parts_e.append(z["emb"]); parts_id.append(z["object_ids"].astype(str))
    E_all = np.concatenate(parts_e); ids_all = np.concatenate(parts_id)
    pos = pd.Index(ids_all).get_indexer(pd.Index(ctx.ids.astype(str)))
    have = pos >= 0
    E = np.full((len(ctx.ids), E_all.shape[1]), np.nan, dtype=np.float32)
    E[have] = E_all[pos[have]]
    info["mapping"] = mapping.to_dict("records")
    info["cache_keys"] = {int(k): v for k, v in keys.items()}
    info["n_with_embedding"] = int(have.sum())
    return E, have, info


# --------------------------------------------------------------------------- annotation and choice
def declared_columns(markers: list) -> dict:
    """label -> column indices into ``ctx.X`` for its declared markers that this run carries."""
    dec = declarations.load()
    out = {}
    for L, spec in dec["labels"].items():
        cols = [markers.index(m) for m in spec["markers"] if m in markers]
        if cols:
            out[L] = cols
    return out


def _cluster_calls(ann: pd.DataFrame, labels: list) -> dict:
    """cluster -> (label, reason, score) from one annotation table."""
    out = {}
    for _, r in ann.iterrows():
        lab, reason = str(r["label"]), str(r.get("reason", "") or "")
        score = np.nan
        if lab in labels:
            score = float(r[f"contrast_{lab}"] - r[f"p95_{lab}"])
        elif lab == AMBIGUOUS:
            top = max(labels, key=lambda L: r[f"contrast_{L}"])
            score = float(r[f"contrast_{top}"] - r[f"p95_{top}"])
        out[r["cluster"]] = (lab, reason, score)
    return out


def annotate_and_choose(Z: np.ndarray, fits: dict, resolutions, seeds, declared: dict, *, primary_seed: int) -> tuple:
    """Annotate every fit on the marker matrix, pick the resolution by seed stability, return the
    primary per-cell ``(label, reason, score)`` arrays and the full details for the sidecar."""
    from sklearn.metrics import adjusted_rand_score
    labels = list(declared)
    per_cell, tables, calls = {}, {}, {}
    for r in resolutions:
        for s in seeds:
            key = kronos_worker.fit_key(r, s)
            ann = gate.null_gated_cluster_annotation(Z, fits[key], declared, seed=int(s), margin_frac=DEFAULT_PARAMS["margin_frac"])
            cc = _cluster_calls(ann, labels)
            lab = np.array([cc[c][0] for c in fits[key]], dtype=object)
            per_cell[key] = lab
            calls[key] = cc
            tables[key] = ann.to_dict("records")
    stability = {}
    for r in resolutions:
        ks = [kronos_worker.fit_key(r, s) for s in seeds]
        aris = [adjusted_rand_score(per_cell[a].astype(str), per_cell[b].astype(str))
                for i, a in enumerate(ks) for b in ks[i + 1:]]
        stability[float(r)] = float(np.mean(aris)) if aris else np.nan
    chosen = sorted(resolutions, key=lambda r: (-(stability[float(r)] if np.isfinite(stability[float(r)]) else -np.inf), float(r)))[0]
    pk = kronos_worker.fit_key(chosen, primary_seed)
    cc = calls[pk]
    lab = np.array([cc[c][0] for c in fits[pk]], dtype=object)
    reason = np.array([cc[c][1] for c in fits[pk]], dtype=object)
    score = np.array([cc[c][2] for c in fits[pk]], dtype=float)
    details = {"stability_mean_pairwise_ari": {f"{k:g}": v for k, v in stability.items()}, "chosen_resolution": float(chosen),
               "primary_fit": pk, "n_clusters": {k: int(len(v)) for k, v in calls.items()},
               "label_counts_per_fit": {k: pd.Series(v.astype(str)).value_counts().to_dict() for k, v in per_cell.items()},
               "annotation_tables": tables}
    return lab, reason, score, details


# --------------------------------------------------------------------------- the arm
def _all(ctx: arm_base.ArmContext, label: str, reason: str) -> pd.DataFrame:
    return pd.DataFrame({"object_id": pd.Series(ctx.ids).astype(str), "label": label, "score": np.nan, "reason": reason})


def _sidecar(ctx: arm_base.ArmContext, payload: dict) -> pathlib.Path:
    d = records.out_dir(ARM, ctx.donor_id, ctx.scope)
    d.mkdir(parents=True, exist_ok=True)
    p = d / f"kronos_details__{ctx.config_id.replace(':', '_')}.json"
    p.write_text(json.dumps(payload, indent=1, default=str))
    return p


def run(ctx: arm_base.ArmContext) -> pd.DataFrame:
    t0 = time.time()
    DEFAULT_PARAMS["run"] = {"donor_id": ctx.donor_id, "config_id": ctx.config_id, "seed": int(ctx.seed),
                             "identity_markers": list(ctx.markers)}
    if ctx.scope != "frame":
        DEFAULT_PARAMS["run"]["deviation"] = REASON_FRAME_ONLY
        return _all(ctx, UNAVAILABLE, REASON_FRAME_ONLY)
    if PURGE_CACHE and CACHE.exists():
        shutil.rmtree(CACHE)
    assets = ensure_assets()
    DEFAULT_PARAMS["run"]["assets"] = assets
    if not assets["ok"]:
        return _all(ctx, UNAVAILABLE, assets["reason"])
    if not KRONOS_VENV_PYTHON.exists():
        raise FileNotFoundError(f"worker interpreter missing: {KRONOS_VENV_PYTHON}")
    DEFAULT_PARAMS["pins"]["worker_versions"] = worker_versions()
    meta = load_marker_metadata(KRONOS_ASSETS / "marker_metadata.csv")
    mapping = map_markers(list(ctx.markers), meta, load_aliases())
    exchange = SCRATCH / f"run__{ctx.donor_id}__{ctx.config_id.replace(':', '_')}__{os.getpid()}"
    try:
        E, have, einfo = embed_frame(ctx, mapping, assets, exchange)
        DEFAULT_PARAMS["run"]["marker_map"] = einfo["mapping"]
        DEFAULT_PARAMS["run"]["dataset_stats"] = einfo["dataset_stats"]
        DEFAULT_PARAMS["run"]["cache_keys"] = einfo["cache_keys"]
        DEFAULT_PARAMS["run"]["n_with_embedding"] = einfo["n_with_embedding"]
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
    declared = declared_columns(list(ctx.markers))
    lab, reason, score, details = annotate_and_choose(Z_all[have], fits, DEFAULT_PARAMS["resolutions"], seeds, declared,
                                                      primary_seed=int(ctx.seed))
    calls = _all(ctx, UNAVAILABLE, REASON_NO_EMBEDDING)
    idx = np.flatnonzero(have)
    calls.loc[idx, "label"] = lab
    calls.loc[idx, "reason"] = reason
    calls.loc[idx, "score"] = score
    DEFAULT_PARAMS["run"].update({"chosen_resolution": details["chosen_resolution"], "primary_fit": details["primary_fit"],
                                  "stability_mean_pairwise_ari": details["stability_mean_pairwise_ari"],
                                  "n_clusters": details["n_clusters"], "graph": {k: v for k, v in cinfo.items() if k != "n_clusters"},
                                  "declared_columns": declared, "seconds": round(time.time() - t0, 1)})
    _sidecar(ctx, {"arm": ARM, "params": DEFAULT_PARAMS, "details": details, "embed": {k: v for k, v in einfo.items() if k != "mapping"}})
    return calls


def _preparse() -> None:
    """Strip this arm's own flags before ``arm_base.cli`` parses the common ones."""
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
