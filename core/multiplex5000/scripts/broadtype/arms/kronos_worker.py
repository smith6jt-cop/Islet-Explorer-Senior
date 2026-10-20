"""KRONOS worker: model forward and graph clustering, run in ``external/.venv-kronos``.

Why a separate process.  KRONOS needs torch 2.6.0+cu124 and xformers 0.0.29.post3, which the
lane's ``.venv`` does not carry (its torch is 2.13/cu130) and must never be re-created for.  The
driver (``kronos_clusters.py``, ``.venv``) therefore does everything that touches the cohort's own
contracts -- the frozen frame, the QuPath id invariant, the registry marker list, the qptiff
window reads (``tifffile`` + ``imagecodecs`` live in ``.venv``; this venv has no LZW codec) -- and
hands this process plain ``.npy`` arrays plus a JSON job.  Three subcommands::

    versions            print the pinned package versions as JSON (the driver records them)
    embed   --job J     per window: normalise planes, cut 64-px cell-centred own-mask patches,
                        forward through KRONOS, write the 384-d CLS embedding per requested cell
    cluster --job J     pynndescent kNN + Leiden (RBConfiguration) over resolutions x seeds

Nothing here reads a marker value table, a QuPath class or a label, and nothing here decides a
cell type.  ``torch``/``kronos`` are imported lazily so :func:`extract_patches` and
:func:`cluster_embeddings` import (and are unit-tested) from ``.venv`` without the GPU stack.

Tutorial idiom followed (KRONOS ``tutorials/utils/cell_phenotyping.py`` at 26f1bc02):
  * intensity / 65535 (uint16) -> (x - marker_mean) / marker_std   (``CellPhenotypingDataset``)
  * a 64-px patch around the cell centre, zero-padded at the image edge (``patch_extraction``)
  * patch * binary OWN-cell mask, applied AFTER normalisation     (``feature_extraction``)
  * ``model(patches, marker_ids=ids)`` -> ``(patch_features, marker_features, token_features)``;
    this arm keeps ``patch_features`` -- the CLS token (``x_norm_clstoken``, 384-d for ViT-S) --
    not the per-marker means the tutorial flattens for its logistic regression.
  * ``token_overlap=True`` (stride 8), the tutorial's cell-phenotyping configuration.

The one difference from the tutorial is where the cell centre comes from: the tutorial reads it
from an annotation table, this worker computes it from the cell's own pixels in the frame mask,
so centre and mask can never disagree.
"""
from __future__ import annotations

import argparse
import json
import pathlib
import sys
import time

import numpy as np

EMBED_DIM = {"vits16": 384, "vitl16": 1024}


# --------------------------------------------------------------------------- patches
def label_centroids(mask: np.ndarray, labels: np.ndarray) -> np.ndarray:
    """Integer (row, col) centre of every requested label from its own pixels.  A label absent
    from the mask is a contract violation (the driver asserted presence) and raises."""
    labels = np.asarray(labels)
    want = np.isin(mask, labels)
    ys, xs = np.nonzero(want)
    if ys.size == 0:
        raise ValueError("none of the requested labels is present in the mask")
    lab = mask[ys, xs]
    uniq, inv = np.unique(lab, return_inverse=True)
    cnt = np.bincount(inv).astype(float)
    cy = np.bincount(inv, weights=ys) / cnt
    cx = np.bincount(inv, weights=xs) / cnt
    pos = np.searchsorted(uniq, labels)
    pos = np.clip(pos, 0, len(uniq) - 1)
    missing = uniq[pos] != labels
    if missing.any():
        raise ValueError(f"{int(missing.sum())} requested labels are absent from the mask, e.g. {labels[missing][:5].tolist()}")
    return np.stack([np.floor(cy[pos] + 0.5), np.floor(cx[pos] + 0.5)], axis=1).astype(np.int64)


def iter_patch_batches(planes: np.ndarray, mask: np.ndarray, labels: np.ndarray, *, size: int = 64,
                       batch: int = 128):
    """Yield ``(row_slice, patches)`` with ``patches`` of shape ``(b, M, size, size)`` float32.

    ``planes`` is ``(M, H, W)`` ALREADY normalised; the crop is ``[c - size//2, c + size//2)`` in
    both axes, zero-padded where it leaves the image (the tutorial pads with zeros), and then
    multiplied by ``mask == own label`` so no neighbour's pixels enter the patch.
    """
    M, H, W = planes.shape
    half = size // 2
    centres = label_centroids(mask, labels)
    n = len(labels)
    for s in range(0, n, batch):
        e = min(n, s + batch)
        out = np.zeros((e - s, M, size, size), dtype=np.float32)
        for j in range(s, e):
            cy, cx = int(centres[j, 0]), int(centres[j, 1])
            y0, y1, x0, x1 = cy - half, cy + half, cx - half, cx + half
            sy0, sy1, sx0, sx1 = max(0, y0), min(H, y1), max(0, x0), min(W, x1)
            crop = planes[:, sy0:sy1, sx0:sx1]
            own = (mask[sy0:sy1, sx0:sx1] == labels[j]).astype(np.float32)
            out[j - s, :, sy0 - y0:sy1 - y0, sx0 - x0:sx1 - x0] = crop * own[None]
        yield slice(s, e), out


def extract_patches(planes: np.ndarray, mask: np.ndarray, labels: np.ndarray, *, size: int = 64) -> np.ndarray:
    """All patches at once (tests and small windows); the embed path streams batches instead."""
    parts = [p for _, p in iter_patch_batches(planes, mask, labels, size=size, batch=max(1, len(labels)))]
    return np.concatenate(parts, axis=0) if parts else np.zeros((0, planes.shape[0], size, size), np.float32)


def normalise_planes(planes: np.ndarray, means: np.ndarray, stds: np.ndarray, max_value: float) -> np.ndarray:
    """``(x / max_value - mean) / std`` per marker, float32 -- the tutorial's dataset transform."""
    x = planes.astype(np.float32) / np.float32(max_value)
    m = np.asarray(means, dtype=np.float32)[:, None, None]
    s = np.asarray(stds, dtype=np.float32)[:, None, None]
    if not np.all(np.isfinite(s)) or np.any(s <= 0):
        raise ValueError("every marker_std must be finite and > 0")
    return (x - m) / s


# --------------------------------------------------------------------------- embed
def _load_model(job: dict):
    import torch
    from kronos import create_model, create_model_from_pretrained
    cfg = {"model_type": job["model_type"], "token_overlap": bool(job["token_overlap"])}
    if job.get("random_init"):
        # Test hook only (shape / device check without weights).  The driver never sets it, and a
        # job that does is stamped ``random_init: true`` in its result so it can never pass as one.
        model, precision, dim = create_model(cfg=cfg)
        model.eval()
    else:
        model, precision, dim = create_model_from_pretrained(checkpoint_path=job["model_path"], cfg=cfg)
    device = torch.device(job.get("device") or ("cuda" if torch.cuda.is_available() else "cpu"))
    return model.to(device), device, int(dim)


def embed_job(job: dict) -> dict:
    """Run every window in ``job["windows"]``; write one ``(n_cells, dim)`` float32 ``.npy`` each."""
    import torch
    t0 = time.time()
    model, device, dim = _load_model(job)
    markers = job["markers"]
    means = np.array([m["mean"] for m in markers], dtype=np.float32)
    stds = np.array([m["std"] for m in markers], dtype=np.float32)
    ids = torch.tensor([int(m["marker_id"]) for m in markers], dtype=torch.int64, device=device)
    size, batch, max_value = int(job["patch_px"]), int(job["batch_size"]), float(job["marker_max_value"])
    windows = []
    with torch.no_grad():
        for w in job["windows"]:
            tw = time.time()
            planes = np.load(w["planes"])
            if planes.ndim != 3 or planes.shape[0] != len(markers):
                raise ValueError(f"{w['planes']}: expected ({len(markers)}, H, W), got {planes.shape}")
            mask = np.load(w["mask"])
            labels = np.load(w["labels"]).astype(np.int64)
            planes = normalise_planes(planes, means, stds, max_value)
            out = np.empty((len(labels), dim), dtype=np.float32)
            for sl, patches in iter_patch_batches(planes, mask, labels, size=size, batch=batch):
                x = torch.from_numpy(patches).to(device=device, dtype=torch.float32)
                cls, _marker_feats, _tokens = model(x, marker_ids=ids.unsqueeze(0).expand(x.shape[0], -1))
                out[sl] = cls.float().cpu().numpy()
            if not np.all(np.isfinite(out)):
                raise FloatingPointError(f"non-finite embedding in window {w['idx']}")
            np.save(w["out"], out)
            windows.append({"idx": w["idx"], "n_cells": int(len(labels)), "seconds": round(time.time() - tw, 2)})
            del planes
    res = {"ok": True, "dim": dim, "device": str(device), "random_init": bool(job.get("random_init", False)),
           "torch": torch.__version__, "windows": windows, "seconds": round(time.time() - t0, 2)}
    if job.get("result"):
        pathlib.Path(job["result"]).write_text(json.dumps(res, indent=1))
    return res


# --------------------------------------------------------------------------- cluster
def knn_edges(E: np.ndarray, *, k: int = 15, seed: int = 0, metric: str = "euclidean") -> np.ndarray:
    """Symmetrised, de-duplicated, self-free kNN edge list ``(n_edges, 2)`` from pynndescent."""
    import pynndescent
    n = E.shape[0]
    if n <= k + 1:
        raise ValueError(f"{n} cells cannot carry a k={k} neighbour graph")
    nnd = pynndescent.NNDescent(np.ascontiguousarray(E, dtype=np.float32), n_neighbors=k + 1, metric=metric,
                                random_state=seed, n_jobs=-1)
    idx, _dist = nnd.neighbor_graph
    rows = np.repeat(np.arange(n), idx.shape[1])
    cols = idx.ravel().astype(np.int64)
    keep = (cols >= 0) & (rows != cols)
    edges = np.stack([rows[keep], cols[keep]], axis=1)
    edges.sort(axis=1)
    return np.unique(edges, axis=0)


def cluster_embeddings(E: np.ndarray, *, resolutions=(0.3, 0.5, 1.0), seeds=(0, 1, 2, 3, 4), k: int = 15,
                       knn_seed: int = 0, metric: str = "euclidean", n_iterations: int = -1) -> tuple:
    """Leiden (RBConfiguration) on one kNN graph at every resolution x seed.

    Returns ``(fits, info)``: ``fits[f"r{res}_s{seed}"]`` is an int32 membership per row of ``E``;
    ``info`` records graph size and per-fit cluster counts.  The number of clusters is whatever
    Leiden returns -- nothing here targets a count.
    """
    import igraph as ig
    import leidenalg as la
    t0 = time.time()
    edges = knn_edges(E, k=k, seed=knn_seed, metric=metric)
    g = ig.Graph(n=E.shape[0], edges=edges.tolist(), directed=False)
    fits, n_clusters = {}, {}
    for r in resolutions:
        for s in seeds:
            part = la.find_partition(g, la.RBConfigurationVertexPartition, resolution_parameter=float(r),
                                     seed=int(s), n_iterations=int(n_iterations))
            key = fit_key(r, s)
            fits[key] = np.asarray(part.membership, dtype=np.int32)
            n_clusters[key] = int(fits[key].max() + 1)
    info = {"n_cells": int(E.shape[0]), "n_edges": int(len(edges)), "k": k, "metric": metric, "knn_seed": knn_seed,
            "n_iterations": n_iterations, "n_clusters": n_clusters, "seconds": round(time.time() - t0, 2)}
    return fits, info


def fit_key(resolution: float, seed: int) -> str:
    return f"r{float(resolution):g}_s{int(seed)}"


def cluster_job(job: dict) -> dict:
    E = np.load(job["embeddings"])
    fits, info = cluster_embeddings(E, resolutions=tuple(job["resolutions"]), seeds=tuple(job["seeds"]),
                                    k=int(job["k"]), knn_seed=int(job["knn_seed"]), metric=job.get("metric", "euclidean"),
                                    n_iterations=int(job.get("n_iterations", -1)))
    np.savez(job["out"], **fits)
    info["ok"] = True
    if job.get("result"):
        pathlib.Path(job["result"]).write_text(json.dumps(info, indent=1))
    return info


# --------------------------------------------------------------------------- versions
def versions() -> dict:
    out = {"python": sys.version.split()[0], "numpy": np.__version__}
    for name in ("torch", "xformers", "kronos", "huggingface_hub", "pynndescent", "leidenalg", "igraph", "sklearn"):
        try:
            mod = __import__(name)
            out[name] = str(getattr(mod, "__version__", "installed"))
        except Exception as exc:  # noqa: BLE001
            out[name] = f"missing ({type(exc).__name__})"
    try:
        from importlib.metadata import version as _v
        out["kronos_dist"] = _v("KRONOS")
    except Exception:  # noqa: BLE001
        pass
    try:
        import torch
        out["cuda_available"] = bool(torch.cuda.is_available())
        out["cuda_device"] = torch.cuda.get_device_name(0) if torch.cuda.is_available() else None
    except Exception:  # noqa: BLE001
        pass
    return out


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description="KRONOS worker (external/.venv-kronos)")
    sub = ap.add_subparsers(dest="cmd", required=True)
    sub.add_parser("versions")
    for name in ("embed", "cluster"):
        p = sub.add_parser(name)
        p.add_argument("--job", required=True)
    a = ap.parse_args(argv)
    if a.cmd == "versions":
        print(json.dumps(versions()))
        return 0
    job = json.loads(pathlib.Path(a.job).read_text())
    res = embed_job(job) if a.cmd == "embed" else cluster_job(job)
    print(json.dumps({k: v for k, v in res.items() if k not in ("windows", "n_clusters")}))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
