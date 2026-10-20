#!/usr/bin/env python
"""The SOM half of arm B, run INSIDE ``external/.venv-som`` -- never in the repo ``.venv``.

    external/.venv-som/bin/python scripts/broadtype/arms/som_kintsugi_worker.py \
        --features <scratch>/features.parquet --out-dir <scratch>/fits \
        --fits 9:0 9:1 ... 49:4 --num-passes 10 --lr-start 0.06 --lr-end 0.001

Why a separate interpreter: pyFlowSOM requires ``numpy < 2`` and the repo ``.venv`` carries
numpy 2.x, so the SOM cannot import there (``scripts/kintsugi/README.md``).  The driver
(``arms/som_kintsugi.py``) writes the feature block once, this worker fits every requested
(nodes, seed) pair on it, and everything else -- annotation, node choice, the record -- happens
back in ``.venv`` where the gate lives.  Parquet is the only thing that crosses the boundary.

Why it streams: the frame is ~48k cells but a whole donor is ~1.9M x 28.  The input is read in
row batches straight into one preallocated float64 array (no pandas copy of the whole table),
and each fit's cluster vector is written as soon as it exists, so a failure on fit 17 leaves
fits 1-16 on disk and the peak footprint is the array plus one fit's temporaries.  pyFlowSOM's
``som`` needs the full training array in memory by design; that is the one copy this cannot
avoid.

Why ``som_adapter.som_cluster`` and not pyFlowSOM directly: the adapter is KINTSUGI's own
``Kanalysis.visualize_som`` statistical core, asserts the pinned checkout
(``KINTSUGI_COMMIT``) before it imports pyFlowSOM, refuses non-finite input rather than
imputing, and relabels nodes to consecutive integers exactly as the notebook does.  Reusing it
keeps arm B's clustering identical to the kintsugi lane's, so the two can be compared.

Every fit's output carries ``object_id`` beside ``cluster``: the driver joins on the UUID, never
on row position (root ``CLAUDE.md``).
"""
from __future__ import annotations

import argparse
import json
import pathlib
import platform
import sys
import time

import numpy as np
import pyarrow as pa
import pyarrow.parquet as pq

HERE = pathlib.Path(__file__).resolve().parent
SCRIPTS = HERE.parents[1]                      # core/multiplex5000/scripts
sys.path.insert(0, str(SCRIPTS / "kintsugi"))
import som_adapter  # noqa: E402

ID_COLUMN = "object_id"


def parse_fit(text: str) -> tuple[int, int]:
    """``'<nodes>:<seed>'`` -> (nodes, seed); anything else is refused loudly."""
    try:
        nodes, seed = text.split(":")
        return int(nodes), int(seed)
    except ValueError as exc:
        raise SystemExit(f"--fits entries are '<nodes>:<seed>', got {text!r}") from exc


def read_features(path: pathlib.Path, batch_rows: int) -> tuple[np.ndarray, np.ndarray, list[str]]:
    """Stream the driver's parquet into one float64 block plus the object_id vector.

    Row batches are filled in order into a preallocated array so the worker never holds a
    second full-size copy of the block (pyarrow -> pandas -> numpy would).
    """
    pf = pq.ParquetFile(path)
    names = [n for n in pf.schema_arrow.names if n != ID_COLUMN]
    if ID_COLUMN not in pf.schema_arrow.names:
        raise SystemExit(f"{path} has no {ID_COLUMN!r} column; the driver must write the join key")
    n = pf.metadata.num_rows
    values = np.empty((n, len(names)), dtype=np.float64)
    ids = np.empty(n, dtype=object)
    row = 0
    for batch in pf.iter_batches(batch_size=batch_rows, columns=[ID_COLUMN, *names]):
        k = batch.num_rows
        for j, name in enumerate(names):
            values[row:row + k, j] = batch.column(name).to_numpy(zero_copy_only=False)
        ids[row:row + k] = batch.column(ID_COLUMN).to_numpy(zero_copy_only=False)
        row += k
    if row != n:
        raise SystemExit(f"read {row} rows, metadata says {n}")
    return values, ids, names


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description="arm B SOM worker (runs in external/.venv-som)")
    ap.add_argument("--features", required=True, type=pathlib.Path)
    ap.add_argument("--out-dir", required=True, type=pathlib.Path)
    ap.add_argument("--fits", nargs="+", required=True, help="'<nodes>:<seed>' per fit")
    ap.add_argument("--num-passes", type=int, default=som_adapter.NOTEBOOK_SOM["num_passes"])
    ap.add_argument("--lr-start", type=float, default=som_adapter.NOTEBOOK_SOM["lr_start"])
    ap.add_argument("--lr-end", type=float, default=som_adapter.NOTEBOOK_SOM["lr_end"])
    ap.add_argument("--batch-rows", type=int, default=100_000)
    a = ap.parse_args(argv)

    # Fail before touching data if the checkout drifted or pyFlowSOM is missing.
    som_adapter.import_pyflowsom()
    fits = [parse_fit(f) for f in a.fits]
    a.out_dir.mkdir(parents=True, exist_ok=True)

    t0 = time.time()
    values, ids, names = read_features(a.features, a.batch_rows)
    manifest = {
        "provenance": som_adapter.provenance(),
        "python": sys.executable,
        "python_version": platform.python_version(),
        "numpy": np.__version__,
        "pyarrow": pa.__version__,
        "n_cells": int(values.shape[0]),
        "features": names,
        "num_passes": a.num_passes, "lr_start": a.lr_start, "lr_end": a.lr_end,
        "read_seconds": round(time.time() - t0, 2),
        "fits": [],
    }
    for nodes, seed in fits:
        t1 = time.time()
        clusters, _profile = som_adapter.som_cluster(
            values, names, num_nodes=nodes, num_passes=a.num_passes,
            lr_start=a.lr_start, lr_end=a.lr_end, seed=seed)
        out = a.out_dir / f"clusters__n{nodes}__s{seed}.parquet"
        pq.write_table(pa.table({ID_COLUMN: pa.array(ids, type=pa.string()),
                                 "cluster": pa.array(clusters.astype(np.int32))}),
                       out, compression="zstd")
        manifest["fits"].append({"nodes": nodes, "seed": seed, "path": str(out),
                                 "n_occupied": int(len(np.unique(clusters))),
                                 "seconds": round(time.time() - t1, 2)})
        print(f"  fit nodes={nodes:>2} seed={seed:>3}: {manifest['fits'][-1]['n_occupied']:>2} occupied, "
              f"{manifest['fits'][-1]['seconds']:.1f}s", flush=True)
    manifest["total_seconds"] = round(time.time() - t0, 2)
    (a.out_dir / "worker_manifest.json").write_text(json.dumps(manifest, indent=1))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
