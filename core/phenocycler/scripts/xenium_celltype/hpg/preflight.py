#!/usr/bin/env python
"""Resolve every cohort donor at this site, size it, and emit the work manifest.

Nothing else in the run discovers work. Every later stage reads ``work_manifest.json`` and
does exactly what is listed there, so "did the run cover the cohort?" is answered by
comparing the manifest to the cohort, not by trusting that a glob found everything.

Fail-closed in three specific ways, each of which has a way of going wrong quietly:

* a donor that resolves at no site is an ERROR, not a skip (``sites.py`` raises);
* a donor centrally excluded in ``phenocycler.cohort`` is refused even if a path exists;
* the manifest is rejected unless every donor in ``qupath_manifest.json['expected_donors']``
  appears in it. A 19-donor run that exits 0 is the failure this guards.

Cell counts come from the file header, not from loading the matrix, so this is seconds even
though it decides the shape of ~1e5 core-hours of work downstream.
"""

from __future__ import annotations

import argparse
import json
import sys
from datetime import datetime, timezone
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(HERE.parents[2]))          # core/phenocycler, for `phenocycler`

import sites  # noqa: E402


def n_cells_of(r: "sites.Resolved") -> tuple[int, int]:
    """(cells in the expression source, cells in cells.parquet) — from headers only."""
    import h5py
    import pyarrow.parquet as pq

    n_spatial = pq.ParquetFile(r.cells_parquet).metadata.num_rows
    if r.expression_source == "h5ad":
        with h5py.File(r.h5ad, "r") as f:
            return int(f["obs/_index"].shape[0]), int(n_spatial)
    with h5py.File(r.bundle / "cell_feature_matrix.h5", "r") as f:
        return int(f["matrix/shape"][1]), int(n_spatial)


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--roi", default="panc")
    ap.add_argument("--site", default=None, help="key in sites.json; default = autodetect")
    ap.add_argument("--repo", type=Path, default=HERE.parents[2],
                    help="core/phenocycler")
    ap.add_argument("--out", type=Path, required=True, help="run directory")
    ap.add_argument("--arms", nargs="+", default=["nonspatial", "spatial"])
    ap.add_argument("--M", type=int, default=96, help="total KODAMA Monte-Carlo iterations")
    ap.add_argument("--m-shard", type=int, default=32, help="iterations per shard task")
    ap.add_argument("--landmarks", type=int, default=2000)
    ap.add_argument("--n-comps", type=int, default=50)
    ap.add_argument("--seed", type=int, default=0)
    a = ap.parse_args(argv)

    if a.M % a.m_shard:
        raise SystemExit(f"--M {a.M} must be a whole multiple of --m-shard {a.m_shard}; "
                         f"a partial shard would silently change the iteration set")

    from phenocycler.cohort import DONOR_EXCLUSIONS
    cohort = sorted(json.loads(
        (a.repo / "data" / "qupath_manifest.json").read_text())["expected_donors"])
    blocked = [d for d in cohort if d in DONOR_EXCLUSIONS]
    if blocked:
        raise SystemExit(f"cohort manifest lists centrally excluded donor(s): {blocked}")

    xpaths = a.repo / "data" / "integration" / "xenium_paths.csv"
    donors, failures = [], []
    for d in cohort:
        try:
            r = sites.resolve(d, xpaths, roi=a.roi, site=a.site)
            n_expr, n_spatial = n_cells_of(r)
        except Exception as exc:                        # noqa: BLE001 — reported, not swallowed
            failures.append(f"{d}: {exc}")
            continue
        donors.append({
            "donor_id": d, "roi": a.roi, "site": r.site,
            "bundle": str(r.bundle), "h5ad": str(r.h5ad) if r.h5ad else None,
            "expression_source": r.expression_source,
            "n_cells_expression": n_expr, "n_cells_spatial": n_spatial,
        })

    if failures:
        print("PREFLIGHT FAILED — these donors did not resolve:", file=sys.stderr)
        for f in failures:
            print(f"  {f}", file=sys.stderr)
        return 2
    missing = set(cohort) - {d["donor_id"] for d in donors}
    if missing:
        print(f"PREFLIGHT FAILED — cohort donors absent from manifest: {sorted(missing)}",
              file=sys.stderr)
        return 2

    n_shards = a.M // a.m_shard
    tasks = [{"donor_id": d["donor_id"], "roi": d["roi"], "arm": arm, "shard": s,
              "seed": a.seed + a.m_shard * s, "m": a.m_shard,
              "n_cells": d["n_cells_expression"]}
             for d in donors for arm in a.arms for s in range(n_shards)]

    a.out.mkdir(parents=True, exist_ok=True)
    manifest = {
        "created_at_utc": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "roi": a.roi, "arms": a.arms,
        "kodama": {"M": a.M, "m_shard": a.m_shard, "n_shards": n_shards,
                   "landmarks": a.landmarks, "seed": a.seed},
        "pca": {"n_comps": a.n_comps},
        "donors": donors, "shard_tasks": tasks,
        "n_donors": len(donors), "n_shard_tasks": len(tasks),
        "n_cells_total": sum(d["n_cells_expression"] for d in donors),
    }
    (a.out / "work_manifest.json").write_text(json.dumps(manifest, indent=2))

    print(f"[preflight] {len(donors)} donors, {manifest['n_cells_total']:,} cells")
    for d in sorted(donors, key=lambda x: -x["n_cells_expression"]):
        drop = d["n_cells_spatial"] - d["n_cells_expression"]
        print(f"  {d['donor_id']:5s} {d['n_cells_expression']:>9,} cells  "
              f"src={d['expression_source']:<7s} site={d['site']:<11s} "
              f"(cells.parquet has {drop:+,} more)")
    print(f"[preflight] arms={a.arms} M={a.M} = {n_shards} x m={a.m_shard} "
          f"-> {len(tasks)} shard tasks; array index 0-{len(tasks)-1}")
    print(f"[preflight] wrote {a.out / 'work_manifest.json'}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
