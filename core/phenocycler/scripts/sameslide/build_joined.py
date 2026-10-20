#!/usr/bin/env python3
"""
Joined per-cell protein + transcript table for a same-slide pair.

One row per Xenium cell. Protein is transferred from PhenoCycler cells by AREA-WEIGHTED
POLYGON OVERLAP in the registered frame: each PhenoCycler cell contributes to each Xenium
cell in proportion to the area they share. That neither forces a bijection nor throws away
extent, which matters because the two segmentations disagree (1,028,140 vs 717,461 cells).

``overlap_fraction`` is the per-row honesty column: a cell assembled from a sliver of one
protein cell is not the same evidence as one covered by a well-matched partner.

Genes go to a sparse sidecar rather than 5,001 parquet columns -- dense would be ~14 GB and
every contract validation would walk all of it.
"""

from __future__ import annotations

import argparse
import json
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd
import shapely

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))

from phenocycler.integration.contract import (  # noqa: E402
    CONTRACT_PROFILE_GEOMETRY, coerce_cell_table, feature_columns, write_cell_table,
)
from phenocycler.integration.transform import Transform  # noqa: E402


def require_unique_ids(table: pd.DataFrame, key: str) -> None:
    if table[key].isna().any() or table[key].astype(str).eq("").any():
        raise ValueError(f"null/empty {key} is not a QuPath UUID")
    if table[key].astype(str).duplicated().any():
        raise ValueError(f"duplicate {key}; UUID joins must be one-to-one")


def overlap_edges(xids, xpolys, pheno_polygons: pd.DataFrame, transform):
    """Positive-area edges, ordered by Xenium ID then immutable QuPath UUID."""
    key = "object_id" if "object_id" in pheno_polygons else "cell_id"
    require_unique_ids(pheno_polygons, key)
    xw = shapely.transform(xpolys, lambda a: transform.apply(a))
    pp = pheno_polygons.sort_values(key, kind="stable").reset_index(drop=True)
    ppolys = shapely.from_wkb(pp["wkb"].to_numpy())
    for arr in (xw, ppolys):
        bad = ~shapely.is_valid(arr)
        arr[bad] = shapely.make_valid(arr[bad])
    xi, pi = shapely.STRtree(ppolys).query(xw, predicate="intersects")
    area = shapely.area(shapely.intersection(xw[xi], ppolys[pi]))
    keep = np.isfinite(area) & (area > 0)
    edges = pd.DataFrame({"cell_id": np.asarray(xids)[xi[keep]],
                          "object_id": pp[key].astype(str).to_numpy()[pi[keep]],
                          "overlap_area_um2": area[keep]})
    edges = edges.sort_values(["cell_id", "object_id"], kind="stable").reset_index(drop=True)
    return edges, shapely.area(xw)


def aggregate_features(xids, edges: pd.DataFrame, pheno: pd.DataFrame, *, include_denominators=False):
    """Finite measurements have their own area denominator; missing stays missing."""
    key = "object_id" if "object_id" in pheno else "cell_id"
    require_unique_ids(pheno, key)
    p = pheno.assign(**{key: pheno[key].astype(str)}).set_index(key)
    e = edges[edges.object_id.isin(p.index)]
    xi = pd.Index(xids).get_indexer(e.cell_id)
    if (xi < 0).any():
        raise ValueError("overlap edges contain unknown Xenium cells")
    weights = e.overlap_area_um2.to_numpy(float)
    out = {"cell_id": np.asarray(xids)}
    for feature in feature_columns(pheno):
        values = p[feature].reindex(e.object_id).to_numpy(float)
        finite = np.isfinite(values)
        numerator = np.bincount(xi[finite], weights=weights[finite] * values[finite], minlength=len(xids))
        denominator = np.bincount(xi[finite], weights=weights[finite], minlength=len(xids))
        out[feature] = np.divide(numerator, denominator, out=np.full(len(xids), np.nan),
                                 where=denominator > 0).astype(np.float32)
        if include_denominators:
            out["support_area_um2__" + feature] = denominator.astype(np.float32)
    return pd.DataFrame(out), e


def xenium_polygons(bundle: Path) -> tuple[np.ndarray, np.ndarray]:
    """Xenium cell outlines (already microns), one polygon per cell_id."""
    b = pd.read_parquet(bundle / "cell_boundaries.parquet")
    b = b.sort_values("cell_id", kind="stable")
    ids, starts = np.unique(b["cell_id"].to_numpy(), return_index=True)
    order = np.argsort(starts)
    ids, starts = ids[order], starts[order]
    coords = b[["vertex_x", "vertex_y"]].to_numpy(dtype=float)
    idx = np.repeat(np.arange(len(starts)), np.diff(np.append(starts, len(b))))
    polys = shapely.polygons(shapely.linearrings(coords, indices=idx))
    return ids, polys


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--pheno", required=True, type=Path)
    ap.add_argument("--pheno-polygons", required=True, type=Path)
    ap.add_argument("--bundle", required=True, type=Path)
    ap.add_argument("--transform-dir", required=True, type=Path)
    ap.add_argument("--transform-stem", required=True)
    ap.add_argument("--out", required=True, type=Path)
    ap.add_argument("--donor", required=True)
    ap.add_argument("--roi", default="panc")
    ap.add_argument("--serial", required=True)
    ap.add_argument("--tissue", default="pancreas")
    args = ap.parse_args(argv)

    t0 = time.time()
    tf = Transform.load(args.transform_dir, stem=args.transform_stem)

    xids, xpolys = xenium_polygons(args.bundle)
    print(f"[joined] xenium polygons {len(xids):,} ({time.time()-t0:.0f}s)", flush=True)
    pp = pd.read_parquet(args.pheno_polygons)
    pheno = pd.read_parquet(args.pheno)
    edges, xarea = overlap_edges(xids, xpolys, pp, tf)
    df, eligible_edges = aggregate_features(xids, edges, pheno)
    xi = pd.Index(xids).get_indexer(eligible_edges.cell_id)
    inter = eligible_edges.overlap_area_um2.to_numpy(float)
    n = len(xids)
    wsum = np.bincount(xi, weights=inter, minlength=n)
    feats = feature_columns(pheno)
    print(f"[joined] aggregated {len(feats)} features ({time.time()-t0:.0f}s)", flush=True)

    cent = shapely.centroid(xpolys)
    df["x_um"] = shapely.get_x(cent)
    df["y_um"] = shapely.get_y(cent)
    df["area_um2"] = shapely.area(xpolys)
    df["n_pheno_cells"] = np.bincount(xi, minlength=n)
    df["overlap_fraction"] = np.where(xarea > 0, wsum / np.maximum(xarea, 1e-12), 0.0)
    df["xenium_sample_key"] = f"{args.serial}__{args.roi}"

    df = coerce_cell_table(df, modality="xenium", donor_id=args.donor,
                           profile=CONTRACT_PROFILE_GEOMETRY, roi=args.roi, tissue=args.tissue)
    args.out.mkdir(parents=True, exist_ok=True)
    path = write_cell_table(df, args.out / "joined_cells.parquet")
    edges.to_parquet(args.out / "overlap_edges.parquet", index=False)

    meta = {"n_xenium_cells": int(n), "n_pheno_cells": int(len(pp)),
            "n_overlapping_pairs": int(len(xi)),
            "transform_stem": args.transform_stem,
            "n_protein_features": len(feats),
            "cells_with_protein": int((df["n_pheno_cells"] > 0).sum()),
            "overlap_fraction_median": float(np.median(df["overlap_fraction"])),
            "n_pheno_cells_median": float(np.median(df["n_pheno_cells"])),
            "seconds": round(time.time() - t0, 1)}
    (args.out / "build_joined.json").write_text(json.dumps(meta, indent=2) + "\n")
    print(f"[joined] wrote {path}")
    print(json.dumps(meta, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
