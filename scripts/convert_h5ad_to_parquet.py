#!/usr/bin/env python3
"""
Convert the enriched Islet Explorer H5AD to Parquet for DuckDB-backed queries.

This is Phase 1.1 of the scaling work. The app currently loads every islet into
an in-memory R dataframe at startup (via ``reticulate`` + ``anndata``), which
works at the 5,214-islet scale but will exhaust per-session RAM when Phenocycler
/ Xenium datasets push us into the millions of rows. Parquet + DuckDB lets us
swap that in-memory table for a lazy view without rewriting the Shiny modules.

The conversion produces:

    data/parquet/islets.parquet              -- main islet-level table (merge
                                                 of .obs, X_scVI_mean, X_umap)
    data/parquet/uns_targets.parquet         -- from .uns (was Islet_Targets)
    data/parquet/uns_markers.parquet         -- from .uns (was Islet_Markers)
    data/parquet/uns_composition.parquet     -- from .uns (was Islet_Composition)
    data/parquet/tissue/case_id=<id>/part-0.parquet  -- partitioned per-donor
                                                       tissue cells (from CSVs)

Design notes
------------
* Column names are preserved EXACTLY as they appear in the H5AD ``.obs`` and the
  reconstructed groovy dataframes. The R layer depends on unchanged names.
* ``.uns`` stores each groovy column under the key
  ``groovy_<sheet>_<column>``. The logic here mirrors
  ``reconstruct_groovy_df_from_list()`` in ``R/data_loading.R``.
* Obsm arrays (``X_scVI_mean``, ``X_umap``) are joined to ``.obs`` before
  writing, so downstream code can query a single table.
* ``tissue/`` is a partitioned dataset keyed on ``case_id`` so DuckDB can prune
  irrelevant partitions when the Spatial tab selects a donor.

Usage
-----
    conda activate scvi-env
    python scripts/convert_h5ad_to_parquet.py
    python scripts/convert_h5ad_to_parquet.py \
        --h5ad data/islet_explorer.h5ad \
        --donors-dir data/donors \
        --output-dir data/parquet
"""
from __future__ import annotations

import argparse
import os
import sys
from pathlib import Path
from typing import Optional

import anndata as ad
import numpy as np
import pandas as pd
import pyarrow as pa
import pyarrow.parquet as pq


# ── Defaults ───────────────────────────────────────────────────────────────

SCRIPT_DIR = Path(__file__).resolve().parent
REPO_ROOT = SCRIPT_DIR.parent
DEFAULT_H5AD = REPO_ROOT / "data" / "islet_explorer.h5ad"
DEFAULT_DONORS_DIR = REPO_ROOT / "data" / "donors"
DEFAULT_OUTPUT_DIR = REPO_ROOT / "data" / "parquet"


# ── Groovy reconstruction (mirrors R/data_loading.R) ───────────────────────

GROOVY_SHEETS = {
    "targets": "uns_targets.parquet",
    "markers": "uns_markers.parquet",
    "composition": "uns_composition.parquet",
}


def reconstruct_groovy_df(uns: dict, sheet: str) -> Optional[pd.DataFrame]:
    """Rebuild a groovy DataFrame from the .uns key-value storage."""
    cols_key = f"groovy_{sheet}_columns"
    nrows_key = f"groovy_{sheet}_n_rows"
    prefix = f"groovy_{sheet}_"

    if cols_key not in uns:
        return None

    cols = list(uns[cols_key])
    n_rows = int(uns.get(nrows_key, 0))
    if n_rows == 0:
        return None

    data = {}
    for col in cols:
        key = prefix + col
        vals = uns.get(key)
        if vals is None:
            continue
        arr = np.asarray(vals)
        # h5ad stores strings as bytes sometimes; normalise to str.
        if arr.dtype.kind in ("S", "O"):
            arr = np.array([
                v.decode("utf-8") if isinstance(v, (bytes, bytearray)) else v
                for v in arr.tolist()
            ], dtype=object)
        data[col] = arr

    if not data:
        return None

    df = pd.DataFrame(data)

    # Standardise "Case ID" to int (matches Excel convention and R-side parsing)
    if "Case ID" in df.columns:
        df["Case ID"] = pd.to_numeric(df["Case ID"], errors="coerce").astype("Int64")

    return df


# ── Obs + obsm flattening ──────────────────────────────────────────────────

def build_islet_table(adata: ad.AnnData) -> pd.DataFrame:
    """Flatten .obs and selected .obsm arrays into a single islet-level table."""
    obs = adata.obs.copy()

    # Promote the obs index so we can treat it like a column in DuckDB / dbplyr.
    obs_index_name = obs.index.name or "islet_id"
    obs = obs.reset_index().rename(columns={"index": obs_index_name})
    if obs_index_name not in obs.columns:
        # reset_index sometimes uses the original name if it existed
        obs.insert(0, obs_index_name, adata.obs.index.astype(str))

    # Attach obsm arrays. Only the columns the Shiny app currently cares about
    # are materialised here. Add more as new panels require them.
    obsm_specs = {
        "X_umap": ("umap_1", "umap_2"),
        "X_scVI_mean": None,  # expanded to scvi_1..scvi_N below
    }
    for key, names in obsm_specs.items():
        if key not in adata.obsm:
            continue
        arr = np.asarray(adata.obsm[key])
        if arr.ndim != 2:
            continue
        if names is None:
            names = tuple(f"scvi_{i+1}" for i in range(arr.shape[1]))
        for i, col in enumerate(names[: arr.shape[1]]):
            if col in obs.columns:
                # Don't overwrite an existing obs column with the same name.
                col = f"{col}_obsm"
            obs[col] = arr[:, i]

    # Normalise any categorical columns to plain Python strings: pyarrow/parquet
    # round-trips strings cleanly, categoricals can lose their dictionary order
    # when consumed by DuckDB.
    for col in obs.columns:
        if pd.api.types.is_categorical_dtype(obs[col]):
            obs[col] = obs[col].astype(str)

    return obs


# ── Tissue CSV partitioning ────────────────────────────────────────────────

def infer_case_id(filename: str) -> Optional[str]:
    """Extract the donor/case id from a tissue CSV filename like ``6505.csv``."""
    stem = Path(filename).stem
    digits = "".join(ch for ch in stem if ch.isdigit())
    return digits or None


def write_tissue_partitions(donors_dir: Path, output_dir: Path) -> int:
    """Concatenate per-donor tissue CSVs into a partitioned Parquet dataset."""
    if not donors_dir.exists():
        print(f"[tissue] No donors directory at {donors_dir}; skipping")
        return 0

    csvs = sorted(donors_dir.glob("*.csv"))
    if not csvs:
        print(f"[tissue] No CSVs under {donors_dir}; skipping")
        return 0

    tissue_root = output_dir / "tissue"
    tissue_root.mkdir(parents=True, exist_ok=True)

    n_written = 0
    for csv in csvs:
        case_id = infer_case_id(csv.name)
        if case_id is None:
            print(f"[tissue] Could not infer case_id from {csv.name}; skipping")
            continue
        df = pd.read_csv(csv)
        if df.empty:
            continue
        df["case_id"] = case_id

        # Partition by case_id so DuckDB can prune when the Spatial tab
        # queries a single donor.
        pq.write_to_dataset(
            pa.Table.from_pandas(df, preserve_index=False),
            root_path=str(tissue_root),
            partition_cols=["case_id"],
        )
        n_written += len(df)
        print(f"[tissue] {csv.name}: wrote {len(df):,} rows (case_id={case_id})")

    print(f"[tissue] Total rows written: {n_written:,}")
    return n_written


# ── Top-level conversion ───────────────────────────────────────────────────

def convert(h5ad_path: Path, donors_dir: Path, output_dir: Path) -> None:
    if not h5ad_path.exists():
        raise FileNotFoundError(f"H5AD not found: {h5ad_path}")

    output_dir.mkdir(parents=True, exist_ok=True)

    print(f"[load] Reading {h5ad_path}")
    adata = ad.read_h5ad(h5ad_path)
    print(f"[load] shape={adata.shape} obs_cols={len(adata.obs.columns)}")

    # 1. Islet-level main table
    islets = build_islet_table(adata)
    islets_path = output_dir / "islets.parquet"
    islets.to_parquet(islets_path, index=False)
    print(f"[write] islets.parquet rows={len(islets):,} cols={len(islets.columns)}")

    # 2. Groovy .uns tables
    uns = dict(adata.uns)  # shallow copy for dict-like access
    for sheet, fname in GROOVY_SHEETS.items():
        df = reconstruct_groovy_df(uns, sheet)
        if df is None or df.empty:
            print(f"[write] {fname} -- no data in .uns for sheet '{sheet}', skipping")
            continue
        path = output_dir / fname
        df.to_parquet(path, index=False)
        print(f"[write] {fname} rows={len(df):,} cols={len(df.columns)}")

    # 3. Tissue partitions (per-donor)
    write_tissue_partitions(donors_dir, output_dir)

    print(f"\n[done] Parquet dataset at {output_dir}")


def main(argv: Optional[list] = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--h5ad", type=Path, default=DEFAULT_H5AD,
                        help=f"Path to islet_explorer.h5ad (default: {DEFAULT_H5AD})")
    parser.add_argument("--donors-dir", type=Path, default=DEFAULT_DONORS_DIR,
                        help=f"Directory of per-donor tissue CSVs (default: {DEFAULT_DONORS_DIR})")
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT_DIR,
                        help=f"Output Parquet root (default: {DEFAULT_OUTPUT_DIR})")
    args = parser.parse_args(argv)

    try:
        convert(args.h5ad, args.donors_dir, args.output_dir)
    except FileNotFoundError as e:
        print(f"ERROR: {e}", file=sys.stderr)
        return 2
    except Exception as e:
        print(f"ERROR: conversion failed: {e}", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
