#!/usr/bin/env python3
"""LEGACY LANE, NOT RUNNABLE HERE -- recorded 2026-09-18.

``scripts/retire/dispositions.py`` records this file as a verified orphan that is also the
NAMED PRODUCER of a live artifact under ``data/app_data`` (code tier 2), which is why it is
kept rather than deleted. It cannot be re-run: the inputs named below are absent from this
checkout.

The current PhenoCycler typing is the subpanel classification under
``core/phenocycler/scripts/typing/`` -- context from QuPath objects, then one question per
context against reference levels drawn from the ordinary tissue, fitted PER DONOR.
``data/development/typing/cells/`` holds 25,230,109 cells over 20 donors (2026-09-15 build,
eight classes plus ``Excluded``). Nothing below describes that method, and ``data/app_data``
is the OTHER app domain: the two datasets share five donors and no vocabulary.

Extract per-islet single-cell data for the Shiny drill-down viewer.

For each qualified islet (≥20 core cells), outputs a CSV file containing all
associated cells (core + peri-islet) with spatial coordinates, phenotype,
region label, cell morphology, and 31 protein marker expressions.

Output: data/cells/{imageid}_Islet_{N}.csv (~1,024 files)

Usage:
    python scripts/extract_per_islet_cells.py
    python scripts/extract_per_islet_cells.py --input path/to/single_cell.h5ad --output-dir data/cells
"""

import argparse
import os
import re
import sys
import warnings

import numpy as np
import pandas as pd

warnings.filterwarnings("ignore", category=FutureWarning)

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
PROJECT_ROOT = os.path.join(SCRIPT_DIR, "..")

DEFAULT_SC_H5AD = os.path.join(
    PROJECT_ROOT, "single_cell_analysis",
    "CODEX_scvi_BioCov_phenotyped_newDuctal.h5ad"
)
DEFAULT_ISLET_H5AD = os.path.join(
    PROJECT_ROOT, "islet_analysis", "islets_core_fixed.h5ad"
)
DEFAULT_OUTPUT_DIR = os.path.join(PROJECT_ROOT, "data", "cells")


def parse_parent(parent_str):
    """Parse Parent column to extract islet name and region.
    Returns (islet_name, region) or (None, None) for non-islet cells.
    """
    if not isinstance(parent_str, str):
        return None, None
    m = re.match(r"^(Islet_\d+)(_exp20um)?$", parent_str)
    if m:
        islet_name = m.group(1)
        region = "peri" if m.group(2) else "core"
        return islet_name, region
    return None, None


def extract_cells(sc_path, islet_path, output_dir):
    """Extract per-islet cell CSVs from single-cell H5AD."""
    import anndata as ad

    print("=" * 60)
    print("Extracting Per-Islet Cell Data")
    print("=" * 60)

    # ── 1. Load qualified islet list ──────────────────────────────────────
    print(f"\n1. Loading qualified islets: {islet_path}")
    islet_adata = ad.read_h5ad(islet_path)
    qualified = islet_adata.obs[["imageid", "base_islet_id"]].copy()
    qualified["imageid"] = qualified["imageid"].astype(str)
    qualified["base_islet_id"] = qualified["base_islet_id"].astype(str)
    qualified["combined_islet_id"] = (
        qualified["imageid"] + "_" + qualified["base_islet_id"]
    )
    qualified_set = set(qualified["combined_islet_id"])
    print(f"   Qualified islets: {len(qualified_set)}")

    # ── 2. Load single-cell H5AD ─────────────────────────────────────────
    print(f"\n2. Loading single-cell H5AD: {sc_path}")
    sc = ad.read_h5ad(sc_path, backed="r")
    var_names = list(sc.var_names)
    print(f"   Total cells: {sc.shape[0]:,}")
    print(f"   Markers ({len(var_names)}): {', '.join(var_names[:10])}...")

    # ── 3. Parse cell-islet assignments ──────────────────────────────────
    print("\n3. Parsing cell-islet assignments...")
    obs = sc.obs[["imageid", "Parent", "phenotype",
                  "X_centroid", "Y_centroid",
                  "Cell Area", "Nucleus Area"]].copy()
    obs["imageid"] = obs["imageid"].astype(str)
    obs["Parent"] = obs["Parent"].astype(str)
    obs["phenotype"] = obs["phenotype"].astype(str)

    parsed = obs["Parent"].apply(parse_parent)
    obs["islet_name"] = [p[0] for p in parsed]
    obs["cell_region"] = [p[1] for p in parsed]

    # Filter to islet cells only
    islet_mask = obs["islet_name"].notna()
    obs_islet = obs[islet_mask].copy()
    obs_islet["combined_islet_id"] = (
        obs_islet["imageid"] + "_" + obs_islet["islet_name"]
    )

    # Filter to qualified islets
    obs_islet = obs_islet[obs_islet["combined_islet_id"].isin(qualified_set)].copy()
    print(f"   Cells in qualified islets: {len(obs_islet):,}")

    # ── 4. Extract expression matrix for these cells ─────────────────────
    print("\n4. Extracting expression data (this may take a minute)...")
    # Get the integer indices of the islet cells
    cell_indices = np.where(islet_mask.values)[0]
    # Further filter to qualified
    qualified_mask = obs.loc[islet_mask, "combined_islet_id_temp" if False else "imageid"].index  # use obs_islet index
    cell_idx_list = obs_islet.index

    # Read X for these cells - chunked to manage memory
    # Since backed mode, slice to get the expression
    chunk_size = 50000
    idx_array = np.array([sc.obs_names.get_loc(i) for i in cell_idx_list])
    idx_array.sort()

    expr_chunks = []
    for start in range(0, len(idx_array), chunk_size):
        end = min(start + chunk_size, len(idx_array))
        chunk_idx = idx_array[start:end]
        chunk_data = sc.X[chunk_idx, :].toarray() if hasattr(sc.X[chunk_idx, :], 'toarray') else np.array(sc.X[chunk_idx, :])
        expr_chunks.append(chunk_data)
        print(f"   Read chunk {start//chunk_size + 1}: cells {start}-{end}")

    expr_matrix = np.vstack(expr_chunks)

    # Build expression DataFrame aligned with obs_islet
    # Re-sort obs_islet to match idx_array order
    idx_to_pos = {idx: pos for pos, idx in enumerate(idx_array)}
    obs_islet_sorted = obs_islet.loc[[sc.obs_names[i] for i in idx_array]].copy()

    expr_df = pd.DataFrame(expr_matrix, columns=var_names, index=obs_islet_sorted.index)

    # ── 5. Write per-islet CSVs ──────────────────────────────────────────
    print(f"\n5. Writing per-islet CSVs to {output_dir}")
    os.makedirs(output_dir, exist_ok=True)

    # Combine obs + expression
    cell_data = pd.concat([
        obs_islet_sorted[["X_centroid", "Y_centroid", "phenotype", "cell_region",
                           "Cell Area", "Nucleus Area"]],
        expr_df
    ], axis=1)
    cell_data["combined_islet_id"] = obs_islet_sorted["combined_islet_id"]

    n_written = 0
    n_skipped = 0
    total_cells = 0
    for cid, group in cell_data.groupby("combined_islet_id"):
        # Use combined_islet_id as filename
        fname = f"{cid}.csv"
        fpath = os.path.join(output_dir, fname)
        out = group.drop(columns=["combined_islet_id"])
        out.to_csv(fpath, index=False)
        n_written += 1
        total_cells += len(out)

    print(f"\n6. Summary:")
    print(f"   Files written: {n_written}")
    print(f"   Total cells extracted: {total_cells:,}")
    print(f"   Avg cells per islet: {total_cells / max(n_written, 1):.0f}")
    print(f"   Output directory: {output_dir}")

    # Check total size
    total_size = sum(
        os.path.getsize(os.path.join(output_dir, f))
        for f in os.listdir(output_dir) if f.endswith(".csv")
    )
    print(f"   Total size: {total_size / (1024*1024):.1f} MB")

    print(f"\n{'='*60}")
    print("DONE")
    print(f"{'='*60}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Extract per-islet single-cell CSVs for drill-down viewer"
    )
    parser.add_argument(
        "--input", default=DEFAULT_SC_H5AD,
        help="Path to single-cell H5AD (phenotyped)"
    )
    parser.add_argument(
        "--islets", default=DEFAULT_ISLET_H5AD,
        help="Path to islets_core_fixed.h5ad (qualified islet list)"
    )
    parser.add_argument(
        "--output-dir", default=DEFAULT_OUTPUT_DIR,
        help="Output directory for per-islet CSVs"
    )
    args = parser.parse_args()

    for path, label in [
        (args.input, "single-cell H5AD"),
        (args.islets, "islets_core_fixed.h5ad"),
    ]:
        if not os.path.exists(path):
            print(f"ERROR: {label} not found: {path}")
            sys.exit(1)

    extract_cells(args.input, args.islets, args.output_dir)
