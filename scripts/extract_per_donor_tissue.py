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

Extract per-donor tissue-wide cell data for the Spatial tab tissue scatter.

For each of the 15 donors, extracts ALL cells (core, peri-islet, and tissue
background) with spatial coordinates, phenotype, region label, and islet name.

Output: data/donors/{imageid}.csv (~15 files, ~225 MB total)

Each CSV has columns:
  X_centroid, Y_centroid, phenotype, cell_region, islet_name

Where cell_region is one of:
  - "core"   : cells inside islet core (Parent = "Islet_N")
  - "peri"   : cells in 20µm expansion zone (Parent = "Islet_N_exp20um")
  - "tissue" : all other cells (tissue background, annotations, etc.)

And islet_name is:
  - "Islet_N" for core/peri cells
  - ""        for tissue background cells

Usage:
    python scripts/extract_per_donor_tissue.py
    python scripts/extract_per_donor_tissue.py --input path/to/single_cell.h5ad --output-dir data/donors
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
DEFAULT_OUTPUT_DIR = os.path.join(PROJECT_ROOT, "data", "donors")


def parse_parent(parent_str):
    """Parse Parent column to extract islet name and region.

    Returns (islet_name, cell_region):
      - ("Islet_N", "core")   for Parent = "Islet_N"
      - ("Islet_N", "peri")   for Parent = "Islet_N_exp20um"
      - ("", "tissue")        for everything else
    """
    if not isinstance(parent_str, str):
        return "", "tissue"
    m = re.match(r"^(Islet_\d+)(_exp20um)?$", parent_str)
    if m:
        islet_name = m.group(1)
        region = "peri" if m.group(2) else "core"
        return islet_name, region
    return "", "tissue"


def extract_donor_tissue(sc_path, output_dir):
    """Extract per-donor tissue CSVs from single-cell H5AD."""
    import anndata as ad

    print("=" * 60)
    print("Extracting Per-Donor Tissue Data")
    print("=" * 60)

    # ── 1. Load single-cell H5AD (backed mode for memory efficiency) ──
    print(f"\n1. Loading single-cell H5AD: {sc_path}")
    sc = ad.read_h5ad(sc_path, backed="r")
    print(f"   Total cells: {sc.shape[0]:,}")

    # ── 2. Extract obs metadata ──
    print("\n2. Extracting cell metadata...")
    obs = sc.obs[["imageid", "Parent", "phenotype",
                  "X_centroid", "Y_centroid"]].copy()
    obs["imageid"] = obs["imageid"].astype(str)
    obs["Parent"] = obs["Parent"].astype(str)
    obs["phenotype"] = obs["phenotype"].astype(str)

    # Parse Parent column to get cell_region and islet_name
    parsed = obs["Parent"].apply(parse_parent)
    obs["islet_name"] = [p[0] for p in parsed]
    obs["cell_region"] = [p[1] for p in parsed]

    # ── 3. Write per-donor CSVs ──
    print(f"\n3. Writing per-donor CSVs to {output_dir}")
    os.makedirs(output_dir, exist_ok=True)

    donors = sorted(obs["imageid"].unique())
    print(f"   Donors: {len(donors)}")

    total_cells = 0
    total_size = 0
    for donor in donors:
        mask = obs["imageid"] == donor
        donor_cells = obs.loc[mask, ["X_centroid", "Y_centroid", "phenotype",
                                      "cell_region", "islet_name"]].copy()
        fpath = os.path.join(output_dir, f"{donor}.csv")
        donor_cells.to_csv(fpath, index=False)
        fsize = os.path.getsize(fpath)
        total_cells += len(donor_cells)
        total_size += fsize

        # Region breakdown
        region_counts = donor_cells["cell_region"].value_counts().to_dict()
        print(f"   {donor}: {len(donor_cells):,} cells "
              f"(core={region_counts.get('core', 0):,}, "
              f"peri={region_counts.get('peri', 0):,}, "
              f"tissue={region_counts.get('tissue', 0):,}) "
              f"-> {fsize / (1024*1024):.1f} MB")

    print(f"\n4. Summary:")
    print(f"   Files written: {len(donors)}")
    print(f"   Total cells: {total_cells:,}")
    print(f"   Total size: {total_size / (1024*1024):.1f} MB")
    print(f"   Output directory: {output_dir}")

    print(f"\n{'='*60}")
    print("DONE")
    print(f"{'='*60}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Extract per-donor tissue CSVs for Spatial tab"
    )
    parser.add_argument(
        "--input", default=DEFAULT_SC_H5AD,
        help="Path to single-cell H5AD (phenotyped)"
    )
    parser.add_argument(
        "--output-dir", default=DEFAULT_OUTPUT_DIR,
        help="Output directory for per-donor CSVs"
    )
    args = parser.parse_args()

    if not os.path.exists(args.input):
        print(f"ERROR: single-cell H5AD not found: {args.input}")
        sys.exit(1)

    extract_donor_tissue(args.input, args.output_dir)
