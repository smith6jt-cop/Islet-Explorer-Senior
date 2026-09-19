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

Compute per-islet peri-islet neighborhood metrics from single-cell H5AD.

Reads the phenotyped single-cell H5AD (2.6M cells) and computes metrics
quantifying the cellular microenvironment around each islet (20µm expansion zone).

Metrics computed (4 categories):
  1. Peri-islet composition — proportion & count of each phenotype in _exp20um zone
  2. Immune infiltration — immune fractions (peri & core), peri/core ratio, ratios
  3. Enrichment z-scores — Poisson z comparing peri-islet vs tissue-wide proportion
  4. Distance metrics — min distance from islet centroid to nearest immune cells

Output: data/neighborhood_metrics.csv (1,015 rows × ~61 columns)

Usage:
    python scripts/compute_neighborhood_metrics.py
    python scripts/compute_neighborhood_metrics.py --input path/to/single_cell.h5ad --output data/neighborhood_metrics.csv
"""

import argparse
import os
import re
import sys
import warnings

import numpy as np
import pandas as pd

warnings.filterwarnings("ignore", category=FutureWarning)

# ── Constants ──────────────────────────────────────────────────────────────

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
PROJECT_ROOT = os.path.join(SCRIPT_DIR, "..")

DEFAULT_SC_H5AD = os.path.join(
    PROJECT_ROOT, "single_cell_analysis",
    "CODEX_scvi_BioCov_phenotyped_newDuctal.h5ad"
)
DEFAULT_ISLET_H5AD = os.path.join(
    PROJECT_ROOT, "islet_analysis", "islets_core_fixed.h5ad"
)
DEFAULT_OUTPUT = os.path.join(PROJECT_ROOT, "data", "neighborhood_metrics.csv")

# Immune phenotypes — drives the AGGREGATE metrics where "what counts as
# immune" must be defined narrowly: immune_frac_peri/core, immune_ratio,
# tcell_density_peri, immune_count_peri/core, and the min_dist_immune_mean
# summary. Do NOT extend this list — it would make those aggregates
# meaningless (e.g. "fraction of peri cells that are 'Beta cell or B cell'").
IMMUNE_TYPES = [
    "CD8a Tcell", "CD4 Tcell", "T cell", "B cell",
    "Macrophage", "APCs", "Immune"
]

# Per-type metrics (enrich_z_*, min_dist_*) are computed for ALL phenotypes
# discovered in the data, minus this noise bucket. This replaces the earlier
# 7-phenotype ENRICH_TYPES subset so the app can expose phenotype-driven
# selectors uniformly across every cell type, not a hand-curated subset.
PER_TYPE_EXCLUDE = {"Unknown"}

PIXEL_SIZE_UM = 0.3774  # micrometers per pixel, matches app constant


def parse_parent(parent_str):
    """Parse Parent column to extract region type and islet ID.

    Returns (islet_name, region) where region is 'core' or 'peri'.
    Returns (None, None) for non-islet cells.
    """
    if not isinstance(parent_str, str):
        return None, None
    m = re.match(r"^(Islet_\d+)(_exp20um)?$", parent_str)
    if m:
        islet_name = m.group(1)
        region = "peri" if m.group(2) else "core"
        return islet_name, region
    return None, None


def compute_metrics(sc_path, islet_path, output_path):
    """Main computation: read single-cell data, compute per-islet neighborhood metrics."""
    import anndata as ad

    print("=" * 60)
    print("Computing Peri-Islet Neighborhood Metrics")
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
    print(f"   Qualified islets: {len(qualified)}")

    # ── 2. Load single-cell H5AD (backed to manage memory) ───────────────
    print(f"\n2. Loading single-cell H5AD: {sc_path}")
    sc = ad.read_h5ad(sc_path, backed="r")
    print(f"   Total cells: {sc.shape[0]:,}")
    print(f"   Markers: {sc.shape[1]}")

    # ── 3. Parse Parent column → islet assignment + region ───────────────
    print("\n3. Parsing cell-islet assignments from Parent column...")
    obs = sc.obs[["imageid", "Parent", "phenotype", "X_centroid", "Y_centroid"]].copy()
    obs["imageid"] = obs["imageid"].astype(str)
    obs["Parent"] = obs["Parent"].astype(str)
    obs["phenotype"] = obs["phenotype"].astype(str)

    parsed = obs["Parent"].apply(parse_parent)
    obs["islet_name"] = [p[0] for p in parsed]
    obs["cell_region"] = [p[1] for p in parsed]

    # Filter to islet-associated cells only
    islet_cells = obs[obs["islet_name"].notna()].copy()
    islet_cells["combined_islet_id"] = (
        islet_cells["imageid"] + "_" + islet_cells["islet_name"]
    )
    print(f"   Islet-associated cells: {len(islet_cells):,}")
    print(f"     Core: {(islet_cells['cell_region'] == 'core').sum():,}")
    print(f"     Peri: {(islet_cells['cell_region'] == 'peri').sum():,}")

    # Filter to qualified islets only
    qualified_set = set(qualified["combined_islet_id"])
    islet_cells = islet_cells[
        islet_cells["combined_islet_id"].isin(qualified_set)
    ].copy()
    print(f"   Cells in qualified islets: {len(islet_cells):,}")

    # ── 4. Compute tissue-wide phenotype proportions (baseline) ──────────
    print("\n4. Computing tissue-wide phenotype baseline...")
    all_phenotypes = sorted(obs["phenotype"].unique())
    tissue_total = len(obs)
    tissue_counts = obs["phenotype"].value_counts()
    tissue_props = (tissue_counts / tissue_total).to_dict()
    print(f"   {len(all_phenotypes)} phenotypes, {tissue_total:,} total cells")
    # Per-type metrics span all phenotypes minus PER_TYPE_EXCLUDE.
    per_type_phenotypes = [p for p in all_phenotypes if p not in PER_TYPE_EXCLUDE]
    print(f"   Per-type metrics (enrich_z, min_dist) computed for "
          f"{len(per_type_phenotypes)} phenotypes")

    # ── 5. Compute per-islet metrics ─────────────────────────────────────
    print("\n5. Computing per-islet neighborhood metrics...")
    results = []

    peri_cells = islet_cells[islet_cells["cell_region"] == "peri"]
    core_cells = islet_cells[islet_cells["cell_region"] == "core"]

    peri_grouped = peri_cells.groupby("combined_islet_id")
    core_grouped = core_cells.groupby("combined_islet_id")

    for _, row in qualified.iterrows():
        cid = row["combined_islet_id"]
        rec = {"combined_islet_id": cid}

        # --- Peri-islet composition ---
        if cid in peri_grouped.groups:
            peri = peri_grouped.get_group(cid)
            peri_total = len(peri)
            peri_pheno_counts = peri["phenotype"].value_counts()
        else:
            peri = pd.DataFrame()
            peri_total = 0
            peri_pheno_counts = pd.Series(dtype=int)

        rec["total_cells_peri"] = peri_total

        for pheno in all_phenotypes:
            safe_name = pheno.replace(" ", "_").replace("+", "plus")
            cnt = int(peri_pheno_counts.get(pheno, 0))
            rec[f"peri_count_{safe_name}"] = cnt
            rec[f"peri_prop_{safe_name}"] = cnt / peri_total if peri_total > 0 else np.nan

        # --- Core composition (for ratio metrics) ---
        if cid in core_grouped.groups:
            core = core_grouped.get_group(cid)
            core_total = len(core)
            core_pheno_counts = core["phenotype"].value_counts()
        else:
            core = pd.DataFrame()
            core_total = 0
            core_pheno_counts = pd.Series(dtype=int)

        rec["total_cells_core"] = core_total

        # --- Immune infiltration metrics ---
        peri_immune = sum(
            int(peri_pheno_counts.get(t, 0)) for t in IMMUNE_TYPES
        )
        core_immune = sum(
            int(core_pheno_counts.get(t, 0)) for t in IMMUNE_TYPES
        )

        rec["immune_count_peri"] = peri_immune
        rec["immune_count_core"] = core_immune
        rec["immune_frac_peri"] = (
            peri_immune / peri_total if peri_total > 0 else np.nan
        )
        rec["immune_frac_core"] = (
            core_immune / core_total if core_total > 0 else np.nan
        )
        # Peri/core immune ratio
        if core_immune > 0 and peri_total > 0:
            rec["immune_ratio"] = (peri_immune / peri_total) / (
                core_immune / core_total
            )
        else:
            rec["immune_ratio"] = np.nan

        # CD8/macrophage ratio in peri-islet
        cd8_peri = int(peri_pheno_counts.get("CD8a Tcell", 0))
        macro_peri = int(peri_pheno_counts.get("Macrophage", 0))
        rec["cd8_to_macro_ratio"] = (
            cd8_peri / macro_peri if macro_peri > 0 else np.nan
        )

        # T-cell density in peri-islet (T cells per 100 peri cells)
        tcell_peri = sum(
            int(peri_pheno_counts.get(t, 0))
            for t in ["CD8a Tcell", "CD4 Tcell", "T cell"]
        )
        rec["tcell_density_peri"] = (
            100 * tcell_peri / peri_total if peri_total > 0 else np.nan
        )

        # --- Enrichment z-scores (Poisson model) — per-phenotype ---
        for etype in per_type_phenotypes:
            safe_name = etype.replace(" ", "_").replace("+", "plus")
            observed = int(peri_pheno_counts.get(etype, 0))
            expected = tissue_props.get(etype, 0) * peri_total
            if expected > 0 and peri_total > 0:
                # Poisson z-score: (observed - expected) / sqrt(expected)
                rec[f"enrich_z_{safe_name}"] = (observed - expected) / np.sqrt(
                    expected
                )
            else:
                rec[f"enrich_z_{safe_name}"] = np.nan

        # --- Distance metrics (min distance from islet core centroid to
        #     nearest cell of each phenotype in the peri-islet zone) ---
        if peri_total > 0 and len(core) > 0:
            # Use core centroid as reference point
            core_cx = core["X_centroid"].astype(float).mean()
            core_cy = core["Y_centroid"].astype(float).mean()

            # Aggregate "any immune cell" min distance — uses IMMUNE_TYPES.
            peri_immune_mask = peri["phenotype"].isin(IMMUNE_TYPES)
            peri_immune_cells = peri[peri_immune_mask]
            if len(peri_immune_cells) > 0:
                dx = peri_immune_cells["X_centroid"].astype(float).values - core_cx
                dy = peri_immune_cells["Y_centroid"].astype(float).values - core_cy
                rec["min_dist_immune_mean"] = float(np.min(np.sqrt(dx**2 + dy**2)))
            else:
                rec["min_dist_immune_mean"] = np.nan

            # Per-phenotype min distances — extends across all phenotypes.
            peri_phenos = peri["phenotype"].values
            peri_x = peri["X_centroid"].astype(float).values
            peri_y = peri["Y_centroid"].astype(float).values
            for dtype in per_type_phenotypes:
                safe_name = dtype.replace(" ", "_").replace("+", "plus")
                type_mask = peri_phenos == dtype
                if type_mask.any():
                    tdx = peri_x[type_mask] - core_cx
                    tdy = peri_y[type_mask] - core_cy
                    rec[f"min_dist_{safe_name}"] = float(
                        np.min(np.sqrt(tdx**2 + tdy**2))
                    )
                else:
                    rec[f"min_dist_{safe_name}"] = np.nan
        else:
            rec["min_dist_immune_mean"] = np.nan
            for dtype in per_type_phenotypes:
                safe_name = dtype.replace(" ", "_").replace("+", "plus")
                rec[f"min_dist_{safe_name}"] = np.nan

        results.append(rec)

    df = pd.DataFrame(results)

    # ── 6. Summary and output ────────────────────────────────────────────
    print(f"\n6. Results summary:")
    print(f"   Total islets: {len(df)}")
    has_peri = df["total_cells_peri"] > 0
    print(f"   Islets with peri-islet data: {has_peri.sum()} ({100*has_peri.mean():.1f}%)")
    print(f"   Islets without peri data: {(~has_peri).sum()}")
    print(f"   Columns: {len(df.columns)}")

    # Quick sanity check
    if "immune_frac_peri" in df.columns:
        ifp = df.loc[has_peri, "immune_frac_peri"]
        print(f"   immune_frac_peri: mean={ifp.mean():.4f}, median={ifp.median():.4f}")

    # Write output
    os.makedirs(os.path.dirname(output_path), exist_ok=True)
    df.to_csv(output_path, index=False)
    print(f"\n   Written to: {output_path}")
    print(f"   Shape: {df.shape}")
    print(f"\n{'='*60}")
    print("DONE")
    print(f"{'='*60}")

    return df


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Compute per-islet peri-islet neighborhood metrics"
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
        "--output", default=DEFAULT_OUTPUT,
        help="Output path for neighborhood_metrics.csv"
    )
    args = parser.parse_args()

    for path, label in [
        (args.input, "single-cell H5AD"),
        (args.islets, "islets_core_fixed.h5ad"),
    ]:
        if not os.path.exists(path):
            print(f"ERROR: {label} not found: {path}")
            sys.exit(1)

    compute_metrics(args.input, args.islets, args.output)
