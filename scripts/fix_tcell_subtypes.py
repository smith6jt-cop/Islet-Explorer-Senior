#!/usr/bin/env python3
"""LEGACY LANE, NOT RUNNABLE HERE -- recorded 2026-09-18.

``scripts/retire/dispositions.py`` records this file as a verified orphan (code tier 1):
nothing in the tree imports it, and the inputs named below are absent from this checkout.
It is kept as the record of how that artifact was built, not as a description of anything
current.

The current PhenoCycler typing is the subpanel classification under
``core/phenocycler/scripts/typing/`` -- context from QuPath objects, then one question per
context against reference levels drawn from the ordinary tissue, fitted PER DONOR.
``data/development/typing/cells/`` holds 25,230,109 cells over 20 donors (2026-09-15 build,
eight classes plus ``Excluded``). Nothing below describes that method.

Fix CD4/CD8 T cell subtype phenotyping in canonical single-cell H5AD.

Reclassifies T cell lineage cells (CD8a Tcell, CD4 Tcell, T cell) based on
rescaled CD8a/CD4 marker expression using a 0.5 gate threshold:
  - CD8a >= 0.5 AND CD4 < 0.5  → CD8a Tcell
  - CD4 >= 0.5 AND CD8a < 0.5  → CD4 Tcell
  - Both positive or both negative → T cell (unresolved)

Usage:
    conda activate scvi-env
    python scripts/fix_tcell_subtypes.py --dry-run   # preview changes
    python scripts/fix_tcell_subtypes.py              # apply and save
"""

import argparse
import shutil
from datetime import datetime
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd


CANONICAL_H5AD = Path("single_cell_analysis/CODEX_scvi_BioCov_phenotyped_newDuctal.h5ad")
TCELL_PHENOTYPES = {"CD8a Tcell", "CD4 Tcell", "T cell"}
GATE = 0.5


def main():
    parser = argparse.ArgumentParser(description="Fix CD4/CD8 T cell subtype phenotyping")
    parser.add_argument("--dry-run", action="store_true", help="Preview changes without writing")
    parser.add_argument("--input", type=Path, default=CANONICAL_H5AD, help="Input H5AD path")
    args = parser.parse_args()

    h5ad_path = args.input
    if not h5ad_path.exists():
        raise FileNotFoundError(f"H5AD not found: {h5ad_path}")

    print(f"Loading {h5ad_path} ...")
    adata = ad.read_h5ad(h5ad_path)
    print(f"  {adata.n_obs:,} cells, {adata.n_vars} markers")

    # Verify required markers exist
    for marker in ["CD8a", "CD4"]:
        if marker not in adata.var_names:
            raise ValueError(f"Marker '{marker}' not found in var_names: {list(adata.var_names)}")

    # Identify T cell lineage cells
    phenotype_col = "phenotype"
    if phenotype_col not in adata.obs.columns:
        raise ValueError(f"Column '{phenotype_col}' not found in obs")

    tcell_mask = adata.obs[phenotype_col].isin(TCELL_PHENOTYPES)
    n_tcells = tcell_mask.sum()
    print(f"\nT cell lineage cells: {n_tcells:,}")

    # Show current distribution
    current_counts = adata.obs.loc[tcell_mask, phenotype_col].value_counts()
    print("\nCurrent phenotype distribution:")
    for pheno, count in current_counts.items():
        print(f"  {pheno}: {count:,}")

    # Get marker expression for T cell lineage
    cd8a_idx = list(adata.var_names).index("CD8a")
    cd4_idx = list(adata.var_names).index("CD4")

    X = adata.X
    if hasattr(X, "toarray"):
        cd8a_vals = np.asarray(X[tcell_mask, cd8a_idx].toarray()).flatten()
        cd4_vals = np.asarray(X[tcell_mask, cd4_idx].toarray()).flatten()
    else:
        cd8a_vals = np.asarray(X[tcell_mask, cd8a_idx]).flatten()
        cd4_vals = np.asarray(X[tcell_mask, cd4_idx]).flatten()

    # Apply gating
    cd8a_pos = cd8a_vals >= GATE
    cd4_pos = cd4_vals >= GATE

    new_phenotypes = np.where(
        cd8a_pos & ~cd4_pos, "CD8a Tcell",
        np.where(
            cd4_pos & ~cd8a_pos, "CD4 Tcell",
            "T cell"
        )
    )

    # Compare old vs new
    old_phenotypes = adata.obs.loc[tcell_mask, phenotype_col].values
    changed_mask = old_phenotypes != new_phenotypes
    n_changed = changed_mask.sum()

    print(f"\nNew phenotype distribution (gate={GATE}):")
    new_series = pd.Series(new_phenotypes)
    for pheno, count in new_series.value_counts().items():
        print(f"  {pheno}: {count:,}")

    print(f"\nCells that would change: {n_changed:,} / {n_tcells:,} ({100*n_changed/n_tcells:.2f}%)")

    if n_changed > 0:
        # Show transition matrix
        transitions = pd.crosstab(
            pd.Series(old_phenotypes[changed_mask], name="Old"),
            pd.Series(new_phenotypes[changed_mask], name="New")
        )
        print("\nTransition matrix (changed cells only):")
        print(transitions.to_string())

        # Show cell regions of changed cells
        tcell_indices = np.where(tcell_mask)[0]
        changed_indices = tcell_indices[changed_mask]
        if "cell_region" in adata.obs.columns:
            regions = adata.obs.iloc[changed_indices]["cell_region"].value_counts()
            print("\nChanged cells by region:")
            for region, count in regions.items():
                print(f"  {region}: {count:,}")

        if "Parent" in adata.obs.columns:
            parents = adata.obs.iloc[changed_indices]["Parent"].value_counts()
            n_islet = sum(c for p, c in parents.items() if "Islet" in str(p))
            n_tissue = sum(c for p, c in parents.items() if "Islet" not in str(p))
            print(f"\nChanged cells: {n_islet} in islets, {n_tissue} in tissue background")

    if args.dry_run:
        print("\n[DRY RUN] No changes written.")
        return

    if n_changed == 0:
        print("\nNo changes needed — phenotypes already correct.")
        return

    # Create backup
    backup_path = h5ad_path.with_suffix(f".backup_{datetime.now():%Y%m%d_%H%M%S}.h5ad")
    print(f"\nCreating backup: {backup_path}")
    shutil.copy2(h5ad_path, backup_path)

    # Apply changes
    adata.obs.loc[tcell_mask, phenotype_col] = new_phenotypes
    print(f"Writing updated H5AD to {h5ad_path} ...")
    adata.write_h5ad(h5ad_path)
    print("Done.")


if __name__ == "__main__":
    main()
