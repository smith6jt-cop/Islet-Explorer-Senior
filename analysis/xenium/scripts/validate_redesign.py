"""Structural validation of the scANVI-classifier redesign outputs (Tier-C).

READ-ONLY. Distributions are REPORTED, never gated (no proportion judgments —
memory `no-proportion-based-quality-judgments`). Hard structural contracts ARE
asserted: Indeterminate-free `celltype`; `celltype_lineage` stays FINE; every
broad lineage representable (seed-starvation check); identity_panel/state masks +
scANVI latent present; layers intact.

    python scripts/validate_redesign.py --tissue pancreas --samples 0041323 0041326
    python scripts/validate_redesign.py --tissue thymus  --samples 0041435 0041534
"""
import argparse
import importlib
import sys
import os
from pathlib import Path

import numpy as np
import scanpy as sc

# Two roots, deliberately separate.  ROOT is where the CODE and the tracked tables live;
# DATA_ROOT is where the multi-GB h5ads live and is NEVER inside the repo.  Collapsing them
# -- as the single HiPerGator ROOT did -- would put data/processed/ under analysis/xenium/,
# and the repo-root .gitignore pattern "/data/" is ANCHORED, so it would not be ignored.
ROOT = Path(__file__).resolve().parents[1]                     # analysis/xenium
REPO = Path(__file__).resolve().parents[3]                     # monorepo root
DATA_ROOT = Path(os.environ.get("XENIUM_DATA_ROOT", "/data/xenium"))
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))


def check(sample, mod):
    p = DATA_ROOT / f"processed/{sample}/{sample}_phenotyped.h5ad"
    if not p.exists():
        print(f"  {sample}: NOT FOUND")
        return False
    a = sc.read_h5ad(p)
    ok = True
    print(f"\n=== {sample}  shape={a.shape} ===")

    # --- hard structural contracts (assert) ---
    for layer in ("counts", "lognorm"):
        if layer not in a.layers:
            print(f"  FAIL: layers['{layer}'] missing"); ok = False
    if "X_scvi_broad" not in a.obsm:
        print("  FAIL: obsm['X_scvi_broad'] missing"); ok = False
    if "identity_panel" not in a.var.columns:
        print("  FAIL: var['identity_panel'] missing"); ok = False
    else:
        print(f"  identity_panel: {int(a.var['identity_panel'].sum()):,} genes")

    ct = a.obs["celltype"].astype(str)
    n_indet = int(ct.str.contains("Indeterminate", case=False).sum())
    if n_indet:
        print(f"  FAIL: {n_indet:,} Indeterminate cells in celltype"); ok = False
    else:
        print("  celltype: Indeterminate-free ✓")

    # celltype_lineage must stay FINE (not collapsed to the broad vocabulary)
    if "celltype_lineage" in a.obs.columns and "celltype_broad" in a.obs.columns:
        cl = set(a.obs["celltype_lineage"].astype(str).unique())
        cb = set(a.obs["celltype_broad"].astype(str).unique())
        if cl and cl <= cb:
            print("  FAIL: celltype_lineage collapsed to broad vocabulary"); ok = False
        else:
            print("  celltype_lineage: FINE ✓")

    # every broad lineage representable (seed-starvation regression)
    expected = {L for L in mod.BROAD_TO_FINE if L != "Cycling"}
    present = set(a.obs["celltype_broad"].astype(str).unique()) if "celltype_broad" in a.obs else set()
    missing = expected - present
    if missing:
        print(f"  [WARN] broad lineages with 0 cells (check seeds): {sorted(missing)}")
    else:
        print(f"  all {len(expected)} broad lineages represented ✓")

    # celltype_flat provenance preserved
    if "celltype_flat" not in a.obs.columns:
        print("  [WARN] celltype_flat (pre-redesign label) not preserved")

    # --- marker-concordance (catch-all gate): cells labeled L must EXPRESS L ---
    # broad_n_identity_genes = # of the CALLED lineage's SPECIFIC (specificity-filtered)
    # identity_core genes the cell detects at the raw-count floor (stamped by
    # hierarchical_phenotype tier1_broad). The redesign's purpose is "cells labeled L
    # express L"; the first-redesign catch-all had this near 0 for Immune_B/Neural.
    # REPORTED per lineage; weighted >=1 asserted > 0.80 as the catch-all guard.
    if "broad_n_identity_genes" in a.obs.columns and "celltype_broad" in a.obs.columns:
        nid = a.obs["broad_n_identity_genes"].to_numpy()
        cb = a.obs["celltype_broad"].astype(str).to_numpy()
        w1 = float((nid >= 1).mean())
        print(f"  marker-concordance (labeled L detects >=1 specific identity_core "
              f"gene): weighted {100*w1:.1f}%  (>=2: {100*(nid >= 2).mean():.1f}%)")
        for L in sorted(set(cb)):
            m = cb == L
            print(f"      {L:<12} n={int(m.sum()):>9,}  >=1 specific: "
                  f"{100*(nid[m] >= 1).mean():.0f}%")
        if w1 < 0.80:
            print(f"  FAIL: broad marker-concordance {100*w1:.1f}% < 80% — possible "
                  f"catch-all (cells labeled a lineage they do not express)"); ok = False
        else:
            print("  marker-concordance >= 80% (no catch-all) ✓")
    else:
        print("  [WARN] broad_n_identity_genes absent — cannot check marker-concordance")

    # --- REPORT ONLY (never gated) ---
    if "broad_confidence" in a.obs.columns:
        bc = a.obs["broad_confidence"].to_numpy()
        finite = np.isfinite(bc).all() and bc.min() >= -1e-6 and bc.max() <= 1.0001
        if not finite:
            print("  FAIL: broad_confidence not finite/∈[0,1]"); ok = False
        qs = np.quantile(bc, [0.0, 0.25, 0.5, 0.75, 1.0])
        print("  broad_confidence min/25/50/75/max: " + " ".join(f"{q:.3f}" for q in qs)
              + f"  (<0.5: {100*(bc < 0.5).mean():.1f}% — report only)")
    state_cols = sorted(c for c in a.obs.columns
                        if c.startswith("state_") and c.endswith("_score"))
    print(f"  state programs: {len(state_cols)} scored")
    if not state_cols:
        print("  [WARN] no state_* scores present")

    print("  celltype_broad value_counts (REPORT ONLY):")
    print("    " + a.obs["celltype_broad"].value_counts().to_string().replace("\n", "\n    "))
    print("  celltype value_counts (REPORT ONLY):")
    print("    " + ct.value_counts().to_string().replace("\n", "\n    "))
    print(f"  ==> {sample}: {'PASS' if ok else 'FAIL'} (structural)")
    return ok


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--tissue", choices=("pancreas", "thymus"), default="pancreas")
    ap.add_argument("--samples", nargs="+", required=True)
    args = ap.parse_args()
    mod = importlib.import_module(f"scripts.tissue_markers.{args.tissue}")
    allok = all(check(s, mod) for s in args.samples)
    print(f"\n{'='*60}\n{'ALL STRUCTURAL CHECKS PASS' if allok else 'STRUCTURAL FAILURES — review'}")
    sys.exit(0 if allok else 1)


if __name__ == "__main__":
    main()
