#!/usr/bin/env python3
"""
Phase 1 (typer B) — Astir probabilistic cell typing (run in the `astir-env` conda env).

Astir learns each cell type's marker distribution from the positive-marker YAML and assigns
per-cell posterior probabilities with NO per-marker gating — sidestepping the empirically
demonstrated unreliability of GMM gates on this 59-plex panel (SOX2/EpCAM/Ker8_18 over-called,
PanCK under-called). It is batch-aware via a per-donor design matrix and outputs a confidence
(max posterior) usable for cluster/QC triage.

Pipeline:
  raw MFI -> arcsinh(x/cofactor) -> Astir(expr, markers.yaml, design=one-hot(imageid))
          -> fit_type -> per-cell posterior + argmax label ('Unknown' below --threshold)

Outputs (per --out-prefix): <prefix>_cells.parquet (cell_id, imageid, phenotype_astir,
astir_prob), <prefix>_probs.parquet (full posterior), <prefix>_assignments.csv summary.

Smoke (single donor):
  python scripts/senior/run_astir.py --donors 6374 --max-cells-per-donor 200000 \
         --markers scripts/senior/astir_markers.yaml --epochs 30 --out-prefix /tmp/astir_6374
"""
from __future__ import annotations
import argparse
from pathlib import Path
import numpy as np
import pandas as pd
import anndata as ad
import yaml

BASE = Path.home() / "IO60panc2nd" / "senior_data"
INP = BASE / "singlecell_protein.h5ad"
YAML_DEFAULT = Path(__file__).resolve().parent / "astir_markers.yaml"


def load_subset(path, donors, max_cells_per_donor, donor_col):
    print(f"[load] reading {path} into memory ...")
    a = ad.read_h5ad(path)
    print(f"[load] full {a.shape}")
    if donors is None and max_cells_per_donor is None:
        return a
    dvals = a.obs[donor_col].astype(str).to_numpy()
    sel = np.ones(dvals.size, dtype=bool) if donors is None else np.isin(dvals, [str(x) for x in donors])
    idx = np.where(sel)[0]
    if max_cells_per_donor is not None:
        keep = []
        for d in pd.unique(dvals[idx]):
            di = idx[dvals[idx] == d]
            if di.size > max_cells_per_donor:
                di = di[np.linspace(0, di.size - 1, max_cells_per_donor).astype(int)]
            keep.append(di)
        idx = np.sort(np.concatenate(keep))
    a = a[idx].copy()
    print(f"[load] subset {a.shape} | donors {pd.unique(a.obs[donor_col]).tolist()}")
    return a


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--input", type=Path, default=INP)
    ap.add_argument("--markers", type=Path, default=YAML_DEFAULT)
    ap.add_argument("--donor-col", default="imageid")
    ap.add_argument("--donors", nargs="*", default=None)
    ap.add_argument("--max-cells-per-donor", type=int, default=None)
    ap.add_argument("--cofactor", type=float, default=5.0)
    ap.add_argument("--epochs", type=int, default=30)
    ap.add_argument("--batch-size", type=int, default=8192)
    ap.add_argument("--lr", type=float, default=2e-3)
    ap.add_argument("--threshold", type=float, default=0.5, help="min max-posterior to assign (else Unknown)")
    ap.add_argument("--no-design", action="store_true", help="disable per-donor batch design matrix")
    ap.add_argument("--out-prefix", type=Path, required=True)
    a = ap.parse_args(argv)

    import torch
    from astir import Astir
    print(f"[astir] torch {torch.__version__} cuda={torch.cuda.is_available()}")

    adata = load_subset(a.input, a.donors, a.max_cells_per_donor, a.donor_col)
    X = adata.X if isinstance(adata.X, np.ndarray) else adata.X.toarray()
    X = np.arcsinh(np.asarray(X, dtype=np.float64) / a.cofactor)
    expr = pd.DataFrame(X, columns=list(adata.var_names), index=adata.obs_names.astype(str))

    with open(a.markers) as f:
        marker_dict = yaml.safe_load(f)

    design = None
    if not a.no_design and adata.obs[a.donor_col].nunique() > 1:
        design = pd.get_dummies(adata.obs[a.donor_col].astype(str)).astype(float)
        design.index = expr.index
        print(f"[astir] batch design: {design.shape[1]} donors")

    ast = Astir(expr, marker_dict, design=design)
    if torch.cuda.is_available():
        try:
            ast = Astir(expr, marker_dict, design=design, dtype=torch.float32)
        except Exception:
            pass
    print(f"[astir] fit_type epochs={a.epochs} batch={a.batch_size} lr={a.lr}")
    ast.fit_type(max_epochs=a.epochs, batch_size=a.batch_size, learning_rate=a.lr, n_init=1)

    probs = ast.get_celltype_probabilities()
    maxp = probs.max(axis=1)
    label = probs.idxmax(axis=1)
    label = label.where(maxp >= a.threshold, other="Unknown")

    print("\n[astir] cell-type distribution:")
    vc = label.value_counts(); tot = len(label)
    for k, v in vc.items():
        print(f"    {k:24s} {v:7d}  {v/tot:6.1%}")
    print(f"\n[astir] mean max-posterior = {float(maxp.mean()):.3f} | "
          f"Unknown (<{a.threshold}) = {(label=='Unknown').mean():.1%}")

    a.out_prefix.parent.mkdir(parents=True, exist_ok=True)
    out = pd.DataFrame({
        "cell_id": expr.index,
        a.donor_col: adata.obs[a.donor_col].astype(str).values,
        "phenotype_astir": label.values,
        "astir_prob": maxp.values,
    })
    try:
        out.to_parquet(f"{a.out_prefix}_cells.parquet")
        probs.to_parquet(f"{a.out_prefix}_probs.parquet")
        print(f"[done] wrote {a.out_prefix}_cells.parquet")
    except Exception as e:
        out.to_csv(f"{a.out_prefix}_cells.csv", index=False)
        probs.to_csv(f"{a.out_prefix}_probs.csv")
        print(f"[done] parquet unavailable ({e}); wrote {a.out_prefix}_cells.csv")
    vc.to_csv(f"{a.out_prefix}_assignments.csv")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
