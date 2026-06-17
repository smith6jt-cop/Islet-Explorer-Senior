#!/usr/bin/env python3
"""
3B(+) Senior: per-islet marker positivity ("markers") from the phenotyped cells.

scimap rescaled adata.X to 0..1 with 0.5 = the gate, so (X > 0.5) == positive.
For each islet (core / peri) and each of the 59 markers, computes the fraction of
cells positive (pos_frac), and stores per-islet pos_frac (islets x markers) in the
islet h5ads' obsm['marker_pos_frac'] + a few convenience obs columns. Replaces the
groovy 'markers' sheet (no groovy). Run in the `scimap` conda env.
"""
from pathlib import Path
import numpy as np
import pandas as pd
import anndata as ad

SD = Path(__file__).resolve().parents[2] / "data"
PHENO = SD / "singlecell_protein_phenotyped.h5ad"


def main():
    a = ad.read_h5ad(PHENO)
    markers = list(a.var_names)
    parent = a.obs["Parent"].astype(str)
    is_peri = parent.str.endswith("_exp20um").to_numpy()
    is_core = (parent.str.startswith("Islet_").to_numpy()) & ~is_peri
    base = parent.str.replace("_exp20um", "", regex=False)
    islet_id = (a.obs["imageid"].astype(str) + "_" + base).to_numpy()
    X = a.X if isinstance(a.X, np.ndarray) else a.X.toarray()

    for region, mask, f in [("core", is_core, "islets_core_senior.h5ad"),
                            ("peri", is_peri, "islets_peri_senior.h5ad")]:
        idx = np.where(mask)[0]
        df = pd.DataFrame((X[idx] > 0.5), columns=markers)
        df["islet_id"] = islet_id[idx]
        posfrac = df.groupby("islet_id")[markers].mean()
        ai = ad.read_h5ad(SD / f)
        pf = posfrac.reindex(ai.obs_names)
        ai.obsm["marker_pos_frac"] = pf.to_numpy(dtype=np.float32)
        ai.uns["marker_pos_frac_names"] = markers
        for m in ["INS", "GCG", "SST", "CD3e", "CD68", "CD8"]:
            if m in markers:
                ai.obs[f"posfrac_{m}"] = pf[m].to_numpy()
        ai.write_h5ad(SD / f)
        print(f"[{region}] marker_pos_frac {ai.obsm['marker_pos_frac'].shape}; "
              f"median INS+ frac across islets = {np.nanmedian(pf['INS']):.3f}")
    # disease sanity on core
    ai = ad.read_h5ad(SD / "islets_core_senior.h5ad")
    print("[core] median INS+ frac by donor status:")
    print(ai.obs.groupby("donor_status", observed=True)["posfrac_INS"].median().round(3).to_string())


if __name__ == "__main__":
    main()
