#!/usr/bin/env python3
"""
3D Senior: per-islet single-cell CSVs for the drill-down overlay, mirroring
Junior's data/cells/<imageid>_<islet_key>.csv. Each file holds the core + peri
cells of one islet with columns the drill-down reads: X_centroid, Y_centroid,
phenotype, cell_region, + 59 raw-MFI markers (for marker coloring).

Streams the phenotyped single-cell h5ad per contiguous donor block; reads raw
MFI from /raw/X. Writes into the app's data dir (app_data/cells).

Run in the `scimap` conda env (~a few minutes; ~9.1k small CSVs).
"""
import itertools
from pathlib import Path

import h5py
import numpy as np
import pandas as pd

SD = Path.home() / "IO60panc2nd" / "senior_data"
PHENO = SD / "singlecell_protein_phenotyped.h5ad"
CELLS_OUT = SD / "app_data" / "cells"


def cat_slice(grp, s, e):
    cats = grp["categories"][()]
    cats = np.array([c.decode() if isinstance(c, (bytes, bytearray)) else c
                     for c in cats], dtype=object)
    codes = grp["codes"][s:e]
    out = np.empty(codes.shape, dtype=object)
    valid = codes >= 0
    out[valid] = cats[codes[valid]]
    out[~valid] = ""
    return out


def main():
    CELLS_OUT.mkdir(parents=True, exist_ok=True)
    f = h5py.File(PHENO, "r")
    obs = f["obs"]
    cats = [c.decode() if isinstance(c, bytes) else c
            for c in obs["imageid"]["categories"][()]]
    codes = obs["imageid"]["codes"][()]
    rawX = f["raw"]["X"]
    rvar = f["raw"]["var"]
    ridx = rvar.attrs.get("_index", "_index")
    if isinstance(ridx, bytes):
        ridx = ridx.decode()
    markers = [v.decode() if isinstance(v, bytes) else v for v in rvar[ridx][()]]

    ranges = []
    start = 0
    for code, grp in itertools.groupby(codes):
        n = sum(1 for _ in grp)
        ranges.append((cats[code], start, start + n))
        start += n

    total = 0
    for donor, s, e in ranges:
        parent = cat_slice(obs["Parent"], s, e)
        region = cat_slice(obs["cell_region"], s, e)
        pheno = cat_slice(obs["phenotype"], s, e)
        islet = np.array([p[:-8] if p.endswith("_exp20um")
                          else ("" if (p == "" or p.startswith("Annotation")) else p)
                          for p in parent], dtype=object)
        keep = islet != ""                      # core + peri cells (drop tissue)
        if keep.sum() == 0:
            continue
        block = rawX[s:e, :]                     # raw MFI, one contiguous read
        df = pd.DataFrame(block[keep], columns=markers).round(2)
        df.insert(0, "X_centroid", obs["X_centroid"][s:e][keep].astype("float32"))
        df.insert(1, "Y_centroid", obs["Y_centroid"][s:e][keep].astype("float32"))
        df.insert(2, "phenotype", pheno[keep])
        df.insert(3, "cell_region", region[keep])
        df["__islet"] = islet[keep]
        n_islets = 0
        for islet_key, grp_df in df.groupby("__islet"):
            grp_df.drop(columns="__islet").to_csv(
                CELLS_OUT / f"{donor}_{islet_key}.csv", index=False)
            n_islets += 1
        total += n_islets
        print(f"[{donor}] {int(keep.sum()):>8,} cells -> {n_islets} islet CSVs")
    f.close()
    print(f"[done] {total} per-islet CSVs in {CELLS_OUT}")


if __name__ == "__main__":
    main()
