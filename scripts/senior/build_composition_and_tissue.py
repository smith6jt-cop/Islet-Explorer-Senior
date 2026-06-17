#!/usr/bin/env python3
"""
3D Senior: stream the phenotyped single-cell h5ad ONCE (per contiguous donor
block) to produce two app inputs at once:

  1. data/islet_composition.parquet
       per-core-islet hormone-coexpression counts (cells_total, Ins/Glu/Stt
       _single/_any, Multi_Pos, Triple_Neg) — the 'composition' groovy sheet
       analog. Hormone positivity = scimap-rescaled X > 0.5 (the gate), same
       convention as compute_islet_markers.py.

  2. data/parquet/tissue/case_id=<donor>/part-0.parquet
       per-donor tissue cells for the Spatial tab scatter, with EXACTLY the
       columns load_donor_tissue() selects: X_centroid, Y_centroid, phenotype,
       cell_region, islet_name.

INS/GCG/SST live at var indices 13/46/35. Each donor is a single contiguous
row block (verified), so X[s:e] is one hyperslab read (~0.6 GB peak).

Run in the `scimap` conda env (anndata 0.11 / h5py).
"""
import itertools
from pathlib import Path

import h5py
import numpy as np
import pandas as pd
import pyarrow as pa
import pyarrow.parquet as pq

SD = Path(__file__).resolve().parents[2] / "data"
PHENO = SD / "singlecell_protein_phenotyped.h5ad"
TISSUE_ROOT = SD / "parquet" / "tissue"


def cat_slice(grp, s, e):
    """Materialize a categorical obs column slice [s:e] as a str object array."""
    cats = grp["categories"][()]
    cats = np.array([c.decode() if isinstance(c, (bytes, bytearray)) else c
                     for c in cats], dtype=object)
    codes = grp["codes"][s:e]
    out = np.empty(codes.shape, dtype=object)
    valid = codes >= 0
    out[valid] = cats[codes[valid]]
    out[~valid] = ""
    return out


def var_index_names(f):
    attr = f["var"].attrs.get("_index", "_index")
    if isinstance(attr, bytes):
        attr = attr.decode()
    return [v.decode() if isinstance(v, bytes) else v for v in f["var"][attr][()]]


def main():
    f = h5py.File(PHENO, "r")
    obs = f["obs"]
    cats = [c.decode() if isinstance(c, bytes) else c
            for c in obs["imageid"]["categories"][()]]
    codes = obs["imageid"]["codes"][()]
    vn = var_index_names(f)
    iINS, iGCG, iSST = vn.index("INS"), vn.index("GCG"), vn.index("SST")
    X = f["X"]

    # Contiguous per-donor ranges (each donor is one run of imageid codes).
    ranges = []
    start = 0
    for code, grp in itertools.groupby(codes):
        n = sum(1 for _ in grp)
        ranges.append((cats[code], start, start + n))
        start += n

    TISSUE_ROOT.mkdir(parents=True, exist_ok=True)
    comp_rows = []

    for donor, s, e in ranges:
        parent = cat_slice(obs["Parent"], s, e)
        region = cat_slice(obs["cell_region"], s, e)
        pheno = cat_slice(obs["phenotype"], s, e)
        xc = obs["X_centroid"][s:e].astype("float32")
        yc = obs["Y_centroid"][s:e].astype("float32")

        # islet a cell belongs to: strip the _exp20um peri suffix; tissue -> "".
        islet_name = np.array([
            p[:-8] if p.endswith("_exp20um")
            else ("" if (p == "" or p.startswith("Annotation")) else p)
            for p in parent
        ], dtype=object)

        # (2) tissue partition — exactly the 5 columns the Spatial loader reads.
        tdf = pd.DataFrame({
            "X_centroid": xc, "Y_centroid": yc,
            "phenotype": pheno, "cell_region": region, "islet_name": islet_name,
        })
        ddir = TISSUE_ROOT / f"case_id={donor}"
        ddir.mkdir(parents=True, exist_ok=True)
        pq.write_table(pa.Table.from_pandas(tdf, preserve_index=False),
                       ddir / "part-0.parquet")

        # (1) composition — core cells only (one row per islet core).
        core = (region == "core") & (islet_name != "")
        n_islets = 0
        if core.sum() > 0:
            block = X[s:e, :]                       # one contiguous read
            ins = block[:, iINS] > 0.5
            glu = block[:, iGCG] > 0.5
            sst = block[:, iSST] > 0.5
            cdf = pd.DataFrame({
                "name": islet_name[core],
                "ins": ins[core], "glu": glu[core], "sst": sst[core],
            })
            npos = cdf[["ins", "glu", "sst"]].sum(axis=1)
            cdf["multi"] = npos >= 2
            cdf["tneg"] = npos == 0
            cdf["ins_s"] = cdf["ins"] & ~cdf["glu"] & ~cdf["sst"]
            cdf["glu_s"] = cdf["glu"] & ~cdf["ins"] & ~cdf["sst"]
            cdf["stt_s"] = cdf["sst"] & ~cdf["ins"] & ~cdf["glu"]
            agg = cdf.groupby("name").agg(
                cells_total=("ins", "size"),
                Ins_any=("ins", "sum"), Glu_any=("glu", "sum"), Stt_any=("sst", "sum"),
                Multi_Pos=("multi", "sum"), Triple_Neg=("tneg", "sum"),
                Ins_single=("ins_s", "sum"), Glu_single=("glu_s", "sum"),
                Stt_single=("stt_s", "sum"),
            ).reset_index()
            agg.insert(0, "Case ID", donor)
            comp_rows.append(agg)
            n_islets = agg["name"].nunique()
        print(f"[{donor}] tissue {len(tdf):>9,} cells | core islets {n_islets}")

    comp = pd.concat(comp_rows, ignore_index=True)
    comp["islet_id"] = comp["name"].str.replace("Islet_", "", regex=False)
    # integer count columns
    int_cols = ["cells_total", "Ins_any", "Glu_any", "Stt_any", "Multi_Pos",
                "Triple_Neg", "Ins_single", "Glu_single", "Stt_single"]
    comp[int_cols] = comp[int_cols].astype("int64")
    comp.to_parquet(SD / "islet_composition.parquet")
    print(f"\n[composition] {comp.shape[0]} islets x {comp.shape[1]} cols "
          f"-> {SD/'islet_composition.parquet'}")
    print("[composition] cells_total: median",
          int(comp['cells_total'].median()), "total core cells",
          int(comp['cells_total'].sum()))
    f.close()


if __name__ == "__main__":
    main()
