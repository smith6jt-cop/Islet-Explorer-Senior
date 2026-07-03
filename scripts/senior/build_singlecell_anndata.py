#!/usr/bin/env python3
"""
Build the single-cell AnnData for scimap phenotyping (run in the `scimap` conda env).

Reads the donor-partitioned cells parquet (~23.3M cells: 59 whole-cell marker
means + morphology/distance + region/islet) and assembles an AnnData mirroring
Junior's scimap input:
  X                   = raw whole-cell marker means (float32); var_names = markers
  obs                 = imageid, Parent, X/Y centroid, cell_region, islet_num,
                        morphology, signed distances, Object ID, Classification,
                        + donor metadata (Donor Status, Age, Sex, AAb)
  uns['all_markers']  = marker list

Output: data/singlecell_protein.h5ad
"""
from __future__ import annotations
import argparse
from pathlib import Path
import numpy as np
import pandas as pd
import anndata as ad
import pyarrow as pa
import pyarrow.compute as pc
import pyarrow.dataset as pads

CELLS = Path(__file__).resolve().parents[2] / "data" / "cells"
DONOR_META = Path.home() / "IO60panc2nd" / "donor_metadata_panc.xlsx"
OUT = Path(__file__).resolve().parents[2] / "data" / "singlecell_protein.h5ad"

# Non-marker columns in the parquet (everything else is a marker mean).
META_COLS = ["donor_id", "image", "object_id", "X_centroid", "Y_centroid",
             "parent_raw", "classification", "cell_region", "islet_num",
             "cell_area", "nucleus_area", "cell_length", "cell_circularity",
             "cell_solidity", "cell_maxdiam", "cell_mindiam", "nucleus_length",
             "nucleus_circularity", "nucleus_solidity", "nucleus_maxdiam",
             "nucleus_mindiam", "nc_area_ratio", "dist_islet", "dist_tissue"]
MORPHO = ["cell_area", "nucleus_area", "cell_length", "cell_circularity",
          "cell_solidity", "cell_maxdiam", "cell_mindiam", "nucleus_length",
          "nucleus_circularity", "nucleus_solidity", "nucleus_maxdiam",
          "nucleus_mindiam", "nc_area_ratio", "dist_islet", "dist_tissue"]
STATUS_MAP = {"ND": "ND", "AAB": "Aab+", "AAb": "Aab+", "AAb+": "Aab+", "T1D": "T1D"}


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", type=Path, default=OUT)
    a = ap.parse_args(argv)

    ds = pads.dataset(str(CELLS), format="parquet", partitioning="hive")
    # Marker INTERSECTION across donor partitions. The batches differ by one marker
    # (batch-1 has CD38, batch-2 has b-Catenin1), so the unified dataset schema is a
    # 60-marker UNION with a null column per cohort — which would NaN-poison X and break
    # scVI. Take only markers present in EVERY partition -> uniform 58-marker panel.
    frag_markers = [set(f.physical_schema.names) - set(META_COLS) for f in ds.get_fragments()]
    common = set.intersection(*frag_markers) if frag_markers else set()
    union = set.union(*frag_markers) if frag_markers else set()
    dropped = sorted(union - common)
    markers = [c for c in ds.schema.names if c in common and c not in META_COLS]
    if dropped:
        print(f"[panel] dropping batch-specific markers {dropped} -> {len(markers)} common markers")
    print(f"[info] {len(markers)} markers; loading parquet ...")
    keep_cols = [c for c in ds.schema.names if c in markers or c in META_COLS]
    tbl = ds.to_table(columns=keep_cols)          # exclude the batch-specific null columns
    # cast markers to float32 in Arrow (halves memory before to_pandas)
    for m in markers:
        i = tbl.schema.get_field_index(m)
        tbl = tbl.set_column(i, m, pc.cast(tbl[m], pa.float32()))
    df = tbl.to_pandas()
    del tbl
    print(f"[info] {len(df):,} cells x {df.shape[1]} cols")

    X = df[markers].to_numpy(dtype=np.float32)
    obs = pd.DataFrame(index=df["object_id"].astype(str).values)
    obs.index.name = None
    obs["imageid"] = df["donor_id"].astype(str).values
    obs["image"] = df["image"].astype(str).values
    obs["Object ID"] = df["object_id"].astype(str).values
    obs["Parent"] = df["parent_raw"].astype(str).values
    obs["Classification"] = df["classification"].astype(str).values
    obs["X_centroid"] = df["X_centroid"].to_numpy()
    obs["Y_centroid"] = df["Y_centroid"].to_numpy()
    obs["cell_region"] = df["cell_region"].astype(str).values
    obs["islet_num"] = df["islet_num"].astype(str).values
    for c in MORPHO:
        obs[c] = df[c].to_numpy()

    meta = pd.read_excel(DONOR_META)
    meta["donor_id"] = meta["donor_id"].astype(str).str.strip()
    mm = meta.set_index("donor_id")
    di = pd.Series(obs["imageid"].values, index=obs.index)
    obs["Donor Status"] = di.map(mm["disease.status"]).astype(str).map(
        lambda v: STATUS_MAP.get(v, v)).values
    for src, dst in [("age", "Age"), ("sex", "Sex"), ("aab", "AAb"),
                     ("aab_count", "AAb_count"), ("duration", "Duration")]:
        obs[dst] = di.map(mm[src]).values if src in mm.columns else np.nan

    adata = ad.AnnData(X=X, obs=obs, var=pd.DataFrame(index=pd.Index(markers)))
    adata.uns["all_markers"] = list(markers)
    a.out.parent.mkdir(parents=True, exist_ok=True)
    adata.write_h5ad(a.out)
    print("[done]", a.out, adata.shape)
    print("Donor Status:", obs["Donor Status"].value_counts().to_dict())
    print("cell_region:", obs["cell_region"].value_counts().to_dict())
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
