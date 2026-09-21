#!/usr/bin/env python3
"""One AnnData per donor for the TissUUmaps viewer, from the CURRENT typing build.

WHY H5AD AND NOT CSV. TissUUmaps 3.2 reads AnnData natively (`tissuumaps/read_h5ad.py`), and one
file then carries three things a CSV cannot carry together: the cell's call, the marker matrix it
was called from, and the per-object overlap the structure rule used. That is what makes the calls
checkable BY EYE against the image rather than against a count -- colour by `cell_class`, switch to
any marker, and the same cells stay selected.

WHAT IS IN IT
  X          the aligned per-cell marker matrix (`celltype/uniform`), var_names = marker names.
             These are the values the classification actually read.
  obsm       `spatial` = (X_centroid, Y_centroid) exactly as `typing/cells` carries them, which
             is QUPATH RESOLUTION #1 (0.5078 um/px) -- NOT micrometres, and NOT the full-res
             pixel space the DZI tiles are in (0.2539 um/px). A consumer must scale by the ratio
             of the two declared pixel sizes, which is 2.0.

             CORRECTED 2026-09-21. This block previously said "IMAGE PIXELS, the space the DZI
             tiles are in", which is wrong by exactly that factor of 2. Measured on five donors
             by comparing the spatial extent against the channel tile dimensions: ratios
             1.97-2.08 in both axes, tightest on 6374 where the tissue nearly fills the slide
             (1.974 / 1.970). A viewer that trusted the old sentence would have placed every
             cell at half its true position -- the whole slide in the top-left quadrant, which
             reads as unremarkable until you zoom to cell level.

             The COORDINATES ARE UNCHANGED and were never wrong; they are what typing/cells
             holds. Only the description of them was. `scripts/senior/build_tissuumaps_project.py`
             derives the factor from the two pixel sizes and falsifies the scaled extent against
             the image bounds before attaching.
  obs        `cell_class` and how it was reached -- `resolved_by`, `object_context`, `islet_route`,
             `islet_stratum`, `excluded_reason` -- plus the three `crossed_frac_*` fractions, so a
             structure call can be inspected against the object it came from and a passthrough
             shows HOW MUCH of the cell that object covers.

WHAT IS NOT IN IT. Nothing is recomputed, nothing is thresholded, no cell is dropped. This is a
view of `typing/cells/`, joined by `object_id` and nothing else.

    .venv/bin/python scripts/senior/build_tissuumaps_h5ad.py [--donors 115 6442]
"""
from __future__ import annotations

import argparse
import pathlib
import sys

import numpy as np
import pandas as pd
import pyarrow.parquet as pq

REPO = pathlib.Path(__file__).resolve().parents[2]
CELLS = REPO / "data/development/typing/cells"
UNIFORM = REPO / "data/development/celltype/uniform"
CROSSED = REPO / "data/development/typing/crossed"
OUT = pathlib.Path("/data/tissuumaps")
#: Where the per-donor projects and their tiled channels already live.
TMAP_DIR = pathlib.Path("/data/tissuumaps")

#: Columns that say WHAT the cell is and HOW that was decided. Carried whenever present, because
#: a call the viewer cannot explain is a call the eye cannot check.
OBS_WANTED = ("cell_class", "resolved_by", "object_context", "cell_region", "cell_region_20um",
              "islet_route", "islet_stratum", "islet_num", "islet_cells", "excluded_reason",
              "immune_definition", "top_expressed_marker", "top_identity_marker",
              "passed_through_Blood_vessel", "passed_through_Nerve", "passed_through_Lymphatic")

#: Never markers: identity, geometry and QC columns that live in `obs`, not in `X`.
NOT_A_MARKER = {"donor_id", "object_id", "qc_analysis_eligible", "X_centroid", "Y_centroid",
                "cell_region", "islet_num", "cell_region_20um", "islet_key", "phenotype"}


def build(donor: str) -> pathlib.Path | None:
    cells_path = CELLS / f"donor_id={donor}/data.parquet"
    uni_path = UNIFORM / f"donor_id={donor}/data.parquet"
    if not cells_path.exists():
        print(f"  {donor}: no typing partition -- skipped", flush=True)
        return None
    import anndata as ad

    have = set(pq.ParquetFile(cells_path).schema_arrow.names)
    cols = ["object_id", "X_centroid", "Y_centroid"] + [c for c in OBS_WANTED if c in have]
    cells = pq.read_table(cells_path, columns=cols).to_pandas()

    # The marker matrix the classification read. Absent for a donor is UNAVAILABLE, not zero: the
    # file is still written, with no X, so the calls remain viewable without inventing values.
    markers: list[str] = []
    if uni_path.exists():
        uhave = [n for n in pq.ParquetFile(uni_path).schema_arrow.names if n not in NOT_A_MARKER]
        uni = pq.read_table(uni_path, columns=["object_id"] + uhave).to_pandas()
        cells = cells.merge(uni, on="object_id", how="left", validate="one_to_one")
        markers = uhave
    else:
        print(f"  {donor}: no aligned marker matrix -- writing calls only", flush=True)

    if CROSSED.joinpath(f"crossed_{donor}.parquet").exists():
        cf = pq.read_table(CROSSED / f"crossed_{donor}.parquet").to_pandas()
        cf = cf.drop(columns=[c for c in ("donor_id",) if c in cf.columns])
        cells = cells.merge(cf, on="object_id", how="left", validate="one_to_one")

    obs = pd.DataFrame(index=pd.Index(cells["object_id"].astype(str), name="object_id"))
    for c in cells.columns:
        if c in ("object_id", "X_centroid", "Y_centroid") or c in markers:
            continue
        # BY VALUE, never by Series. `obs` is indexed by object_id and `cells` by position, so
        # assigning the Series lets pandas align the two indexes and silently write all-NaN --
        # a file that loads, shows every cell, and has no labels on any of them.
        v = cells[c].to_numpy()
        # Categorical for the label columns TissUUmaps groups by; it reads them from
        # /obs/<col>/categories, and a plain object column there is not selectable as a grouping.
        obs[c] = pd.Categorical(v) if cells[c].dtype == object else v

    X = (cells[markers].to_numpy(dtype=np.float32) if markers
         else np.zeros((len(cells), 0), dtype=np.float32))
    adata = ad.AnnData(X=X, obs=obs,
                       var=pd.DataFrame(index=pd.Index(markers, name="marker")))
    adata.obsm["spatial"] = cells[["X_centroid", "Y_centroid"]].to_numpy(dtype=np.float64)
    # The donor's EXISTING tmap project -- its 59 image layers, opacities and filters -- carried
    # inside the file. TissUUmaps `read_h5ad.h5ad_to_tmap` uses `/uns/tmap` as the base project
    # and appends the marker files it generates, so one h5ad opens with the images AND the cells
    # and there is no second file to keep in step. The layer tileSources are relative to the tmap
    # directory, which is why the h5ad is written beside them.
    tmap_src = TMAP_DIR / f"{donor}.tmap"
    if tmap_src.exists():
        import json as _json
        base = _json.loads(tmap_src.read_text())
        base.pop("markerFiles", None)          # regenerated from this file's own obs
        base["filename"] = donor
        adata.uns["tmap"] = _json.dumps(base)
        # How the observation menu is grouped, so the call and the reasoning are not buried in one
        # alphabetical list with 58 markers.
        adata.uns["tmap_obsgroups"] = _json.dumps({
            "The call": [c for c in ("cell_class", "immune_definition") if c in obs],
            "How it was decided": [c for c in ("resolved_by", "islet_route", "islet_stratum",
                                               "excluded_reason", "top_expressed_marker",
                                               "top_identity_marker") if c in obs],
            "Where the cell is": [c for c in ("object_context", "cell_region", "cell_region_20um",
                                              "islet_num", "islet_cells") if c in obs],
            "How much an object covers it": [
                c for c in obs if c.startswith(("crossed_frac_", "passed_through_"))],
        })
    adata.uns["source"] = {
        "typing": str(cells_path.relative_to(REPO)),
        "markers": str(uni_path.relative_to(REPO)) if markers else "UNAVAILABLE",
        "coordinates": "QuPath resolution #1 (0.5078 um/px), from typing/cells X_centroid/Y_centroid; multiply by 2.0 to reach the 0.2539 um/px full-res pixel space the DZI tiles use",
        "note": "a view of the typing build; nothing recomputed, nothing thresholded, no cell "
                "dropped. Joined by object_id only.",
    }
    OUT.mkdir(parents=True, exist_ok=True)
    target = OUT / f"{donor}.h5ad"
    adata.write_h5ad(target, compression="gzip")
    _flatten_string_arrays(target)
    n = obs["cell_class"].value_counts() if "cell_class" in obs else pd.Series(dtype=int)
    print(f"  {donor}: {adata.n_obs:>9,} cells x {adata.n_vars} markers -> {target.name}"
          + ("   " + "  ".join(f"{k} {v:,}" for k, v in n.head(4).items()) if len(n) else ""),
          flush=True)
    return target


def _flatten_string_arrays(path: pathlib.Path) -> list[str]:
    """Write string columns as PLAIN DATASETS, which is what the viewer can read.

    anndata 0.13 encodes a string column as a GROUP -- ``encoding-type:
    nullable-string-array``, holding ``values`` and ``mask``.  TissUUmaps 3.2.1
    (`read_h5ad.getVarList`, line 66) does ``[x.decode("utf-8") for x in adata.get("/var/<index>")]``
    and an h5py Group iterates its MEMBER NAMES, so it silently reads ``["values", "mask"]`` as if
    those were the markers, and then fails on ``str.decode``.

    Rewriting the group as the plain dataset the older encoding used loses nothing here: the mask
    is all-False on every column this builder writes, which is checked before the rewrite and
    raises rather than dropping a null on the floor.  Categoricals are untouched -- they are
    already plain ``categories``/``codes`` datasets and the viewer reads them correctly.
    """
    import h5py

    rewritten: list[str] = []
    with h5py.File(path, "r+") as f:
        def walk(group, prefix=""):
            for key in list(group.keys()):
                item = group[key]
                if not isinstance(item, h5py.Group):
                    continue
                if item.attrs.get("encoding-type") == "nullable-string-array":
                    mask = item["mask"][...]
                    if mask.any():
                        raise ValueError(
                            f"{path.name}: {prefix}/{key} has {int(mask.sum())} null strings; "
                            f"flattening would turn a missing value into the na-value, and a "
                            f"missing measurement is unavailable, not a label")
                    values = item["values"][...]
                    del group[key]
                    ds = group.create_dataset(key, data=values)
                    ds.attrs["encoding-type"] = "string-array"
                    ds.attrs["encoding-version"] = "0.2.0"
                    rewritten.append(f"{prefix}/{key}")
                else:
                    walk(item, f"{prefix}/{key}")
        walk(f)
    return rewritten

def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--donors", nargs="*", default=None)
    args = ap.parse_args()
    donors = args.donors or sorted(p.name.split("=")[1] for p in CELLS.glob("donor_id=*"))
    for d in donors:
        build(d)
    print(f"\n-> {OUT}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
