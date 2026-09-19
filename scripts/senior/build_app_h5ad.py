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

3D Senior: assemble the single enriched app H5AD that the Shiny app loads,
mirroring Junior's data/islet_explorer.h5ad schema (groovy .uns sheets +
enriched .obs) but built from the Senior islet-level datasets (NO groovy export).

Inputs (data/):
  islets_core_senior.h5ad   -- 9,130 core islets: obs prop_*/count_* (25 phenos),
                               demographics, areas/diam, dpt/DC/leiden, posfrac_*;
                               obsm X_umap / marker_pos_frac (59) / phenotype_*
  islets_peri_senior.h5ad   -- 9,130 peri islets (same index): marker_pos_frac, total_cells
  islet_targets.parquet     -- core/peri area_um2, counts, diameter (by islet_id)
  islet_composition.parquet -- per-core-islet hormone coexpression counts
  neighborhood_metrics.csv  -- 108 peri-islet metrics keyed by combined_islet_id

Output:
  data/islet_explorer_senior.h5ad

Groovy sheets reproduced (consumed by prep_data() via reconstruct_groovy_df):
  targets     : per islet x {core,band} x class{Nerve,Capillary,Lymphatic}
                region_um2 = QuPath annotation area (drives islet_diam_um);
                area_um2/area_density/count = NaN (structural intersection is a
                deferred v1.1 step — see results/*.geojson).
  markers     : per islet x {core,band} x 59 markers — n_cells, pos_count, pos_frac
                (from obsm marker_pos_frac; richer than Junior's 4 markers).
  composition : per core islet — cells_total + hormone single/any/multi/triple-neg.

Run in the `scimap` conda env (anndata 0.11).
"""
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd

SD = Path(__file__).resolve().parents[2] / "data"
OUT = SD / "islet_explorer_senior.h5ad"
AAB_COLS = ["AAb_GADA", "AAb_IA2A", "AAb_ZnT8A", "AAb_mIAA"]
TARGET_CLASSES = ["Nerve", "Capillary", "Lymphatic"]  # structural classes (areas deferred)


def store_groovy(adata, df, sheet):
    """Serialize a DataFrame into adata.uns as groovy_<sheet>_<col> arrays
    (mirrors scripts/build_h5ad_for_app.py so reconstruct_groovy_df reads it)."""
    adata.uns[f"groovy_{sheet}_columns"] = list(df.columns)
    adata.uns[f"groovy_{sheet}_n_rows"] = len(df)
    for col in df.columns:
        vals = df[col].values
        key = f"groovy_{sheet}_{col}"
        if pd.api.types.is_bool_dtype(vals):
            adata.uns[key] = vals.astype(bool)
        elif pd.api.types.is_numeric_dtype(vals):
            adata.uns[key] = np.asarray(vals, dtype=float)
        else:
            adata.uns[key] = np.asarray(vals, dtype=str)
    print(f"   groovy[{sheet}]: {len(df):,} rows x {len(df.columns)} cols")


def main():
    print("1. Loading Senior islet datasets")
    core = ad.read_h5ad(SD / "islets_core_senior.h5ad")
    peri = ad.read_h5ad(SD / "islets_peri_senior.h5ad")
    peri = peri[core.obs_names].copy()                      # align row order to core
    tgt = pd.read_parquet(SD / "islet_targets.parquet").reindex(core.obs_names)
    comp_src = pd.read_parquet(SD / "islet_composition.parquet")
    nbr = pd.read_csv(SD / "neighborhood_metrics.csv").set_index("combined_islet_id")
    print(f"   core={core.n_obs} peri={peri.n_obs} targets={len(tgt)} "
          f"comp={len(comp_src)} nbr={nbr.shape}")

    caseid = core.obs["imageid"].astype(str).to_numpy()      # donor id "115"
    status = core.obs["donor_status"].astype(str).to_numpy()  # ND / Aab+ / T1D
    base = core.obs["base_islet_id"].astype(str).to_numpy()   # "Islet_755"
    donor_status_map = dict(zip(caseid, status))

    # ── targets sheet ────────────────────────────────────────────────────────
    print("2. Building groovy 'targets'")
    region_specs = [
        ("islet_core", "_core", tgt["core_area_um2"].to_numpy()),
        ("islet_band", "_band", tgt["peri_area_um2"].to_numpy()),
    ]
    tgt_frames = []
    for rtype, suffix, reg_um2 in region_specs:
        b = pd.DataFrame({
            "Case ID": caseid,
            "Donor Status": status,
            "region": pd.Series(base).astype(str) + suffix,
            "type": rtype,
            "region_um2": reg_um2.astype(float),
        })
        for cls in TARGET_CLASSES:
            r = b.copy()
            r["class"] = cls
            r["area_um2"] = np.nan          # structural area within region (deferred)
            r["area_density"] = np.nan
            r["count"] = np.nan
            tgt_frames.append(r)
    targets_df = pd.concat(tgt_frames, ignore_index=True)
    for a in AAB_COLS:
        targets_df[a] = False
    targets_df = targets_df[["Case ID", "Donor Status", "region", "type", "class",
                             "area_um2", "region_um2", "area_density", "count"] + AAB_COLS]

    # ── markers sheet (59 markers x core/band) ───────────────────────────────
    print("3. Building groovy 'markers' (59 markers x 2 regions)")
    mk_names = list(core.uns["marker_pos_frac_names"])

    def marker_long(pf, ncells, rtype, suffix):
        df = pd.DataFrame(np.asarray(pf), columns=mk_names)
        df.insert(0, "Case ID", caseid)
        df.insert(1, "Donor Status", status)
        df.insert(2, "region", pd.Series(base).astype(str) + suffix)
        df.insert(3, "n_cells", np.asarray(ncells, dtype=float))
        long = df.melt(id_vars=["Case ID", "Donor Status", "region", "n_cells"],
                       var_name="marker", value_name="pos_frac")
        long["region_type"] = rtype
        long["pos_count"] = np.round(long["pos_frac"].astype(float)
                                     * long["n_cells"].astype(float))
        return long

    markers_df = pd.concat([
        marker_long(core.obsm["marker_pos_frac"], core.obs["total_cells"].to_numpy(),
                    "islet_core", "_core"),
        marker_long(peri.obsm["marker_pos_frac"], peri.obs["total_cells"].to_numpy(),
                    "islet_band", "_band"),
    ], ignore_index=True)
    for a in AAB_COLS:
        markers_df[a] = False
    markers_df = markers_df[["Case ID", "Donor Status", "region", "region_type",
                             "marker", "n_cells", "pos_count", "pos_frac"] + AAB_COLS]

    # ── composition sheet ────────────────────────────────────────────────────
    print("4. Building groovy 'composition'")
    comp_df = comp_src.copy()
    comp_df["Donor Status"] = comp_df["Case ID"].astype(str).map(donor_status_map)
    for a in AAB_COLS:
        comp_df[a] = False
    comp_cols = ["Case ID", "Donor Status", "name", "islet_id", "cells_total",
                 "Ins_single", "Glu_single", "Stt_single", "Multi_Pos", "Triple_Neg",
                 "Ins_any", "Glu_any", "Stt_any"] + AAB_COLS
    comp_df = comp_df[comp_cols]

    # ── enrich obs: neighborhood + leiden_umap + AAb flags ───────────────────
    print("5. Enriching obs (neighborhood + leiden_umap + AAb)")
    core.obs = core.obs.join(nbr)                            # by index == combined_islet_id
    n_peri = int(core.obs["total_cells_peri"].notna().sum()) if "total_cells_peri" in core.obs else 0
    print(f"   merged {nbr.shape[1]} neighborhood cols; islets with peri data: {n_peri}/{core.n_obs}")
    if "X_umap" in core.obsm:
        core.obs["leiden_umap_1"] = np.asarray(core.obsm["X_umap"])[:, 0]
        core.obs["leiden_umap_2"] = np.asarray(core.obsm["X_umap"])[:, 1]
    for a in AAB_COLS:
        core.obs[a] = False

    # Strip NaN from string/categorical obs cols so h5py can write them.
    for col in core.obs.columns:
        if isinstance(core.obs[col].dtype, pd.CategoricalDtype):
            core.obs[col] = core.obs[col].astype(str)
        elif core.obs[col].dtype == object:
            core.obs[col] = core.obs[col].fillna("").astype(str)

    # ── store + write ────────────────────────────────────────────────────────
    print("6. Storing groovy sheets in uns")
    store_groovy(core, targets_df, "targets")
    store_groovy(core, markers_df, "markers")
    store_groovy(core, comp_df, "composition")
    core.uns["data_provenance"] = {
        "pipeline": "Islet Explorer Senior v1 (3D app integration)",
        "n_target_rows": len(targets_df), "n_marker_rows": len(markers_df),
        "n_composition_rows": len(comp_df),
    }
    # Drop heavy obsm the app doesn't need from this file (keep X_umap for converter).
    for k in ["X_diffmap", "phenotype_counts", "marker_pos_frac", "morphology"]:
        core.obsm.pop(k, None)

    print(f"7. Writing {OUT}")
    core.write_h5ad(OUT)
    sz = OUT.stat().st_size / 1e6
    print(f"   wrote {sz:.1f} MB | shape {core.shape} | obs cols {core.obs.shape[1]} "
          f"| obsm {list(core.obsm.keys())}")
    print("[done]")


if __name__ == "__main__":
    main()
