#!/usr/bin/env python3
"""
3B(+) Senior: per-islet region areas / "targets" from AnnotationsMeasurements.csv
(NO groovy). QuPath annotation export gives Area µm^2 + Num Detections per
annotation; Classification Islet -> core, IsletExpanded -> peri. Computes
islet_diam_um from core area and merges into the islet-level h5ads.

islet_id = "<donorDigits>_<Name>" matches the islet h5ad obs_names
(imageid_Islet_N from create_islet_dataset_fixed).
"""
import re
from pathlib import Path
import numpy as np
import pandas as pd
import anndata as ad

BASE = Path.home() / "IO60panc2nd"
ANN = BASE / "AnnotationsMeasurements.csv"
SD = Path(__file__).resolve().parents[2] / "data"


def donor_of(img):
    m = re.search(r"[0-9]+", str(img))
    return m.group(0) if m else ""


def main():
    df = pd.read_csv(ANN)
    df["donor_id"] = df["Image"].map(donor_of)
    df["Classification"] = df["Classification"].astype(str)
    df["Name"] = df["Name"].astype(str)

    core = df[df["Classification"] == "Islet"].copy()
    peri = df[df["Classification"] == "IsletExpanded"].copy()
    core["islet_id"] = core["donor_id"] + "_" + core["Name"]
    peri["islet_id"] = peri["donor_id"] + "_" + peri["Name"].str.replace("_exp20um", "", regex=False)

    core_t = core.set_index("islet_id")[
        ["Area µm^2", "Num Detections", "Perimeter µm", "Max diameter µm", "Min diameter µm"]]
    core_t.columns = ["core_area_um2", "core_count", "core_perimeter_um", "core_maxdiam_um", "core_mindiam_um"]
    peri_t = peri.set_index("islet_id")[["Area µm^2", "Num Detections"]]
    peri_t.columns = ["peri_area_um2", "peri_count"]
    core_t = core_t[~core_t.index.duplicated()]
    peri_t = peri_t[~peri_t.index.duplicated()]

    tgt = core_t.join(peri_t, how="outer")
    tgt["islet_diam_um"] = 2 * np.sqrt(tgt["core_area_um2"].clip(lower=0) / np.pi)
    tgt["core_density"] = tgt["core_count"] / tgt["core_area_um2"]
    print(f"[targets] {tgt.shape[0]} islets; with core area: {tgt['core_area_um2'].notna().sum()}")
    tgt.to_parquet(SD / "islet_targets.parquet")

    for f in ["islets_core_senior.h5ad", "islets_peri_senior.h5ad"]:
        a = ad.read_h5ad(SD / f)
        for col in tgt.columns:
            a.obs[col] = tgt[col].reindex(a.obs_names).to_numpy()
        a.write_h5ad(SD / f)
        n = int(a.obs["core_area_um2"].notna().sum())
        print(f"[merge] {f}: {n}/{a.n_obs} islets matched; "
              f"islet_diam_um median={np.nanmedian(a.obs['islet_diam_um']):.1f} µm")


if __name__ == "__main__":
    main()
