#!/usr/bin/env python
"""Validate the cell-derived region layer against QuPath's own parenting.

The two are independent: QuPath assigns a vessel to a region by hierarchy at
segmentation time, the region layer by which cell-derived polygon contains the
vessel centroid. Agreement is therefore evidence about both.

The headline case is 1901HBMP004, whose corrected export re-parents 850
SmallVessels to Follicle where the polygon method finds ~881 — 3.5% apart. That
number only means something against the spread across donors that have both,
which is what this script measures.

Outputs
-------
analysis/tables/region_layer_validation.csv     per donor per region
analysis/tables/region_layer_identifiability.csv  platform/genotype design checks

Usage
-----
    python scripts/validate_region_layer.py
"""
import sys
import warnings
from pathlib import Path

warnings.filterwarnings("ignore")
import numpy as np
import pandas as pd

CODE_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(CODE_ROOT / "analysis"))
import data_utils as du  # noqa: E402

FIXED_1901 = du.DATA_ROOT / "Measurements" / "CellsFixed1901HBMP004.csv"


def qupath_counts(df):
    """SmallVessel counts per image per region from the QuPath Parent string.

    `from_parent=True` is load-bearing, not decoration. get_vessels defaults to
    the region layer now, so without it this function returns the very numbers
    it is supposed to be an independent check on, and the comparison below
    reports 0.0% disagreement everywhere — a passing validation that validates
    nothing.
    """
    v = du.get_vessels(df, from_parent=True)
    return (v[v["Region"].isin(du.MAIN_REGIONS)]
            .groupby(["Image", "Region"], observed=True).size()
            .rename("QuPath_Count").reset_index())


def main():
    layer = du.load_region_layer()
    df = du.load_all_data()

    # The corrected export repairs 1,014 vessels orphaned at the image root.
    # Under the pivot it is no longer the primary input — the region layer never
    # reads Parent — so it is used here purely as the reference to compare to.
    if FIXED_1901.exists():
        fix = pd.read_csv(FIXED_1901)
        fix["Sample"] = fix["Image"].apply(du.extract_sample_id)
        fix["Genotype"] = fix["Sample"].map(du.GENOTYPE_MAP)
        fix["Platform"] = "CODEX"
        imgs = set(fix["Image"])
        df = pd.concat([df[~df["Image"].isin(imgs)], fix], ignore_index=True)
        print(f"applied corrected export for {sorted(imgs)}")

    qp = qupath_counts(df)
    cell = du.compute_density_from_regions(df, region_layer=layer,
                                           regions=du.MAIN_REGIONS)
    if cell.empty:
        print("no region layer found — run scripts/build_region_layer_qupath.py")
        return

    m = cell.merge(qp, on=["Image", "Region"], how="left")
    m["QuPath_Count"] = m["QuPath_Count"].fillna(0).astype(int)
    denom = m[["Vessel_Count", "QuPath_Count"]].max(axis=1).replace(0, np.nan)
    m["Abs_Pct_Diff"] = 100 * (m["Vessel_Count"] - m["QuPath_Count"]).abs() / denom

    out = du.TABLES_DIR / "region_layer_validation.csv"
    m.to_csv(out, index=False)

    pd.set_option("display.width", 200)
    print("\n=== Vessel counts: cell-derived polygons vs QuPath Parent ===")
    show = m[["Sample", "Region", "Vessel_Count", "QuPath_Count", "Abs_Pct_Diff"]]
    print(show.to_string(index=False, float_format=lambda v: f"{v:.1f}"))

    fol = m[m["Region"] == "Follicle"].dropna(subset=["Abs_Pct_Diff"])
    if len(fol):
        print(f"\nFollicle agreement across {len(fol)} donors: "
              f"median {fol['Abs_Pct_Diff'].median():.1f}%, "
              f"IQR {fol['Abs_Pct_Diff'].quantile(.25):.1f}-"
              f"{fol['Abs_Pct_Diff'].quantile(.75):.1f}%, "
              f"range {fol['Abs_Pct_Diff'].min():.1f}-{fol['Abs_Pct_Diff'].max():.1f}%")
        one = fol[fol["Sample"] == "1901HBMP004"]
        if len(one):
            d = one["Abs_Pct_Diff"].iloc[0]
            pct = 100 * (fol["Abs_Pct_Diff"] <= d).mean()
            print(f"  1901HBMP004 sits at {d:.1f}% — {pct:.0f}th percentile of that spread")

    # ---- identifiability -------------------------------------------------
    donors = (df.drop_duplicates("Sample")[["Sample", "Platform", "Genotype"]]
              .dropna(subset=["Genotype"]))
    donors = donors[~donors["Sample"].isin(du.EXCLUDE_SAMPLES)]
    ct = pd.crosstab(donors["Platform"], donors["Genotype"])
    print("\n=== Identifiability ===")
    print(ct.to_string())
    dosage = donors["Genotype"].map({"C/C": 0, "C/T": 1, "T/T": 2}).astype(float)
    plat = (donors["Platform"] == "CODEX").astype(float)
    r = np.corrcoef(plat, dosage)[0, 1]
    # VIF of platform against dosage in a two-term design
    vif = 1 / (1 - r ** 2) if abs(r) < 1 else np.inf
    X = np.column_stack([np.ones(len(donors)), plat, dosage])
    print(f"  n donors                 = {len(donors)}")
    print(f"  corr(platform, dosage)   = {r:+.3f}   (was +0.638 with CODEX T/T-only)")
    print(f"  VIF                      = {vif:.2f}")
    print(f"  design matrix rank       = {np.linalg.matrix_rank(X)} of {X.shape[1]}")
    print(f"  CODEX genotypes present  = {sorted(donors[plat==1]['Genotype'].unique())}")
    pd.DataFrame([{"n_donors": len(donors), "corr_platform_dosage": r, "vif": vif,
                   "rank": int(np.linalg.matrix_rank(X)), "n_terms": X.shape[1]}]
                 ).to_csv(du.TABLES_DIR / "region_layer_identifiability.csv", index=False)
    print(f"\nwrote {out}")


if __name__ == "__main__":
    main()
