"""Simplest direct check: are follicle area distributions visibly different by rs3184504?

Pool every follicle from both cohorts, color by genotype, plot. No aggregation,
no per-test framework. If the answer isn't obvious here, it isn't there.
"""
from pathlib import Path
import sys

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

CODE_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(CODE_ROOT / "analysis"))
import data_utils as du

GENO_ORDER = ["C/C", "C/T", "T/T"]
PALETTE = du.GENO_PALETTE

# --- H&E cohort (per-follicle) ---------------------------------------------
he = pd.read_csv(CODE_ROOT / "27_donor_spleen_measurements1.csv")
he = he[he["Classification"] == "Follicle"].copy()
he = he[he["Area µm^2"] <= 5_000_000]  # H18's outlier filter
he["Sample"] = he["Image"].apply(du.extract_sample_id)

he_geno = pd.read_excel(CODE_ROOT / "Spleen_rs3184504_Genotypes.xlsx")
geno_col = next(c for c in he_geno.columns if "rs3184504" in c)
he_geno_map = dict(zip(he_geno["Sample ID"].astype(str).str.strip(),
                        he_geno[geno_col].apply(du._normalize_genotype)))
he["Genotype"] = he["Sample"].map(he_geno_map)
he = he.dropna(subset=["Genotype"])
he["Cohort"] = "H&E"

# --- Fluorescent cohort (per-follicle) -------------------------------------
flu = pd.read_csv(du.DATA_ROOT / "Measurements" / "AnnotationsFinal.csv")
flu = flu[flu["Classification"] == "Follicle"].copy()
flu["Sample"] = flu["Image"].apply(du.extract_sample_id)
flu_geno_map = du._build_genotype_map()
flu["Genotype"] = flu["Sample"].map(flu_geno_map)
flu = flu.dropna(subset=["Genotype"])
# Use the area column from this CSV
area_col = next(c for c in flu.columns if c.startswith("Area"))
flu = flu.rename(columns={area_col: "Area µm^2"})
flu["Cohort"] = "Fluorescent"

# Combine
keep = ["Sample", "Genotype", "Cohort", "Area µm^2"]
combined = pd.concat([he[keep], flu[keep]], ignore_index=True)
combined["Genotype"] = pd.Categorical(combined["Genotype"], categories=GENO_ORDER, ordered=True)

print(f"Total follicles pooled: {len(combined):,}")
print(combined.groupby(["Cohort", "Genotype"], observed=True).agg(
    n_follicles=("Area µm^2", "size"),
    n_donors=("Sample", "nunique"),
    median_area=("Area µm^2", "median"),
    mean_area=("Area µm^2", "mean"),
).round(0))
print()

# --- The basic plot --------------------------------------------------------
du.setup_style()
fig, axes = plt.subplots(1, 2, figsize=(14, 5), sharey=True)

logarea_all = np.log10(combined["Area µm^2"])
x_grid = np.linspace(logarea_all.min(), logarea_all.max(), 400)

for ax, cohort in zip(axes, ["H&E", "Fluorescent"]):
    sub = combined[combined["Cohort"] == cohort]
    for geno in GENO_ORDER:
        s = sub[sub["Genotype"] == geno]["Area µm^2"]
        if len(s) < 5:
            continue
        ax.hist(np.log10(s), bins=60, density=True, color=PALETTE[geno],
                alpha=0.45, label=f"{geno} (n={len(s)} foll, "
                                  f"{sub[sub['Genotype']==geno]['Sample'].nunique()} donors)")
        ax.axvline(np.log10(s.median()), color=PALETTE[geno], lw=2.2, ls="--")
    ax.set_title(f"{cohort} cohort")
    ax.set_xlabel("log10(Follicle Area, µm²)")
    ax.set_ylabel("Density")
    ax.legend(fontsize=8, title="Genotype (median = dashed)")

fig.suptitle("Basic check — pooled follicle area distributions by rs3184504",
             y=1.02, fontsize=13)
plt.tight_layout()
du.FIGURES_DIR.mkdir(parents=True, exist_ok=True)
out = du.FIGURES_DIR / "basic_follicle_area_by_genotype.png"
fig.savefig(out, dpi=150, bbox_inches="tight", facecolor="white")
print(f"Saved: {out}")
plt.close(fig)

# --- Honest one-line ratio check ------------------------------------------
print("\nMedian follicle area ratio (T/T vs C/C, by cohort):")
for cohort in ["H&E", "Fluorescent"]:
    sub = combined[combined["Cohort"] == cohort]
    cc = sub[sub["Genotype"] == "C/C"]["Area µm^2"].median()
    tt = sub[sub["Genotype"] == "T/T"]["Area µm^2"].median()
    print(f"  {cohort:12s}  C/C median = {cc:>9,.0f} µm²   "
          f"T/T median = {tt:>9,.0f} µm²   ratio T/T / C/C = {tt/cc:.2f}")
