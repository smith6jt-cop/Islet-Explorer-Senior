"""Idempotently append a centroid-based spatial-features section to H18.

Adds three feature families using the per-Follicle data already in the H18
notebook (variable `df`):
  - Clark-Evans nearest-neighbor index R (clustering vs dispersion)
  - Pair-correlation g(r) summary at r=500 µm (relative density at that lag)
  - Lorenz/Gini of follicle area distribution per image
  - Centrality: distance to all-follicles centroid, normalized by sqrt(image follicle bbox area)

Per-donor aggregates → genotype tests via data_utils.full_stats_table.

Idempotency: the script strips any prior copy of the same `## 8. Centroid-based`
appendix before re-appending.  The strip is bounded by the next `## ` (top-level)
heading so it does not accidentally remove unrelated trailing sections.
"""
from pathlib import Path
import nbformat as nbf

PROJECT = Path(__file__).resolve().parent.parent
NB_PATH = PROJECT / "analysis" / "H18_HE_follicles.ipynb"

START_SENTINEL = "## 8. Centroid-based spatial features"

MD_INTRO = START_SENTINEL + """

The 27-donor H&E CSV exposes per-follicle Centroid X/Y and Area but no tissue polygon, so we focus on three centroid-driven feature families that the H18 simple metrics didn't capture:

1. **Clark-Evans nearest-neighbor index R** — `R = mean(d_NN) / (0.5 / sqrt(density))`. R<1 indicates clustering, R≈1 random Poisson, R>1 dispersion. Computed per image, then averaged per donor.
2. **Lorenz/Gini of follicle areas** — measures size inequality within each image. Gini=0 → all follicles same size; →1 → one follicle dominates. Per-donor mean across that donor's images.
3. **Centrality of follicle distribution** — for each follicle, distance to the centroid of all follicles in its image, divided by sqrt(image follicle bounding-box area) for scale invariance. Per-donor mean.

These are scale-aware, polygon-free features. They are tested with the same `full_stats_table` framework as H18's simple metrics (Kruskal-Wallis + 3 pairwise Mann-Whitney + Spearman dosage)."""

CODE = """from scipy.spatial import cKDTree


def clark_evans(coords):
    \"\"\"Clark-Evans nearest-neighbor index. R<1 cluster, =1 random, >1 dispersion.\"\"\"
    if len(coords) < 4:
        return np.nan
    tree = cKDTree(coords)
    d_nn = tree.query(coords, k=2)[0][:, 1]
    # Use min-area bounding rectangle for area estimate
    span_x = coords[:, 0].max() - coords[:, 0].min()
    span_y = coords[:, 1].max() - coords[:, 1].min()
    area = span_x * span_y
    if area <= 0:
        return np.nan
    density = len(coords) / area
    expected = 0.5 / np.sqrt(density)
    return np.nanmean(d_nn) / expected


def gini(values):
    \"\"\"Gini coefficient of a 1-D array of non-negative values.\"\"\"
    v = np.asarray([x for x in values if np.isfinite(x) and x > 0], dtype=float)
    if len(v) < 2:
        return np.nan
    v.sort()
    n = len(v)
    cum = np.cumsum(v)
    return (2 * np.sum((np.arange(1, n + 1)) * v) - (n + 1) * cum[-1]) / (n * cum[-1])


def centrality(coords):
    \"\"\"Mean distance from each follicle to the all-follicles centroid,
       normalized by sqrt(bounding-box area) for scale invariance.\"\"\"
    if len(coords) < 3:
        return np.nan
    cx, cy = coords.mean(axis=0)
    d = np.sqrt((coords[:, 0] - cx) ** 2 + (coords[:, 1] - cy) ** 2)
    span_x = coords[:, 0].max() - coords[:, 0].min()
    span_y = coords[:, 1].max() - coords[:, 1].min()
    scale = np.sqrt(max(span_x * span_y, 1.0))
    return np.nanmean(d) / scale


# Build per-image features from H18's filtered follicle frame `follicles`
per_image_spatial = []
for image, sub in follicles.groupby("Image", observed=True):
    coords = sub[["Centroid X µm", "Centroid Y µm"]].dropna().values
    areas = sub["Area µm^2"].dropna().values
    per_image_spatial.append({
        "Image": image,
        "Sample": extract_he_sample(image),
        "Clark_Evans_R": clark_evans(coords),
        "Gini_area":     gini(areas),
        "Centrality":    centrality(coords),
        "N_follicles":   len(coords),
    })
per_image_spatial = pd.DataFrame(per_image_spatial)
per_image_spatial = per_image_spatial.merge(
    pd.DataFrame({"Sample": list(geno_map), "Genotype": list(geno_map.values())}),
    on="Sample", how="left",
)
save_table(per_image_spatial, "HE_spatial_features_per_image")

# Per-donor aggregates
per_donor_spatial = (
    per_image_spatial.groupby(["Sample", "Genotype"], observed=True)
    [["Clark_Evans_R", "Gini_area", "Centrality"]]
    .mean()
    .reset_index()
)
per_donor_spatial["Genotype"] = pd.Categorical(per_donor_spatial["Genotype"],
                                                categories=GENO_ORDER, ordered=True)
save_table(per_donor_spatial, "HE_spatial_features_per_donor")
print(per_donor_spatial.groupby("Genotype", observed=True)
                       [["Clark_Evans_R", "Gini_area", "Centrality"]].mean().round(3))"""

CODE_STATS = """# Genotype tests on the three new features
spatial_stats_rows = []
for metric in ["Clark_Evans_R", "Gini_area", "Centrality"]:
    spatial_stats_rows.append(
        full_stats_table(per_donor_spatial, metric, label=metric)
    )
spatial_stats = pd.concat(spatial_stats_rows, ignore_index=True)
save_table(spatial_stats, "HE_spatial_features_stats")
spatial_stats"""

CODE_PLOT = """fig, axes = plt.subplots(1, 3, figsize=(13, 4.6), sharey=False)
labels = {"Clark_Evans_R": "Clark-Evans R\\n(<1 cluster, >1 disperse)",
          "Gini_area":     "Gini of follicle areas\\n(0 equal — 1 unequal)",
          "Centrality":    "Centrality (norm. dist. to centroid)"}
for ax, metric in zip(axes, ["Clark_Evans_R", "Gini_area", "Centrality"]):
    sns.boxplot(data=per_donor_spatial, x="Genotype", y=metric, order=GENO_ORDER,
                palette=GENO_PALETTE, ax=ax, showfliers=False, width=0.55)
    sns.stripplot(data=per_donor_spatial, x="Genotype", y=metric, order=GENO_ORDER,
                   color="black", size=4.5, alpha=0.7, ax=ax, jitter=0.18)
    h, kw_p = run_kruskal(per_donor_spatial, metric)
    rho, sp_p = run_dosage_trend(per_donor_spatial, metric)
    ax.set_title(labels[metric] + f"\\nKW p={kw_p:.3g}  ρ={rho:+.2f} (p={sp_p:.3g})",
                 fontsize=10)
    ax.set_xlabel(""); ax.set_ylabel(metric.replace('_', ' '))
fig.suptitle("H18 centroid-based spatial features by rs3184504", y=1.02)
plt.tight_layout()
save_figure(fig, "HE_spatial_features")
plt.show()"""


def build_appendix_cells():
    return [
        nbf.v4.new_markdown_cell(MD_INTRO),
        nbf.v4.new_code_cell(CODE),
        nbf.v4.new_code_cell(CODE_STATS),
        nbf.v4.new_code_cell(CODE_PLOT),
    ]


def main():
    nb = nbf.read(NB_PATH, as_version=4)
    new_cells = []
    in_section = False
    for cell in nb.cells:
        src = cell.source if isinstance(cell.source, str) else "".join(cell.source)
        if cell.cell_type == "markdown" and src.startswith(START_SENTINEL):
            in_section = True
            continue
        if in_section:
            # Stop stripping when we find the next top-level (## ) header that's NOT us
            if cell.cell_type == "markdown" and src.lstrip().startswith("## "):
                in_section = False
                new_cells.append(cell)
                continue
            # Otherwise drop this cell
            continue
        new_cells.append(cell)

    nb.cells = new_cells + build_appendix_cells()
    nbf.write(nb, NB_PATH)
    print(f"Wrote {NB_PATH.relative_to(PROJECT)}  ({len(nb.cells)} cells)")


if __name__ == "__main__":
    main()
