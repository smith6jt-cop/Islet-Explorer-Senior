#!/usr/bin/env python
"""Primary endpoint: SmallVessel density per region per donor, by rs3184504.

This is the headline result, and until now it had no producer.
`plot_region_layer_results.py` *reads* both tables written here and draws them;
nothing wrote them. They were therefore whatever the last interactive session
happened to leave on disk, and after the region layer changed underneath them
the figures kept plotting the old numbers without anything failing. That is the
same defect class as `radial_pfz_width.csv`, which had a docstring claiming a
producer that did not exist.

Method
------
Vessel region is decided by containment in the region-layer polygons, region
area by the union of those polygons. Neither reads QuPath's `Parent` string,
which strands 1,014 SmallVessels at the image root — 21.5% of them in
`1901HBMP004` against under 0.3% in every other image — nor sums annotation
areas, which double-counts the overlapping CODEX `Red_Pulp`/`Sinusoid` pair.

Unit of analysis is the donor, not the annotation. Tests are the project
standard: Kruskal-Wallis across the three genotypes, Mann-Whitney with
rank-biserial r pairwise, and a Spearman dosage trend on C/C=0, C/T=1, T/T=2.

    python scripts/compute_primary_endpoint.py

Writes analysis/tables/P1_vessel_density_cell_regions.csv (per donor x region)
and analysis/tables/P1_vessel_density_stats.csv (the tests). Run
plot_region_layer_results.py afterwards to redraw the figures.
"""
import argparse
import sys
import warnings
from pathlib import Path

warnings.filterwarnings("ignore")
import pandas as pd

CODE_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(CODE_ROOT / "analysis"))
import data_utils as du  # noqa: E402


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--regions", nargs="*",
                    default=du.MAIN_REGIONS + ["Perifollicular"])
    args = ap.parse_args()

    layer = du.load_region_layer()
    df = du.load_all_data()

    density = du.compute_density_from_regions(df, region_layer=layer,
                                              regions=args.regions)
    if density.empty:
        raise SystemExit("no region layer — run scripts/build_region_layer_qupath.py")

    density["Platform"] = density["Sample"].apply(
        lambda s: "CODEX" if s in du.CODEX_SAMPLES else "Phenocycler")
    du.save_table(density, "P1_vessel_density_cell_regions")

    rows = []
    for region in args.regions:
        sub = density[density["Region"] == region].dropna(
            subset=["Density_per_mm2"])
        if sub["Sample"].nunique() < 3:
            print(f"{region}: only {sub['Sample'].nunique()} donors — skipped")
            continue
        stats = du.full_stats_table(sub, "Density_per_mm2", label=region)
        stats.insert(0, "Region", region)
        stats["n"] = sub["Sample"].nunique()
        rows.append(stats)

        kw = stats[stats["Test"] == "Kruskal-Wallis"].iloc[0]
        sp = stats[stats["Test"] == "Spearman dosage"].iloc[0]
        print(f"{region:<16} n={sub['Sample'].nunique():<3} "
              f"KW p={kw['p']:.3f}   dosage {sp['Effect_Size']}, p={sp['p']:.3f}")

    du.save_table(pd.concat(rows, ignore_index=True), "P1_vessel_density_stats")
    return density


if __name__ == "__main__":
    main()
