#!/usr/bin/env python
"""Build the region layer from the QuPath pixel-classifier exports.

The cell-derived route was tried and abandoned. Deriving Follicle/PALS from two
markers cannot match a classifier trained on 25-33 channels: comparing z-scored
CD20 and CD3e fields labels anything sitting below CD20's mean as PALS whether
or not CD3e is elevated, which gave HDL098 a "PALS" whose CD3e median (40) is
*lower* than its red pulp (71). These annotations are the reference.

Region sources
--------------
Follicle, PALS   QuPath geometry from analysis/geojson_regions/ (pixel coords,
                 scaled per donor by data_utils.get_pixel_size), falling back
                 per file to the older Follicle+PALS-only analysis/geojson/.
                 Peripheral_White_Pulp maps to Perifollicular, NOT to PALS —
                 see CLASS_TO_REGION below for why _CODEX_REGION_MAP is wrong
                 for this purpose.
Trabeculae       QuPath. Every donor in project.qpproj now has it, via
                 export_all_regions_geojson.groovy; before that export only the
                 three ForSH2B3 CODEX sections did, and everyone else's
                 trabecular area sat inside RedPulp.
RedPulp          tissue minus everything above. Tissue comes from cell density
                 (data_utils.tissue_extent), which reproduces QuPath's total
                 annotated area at 0.92-1.12x, median 1.00.
Perifollicular   ring immediately adjacent to the follicle, width measured per
                 donor by radial_follicle_profiles.py (the SMA peak just
                 outside the boundary, median 75 µm), unioned with the
                 classifier's Peripheral_White_Pulp where it exists. Anchored
                 on the follicle; PALS wins where the two overlap. PALS is
                 validated against the CSV before that trim.

Small QuPath follicle objects are fragments, not follicles: their area
distribution is bimodal and a 2-component GMM on log10 area puts the boundary
near 8,000 µm2 (the filter H6 applies). Below it sits about half the objects
but only ~2% of the follicle area. scripts/fit_area_floors.py re-derives both
this and the PALS floor — note the PALS one is a geometric floor rather than a
GMM crossover, for the reason recorded beside it in data_utils.

Every donor's Follicle and PALS area is checked against the annotation CSV, so
a pixel-size or coordinate error cannot pass silently.

    python scripts/build_region_layer_qupath.py
"""
import argparse
import json
import sys
import warnings
from pathlib import Path

warnings.filterwarnings("ignore")
import numpy as np
import pandas as pd
from shapely.affinity import scale as shapely_scale
from shapely.geometry import mapping, shape
from shapely.ops import unary_union

CODE_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(CODE_ROOT / "analysis"))
import data_utils as du  # noqa: E402

GEOJSON_DIRS = [du.OUTPUT_ROOT / "analysis" / "geojson_regions",  # preferred
                du.OUTPUT_ROOT / "analysis" / "geojson"]  # Follicle+PALS only
# QuPath class name -> canonical region.
#
# Peripheral_White_Pulp maps to Perifollicular, NOT to PALS. data_utils'
# _CODEX_REGION_MAP folds it into PALS, which is wrong for this purpose: the
# peripheral white pulp is the zone around the follicle, while PALS is the
# sheath around the central arteriole. Folding it in gave the CODEX donors a
# "PALS" collar hugging every follicle that then consumed the perifollicular
# band — HBMP006 had 2 of 7 follicles with a band, HDL098 1 of 8 — and it
# doubled their PALS area (HBMP006: 0.90 mm2 of PALS reported as 1.76).
CLASS_TO_REGION = {
    "Follicle": "Follicle",
    "PALS": "PALS",
    "Peripheral_White_Pulp": "Perifollicular",
    "Trabecula": "Trabeculae",
    "Trabeculae": "Trabeculae",
}


def find_geojson(sample):
    for d in GEOJSON_DIRS:
        if not d.exists():
            continue
        for path in sorted(d.glob("*")):
            if path.suffix not in (".geojson", ".qptiff") and "qptiff" not in path.name:
                continue
            if du.extract_sample_id(path.name.replace(".geojson", "")) == sample:
                return path
    return None


def load_regions(path, sample):
    """QuPath polygons in µm, keyed by canonical region name."""
    px = du.get_pixel_size(sample)
    data = json.load(open(path))
    feats = data["features"] if isinstance(data, dict) else data
    out = {}
    for ft in feats:
        cls = (ft.get("properties", {}) or {}).get("classification") or {}
        name = cls.get("name") if isinstance(cls, dict) else cls
        region = CLASS_TO_REGION.get(name)
        if region is None:
            continue
        g = shapely_scale(shape(ft["geometry"]), xfact=px, yfact=px, origin=(0, 0))
        out.setdefault(region, []).extend(du._as_polygon_list(g, min_area=0.0))
    return out


def annotation_areas():
    """Per-sample QuPath region area in mm², straight from the CSVs."""
    frames = []
    for csv, colmap in ((du.DATA_ROOT / "Measurements" / "AnnotationsFinal.csv", None),
                        (du.DATA_ROOT / "Measurements" / "ForSH2B3.csv", None)):
        if not csv.exists():
            continue
        df = pd.read_csv(csv, usecols=["Image", "Classification", "Area µm^2"])
        if colmap:
            df["Classification"] = df["Classification"].map(lambda c: colmap.get(c, c))
        df["Sample"] = df["Image"].apply(du.extract_sample_id)
        frames.append(df)
    df = pd.concat(frames)
    df["Region"] = df["Classification"].map(lambda c: CLASS_TO_REGION.get(c, c))
    return (df.groupby(["Sample", "Region"])["Area µm^2"].sum() / 1e6)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out-dir", default=str(du.OUTPUT_ROOT / "analysis" / "regions"))
    ap.add_argument("--table", default=str(du.TABLES_DIR / "region_layer_qc.csv"))
    ap.add_argument("--images", nargs="*", default=None,
                    help="rebuild only these images or samples (substring match). "
                         "The QC table is then written for this subset only, so "
                         "use it for iteration, not to produce the cohort table.")
    ap.add_argument("--include-excluded", action="store_true",
                    help="also build HDL018/HDL021, which are excluded for "
                         "quality but which H12/H13 deliberately re-include. "
                         "Pair with --images and a separate --table so the "
                         "cohort QC table stays at the 14 analysis donors.")
    args = ap.parse_args()
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    ref = annotation_areas()
    wpath = du.TABLES_DIR / "radial_pfz_width.csv"
    PFZ_WIDTH, PFZ_WIDTH_DEFAULT = {}, 75.0
    if wpath.exists():
        w = pd.read_csv(wpath)
        PFZ_WIDTH_DEFAULT = float(w["PFZ_width_um"].median())
        PFZ_WIDTH = {r.Sample: float(r.PFZ_width_um) for r in w.itertuples()
                     if r.PFZ_width_um == r.PFZ_width_um and r.PFZ_width_um <= 200}
        print(f"perifollicular width: {len(PFZ_WIDTH)} donors measured, "
              f"default {PFZ_WIDTH_DEFAULT:.0f} µm", flush=True)
    print("loading cells (tissue mask + perifollicular markers) …", flush=True)
    PFZ_MARKERS = [*du.PFZ_MARKERS, *du.PFZ_VASCULAR_MARKERS]
    cells = du.load_single_cells(markers=PFZ_MARKERS,
                                 include_excluded=args.include_excluded)
    print(f"  {len(cells):,} cells, {cells['Sample'].nunique()} donors", flush=True)

    rows = []
    for image, g in cells.groupby("Image", observed=True):
        sample = g["Sample"].iloc[0]
        if args.images and not any(w in image or w in sample for w in args.images):
            continue
        path = find_geojson(sample)
        if path is None:
            print(f"  {sample:<14} NO GEOJSON — skipped", flush=True)
            continue
        qp = load_regions(path, sample)

        tissue = du.tissue_extent(g[["x", "y"]].values.astype(float))
        tissue_u = unary_union(tissue) if tissue else None

        regions = {}
        for name in ("Follicle", "PALS", "Trabeculae"):
            polys = qp.get(name, [])
            if tissue_u is not None and polys:
                polys = [q for p in polys
                         for q in du._as_polygon_list(p.intersection(tissue_u), 0.0)]
            if name == "Follicle":
                polys = [p for p in polys if p.area >= du.MIN_QUPATH_FOLLICLE_AREA]
            elif name == "PALS":
                polys = [p for p in polys if p.area >= du.MIN_QUPATH_PALS_AREA]
            regions[name] = polys

        # Perifollicular is the ring immediately adjacent to the follicle, and
        # its width is measured rather than assumed. The radial profiles
        # (scripts/radial_follicle_profiles.py) show SMA low inside the
        # follicle, peaking at +15 µm outside it, and back to red-pulp baseline
        # by ~135 µm; the half-max-above-baseline crossing gives a per-donor
        # width of 45-105 µm, median 75. Donors with no SMA peak (HDL070,
        # HDL073) take the cohort median.
        #
        # The band is an exact geometric buffer of the follicle, not a raster,
        # so the boundary sits wherever the width places it rather than on a
        # 20 µm grid, and point-in-polygon assignment of cells and vessels is
        # exact.
        width = PFZ_WIDTH.get(sample, PFZ_WIDTH_DEFAULT)
        fol_u = unary_union(regions["Follicle"]) if regions["Follicle"] else None
        if fol_u is not None:
            band = fol_u.buffer(width).difference(fol_u)
            if tissue_u is not None:
                band = band.intersection(tissue_u)
            if regions["Trabeculae"]:
                band = band.difference(unary_union(regions["Trabeculae"]))
            regions["Perifollicular"] = du._as_polygon_list(band, 0.0)
        else:
            regions["Perifollicular"] = []
        # Where the classifier traced peripheral white pulp explicitly, use it
        # as well as the measured band — it is the same structure, drawn from
        # the full marker panel.
        if qp.get("Perifollicular"):
            pwp = qp["Perifollicular"]
            if tissue_u is not None:
                pwp = [q for p in pwp
                       for q in du._as_polygon_list(p.intersection(tissue_u), 0.0)]
            regions["Perifollicular"] = du._as_polygon_list(
                unary_union(regions["Perifollicular"] + pwp), 0.0)
        if regions["Perifollicular"] and regions["Follicle"]:
            regions["Perifollicular"] = [
                q for p in regions["Perifollicular"]
                for q in du._as_polygon_list(p.difference(fol_u), 0.0)]
        row_width = width

        pals_raw_mm2 = (unary_union(regions["PALS"]).area / 1e6
                        if regions["PALS"] else 0.0)
        # PALS keeps its territory where the two meet. Letting the band win
        # erased PALS for donors whose sheath hugs the follicle — HDL073 went to
        # 0.069 mm2 of PALS from 0.20. The band is the part of the ring that is
        # not already classified sheath.
        if regions["Perifollicular"] and regions["PALS"]:
            pals_u = unary_union(regions["PALS"])
            regions["Perifollicular"] = [
                q for p in regions["Perifollicular"]
                for q in du._as_polygon_list(p.difference(pals_u), 0.0)]
        # A perifollicular piece with no follicle to be peri- to is an artifact
        # of differencing, not a region.
        if regions["Perifollicular"] and fol_u is not None:
            near = fol_u.buffer(1.0)
            regions["Perifollicular"] = [p for p in regions["Perifollicular"]
                                         if p.intersects(near)]

        claimed = [p for k in ("Follicle", "PALS", "Trabeculae", "Perifollicular")
                   for p in regions[k]]
        regions["RedPulp"] = (
            du._as_polygon_list(tissue_u.difference(unary_union(claimed)))
            if tissue_u is not None and claimed else (tissue or []))

        regions = {k: [p for p in v if p.area >= du.SLIVER_AREA]
                   for k, v in regions.items()}

        stem = image.replace(".ome.tiff", "").replace(".ome.tif", "").replace(".tif", "")
        feats = [{"type": "Feature",
                  "properties": {"region": name, "sample": sample, "image": image},
                  "geometry": mapping(poly)}
                 for name, polys in regions.items() for poly in polys]
        (out_dir / f"{stem}.geojson").write_text(
            json.dumps({"type": "FeatureCollection", "features": feats}))

        row = {"Image": image, "Sample": sample,
               "Genotype": str(g["Genotype"].iloc[0]),
               "Platform": g["Platform"].iloc[0], "n_cells": len(g),
               "tissue_mm2": (tissue_u.area / 1e6) if tissue_u is not None else 0.0,
               "geojson": path.name,
               "trabeculae_source": "qupath" if regions["Trabeculae"] else "missing",
               "flagged": "", "pfz_width_um": row_width}
        checks = []
        for name in du.REGION_ORDER:
            area = unary_union(regions[name]).area / 1e6 if regions[name] else 0.0
            row[f"{name}_mm2"] = area
            row[f"{name}_n"] = len(regions[name])
            if name == "Perifollicular":
                row["Perifollicular_follicles_with_band"] = (
                    sum(1 for f in regions["Follicle"]
                        if any(b.intersects(f.buffer(1.0))
                               for b in regions["Perifollicular"])))
            if name in ("Follicle", "PALS"):
                exp = ref.get((sample, name), np.nan)
                if name == "PALS":
                    area = pals_raw_mm2      # check before the PFZ trim
                row[f"{name}_csv_mm2"] = exp
                row[f"{name}_ratio"] = area / exp if exp and exp == exp else np.nan
                checks.append(f"{name} {area:.2f}/{exp:.2f}={row[f'{name}_ratio']:.2f}x")
        r_f = row.get("Follicle_ratio", np.nan)
        if r_f == r_f and not (0.9 <= r_f <= 1.1):
            row["flagged"] = "Follicle"
        rows.append(row)
        print(f"  {sample:<14} tissue {row['tissue_mm2']:7.2f}  "
              f"Fol {row['Follicle_mm2']:7.3f} ({row['Follicle_n']:4d})  "
              f"PALS {row['PALS_mm2']:6.3f} ({row['PALS_n']:4d})  "
              f"RP {row['RedPulp_mm2']:8.3f}   [{'; '.join(checks)}]", flush=True)

    qc = pd.DataFrame(rows).sort_values("Sample")
    qc.to_csv(args.table, index=False)
    print(f"\nwrote {len(qc)} images → {out_dir}\nwrote QC table → {args.table}")

    bad = qc[(qc.Follicle_ratio.notna()) &
             ((qc.Follicle_ratio < 0.9) | (qc.Follicle_ratio > 1.1))]
    if len(bad):
        print("\nGeometry does NOT reproduce the annotation CSV areas:")
        print(bad[["Sample", "Follicle_mm2", "Follicle_csv_mm2",
                   "Follicle_ratio"]].to_string(index=False))
    else:
        print("\nall donors reproduce the CSV follicle area within 10%")
    miss = qc[qc.trabeculae_source == "missing"]["Sample"].tolist()
    if miss:
        print(f"\nno Trabeculae geometry for {len(miss)} donors "
              f"(area sits in RedPulp): {', '.join(miss)}")


if __name__ == "__main__":
    main()
