#!/usr/bin/env python3
"""
3D Senior: build the segmentation-viewer assets the Shiny app reads from data/:

  1. islet_spatial_lookup.csv  (case_id, islet_key, centroid_x_um, centroid_y_um,
     area_um2) -- per-islet core-cell centroid (µm) from the tissue partitions +
     core area from islet_targets. select_islet_from_key() needs this to resolve a
     clicked islet to a centroid (segmentation + trajectory clicks).

  2. json/<case_id>.geojson  -- one combined FeatureCollection per donor merging the
     per-class QuPath exports in results/: Islets (class 'Islet', keeps name Islet_N),
     Blood_vessel -> 'Capillary', Nerve -> 'Nerve', Lymphatic -> 'Lymphatic'. Coords
     stay in pixels (resolution #1; PIXEL_SIZE_UM=0.5078 in 00_globals maps µm->px).
     IsletExpanded (peri) boundaries are not exported by QuPath -> omitted (viewer
     shows core + structures only).

Writes into the app's data dir (the app_data symlink target).

Run in the `scimap` conda env.
"""
import json
import re
from pathlib import Path

import numpy as np
import pandas as pd

SD = Path.home() / "IO60panc2nd" / "senior_data"
RESULTS = Path.home() / "IO60panc2nd" / "results"
APPDATA = SD / "app_data"
JSON_OUT = APPDATA / "json"

# QuPath source class file-prefix -> viewer classification name.
# v1: Islet boundaries ONLY. The Blood_vessel/Nerve/Lymphatic exports carry
# 200k+ features per donor (~300 MB combined geojson), which makes sf::st_read
# unusably slow on first click. The structure OVERLAY in the viewer is deferred
# to v1.1 (needs spatial pre-indexing / per-islet clipping). Structural class
# densities in the Microenvironment tab are separately deferred (targets area_um2
# = NaN). To re-enable, add the structure prefixes back here.
CLASS_MAP = {"Islets": "Islet"}
STRUCTURE_MAP = {"Blood_vessel": "Capillary", "Nerve": "Nerve", "Lymphatic": "Lymphatic"}  # v1.1

DONORS = ["115", "6374", "6380", "6436", "6442", "6457", "6461", "6476",
          "6482", "6493", "6521", "6539", "6551", "6579", "6593"]


def token_to_case(token):
    """Map a filename donor token to the data case_id."""
    if token.isdigit():
        return token
    if "HDL" in token:           # HDLNpancLN pancreatic lymph node == donor 115
        return "115"
    return None


def find_class_file(cls_prefix, case_id):
    """Locate the results/ geojson for a (class, donor)."""
    for f in RESULTS.glob(f"{cls_prefix}_*.geojson"):
        m = re.search(rf"{cls_prefix}_+(.+?)_Scan", f.name)
        if m and token_to_case(m.group(1)) == case_id:
            return f
    return None


def build_lookup():
    rows = []
    tgt = pd.read_parquet(SD / "islet_targets.parquet")  # index "donor_Islet_N"
    for donor in DONORS:
        tp = SD / "parquet" / "tissue" / f"case_id={donor}" / "part-0.parquet"
        if not tp.exists():
            continue
        t = pd.read_parquet(tp, columns=["X_centroid", "Y_centroid",
                                         "cell_region", "islet_name"])
        core = t[(t["cell_region"] == "core") & (t["islet_name"] != "")]
        g = core.groupby("islet_name")[["X_centroid", "Y_centroid"]].mean()
        for islet_key, (cx, cy) in g.iterrows():
            area = tgt["core_area_um2"].get(f"{donor}_{islet_key}", np.nan)
            rows.append((donor, islet_key, round(cx, 2), round(cy, 2),
                         round(area, 2) if pd.notna(area) else np.nan))
    df = pd.DataFrame(rows, columns=["case_id", "islet_key",
                                     "centroid_x_um", "centroid_y_um", "area_um2"])
    df.to_csv(APPDATA / "islet_spatial_lookup.csv", index=False)
    print(f"[lookup] {len(df)} islets -> islet_spatial_lookup.csv "
          f"({df['case_id'].nunique()} donors)")


def build_geojsons():
    JSON_OUT.mkdir(parents=True, exist_ok=True)
    for donor in DONORS:
        feats = []
        counts = {}
        for cls_prefix, cls_name in CLASS_MAP.items():
            f = find_class_file(cls_prefix, donor)
            if f is None:
                continue
            gj = json.load(open(f))
            for ft in gj.get("features", []):
                props = ft.get("properties", {}) or {}
                out = {
                    "type": "Feature",
                    "geometry": ft["geometry"],
                    "properties": {
                        "name": props.get("name"),
                        "classification": {"name": cls_name},
                    },
                }
                feats.append(out)
            counts[cls_name] = counts.get(cls_name, 0) + len(gj.get("features", []))
        if not feats:
            print(f"[geojson] {donor}: NO source files found -- skipped")
            continue
        fc = {"type": "FeatureCollection", "features": feats}
        with open(JSON_OUT / f"{donor}.geojson", "w") as fh:
            json.dump(fc, fh)
        print(f"[geojson] {donor}: {len(feats)} features {counts}")


if __name__ == "__main__":
    build_lookup()
    build_geojsons()
    print("[done] viewer assets in", APPDATA)
