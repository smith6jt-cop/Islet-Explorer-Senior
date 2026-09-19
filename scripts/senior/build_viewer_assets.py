#!/usr/bin/env python3
"""LEGACY LANE, NOT RUNNABLE HERE -- recorded 2026-09-18.

``scripts/retire/dispositions.py`` records this file as a verified orphan that is also the
NAMED PRODUCER of a live artifact under ``data/app_data`` (code tier 2), which is why it is
kept rather than deleted. It cannot be re-run: the inputs named below are absent from this
checkout.

The current PhenoCycler typing is the subpanel classification under
``core/phenocycler/scripts/typing/`` -- context from QuPath objects, then one question per
context against reference levels drawn from the ordinary tissue, fitted PER DONOR.
``data/development/typing/cells/`` holds 25,230,109 cells over 20 donors (2026-09-15 build,
eight classes plus ``Excluded``). Nothing below describes that method, and ``data/app_data``
is the OTHER app domain: the two datasets share five donors and no vocabulary.

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

SD = Path(__file__).resolve().parents[2] / "data"
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
# Used ONLY by build_structure_geojsons() (--structures), which writes a separate
# file for the TissUUmaps viewer. Still absent from CLASS_MAP on purpose: the sf
# deferral above is about the Shiny app's drill-down and has not changed.
STRUCTURE_MAP = {"Blood_vessel": "Capillary", "Nerve": "Nerve", "Lymphatic": "Lymphatic"}

def discover_donors():
    """Auto-discover the cohort from the tissue partitions (case_id=<donor>), so adding
    donors needs no code edit. Falls back to the cells partitions if tissue isn't built."""
    for sub in ("parquet/tissue", "cells"):
        base = SD / sub
        if base.exists():
            key = "case_id" if sub.endswith("tissue") else "donor_id"
            ds = sorted(p.name.split("=")[1] for p in base.glob(f"{key}=*"))
            if ds:
                return ds
    return []


DONORS = discover_donors()


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
    _write_class_geojsons(CLASS_MAP, JSON_OUT, "{donor}.geojson", "geojson")


def build_structure_geojsons():
    """Vessel/nerve/lymphatic outlines for the TissUUmaps viewer ONLY.

    Deliberately a SEPARATE file from the per-donor `<donor>.geojson` above, and
    STRUCTURE_MAP is deliberately still absent from CLASS_MAP. Merging these into
    the app's own GeoJSON is what the v1 deferral at the top of this file rules
    out: ~200k features per donor makes `sf::st_read` unusably slow on the first
    drill-down click, and that is still true. The viewer reads its own file, as
    its own region layer, which it does not auto-load.

    Written to a `structures/` subdirectory so a case-id search can never match
    both this and the islet file for the same donor -- the builder's region
    lookup refuses to guess between two matches and would attach neither.
    """
    out_dir = JSON_OUT / "structures"
    _write_class_geojsons(STRUCTURE_MAP, out_dir, "{donor}_structures.geojson", "structures")


def _write_class_geojsons(class_map, out_dir, name_tmpl, tag):
    """Merge each donor's per-class QuPath exports into one FeatureCollection.

    Properties are stripped to name + classification on purpose: QuPath's
    measurement dictionaries travel with every feature and dominate the file
    size, and `regionUtils.geoJSON2regions()` keeps whatever it is given on each
    region in memory.
    """
    out_dir.mkdir(parents=True, exist_ok=True)
    for donor in DONORS:
        feats = []
        counts = {}
        for cls_prefix, cls_name in class_map.items():
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
            print(f"[{tag}] {donor}: NO source files found -- skipped")
            continue
        fc = {"type": "FeatureCollection", "features": feats}
        dst = out_dir / name_tmpl.format(donor=donor)
        with open(dst, "w") as fh:
            json.dump(fc, fh)
        print(f"[{tag}] {donor}: {len(feats)} features {counts} "
              f"({dst.stat().st_size / 1e6:.1f} MB)")


if __name__ == "__main__":
    import argparse

    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--structures", action="store_true",
                    help="Also write structures/<donor>_structures.geojson for the "
                         "TissUUmaps viewer. Off by default: these are large and "
                         "nothing in the Shiny app reads them.")
    a = ap.parse_args()

    build_lookup()
    build_geojsons()
    if a.structures:
        build_structure_geojsons()
    print("[done] viewer assets in", APPDATA)
