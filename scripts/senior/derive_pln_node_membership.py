#!/usr/bin/env python
"""Derive pLN node membership geometrically, from the QuPath annotations.

The node is the pLN's structural unit. Which annotation delimits it is recorded by hand in
core/phenocycler/data/integration/pln_node_annotations.csv and is never inferred here.

Membership is a point-in-polygon test of `Centroid X µm` / `Centroid Y µm` against the node
annotation ROI. It is never taken from the TSV's `Parent` column, which reads
"Annotation (Tissue)" on every row, and never from centroid matching between modalities.

Inputs
  core/phenocycler/data/integration/pln_node_annotations.csv   which annotation is the node
  /data/IO60pLN/results/node_geojson/entry<N>.geojson          annotations, exported from the
      QuPath .qpdata by core/phenocycler/scripts/groovy/export_pln_node_annotations.groovy
  /data/IO60pLN/data/<N>/server.json                           pixel size, to put the pixel-space
      ROI into the micron space the centroids are recorded in
  a 5-column cut of /data/IO60pLN/results/IO60pLN_cells.tsv    Image, Object ID, Parent,
      Centroid X µm, Centroid Y µm. The cut is not kept in the checkout; remake it in ~2 min.
      The micron columns are spelled with U+00B5 MICRO SIGN, not "um" and not U+03BC GREEK
      SMALL LETTER MU -- quote them exactly as written here or duckdb raises BinderException:

          COPY (SELECT "Image", "Object ID", "Parent", "Centroid X µm", "Centroid Y µm"
                FROM read_csv('/data/IO60pLN/results/IO60pLN_cells.tsv', header=true,
                              delim='\t', all_varchar=true))
          TO '/data/tmp/duckdb_pln/pln_centroids.tsv' (FORMAT CSV, DELIMITER '\t', HEADER)

      (10,192,986 rows, measured 2026-09-17; the SELECT above was run against the export on
      2026-09-17 and returns rows -- it is pasteable as written)

Outputs (under data/development/pln/)
  node_membership/pln_node_membership_<donor>.parquet
      object_id, donor_id, image_file, node_id, x_um, y_um
  node_membership_summary.csv       one row per node: cells, area, how it is delimited
  node_geometry_audit.csv           one row per annotation, including the ones refused

The parquet files are NOT hive-partitioned and the CSVs sit outside their directory, both on
purpose. `donor_id` is a string column inside each file, so a hive key of the same name would
be inferred as int32 and collide -- the failure docs/ISLET_LAYERS_SPEC.md records for expression,
cells_xen, cells_pheno and typing/immune. Read the tree with
`pyarrow.dataset("data/development/pln/node_membership")`: the default partitioning is None and
donor_id comes back as the string it was written as.

A node whose annotation carries no area cannot deliver membership. It is reported as
unavailable with the measured reason, never as zero cells.

This is the measurement, and it is a development lane. The published layer artifact --
$ISLET_LAYERS_DIR (default /data/islet_layers) then /pln, with its MANIFEST and the QuPath
cross-check -- is built from this output by core/phenocycler/scripts/pln/build_pln_layer.py.
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pyarrow as pa
import pyarrow.parquet as pq
from shapely.geometry import shape
from shapely.strtree import STRtree

REPO = Path(__file__).resolve().parents[2]
ANNOT_CSV = REPO / "core/phenocycler/data/integration/pln_node_annotations.csv"
QP_PROJECT = Path("/data/IO60pLN")
GEOJSON_DIR = QP_PROJECT / "results/node_geojson"
OUT_ROOT = REPO / "data/development/pln"

# 6380 carries no pLN_* region; the maintainer's ruling is that its single Tissue annotation IS
# the node. The ruling is keyed on the donor so it cannot silently widen to another slide.
TISSUE_IS_NODE_DONOR = "6380"


def read_annot_csv() -> pd.DataFrame:
    df = pd.read_csv(ANNOT_CSV, comment="#", dtype={"donor_id": str, "qupath_entry": int})
    df["node_classes"] = df["node_classes"].str.split(";")
    df["excluded"] = df["excluded"].astype(str).str.lower() == "true"
    return df


def pixel_size_um(entry: int) -> tuple[float, float]:
    meta = json.loads((QP_PROJECT / f"data/{entry}/server.json").read_text())
    cal = meta["metadata"]["pixelCalibration"]
    return float(cal["pixelWidth"]["value"]), float(cal["pixelHeight"]["value"])


def is_axis_aligned_rectangle(geom) -> bool:
    """True when the ring is exactly its own bounding box.

    It matters. A rectangle is a separator drawn to say which piece of tissue is which node --
    it is not an outline of the node, and its area is not the node's area.
    """
    if geom.geom_type != "Polygon" or len(geom.interiors) > 0:
        return False
    coords = list(geom.exterior.coords)
    if len(coords) != 5:
        return False
    xs = {round(x, 6) for x, _ in coords}
    ys = {round(y, 6) for _, y in coords}
    return len(xs) == 2 and len(ys) == 2


def load_annotations(entry: int) -> list[dict]:
    """Every annotation on one image, as {classification, geom_px, area_px, is_area}."""
    gj = json.loads((GEOJSON_DIR / f"entry{entry}.geojson").read_text())
    feats = gj["features"] if gj.get("type") == "FeatureCollection" else gj
    out = []
    for f in feats:
        props = f.get("properties") or {}
        cls = (props.get("classification") or {}).get("name")
        geom = shape(f["geometry"])
        out.append(
            {
                "object_id": f.get("id"),
                "classification": cls,
                "geom_type": f["geometry"]["type"],
                "geom_px": geom,
                "area_px": float(geom.area),
                "is_area": geom.geom_type in ("Polygon", "MultiPolygon") and geom.area > 0,
                "is_rectangle": is_axis_aligned_rectangle(geom),
                "n_vertices": len(geom.exterior.coords) if geom.geom_type == "Polygon" else None,
            }
        )
    return out


def node_annotations(row, annots: list[dict]) -> tuple[list[dict], list[dict]]:
    """Split the image's annotations into the ones that delimit a node and the ones refused."""
    recorded = list(row["node_classes"])
    keep, refused = [], []
    if row["donor_id"] == TISSUE_IS_NODE_DONOR:
        if recorded != ["Tissue"]:
            raise SystemExit(f"{row['donor_id']}: expected node_classes Tissue, got {recorded}")
        cand = [a for a in annots if a["classification"] == "Tissue"]
        if len(cand) != 1:
            raise SystemExit(f"{row['donor_id']}: Tissue-is-node ruling needs exactly 1 Tissue, "
                             f"found {len(cand)}")
    else:
        cand = [a for a in annots if a["classification"] in recorded]
    found = sorted({a["classification"] for a in cand})
    if found != sorted(recorded):
        raise SystemExit(f"{row['donor_id']}: recorded node_classes {sorted(recorded)} but the "
                         f"annotations carry {found}")
    for a in cand:
        (keep if a["is_area"] else refused).append(a)
    return keep, refused


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--centroids", required=True,
                    help="TSV cut of IO60pLN_cells.tsv: Image, Object ID, Parent, "
                         "Centroid X \u00b5m, Centroid Y \u00b5m (U+00B5, as the export "
                         "spells them)")
    ap.add_argument("--out", default=str(OUT_ROOT))
    args = ap.parse_args()
    out_root = Path(args.out)
    out_root.mkdir(parents=True, exist_ok=True)

    roster = read_annot_csv()
    eligible = roster[~roster["excluded"]].copy()
    print(f"roster: {len(roster)} images, {len(eligible)} eligible, "
          f"{int(eligible['n_nodes'].sum())} nodes recorded", flush=True)

    # geometry first, so a refusal happens before the 9M-row read
    per_image: dict[str, dict] = {}
    audit_rows: list[dict] = []
    for _, row in eligible.iterrows():
        entry = int(row["qupath_entry"])
        annots = load_annotations(entry)
        px_w, px_h = pixel_size_um(entry)
        keep, refused = node_annotations(row, annots)
        tissue = [a for a in annots if a["classification"] == "Tissue" and a["is_area"]]
        for a in annots:
            is_node = a in keep
            audit_rows.append(
                {
                    "donor_id": row["donor_id"], "qupath_entry": entry,
                    "image_file": row["image_file"], "annotation_id": a["object_id"],
                    "classification": a["classification"], "geom_type": a["geom_type"],
                    "area_px": a["area_px"],
                    "area_um2": a["area_px"] * px_w * px_h,
                    "is_axis_aligned_rectangle": a["is_rectangle"],
                    "n_vertices": a["n_vertices"],
                    "role": "node" if is_node else (
                        "refused_no_area" if a in refused else "not_a_node"),
                }
            )
        per_image[row["image_file"]] = {
            "row": row, "entry": entry, "px_w": px_w, "px_h": px_h,
            "nodes": keep, "refused": refused, "tissue": tissue,
        }

    (out_root / "node_membership").mkdir(parents=True, exist_ok=True)
    pd.DataFrame(audit_rows).to_csv(out_root / "node_geometry_audit.csv", index=False)

    # centroids, one image at a time
    print(f"reading {args.centroids}", flush=True)
    cent = pd.read_csv(args.centroids, sep="\t", dtype=str)
    cent.columns = [c.strip() for c in cent.columns]
    xcol = [c for c in cent.columns if c.startswith("Centroid X")][0]
    ycol = [c for c in cent.columns if c.startswith("Centroid Y")][0]
    cent["image_file"] = cent["Image"].str.split(" - resolution", regex=False).str[0]
    print(f"centroids: {len(cent):,} rows, {cent['image_file'].nunique()} images", flush=True)

    summary_rows: list[dict] = []
    for image_file, info in per_image.items():
        sub = cent[cent["image_file"] == image_file]
        row, px_w, px_h = info["row"], info["px_w"], info["px_h"]
        donor = row["donor_id"]
        if len(sub) == 0:
            print(f"  {donor} {image_file}: NO ROWS in the centroid table", flush=True)
            continue
        x_um = sub[xcol].astype("float64").to_numpy()
        y_um = sub[ycol].astype("float64").to_numpy()
        # the ROI is in pixels; put it in the micron space the centroids are recorded in
        from shapely import points as sh_points
        P = sh_points(np.column_stack([x_um / px_w, y_um / px_h]))

        def assign(geoms: list, labels: list[str]) -> tuple[np.ndarray, np.ndarray]:
            """label index per point, and how many geometries contain it."""
            idx = np.full(len(P), -1, dtype=np.int64)
            hits = np.zeros(len(P), dtype=np.int64)
            tree = STRtree(geoms)
            pairs = tree.query(P, predicate="within")  # (2, n) point-idx, geom-idx
            for p_i, g_i in zip(pairs[0], pairs[1]):
                hits[p_i] += 1
                if idx[p_i] == -1:
                    idx[p_i] = g_i
            return idx, hits

        node_geoms = [a["geom_px"] for a in info["nodes"]]
        node_labels = [a["classification"] for a in info["nodes"]]
        node_idx, node_hits = assign(node_geoms, node_labels)

        tis_geoms = [a["geom_px"] for a in info["tissue"]]
        tis_idx, tis_hits = assign(tis_geoms, [a["object_id"] for a in info["tissue"]])

        node_id = np.array(["unassigned"] * len(sub), dtype=object)
        ok = node_idx >= 0
        node_id[ok] = [node_labels[i] for i in node_idx[ok]]
        ambiguous = int((node_hits > 1).sum())

        out_file = out_root / "node_membership" / f"pln_node_membership_{donor}.parquet"
        tbl = pa.table(
            {
                "object_id": pa.array(sub["Object ID"].to_numpy(), type=pa.string()),
                "donor_id": pa.array([donor] * len(sub), type=pa.string()),
                "image_file": pa.array([image_file] * len(sub), type=pa.string()),
                "node_id": pa.array(node_id.astype(str), type=pa.string()),
                "x_um": pa.array(x_um, type=pa.float64()),
                "y_um": pa.array(y_um, type=pa.float64()),
                "n_node_regions_containing": pa.array(node_hits, type=pa.int64()),
                "tissue_annotation_index": pa.array(tis_idx, type=pa.int64()),
            }
        )
        pq.write_table(tbl, out_file, compression="zstd")

        counts = pd.Series(node_id.astype(str)).value_counts()
        for label in node_labels:
            summary_rows.append(
                {
                    "donor_id": donor, "image_file": image_file,
                    "qupath_entry": info["entry"], "node_id": label,
                    "annotation_geom": next(a["geom_type"] for a in info["nodes"]
                                            if a["classification"] == label),
                    "delimiter_kind": (
                        "rectangle" if next(a["is_rectangle"] for a in info["nodes"]
                                            if a["classification"] == label) else "outline"),
                    "area_um2": next(a["area_px"] for a in info["nodes"]
                                     if a["classification"] == label) * px_w * px_h,
                    "n_cells": int(counts.get(label, 0)),
                    "status": "derived",
                }
            )
        for a in info["refused"]:
            summary_rows.append(
                {
                    "donor_id": donor, "image_file": image_file,
                    "qupath_entry": info["entry"], "node_id": a["classification"],
                    "annotation_geom": a["geom_type"], "delimiter_kind": "none",
                    "area_um2": np.nan,
                    "n_cells": pd.NA,
                    "status": f"unavailable: annotation is {a['geom_type']}, zero area -- "
                              f"it delimits no region, so no cell can be tested against it",
                }
            )
        print(f"  {donor} {image_file}: {len(sub):,} cells, "
              f"{int(counts.drop('unassigned', errors='ignore').sum()):,} in a node, "
              f"{int(counts.get('unassigned', 0)):,} unassigned, {ambiguous:,} in >1 node, "
              f"{len(info['refused'])} node(s) refused", flush=True)

    summary = pd.DataFrame(summary_rows)
    summary.to_csv(out_root / "node_membership_summary.csv", index=False)
    print(f"\nwrote {out_root}/node_membership_summary.csv")
    return 0


if __name__ == "__main__":
    sys.exit(main())
