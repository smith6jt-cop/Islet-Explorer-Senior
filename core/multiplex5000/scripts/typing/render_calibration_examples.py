#!/usr/bin/env python3
"""Render illustrative raw-image crops for matched REDSEA/all-raw call changes.

The same source pixels appear for both arms: only the marker measurement/call policy changes.
These crops are examples for review, not independent truth labels.
"""
from __future__ import annotations

import argparse
import json
import mmap
import re
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pyarrow as pa
import pyarrow.compute as pc
import pyarrow.parquet as pq
import pyvips
import tifffile
from matplotlib.patches import Polygon


REPO = Path(__file__).resolve().parents[4]
DEFAULT_REVIEW = REPO / "results/v3/release_20261003/calibration_review"
EXAMPLE_TRANSITIONS = [
    ("Immune", "Exocrine"),
    ("Exocrine", "Immune"),
    ("Blood vessel", "Exocrine"),
]


def cell_feature(geojson: Path, object_id: str) -> dict:
    needle = b'"id":"' + object_id.encode() + b'"'
    with geojson.open("rb") as fh:
        mm = mmap.mmap(fh.fileno(), 0, access=mmap.ACCESS_READ)
        at = mm.find(needle)
        if at < 0:
            mm.close()
            raise KeyError(f"{object_id} not found in {geojson}")
        start = mm.rfind(b'{"type":"Feature","id":', 0, at + 1)
        if start < 0:
            mm.close()
            raise ValueError(f"Could not find feature boundary for {object_id}")
        depth, in_string, escaped, end = 0, False, False, -1
        for pos in range(start, len(mm)):
            c = mm[pos]
            if in_string:
                if escaped:
                    escaped = False
                elif c == 92:
                    escaped = True
                elif c == 34:
                    in_string = False
            elif c == 34:
                in_string = True
            elif c == 123:
                depth += 1
            elif c == 125:
                depth -= 1
                if depth == 0:
                    end = pos + 1
                    break
        if end < 0:
            mm.close()
            raise ValueError(f"Could not parse feature boundary for {object_id}")
        feature = json.loads(mm[start:end])
        mm.close()
    return feature


def polygon_xy(feature: dict) -> np.ndarray:
    geom = feature["geometry"]
    if geom["type"] == "Polygon":
        ring = geom["coordinates"][0]
    elif geom["type"] == "MultiPolygon":
        ring = max(geom["coordinates"], key=lambda item: len(item[0]))[0]
    else:
        raise ValueError(f"Expected cell Polygon geometry, got {geom['type']}")
    return np.asarray(ring, dtype=float)


def qptiff_page_map(path: Path) -> dict[str, int]:
    mapping = {}
    with tifffile.TiffFile(path) as tf:
        for index, page in enumerate(tf.pages[:70]):
            desc = page.description or ""
            match = re.search(r"<Biomarker>(.*?)</Biomarker>", desc, flags=re.S)
            if match:
                name = match.group(1).strip()
                mapping.setdefault(name, index)
    return mapping


def crop_channel(path: Path, page: int, x: int, y: int, half: int) -> np.ndarray:
    image = pyvips.Image.new_from_file(str(path), page=page)
    x0 = max(0, min(image.width - 1, x - half))
    y0 = max(0, min(image.height - 1, y - half))
    width = min(2 * half, image.width - x0)
    height = min(2 * half, image.height - y0)
    arr = np.ndarray(buffer=image.crop(x0, y0, width, height).write_to_memory(),
                     dtype=np.uint16, shape=(height, width))
    return arr.copy()


def stretch(arr: np.ndarray) -> np.ndarray:
    valid = arr[np.isfinite(arr)]
    hi = float(np.quantile(valid, .997)) if valid.size else 1.0
    if hi <= 0:
        hi = 1.0
    return np.clip(arr.astype(np.float32) / hi, 0, 1) ** 0.65


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--review-dir", type=Path, default=DEFAULT_REVIEW)
    ap.add_argument("--donor", default="115")
    ap.add_argument("--half-width-px", type=int, default=95)
    args = ap.parse_args()
    root = args.review_dir.resolve()
    examples = pd.read_csv(root / "candidate_image_examples.csv", dtype={"donor_id": str})
    candidates = examples[(examples.donor_id == args.donor) &
                          (examples.arm == "raw_shifted_frozen")]
    chosen = []
    for old, new in EXAMPLE_TRANSITIONS:
        subset = candidates[(candidates.from_cell_class == old) &
                            (candidates.to_cell_class == new)]
        if not subset.empty:
            chosen.append(subset.sort_values(["transition_n", "object_id"],
                                             ascending=[False, True]).iloc[0])
    if not chosen:
        raise RuntimeError(f"No changed-cell examples found for donor {args.donor}")
    donor_record = next(r for r in json.loads((root.parent / "current_input_manifest.json").read_text())["records"]
                        if r.get("donor_id") == args.donor and r.get("tissue") == "pancreas"
                        and r.get("modality") == "protein_if" and r.get("source_path"))
    qptiff = Path(donor_record["source_path"])
    if not qptiff.is_file():
        raise FileNotFoundError(qptiff)
    stem = qptiff.name + "_-_resolution_1.geojson"
    geojson = Path("/home/smith6jt/IO60panc2nd/results/v3/export/cell_geojson") / stem
    if not geojson.is_file():
        matches = list(Path("/home/smith6jt/IO60panc2nd/results/v3/export/cell_geojson").glob(f"*{qptiff.name}*resolution_1.geojson"))
        if len(matches) != 1:
            raise FileNotFoundError(f"could not resolve one cell GeoJSON for {qptiff.name}: {matches}")
        geojson = matches[0]
    pages = qptiff_page_map(qptiff)
    run_red = json.loads(Path("/data/io60-calibration/arms/redsea_shifted_frozen/run.json").read_text())
    run_raw = json.loads(Path("/data/io60-calibration/arms/raw_shifted_frozen/run.json").read_text())
    red_path = Path(run_red["aligned_root"]) / f"donor_id={args.donor}/data.parquet"
    raw_path = Path(run_raw["aligned_root"]) / f"donor_id={args.donor}/data.parquet"
    ids = [str(row.object_id) for row in chosen]
    value_rows = {}
    for policy, path in [("REDSEA", red_path), ("all raw", raw_path)]:
        pf = pq.ParquetFile(path)
        marker_names = sorted({str(row.top_expressed_marker) for row in chosen})
        cols = ["object_id"] + [m for m in marker_names if m in pf.schema_arrow.names]
        table = pq.read_table(path, columns=cols)
        matched = table.filter(pc.is_in(table["object_id"], value_set=pa.array(ids)))
        for row in matched.to_pylist():
            oid = row.pop("object_id")
            value_rows.setdefault(oid, {})[policy] = row

    fig, axes = plt.subplots(1, len(chosen), figsize=(5.2 * len(chosen), 5.8), squeeze=False,
                             constrained_layout=True)
    example_records = []
    for ax, row in zip(axes[0], chosen):
        oid, marker = str(row.object_id), str(row.top_expressed_marker)
        if marker not in pages:
            raise KeyError(f"{marker} not found among QPTIFF biomarker pages")
        feature = cell_feature(geojson, oid)
        poly = polygon_xy(feature)
        # The typing table stores physical micrometres; the registered QuPath cell
        # GeoJSON stores full-resolution image pixels. Use the exact-ID polygon for
        # crop placement and keep the centroid conversion discrepancy as a QA field.
        spacing_x = float(donor_record["pixel_size_um_x"])
        spacing_y = float(donor_record["pixel_size_um_y"])
        median_xy = np.median(poly, axis=0)
        x, y = int(round(median_xy[0])), int(round(median_xy[1]))
        centroid_delta_px = float(np.linalg.norm(
            median_xy - np.array([float(row.X_centroid) / spacing_x,
                                  float(row.Y_centroid) / spacing_y])))
        half = args.half_width_px
        dapi = crop_channel(qptiff, pages["DAPI"], x, y, half)
        marker_img = crop_channel(qptiff, pages[marker], x, y, half)
        # DAPI blue, selected marker magenta; colors are display-only.
        blue, red = stretch(dapi), stretch(marker_img)
        rgb = np.stack([red, red * .14, blue], axis=-1)
        ax.imshow(rgb, origin="upper", interpolation="nearest")
        ax.add_patch(Polygon(poly - np.array([x - half, y - half]), closed=True,
                             fill=False, edgecolor="#ffe86b", linewidth=1.1))
        ax.scatter([half], [half], s=18, c="#ffffff", edgecolors="#111111", linewidths=.5)
        red_val = value_rows.get(oid, {}).get("REDSEA", {}).get(marker, np.nan)
        raw_val = value_rows.get(oid, {}).get("all raw", {}).get(marker, np.nan)
        title = (f"{row.from_cell_class} → {row.to_cell_class}\n"
                 f"{marker}: REDSEA {red_val:.3g} · raw {raw_val:.3g}\n"
                 f"transition n={int(row.transition_n):,}; cell {oid[:8]}")
        ax.set_title(title, fontsize=10, pad=9)
        ax.set_axis_off()
        example_records.append({"donor_id": args.donor, "object_id": oid,
            "redsea_cell_class": row.from_cell_class, "all_raw_cell_class": row.to_cell_class,
            "transition_n_donor": int(row.transition_n), "marker": marker,
            "redsea_marker_value": red_val, "all_raw_marker_value": raw_val,
            "X_um": float(row.X_centroid), "Y_um": float(row.Y_centroid),
            "X_pixel": x, "Y_pixel": y, "pixel_size_um_x": spacing_x,
            "pixel_size_um_y": spacing_y, "centroid_vs_polygon_delta_px": centroid_delta_px,
            "qupath_image": str(qptiff),
            "cell_geojson": str(geojson), "crop_file": "calibration_image_examples.png"})
    fig.suptitle("Candidate REDSEA vs all-raw call changes — same source image in both arms\n"
                 "Raw IF intensity display with DAPI (blue), selected marker (magenta), cell outline (yellow)",
                 fontsize=12)
    fig.savefig(root / "calibration_image_examples.png", dpi=220, facecolor="white")
    plt.close(fig)
    pd.DataFrame(example_records).to_csv(root / "calibration_image_examples.csv", index=False)
    print(f"wrote {root / 'calibration_image_examples.png'} with {len(example_records)} candidates")


if __name__ == "__main__":
    main()
