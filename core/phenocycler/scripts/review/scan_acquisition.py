"""Extract the instrument-level acquisition parameters from every qptiff.

These are properties of the SCANNER, not of the biology: exposure time, the
calibrated filter responsivity and its calibration date, the signal units, and
which corrections the acquisition software already applied. If they vary between
donors, raw counts are not on a common physical scale and any cross-donor
comparison of intensity is confounded by acquisition, not biology.
"""

from __future__ import annotations

import json
import xml.etree.ElementTree as ET
from pathlib import Path

import pandas as pd
import tifffile

REPO = Path(__file__).resolve().parents[2]
MANIFEST = REPO / "data" / "qupath_manifest.json"


def channel_records(path: str):
    rows = []
    with tifffile.TiffFile(path) as tf:
        for idx, page in enumerate(tf.series[0].pages):
            desc = page.description or ""
            if not desc.startswith("<?xml"):
                continue
            try:
                root = ET.fromstring(desc)
            except ET.ParseError:
                continue

            def text(tag):
                el = root.find(tag)
                return el.text if el is not None else None

            resp = root.find("Responsivity/Filter/Response")
            date = root.find("Responsivity/Filter/Date")
            fid = root.find("Responsivity/Filter/FilterID")
            rows.append({
                "slide": text("SlideID"),
                "channel_index": idx,
                "marker": text("Name") or text("Biomarker"),
                "exposure": float(text("ExposureTime")) if text("ExposureTime") else None,
                "signal_units": text("SignalUnits"),
                "responsivity": float(resp.text) if resp is not None else None,
                "calib_date": (date.text[:10] if date is not None and date.text else None),
                "filter_id": fid.text if fid is not None else None,
                "processing": text("ImageProcessingApplied"),
                "af_subtracted": text("AutofluorescenceSubtracted"),
                "software": text("AcquisitionSoftware"),
                "objective": text("Objective"),
            })
    return rows


def main() -> None:
    # Iterate the cohort manifest, never a directory glob: the manifest is the
    # only donor authority. 6523panc_Scan1.er.qptiff carries the internal
    # SlideID 6623panc (a confirmed acquisition typo), so anything that keys a
    # donor on the image itself misattributes a T1D donor's cells.
    manifest = json.loads(MANIFEST.read_text())
    rows = []
    for image in manifest["images"]:
        path = image["qptiff"]["path"]
        try:
            for record in channel_records(path):
                record["donor_id"] = str(image["donor_id"])
                rows.append(record)
        except Exception as exc:  # noqa: BLE001
            print(f"SKIP {Path(path).name}: {exc}")
    df = pd.DataFrame(rows)
    out = REPO / "data" / "pilots" / "acquisition_parameters.parquet"
    out.parent.mkdir(parents=True, exist_ok=True)
    df.to_parquet(out, index=False)
    print(f"{len(df)} channel records from {df.slide.nunique()} slides -> {out}\n")

    pd.set_option("display.width", 200)
    print("=== acquisition software / corrections already applied ===")
    print(df.groupby(["software", "processing", "af_subtracted", "objective",
                      "signal_units"]).slide.nunique().to_string())
    print("\n=== exposure time: distinct values across the cohort ===")
    print(df.exposure.value_counts().to_string())
    print("\n=== responsivity calibration dates (per slide) ===")
    print(df.groupby("calib_date").slide.nunique().to_string())
    print("\n=== per-marker responsivity spread across slides ===")
    g = df.dropna(subset=["responsivity"]).groupby("marker").responsivity.agg(
        ["nunique", "min", "max", "mean"])
    g["fold"] = g["max"] / g["min"]
    print(g.sort_values("fold", ascending=False).head(20).round(3).to_string())


if __name__ == "__main__":
    main()
