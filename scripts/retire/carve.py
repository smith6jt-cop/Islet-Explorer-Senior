#!/usr/bin/env python3
"""Rescue the columns and partitions that Tier 2's replacements do not carry.

Run BEFORE retire.py --tier2. Each carve-out is written under
data/ARCHIVE_2026-09-15/carve/ with a sidecar recording what it rescued and why.

Measured, not assumed (2026-09-16):
  data/cells, data/cells_preqc  14 columns absent from runs/.../cells
                                cell_circularity cell_length cell_maxdiam cell_mindiam
                                classification dist_islet dist_tissue nc_area_ratio
                                nucleus_area nucleus_circularity nucleus_length
                                nucleus_maxdiam nucleus_mindiam nucleus_solidity
  data/cells                    donor_id=6457 and 6579 -- in NO replacement (the run,
                                typing/cells, qupath_fundamentals and kintsugi/distance
                                are all 20-donor)
  data/cell_qc                  9 columns, + the retired gate's own parameters
  subpanel_context              in_Nerve / in_Lymphatic -- reconstruction from cells/
                                loses 18.4% of in_Nerve
  data/cells_redsea             NOTHING -- its only loss, cell_area_px, is present in
                                runs/.../redsea. No carve-out; the manifest records the
                                wider superseded_by instead.
"""
from __future__ import annotations

import argparse
import glob
import json
import os
import shutil
from datetime import datetime, timezone
from pathlib import Path

import pyarrow.parquet as pq

REPO = Path("/home/smith6jt/IO60panc2nd/Islet-Explorer-Senior")
CARVE = REPO / "data" / "ARCHIVE_2026-09-15" / "carve"

MORPH = ["cell_circularity", "cell_length", "cell_maxdiam", "cell_mindiam", "classification",
         "dist_islet", "dist_tissue", "nc_area_ratio", "nucleus_area", "nucleus_circularity",
         "nucleus_length", "nucleus_maxdiam", "nucleus_mindiam", "nucleus_solidity"]
QC = ["cell_area", "cell_area_px", "cell_solidity", "nc_area_ratio", "nucleus_area",
      "qc_estimation_eligible", "qc_exclude_reason", "qc_no_cytoplasm", "qc_pass"]


def carve_columns(src: Path, out: Path, keep: list[str], why: str) -> dict:
    """Copy object_id + `keep` from every donor partition of a hive root."""
    out.mkdir(parents=True, exist_ok=True)
    rows = donors = 0
    kept_cols: list[str] = []
    for part in sorted(glob.glob(str(src / "donor_id=*"))):
        d = os.path.basename(part)
        files = sorted(glob.glob(os.path.join(part, "*.parquet")))
        if not files:
            continue
        have = set(pq.ParquetFile(files[0]).schema_arrow.names)
        cols = ["object_id"] + [c for c in keep if c in have]
        kept_cols = cols
        t = pq.read_table(files[0], columns=cols)
        (out / d).mkdir(exist_ok=True)
        pq.write_table(t, out / d / "data.parquet", compression="zstd")
        rows += t.num_rows
        donors += 1
    return dict(kind="columns", source=str(src), rows=rows, donors=donors,
                columns=kept_cols, why=why)


def carve_partitions(src: Path, out: Path, names: list[str], why: str) -> dict:
    out.mkdir(parents=True, exist_ok=True)
    rows = 0
    got = []
    for n in names:
        p = src / f"donor_id={n}"
        if not p.exists():
            continue
        shutil.copytree(p, out / f"donor_id={n}", dirs_exist_ok=True)
        for f in glob.glob(str(out / f"donor_id={n}" / "*.parquet")):
            rows += pq.ParquetFile(f).metadata.num_rows
        got.append(n)
    return dict(kind="partitions", source=str(src), rows=rows, donors=len(got),
                partitions=got, why=why)


def carve_files(files: list[Path], out: Path, why: str) -> dict:
    out.mkdir(parents=True, exist_ok=True)
    got = []
    for f in files:
        if f.exists():
            shutil.copy2(f, out / f.name)
            got.append(f.name)
    return dict(kind="files", files=got, why=why)


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--execute", action="store_true", help="without this, only reports the plan")
    a = ap.parse_args()

    jobs = [
        ("cells_excluded_donors", carve_partitions, (REPO / "data/cells", ["6457", "6579"]),
         "donor_id=6457 and 6579 exist in data/cells and in NO replacement -- the run, "
         "typing/cells, qupath_fundamentals and kintsugi/distance are all 20-donor. Centrally "
         "excluded from analysis, but this is the only copy of their measurements."),
        ("cells_morphometry", carve_columns, (REPO / "data/cells", MORPH),
         "14 QuPath morphometry + signed-distance columns absent from runs/.../cells. "
         "signed distances survive for the 20 cohort donors in kintsugi/distance, but not for "
         "6457/6579 and not under these names."),
        ("cells_preqc_morphometry", carve_columns, (REPO / "data/cells_preqc", MORPH),
         "Same 14 columns on the pre-QC row set, which is a different population from data/cells."),
        ("cell_qc_verdict", carve_columns, (REPO / "data/cell_qc", QC),
         "9 columns with no like-named successor. qc_pass / qc_exclude_reason / qc_no_cytoplasm "
         "are the retired gate's verdict -- geometry_qc re-derives a different three-axis "
         "eligibility rather than reproducing it."),
        ("subpanel_context_objects", "flat", (REPO / "data/development/typing/scratch/subpanel_context",
                                              ["in_islet", "in_Blood vessel", "in_Nerve", "in_Lymphatic"]),
         "in_Nerve and in_Lymphatic are not reproducible from cells/: measured reconstruction "
         "loses 18.4% of in_Nerve and 325 cells of in_Lymphatic."),
    ]

    if not a.execute:
        print("PLAN (nothing written; pass --execute)\n")
        for name, _fn, (src, arg), why in jobs:
            print(f"  {name}\n    from {src}\n    {arg if isinstance(arg, list) else arg}\n    {why}\n")
        print(f"  data/cells_redsea -> NO CARVE (cell_area_px present in runs/.../redsea)")
        return 0

    CARVE.mkdir(parents=True, exist_ok=True)
    report = []
    for name, fn, (src, arg), why in jobs:
        out = CARVE / name
        print(f"[carve] {name} ...", flush=True)
        if fn == "flat":  # flat <donor>__geometry.parquet files, not a hive
            out.mkdir(parents=True, exist_ok=True)
            rows = 0
            n = 0
            for f in sorted(glob.glob(str(src / "*__geometry.parquet"))):
                have = set(pq.ParquetFile(f).schema_arrow.names)
                cols = ["object_id"] + [c for c in arg if c in have]
                t = pq.read_table(f, columns=cols)
                pq.write_table(t, out / os.path.basename(f), compression="zstd")
                rows += t.num_rows
                n += 1
            rec = dict(kind="columns", source=str(src), rows=rows, donors=n, columns=cols, why=why)
        elif fn is carve_partitions:
            rec = carve_partitions(src, out, arg, why)
        else:
            rec = carve_columns(src, out, arg, why)
        rec["name"] = name
        rec["bytes"] = sum(p.stat().st_size for p in out.rglob("*") if p.is_file())
        report.append(rec)
        print(f"         {rec['rows']:,} rows / {rec.get('donors', 0)} parts -> "
              f"{rec['bytes'] / 1e6:,.1f} MB")

    report.append(dict(name="cells_redsea", kind="none", why=(
        "No carve-out. Its only column absent from runs/.../expression is cell_area_px, which is "
        "present in runs/.../redsea. Recorded as a wider superseded_by instead.")))
    (CARVE / "CARVE.json").write_text(json.dumps(
        dict(created_utc=datetime.now(timezone.utc).isoformat(timespec="seconds"),
             carves=report), indent=2))
    print(f"\ntotal carved: {sum(r.get('bytes', 0) for r in report) / 1e9:.2f} GB")
    print(f"wrote {CARVE / 'CARVE.json'}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
