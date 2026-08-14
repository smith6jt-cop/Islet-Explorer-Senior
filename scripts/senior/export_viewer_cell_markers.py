#!/usr/bin/env python3
"""Export one whole-slide per-cell CSV per donor, for the TissUUmaps viewer.

``build_tissuumaps_project.py --markers`` wants a single CSV of every cell on a
slide. Every usage example in this repo names ``<case>_tissue.csv`` — and no
script has ever written it. Whole-slide per-cell data exists only as the Parquet
tissue partitions (`build_composition_and_tissue.py`) and, for the older 15-donor
set, as `data/app_data/donors/<case>.csv`. This closes that gap.

Output goes to ``data/tissuumaps_inputs/markers/<case>.csv``, NOT to
``data/app_data/cells/``. That directory is the per-ISLET drill-down namespace,
and `drilldown_helpers.R::drilldown_available()` decides the whole drill-down
feature is available by globbing it for any ``*.csv``. Dropping a whole-slide
file in there would switch the drill-down on for a deployment that has no
per-islet data at all.

    python scripts/senior/export_viewer_cell_markers.py
    python scripts/senior/export_viewer_cell_markers.py --cases 6539 6450

Writes ``coverage.json`` alongside the CSVs recording what each donor got and
from where. The batch driver reads it to decide whether to pass ``--markers``, so
a donor with no cell data is a recorded fact rather than a build failure — seven
batch-2 donors have slides but no cells.
"""

from __future__ import annotations

import argparse
import csv
import json
import math
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parents[2]

# Exactly the columns build_tissuumaps_project.py expects, in order. The first
# three are what the viewer plots; cell_region and islet_name ride along because
# TissUUmaps lets the user re-pick the grouping column in the Markers tab at zero
# build cost, and "colour by core/peri/tissue" is worth having for free.
MARKER_COLS = ["X_centroid", "Y_centroid", "phenotype", "cell_region", "islet_name"]

# Centroids are float32 µm. Full repr roughly doubles the file for precision far
# below a pixel (0.2539 µm/px), which nobody can see and every viewer must parse.
COORD_DP = 2


def log(msg: str) -> None:
    print(f"[markers] {msg}", flush=True)


def read_parquet_tissue(parquet_root: Path, case_id: str):
    """Rows from the canonical Parquet tissue partition, or None if absent.

    This is the same source `spatial_helpers.R::load_donor_tissue()` prefers, and
    it is written with exactly these columns by
    `scripts/senior/build_composition_and_tissue.py`. pandas is imported lazily so
    the legacy path and the coverage logic stay usable without it.
    """
    part = parquet_root / f"case_id={case_id}"
    if not part.is_dir() or not any(part.glob("*.parquet")):
        return None
    import pandas as pd

    frames = [pd.read_parquet(f, columns=MARKER_COLS) for f in sorted(part.glob("*.parquet"))]
    df = pd.concat(frames, ignore_index=True) if len(frames) > 1 else frames[0]
    return df.itertuples(index=False, name=None)


def read_legacy_csv(donors_dir: Path, case_id: str):
    """Rows from the older per-donor CSV, or None. Pure stdlib on purpose."""
    src = donors_dir / f"{case_id}.csv"
    if not src.is_file():
        return None

    def rows():
        with open(src, newline="") as fh:
            for rec in csv.DictReader(fh):
                yield tuple(rec.get(c, "") for c in MARKER_COLS)

    return rows()


def write_markers(rows, dst: Path) -> dict:
    """Write the viewer CSV, returning what was kept and what was dropped.

    Non-finite coordinates are dropped rather than passed through: TissUUmaps
    would place them at NaN, which renders as nothing or as a pile at the origin,
    and either way they are cells the viewer silently lies about. The count is
    reported so a data problem upstream is visible here rather than inferred from
    a marker layer that looks slightly wrong.
    """
    dst.parent.mkdir(parents=True, exist_ok=True)
    kept = dropped = 0
    x_max = y_max = float("-inf")
    tmp = dst.with_suffix(".csv.part")
    with open(tmp, "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(MARKER_COLS)
        for row in rows:
            try:
                x, y = float(row[0]), float(row[1])
            except (TypeError, ValueError):
                dropped += 1
                continue
            if not (math.isfinite(x) and math.isfinite(y)):
                dropped += 1
                continue
            x_max, y_max = max(x_max, x), max(y_max, y)
            w.writerow([round(x, COORD_DP), round(y, COORD_DP)] + list(row[2:]))
            kept += 1
    tmp.replace(dst)
    return {
        "n_cells": kept,
        "n_dropped": dropped,
        "max_x_um": None if x_max == float("-inf") else round(x_max, COORD_DP),
        "max_y_um": None if y_max == float("-inf") else round(y_max, COORD_DP),
        "bytes": dst.stat().st_size,
    }


def discover_cases(parquet_root: Path, donors_dir: Path) -> list[str]:
    """Every donor with cell data, from either source. Discovered, not listed —
    the same reasoning as `build_viewer_assets.py::discover_donors()`: adding a
    donor should not need a code edit."""
    cases = {p.name.split("=", 1)[1] for p in parquet_root.glob("case_id=*") if p.is_dir()}
    cases |= {p.stem for p in donors_dir.glob("*.csv")}
    return sorted(cases)


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--data-dir", type=Path, default=REPO / "data",
                    help="Repository data directory (default: <repo>/data)")
    ap.add_argument("--out-dir", type=Path, default=None,
                    help="Default: <data-dir>/tissuumaps_inputs/markers")
    ap.add_argument("--cases", nargs="*", default=None,
                    help="Restrict to these case ids (default: every donor with cell data)")
    ap.add_argument("--force", action="store_true", help="Rewrite CSVs that already exist")
    args = ap.parse_args()

    data_dir = args.data_dir.resolve()
    parquet_root = data_dir / "parquet" / "tissue"
    donors_dir = data_dir / "app_data" / "donors"
    out_dir = (args.out_dir or data_dir / "tissuumaps_inputs" / "markers").resolve()

    cases = args.cases or discover_cases(parquet_root, donors_dir)
    if not cases:
        log(f"ERROR: no donors found under {parquet_root} or {donors_dir}")
        return 1
    log(f"{len(cases)} donor(s) -> {out_dir}")

    coverage: dict[str, dict] = {}
    for case_id in cases:
        dst = out_dir / f"{case_id}.csv"
        if dst.exists() and not args.force:
            coverage[case_id] = {"status": "cached", "source": None,
                                 "bytes": dst.stat().st_size}
            log(f"  {case_id:<8} cached ({dst.stat().st_size / 1e6:.1f} MB)")
            continue

        rows, source = read_parquet_tissue(parquet_root, case_id), "parquet"
        if rows is None:
            rows, source = read_legacy_csv(donors_dir, case_id), "legacy_csv"
        if rows is None:
            # Not an error. Seven batch-2 donors have slides and no cells; the
            # viewer still shows their channels and outlines.
            coverage[case_id] = {"status": "no_cell_data", "source": None}
            log(f"  {case_id:<8} SKIPPED — no tissue partition and no legacy CSV")
            continue

        stats = write_markers(rows, dst)
        coverage[case_id] = {"status": "written", "source": source, **stats}
        note = f", {stats['n_dropped']} dropped (non-finite coords)" if stats["n_dropped"] else ""
        log(f"  {case_id:<8} {stats['n_cells']:>8,} cells from {source} "
            f"({stats['bytes'] / 1e6:.1f} MB){note}")

    out_dir.mkdir(parents=True, exist_ok=True)
    (out_dir / "coverage.json").write_text(json.dumps(coverage, indent=2, sort_keys=True))
    written = sum(1 for c in coverage.values() if c["status"] == "written")
    missing = [k for k, c in coverage.items() if c["status"] == "no_cell_data"]
    log(f"wrote {written}, cached {sum(1 for c in coverage.values() if c['status'] == 'cached')}, "
        f"no data {len(missing)}" + (f": {', '.join(missing)}" if missing else ""))
    log(f"coverage -> {out_dir / 'coverage.json'}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
