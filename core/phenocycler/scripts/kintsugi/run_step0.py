#!/usr/bin/env python
"""Step 0 driver: extract QuPath's signed islet distance, then check it against cell_region.

    python scripts/kintsugi/run_step0.py --smoke          # a few rows, proves the SQL
    python scripts/kintsugi/run_step0.py                  # whole cohort, one pass per export
    python scripts/kintsugi/run_step0.py --validate-only  # re-check what is already written

The two exports are 147 GB and 64 GB, so each is scanned exactly once and only four columns
are converted.  Nothing here writes a ``StageManifest`` or touches ``data/runs/``.
"""
from __future__ import annotations

import argparse
import json
import pathlib
import sys
import time

HERE = pathlib.Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(HERE.parents[1]))

import distance                                    # noqa: E402
from phenocycler.artifacts import QuPathCohortManifest   # noqa: E402


def format_validation_line(r: dict) -> str:
    """One donor's validation row, as the operator sees it.

    Split out of ``main`` so the non-nominal-ring branch is reachable by a test.  It fires
    only where ``ring_matches_nominal`` is false -- i.e. only on donor 6482 -- so a defect
    there is invisible on 19 of 20 donors and raises on the one donor the whole ring-width
    finding is about.  It did: the key is ``ring_width_um_median``, not ``ring_width_um``.
    """
    flag = ""
    if r["unmatched_object_id"]:
        flag = "  <-- UNMATCHED object_id"
    elif max(r["core_positive_rate"], r["peri_beyond_ring_rate"]) > distance.MAX_VIOLATION_RATE:
        flag = "  <-- OVER BAR"
    ring = ("" if r["ring_matches_nominal"]
            else f"  ring={r['ring_width_um_median']:.1f}um!")
    return (f"  {r['donor_id']:>6}  n={r['n_expression']:>9,}"
            f"  unmatched={r['unmatched_object_id']:>6,}"
            f"  core+ {r['core_positive_rate']:.4f}"
            f"  peri<=0 {r['peri_nonpositive_rate']:.4f}"
            f"  peri>ring {r['peri_beyond_ring_rate']:.4f}{ring}{flag}")


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--donors", nargs="*", default=None,
                    help="restrict to these donors (the export is still scanned once)")
    ap.add_argument("--threads", type=int, default=18)
    ap.add_argument("--memory-limit", default="48GB")
    ap.add_argument("--smoke", action="store_true",
                    help="read a bounded slice and print it; writes nothing")
    ap.add_argument("--validate-only", action="store_true")
    args = ap.parse_args()

    import duckdb

    manifest = distance.panel.PKG / "data" / "qupath_manifest.json"
    cohort = QuPathCohortManifest.read_json(manifest)
    groups = distance.csv_groups(cohort)

    if args.donors:
        keep = set(args.donors)
        groups = {csv: {i: d for i, d in m.items() if d in keep} for csv, m in groups.items()}
        groups = {csv: m for csv, m in groups.items() if m}
        if not groups:
            raise SystemExit(f"no manifest image maps to donors {sorted(keep)}")

    if not args.validate_only:
        con = duckdb.connect()
        con.execute(f"PRAGMA threads={args.threads}")
        con.execute(f"SET memory_limit='{args.memory_limit}'")
        con.execute("SET preserve_insertion_order=false")

        for csv_path, image_to_donor in sorted(groups.items()):
            donors = sorted(set(image_to_donor.values()))
            print(f"\n=== {csv_path.name}  ({len(donors)} donors: {', '.join(donors)}) ===",
                  flush=True)
            distance.assert_columns_present(con, csv_path)
            print("  distance columns present", flush=True)

            if args.smoke:
                escaped = str(csv_path).replace("'", "''")
                allowed = ", ".join(distance._sql_string(i) for i in sorted(image_to_donor))
                rows = con.execute(f"""
                    SELECT "Image", "Object ID",
                           TRY_CAST("{distance.ISLET_COL}" AS DOUBLE) AS islet_um,
                           TRY_CAST("{distance.TISSUE_COL}" AS DOUBLE) AS tissue_um
                    FROM read_csv('{escaped}', header=true, all_varchar=true,
                                  ignore_errors=true)
                    WHERE "Image" IN ({allowed}) LIMIT 5
                """).fetchall()
                for r in rows:
                    print("   ", r, flush=True)
                continue

            distance.OUT.mkdir(parents=True, exist_ok=True)
            t0 = time.time()
            con.execute(distance.extract_sql(csv_path, image_to_donor, distance.OUT))
            print(f"  wrote in {time.time() - t0:.0f}s", flush=True)

    if args.smoke:
        return 0

    report = []
    for csv_path, image_to_donor in sorted(groups.items()):
        for donor in sorted(set(image_to_donor.values())):
            r = distance.validate(donor)
            report.append(r)
            print(format_validation_line(r), flush=True)

    distance.LANE_OUT.mkdir(parents=True, exist_ok=True)
    out = distance.LANE_OUT / "distance_validation.json"
    out.write_text(json.dumps({"run_id": distance.panel.RUN_ID,
                               "max_violation_rate_bar": distance.MAX_VIOLATION_RATE,
                               "donors": report}, indent=2))
    print(f"\nwrote {out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
