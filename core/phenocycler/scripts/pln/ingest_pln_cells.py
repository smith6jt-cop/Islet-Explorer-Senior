#!/usr/bin/env python3
"""Reduce the PhenoCycler pLN QuPath export to a per-cell table the pLN arm can use.

WHY A REDUCED TABLE AND NOT A FULL INGEST INTO THE SOURCE WORKFLOW
------------------------------------------------------------------
`/data/IO60pLN/results/IO60pLN_cells.tsv` is a QuPath measurement export with exactly the
same shape as a pancreas export -- 4 compartments x 4 statistics per marker -- so
`phenocycler.cells_parquet` would parse it happily.  It must not be pointed at it:

1.  **Donor partitions would collide.**  18 of the 19 pLN images are donors that already own
    a pancreas partition under `core/phenocycler/data/cells/donor_id=<id>/`.  A pLN ingest
    into that root puts two tissues, two panels and two segmentation runs into one partition
    with no column that tells them apart.  Every downstream stage reads that root by donor.
2.  **The panel is different.**  The pancreas cohort measures 70 markers, the pLN export 56.
    19 pancreas markers are absent here (every endocrine one: GCG, INS, SST, IAPP) and 3 pLN
    markers are absent from `marker_registry.json` (IFNG, GP100, CD14).  `reconcile_headers`
    would take the union and mark the rest unavailable -- correct for two pancreas batches,
    wrong for two tissues, because it would present a cohort-wide 73-marker panel that no
    slide carries.
3.  **There is no islet.**  `cells_parquet.build_sql` derives `cell_region` and `islet_num`
    from `Parent`.  On this export `Parent` is `Annotation (Tissue)` on 100% of rows, so
    every pLN cell would be labelled `tissue` / `islet_num = NULL` -- true but useless, and
    it would look like a measurement rather than the absence of one.

So this writes its own root, `core/phenocycler/data/integration/cells_pln/`, beside the
`cells_pheno` and `cells_xen` roots the integration lane already keeps.

WHAT THIS DOES *NOT* DO
-----------------------
It does not assign node membership.  The node is the pLN's structural unit and membership is
geometric -- the `pLN_*` (and, for 6380 only, `Tissue`) annotation polygons against
`Centroid X/Y um`.  This script carries `parent_raw` for audit only; it is
`Annotation (Tissue)` on every row and can never give membership.  Deriving it belongs to
`build_pln_layer.py`, which needs the annotation polygons this export does not contain.

It does not apply the raster or duplicate rules from `[geometry_qc]`.  Those need the
whole-slide label raster and the cell outline geojson, neither of which exists for pLN in
this checkout.  `analysis_reason` records which rules were actually applied, so a reader is
never left to assume the full QC ran.

DONOR IDENTITY
--------------
`Image` -> `donor_id` comes from `core/phenocycler/data/integration/pln_node_annotations.csv`,
a tracked hand-recorded manifest, exactly as the pancreas lane takes it from
`qupath_manifest.json`.  Never from parsing the image name.  An `Image` value the manifest
does not name lands in a `donor_id=UNMAPPED` partition rather than being dropped.

WHO MAY BE WRITTEN
------------------
Two authorities have to agree.  `phenocycler.cohort` holds the cohort-wide exclusions and is
asked, never copied -- no lane restates that list.  `pln_node_annotations.csv` separately
carries a hand-recorded `excluded` flag per image, which is the maintainer's note about these
particular slides.  A donor is written only when both admit it; one that either withholds is
counted in the census and refused at write, so the refusal is a measured number rather than an
absence.

As of 2026-09-17 they disagree about 6579: the working tree reinstated it in `cohort.py`
while the tracked pLN manifest still records `excluded=True`.  This script withholds it and
says so on stdout and in the report.  Resolving that belongs in whichever file should carry
the decision, not in a builder that happens to read both.

    core/phenocycler/scripts/pln/ingest_pln_cells.py            # full ingest
    core/phenocycler/scripts/pln/ingest_pln_cells.py --census-only
    core/phenocycler/scripts/pln/ingest_pln_cells.py --limit 200000 --out /tmp/probe
"""
from __future__ import annotations

import argparse
import configparser
import csv
import json
import pathlib
import re
import sys
import time

REPO = pathlib.Path(__file__).resolve().parents[4]
sys.path.insert(0, str(REPO / "core" / "phenocycler"))

from phenocycler.cells_parquet import sanitize  # noqa: E402
from phenocycler.cohort import DONOR_EXCLUSIONS, filter_eligible_donors  # noqa: E402
from phenocycler.marker_registry import load_registry  # noqa: E402

SOURCE = pathlib.Path("/data/IO60pLN/results/IO60pLN_cells.tsv")
NODES = REPO / "core/phenocycler/data/integration/pln_node_annotations.csv"
OUT_ROOT = REPO / "core/phenocycler/data/integration/cells_pln"
REPORT = REPO / "core/phenocycler/data/integration/cells_pln_ingest_report.json"
CONFIG = REPO / "core/phenocycler/config.ini"

#: QuPath appends this to the project entry name in the `Image` column.  Verified present on
#: all 19 distinct values before it was relied on; the mapping still matches the bare file
#: name too, so a re-export without the suffix is read rather than sent to UNMAPPED.
RESOLUTION_SUFFIX = " - resolution #1"

#: Context columns, QuPath name -> output name.  Same spellings the pancreas lane uses
#: (`phenocycler.cells_parquet.build_sql`, `scripts/typing/rebuild_donor_from_revised_export`)
#: so a pLN frame and a pancreas frame can be put side by side without a rename step.
CONTEXT = {
    "Image": "image",
    "Object ID": "object_id",
    "Centroid X µm": "X_centroid",
    "Centroid Y µm": "Y_centroid",
    "Parent": "parent_raw",
    "Classification": "qupath_classification",
}
#: Geometry-QC inputs.  Carried so every `[geometry_qc]` rule that *can* be applied without
#: the label raster is reproducible from this table alone -- the 43.8 GB export is never read
#: a second time to re-decide eligibility.
SHAPE = {
    "Nucleus: Area µm^2": "nucleus_area_um2",
    "Cell: Area µm^2": "cell_area",
    "Cell: Solidity": "cell_solidity",
    "Nucleus/Cell area ratio": "nucleus_cell_area_ratio",
}
NUMERIC_CONTEXT = ("X_centroid", "Y_centroid")


def read_header(path: pathlib.Path) -> list[str]:
    with path.open(newline="", encoding="utf-8", errors="replace") as fh:
        return next(csv.reader(fh, delimiter="\t"))


def compartment_means(header: list[str], compartment: str) -> dict[str, str]:
    """`<compartment>: <marker>: Mean` -> sanitized marker name, in header order."""
    pattern = re.compile(rf"^{re.escape(compartment)}: (.+): Mean$")
    out: dict[str, str] = {}
    for column in header:
        match = pattern.match(column)
        if match:
            out[column] = sanitize(match.group(1))
    return out


def donor_map() -> tuple[dict[str, str], dict[str, str], dict[str, dict]]:
    """`Image` value -> donor_id, `Image` value -> roi, and the manifest row behind each image.

    Both the bare QuPath entry name and the `- resolution #1` form map to the same donor, so
    the mapping does not depend on which one a given export writes.
    """
    rows = {}
    with NODES.open(newline="", encoding="utf-8") as fh:
        for row in csv.DictReader(line for line in fh if not line.startswith("#")):
            rows[row["image_file"]] = row
    mapping = {}
    roi = {}
    for image_file, row in rows.items():
        for key in (image_file, image_file + RESOLUTION_SUFFIX):
            mapping[key] = row["donor_id"]
            # The manifest records the roi per image and every row reads `pln_1` today. It is
            # read rather than written as a literal so a donor who later contributes a second
            # node gets its own partition instead of being merged into the first one.
            roi[key] = row["roi"]
    return mapping, roi, rows


def thresholds() -> dict[str, float]:
    cp = configparser.ConfigParser()
    cp.read(CONFIG)
    g = cp["geometry_qc"]
    return {
        "min_cell_area": g.getfloat("min_cell_area"),
        "max_area_median_multiple": g.getfloat("max_area_median_multiple"),
        "min_cell_solidity": g.getfloat("min_cell_solidity"),
        "no_cytoplasm_nc_ratio": g.getfloat("no_cytoplasm_nc_ratio"),
    }


def sql_string(value: str) -> str:
    return "'" + str(value).replace("'", "''") + "'"


def image_case(mapping: dict[str, str], *, column: str = '"Image"', default: str = "UNMAPPED") -> str:
    """CASE over the exact `Image` value.  ELSE is a visible partition, not a silent drop."""
    arms = " ".join(
        f"WHEN {column} = {sql_string(image)} THEN {sql_string(value)}"
        for image, value in sorted(mapping.items())
    )
    return f"CASE {arms} ELSE {sql_string(default)} END"


def absent(column: str) -> str:
    """SQL for "this measurement is unavailable".

    QuPath writes an empty field for a measurement it did not make and the literal ``NaN``
    for one it attempted and could not compute; ``TRY_CAST`` turns the first into NULL and
    the second into NaN. Both mean *no value*, and only testing NULL would read a NaN as a
    number and let it through every threshold comparison as false.
    """
    return f"({column} IS NULL OR NOT isfinite({column}))"


def reason_sql(thr: dict[str, float], *, estimation: bool) -> str:
    """The exclusion ladder, as one SQL CASE. The first rule a cell fails is its reason.

    Same rules, same order and same reason names as
    ``scripts/typing/rebuild_donor_from_revised_export.eligibility()``, plus one this export
    forced into existence.

    ``no_measurements_exported`` comes first and is NOT a variant of
    ``missing_nucleus_geometry``. Four images -- 6551, 6566, 6593, 6623 -- carry rows with all
    918 fields present and only seven of them filled: Image, Object ID, Object type, Parent,
    ROI and the two centroids. No area, no solidity, no intensity, in any compartment, on any
    marker, on 100% of their rows. Those are detections with no measurements, not cells whose
    nucleus was missed, and folding them into the nucleus rule would report a segmentation
    outcome for what is an export failure. The distinction is the whole reason a card can say
    ``unavailable(no_measurements_in_pln_export)`` instead of quietly dropping 2.79M cells.

    There is deliberately no rule that infers a value for any of them. A missing measurement
    is unavailable -- not zero, not negative.
    """
    nothing = (
        f"({absent('cell_area')} AND {absent('nucleus_area_um2')} AND {absent('DAPI')})"
    )
    miss = absent("nucleus_area_um2")
    below = f"({absent('cell_area')} OR cell_area < {thr['min_cell_area']})"
    above = f"(cell_area > {thr['max_area_median_multiple']} * med_cell_area)"
    low = f"(cell_solidity < {thr['min_cell_solidity']})"
    arms = [
        (nothing, "no_measurements_exported"),
        (miss, "missing_nucleus_geometry"),
        (below, "cell_area_below_min"),
        (above, "cell_area_above_max"),
        (low, "low_solidity"),
    ]
    if estimation:
        arms.append(
            (f"(nucleus_cell_area_ratio > {thr['no_cytoplasm_nc_ratio']})", "no_cytoplasm")
        )
    return "CASE " + " ".join(
        f"WHEN {test} THEN {sql_string(name)}" for test, name in arms
    ) + " ELSE '' END"


def build(args) -> int:
    import duckdb

    if not SOURCE.exists():
        raise SystemExit(f"[err] pLN export not found: {SOURCE}")
    header = read_header(SOURCE)
    missing = [c for c in list(CONTEXT) + list(SHAPE) if c not in header]
    if missing:
        raise SystemExit(f"[err] pLN export is missing required columns: {missing}")

    cell_mm = compartment_means(header, "Cell")
    nucleus_mm = compartment_means(header, "Nucleus")
    if sorted(cell_mm.values()) != sorted(nucleus_mm.values()):
        raise SystemExit(
            "[err] Cell and Nucleus marker means do not cover the same markers: "
            f"{sorted(set(cell_mm.values()) ^ set(nucleus_mm.values()))}"
        )
    registry = {marker.name for marker in load_registry().markers}
    in_registry = sorted(m for m in cell_mm.values() if m in registry)
    off_registry = sorted(m for m in cell_mm.values() if m not in registry)

    mapping, roi_map, manifest_rows = donor_map()
    all_donors = sorted({row["donor_id"] for row in manifest_rows.values()})

    # TWO AUTHORITIES, AND THIS LANE REFUSES TO PICK BETWEEN THEM SILENTLY.
    #
    # `phenocycler.cohort` holds the cohort-wide exclusions and no lane may restate that list,
    # so it is asked rather than copied.  `pln_node_annotations.csv` separately carries a
    # hand-recorded `excluded` flag per image -- the maintainer's own note about *these* slides,
    # written when the node annotations were recorded.  They are not the same statement and
    # they can disagree.
    #
    # When they do, the safe direction is the strict one: a donor either authority withholds is
    # not written.  Publishing on the permissive reading would put cells into a stable root on
    # the strength of whichever file happened to be edited last, and a reader of the parquet
    # cannot see which that was.  The disagreement is named in the report and on stdout instead,
    # for a human to settle in the file that should carry the decision.
    by_cohort = set(filter_eligible_donors(all_donors))
    by_manifest = {
        row["donor_id"] for row in manifest_rows.values()
        if row["excluded"].strip().lower() not in {"true", "1", "yes"}
    }
    disagreement = sorted(by_cohort ^ by_manifest)
    if disagreement:
        for donor in disagreement:
            where = "cohort.py" if donor in by_cohort else "pln_node_annotations.csv"
            print(f"[warn] {donor}: admitted by {where} and withheld by the other authority; "
                  "withheld here until they agree")
    eligible_donors = sorted(by_cohort & by_manifest)
    refused = sorted(set(all_donors) - set(eligible_donors))
    eligible_images = sorted(
        image for image, donor in mapping.items() if donor in eligible_donors
    )
    thr = thresholds()

    out_root = pathlib.Path(args.out) if args.out else OUT_ROOT
    # A probe writes its report beside its own output. The canonical report must describe the
    # canonical root and nothing else, or a --limit run silently replaces the real measurement.
    report_path = (
        REPORT if out_root == OUT_ROOT
        else out_root.with_name(out_root.name + "_ingest_report.json")
    )
    selects = [f'"{orig}" AS "{name}"' for orig, name in CONTEXT.items() if name not in NUMERIC_CONTEXT]
    selects += [f'TRY_CAST("{orig}" AS DOUBLE) AS "{CONTEXT[orig]}"'
                for orig in CONTEXT if CONTEXT[orig] in NUMERIC_CONTEXT]
    selects += [f'TRY_CAST("{orig}" AS DOUBLE) AS "{name}"' for orig, name in SHAPE.items()]
    selects += [f'TRY_CAST("{orig}" AS DOUBLE) AS "{name}"' for orig, name in cell_mm.items()]
    selects += [f'TRY_CAST("{orig}" AS DOUBLE) AS "Nucleus__{name}"'
                for orig, name in nucleus_mm.items()]

    read = (
        f"read_csv({sql_string(str(SOURCE))}, header=true, delim='\\t', all_varchar=true, "
        "ignore_errors=true, store_rejects=true, rejects_table='pln_rejects', "
        "rejects_scan='pln_reject_scans', rejects_limit=10000)"
    )
    inner = "SELECT\n  " + ",\n  ".join(selects) + f"\nFROM {read}"
    if args.limit:
        inner += f"\nLIMIT {args.limit}"

    con = duckdb.connect()
    con.execute(f"PRAGMA threads={args.threads}")
    con.execute(f"PRAGMA memory_limit='{args.memory_limit}'")
    con.execute(f"PRAGMA temp_directory={sql_string(args.temp_dir)}")

    started = time.time()
    print(f"[info] scanning {SOURCE} ({SOURCE.stat().st_size / 1e9:.2f} GB, {len(header)} columns)")
    print(f"[info] keeping {len(selects)} columns: {len(cell_mm)} Cell means + "
          f"{len(nucleus_mm)} Nucleus means + {len(CONTEXT) + len(SHAPE)} context/geometry")
    con.execute(f"CREATE TABLE pln_cells AS {inner}")
    scanned = time.time() - started
    n_rows = con.execute("SELECT count(*) FROM pln_cells").fetchone()[0]
    print(f"[info] {n_rows:,} rows in {scanned / 60:.1f} min")

    rejects = con.execute("SELECT count(*) FROM pln_rejects").fetchone()[0]
    if rejects:
        print(f"[warn] CSV reader rejected {rejects} row(s); see the report")

    census = con.execute(
        'SELECT image, count(*) AS n_rows, '
        'count(*) FILTER (WHERE nucleus_area_um2 IS NULL OR NOT isfinite(nucleus_area_um2)) '
        '  AS n_no_nucleus, '
        # A row with no cell area, no nucleus area and no Cell DAPI mean carries no
        # measurement at all. Counted separately from the nucleus rule because the two have
        # different causes and different remedies -- see reason_sql().
        'count(*) FILTER (WHERE (cell_area IS NULL OR NOT isfinite(cell_area)) '
        '  AND (nucleus_area_um2 IS NULL OR NOT isfinite(nucleus_area_um2)) '
        '  AND (DAPI IS NULL OR NOT isfinite(DAPI))) AS n_no_measurements, '
        'median(cell_area) AS median_cell_area '
        'FROM pln_cells GROUP BY image ORDER BY image'
    ).fetchall()
    unmapped = [row for row in census if row[0] not in mapping]
    print(f"[info] {len(census)} distinct Image values; {len(unmapped)} unmapped")

    # Eligibility.  Same rules, same order, same reason names as
    # scripts/typing/rebuild_donor_from_revised_export.eligibility(), minus the raster and
    # duplicate rules that need artifacts pLN does not have in this checkout.
    con.execute(
        "CREATE TABLE image_median AS "
        "SELECT image, median(cell_area) AS med_cell_area FROM pln_cells "
        "WHERE isfinite(cell_area) GROUP BY image"
    )
    analysis_reason = reason_sql(thr, estimation=False)
    estimation_reason = reason_sql(thr, estimation=True)

    allowed = ", ".join(sql_string(image) for image in eligible_images)
    write_sql = (
        "SELECT "
        f"{image_case(mapping, column='c.image')} AS donor_id, "
        f"{image_case(roi_map, column='c.image')} AS roi, "
        "c.image AS image, "
        "c.* EXCLUDE (image), "
        f"{analysis_reason} AS analysis_reason, "
        f"({analysis_reason}) = '' AS qc_analysis_eligible, "
        f"{estimation_reason} AS estimation_reason, "
        f"({estimation_reason}) = '' AS qc_estimation_eligible "
        "FROM pln_cells c LEFT JOIN image_median m USING (image) "
        f"WHERE c.image IN ({allowed})"
    )
    if args.census_only:
        print("[info] --census-only: nothing written")
    else:
        out_root.parent.mkdir(parents=True, exist_ok=True)
        con.execute(
            f"COPY ({write_sql}) TO {sql_string(str(out_root))} "
            "(FORMAT PARQUET, PARTITION_BY (donor_id, roi), OVERWRITE_OR_IGNORE, "
            "COMPRESSION ZSTD)"
        )
        print(f"[info] wrote {out_root}")

    written = con.execute(
        f"SELECT count(*), count(*) FILTER (WHERE qc_analysis_eligible) FROM ({write_sql})"
    ).fetchone()

    report_per_image = [
        {
            "image": row[0],
            "donor_id": mapping.get(row[0], "UNMAPPED"),
            "n_rows_in_export": int(row[1]),
            "n_no_nucleus": int(row[2]),
            "n_no_measurements_exported": int(row[3]),
            "measurements_state": (
                "unavailable(no_measurements_in_pln_export)" if row[3] == row[1]
                else "measured" if row[3] == 0
                else "partial"
            ),
            "median_cell_area_um2": None if row[4] is None else round(float(row[4]), 4),
            "n_cells_recorded_in_manifest": (
                int(manifest_rows[row[0].removesuffix(RESOLUTION_SUFFIX)]["n_cells"])
                if row[0].removesuffix(RESOLUTION_SUFFIX) in manifest_rows else None
            ),
            "excluded_donor": mapping.get(row[0]) in DONOR_EXCLUSIONS,
        }
        for row in census
    ]

    report = {
        "probe_only_limit": args.limit,
        "generated_at_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
        "source": {
            "path": str(SOURCE),
            "size_bytes": SOURCE.stat().st_size,
            "columns": len(header),
            "rows": int(n_rows),
            "csv_rejected_rows": int(rejects),
            "scan_seconds": round(scanned, 1),
        },
        "output": {
            "root": str(out_root),
            "layout": "donor_id=<id>/roi=pln_1/data_*.parquet",
            "donor_id_is_a_partition_key_only": (
                "no donor_id column is written inside the files, so "
                "pyarrow.dataset(root, partitioning='hive') reads this root without the "
                "large_string collision documented for expression/cells_xen/cells_pheno"
            ),
            "read_donor_id_as_a_string": (
                "bare partitioning='hive' INFERS donor_id as int32 here, because every pLN "
                "donor id happens to be all digits. That silently turns '0112'-style ids into "
                "integers and will not join against the string donor_id in cells_pheno or "
                "typing/cells. Read it with an explicit schema: "
                "ds.partitioning(pa.schema([('donor_id', pa.string()), ('roi', pa.string())]), "
                "flavor='hive')"
            ),
            "columns_kept": len(selects) + 6,
            "rows_written": int(written[0]),
            "rows_analysis_eligible": int(written[1]),
        },
        "donor_mapping": {
            "source": str(NODES),
            "method": "exact Image match against a tracked hand-recorded manifest",
            "images_in_manifest": len(manifest_rows),
            "donors_in_manifest": all_donors,
            "donors_written": eligible_donors,
            "donors_refused": {
                d: (
                    DONOR_EXCLUSIONS.get(d)
                    or "withheld by pln_node_annotations.csv (excluded=True)"
                )
                for d in refused
            },
            "authorities_disagree_about": {
                d: (
                    "admitted by phenocycler.cohort, withheld by pln_node_annotations.csv"
                    if d in by_cohort else
                    "admitted by pln_node_annotations.csv, withheld by phenocycler.cohort"
                )
                for d in disagreement
            },
            "exclusion_rule": (
                "a donor is written only when phenocycler.cohort AND the tracked pLN manifest "
                "both admit it. The exclusion list is never restated here -- both authorities "
                "are read"
            ),
            "unmapped_images": [{"image": r[0], "n_rows": int(r[1])} for r in unmapped],
        },
        "panel": {
            "markers_in_export": len(cell_mm),
            "compartments_in_export": ["Cell", "Cytoplasm", "Membrane", "Nucleus"],
            "compartments_kept": ["Cell", "Nucleus"],
            "statistics_in_export": ["Mean", "Min", "Max", "Std.Dev."],
            "statistics_kept": ["Mean"],
            "in_marker_registry": in_registry,
            "not_in_marker_registry": off_registry,
            "registry_markers_absent_from_this_export": sorted(registry - set(cell_mm.values())),
            "name_collisions_not_asserted": (
                "the export writes PDL1; the registry writes PD_L1. They are not mapped onto "
                "each other here -- that is a panel-identity question for the maintainer, and "
                "guessing it would put an unverified marker identity into a stable root"
            ),
        },
        "geometry_qc": {
            "thresholds": thr,
            "thresholds_source": str(CONFIG) + " [geometry_qc]",
            "rules_applied": [
                "missing_nucleus_geometry", "cell_area_below_min",
                "cell_area_above_max", "low_solidity", "no_cytoplasm (estimation only)",
            ],
            "rules_not_applied": {
                "raster_area_mismatch/raster_empty/nucleus_raster_empty":
                    "unavailable(no_whole_slide_label_raster_for_pln)",
                "duplicate_policy":
                    "unavailable(no_cell_outline_geojson_for_pln)",
            },
            "median_cell_area_is_per_image": True,
        },
        "per_image": report_per_image,
        "measurement_availability": {
            "rule": (
                "a row whose cell area, nucleus area and Cell DAPI mean are all absent carries "
                "no measurement of any kind; it is a detection, not a measured cell"
            ),
            "images_with_no_measurements": [
                {"image": e["image"], "donor_id": e["donor_id"], "n_rows": e["n_rows_in_export"]}
                for e in report_per_image if e["measurements_state"].startswith("unavailable")
            ],
            "n_cells_with_no_measurements": sum(
                e["n_no_measurements_exported"] for e in report_per_image
            ),
            "for_the_availability_layer": (
                "phenocycler_pln is unavailable(no_measurements_in_pln_export) for those donors. "
                "Their cells are written with their ids and centroids so the count is visible and "
                "a re-export can be joined on object_id; no value is imputed for any of them"
            ),
        },
        "parent_column": (
            "carried as parent_raw for audit only. It is 'Annotation (Tissue)' on every row, "
            "so it cannot give node membership; membership is geometric and is not derived here"
        ),
    }
    report_path.parent.mkdir(parents=True, exist_ok=True)
    report_path.write_text(json.dumps(report, indent=2) + "\n")
    print(f"[info] report -> {report_path}")
    for entry in report["per_image"]:
        delta = (
            "" if entry["n_cells_recorded_in_manifest"] is None
            else f"  manifest={entry['n_cells_recorded_in_manifest']:,}"
                 f"  delta={entry['n_rows_in_export'] - entry['n_cells_recorded_in_manifest']:+d}"
        )
        flag = "  EXCLUDED" if entry["excluded_donor"] else ""
        print(f"  {entry['donor_id']:>8}  {entry['n_rows_in_export']:>10,} rows"
              f"  no_nucleus={entry['n_no_nucleus']:>8,}"
              f"  no_measurements={entry['n_no_measurements_exported']:>9,}{delta}{flag}")
    return 0


def main(argv=None) -> int:
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--out", default=None, help="output root (default: cells_pln beside cells_pheno)")
    p.add_argument("--limit", type=int, default=None, help="stop after N source rows (probe)")
    p.add_argument("--census-only", action="store_true", help="measure, write no parquet")
    p.add_argument("--threads", type=int, default=16)
    p.add_argument("--memory-limit", default="96GB")
    p.add_argument("--temp-dir", default="/data/tmp/duckdb_pln")
    return build(p.parse_args(argv))


if __name__ == "__main__":
    raise SystemExit(main())
