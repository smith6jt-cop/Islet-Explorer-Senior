#!/usr/bin/env python3
"""Reduce the PhenoCycler pLN QuPath export to a per-cell table the pLN arm can use.

WHY A REDUCED TABLE AND NOT A FULL INGEST INTO THE SOURCE WORKFLOW
------------------------------------------------------------------
The pLN export (`data/pLN/IO60pLN_cells.tsv`) is a QuPath measurement export with exactly the
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

They disagreed about 6579 on 2026-09-17 -- `cohort.py` reinstated it at 14:16, the tracked pLN
manifest still said `excluded=True`, and the 13:31 run withheld it. The manifest was corrected
the same day, so both now admit 6579 and a re-run writes its 1,231,818 cells. The rule that
produced the withholding is unchanged and deliberately so; what changed is that the withholding
is now printed as a named donor in a ledger rather than being visible only as a partition that
is not there.

A DONOR THE AUTHORITIES ADMIT AND THE EXPORT DOES NOT CONTAIN IS ALSO NAMED. The per-image
census can only report images the scan saw, so a donor whose `Image` value never appears makes
no row anywhere -- indistinguishable, in a per-image listing, from a donor nobody asked about.
The ledger is built from the manifest's donor list instead, and prints one line per donor:
ingested with its counts, refused with which authority withheld it, or absent with the fact
that it was admitted and still produced nothing.

WHAT "NO INTENSITY MEASUREMENT" MEANS, AND WHY IT IS RE-MEASURED EVERY RUN
-------------------------------------------------------------------------
Four images -- 6551, 6566, 6593, 6623 -- carried no measurement of any kind when this was last
measured (2026-09-17 13:31, against the 43.76 GB export of 2026-09-16 13:29: `n_no_measurements`
equal to `n_rows` on all four, 2,786,845 cells). That verdict is about COLUMNS, and repeating it
would be inheriting it. So each run counts, per image, how many of the kept measurement columns
have no finite value anywhere, names them, and reports `columns_state` beside the row-level
`measurements_state`. When the two disagree the run says so and neither may be quoted. Nothing
is imputed for an empty column: a missing measurement is unavailable, not zero.

THERE ARE TWO COPIES OF THE EXPORT on this machine, `data/pLN/IO60pLN_cells.tsv` and
`/data/IO60pLN/results/IO60pLN_cells.tsv` -- different inodes, identical size and mtime. The run
picks the newest, prints every candidate, and records the choice in the report, because "which
file did that number come from" is otherwise unanswerable.

    core/multiplex5000/scripts/pln/ingest_pln_cells.py            # full ingest
    core/multiplex5000/scripts/pln/ingest_pln_cells.py --census-only
    core/multiplex5000/scripts/pln/ingest_pln_cells.py --source data/pLN/IO60pLN_cells.tsv
    core/multiplex5000/scripts/pln/ingest_pln_cells.py --limit 200000 --out /tmp/probe

A DIRECTORY OF PER-SLIDE PARQUET TABLES IS A SOURCE TOO (the v3 QuPath export, 2026-10-01:
`~/IO60pLN/results/v3/export/cells/`, one table per slide with QuPath's own column names). It is
read only when `--source` names it -- never picked by the newest-candidate rule -- and it goes
through the same SELECT, census, ledger and report as the TSV. A column a slide's table lacks is
NULL for that slide (unavailable, not zero), and the report lists which slide lacks which column.

    core/multiplex5000/scripts/pln/ingest_pln_cells.py --source ~/IO60pLN/results/v3/export/cells
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
sys.path.insert(0, str(REPO / "core" / "multiplex5000"))

from phenocycler.cells_parquet import sanitize  # noqa: E402
from phenocycler.cohort import DONOR_EXCLUSIONS, filter_eligible_donors  # noqa: E402
from phenocycler.marker_registry import load_registry  # noqa: E402

#: Where the pLN export may be, newest first at equal standing. Two copies of the 43.76 GB
#: re-export exist on this machine -- `data/pLN/IO60pLN_cells.tsv` and
#: `/data/IO60pLN/results/IO60pLN_cells.tsv`, different inodes, identical size and mtime
#: (2026-09-16 13:29). Which one a run read is not recoverable from the report unless the run
#: says so, so `resolve_source()` picks the newest, prints every candidate it saw with its size
#: and mtime, and records the choice. `--source` overrides it.
SOURCE_CANDIDATES = (
    REPO / "data/pLN/IO60pLN_cells.tsv",
    pathlib.Path("/data/IO60pLN/results/IO60pLN_cells.tsv"),
)
NODES = REPO / "core/phenocycler/data/integration/pln_node_annotations.csv"
OUT_ROOT = REPO / "core/phenocycler/data/integration/cells_pln"
REPORT = REPO / "core/phenocycler/data/integration/cells_pln_ingest_report.json"
CONFIG = REPO / "core/multiplex5000/config.ini"

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


def resolve_source(override: str | None = None) -> tuple[pathlib.Path, list[dict]]:
    """(the export this run reads, every candidate that exists).

    The newest candidate wins; ties keep SOURCE_CANDIDATES order. Nothing is inferred from a
    path looking canonical -- the choice and the alternatives are both printed and both recorded,
    because a second copy of a 43.76 GB file is exactly the kind of thing that makes two runs
    disagree with no way to tell which read what.
    """
    seen = []
    for candidate in ((pathlib.Path(override).expanduser(),) if override else SOURCE_CANDIDATES):
        if candidate.is_dir():
            tables = parquet_tables(candidate)
            if not tables:
                raise SystemExit(f"[err] {candidate} is a directory with no *.parquet table in it")
            newest = max(t.stat().st_mtime for t in tables)
            seen.append({"path": str(candidate), "size_bytes": sum(t.stat().st_size for t in tables),
                         "mtime_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime(newest)),
                         "mtime": newest, "parquet_tables": len(tables)})
        elif candidate.exists():
            stat = candidate.stat()
            seen.append({"path": str(candidate), "size_bytes": stat.st_size,
                         "mtime_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ",
                                                    time.gmtime(stat.st_mtime)),
                         "mtime": stat.st_mtime})
    if not seen:
        tried = ", ".join(str(c) for c in ((override,) if override else SOURCE_CANDIDATES))
        raise SystemExit(f"[err] pLN export not found; tried {tried}")
    chosen = max(seen, key=lambda entry: entry["mtime"])
    return pathlib.Path(chosen["path"]), seen


def parquet_tables(directory: pathlib.Path) -> list[pathlib.Path]:
    return sorted(directory.glob("*.parquet"))


def parquet_columns(directory: pathlib.Path) -> tuple[list[str], dict[str, list[str]]]:
    """(every column any table carries, in first-seen order; table name -> columns it lacks)."""
    import pyarrow.parquet as pq

    per_table = {t.name: pq.ParquetFile(t).schema_arrow.names for t in parquet_tables(directory)}
    union: list[str] = []
    for names in per_table.values():
        union.extend(name for name in names if name not in union)
    lacking = {table: [c for c in union if c not in set(names)] for table, names in per_table.items()}
    return union, {table: cols for table, cols in lacking.items() if cols}


def read_header(path: pathlib.Path) -> list[str]:
    if path.is_dir():
        return parquet_columns(path)[0]
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

    source, candidates = resolve_source(args.source)
    for entry in candidates:
        mark = "->" if entry["path"] == str(source) else "  "
        print(f"[info] {mark} candidate {entry['path']}  "
              f"{entry['size_bytes'] / 1e9:.2f} GB  {entry['mtime_utc']}")
    if len(candidates) > 1:
        print(f"[warn] {len(candidates)} copies of the pLN export exist; reading the newest")
    header = read_header(source)
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

    from_parquet = source.is_dir()
    source_bytes = (sum(t.stat().st_size for t in parquet_tables(source)) if from_parquet
                    else source.stat().st_size)
    columns_lacking = parquet_columns(source)[1] if from_parquet else {}
    for table, cols in columns_lacking.items():
        print(f"[warn] {table} lacks {len(cols)} column(s) other tables carry; they are NULL "
              f"for that slide: {cols[:6]}{' ...' if len(cols) > 6 else ''}")
    if from_parquet:
        # union_by_name: a column one slide's table lacks is NULL for that slide, never shifted
        read = (f"read_parquet({sql_string(str(source / '*.parquet'))}, union_by_name=true)")
    else:
        read = (
            f"read_csv({sql_string(str(source))}, header=true, delim='\\t', all_varchar=true, "
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
    print(f"[info] scanning {source} ({source_bytes / 1e9:.2f} GB, {len(header)} columns"
          + (f", {len(parquet_tables(source))} Parquet tables)" if from_parquet else ")"))
    print(f"[info] keeping {len(selects)} columns: {len(cell_mm)} Cell means + "
          f"{len(nucleus_mm)} Nucleus means + {len(CONTEXT) + len(SHAPE)} context/geometry")
    con.execute(f"CREATE TABLE pln_cells AS {inner}")
    scanned = time.time() - started
    n_rows = con.execute("SELECT count(*) FROM pln_cells").fetchone()[0]
    print(f"[info] {n_rows:,} rows in {scanned / 60:.1f} min")

    # a Parquet table has no malformed line to reject; the count stays in the report as a zero that
    # was measured on a format where it cannot be anything else
    rejects = 0 if from_parquet else con.execute("SELECT count(*) FROM pln_rejects").fetchone()[0]
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

    # WHICH COLUMNS THIS EXPORT ACTUALLY FILLED, PER IMAGE, MEASURED HERE.
    #
    # The row-level rule above answers "does this row carry any measurement". It cannot answer
    # "did this re-export fix the donors that carried none", because a verdict like
    # `no intensity measurement` is a statement about columns, and repeating an earlier run's
    # statement is inheriting it rather than measuring it. So every kept measurement column is
    # counted per image, and a column with no finite value anywhere on that image is named.
    # Nothing is imputed for an empty column: absent is unavailable, not zero.
    measurement_columns = (
        list(SHAPE.values()) + list(cell_mm.values())
        + [f"Nucleus__{name}" for name in nucleus_mm.values()]
    )
    finite_sql = ", ".join(
        f'count(*) FILTER (WHERE isfinite("{column}")) AS "n_{index}"'
        for index, column in enumerate(measurement_columns)
    )
    finite_rows = con.execute(
        f"SELECT image, {finite_sql} FROM pln_cells GROUP BY image ORDER BY image"
    ).fetchall()
    all_null: dict[str, list[str]] = {}
    for row in finite_rows:
        all_null[row[0]] = [
            column for column, count in zip(measurement_columns, row[1:], strict=True)
            if not count]

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

    # `IN ()` is a syntax error, and an empty eligible set is a real state (every donor
    # withheld). It becomes a predicate that selects nothing, so the run still reports.
    allowed = ", ".join(sql_string(image) for image in eligible_images) or "NULL"
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

    # One pass for both the totals and the per-donor counts: write_sql re-derives the whole
    # eligibility ladder, so it is evaluated once rather than once per question.
    per_donor_written = con.execute(
        "SELECT donor_id, count(*), count(*) FILTER (WHERE qc_analysis_eligible) "
        f"FROM ({write_sql}) GROUP BY donor_id"
    ).fetchall()
    rows_by_donor = {row[0]: int(row[1]) for row in per_donor_written}
    written = (sum(int(row[1]) for row in per_donor_written),
               sum(int(row[2]) for row in per_donor_written))

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
            # Measured on this export, column by column. `columns_state` and
            # `measurements_state` answer different questions and are reported separately:
            # one is "every row is empty", the other is "this column is empty on every row".
            "n_measurement_columns": len(measurement_columns),
            "n_all_null_measurement_columns": len(all_null.get(row[0], [])),
            "all_null_measurement_columns": all_null.get(row[0], []),
            "columns_state": (
                "unavailable(every_measurement_column_all_null)"
                if len(all_null.get(row[0], [])) == len(measurement_columns)
                else "measured" if not all_null.get(row[0])
                else f"partial({len(all_null.get(row[0], []))}"
                     f"_of_{len(measurement_columns)}_columns_all_null)"
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

    # The row rule and the column census must agree about a donor that carries nothing. When
    # they do not, the export changed under one of them and neither verdict may be repeated.
    for entry in report_per_image:
        row_empty = entry["measurements_state"].startswith("unavailable")
        column_empty = entry["columns_state"].startswith("unavailable")
        if row_empty != column_empty:
            print(f"[warn] {entry['donor_id']}: rows say "
                  f"{entry['measurements_state']} and columns say {entry['columns_state']}; "
                  "re-measure before either is quoted")

    # ONE LINE PER DONOR THE MANIFEST NAMES, INGESTED OR NOT. A donor that never reaches the
    # export produces no census row at all, so a per-image loop cannot say it is missing -- and
    # that silence is how 6579 went unnoticed for a day. The ledger is built from the manifest's
    # donors, not from what happened to be scanned.
    by_donor_image = {}
    for entry in report_per_image:
        by_donor_image.setdefault(entry["donor_id"], []).append(entry)
    ledger = []
    for donor in all_donors:
        images = by_donor_image.get(donor, [])
        if donor not in by_cohort:
            outcome, why = "refused", DONOR_EXCLUSIONS.get(
                donor, "withheld by phenocycler.cohort")
        elif donor not in by_manifest:
            outcome, why = "refused", f"withheld by {NODES.name} (excluded=True)"
        elif not images:
            outcome, why = "absent", (
                "admitted by both authorities and carries no row in this export; its Image "
                f"value is not one of the {len(census)} the scan saw")
        elif not rows_by_donor.get(donor):
            outcome, why = "absent", "admitted, matched an image, and wrote no row"
        else:
            outcome = "ingested"
            why = "; ".join(f"{e['image']}: {e['columns_state']}" for e in images)
        ledger.append({
            "donor_id": donor,
            "outcome": outcome,
            "reason": why,
            "n_rows_in_export": sum(e["n_rows_in_export"] for e in images) or None,
            "n_rows_written": int(rows_by_donor.get(donor, 0)) or None,
            "n_cells_recorded_in_manifest": sum(
                e["n_cells_recorded_in_manifest"] or 0 for e in images) or None,
        })

    report = {
        "probe_only_limit": args.limit,
        "generated_at_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
        "source": {
            "path": str(source),
            "format": "parquet directory" if from_parquet else "tsv",
            "size_bytes": source_bytes,
            "columns_lacking_by_table": columns_lacking,
            "candidates_seen": [
                {k: v for k, v in entry.items() if k != "mtime"} for entry in candidates
            ],
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
        "per_donor": ledger,
        "columns_measured_on_this_export": {
            "rule": (
                "a measurement column with no finite value on an image is named for that image. "
                "It is measured here, on this export, and never carried over from an earlier "
                "run's verdict -- a re-export can fill a column that was empty before, and the "
                "only way to know is to look"
            ),
            "n_measurement_columns": len(measurement_columns),
            "donors_with_every_column_all_null": sorted(
                e["donor_id"] for e in report_per_image
                if e["columns_state"].startswith("unavailable")
            ),
            "donors_with_some_columns_all_null": {
                e["donor_id"]: e["all_null_measurement_columns"]
                for e in report_per_image
                if e["all_null_measurement_columns"]
                and not e["columns_state"].startswith("unavailable")
            },
        },
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
              f"  no_measurements={entry['n_no_measurements_exported']:>9,}"
              f"  all_null_cols={entry['n_all_null_measurement_columns']:>3}"
              f"/{entry['n_measurement_columns']}{delta}{flag}")

    print(f"\n[info] per donor ({len(ledger)} in {NODES.name}) -- ingested, or why not:")
    for entry in ledger:
        wrote = "-" if entry["n_rows_written"] is None else f"{entry['n_rows_written']:,}"
        print(f"  {entry['donor_id']:>8}  {entry['outcome']:<9} wrote={wrote:>11}  "
              f"{entry['reason']}")
    for outcome in ("refused", "absent"):
        named = [e["donor_id"] for e in ledger if e["outcome"] == outcome]
        if named:
            print(f"[warn] {outcome}: {' '.join(named)}")
    return 0


def main(argv=None) -> int:
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--source", default=None,
                   help="the pLN export to read; default is the newest of "
                        + ", ".join(str(c) for c in SOURCE_CANDIDATES))
    p.add_argument("--out", default=None, help="output root (default: cells_pln beside cells_pheno)")
    p.add_argument("--limit", type=int, default=None, help="stop after N source rows (probe)")
    p.add_argument("--census-only", action="store_true", help="measure, write no parquet")
    p.add_argument("--threads", type=int, default=16)
    p.add_argument("--memory-limit", default="96GB")
    p.add_argument("--temp-dir", default="/data/tmp/duckdb_pln")
    return build(p.parse_args(argv))


if __name__ == "__main__":
    raise SystemExit(main())
