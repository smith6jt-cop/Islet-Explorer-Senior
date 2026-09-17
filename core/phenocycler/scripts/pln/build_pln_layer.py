#!/usr/bin/env python3
"""Publish the pLN node-membership layer artifact, and reconcile it against QuPath.

WHAT THIS IS
------------
`scripts/senior/derive_pln_node_membership.py` does the measurement: it tests every pLN cell's
`Centroid X µm` / `Centroid Y µm` against the node annotation ROIs exported from the QuPath
`.qpdata`.  This builder does no geometry of its own except for the handful of cells that fell
in no node.  It takes that development-lane output and publishes it as the layer artifact
`docs/ISLET_LAYERS_SPEC.md` section 1.7 names, adds the node-level table, and writes the
MANIFEST that records what each number was checked against.

    derive_pln_node_membership.py  ->  data/development/pln/        (measurement, development)
    build_pln_layer.py             ->  /data/islet_layers/pln/      (artifact, what a card reads)

THE CROSS-CHECK, AND WHY IT IS NOT THE SOURCE
---------------------------------------------
QuPath's own hierarchy records, per annotation, how many detections it holds (`nChildren`).
On these slides the cells are children of the `Tissue` outlines, never of the `pLN_*`
annotations -- every `pLN_*` line in the audit log reads `nChildren=0`, which is also why the
export's `Parent` column says `Annotation (Tissue)` on 100% of rows and can never give node
membership.  Each node region contains exactly one `Tissue` outline, so QuPath's per-`Tissue`
child count is an independent count of the same cells.

It is a **cross-check, never the source**.  Membership is the geometric test; a disagreement is
reported with its size, never reconciled by adopting QuPath's number.

The audit log is the stdout of
`core/phenocycler/scripts/groovy/export_pln_node_annotations.groovy`, captured per entry.  It is
not in this checkout -- it lives beside the artifact, under `<out>/provenance/`.  Regenerate it
with the roster's own entry list (`qupath_entry` column, one per image):

    G=core/phenocycler/scripts/groovy/export_pln_node_annotations.groovy
    for e in 1 2 4 5 6 7 9 10 11 12 13 14 15 16 17 18 19 20 21; do
      echo "ENTRY=$e"
      JAVA_TOOL_OPTIONS="-Xmx48g" /home/smith6jt/QuPath/bin/QuPath script \
        --args /data/IO60pLN/data/$e/data.qpdata \
        --args /data/IO60pLN/results/node_geojson/entry$e.geojson "$G" \
        | grep -E '^(ANNOT|READ_MS)'
    done > "${ISLET_LAYERS_DIR:-/data/islet_layers}/pln/provenance/qupath_annotation_audit.log"

Without it every node's `qupath_tissue_child_count` is written as unavailable with that reason.
It is never defaulted to zero, and a node is never called "reproduced" on a check that did not run.

The MANIFEST counts what that log holds rather than repeating a number: annotations are the
`ANNOT` records alone, and the count is cross-checked against the exporter's own
`N_ANNOTATIONS` per entry and against the geojson feature counts.  Adding the one `READ_MS` line
per image to the `ANNOT` records gives a LINE count, not an annotation count -- 90 annotations
over 19 images make 109 such lines -- and reporting that sum as the number of annotations is a
mistake the derived count now makes impossible.

ABSENCE
-------
Three of donor 6534's four recorded nodes are `PolylineROI` with zero area.  A line bounds no
region, so no cell can be tested against it: those nodes are written into `nodes.parquet` with
`state='unavailable'` and the measured reason, keeping their `node_index` from the roster.  They
are not 0-cell nodes and they are not absent.

A cell that fell in no node keeps its row with `node_name = NULL` and a reason.  It is never
dropped and never labelled with a node it is not in.

NODE SCOPE -- WHY A PER-NODE NUMBER IS NOT ALWAYS A PER-NODE NUMBER
-------------------------------------------------------------------
On some images one annotation is the whole slide wearing a node's name.  6534 is the clearest
case: its one surviving node, `pLN_bottom`, outlines 10,742,163 um2 against a `Tissue` outline of
9,863,458 um2 and holds all 147,959 of the image's cells -- the "node" IS the slide.  Pooling
that with the two- and three-node images would pool a per-slide number with per-node ones.

So every measured node carries a `node_scope`, measured on each build from three facts and keyed
on no donor id: it is the only usable node on its image, it holds every cell assigned on that
image, and its annotation area is at least the total `Tissue` area on that image.  All three ->
`slide`; otherwise `node`.  The scope is written PER ROW in the membership file as well as in
`nodes.parquet`, so a consumer that pools nodes across donors cannot do it unknowingly.

EXCLUSIONS
----------
Two authorities have to agree, exactly as `ingest_pln_cells.py` does it: `phenocycler.cohort`
(asked, never restated) and the tracked `pln_node_annotations.csv` `excluded` flag.  A donor
either one withholds is not written, and the disagreement is named in the MANIFEST.

That second flag is for withholding ONE IMAGE for a reason about that image, and it must say
which in `excluded_reason`.  It is not a place to restate a central cohort exclusion, and this
builder refuses a row that tries: an `excluded=True` with no reason, or on a donor
`phenocycler.cohort` has REINSTATED, stops the build.  6579 carried exactly that stale echo --
excluded=True on the central pancreatitis-outlier exclusion, which cohort.py lifted on
2026-09-17 -- and withheld 1,231,818 cells from a decision that had lapsed.

A PARTIAL RUN MAY NOT SHRINK THE PUBLISHED ROOT
-----------------------------------------------
The shared rule is in `partial_run_guard.py`.  This lane takes its REFUSE branch: a donor the
output root already holds and this run does not write stops the build when the reason is a
MISSING INPUT (its derivation has not been run), and is removed as before when the reason is a
DECISION (`phenocycler.cohort`, or the roster's `excluded` flag).  There is no such thing as a
one-donor node-membership layer, so there is no scoped copy worth writing -- build a partial
run into its own root with `--out` instead, which leaves the published root untouched.

THE OUTPUT ROOT IS RECONCILED, NOT APPENDED TO
-----------------------------------------------
A rebuild into a populated root used to leave behind the partition of any donor now withheld or
unavailable, while the MANIFEST reported the smaller total -- the MANIFEST reconciled against
this builder's own in-memory counter rather than against what was on disk.  Now the run removes
every `donor_id=*/pln_node_membership.parquet` it did not itself write (refusing to touch any
other file it finds there), and then MEASURES the totals back off the files with
`pq.ParquetFile(...).metadata.num_rows`.  The in-memory counter is still recorded, next to the
measured number and a boolean saying whether they agree.

Usage
-----
    .venv/bin/python core/phenocycler/scripts/pln/build_pln_layer.py
    .venv/bin/python core/phenocycler/scripts/pln/build_pln_layer.py --out /tmp/probe_layer
"""
from __future__ import annotations

import argparse
import collections
import glob as globmod
import json
import os
import pathlib
import re
import sys
import time

import numpy as np
import pandas as pd
import pyarrow as pa
import pyarrow.dataset as ds
import pyarrow.parquet as pq

REPO = pathlib.Path(__file__).resolve().parents[4]
sys.path.insert(0, str(REPO / "core" / "phenocycler"))
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))

from partial_run_guard import refuse_shrink  # noqa: E402

from phenocycler.cohort import (  # noqa: E402
    DONOR_EXCLUSIONS,
    REINSTATED,
    filter_eligible_donors,
)

ROSTER = REPO / "core/phenocycler/data/integration/pln_node_annotations.csv"
DERIVED = REPO / "data/development/pln"
QP_PROJECT = pathlib.Path("/data/IO60pLN")
GEOJSON_DIR = QP_PROJECT / "results/node_geojson"
INGEST_REPORT = REPO / "core/phenocycler/data/integration/cells_pln_ingest_report.json"
CELLS_PLN = REPO / "core/phenocycler/data/integration/cells_pln"

#: Spec 0.5: the default layers root is $ISLET_LAYERS_DIR, /data/islet_layers, and this lane must
#: use the same default as its readers.  It briefly defaulted to the in-repo `data/islet_layers`
#: instead, on a mistaken instruction, and the artifact was then invisible: every consumer --
#: ie_layer_io.R::ie_layers_dir and cards x1/x2/s1/s2 among them -- resolves /data/islet_layers,
#: where the other six layers already sit.  A layer nothing can read is the defect this whole
#: lane exists to fix, so the default is pinned here rather than left to agree by accident.
DEFAULT_OUT = pathlib.Path(
    os.environ.get("ISLET_LAYERS_DIR", "/data/islet_layers")
) / "pln"
#: The provenance log lives beside the artifact, so --audit-log follows --out rather than
#: being pinned here: a probe build into another root must not read the published provenance.
AUDIT_LOG_NAME = "provenance/qupath_annotation_audit.log"

#: Spec 0.5.  kind:donor:roi; every pLN image carries one roi and the roster records it.
FRAME_KIND = "phenocycler_um"

UNASSIGNED_REASON = "centroid_outside_every_node_annotation_on_this_image"
NO_AREA_REASON = "annotation_is_{geom}_with_zero_area_it_bounds_no_region"

#: The only file name this builder writes into a `donor_id=` partition, and so the only one it
#: will remove when reconciling the root.  Anything else there stops the build.
PARTITION_FILE = "pln_node_membership.parquet"


def norm_donor(value: str) -> str:
    """'0112' -> '112'.  Spec 0.3.  Never int(): a float column would ruin 6533."""
    return str(value).strip().lstrip("0") or "0"


def read_roster(path: pathlib.Path = ROSTER) -> pd.DataFrame:
    df = pd.read_csv(path, comment="#", dtype={"donor_id": str, "qupath_entry": int})
    df["donor_id"] = df["donor_id"].map(norm_donor)
    df["node_classes"] = df["node_classes"].str.split(";")
    df["excluded"] = df["excluded"].astype(str).str.strip().str.lower() == "true"
    if "excluded_reason" not in df.columns:
        df["excluded_reason"] = ""
    df["excluded_reason"] = df["excluded_reason"].fillna("").astype(str).str.strip()
    return df


def check_roster_exclusions(roster: pd.DataFrame, path: pathlib.Path) -> None:
    """Refuse an `excluded=True` that says nothing, or one a lapsed central ruling left behind.

    The central exclusions live once, in `phenocycler.cohort`.  A row here that merely echoed
    them went stale the moment the maintainer lifted one -- which is how 6579 kept withholding
    1,231,818 cells after 2026-09-17.  Stopping the build is how this file learns the truth.
    """
    for _, row in roster[roster["excluded"]].iterrows():
        donor = row["donor_id"]
        if not row["excluded_reason"]:
            raise SystemExit(
                f"[err] {path}: {donor} is excluded=True with an empty excluded_reason. "
                "Say what about THIS IMAGE withholds it; a central cohort exclusion belongs in "
                "core/phenocycler/phenocycler/cohort.py and is asked for there, never restated "
                "here."
            )
        if donor in REINSTATED:
            raise SystemExit(
                f"[err] {path}: {donor} is excluded=True, but phenocycler.cohort has REINSTATED "
                f"it -- {REINSTATED[donor]}. The flag has lapsed. Set excluded=False, or give an "
                "excluded_reason about this image that the reinstatement does not answer."
            )
        if donor in DONOR_EXCLUSIONS:
            print(f"[info] {donor}: excluded here AND by phenocycler.cohort -- "
                  f"{row['excluded_reason']}")


def parse_audit_log(path: pathlib.Path) -> dict[str, dict]:
    """annotation uuid -> {n_children, roi, n_points, area_px}, from the groovy stdout.

    Absent file -> empty mapping.  The caller then writes the check as unavailable; it never
    substitutes a zero.
    """
    if not path.exists():
        return {}
    out: dict[str, dict] = {}
    field = re.compile(r"(\w+)=(.*)")
    for line in path.read_text().splitlines():
        if not line.startswith("ANNOT\t"):
            continue
        kv = {}
        for part in line.split("\t")[1:]:
            m = field.match(part)
            if m:
                kv[m.group(1)] = m.group(2)
        out[kv["id"]] = {
            "n_children": int(kv["nChildren"]),
            "roi": kv.get("roi"),
            "n_points": int(kv.get("nPoints", -1)),
            "area_px": float(kv.get("areaPx", "nan")),
            "classification": kv.get("class"),
        }
    return out


def measure_audit_log(path: pathlib.Path) -> dict:
    """Count what the provenance log actually holds, rather than transcribing a number.

    The log interleaves three line kinds: one `ENTRY=` echo and one `READ_MS=` summary per image,
    and one `ANNOT` record per annotation.  The annotation count is the ANNOT records alone --
    `grep -c '^ANNOT'`.  Adding the READ_MS lines gives a LINE count, not an annotation count,
    and a log of 90 annotations across 19 images has 109 such lines; recording that sum as the
    number of annotations is the transcription error this function exists to make impossible.

    Each `READ_MS=` line carries the exporter's own `N_ANNOTATIONS` for that image, so the count
    is checked against the exporter per entry and every disagreement is reported.
    """
    if not path.exists():
        return {"state": "unavailable(audit_log_absent)", "path": str(path)}
    lines = path.read_text().splitlines()
    declared: dict[str, int] = {}
    counted: dict[str, int] = collections.defaultdict(int)
    entry = None
    order: list[str] = []
    for line in lines:
        if line.startswith("ENTRY="):
            entry = line.split("=", 1)[1].strip()
            order.append(entry)
            counted.setdefault(entry, 0)
        elif line.startswith("READ_MS="):
            m = re.search(r"N_ANNOTATIONS=(\d+)", line)
            if m and entry is not None:
                declared[entry] = int(m.group(1))
        elif line.startswith("ANNOT\t") and entry is not None:
            counted[entry] += 1
    n_annot = sum(1 for line in lines if line.startswith("ANNOT\t"))
    n_read_ms = sum(1 for line in lines if line.startswith("READ_MS="))
    mismatch = [
        {"entry": e, "exporter_N_ANNOTATIONS": declared.get(e), "annot_records_counted": counted[e]}
        for e in order if declared.get(e) != counted[e]
    ]
    return {
        "state": "measured",
        "path": str(path),
        "how_counted": "len([line for line in log if line.startswith('ANNOT\\t')])",
        "annotation_records": n_annot,
        "entries": len(order),
        "read_ms_lines": n_read_ms,
        "lines_total": len(lines),
        "annot_plus_read_ms_lines": n_annot + n_read_ms,
        "not_the_annotation_count": (
            f"{n_annot + n_read_ms} is ANNOT + READ_MS lines, a line count. The annotations "
            f"number {n_annot}"
        ),
        "exporter_declared_total": sum(declared.values()) if declared else None,
        "entries_declaring_a_different_count": mismatch,
        "agrees_with_exporter": not mismatch and len(declared) == len(order),
    }


def measure_geojson_annotations(entries: list[int], geojson_dir: pathlib.Path) -> dict:
    """Feature counts of the geojson exports, as a second count of the same annotations."""
    per_entry, missing = {}, []
    for entry in entries:
        path = geojson_dir / f"entry{entry}.geojson"
        if not path.exists():
            missing.append(entry)
            continue
        gj = json.loads(path.read_text())
        feats = gj["features"] if isinstance(gj, dict) and gj.get("type") == "FeatureCollection" \
            else gj
        per_entry[str(entry)] = len(feats)
    return {
        "features_per_entry": per_entry,
        "features_total": sum(per_entry.values()),
        "entries_with_no_geojson": missing,
    }


def reconcile_out_root(out_root: pathlib.Path, written: set[str]) -> list[dict]:
    """Remove the partitions this run did not write, so the root IS what the MANIFEST says.

    Only `donor_id=<d>/pln_node_membership.parquet` is ever removed, and only for a donor this
    run did not write.  Any other file in a `donor_id=` directory stops the build rather than
    being deleted: this function must never be a way to lose something it did not create.
    """
    removed: list[dict] = []
    for part in sorted(out_root.glob("donor_id=*")):
        if not part.is_dir():
            continue
        donor = part.name.split("=", 1)[1]
        if donor in written:
            continue
        for found in sorted(part.iterdir()):
            if found.name != PARTITION_FILE:
                raise SystemExit(
                    f"[err] {found}: a stale partition holds a file this builder never writes. "
                    "Refusing to remove it -- move it aside by hand and rerun."
                )
            rows = pq.ParquetFile(found).metadata.num_rows
            found.unlink()
            removed.append({"donor_id": donor, "path": str(found), "rows_removed": int(rows)})
            print(f"[info] {donor}: removed a stale partition of {rows:,} rows -- this run does "
                  "not write this donor")
        part.rmdir()
    return removed


def measure_partitions(out_root: pathlib.Path) -> list[dict]:
    """Rows per donor, read back off the files.  This, not a counter, is the artifact's total."""
    out = []
    for path in sorted(out_root.glob(f"donor_id=*/{PARTITION_FILE}")):
        out.append(
            {
                "donor_id": path.parent.name.split("=", 1)[1],
                "path": str(path.relative_to(out_root)),
                "rows": int(pq.ParquetFile(path).metadata.num_rows),
                "bytes": int(path.stat().st_size),
            }
        )
    return out


def probe_read_routes(out_root: pathlib.Path) -> dict:
    """Run the three candidate read routes and record what each one actually did, this build.

    These used to be sentences with a row count baked into them; a rebuild with a different
    cohort re-emitted the old number unchanged.  Now each route is attempted here and the
    outcome -- the count, or the exception type and message -- is whatever happened.
    """
    routes: dict[str, str] = {}
    files = sorted(globmod.glob(str(out_root / "donor_id=*" / "*.parquet")))
    if not files:
        return {
            "state": "unavailable(no_partition_file_in_this_root)",
            "measured_on": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
            "measured_files": 0,
        }

    def attempt(label: str, fn) -> None:
        try:
            routes[label] = fn()
        except Exception as exc:                                   # noqa: BLE001 -- recorded
            first = str(exc).splitlines()[0] if str(exc) else ""
            routes[label] = f"FAILS {type(exc).__name__}: {first[:220]}"

    def by_glob() -> str:
        dataset = ds.dataset(files)
        return (f"OK, {dataset.count_rows():,} rows, donor_id "
                f"{dataset.schema.field('donor_id').type} -- this is the route to use")

    attempt("ds.dataset(glob.glob(root + '/donor_id=*/*.parquet'))", by_glob)
    attempt("ds.dataset(root)",
            lambda: f"OK, {ds.dataset(str(out_root)).count_rows():,} rows")
    attempt(
        "ds.dataset(root, ignore_prefixes=[...], partitioning='hive')",
        lambda: "OK, {:,} rows".format(
            ds.dataset(
                str(out_root), partitioning="hive",
                ignore_prefixes=["MANIFEST", "nodes", "unassigned", "provenance", ".", "_"],
            ).count_rows()
        ),
    )
    routes["measured_on"] = time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())
    routes["measured_files"] = len(files)
    return routes


def geojson_annotations(entry: int, geojson_dir: pathlib.Path = GEOJSON_DIR) -> list[dict]:
    """Every annotation on one image, in the order derive_pln_node_membership.py reads them."""
    from shapely.geometry import shape

    gj = json.loads((geojson_dir / f"entry{entry}.geojson").read_text())
    feats = gj["features"] if gj.get("type") == "FeatureCollection" else gj
    out = []
    for f in feats:
        props = f.get("properties") or {}
        geom = shape(f["geometry"])
        out.append(
            {
                "annotation_id": f.get("id"),
                "classification": (props.get("classification") or {}).get("name"),
                "geom_type": f["geometry"]["type"],
                "geom": geom,
                "is_area": geom.geom_type in ("Polygon", "MultiPolygon") and geom.area > 0,
            }
        )
    return out


def pixel_size_um(entry: int, project: pathlib.Path = QP_PROJECT) -> float:
    meta = json.loads((project / f"data/{entry}/server.json").read_text())
    cal = meta["metadata"]["pixelCalibration"]
    return float(cal["pixelWidth"]["value"])


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--derived", default=str(DERIVED),
                    help="output root of derive_pln_node_membership.py")
    ap.add_argument("--out", default=str(DEFAULT_OUT), help="layer artifact root")
    ap.add_argument("--audit-log", default=None,
                    help="captured stdout of export_pln_node_annotations.groovy; the default "
                         "follows --out, so a probe build into another root looks for its own "
                         "provenance rather than the published one")
    ap.add_argument("--roster", default=str(ROSTER), help="pln_node_annotations.csv")
    ap.add_argument("--geojson-dir", default=str(GEOJSON_DIR),
                    help="directory of entry<N>.geojson exports")
    ap.add_argument("--qupath-project", default=str(QP_PROJECT),
                    help="QuPath project root, for data/<entry>/server.json")
    ap.add_argument("--ingest-report", default=str(INGEST_REPORT),
                    help="cells_pln_ingest_report.json, for the no-measurement images")
    ap.add_argument("--cells-pln", default=str(CELLS_PLN),
                    help="the ingested pLN cells root, for the donor-coverage cross-check")
    args = ap.parse_args(argv)

    derived = pathlib.Path(args.derived)
    out_root = pathlib.Path(args.out)
    audit_log = (pathlib.Path(args.audit_log) if args.audit_log
                 else out_root / AUDIT_LOG_NAME)
    geojson_dir = pathlib.Path(args.geojson_dir)
    project = pathlib.Path(args.qupath_project)
    ingest_report = pathlib.Path(args.ingest_report)
    cells_pln = pathlib.Path(args.cells_pln)
    out_root.mkdir(parents=True, exist_ok=True)

    # What is in the root BEFORE this run touches it.  A layer must not gain or lose a donor
    # silently, so the before state is measured rather than assumed empty.
    before_partitions = measure_partitions(out_root)
    before_manifest_totals = None
    if (out_root / "MANIFEST.json").exists():
        try:
            before_manifest_totals = json.loads(
                (out_root / "MANIFEST.json").read_text())["totals"]
        except (json.JSONDecodeError, KeyError) as exc:                # noqa: BLE001 -- recorded
            before_manifest_totals = {"state": f"unreadable({type(exc).__name__})"}

    roster = read_roster(pathlib.Path(args.roster))
    check_roster_exclusions(roster, pathlib.Path(args.roster))
    audit_geom = pd.read_csv(derived / "node_geometry_audit.csv", dtype={"donor_id": str})
    audit_geom["donor_id"] = audit_geom["donor_id"].map(norm_donor)
    qp = parse_audit_log(audit_log)
    if not qp:
        print(f"[warn] no QuPath annotation audit log at {audit_log}; every child-count check "
              "is written unavailable")

    # -- who may be written.  Both authorities, neither restated. ------------------------
    all_donors = sorted(set(roster["donor_id"]))
    by_cohort = set(filter_eligible_donors(all_donors))
    by_manifest = set(roster.loc[~roster["excluded"], "donor_id"])
    disagreement = sorted(by_cohort ^ by_manifest)
    eligible = sorted(by_cohort & by_manifest)
    roster_reason = dict(zip(roster["donor_id"], roster["excluded_reason"]))
    withheld = {}
    for donor in sorted(set(all_donors) - set(eligible)):
        if donor in DONOR_EXCLUSIONS:
            withheld[donor] = f"phenocycler.cohort: {DONOR_EXCLUSIONS[donor]}"
        else:
            withheld[donor] = (
                "pln_node_annotations.csv excluded=True: "
                f"{roster_reason.get(donor) or 'no excluded_reason recorded'}"
            )
        print(f"[info] {donor}: withheld -- {withheld[donor]}")

    cell_rows_counted_in_memory = 0
    donors_written: set[str] = set()
    donors_manifest: list[dict] = []
    node_rows: list[dict] = []
    unassigned_detail: list[dict] = []

    for _, row in roster.iterrows():
        donor, entry = row["donor_id"], int(row["qupath_entry"])
        image = row["image_file"]
        recorded = list(row["node_classes"])
        src = derived / "node_membership" / f"pln_node_membership_{donor}.parquet"

        # A donor this lane cannot write is RECORDED, with its state and the reason.  It is
        # never a 0-cell donor and it is never simply missing from the artifact: a card that
        # reads 0 cells for a donor whose export never arrived would report a measurement.
        state, reason = "measured", ""
        if donor in withheld:
            state, reason = "withheld", withheld[donor]
        elif not src.exists():
            state = "unavailable"
            reason = (f"no membership artifact under {derived}/node_membership -- "
                      "scripts/senior/derive_pln_node_membership.py has not been run for "
                      "this image")
        if state != "measured":
            print(f"[info] {donor}: {state} -- {reason}")
            donors_manifest.append(
                {
                    "donor_id": donor, "image": image, "qupath_entry": entry,
                    "roi": row["roi"], "frame": None, "state": state, "reason": reason,
                    "n_nodes_recorded": int(row["n_nodes"]),
                    "n_nodes_usable": None, "n_nodes_unavailable": None,
                    "cells_assigned": None, "cells_unassigned": None, "cells_total": None,
                    "roster_n_cells": int(row["n_cells"]),
                    "reconciles_with_roster": None,
                    "qupath_tissue_child_count_total": None,
                    "nodes_checked_against_qupath": 0,
                    "nodes_reproducing_qupath_exactly": 0,
                }
            )
            for name in recorded:
                node_rows.append(
                    {
                        "donor_id": donor, "image": image, "qupath_entry": entry,
                        "node_index": recorded.index(name), "node_name": name,
                        "annotation_id": None, "annotation_geom": None,
                        "delimiter_kind": None, "area_um2": None,
                        "state": state, "reason": reason, "n_cells": None,
                        "image_tissue_area_um2": None,
                        "area_ratio_node_over_image_tissue": None,
                        "n_usable_nodes_on_image": None,
                        "holds_every_assigned_cell_on_image": None,
                        "node_scope": None,
                        "node_scope_reason": f"not_applicable_donor_is_{state}",
                        "qupath_tissue_annotation_id": None,
                        "qupath_tissue_child_count": None,
                        "qupath_check": f"not_applicable_donor_is_{state}",
                        "delta_vs_qupath": None,
                    }
                )
            continue

        annots = geojson_annotations(entry, geojson_dir)
        # derive_pln_node_membership.py indexes its tissue list this way; keep it identical or
        # tissue_annotation_index points at the wrong annotation.
        tissue_list = [a for a in annots if a["classification"] == "Tissue" and a["is_area"]]
        px_um = pixel_size_um(entry, project)

        table = pq.read_table(src)
        node_id = np.asarray(table.column("node_id").to_pylist(), dtype=object)
        assigned = node_id != "unassigned"
        index_of = {name: i for i, name in enumerate(recorded)}
        unknown = sorted(set(node_id[assigned]) - set(index_of))
        if unknown:
            raise SystemExit(f"[err] {donor}: membership names nodes the roster does not "
                             f"record: {unknown}")

        node_index = np.array(
            [index_of[n] if a else None for n, a in zip(node_id, assigned)], dtype=object
        )
        n_cells = table.num_rows
        n_assigned = int(assigned.sum())
        n_unassigned = n_cells - n_assigned
        object_id = table.column("object_id").to_pylist()
        if len(set(object_id)) != n_cells:
            raise SystemExit(f"[err] {donor}: object_id is not unique within the donor "
                             f"({len(set(object_id)):,} distinct of {n_cells:,})")

        frame = f"{FRAME_KIND}:{donor}:{row['roi']}"

        # -- the cells that fell in no node: how far outside, measured -------------------
        if n_unassigned:
            from shapely.geometry import Point

            node_geoms = {
                a["classification"]: a["geom"] for a in annots
                if a["classification"] in recorded and a["is_area"]
            }
            xs = table.column("x_um").to_pylist()
            ys = table.column("y_um").to_pylist()
            for i in np.flatnonzero(~assigned):
                p = Point(xs[i] / px_um, ys[i] / px_um)
                near, dist = min(
                    ((name, g.distance(p) * px_um) for name, g in node_geoms.items()),
                    key=lambda kv: kv[1],
                )
                unassigned_detail.append(
                    {
                        "donor_id": donor, "image": image, "object_id": object_id[i],
                        "x_um": round(float(xs[i]), 3), "y_um": round(float(ys[i]), 3),
                        "nearest_node": near,
                        "distance_outside_um": round(float(dist), 4),
                        "reason": UNASSIGNED_REASON,
                    }
                )

        # -- node-level table, including the nodes that bound nothing --------------------
        per_node_counts = collections.Counter(node_id[assigned].tolist())
        tissue_idx = np.asarray(table.column("tissue_annotation_index").to_pylist())
        donor_checked = donor_matched = 0

        # Scope: is this a node, or the slide wearing a node's name?  Measured per image, from
        # the image's own geometry -- never keyed on a donor id.
        image_geom = audit_geom[audit_geom["donor_id"] == donor]
        tissue_area_um2 = float(
            image_geom.loc[(image_geom["classification"] == "Tissue")
                           & (image_geom["area_um2"] > 0), "area_um2"].sum()
        )
        n_usable_nodes = int((image_geom["role"] == "node").sum())
        donor_node_scope: dict[str, str] = {}

        for name in recorded:
            geom_row = audit_geom[(audit_geom["donor_id"] == donor)
                                  & (audit_geom["classification"] == name)]
            if len(geom_row) != 1:
                raise SystemExit(f"[err] {donor}/{name}: {len(geom_row)} geometry audit rows")
            geom_row = geom_row.iloc[0]
            usable = bool(geom_row["role"] == "node")
            record = {
                "donor_id": donor, "image": image, "qupath_entry": entry,
                "node_index": index_of[name], "node_name": name,
                "annotation_id": geom_row["annotation_id"],
                "annotation_geom": geom_row["geom_type"],
                "delimiter_kind": ("rectangle" if bool(geom_row["is_axis_aligned_rectangle"])
                                   else "outline") if usable else "none",
                "area_um2": float(geom_row["area_um2"]) if usable else None,
                "state": "measured" if usable else "unavailable",
                "reason": "" if usable else NO_AREA_REASON.format(
                    geom=str(geom_row["geom_type"]).lower()),
                "n_cells": int(per_node_counts.get(name, 0)) if usable else None,
                "image_tissue_area_um2": tissue_area_um2 if usable else None,
                "area_ratio_node_over_image_tissue": None,
                "n_usable_nodes_on_image": n_usable_nodes if usable else None,
                "holds_every_assigned_cell_on_image": None,
                "node_scope": None,
                "node_scope_reason": "not_applicable_node_bounds_no_region",
                "qupath_tissue_annotation_id": None,
                "qupath_tissue_child_count": None,
                "qupath_check": "not_applicable_node_bounds_no_region",
                "delta_vs_qupath": None,
            }
            if usable:
                # -- node or slide?  Three measured facts, no threshold chosen by hand.
                holds_all = bool(record["n_cells"] == n_assigned)
                ratio = (float(record["area_um2"]) / tissue_area_um2
                         if tissue_area_um2 > 0 else None)
                record["holds_every_assigned_cell_on_image"] = holds_all
                record["area_ratio_node_over_image_tissue"] = ratio
                only_node = n_usable_nodes == 1
                covers = ratio is not None and ratio >= 1.0
                if only_node and holds_all and covers:
                    record["node_scope"] = "slide"
                    record["node_scope_reason"] = (
                        f"the only usable node on this image, holds all {n_assigned:,} of its "
                        f"assigned cells, and outlines {record['area_um2']:,.0f} um2 against "
                        f"{tissue_area_um2:,.0f} um2 of Tissue (ratio {ratio:.3f}) -- this is a "
                        "per-slide region wearing a node's name, and a per-node number for it is "
                        "a per-slide number"
                    )
                else:
                    record["node_scope"] = "node"
                    record["node_scope_reason"] = (
                        f"one of {n_usable_nodes} usable node(s) on this image; holds "
                        f"{record['n_cells']:,} of {n_assigned:,} assigned cells; area ratio to "
                        f"the image's Tissue outlines "
                        + ("unavailable" if ratio is None else f"{ratio:.3f}")
                    )
                donor_node_scope[name] = record["node_scope"]

                # which Tissue outline sits in this node: the one this node's own cells are in.
                in_node = tissue_idx[(node_id == name) & assigned]
                seen = collections.Counter(int(v) for v in in_node if v >= 0)
                if len(seen) == 1:
                    tissue = tissue_list[next(iter(seen))]
                    record["qupath_tissue_annotation_id"] = tissue["annotation_id"]
                    hit = qp.get(tissue["annotation_id"])
                    if hit is None:
                        record["qupath_check"] = (
                            "unavailable(qupath_annotation_audit_log_absent)" if not qp
                            else "unavailable(annotation_absent_from_audit_log)"
                        )
                    else:
                        record["qupath_tissue_child_count"] = int(hit["n_children"])
                        delta = record["n_cells"] - int(hit["n_children"])
                        record["delta_vs_qupath"] = int(delta)
                        record["qupath_check"] = "reproduced" if delta == 0 else "differs"
                        donor_checked += 1
                        donor_matched += int(delta == 0)
                else:
                    record["qupath_check"] = (
                        f"unavailable(node_spans_{len(seen)}_tissue_annotations)"
                    )
            node_rows.append(record)

        # -- the membership file, written now that every node's scope is measured ---------
        scopeless = sorted(set(node_id[assigned].tolist()) - set(donor_node_scope))
        if scopeless:
            raise SystemExit(f"[err] {donor}: cells are assigned to node(s) with no measured "
                             f"scope: {scopeless}")
        out_dir = out_root / f"donor_id={donor}"
        out_dir.mkdir(parents=True, exist_ok=True)
        pq.write_table(
            pa.table(
                {
                    "object_id": pa.array(object_id, type=pa.string()),
                    "donor_id": pa.array([donor] * n_cells, type=pa.string()),
                    "image": pa.array([image] * n_cells, type=pa.string()),
                    "node_name": pa.array(
                        [n if a else None for n, a in zip(node_id, assigned)], type=pa.string()
                    ),
                    "node_index": pa.array(node_index.tolist(), type=pa.int32()),
                    "node_state": pa.array(
                        ["measured" if a else "unassigned" for a in assigned], type=pa.string()
                    ),
                    # Carried per row on purpose: a consumer that pools nodes across donors sees
                    # the slide-scope ones here, not only in nodes.parquet or in the MANIFEST.
                    "node_scope": pa.array(
                        [donor_node_scope[n] if a else None for n, a in zip(node_id, assigned)],
                        type=pa.string(),
                    ),
                    "reason": pa.array(
                        ["" if a else UNASSIGNED_REASON for a in assigned], type=pa.string()
                    ),
                    "x_um": table.column("x_um").cast(pa.float64()),
                    "y_um": table.column("y_um").cast(pa.float64()),
                    "frame": pa.array([frame] * n_cells, type=pa.string()),
                }
            ),
            out_dir / PARTITION_FILE,
            compression="zstd",
        )
        cell_rows_counted_in_memory += n_cells
        donors_written.add(donor)

        donors_manifest.append(
            {
                "donor_id": donor,
                "image": image,
                "qupath_entry": entry,
                "roi": row["roi"],
                "frame": frame,
                "state": "measured",
                "reason": "",
                "n_nodes_recorded": int(row["n_nodes"]),
                "n_nodes_usable": int(sum(1 for r in node_rows
                                          if r["donor_id"] == donor and r["state"] == "measured")),
                "n_nodes_unavailable": int(sum(1 for r in node_rows
                                               if r["donor_id"] == donor
                                               and r["state"] == "unavailable")),
                "cells_assigned": n_assigned,
                "cells_unassigned": n_unassigned,
                "cells_total": n_cells,
                "roster_n_cells": int(row["n_cells"]),
                "reconciles_with_roster": bool(n_cells == int(row["n_cells"])),
                "qupath_tissue_child_count_total": (
                    int(sum(r["qupath_tissue_child_count"] for r in node_rows
                            if r["donor_id"] == donor
                            and r["qupath_tissue_child_count"] is not None))
                    if donor_checked else None
                ),
                "nodes_checked_against_qupath": donor_checked,
                "nodes_reproducing_qupath_exactly": donor_matched,
            }
        )
        print(f"  {donor:>5} {image}: {n_cells:>9,} cells, {n_assigned:>9,} in a node, "
              f"{n_unassigned} unassigned, {donor_matched}/{donor_checked} nodes reproduce QuPath")

    # -- a partial run may not shrink a published artifact --------------------------------
    # The rule is `partial_run_guard.py`'s, and this lane takes its REFUSE branch: there is no
    # such thing as a one-donor node-membership layer, so a scoped copy would be a thing no
    # reader wants.  A donor the root already holds and this run does not write is allowed to
    # go ONLY when a decision withheld it -- phenocycler.cohort, or the tracked roster's
    # `excluded` flag.  A donor missing because its derivation has not been run is an absent
    # input, and removing it would republish a smaller layer under the same name with nothing
    # saying so.  `build_pln_immune_pheno.py --donors 6539 --stage standout` did exactly that
    # to three cohort files on 2026-09-17; this is the same defect wearing this lane's clothes.
    state_by_donor = {d["donor_id"]: (d["state"], d["reason"]) for d in donors_manifest}
    held_before = {p["donor_id"] for p in before_partitions}
    would_lose_an_input = {}
    for donor in sorted(held_before - donors_written):
        state, reason = state_by_donor.get(donor, ("absent_from_the_roster",
                                                   "this run's roster has no row for it"))
        if state in ("withheld",):
            continue                                   # a decision; the reconciliation says so
        would_lose_an_input[donor] = f"{state}: {reason}"
    refuse_shrink(
        artifact=f"the pLN node-membership layer at {out_root}",
        unit_kind="donor",
        losing=would_lose_an_input,
        hint=("run scripts/senior/derive_pln_node_membership.py for those donors and build "
              f"again, or build this partial run into its own root with --out (the published "
              f"root {out_root} is then untouched)."),
    )

    # -- the root is reconciled against this run, then MEASURED back off the files ---------
    stale_removed = reconcile_out_root(out_root, donors_written)
    partitions = measure_partitions(out_root)
    cells_on_disk = int(sum(p["rows"] for p in partitions))
    donors_on_disk = sorted(p["donor_id"] for p in partitions)

    measured_donors = [d for d in donors_manifest if d["state"] == "measured"]
    roster_measured_cells = int(sum(d["roster_n_cells"] for d in measured_donors))
    rows_by_donor = {p["donor_id"]: p["rows"] for p in partitions}
    for d in measured_donors:
        d["rows_on_disk"] = rows_by_donor.get(d["donor_id"])
        d["file_matches_cells_total"] = bool(d["rows_on_disk"] == d["cells_total"])
    nodes = pd.DataFrame(node_rows)
    # Explicit types.  Through pandas, a None in an integer column becomes NaN in a float column
    # and an unavailable node starts reading like a measured 0.0.
    node_schema = pa.schema(
        [
            ("donor_id", pa.string()), ("image", pa.string()), ("qupath_entry", pa.int32()),
            ("node_index", pa.int32()), ("node_name", pa.string()),
            ("annotation_id", pa.string()), ("annotation_geom", pa.string()),
            ("delimiter_kind", pa.string()), ("area_um2", pa.float64()),
            ("state", pa.string()), ("reason", pa.string()), ("n_cells", pa.int64()),
            ("image_tissue_area_um2", pa.float64()),
            ("area_ratio_node_over_image_tissue", pa.float64()),
            ("n_usable_nodes_on_image", pa.int32()),
            ("holds_every_assigned_cell_on_image", pa.bool_()),
            ("node_scope", pa.string()), ("node_scope_reason", pa.string()),
            ("qupath_tissue_annotation_id", pa.string()),
            ("qupath_tissue_child_count", pa.int64()), ("qupath_check", pa.string()),
            ("delta_vs_qupath", pa.int64()),
        ]
    )
    pq.write_table(
        pa.table({f.name: pa.array([r[f.name] for r in node_rows], type=f.type)
                  for f in node_schema}),
        out_root / "nodes.parquet", compression="zstd",
    )
    if unassigned_detail:
        pd.DataFrame(unassigned_detail).to_csv(out_root / "unassigned_cells.csv", index=False)

    # -- what the export itself cannot deliver, read from the ingest's own measurement ----
    no_measurement = {"state": "unavailable(cells_pln_ingest_report_absent)", "images": []}
    if ingest_report.exists():
        rep = json.loads(ingest_report.read_text())
        no_measurement = {
            "measured_by": str(ingest_report),
            "measured_at_utc": rep["generated_at_utc"],
            "rule": rep["measurement_availability"]["rule"],
            "n_cells_with_no_measurements": rep["measurement_availability"][
                "n_cells_with_no_measurements"],
            "images": [
                {
                    "donor_id": norm_donor(e["donor_id"]),
                    "image": e["image"],
                    "n_rows": e["n_rows_in_export"],
                    "n_rows_with_no_measurement": e["n_no_measurements_exported"],
                    "state": e["measurements_state"],
                }
                for e in rep["per_image"] if e["measurements_state"].startswith("unavailable")
            ],
        }

    # -- does the ingested cell root carry the same donors this layer does? ---------------
    if cells_pln.exists():
        ingested = sorted(d.name.split("=", 1)[1] for d in cells_pln.glob("donor_id=*")
                          if d.is_dir())
        coverage = {
            "root": str(cells_pln),
            "donors_in_that_root": ingested,
            "in_this_layer_but_not_ingested": sorted(set(donors_on_disk) - set(ingested)),
            "ingested_but_not_in_this_layer": sorted(set(ingested) - set(donors_on_disk)),
        }
        coverage["agrees"] = not (coverage["in_this_layer_but_not_ingested"]
                                  or coverage["ingested_but_not_in_this_layer"])
        coverage["what_a_disagreement_means"] = (
            "a donor this layer carries but that root does not cannot be joined to any "
            "measurement: an object_id join drops it silently. Reported, not worked around -- "
            "re-running the ingest is core/phenocycler/scripts/pln/ingest_pln_cells.py's job, "
            "not this builder's"
        )
    else:
        coverage = {"state": "unavailable(cells_pln_root_absent)", "root": str(cells_pln)}

    checked = [r for r in node_rows if r["qupath_check"] in ("reproduced", "differs")]
    differs = [r for r in node_rows if r["qupath_check"] == "differs"]
    slide_scope = [r for r in node_rows if r["node_scope"] == "slide"]
    log_measured = measure_audit_log(audit_log)
    geojson_measured = measure_geojson_annotations(
        sorted(int(e) for e in roster["qupath_entry"]), geojson_dir)
    if unassigned_detail:
        distances = [u["distance_outside_um"] for u in unassigned_detail]
        unassigned_note = (
            f"{len(unassigned_detail)} cell(s) fell in no node on this build. Measured here: "
            f"each is between {min(distances):g} and {max(distances):g} um outside its nearest "
            "node boundary. Strict containment excludes a centroid on the boundary and QuPath's "
            "hierarchy does not. The difference is reported, not reconciled"
        )
    else:
        unassigned_note = (
            "no cell fell outside every node annotation on this build, over the "
            f"{cells_on_disk:,} rows measured off the files"
        )
    manifest = {
        "artifact": "pln node membership",
        "generated_at_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
        "built_by": "core/phenocycler/scripts/pln/build_pln_layer.py",
        "root": str(out_root),
        "layout": f"donor_id=<donor>/{PARTITION_FILE}, plus nodes.parquet",
        "how_to_read": {
            "rule": (
                "donor_id is BOTH the directory key and a string column in the file. Take it "
                "from the column, with the pyarrow default partitioning=None. Each route below "
                "was attempted against THIS root during THIS build; the outcomes are what "
                "happened, not a transcript of an earlier one."
            ),
            **probe_read_routes(out_root),
        },
        "measurement": {
            "by": "scripts/senior/derive_pln_node_membership.py",
            "consumed_from": str(derived),
            "method": (
                "point-in-polygon of Centroid X/Y um against the node annotation ROI exported "
                "from the QuPath .qpdata; strict containment (shapely 'within')"
            ),
            "object_id": (
                "the QuPath UUID, carried verbatim from the export; no centroid matching, no "
                "row order, no filename parsing"
            ),
            "never_from": (
                "the export's Parent column, which reads 'Annotation (Tissue)' on 100% of rows"
            ),
        },
        "totals": {
            "roster_images": int(len(roster)),
            "donors_measured": len(measured_donors),
            "donors_not_written": [
                {"donor_id": d["donor_id"], "state": d["state"], "reason": d["reason"]}
                for d in donors_manifest if d["state"] != "measured"
            ],
            #: MEASURED off the files with pq.ParquetFile(...).metadata.num_rows, after the root
            #: was reconciled.  Not a counter this process kept: a counter cannot see a stale
            #: partition left by an earlier build, and that is exactly how this total went wrong.
            "cells_written": cells_on_disk,
            "cells_written_measured_by": (
                "sum(pq.ParquetFile(f).metadata.num_rows for f in "
                f"root/donor_id=*/{PARTITION_FILE})"
            ),
            "cells_counted_in_memory_by_this_run": cell_rows_counted_in_memory,
            "files_match_the_in_memory_counter":
                bool(cells_on_disk == cell_rows_counted_in_memory),
            "cells_assigned_to_a_node": sum(d["cells_assigned"] for d in measured_donors),
            "cells_unassigned": sum(d["cells_unassigned"] for d in measured_donors),
            "roster_total_cells_for_measured_donors": roster_measured_cells,
            "reconciles_with_roster": bool(cells_on_disk == roster_measured_cells),
            "nodes_recorded_all_images": int(roster["n_nodes"].sum()),
            "nodes_recorded_for_measured_donors": int(
                sum(d["n_nodes_recorded"] for d in measured_donors)),
            "nodes_measured": int((nodes["state"] == "measured").sum()),
            "nodes_unavailable": int((nodes["state"] == "unavailable").sum()),
            "nodes_withheld": int((nodes["state"] == "withheld").sum()),
            "nodes_slide_scope": len(slide_scope),
            "nodes_node_scope": int((nodes["node_scope"] == "node").sum()),
        },
        "output_root_reconciliation": {
            "why": (
                "the builder used to append to whatever was already in the root, so a donor "
                "withheld since the last build kept its partition while the MANIFEST reported "
                "the smaller total. The MANIFEST then reconciled against this builder's "
                "bookkeeping rather than against the files"
            ),
            "rule": (
                f"after writing, every donor_id=*/{PARTITION_FILE} whose donor this run did not "
                "write is removed; any other file found in a donor_id= directory stops the build "
                "rather than being removed"
            ),
            "donors_written_this_run": sorted(donors_written),
            "donors_on_disk_after": donors_on_disk,
            "root_matches_this_run": bool(donors_on_disk == sorted(donors_written)),
            "stale_partitions_removed": stale_removed,
            "partitions": partitions,
            #: Before/after, both measured off the files.  A donor may be gained or lost, but
            #: never silently: the change is named here on every build.
            "before_this_run": {
                "donors": sorted(p["donor_id"] for p in before_partitions),
                "cells": int(sum(p["rows"] for p in before_partitions)),
                "previous_manifest_totals": before_manifest_totals,
            },
            "donors_gained": sorted(set(donors_on_disk)
                                    - {p["donor_id"] for p in before_partitions}),
            "donors_lost": sorted({p["donor_id"] for p in before_partitions}
                                  - set(donors_on_disk)),
            "cells_delta": int(cells_on_disk - sum(p["rows"] for p in before_partitions)),
            #: The partial-run rule (partial_run_guard.py), as this lane applies it.  Reaching
            #: the MANIFEST at all means no donor was about to be removed for want of an input.
            "partial_run_rule": {
                "rule": ("a run that would take a donor OUT of this published root because its "
                         "derivation is missing is refused; a donor removed by a DECISION -- "
                         "phenocycler.cohort, or the roster's excluded flag -- is removed and "
                         "named above"),
                "donors_the_root_held_before_this_run": sorted(
                    p["donor_id"] for p in before_partitions),
                "donors_removed_by_a_decision": sorted(
                    r["donor_id"] for r in stale_removed),
                "donors_removed_for_a_missing_input": [],
            },
        },
        "node_scope": {
            "what": (
                "on some images one annotation is the slide wearing a node's name, so a "
                "per-node number for it is a per-slide number. Pooling nodes across donors "
                "without looking at this mixes the two"
            ),
            "rule": (
                "slide when all three hold: it is the only usable node on its image, it holds "
                "every cell assigned on that image, and its annotation area is at least the "
                "total area of the Tissue outlines on that image. Otherwise node. Measured each "
                "build from the image's own geometry; keyed on no donor id"
            ),
            "carried": (
                "nodes.parquet columns node_scope / node_scope_reason / "
                "area_ratio_node_over_image_tissue / n_usable_nodes_on_image / "
                "holds_every_assigned_cell_on_image, AND per row in "
                f"donor_id=*/{PARTITION_FILE} as node_scope"
            ),
            "slide_scope_nodes": [
                {
                    "donor_id": r["donor_id"], "image": r["image"],
                    "node_name": r["node_name"], "n_cells": r["n_cells"],
                    "area_um2": r["area_um2"],
                    "image_tissue_area_um2": r["image_tissue_area_um2"],
                    "area_ratio_node_over_image_tissue":
                        r["area_ratio_node_over_image_tissue"],
                    "n_usable_nodes_on_image": r["n_usable_nodes_on_image"],
                    "reason": r["node_scope_reason"],
                }
                for r in slide_scope
            ],
        },
        "qupath_cross_check": {
            "source": str(audit_log),
            #: Counted here, never transcribed.  ANNOT records are the annotations; adding the
            #: one READ_MS line per image gives a line count, and reporting that sum as the
            #: number of annotations is the error this block replaces.
            "audit_log_measured": log_measured,
            "geojson_annotations_measured": geojson_measured,
            "audit_log_agrees_with_geojson": (
                bool(log_measured.get("annotation_records")
                     == geojson_measured["features_total"])
                if log_measured.get("state") == "measured" else None
            ),
            "what_it_is": (
                "nChildren of the Tissue annotation inside each node region, from QuPath's own "
                "hierarchy. Every pLN_* annotation reads nChildren=0 -- the cells are children "
                "of Tissue -- so the per-Tissue count is the independent count of the same cells"
            ),
            "status": "a cross-check, never the source of membership",
            "nodes_checked": len(checked),
            "nodes_reproducing_exactly": len(checked) - len(differs),
            "nodes_differing": [
                {
                    "donor_id": r["donor_id"], "node_name": r["node_name"],
                    "derived_n_cells": r["n_cells"],
                    "qupath_tissue_child_count": r["qupath_tissue_child_count"],
                    "delta": r["delta_vs_qupath"],
                }
                for r in differs
            ],
            "unassigned_cells": unassigned_detail,
            #: Derived from this build's own unassigned_detail.  It was a sentence with a
            #: distance baked into it, and a rebuild over a different cohort re-emitted it.
            "unassigned_cells_note": unassigned_note,
        },
        "nodes_unavailable": [
            {"donor_id": r["donor_id"], "node_name": r["node_name"],
             "node_index": r["node_index"], "annotation_geom": r["annotation_geom"],
             "state": r["state"], "reason": r["reason"]}
            for r in node_rows if r["state"] == "unavailable"
        ],
        "ingested_cells_coverage": coverage,
        "donors": donors_manifest,
        "donors_withheld": withheld,
        "authorities_disagree_about": {
            d: ("admitted by phenocycler.cohort, withheld by pln_node_annotations.csv"
                if d in by_cohort else
                "admitted by pln_node_annotations.csv, withheld by phenocycler.cohort")
            for d in disagreement
        },
        "maintainer_action_required": {
            "no_intensity_measurements_in_the_pln_export": {
                "what": (
                    "four images carry rows with all 918 export fields present and only the "
                    "seven context fields filled -- no area, no solidity, no intensity, in any "
                    "compartment, on any marker, on 100% of their rows. Node membership is "
                    "unaffected: it needs only the centroids, which are present. Anything that "
                    "needs an intensity on these donors is unavailable until QuPath re-runs the "
                    "measurements. Not worked around here, and no value is imputed"
                ),
                **no_measurement,
            },
        },
    }
    (out_root / "MANIFEST.json").write_text(json.dumps(manifest, indent=2) + "\n")
    print(f"\n[info] {cells_on_disk:,} cells measured off {len(partitions)} partition file(s), "
          f"{manifest['totals']['cells_assigned_to_a_node']:,} assigned, "
          f"{manifest['totals']['cells_unassigned']:,} unassigned; "
          f"reconciles_with_roster={manifest['totals']['reconciles_with_roster']}")
    recon = manifest["output_root_reconciliation"]
    if recon["donors_gained"] or recon["donors_lost"] or recon["cells_delta"]:
        print(f"[info] the root changed: gained {recon['donors_gained'] or 'no donor'}, lost "
              f"{recon['donors_lost'] or 'no donor'}, {recon['cells_delta']:+,} cells against "
              f"{recon['before_this_run']['cells']:,} before this run")
    if cells_on_disk != cell_rows_counted_in_memory:
        print(f"[warn] the files hold {cells_on_disk:,} rows, this run wrote "
              f"{cell_rows_counted_in_memory:,} -- the MANIFEST reports the files")
    print(f"[info] {len(checked) - len(differs)}/{len(checked)} nodes reproduce QuPath exactly")
    if slide_scope:
        print("[info] slide-scope node(s) -- a per-node number here is a per-slide number: "
              + ", ".join(f"{r['donor_id']}/{r['node_name']}" for r in slide_scope))
    if log_measured.get("state") == "measured":
        print(f"[info] provenance log: {log_measured['annotation_records']} annotation records "
              f"over {log_measured['entries']} entries, geojson carries "
              f"{geojson_measured['features_total']}")
    if coverage.get("agrees") is False:
        print(f"[warn] donor coverage differs from {coverage['root']}: this layer has "
              f"{coverage['in_this_layer_but_not_ingested'] or 'no extra donor'}, that root has "
              f"{coverage['ingested_but_not_in_this_layer'] or 'no extra donor'}")
    print(f"[info] MANIFEST -> {out_root / 'MANIFEST.json'}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
