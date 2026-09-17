#!/usr/bin/env python3
"""Aggregate current typing against UUID-keyed islet membership, never QuPath Parent.

The output is ``$ISLET_LAYERS_DIR/distance_profiles/data.parquet``. In addition
to the columns in ISLET_LAYERS_SPEC 1.3, ``state``/``reason`` distinguish absent
measurements and ``compartment`` preserves outline membership for the three
compartment card. Core membership follows ``in_islet``; the peri-islet expansion
is outside the core with 0 < signed distance <= 20 um. Other pancreas cells
remain in the gradient but do not enter either compartment.

``n_cells`` is the eligible, non-excluded denominator, repeated for each immune
subtype/provenance pair. Do not sum denominators across those pairs. ``Immune``
with ``resolved_by=all`` is the total; other rows preserve the existing calls,
including supporting-tier cells whose subtype was not resolved. No typing or
marker threshold is fitted here. Missing UUID coverage withholds the whole
donor's profile instead of silently changing its denominator.

Cell files have no donor column: only their explicit donor_id= partition is
read. Membership files are read with partitioning=None and joined on object_id
alone, including when their optional donor column is absent. Root-level audit
files such as membership_diff.parquet are not cell partitions.
"""

from __future__ import annotations

import argparse
import csv
import json
import math
import os
import sys
from pathlib import Path

import duckdb
import pyarrow as pa
import pyarrow.dataset as ds
import pyarrow.parquet as pq

REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO))
from scripts.senior.build_sections_layer import filter_eligible_donors, norm_donor

TYPING = REPO / "data" / "development" / "typing"
SUPERSEDED = (
    "Refusing superseded typing/subpanel/subpanel_calls.parquet: "
    "Beta 191,485 / Alpha 246,055 / Delta 62,797; "
    "current typing/cells has Beta 226,183 / Alpha 279,686 / Delta 90,658. "
    "Use typing/cells and typing/islet_membership."
)
SCHEMA = pa.schema([
    ("donor_id", pa.string()), ("donor_status", pa.string()),
    ("distance_bin_um", pa.float64()), ("n_cells", pa.int64()),
    ("n_immune", pa.int64()), ("immune_fraction", pa.float64()),
    ("immune_subtype", pa.string()), ("resolved_by", pa.string()),
    ("frame", pa.string()), ("compartment", pa.string()),
    ("state", pa.string()), ("reason", pa.string()),
])


def refuse_superseded(path: str | Path) -> None:
    """Reject the old input before any existence, schema, or output checks."""
    path = Path(path)
    if any(p.name == "subpanel_calls.parquet" for p in (path, path.resolve())):
        raise ValueError(SUPERSEDED)


def _partitions(root: Path) -> list[Path]:
    refuse_superseded(root)
    if root.is_file():
        return [root]
    return sorted(root.glob("donor_id=*/*.parquet"))


def _dataset(files: list[Path], required: set[str]) -> ds.Dataset:
    for path in files:
        refuse_superseded(path)
        missing = required - set(pq.ParquetFile(path).schema_arrow.names)
        if missing:
            raise ValueError(f"{path}: missing required columns {sorted(missing)}")
    # In particular, never infer a hive donor key from the membership tree.
    return ds.dataset([str(p) for p in files], format="parquet", partitioning=None)


def _status_by_donor(path: Path) -> dict[str, str]:
    if not path.is_file():
        return {}
    with path.open(newline="") as handle:
        reader = csv.DictReader(handle)
        fields = reader.fieldnames or []
        key = "donor_id" if "donor_id" in fields else "case_id"
        if key not in fields or "donor_status" not in fields:
            raise ValueError(f"{path}: requires case_id (or donor_id) and donor_status")
        statuses = {}
        for row in reader:
            if not row[key] or not row["donor_status"]:
                continue
            donor, status = norm_donor(row[key]), row["donor_status"].strip()
            if donor in statuses and statuses[donor] != status:
                raise ValueError(f"{path}: conflicting donor_status for {donor}")
            statuses[donor] = status
        return statuses


def _absence(donor: str, status: str | None, reason: str) -> dict:
    return {"donor_id": donor, "donor_status": status, "distance_bin_um": None,
            "n_cells": None, "n_immune": None, "immune_fraction": None,
            "immune_subtype": "Immune", "resolved_by": "all",
            "frame": "signed_distance_um", "compartment": None,
            "state": "unavailable", "reason": reason}


def build_distance_profiles(
    cells_root: str | Path = TYPING / "cells",
    membership_root: str | Path = TYPING / "islet_membership",
    donor_metadata: str | Path = TYPING / "app" / "donor_metadata.csv",
    *,
    bin_width_um: float = 10.0,
    donors: list[str] | None = None,
) -> pa.Table:
    """Build a small artifact in memory; no source or output file is modified."""
    for path in (cells_root, membership_root, donor_metadata):
        refuse_superseded(path)
    if not math.isfinite(bin_width_um) or bin_width_um <= 0:
        raise ValueError("bin_width_um must be finite and positive")
    cells_root, membership_root, donor_metadata = map(
        Path, (cells_root, membership_root, donor_metadata)
    )
    files_by_donor: dict[str, list[Path]] = {}
    for path in _partitions(cells_root):
        if not path.parent.name.startswith("donor_id="):
            raise ValueError(f"{path}: typing/cells requires an explicit donor_id= partition")
        donor = norm_donor(path.parent.name.removeprefix("donor_id="))
        files_by_donor.setdefault(donor, []).append(path)
    statuses = _status_by_donor(donor_metadata)
    requested = {norm_donor(d) for d in donors} if donors is not None else set(files_by_donor)
    if not requested:
        requested = set(statuses)
    eligible = sorted(filter_eligible_donors(requested))
    if not eligible:
        raise ValueError("unavailable(no_eligible_donors_in_cells_or_metadata)")

    member_files = _partitions(membership_root)
    rows, audit = [], {}
    with duckdb.connect(config={"threads": 2, "memory_limit": "2GB"}) as con:
        if member_files:
            con.register("members", _dataset(member_files, {
                "object_id", "in_islet", "signed_dist_islet_um",
            }))
            n, unique, missing = con.execute("""
                SELECT count(*), count(DISTINCT object_id),
                       count(*) FILTER (WHERE object_id IS NULL OR trim(object_id) = '')
                FROM members
            """).fetchone()
            if missing or n != unique:
                raise ValueError("islet_membership has missing or duplicate object_id values")

        for donor in eligible:
            status = statuses.get(donor)
            if not status:
                rows.append(_absence(donor, None, "donor_status_not_recorded"))
                continue
            files = files_by_donor.get(donor, [])
            if not files or not member_files:
                reason = "typing_cells_absent" if not files else "islet_membership_export_absent"
                rows.append(_absence(donor, status, reason))
                continue
            con.register("cells", _dataset(files, {
                "object_id", "cell_class", "excluded_reason", "immune_definition",
                "resolved_by", "qc_analysis_eligible",
            }))
            n, unique, missing = con.execute("""
                SELECT count(*), count(DISTINCT object_id),
                       count(*) FILTER (WHERE object_id IS NULL OR trim(object_id) = '')
                FROM cells
            """).fetchone()
            if missing or n != unique:
                raise ValueError(f"{donor}: typing/cells has missing or duplicate object_id values")
            con.execute("""
                CREATE OR REPLACE TEMP TABLE joined AS
                SELECT c.object_id, c.cell_class, c.immune_definition, c.resolved_by,
                       m.object_id AS member_id, m.in_islet, m.signed_dist_islet_um
                FROM cells c LEFT JOIN members m ON c.object_id = m.object_id
                WHERE c.qc_analysis_eligible IS TRUE
                  AND c.cell_class IS DISTINCT FROM 'Excluded'
                  AND coalesce(c.excluded_reason, '') = ''
            """)
            total, unmatched, invalid, untyped, disagree = con.execute("""
                SELECT count(*), count(*) FILTER (WHERE member_id IS NULL),
                       count(*) FILTER (WHERE member_id IS NOT NULL AND
                           (in_islet IS NULL OR signed_dist_islet_um IS NULL
                            OR NOT isfinite(signed_dist_islet_um))),
                       count(*) FILTER (WHERE cell_class IS NULL OR trim(cell_class) = ''
                           OR resolved_by IS NULL OR trim(resolved_by) = ''),
                       count(*) FILTER (WHERE in_islet != (signed_dist_islet_um <= 0))
                FROM joined
            """).fetchone()
            audit[donor] = {"source_cells": n, "eligible_cells": total,
                            "omitted_qc_or_excluded": n - total,
                            "missing_membership": unmatched, "invalid_membership": invalid,
                            "missing_typing": untyped, "membership_distance_disagreement": disagree}
            reason = None
            if not total:
                reason = "no_analysis_eligible_cells"
            elif unmatched:
                reason = f"islet_membership_missing_uuid:{unmatched}_of_{total}"
            elif invalid:
                reason = f"islet_membership_invalid_measurement:{invalid}_of_{total}"
            elif untyped:
                reason = f"typing_or_provenance_missing:{untyped}_of_{total}"
            if reason:
                rows.append(_absence(donor, status, reason))
                continue

            # Membership, rather than the sign of a possibly disagreeing distance,
            # decides the core. Exact zeros retain their own visible bin.
            con.execute("""
                CREATE OR REPLACE TEMP TABLE binned AS
                SELECT *, sign(signed_dist_islet_um) *
                          ceil(abs(signed_dist_islet_um) / ?) * ? AS distance_bin_um,
                       CASE WHEN in_islet THEN 'islet_core'
                            WHEN signed_dist_islet_um > 0 AND signed_dist_islet_um <= 20
                                THEN 'peri_islet'
                            ELSE 'other_pancreas' END AS compartment,
                       CASE WHEN coalesce(trim(immune_definition), '') = ''
                            THEN 'Immune (subtype unresolved)'
                            ELSE immune_definition END AS immune_subtype
                FROM joined
            """, [bin_width_um, bin_width_um])
            table = con.execute("""
                WITH denominators AS (
                    SELECT distance_bin_um, compartment, count(*) AS n_cells,
                           count(*) FILTER (WHERE cell_class = 'Immune') AS n_immune
                    FROM binned GROUP BY ALL
                ), pairs AS (
                    SELECT DISTINCT immune_subtype, resolved_by FROM binned
                    WHERE cell_class = 'Immune'
                ), numerators AS (
                    SELECT distance_bin_um, compartment, immune_subtype, resolved_by,
                           count(*) AS n_immune
                    FROM binned WHERE cell_class = 'Immune' GROUP BY ALL
                ), counts AS (
                    SELECT distance_bin_um, compartment, n_cells, n_immune,
                           'Immune' AS immune_subtype, 'all' AS resolved_by FROM denominators
                    UNION ALL
                    SELECT d.distance_bin_um, d.compartment, d.n_cells,
                           coalesce(n.n_immune, 0), p.immune_subtype, p.resolved_by
                    FROM denominators d CROSS JOIN pairs p LEFT JOIN numerators n
                      ON d.distance_bin_um = n.distance_bin_um AND d.compartment = n.compartment
                     AND p.immune_subtype = n.immune_subtype AND p.resolved_by = n.resolved_by
                )
                SELECT ? AS donor_id, ? AS donor_status, distance_bin_um, n_cells, n_immune,
                       n_immune::DOUBLE / n_cells AS immune_fraction, immune_subtype,
                       resolved_by, 'signed_distance_um' AS frame, compartment,
                       'measured' AS state, '' AS reason
                FROM counts ORDER BY immune_subtype, resolved_by, compartment, distance_bin_um
            """, [donor, status]).to_arrow_table()
            rows.extend(table.to_pylist())

    metadata = {
        "source_cells": str(cells_root.resolve()),
        "source_membership": str(membership_root.resolve()),
        "source_domain_root": str(donor_metadata.parent.resolve()),
        "bin_width_um": str(bin_width_um),
        "binning": "ceil(abs(distance)/width)*width with original sign; zero is its own bin",
        "denominator": "QC-eligible, non-excluded cells; repeated across subtype/provenance pairs",
        "audit": json.dumps(audit, sort_keys=True),
    }
    return pa.Table.from_pylist(rows, schema=SCHEMA).replace_schema_metadata(metadata)


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--cells-root", "--cells", type=Path, default=TYPING / "cells")
    parser.add_argument("--membership-root", type=Path, default=TYPING / "islet_membership")
    parser.add_argument("--donor-metadata", type=Path, default=TYPING / "app" / "donor_metadata.csv")
    parser.add_argument("--layers-dir", type=Path,
                        default=Path(os.environ.get("ISLET_LAYERS_DIR", "/data/islet_layers")))
    parser.add_argument("--bin-width-um", type=float, default=10.0)
    parser.add_argument("--donors", nargs="+")
    parser.add_argument("--dry-run", action="store_true", help="Measure without writing an artifact")
    args = parser.parse_args(argv)
    table = build_distance_profiles(args.cells_root, args.membership_root, args.donor_metadata,
                                    bin_width_um=args.bin_width_um, donors=args.donors)
    measured = table.filter(pa.compute.equal(table["state"], "measured"))
    totals = measured.filter(pa.compute.equal(measured["resolved_by"], "all"))
    print(f"{len(set(measured['donor_id'].to_pylist()))} measured donors; "
          f"{sum(totals['n_cells'].to_pylist()):,} cells; "
          f"{sum(totals['n_immune'].to_pylist()):,} immune; {table.num_rows:,} artifact rows")
    for row in table.filter(pa.compute.not_equal(table["state"], "measured")).to_pylist():
        print(f"{row['donor_id']}: {row['state']}({row['reason']})")
    if not args.dry_run:
        target = args.layers_dir / "distance_profiles" / "data.parquet"
        target.parent.mkdir(parents=True, exist_ok=True)
        temporary = target.with_suffix(".parquet.tmp")
        pq.write_table(table, temporary, compression="zstd")
        temporary.replace(target)
        print(target)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
