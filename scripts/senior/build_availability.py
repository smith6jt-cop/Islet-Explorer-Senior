#!/usr/bin/env python3
"""Measure donor x capability availability without turning missing inputs into zero.

Counts have explicit units: eligible current typing cells for phenocycler_panc;
inventoried sections for phenocycler_pln; typed cells (lineage_native only) for
Xenium; paired IF cells whose recorded hormone release check passes for restain_if;
nonempty exported islet outlines; finite core-or-peri islet trajectory summaries;
and maintainer-recorded pLN nodes. pLN export cell counts are not used as a
nucleus-QC denominator. Availability does not certify cross-donor comparability.

Source files are read without hive inference. Only typing/cells gets donor_id
from its explicit donor_id= partition, because its files lack that column.
No cell tables are joined, and no identity is inferred from coordinates or names.
Recorded failures/withholding outrank the existence of files. All non-measured
states carry a reason and a null int64 count.
"""

from __future__ import annotations

import argparse
import json
import math
import sys
from collections import Counter, defaultdict
from pathlib import Path

import pyarrow as pa
import pyarrow.compute as pc
import pyarrow.parquet as pq

REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO))

from scripts.senior.build_sections_layer import (
    add_common_arguments,
    filter_eligible_donors,
    norm_donor,
    read_inventory,
    section_inventories,
    write_layer,
)

STATES = ("measured", "blank", "unavailable", "withheld", "failed")
CAPABILITIES = (
    "phenocycler_panc", "phenocycler_pln", "xenium_panc", "xenium_pln",
    "restain_if", "islet_outline", "trajectory", "pln_nodes",
)
AVAILABILITY_SCHEMA = pa.schema([
    ("donor_id", pa.string()), ("capability", pa.string()), ("state", pa.string()),
    ("reason", pa.string()), ("n", pa.int64()), ("source", pa.string()),
])


def availability(donor: str, capability: str, state: str, source: Path,
                 *, reason: str = "", n: int | None = None) -> dict:
    if state not in STATES or capability not in CAPABILITIES:
        raise ValueError(f"Unknown capability/state: {capability}/{state}")
    if state == "measured":
        if not isinstance(n, int) or isinstance(n, bool) or n < 0:
            raise ValueError("A measured state requires a nonnegative integer count")
    elif not reason.strip() or n is not None:
        raise ValueError("Every non-measured state requires a reason and n=None")
    return {"donor_id": norm_donor(donor), "capability": capability, "state": state,
            "reason": reason, "n": n, "source": str(source)}


def parquet_files(root: Path) -> list[Path]:
    return [root] if root.is_file() else sorted(root.rglob("*.parquet"))


def batches(path: Path, columns: tuple[str, ...]):
    """ParquetFile reads physical columns with partitioning=None semantics."""
    parquet = pq.ParquetFile(path)
    missing = set(columns) - set(parquet.schema_arrow.names)
    if missing:
        raise ValueError(f"{path}: missing required columns {sorted(missing)}")
    return parquet.iter_batches(batch_size=131_072, columns=list(columns))


def donor_id(value) -> str:
    if value is None or not str(value).strip().isdigit():
        raise ValueError(f"Missing or invalid source donor_id: {value!r}")
    return norm_donor(value)


def typing_availability(root: Path) -> dict[str, dict]:
    counts = Counter()
    for path in sorted(root.glob("donor_id=*/*.parquet")):
        donor = donor_id(path.parent.name.removeprefix("donor_id="))
        if not filter_eligible_donors([donor]):
            continue
        counts[donor] += 0
        for batch in batches(path, ("qc_analysis_eligible",)):
            counts[donor] += pc.sum(batch["qc_analysis_eligible"]).as_py() or 0
    return {donor: availability(
        donor, "phenocycler_panc", "measured" if n else "unavailable", root,
        n=n if n else None,
        reason="qc_analysis_eligible_cells" if n else "no_qc_eligible_typing_cells",
    ) for donor, n in counts.items()}


def xenium_availability(root: Path) -> tuple[dict[tuple[str, str], dict], Counter]:
    counts = defaultdict(Counter)
    paired_if = Counter()
    for path in parquet_files(root):
        schema = pq.ParquetFile(path).schema_arrow
        columns = ("donor_id", "roi", "lineage_native")
        features = ("feat__if_INS", "feat__if_GCG")
        has_if = set(features).issubset(schema.names)
        for batch in batches(path, columns + (features if has_if else ())):
            table = pa.Table.from_batches([batch])
            # Count physical rows, including null native identities, rather than
            # count(lineage_native), which would discard untyped nulls.
            table = table.append_column("row_count", pa.repeat(pa.scalar(1, pa.int64()), len(table)))
            grouped = table.group_by(list(columns)).aggregate([("row_count", "sum")])
            for row in grouped.to_pylist():
                donor = donor_id(row["donor_id"])
                if not filter_eligible_donors([donor]):
                    continue
                roi = row["roi"]
                if roi != "panc" and not (isinstance(roi, str) and roi.startswith("pln_")):
                    raise ValueError(f"{path}: unrecognised recorded ROI {roi!r}")
                capability = "xenium_panc" if roi == "panc" else "xenium_pln"
                native = row["lineage_native"]
                state = "blank" if native is None or native == "" else (
                    "unavailable" if native == "Unavailable" else "measured")
                counts[donor, capability][state] += row["row_count_sum"]
            if has_if:
                valid = pc.fill_null(pc.and_(pc.is_finite(table[features[0]]),
                                            pc.is_finite(table[features[1]])), False)
                paired = table.filter(valid).group_by(["donor_id", "roi"]).aggregate([
                    ("row_count", "sum")])
                for row in paired.to_pylist():
                    donor = donor_id(row["donor_id"])
                    if row["roi"] == "panc" and filter_eligible_donors([donor]):
                        paired_if[donor] += row["row_count_sum"]
    result = {}
    for (donor, capability), count in counts.items():
        reason = (f"lineage_native: measured={count['measured']}; "
                  f"blank={count['blank']}; qc_unavailable={count['unavailable']}")
        state = "measured" if count["measured"] else (
            "unavailable" if count["unavailable"] else "blank")
        result[donor, capability] = availability(
            donor, capability, state, root, reason=reason,
            n=count["measured"] if state == "measured" else None,
        )
    return result, paired_if


def outline_availability(root: Path) -> dict[str, dict]:
    outlines = defaultdict(set)
    for path in parquet_files(root):
        for batch in batches(path, ("donor_id", "islet_key", "wkb")):
            for row in batch.to_pylist():
                donor = donor_id(row["donor_id"])
                if filter_eligible_donors([donor]):
                    outlines[donor]  # Keep donors whose exported outlines are all missing.
                    if row["wkb"] and row["islet_key"]:
                        outlines[donor].add(row["islet_key"])
    return {donor: availability(
        donor, "islet_outline", "measured" if keys else "unavailable", root,
        reason="" if keys else "no_islet_outline_exported", n=len(keys) if keys else None,
    ) for donor, keys in outlines.items()}


def trajectory_availability(path: Path) -> dict[str, dict]:
    islets = defaultdict(set)
    if path.is_file():
        for batch in batches(path, (
            "case_id", "islet_key", "sc_pseudotime_core", "sc_pseudotime_peri",
        )):
            for row in batch.to_pylist():
                donor = donor_id(row["case_id"])
                if not filter_eligible_donors([donor]):
                    continue
                islets[donor]
                finite = any(row[column] is not None and math.isfinite(row[column])
                             for column in ("sc_pseudotime_core", "sc_pseudotime_peri"))
                if finite and row["islet_key"]:
                    islets[donor].add(row["islet_key"])
    return {donor: availability(
        donor, "trajectory", "measured" if keys else "unavailable", path,
        reason="finite_core_or_peri_islet_summaries" if keys else "no_finite_islet_pseudotime",
        n=len(keys) if keys else None,
    ) for donor, keys in islets.items()}


def pln_availability(path: Path) -> dict[tuple[str, str], dict]:
    if not path.is_file():
        return {}
    records = read_inventory(path, ("donor_id", "roi", "image_file", "n_nodes"))
    nodes, sections = Counter(), defaultdict(set)
    seen = set()
    for row in records:
        key = row["donor_id"], row["roi"], row["image_file"]
        if key in seen:
            raise ValueError(f"{path}: repeated pLN section {key}")
        seen.add(key)
        n = int(row["n_nodes"])
        if n < 0:
            raise ValueError(f"{path}: negative node count")
        nodes[row["donor_id"]] += n
        sections[row["donor_id"]].add((row["roi"], row["image_file"]))
    result = {}
    for donor, n in nodes.items():
        result[donor, "pln_nodes"] = availability(
            donor, "pln_nodes", "measured", path, n=n, reason="recorded_node_annotations")
        result[donor, "phenocycler_pln"] = availability(
            donor, "phenocycler_pln", "measured", path, n=len(sections[donor]),
            reason="recorded_pheno_pln_sections")
    return result


def restain_availability(root: Path, sections: list[dict], paired_if: Counter,
                         cells_root: Path) -> dict[str, dict]:
    hormone_file = root / "hormone_panc_refined.json"
    hormones = {}
    if hormone_file.is_file():
        for record in json.loads(hormone_file.read_text()):
            donor = donor_id(record["donor_id"])
            if record["roi"] == "panc":
                if donor in hormones:
                    raise ValueError(f"{hormone_file}: duplicate donor {donor}")
                hormones[donor] = record["hormone"]
    gates = {}
    for path in sorted((root / "valis").glob("donor_id=*/alignment_refined.json")):
        record = json.loads(path.read_text())
        if record["roi"] == "panc":
            donor = donor_id(record["donor_id"])
            if donor in gates:
                raise ValueError(f"{root}: duplicate refined alignment for {donor}")
            gates[donor] = record["alignment"], path
    result = {}
    for donor in {row["donor_id"] for row in sections if row["roi"] == "panc"}:
        if donor in gates and gates[donor][0]["verdict"] == "FAIL":
            gate, path = gates[donor]
            result[donor] = availability(
                donor, "restain_if", "failed", path,
                reason=f"restain_registration_failed: median_offset_um={gate['median_offset_um']}; "
                       f"{gate['verdict_why']}",
            )
        elif donor in gates and gates[donor][0]["verdict"] != "PASS":
            gate, path = gates[donor]
            result[donor] = availability(
                donor, "restain_if", "unavailable", path,
                reason=f"restain_registration_{gate['verdict'].lower()}: {gate['verdict_why']}")
        elif paired_if[donor]:
            hormone = hormones.get(donor)
            if hormone is None:
                result[donor] = availability(
                    donor, "restain_if", "unavailable", hormone_file,
                    reason="restain_hormone_release_not_recorded")
                continue
            withheld = [marker for marker in ("INS", "GCG")
                        if hormone[marker]["clears_preregistered_bar"] is not True]
            if withheld:
                values = "; ".join(f"{marker}_auroc={hormone[marker]['auroc']}" for marker in withheld)
                result[donor] = availability(
                    donor, "restain_if", "withheld", hormone_file,
                    reason=f"hormone_calls_withheld: {values}; "
                           f"recorded_auroc_bar={hormone['preregistered_auroc_bar']}; "
                           f"paired_intensity_cells={paired_if[donor]}")
            else:
                result[donor] = availability(
                    donor, "restain_if", "measured", cells_root, n=paired_if[donor],
                    reason=f"paired_if_cells; hormone_release_source={hormone_file}")
    return result


def default_phenocycler_data(repo: Path) -> Path:
    root = repo / "core/phenocycler/data"
    # Some worktrees retain tracked CSVs here and mount the full tree at data/data.
    return root / "data" if (root / "data").is_symlink() else root


def build_availability(repo_root: Path = REPO, *, inventory_dir: Path | None = None,
                       phenocycler_data: Path | None = None) -> pa.Table:
    repo_root = Path(repo_root)
    inventory = inventory_dir or repo_root / "core/phenocycler/data/integration"
    core_data = phenocycler_data or default_phenocycler_data(repo_root)
    restain, xenium = section_inventories(inventory)
    typing = repo_root / "data/development/typing"
    cells = core_data / "integration/cells_xen"
    nodes = inventory / "pln_node_annotations.csv"
    summaries = typing / "trajectory/singlecell_islet_summary.parquet"
    observed, paired_if = xenium_availability(cells)
    observed.update(pln_availability(nodes))
    for capability, records in (
        ("phenocycler_panc", typing_availability(typing / "cells")),
        ("islet_outline", outline_availability(typing / "geometry")),
        ("trajectory", trajectory_availability(summaries)),
        ("restain_if", restain_availability(repo_root / "data/development/restain",
                                          restain, paired_if, cells)),
    ):
        observed.update({(donor, capability): row for donor, row in records.items()})
    donors = set(filter_eligible_donors(
        {r["donor_id"] for r in restain + xenium} | {key[0] for key in observed}))
    if not donors:
        raise ValueError("No eligible donors discovered")
    missing = {
        "phenocycler_panc": (typing / "cells", "no_typing_cells_exported"),
        "phenocycler_pln": (nodes, "no_pheno_pln_section_recorded"),
        "xenium_panc": (cells, "no_xenium_cells_exported_for_roi"),
        "xenium_pln": (cells, "no_xenium_cells_exported_for_roi"),
        "restain_if": (cells, "no_restain_pairing"),
        "islet_outline": (typing / "geometry", "no_islet_outline_exported"),
        "trajectory": (summaries, "no_islet_trajectory_exported"),
        "pln_nodes": (nodes, "no_pln_node_annotations_recorded"),
    }
    rows = []
    for donor in sorted(donors):
        for capability in CAPABILITIES:
            row = observed.get((donor, capability))
            if row is None:
                source, reason = missing[capability]
                if capability.startswith("xenium_"):
                    records = [r for r in xenium if r["donor_id"] == donor and
                               (r["roi"] == "panc") == (capability == "xenium_panc")]
                    if records and all(r["source_path"].upper() == "NA" for r in records):
                        source = inventory / "xenium_paths.csv"
                        reason = "xenium_source_path_recorded_unavailable"
                row = availability(donor, capability, "unavailable", source, reason=reason)
            rows.append(row)
    return pa.Table.from_pylist(rows, schema=AVAILABILITY_SCHEMA)


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    add_common_arguments(parser)
    parser.add_argument("--phenocycler-data", type=Path,
                        help="Source data tree, separate from the tracked section inventories")
    args = parser.parse_args(argv)
    table = build_availability(args.repo_root, inventory_dir=args.inventory_dir,
                               phenocycler_data=args.phenocycler_data)
    path = write_layer(table, args.layers_dir, "availability")
    print(json.dumps({"artifact": str(path), "rows": table.num_rows,
                      "donors": len(set(table["donor_id"].to_pylist())),
                      "states": dict(Counter(table["state"].to_pylist()))}))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
