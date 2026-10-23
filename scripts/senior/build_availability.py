#!/usr/bin/env python3
"""Measure donor x capability availability without turning missing inputs into zero.

Four capabilities, each with a card that reads it, and counts with explicit
units: ``qc_analysis_eligible`` cells of the current
``data/development/typing/cells`` build for phenocycler_panc and nonempty
exported islet outlines for islet_outline (both read by card_x1_gradient.R:284);
paired IF cells whose recorded hormone release check passes for restain_if
(card_x2_overlay.R:47; the untyped ``cells_xen`` ingest is read for that pairing
only, being the one tree that carries the IF columns); and finite core-or-peri
islet trajectory summaries for trajectory (card_s2_trajectory.R:56).
Availability does not certify cross-donor comparability.

Until 2026-09-23 the layer carried four more -- phenocycler_pln, xenium_panc,
xenium_pln and pln_nodes -- and no card read any of them. They were dropped,
with the two reads that fed only them: the Xenium releases (resolved through
``xenium_release``) and ``pln_node_annotations.csv``. The releases are still
read by distance_profiles, compartments and grid_composition. The donor rows
are discovered from the two section inventories and the four sources above.

Source files are read without hive inference. ``typing/cells`` takes donor_id
from its explicit ``donor_id=`` path segment, read as text, because its files
lack that column (inferred hive partitioning would type an all-digit donor as
int32 and stop matching the donor it names). No cell tables are joined, and no
identity is inferred from coordinates or names. Recorded failures/withholding
outrank the existence of files. All non-measured states carry a reason and a
null int64 count.

Exclusions are asked per capability, of ``cohort.excluded_on`` for the platforms
that capability reads (``CAPABILITY_PLATFORMS``). A donor whose exclusion covers
one of them -- 6457's covers its PhenoCycler data only -- keeps the capabilities
it is admitted to and is written ``withheld``, with the exclusion's reason, on
the rest.

The PhenoCycler pancreas count is whatever typing/cells holds when this runs --
since 2026-09-19, the islet-chain build of all 20 donors. This builder neither
types nor re-types anything; NO_VERIFIED_TYPING remains latched.
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

sys.path.insert(0, str(REPO / "core" / "multiplex5000"))

from phenocycler.cohort import DONOR_EXCLUSIONS, donor_of, excluded_on  # noqa: E402

from scripts.senior.availability_states import AVAILABILITY_STATES  # noqa: E402
from scripts.senior.build_sections_layer import (
    add_common_arguments,
    norm_donor,
    section_inventories,
    write_layer,
)

#: Derived, never restated. `scripts/senior/availability_states.py` is the one authority and
#: mirrors `IE_STATES` in `apps/R/ieCore/R/ie_availability.R`, order included. The five-element
#: copy that stood here could not carry `uncalled`, so a capability the producer had looked at
#: and resolved nothing for had to be reported as one of the other five.
STATES = AVAILABILITY_STATES
#: Each is read by a card (see the module docstring); a capability no card reads is not written.
CAPABILITIES = ("phenocycler_panc", "restain_if", "islet_outline", "trajectory")
#: The platforms each capability's data comes from, for ``cohort.excluded_on``. An exclusion can
#: cover one platform only -- 6457's covers its PhenoCycler data (maintainer, 2026-10-22) -- so a
#: donor is gated per capability: the IF pairing is Xenium cells carrying the post-Xenium IF
#: scanned on the same slide, and the other three are the PhenoCycler pancreas.
CAPABILITY_PLATFORMS = {
    "phenocycler_panc": ("phenocycler",),
    "restain_if": ("xenium", "if"),
    "islet_outline": ("phenocycler",),
    "trajectory": ("phenocycler",),
}
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


def admitted(donor: str, capability: str) -> bool:
    """Whether no central exclusion covers a platform this capability reads for the donor."""
    return not any(excluded_on(donor, p) for p in CAPABILITY_PLATFORMS[capability])


def withheld_by_cohort(donor: str, capability: str, source: Path) -> dict:
    """The row for a donor whose exclusion covers this capability: named, never a zero."""
    platforms = [p for p in CAPABILITY_PLATFORMS[capability] if excluded_on(donor, p)]
    return availability(
        donor, capability, "withheld", source,
        reason=f"centrally_excluded_on_{'_and_'.join(platforms)}: "
               f"{DONOR_EXCLUSIONS[donor_of(donor)]}")


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
    counts, n_rows, n_null = Counter(), Counter(), Counter()
    refused: dict[str, str] = {}
    for path in sorted(root.glob("donor_id=*/*.parquet")):
        donor = donor_id(path.parent.name.removeprefix("donor_id="))
        if not admitted(donor, "phenocycler_panc"):
            continue
        counts[donor] += 0
        if donor in refused:
            continue
        for batch in batches(path, ("qc_analysis_eligible",)):
            column = batch["qc_analysis_eligible"]
            # The boundary: only a boolean column (nulls allowed) can be summed
            # as a count of eligible cells.  pc.sum over int8 [1, -1] is 0 and
            # read as a measured zero; over a string column it raised out of the
            # build.  Neither is a measurement, so the donor is unavailable with
            # the dtype named, and nothing below is hardened against it.
            if not pa.types.is_boolean(column.type):
                refused[donor] = f"qc_analysis_eligible_not_boolean; dtype={column.type}"
                break
            n_rows[donor] += len(column)
            n_null[donor] += column.null_count
            # pc.sum skips nulls and returns null over a column that is all null;
            # ``or 0`` would read that as a measured zero.  The null tally above
            # is what keeps an unmeasured column apart from zero eligible cells.
            counts[donor] += pc.sum(column).as_py() or 0
    result = {}
    for donor, n in counts.items():
        if donor in refused:
            result[donor] = availability(donor, "phenocycler_panc", "unavailable", root,
                                         reason=refused[donor])
            continue
        # Rows with no measured value are a quantity of their own, carried beside
        # n on EVERY typing row and never folded into it: [True, None, None] is
        # one eligible cell and two unmeasured rows, [False, None, None, None] is
        # zero eligible among one measured and three unmeasured.  It travels in
        # the reason text with its number because availability() forbids n on a
        # non-measured row and the layer's six columns are pinned
        # (test_availability_layer.py:29); a reader splits on "; n_null=".
        null_note = f"n_null={n_null[donor]}"
        if n:
            result[donor] = availability(donor, "phenocycler_panc", "measured", root,
                                         n=n, reason=f"qc_analysis_eligible_cells; {null_note}")
        elif n_rows[donor] and n_null[donor] == n_rows[donor]:
            # Rows exist and none carries a value: the column was never measured
            # for this donor.  That is a missing measurement, not zero cells.
            result[donor] = availability(donor, "phenocycler_panc", "unavailable", root,
                                         reason=f"qc_analysis_eligible_unmeasured; {null_note}")
        else:
            result[donor] = availability(donor, "phenocycler_panc", "unavailable", root,
                                         reason=f"no_qc_eligible_typing_cells; {null_note}")
    return result


def paired_if_counts(root: Path) -> Counter:
    """Cells of a pancreas section carrying BOTH post-Xenium IF intensities.

    This stays on the integration ingest, which is the only tree that carries the
    ``feat__if_*`` columns at all: a Xenium typing release types transcripts and
    has no IF intensity in it.  A cell is paired only where both intensities are
    finite -- a missing intensity is unavailable, never zero and never negative.
    """
    paired_if = Counter()
    features = ("feat__if_INS", "feat__if_GCG")
    for path in parquet_files(root):
        schema = pq.ParquetFile(path).schema_arrow
        if not set(features).issubset(schema.names):
            continue
        for batch in batches(path, ("donor_id", "roi") + features):
            table = pa.Table.from_batches([batch])
            table = table.append_column("row_count", pa.repeat(pa.scalar(1, pa.int64()), len(table)))
            valid = pc.fill_null(pc.and_(pc.is_finite(table[features[0]]),
                                         pc.is_finite(table[features[1]])), False)
            paired = table.filter(valid).group_by(["donor_id", "roi"]).aggregate([
                ("row_count", "sum")])
            for row in paired.to_pylist():
                donor = donor_id(row["donor_id"])
                if row["roi"] == "panc" and admitted(donor, "restain_if"):
                    paired_if[donor] += row["row_count_sum"]
    return paired_if


def outline_availability(root: Path) -> dict[str, dict]:
    outlines = defaultdict(set)
    for path in parquet_files(root):
        for batch in batches(path, ("donor_id", "islet_key", "wkb")):
            for row in batch.to_pylist():
                donor = donor_id(row["donor_id"])
                if admitted(donor, "islet_outline"):
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
                if not admitted(donor, "trajectory"):
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
        if not admitted(donor, "restain_if"):
            continue
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
    summaries = typing / "trajectory/singlecell_islet_summary.parquet"
    # The ingest answers only the IF pairing, because it is the only tree that carries an IF
    # column; its lineage columns are empty since config.ini:153 was blanked (2026-09-23) and
    # are not read here.
    paired_if = paired_if_counts(cells)
    observed = {}
    for capability, records in (
        ("phenocycler_panc", typing_availability(typing / "cells")),
        ("islet_outline", outline_availability(typing / "geometry")),
        ("trajectory", trajectory_availability(summaries)),
        ("restain_if", restain_availability(repo_root / "data/development/restain",
                                          restain, paired_if, cells)),
    ):
        observed.update({(donor, capability): row for donor, row in records.items()})
    # A donor is in the layer when some capability admits it: the Xenium and IF inventories are
    # already gated on their own platform, so 6457 arrives through its Xenium sections, and its
    # PhenoCycler capabilities are then written withheld with the exclusion's reason below.
    donors = {d for d in {r["donor_id"] for r in restain + xenium} | {key[0] for key in observed}
              if any(admitted(d, capability) for capability in CAPABILITIES)}
    if not donors:
        raise ValueError("No eligible donors discovered")
    # A capability with no observed row still names WHICH source was asked and what it said.
    missing = {
        "phenocycler_panc": (typing / "cells", "no_typing_cells_exported"),
        "restain_if": (cells, "no_restain_pairing"),
        "islet_outline": (typing / "geometry", "no_islet_outline_exported"),
        "trajectory": (summaries, "no_islet_trajectory_exported"),
    }
    rows = []
    for donor in sorted(donors):
        for capability in CAPABILITIES:
            row = observed.get((donor, capability))
            if not admitted(donor, capability):
                row = withheld_by_cohort(donor, capability, missing[capability][0])
            elif row is None:
                source, reason = missing[capability]
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
