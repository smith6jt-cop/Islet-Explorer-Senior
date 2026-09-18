#!/usr/bin/env python3
"""Build sections/data.parquet from the two maintainer-recorded inventories.

Each Xenium inventory row names a section, not a bundle-selection rule. Regions
sharing a serial retain different frames; different serials for one donor/ROI
retain different rows. ``image_file`` preserves the recorded Xenium source path
(a bundle directory), or the recorded IF image filename. The explicit NA path
remains a row with null image/serial metadata; availability records its reason.

Only restain_sections.csv supplies pixel_size_um. Xenium has no scale in these
inventories, so its scale stays null. No image is opened by either builder.

The helpers here are shared with build_availability.py so both builders use the
exact donor normalisation in docs/ISLET_LAYERS_SPEC.md, section 0.3.

``section_relationship`` -- A FILTERING AND SORTING ATTRIBUTE, AND NOTHING ELSE
------------------------------------------------------------------------------
Some donors carry pancreas and pancreatic lymph node in one image; most do not. That is
worth being able to filter and sort on, and it is derived from evidence rather than
declared:

``same_block``   the donor's PhenoCycler pancreas image carries a pLN annotation beside the
                 pancreas tissue annotation. Read from the per-image annotation census
                 ``data/development/typing/revised_cells/section_parents.json``, which
                 ``scripts/typing/rebuild_donor_from_revised_export.py`` writes from the
                 export's own ``Parent`` column. Measured on 6533: ``Annotation (pLN)``
                 34,841 cells beside ``Annotation (Tissue)`` 931,910 in
                 ``0080870_6533_pLNpanc_Scan3.er.qptiff``.
``same_slide``   this section is one physical section imaged twice: a ``restain_if`` row and
                 a ``xenium`` row sharing ``(donor_id, roi, xenium_serial)``. The post-Xenium
                 IF is a restain of the section Xenium ran on, and the recorded serial is
                 what says so.
``different_block``  reserved for positive evidence that a donor's two tissues came from
                 different blocks. NOTHING IN THIS CHECKOUT RECORDS BLOCK IDENTITY, so it is
                 declared and never assigned. It is deliberately not inferred from the
                 absence of ``same_block`` evidence: not measured is not measured.
``null``         no evidence either way; ``section_relationship_basis`` names the reason as
                 ``unavailable(...)``.

Precedence when a row has both, and why: ``same_block`` wins over ``same_slide``, because
``same_slide`` is already fully recoverable from this table (``modality`` plus a shared
``xenium_serial``) and ``same_block`` is recoverable from nothing else in it. Neither fact
is lost -- ``section_relationship_basis`` carries every relationship measured for the row.

THIS COLUMN LICENSES NO DISTANCE. Islet-to-lymphoid distance is not wanted, and two tissues
in one image share a pixel grid and no anatomy. ``refuse_cross_tissue_distance()`` in
``build_distance_profiles.py`` stands exactly as it is; nothing here weakens it.

Note for a reader who goes looking: ``docs/ISLET_LAYERS_SPEC.md`` section 1.3 and
``apps/R/ieCore/R/ie_frames.R`` both state that pancreas and pLN mounted on one slide are cut
from *different* blocks. That is a declaration with no cited measurement, and it is about
the Xenium slide layout rather than about a QuPath annotation census. This builder does not
restate it and does not contradict it: it records only what the export measured.
"""

from __future__ import annotations

import argparse
import csv
import json
import math
import os
import sys
from pathlib import Path

import pyarrow as pa
import pyarrow.parquet as pq

REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO / "core" / "phenocycler"))

from phenocycler.cohort import filter_eligible_donors
from phenocycler.integration.manifest import parse_bundle_path

SECTIONS_SCHEMA = pa.schema([
    ("donor_id", pa.string()),
    ("roi", pa.string()),
    ("modality", pa.string()),
    ("image_file", pa.string()),
    ("xenium_serial", pa.string()),
    ("pixel_size_um", pa.float64()),
    ("frame", pa.string()),
    ("section_relationship", pa.string()),
    ("section_relationship_basis", pa.string()),
])

#: Strongest measured co-location first. `different_block` is in the vocabulary and is never
#: assigned here -- see the module docstring.
RELATIONSHIP_PRECEDENCE = ("same_block", "same_slide", "different_block")

#: Where rebuild_donor_from_revised_export.py leaves the per-image annotation census.
PARENT_EVIDENCE = REPO / "data/development/typing/revised_cells/section_parents.json"


def norm_donor(value: str) -> str:
    """'0112' -> '112', using the shared layer contract without numeric coercion."""
    return str(value).strip().lstrip("0") or "0"


def read_inventory(path: Path, required: tuple[str, ...]) -> list[dict[str, str]]:
    """Read text identifiers, refusing blanks instead of inferring missing facts."""
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(line for line in handle if not line.lstrip().startswith("#"))
        if not set(required).issubset(reader.fieldnames or []):
            raise ValueError(f"{path}: required columns are {required}")
        rows = []
        for row in reader:
            if None in row:
                raise ValueError(f"{path}: too many CSV fields")
            row = {key: (value or "").strip() for key, value in row.items()}
            for key in required:
                if not row[key] or (row[key].upper() == "NA" and key != "source_path"):
                    raise ValueError(f"{path}: blank or NA required field {key}")
            if not row["donor_id"].isdigit():
                raise ValueError(f"{path}: donor_id must be a string of digits")
            row["donor_id"] = norm_donor(row["donor_id"])
            if filter_eligible_donors([row["donor_id"]]):
                rows.append(row)
    return rows


def section_inventories(inventory_dir: Path) -> tuple[list[dict], list[dict]]:
    restain = read_inventory(inventory_dir / "restain_sections.csv", (
        "donor_id", "roi", "image_file", "xenium_serial", "pixel_size_um",
    ))
    xenium = read_inventory(inventory_dir / "xenium_paths.csv", (
        "donor_id", "roi", "source_path",
    ))
    return restain, xenium


def read_parent_evidence(path: Path | None) -> dict[str, dict]:
    """donor_id -> that donor's per-image annotation census, keyed by the layer's donor spelling.

    Absent, unreadable or malformed is an empty mapping rather than an exception: a layer whose
    ``section_relationship`` is unknown still has to build, and every row then carries a named
    ``unavailable`` saying which evidence was missing.
    """
    if path is None or not Path(path).exists():
        return {}
    try:
        payload = json.loads(Path(path).read_text())
    except (json.JSONDecodeError, OSError):
        return {}
    donors = payload.get("donors")
    if not isinstance(donors, dict):
        return {}
    out = {}
    for donor, entry in donors.items():
        if isinstance(entry, dict):
            out[norm_donor(donor)] = entry
    return out


def same_block_evidence(entry: dict) -> str | None:
    """The basis string for a donor whose image carries a pLN annotation, else None."""
    classes = entry.get("pln_annotation_classes") or []
    counts = entry.get("annotation_parent_counts") or {}
    named = [name for name in classes if counts.get(name)]
    if not named:
        return None
    detail = ",".join(f"{name}:{counts[name]}" for name in sorted(named))
    return (f"same_block:pln_annotation_in_pancreas_image="
            f"{entry.get('image') or entry.get('source_tsv') or 'unknown'}({detail})")


def assign_section_relationship(rows: list[dict], evidence: dict[str, dict]) -> None:
    """Fill ``section_relationship`` / ``section_relationship_basis`` in place.

    Every relationship measured for a row goes into the basis; the enum takes the first of
    RELATIONSHIP_PRECEDENCE that was measured, so nothing is hidden by the choice.
    """
    same_slide_keys = {
        (row["donor_id"], row["roi"], row["xenium_serial"])
        for row in rows if row["modality"] == "restain_if" and row["xenium_serial"]
    }
    for row in rows:
        measured: dict[str, str] = {}
        key = (row["donor_id"], row["roi"], row["xenium_serial"])
        if row["xenium_serial"] and key in same_slide_keys:
            measured["same_slide"] = (
                f"same_slide:restain_if_and_xenium_share_serial={row['xenium_serial']}")
        entry = evidence.get(row["donor_id"])
        if entry is not None:
            basis = same_block_evidence(entry)
            if basis:
                measured["same_block"] = basis
        if measured:
            row["section_relationship"] = next(
                name for name in RELATIONSHIP_PRECEDENCE if name in measured)
            row["section_relationship_basis"] = "; ".join(
                measured[name] for name in RELATIONSHIP_PRECEDENCE if name in measured)
        else:
            row["section_relationship"] = None
            reason = ("annotation_census_not_built" if not evidence
                      else "donor_absent_from_annotation_census" if entry is None
                      else "no_pln_annotation_in_this_donors_image")
            row["section_relationship_basis"] = f"unavailable({reason})"


def build_sections(inventory_dir: Path, parent_evidence: Path | None = PARENT_EVIDENCE) -> pa.Table:
    restain, xenium = section_inventories(Path(inventory_dir))
    rows = []
    for row in xenium:
        source = row["source_path"]
        serial = None
        if source.upper() != "NA":
            serial, _ = parse_bundle_path(source)
            if not serial:
                raise ValueError(f"xenium_paths.csv: unrecognised bundle path {source!r}")
        rows.append({
            "donor_id": row["donor_id"], "roi": row["roi"], "modality": "xenium",
            "image_file": source if serial else None, "xenium_serial": serial,
            "pixel_size_um": None,
            "frame": f"xenium_um:{row['donor_id']}:{row['roi']}",
        })
    recorded_serials = {(r["donor_id"], r["roi"], r["xenium_serial"]) for r in rows}
    for row in restain:
        scale = float(row["pixel_size_um"])
        if not math.isfinite(scale) or scale <= 0:
            raise ValueError(f"restain_sections.csv: invalid pixel_size_um for {row['image_file']}")
        if (row["donor_id"], row["roi"], row["xenium_serial"]) not in recorded_serials:
            raise ValueError(f"restain_sections.csv: section absent from xenium_paths.csv: {row}")
        rows.append({
            "donor_id": row["donor_id"], "roi": row["roi"], "modality": "restain_if",
            "image_file": row["image_file"], "xenium_serial": row["xenium_serial"],
            "pixel_size_um": scale, "frame": f"if_um:{row['donor_id']}:{row['roi']}",
        })
    assign_section_relationship(rows, read_parent_evidence(parent_evidence))
    unique = {}
    for row in rows:
        key = tuple(row[k] for k in ("donor_id", "roi", "modality", "xenium_serial"))
        if key in unique and unique[key] != row:
            raise ValueError(f"Conflicting records for section {key}")
        unique[key] = row
    if not unique:
        raise ValueError("The section inventories contain no eligible sections")
    return pa.Table.from_pylist(
        sorted(unique.values(), key=lambda r: tuple(str(r[k] or "") for k in (
            "donor_id", "roi", "modality", "xenium_serial"))),
        schema=SECTIONS_SCHEMA,
    )


def write_layer(table: pa.Table, layers_dir: Path, name: str) -> Path:
    """Replace this builder's single artifact only after the parquet write succeeds."""
    target = Path(layers_dir) / name / "data.parquet"
    target.parent.mkdir(parents=True, exist_ok=True)
    temporary = target.with_suffix(".parquet.tmp")
    try:
        pq.write_table(table, temporary, compression="zstd")
        temporary.replace(target)
    finally:
        temporary.unlink(missing_ok=True)
    return target


def add_common_arguments(parser: argparse.ArgumentParser) -> None:
    parser.add_argument("--repo-root", type=Path, default=REPO)
    parser.add_argument("--inventory-dir", type=Path,
                        help="Default: <repo-root>/core/phenocycler/data/integration")
    parser.add_argument("--layers-dir", type=Path,
                        default=Path(os.environ.get("ISLET_LAYERS_DIR", "/data/islet_layers")))


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    add_common_arguments(parser)
    # Only this builder reads the annotation census, so the flag stays out of
    # add_common_arguments(), which build_availability.py shares.
    parser.add_argument("--parent-evidence", type=Path, default=PARENT_EVIDENCE,
                        help="per-image annotation census written by "
                             "scripts/typing/rebuild_donor_from_revised_export.py; absent means "
                             "section_relationship is unavailable, never guessed")
    args = parser.parse_args(argv)
    inventory = args.inventory_dir or args.repo_root / "core/phenocycler/data/integration"
    table = build_sections(inventory, args.parent_evidence)
    path = write_layer(table, args.layers_dir, "sections")
    relationship: dict[str, int] = {}
    for value in table["section_relationship"].to_pylist():
        key = value or "unavailable"
        relationship[key] = relationship.get(key, 0) + 1
    print(json.dumps({"artifact": str(path), "sections": table.num_rows,
                      "donors": len(set(table["donor_id"].to_pylist())),
                      "parent_evidence": str(args.parent_evidence),
                      "parent_evidence_present": Path(args.parent_evidence).exists(),
                      "section_relationship": dict(sorted(relationship.items()))}))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
