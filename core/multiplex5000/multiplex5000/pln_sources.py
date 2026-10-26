"""Bind dedicated pLN QuPath exports to recorded images and physical sections.

The pLN annotation roster supplies donor identity. A site index supplies one
exact Parquet path and segmentation ID for each image. No newest-file search,
filename-derived donor, or substitution of another export is allowed. The
roster's 2026-09-16 cell count is comparison evidence, not an expected count
for the later v3 segmentation.
"""

from __future__ import annotations

import csv
from dataclasses import dataclass
from pathlib import Path

import pyarrow.parquet as pq

from .cohort import excluded_on

INDEX_COLUMNS = ["image_file", "source_path", "segmentation_id"]
REQUIRED_COLUMNS = {"Image", "Object ID", "Centroid X µm", "Centroid Y µm",
                    "Nucleus: Area µm^2"}


class PLNSourceError(ValueError):
    """A dedicated pLN export is missing or its identity is ambiguous."""


@dataclass(frozen=True)
class PLNExport:
    image_file: str
    donor_id: str
    roi: str
    source_path: Path
    segmentation_id: str
    roster_cells: int
    export_cells: int | None


def _rows(path: Path) -> tuple[list[str], list[dict[str, str]]]:
    with Path(path).open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(line for line in handle if not line.startswith("#"))
        return list(reader.fieldnames or ()), list(reader)


def load_exports(roster_path: Path, index_path: Path, *, source_root: Path,
                 selected_images: set[str] | None = None,
                 require_sources: bool = False) -> tuple[PLNExport, ...]:
    """Return an exact, complete image roster; inspect selected files when asked.

    The index must name every roster image, even if a bounded pilot processes
    only ``selected_images``. Unselected files need not be present on this host.
    A selected image withheld by either recorded authority is refused. Paths
    in the tracked index are relative to the explicit site source root, so
    the same roster and index can be installed on another computer.
    """
    roster_fields, roster_rows = _rows(roster_path)
    if not {"image_file", "donor_id", "roi", "n_cells", "excluded",
            "excluded_reason"} <= set(roster_fields):
        raise PLNSourceError("pLN roster lacks image, donor, ROI or exclusion columns")
    index_fields, index_rows = _rows(index_path)
    if index_fields != INDEX_COLUMNS:
        raise PLNSourceError("pLN site index needs image_file,source_path,segmentation_id")
    source_root = Path(source_root)
    if not source_root.is_absolute():
        raise PLNSourceError("pLN export source root must be absolute")
    roster = {row["image_file"]: row for row in roster_rows}
    indexed = {row["image_file"]: row for row in index_rows}
    if (not roster or len(roster) != len(roster_rows) or len(indexed) != len(index_rows)
            or set(indexed) != set(roster)):
        raise PLNSourceError("pLN index and recorded image roster differ or repeat an image")
    indexed_paths = [row["source_path"] for row in index_rows]
    if len(set(indexed_paths)) != len(indexed_paths):
        raise PLNSourceError("pLN site index assigns one source to multiple images")
    chosen = set(roster) if selected_images is None else set(selected_images)
    if not chosen or not chosen <= set(roster):
        raise PLNSourceError(f"selected pLN images are outside the roster: {sorted(chosen-set(roster))}")
    paths: set[Path] = set()
    exports: list[PLNExport] = []
    for image_file in sorted(chosen):
        row, location = roster[image_file], indexed[image_file]
        donor = row["donor_id"].strip()
        if not donor or not row["roi"].strip():
            raise PLNSourceError(f"{image_file}: invalid donor or pLN ROI")
        if row["excluded"] not in {"True", "False"}:
            raise PLNSourceError(f"{image_file}: invalid recorded exclusion state")
        if row["excluded"] == "False" and row["excluded_reason"].strip():
            raise PLNSourceError(f"{image_file}: stale exclusion reason on admitted image")
        try:
            roster_cells = int(row["n_cells"])
        except ValueError as error:
            raise PLNSourceError(f"{image_file}: invalid recorded cell count") from error
        if roster_cells < 0:
            raise PLNSourceError(f"{image_file}: negative recorded cell count")
        relative = Path(location["source_path"])
        if (not relative.parts or relative.is_absolute()
                or ".." in relative.parts or relative == Path(".")):
            raise PLNSourceError(f"{image_file}: source must be relative to the site root")
        path = source_root / relative
        segmentation = location["segmentation_id"].strip()
        if not segmentation or path in paths:
            raise PLNSourceError(f"{image_file}: missing or repeated exact source")
        paths.add(path)
        reason = None
        if row["excluded"] == "True":
            reason = row["excluded_reason"].strip() or "recorded image exclusion"
        if excluded_on(donor, "phenocycler"):
            reason = "cohort.DONOR_EXCLUSIONS"
        if reason:
            raise PLNSourceError(f"{image_file}: selected image withheld by {reason}")
        export_cells = None
        if require_sources:
            if not path.is_file():
                raise PLNSourceError(f"{image_file}: exact source is missing: {path}")
            source = pq.ParquetFile(path)
            fields = set(source.schema_arrow.names)
            if not REQUIRED_COLUMNS <= fields or not any(
                    name.startswith("Cell: ") and name.endswith(": Mean") for name in fields):
                raise PLNSourceError(f"{image_file}: not a full-cell QuPath export")
            export_cells = source.metadata.num_rows
            if export_cells <= 0:
                raise PLNSourceError(f"{image_file}: empty QuPath export")
            first = next(source.iter_batches(batch_size=1, columns=["Image"]))
            declared = first.column(0)[0].as_py()
            if declared not in {image_file, f"{image_file} - resolution #1"}:
                raise PLNSourceError(f"{image_file}: source Image column names {declared!r}")
        exports.append(PLNExport(image_file, donor, row["roi"], path, segmentation,
                                 roster_cells, export_cells))
    return tuple(exports)
