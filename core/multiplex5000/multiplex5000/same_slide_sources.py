"""Read exact same-slide sources without inferring cell links from donor names.

The site manifest records a physical-section relationship. On the data host,
``require_sources=True`` also verifies the section table, export headers and raw
bundle presence. Registration and cell overlap remain separate Stage 07/08 work.
"""

from __future__ import annotations

import csv
import gzip
import hashlib
from dataclasses import dataclass
from pathlib import Path

import yaml

from .integration.manifest import parse_bundle_path
from .xenium_sources import CELL_FILES, INSTRUMENT_FILES


class SameSlideSourceError(ValueError):
    """A recorded same-slide pair is missing, ambiguous or inconsistent."""


@dataclass(frozen=True)
class SameSlidePair:
    tissue: str
    donor_id: str
    protein_entry: str
    protein_export: Path
    protein_image: Path
    xenium_acquisition_id: str
    xenium_output: Path
    relationship: str
    protein_frame: str
    xenium_frame: str
    crop_origin_px: tuple[int, int] | None = None


def load_pairs(manifest_path: Path, *, require_sources: bool = False) -> tuple[SameSlidePair, ...]:
    """Load explicit pairs; only verified sources receive a crop-frame offset."""
    manifest = yaml.safe_load(Path(manifest_path).read_text())
    if not isinstance(manifest, dict) or manifest.get("schema") != "multiplex5000-same-slide-sources/1":
        raise SameSlideSourceError("unknown same-slide source manifest schema")
    if manifest.get("relationship") != "same_section_after_phenocycler":
        raise SameSlideSourceError("pair relationship must be recorded as same-section after PC")
    tissue = manifest.get("tissue")
    if not isinstance(tissue, str) or not tissue.strip():
        raise SameSlideSourceError("pair tissue must be declared")
    if (manifest.get("protein_frame"), manifest.get("xenium_frame")) != (
            "crop_local_um", "region_local_um"):
        raise SameSlideSourceError("this manifest requires declared crop and region frames")
    rows = manifest.get("pairs")
    if not isinstance(rows, list) or not rows or len(rows) != manifest.get("expected_pairs"):
        raise SameSlideSourceError("same-slide source count differs from expected_pairs")
    section_spec = manifest.get("section_table")
    if not isinstance(section_spec, dict) or not isinstance(section_spec.get("sha256"), str):
        raise SameSlideSourceError("section table and its digest are required")
    section_path = Path(section_spec["path"])
    if not section_path.is_absolute():
        raise SameSlideSourceError("section table path must be absolute")
    sections = {}
    if require_sources:
        if not section_path.is_file():
            raise SameSlideSourceError(f"missing section table: {section_path}")
        if hashlib.sha256(section_path.read_bytes()).hexdigest() != section_spec["sha256"]:
            raise SameSlideSourceError("section table changed from its recorded digest")
        with section_path.open(newline="") as handle:
            for section in csv.DictReader(handle):
                entry = section["entry_name"]
                if entry in sections:
                    raise SameSlideSourceError(f"repeated protein section {entry}")
                sections[entry] = section
    pairs, entries, acquisitions, outputs = [], set(), set(), set()
    for row in rows:
        if not isinstance(row, dict) or set(row) != {
                "donor_id", "protein_entry", "protein_export", "protein_image",
                "xenium_acquisition_id", "xenium_output"}:
            raise SameSlideSourceError("pair must name its exact donor, entry and two sources")
        donor, entry, acquisition = (str(row[key]) for key in (
            "donor_id", "protein_entry", "xenium_acquisition_id"))
        export, image, output = (Path(row[key]) for key in (
            "protein_export", "protein_image", "xenium_output"))
        if not all(path.is_absolute() for path in (export, image, output)):
            raise SameSlideSourceError(f"{entry}: source paths must be absolute")
        serial, region = parse_bundle_path(output.name)
        if not serial or acquisition != f"xenium:{serial}:{region}":
            raise SameSlideSourceError(f"{entry}: Xenium acquisition differs from output identity")
        if export.name != f"{entry}.csv.gz":
            raise SameSlideSourceError(f"{entry}: protein export names another section")
        if entry in entries or acquisition in acquisitions or output in outputs:
            raise SameSlideSourceError(f"{entry}: repeated section or Xenium region")
        entries.add(entry)
        acquisitions.add(acquisition)
        outputs.add(output)
        origin = None
        if require_sources:
            section = sections.get(entry)
            if not section or section["donor_id"] != donor or section["slide"] != serial:
                raise SameSlideSourceError(f"{entry}: section table disagrees with the pair")
            if serial not in image.name:
                raise SameSlideSourceError(f"{entry}: protein image names another slide")
            missing = [str(path) for path in (export, image)
                       if not path.is_file() or path.stat().st_size == 0]
            missing += [str(output / name) for name in (*INSTRUMENT_FILES, *CELL_FILES)
                        if not (output / name).is_file() or (output / name).stat().st_size == 0]
            if missing:
                raise SameSlideSourceError(f"{entry}: missing sources: {missing}")
            with gzip.open(export, "rb") as handle:
                header = handle.readline()
            if hashlib.sha256(header).hexdigest() != manifest.get("protein_header_sha256"):
                raise SameSlideSourceError(f"{entry}: protein export header differs from PC67")
            origin = (int(section["bbox_x"]), int(section["bbox_y"]))
        pairs.append(SameSlidePair(tissue, donor, entry, export, image, acquisition, output,
                                   manifest["relationship"], manifest["protein_frame"],
                                   manifest["xenium_frame"], origin))
    return tuple(pairs)
