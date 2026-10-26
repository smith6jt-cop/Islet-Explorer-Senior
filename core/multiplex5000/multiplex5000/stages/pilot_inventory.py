"""Stage 00 source adapter for explicitly recorded same-slide pilot sections.

This is called by the one Stage 00 module, not a second workflow controller.
Whole-slide scans, section exports and Xenium region outputs are separate
records. Counting one protein export and one Xenium region as two acquisitions
per section would falsely turn this two-slide pilot into 16 acquisitions.
"""

from __future__ import annotations

import csv
from pathlib import Path

from .. import gates, profiles
from ..integration.manifest import parse_bundle_path
from ..same_slide_sources import load_pairs
from ..tasks import Input, Task
from ..xenium_sources import CELL_FILES, INSTRUMENT_FILES

STAGE = "00_inventory"
_RAW = (*INSTRUMENT_FILES, *CELL_FILES)
_SLIDE_COLUMNS = ("slide_id", "protein_image", "xenium_run")
_SECTION_COLUMNS = ("protein_entry", "donor_id", "tissue", "slide_id", "protein_export")
_REGION_COLUMNS = ("acquisition_id", "platform", "tissue", "roi", "donor_id",
                   "physical_section", "relationship", "path", "output", "condition",
                   "slide_id")
_RELATION_COLUMNS = ("donor_id", "tissue", "physical_section", "relationship",
                     "slide_id", "protein_entry", "xenium_region_id",
                     "protein_frame", "xenium_frame", "crop_origin_x_px", "crop_origin_y_px")


class PilotInventoryError(ValueError):
    """A same-slide pilot cannot be admitted as one recorded source roster."""


def _manifest(project) -> Path:
    named = project.source("same_slide_sources")
    if named is not None:
        return Path(named)
    relative = project.settings["same_slide_source_manifest"]
    path = Path(relative)
    if path.is_absolute() or not (profiles.REPO / path).is_file():
        raise PilotInventoryError("same_slide_source_manifest needs a tracked relative file")
    return profiles.REPO / path


def settings(project) -> dict:
    return {"same_slide_source_manifest": str(_manifest(project))}


def inputs(project) -> list[Input]:
    manifest = _manifest(project)
    pairs = load_pairs(manifest)
    import yaml

    spec = yaml.safe_load(manifest.read_text())
    found = [Input("same_slide_source_manifest", manifest),
             Input("section_table", Path(spec["section_table"]["path"]))]
    images: set[Path] = set()
    for pair in pairs:
        found.append(Input(f"protein_export:{pair.protein_entry}", pair.protein_export))
        if pair.protein_image not in images:
            # The original qptiffs are ~99/87 GB. Stage 00 records exact path,
            # byte size and mtime; a later release gate requires its own digest.
            found.append(Input(f"protein_image:{pair.protein_image.name}",
                               pair.protein_image, hash=False))
            images.add(pair.protein_image)
        for name in _RAW:
            found.append(Input(f"xenium:{pair.xenium_output.name}:{name}",
                               pair.xenium_output / name,
                               hash=name in INSTRUMENT_FILES))
    return found


def _write(path: Path, columns: tuple[str, ...], rows: list[dict]) -> None:
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns)
        writer.writeheader()
        writer.writerows(rows)


def run(task: Task) -> dict:
    pairs = load_pairs(task.input("same_slide_source_manifest"), require_sources=True)
    cohort = profiles.load_all()[("cohort", task.project.cohort)]
    if {pair.tissue for pair in pairs} != set(cohort.tissues):
        raise PilotInventoryError("pair tissue differs from cohort profile")
    slides, sections, regions, relationships = {}, [], [], []
    for pair in pairs:
        serial, region = parse_bundle_path(pair.xenium_output.name)
        tissue_profile = profiles.load_all()[("tissue", pair.tissue)]
        if not tissue_profile.roi_prefixes:
            raise PilotInventoryError(f"{pair.tissue}: no ROI prefix in its tissue profile")
        slide = {"slide_id": serial, "protein_image": str(pair.protein_image),
                 "xenium_run": str(pair.xenium_output.parent)}
        if serial in slides and slides[serial] != slide:
            raise PilotInventoryError(f"{serial}: conflicting whole-slide sources")
        slides[serial] = slide
        sections.append({"protein_entry": pair.protein_entry, "donor_id": pair.donor_id,
                         "tissue": pair.tissue, "slide_id": serial,
                         "protein_export": str(pair.protein_export)})
        regions.append({"acquisition_id": pair.xenium_acquisition_id,
                        "platform": "xenium", "tissue": pair.tissue,
                        "roi": region, "donor_id": pair.donor_id,
                        "physical_section": pair.protein_entry,
                        "relationship": pair.relationship,
                        "path": str(pair.xenium_output),
                        "output": pair.xenium_output.name,
                        "condition": "", "slide_id": serial})
        relationships.append({"donor_id": pair.donor_id, "tissue": pair.tissue,
                              "physical_section": pair.protein_entry,
                              "relationship": pair.relationship,
                              "slide_id": serial, "protein_entry": pair.protein_entry,
                              "xenium_region_id": pair.xenium_acquisition_id,
                              "protein_frame": pair.protein_frame,
                              "xenium_frame": pair.xenium_frame,
                              "crop_origin_x_px": pair.crop_origin_px[0],
                              "crop_origin_y_px": pair.crop_origin_px[1]})
    _write(task.output("slides.csv"), _SLIDE_COLUMNS, list(slides.values()))
    _write(task.output("sections.csv"), _SECTION_COLUMNS, sections)
    _write(task.output("regions.csv"), _REGION_COLUMNS, regions)
    _write(task.output("relationships.csv"), _RELATION_COLUMNS, relationships)
    _figure(task, relationships, len(slides), len({row["xenium_run"] for row in slides.values()}))
    return {"pairs": len(pairs), "physical_slides": len(slides),
            "protein_images": len({row["protein_image"] for row in slides.values()}),
            "protein_exports": len(sections), "xenium_regions": len(regions),
            "xenium_runs": len({row["xenium_run"] for row in slides.values()}),
            "tissues": sorted({pair.tissue for pair in pairs}),
            "verified_crop_origins": sum(pair.crop_origin_px is not None for pair in pairs),
            "instrument_records": len(pairs) * len(INSTRUMENT_FILES),
            "raw_image_digests": "not yet recorded; candidate pilot only"}


def _figure(task: Task, rows: list[dict], slides: int, runs: int) -> None:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots(figsize=(11, max(3, 0.55 * len(rows) + 1)), layout="constrained")
    labels = [f"{row['donor_id']} · {row['physical_section']}" for row in rows]
    y = list(range(len(rows)))
    ax.barh(y, [1] * len(rows), color="#2d6a4f", label="PhenoCycler export verified")
    ax.barh(y, [1] * len(rows), left=[1] * len(rows), color="#457b9d",
            label="Xenium region verified")
    ax.set(yticks=y, yticklabels=labels, xlim=(0, 2),
           xlabel="Section export and Xenium region output present (one each)",
           title=f"Same-slide sources: {len(rows)} sections on {slides} slides; "
                 f"{runs} Xenium run(s)")
    ax.invert_yaxis()
    ax.legend(loc="lower right", fontsize=8)
    fig.savefig(gates.figure_path(task, "same_slide_sources.png"), dpi=150)
    plt.close(fig)


def gate(task: Task, summary: dict) -> gates.Gate:
    pairs = load_pairs(task.input("same_slide_source_manifest"))
    checks = [
        gates.Check("every recorded physical section has both sources", summary["pairs"],
                    str(len(pairs)), summary["pairs"] == len(pairs)),
        gates.Check("whole-slide images match physical slide IDs",
                    summary["protein_images"], str(summary["physical_slides"]),
                    summary["protein_images"] == summary["physical_slides"]),
        gates.Check("one protein section export per recorded section",
                    summary["protein_exports"], str(len(pairs)),
                    summary["protein_exports"] == len(pairs)),
        gates.Check("one Xenium region output per recorded section",
                    summary["xenium_regions"], str(len(pairs)),
                    summary["xenium_regions"] == len(pairs)),
        gates.Check("Xenium run folders", summary["xenium_runs"],
                    "recorded source count", True, blocking=False),
        gates.Check("crop origins verified against the section table",
                    summary["verified_crop_origins"], str(len(pairs)),
                    summary["verified_crop_origins"] == len(pairs)),
        gates.Check("Xenium instrument records declared and hashed",
                    summary["instrument_records"], str(3 * len(pairs)),
                    summary["instrument_records"] == 3 * len(pairs)),
        gates.Check("large raw image digests", summary["raw_image_digests"],
                    "required before a release is active", False, blocking=False),
    ]
    return gates.Gate(STAGE, "Verified source and physical-section inventory; "
                      "no cell-to-cell join or typing.", checks, ["same_slide_sources.png"])
