"""Stage 00 source adapter for explicitly recorded same-slide pilot sections.

This is called by the one Stage 00 module, not a second workflow controller.
The accepted acquisition table can feed the existing Stage 04 Xenium reader.
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
_ACQUISITION_COLUMNS = ("acquisition_id", "platform", "tissue", "roi", "donor_id",
                        "physical_section", "relationship", "path", "output", "condition")
_RELATION_COLUMNS = ("donor_id", "tissue", "physical_section", "relationship",
                     "phenocycler_acquisition_id", "xenium_acquisition_id",
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
    acquisitions, relationships = [], []
    for pair in pairs:
        serial, region = parse_bundle_path(pair.xenium_output.name)
        tissue_profile = profiles.load_all()[("tissue", pair.tissue)]
        if not tissue_profile.roi_prefixes:
            raise PilotInventoryError(f"{pair.tissue}: no ROI prefix in its tissue profile")
        pc_id = f"phenocycler:{pair.protein_entry}"
        base = {"tissue": pair.tissue, "donor_id": pair.donor_id,
                "physical_section": pair.protein_entry,
                "relationship": pair.relationship, "condition": ""}
        acquisitions.append({**base, "acquisition_id": pc_id, "platform": "phenocycler",
                             "roi": tissue_profile.roi_prefixes[0],
                             "path": str(pair.protein_image), "output": ""})
        acquisitions.append({**base, "acquisition_id": pair.xenium_acquisition_id,
                             "platform": "xenium", "roi": region,
                             "path": str(pair.xenium_output), "output": pair.xenium_output.name})
        relationships.append({"donor_id": pair.donor_id, "tissue": pair.tissue,
                              "physical_section": pair.protein_entry,
                              "relationship": pair.relationship,
                              "phenocycler_acquisition_id": pc_id,
                              "xenium_acquisition_id": pair.xenium_acquisition_id,
                              "protein_frame": pair.protein_frame,
                              "xenium_frame": pair.xenium_frame,
                              "crop_origin_x_px": pair.crop_origin_px[0],
                              "crop_origin_y_px": pair.crop_origin_px[1]})
    _write(task.output("acquisitions.csv"), _ACQUISITION_COLUMNS, acquisitions)
    _write(task.output("relationships.csv"), _RELATION_COLUMNS, relationships)
    _figure(task, relationships)
    return {"pairs": len(pairs), "acquisitions": len(acquisitions),
            "regions": len(relationships), "tissues": sorted({pair.tissue for pair in pairs}),
            "verified_crop_origins": sum(pair.crop_origin_px is not None for pair in pairs),
            "instrument_records": len(pairs) * len(INSTRUMENT_FILES),
            "raw_image_digests": "not yet recorded; candidate pilot only"}


def _figure(task: Task, rows: list[dict]) -> None:
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
           xlabel="Source acquisitions verified per physical section (n)",
           title=f"Same-slide source coverage: {len(rows)} sections, {len(rows) * 2} acquisitions")
    ax.invert_yaxis()
    ax.legend(loc="lower right", fontsize=8)
    fig.savefig(gates.figure_path(task, "same_slide_sources.png"), dpi=150)
    plt.close(fig)


def gate(task: Task, summary: dict) -> gates.Gate:
    pairs = load_pairs(task.input("same_slide_source_manifest"))
    checks = [
        gates.Check("every recorded physical section has both sources", summary["pairs"],
                    str(len(pairs)), summary["pairs"] == len(pairs)),
        gates.Check("two distinct assay acquisitions per section", summary["acquisitions"],
                    str(2 * len(pairs)), summary["acquisitions"] == 2 * len(pairs)),
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
