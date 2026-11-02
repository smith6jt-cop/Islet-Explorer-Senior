#!/usr/bin/env python3
"""Read MALDI source records and recorded Xenium blocks; never infer section links."""
import csv
import hashlib
import json
import xml.etree.ElementTree as ET
import zipfile
from pathlib import Path

BASE = Path(__file__).resolve().parents[2]
OUT = BASE / "outputs/maldi_completion_20261102"
SOURCE = Path("/data/io60-current/sources/maldi/wilkes_2026-10")
INVENTORY = Path("/data/io60-current/releases/io60-maldi-candidate-1/00_inventory")
NS = {"s": "http://schemas.openxmlformats.org/spreadsheetml/2006/main"}


def sha(path):
    with path.open("rb") as handle:
        return hashlib.file_digest(handle, "sha256").hexdigest()


def ref(path):
    return {"path": str(path), "sha256": sha(path), "bytes": path.stat().st_size}


def rows(path):
    with path.open(newline="") as handle:
        return list(csv.DictReader(line for line in handle if not line.startswith("#")))


def workbook_rows(path):
    """Extract only section metadata sheets, retaining exact source cell addresses."""
    with zipfile.ZipFile(path) as archive:
        strings = ["".join(node.itertext()) for node in ET.fromstring(
            archive.read("xl/sharedStrings.xml")).findall("s:si", NS)]
        result = []
        for sheet, filename, keep in (
                ("Xenium 5k Donor Cohort", "xl/worksheets/sheet2.xml", {"A", "B", "M", "N", "O", "P", "Q"}),
                ("PanTissue5k", "xl/worksheets/sheet3.xml", {"A", "B", "C", "D", "E"})):
            for row in ET.fromstring(archive.read(filename)).findall("s:sheetData/s:row", NS):
                values = {}
                for cell in row.findall("s:c", NS):
                    address = cell.get("r")
                    if "".join(c for c in address if c.isalpha()) not in keep:
                        continue
                    value = cell.find("s:v", NS)
                    text = value.text if value is not None else ""
                    if cell.get("t") == "s":
                        text = strings[int(text)]
                    if text:
                        values[address] = text
                if values:
                    result.append({"sheet": sheet, "row": int(row.get("r")), "cells": values})
        return result


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    acquisitions = rows(INVENTORY / "acquisitions.csv")
    maldi = rows(INVENTORY / "maldi_acquisitions.csv")
    donors = {r["donor_id"] for r in maldi}
    source_zip = Path.home() / "Downloads/wilkes_pancreas_glyc.zip"
    with zipfile.ZipFile(source_zip) as archive:
        members = [{"name": r.filename, "bytes": r.file_size} for r in archive.infolist()]
    spreadsheet = BASE / "Xenium_Visium request.xlsx"
    sheet_rows = workbook_rows(spreadsheet)
    retained = []
    for row in sheet_rows:
        if row["row"] == 1:
            retained.append(row)
            continue
        label = row["cells"].get(f"A{row['row']}", "")
        # This spelling conversion only selects the donor's recorded row. It creates no
        # MALDI block, tissue, nearest section, or acquisition relationship.
        donor = label.removeprefix("HDL").removesuffix(".0").lstrip("0")
        if donor in donors:
            retained.append(row)
    candidates = []
    for section in maldi:
        candidates.append({"tissue_id": section["tissue_id"], "donor_id": section["donor_id"],
                           "neighbour_acquisition": None, "tissue": None,
                           "reference_candidates": [
                               {k: a.get(k) for k in ("acquisition_id", "platform", "tissue", "path", "qptiff", "roi", "recorded_slide")}
                               for a in acquisitions if a["donor_id"] == section["donor_id"]
                               and a["platform"] in {"phenocycler", "xenium"}],
                           "interpretation": "same-donor inventory candidates only; not established serial neighbours"})
    report = {
        "schema": "maldi-source-evidence/1", "scope": "read-only source audit before new visual inspection",
        "sources": [ref(p) for p in (SOURCE / "README.txt", SOURCE / "SHA256SUMS", source_zip,
                     spreadsheet, INVENTORY / "receipt.json", INVENTORY / "acquisitions.csv",
                     INVENTORY / "maldi_acquisitions.csv", BASE / "docs/PLAN_2026-10-31_MALDI_SEQUENTIAL.md")],
        "supplier_archive_members": members,
        "supplier_archive_has_optical_images_or_section_manifest": False,
        "raw_metadata_columns": ["spotId", "raster", "x", "y", "z", "Tissue_id", "Date"],
        "source_directory_members": sorted(p.name for p in SOURCE.iterdir()),
        "recorded_serial_scope": "Maintainer declaration in the dated plan: MALDI sections are sequential slides from the same series as PhenoCycler and Xenium; different cells.",
        "inspection_authority": "The dated plan specifies visual comparison of ion images with candidate reference slides, never Tissue_id parsing.",
        "unrecorded": ["MALDI physical slide identifiers", "exact cutting order or section-to-section distances",
                       "MALDI block-to-acquisition manifest", "independent MALDI optical scan"],
        "comparison_figure": ref(INVENTORY / "gate/maldi_sections.png"),
        "xenium_block_records": retained,
        "xenium_block_record_limit": "Run-sheet rows record Xenium blocks and serials, not MALDI block identity or physical adjacency. Some proposal/current run-sheet block values differ; preserve source sheet and cell.",
        "sections": candidates,
        "decision": "The source files support spot identity and same-series scope. Exact tissue/reference correspondence requires the specified image inspection; visual ambiguity cannot be resolved by donor or Tissue_id text."}
    (OUT / "source_evidence.json").write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps({"sections": len(candidates), "donors": len(donors), "archive_members": members,
                      "workbook_metadata_rows": len(retained), "report": str(OUT / "source_evidence.json")}, indent=2))


if __name__ == "__main__":
    main()
