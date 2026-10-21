"""A small synthetic MultiPlex5000 project: every source stage 00 reads, for two donors.

``synthetic_project(**replace)`` writes the sources and a project file under ``tmp_path`` and
returns the project file's path. ``replace`` swaps one source's text (by its file name), so a test
can break exactly one thing. Donor 0112 is spelled ``HDL112`` and ``112`` in the sources, as in
the real ones, so the cohort profile's declared aliases are exercised.

What the instruments recorded is written the way they write it: each qptiff is a real (tiny) TIFF
whose first page carries the slide, experiment, dates and markers; each Xenium output carries
``experiment.xenium``, ``gene_panel.json`` and ``metrics_summary.csv``. The two pancreas images
were scanned under one experiment on one day (one run); the pLN slide was run alone. Both Xenium
slides share one run. The pLN output is held here only in part, so its records come from the
copy under ``xenium_records``, as on the workstation.
"""

from __future__ import annotations

import datetime
import html
import json
from pathlib import Path

import pytest
import yaml

NOTES = ("IO60panc2nd Notes\n"
         "6374: pancTail no pLN, sequential Xenium\n"
         "HDL112: pancTail and 1 pLN different-block, same-slide post-phenocycler Xenium\n")

SOURCES = {
    "xenium_paths.csv": "donor_id,roi,batch,source_path\n"
                        "6374,panc,B1,bundles/output-XETG00298__0041000__Pan__20240101__000000\n"
                        "112,panc,B1,bundles/output-XETG00298__0041001__Pancreas__20240101__000000\n"
                        "112,pln_1,B1,bundles/output-XETG00298__0041001__pLN__20240101__000000\n",
    "section_modes.csv": "donor_id,roi,xenium_serial,mode,evidence,note\n"
                         "112,panc,0041001,same_slide,same slide after PhenoCycler,\n",
    "restain_sections.csv": "image_file,donor_id,roi,xenium_serial\n"
                            "if_6374.qptiff,6374,panc,0041000\n",
    "restain_registration_inputs.csv": "image_file,donor_id,roi,xenium_serial,provenance_source\n"
                                       "if_6374.qptiff,6374,panc,0041000,restain_sections.csv\n",
    "pln_node_annotations.csv": "image_file,donor_id,roi,qupath_entry,node_classes,n_nodes,n_cells,"
                                "excluded,recorded_by,recorded_on,excluded_reason\n"
                                "6374pLN_Scan1.qptiff,6374,pln_1,1,pLN,1,100,False,t,2026-01-01,\n",
    "notes.txt": NOTES,
    "notes_copy.md": NOTES,
    "tissue_condition.csv": "# synthetic\n"
                            "acquisition_id,donor_id,tissue,condition,when,detail,recorded_by,"
                            "recorded_on\n"
                            "xenium:0041001:Pancreas,112,pancreas,damaged,before Xenium,x,t,"
                            "2026-01-01\n",
    "pln_blocks.csv": "image,entry_id,donor,pln_annotation,object_id,pieces,area_mm2,"
                      "overlap_with_tissue_um2,gap_to_tissue_um,measured_relation,"
                      "block_from_notes,notes\n"
                      "img112,1,HDL112,1,uuid-1,1,1.0,0.0,664.7,separate,different-block,x\n",
}

#: The project's name for each source file.
NAMES = {"xenium_paths.csv": "xenium_paths", "section_modes.csv": "section_modes",
         "restain_sections.csv": "restain_sections",
         "restain_registration_inputs.csv": "restain_registration_inputs",
         "pln_node_annotations.csv": "pln_node_annotations", "notes.txt": "donor_notes",
         "notes_copy.md": "donor_notes_copy", "pln_blocks.csv": "pln_blocks",
         "tissue_condition.csv": "tissue_condition",
         "cohort_spec.json": "cohort_spec", "qupath_manifest.json": "qupath_manifest"}

#: The synthetic panels: "2+ CD38" (CD3e and CD38, DAPI not counted) and the pLN slides'
#: "2 markers".
PANCREAS_MARKERS = ("DAPI", "CD3e", "CD38")
PLN_MARKERS = ("DAPI", "CD3e", "CD20")


def write_qptiff(path: Path, *, slide: str, experiment: str, position: str,
                 experiment_date: str, scanned: datetime.datetime, markers) -> Path:
    """A tiny TIFF whose first page says what a PhenoCycler qptiff's does."""
    import numpy as np
    import tifffile

    profile = {"experimentDescription": {"name": f"{experiment}_{position}", "wells": [
        {"wellName": "H1", "items": [{"markerName": m} for m in markers]}]}}
    description = ("<PerkinElmer-QPI-ImageDescription>"
                   f"<SlideID>{slide}</SlideID><Date>{experiment_date}</Date>"
                   "<AcquisitionSoftware>Fusion 2.4.1</AcquisitionSoftware>"
                   f"<ScanProfile>{html.escape(json.dumps(profile))}</ScanProfile>"
                   "</PerkinElmer-QPI-ImageDescription>")
    tifffile.imwrite(path, np.zeros((2, 2), np.uint8), description=description,
                     datetime=scanned, metadata=None)
    return path


def write_xenium_records(folder: Path, *, run: str, start: str, slide: str, region: str,
                         custom: int = 100, design: str = "6KTRAF") -> None:
    folder.mkdir(parents=True, exist_ok=True)
    (folder / "experiment.xenium").write_text(json.dumps({
        "run_name": run, "run_start_time": start, "slide_id": slide, "region_name": region,
        "panel_type": "add_on" if custom else "predesigned", "panel_design_id": design,
        "panel_num_targets_predesigned": 5001, "panel_num_targets_custom": custom}))
    (folder / "gene_panel.json").write_text(json.dumps(
        {"payload": {"panel": {"identity": {"design_id": design}, "num_gene_targets": custom}}}))
    (folder / "metrics_summary.csv").write_text(
        "run_name,region_name,panel_design_id\n" f"{run},{region},{design}\n")


@pytest.fixture
def synthetic_project(tmp_path):
    def build(*, release="r1", automatic=(), stages=("00_inventory",), **replace) -> Path:
        source = tmp_path / "sources"
        raw = tmp_path / "raw"
        pln = tmp_path / "pln_images"
        bundles = tmp_path / "xenium" / "bundles"
        copies = tmp_path / "xenium_records"
        for folder in (source, raw, pln, bundles, copies):
            folder.mkdir(parents=True, exist_ok=True)
        run, start = "R1_synthetic", "2024-01-01T10:00:00Z"
        write_xenium_records(bundles / "output-XETG00298__0041000__Pan__20240101__000000",
                             run=run, start=start, slide="0041000", region="Pan")
        write_xenium_records(bundles / "output-XETG00298__0041001__Pancreas__20240101__000000",
                             run=run, start=start, slide="0041001", region="Pancreas")
        # held here only in part: its records come from the copy of the complete output
        partial = bundles / "output-XETG00298__0041001__pLN__20240101__000000"
        partial.mkdir(exist_ok=True)
        (partial / "cells.parquet").write_bytes(b"")
        write_xenium_records(copies / partial.name, run=run, start=start, slide="0041001",
                             region="pLN")
        (copies / "MANIFEST.tsv").write_text(
            "# synthetic\noutput\tfile\tsource\tsha256\tbytes\tcopied_utc\n" + "".join(
                f"{partial.name}\t{f}\thpg:/complete/{partial.name}/{f}\t0\t0\t2026-01-01\n"
                for f in ("experiment.xenium", "gene_panel.json", "metrics_summary.csv")))
        images = []
        for donor, image, position, minute in (("6374", "6374_Scan1 - resolution #1", "Slide1", 0),
                                               ("0112", "HDL112_Scan3 - resolution #1", "Slide2",
                                                20)):
            path = write_qptiff(raw / f"{donor}.qptiff", slide=donor,
                                experiment="Synthetic_Pancreas", position=position,
                                experiment_date="2025-12-01T09:00:00Z",
                                scanned=datetime.datetime(2026, 1, 2, 14, minute),
                                markers=PANCREAS_MARKERS)
            stat = path.stat()
            images.append({"donor_id": donor, "image_id": image, "panel_id": "panel",
                           "segmentation_version": "seg",
                           "qptiff": {"path": str(path), "size_bytes": stat.st_size,
                                      "mtime_ns": stat.st_mtime_ns, "sha256": "0" * 64,
                                      "hash_mode": "sampled"}})
        write_qptiff(pln / "6374pLN_Scan1.qptiff", slide="6374pLN", experiment="Synthetic_pLN",
                     position="Slide1", experiment_date="2025-12-05T09:00:00Z",
                     scanned=datetime.datetime(2026, 1, 9, 12, 0), markers=PLN_MARKERS)
        files = dict(SOURCES)
        files["cohort_spec.json"] = json.dumps({
            "spec_version": 1, "cohort_name": "synthetic", "expected_donors": ["6374", "0112"],
            "defaults": {}, "images": [{"donor_id": i["donor_id"], "image_id": i["image_id"]}
                                       for i in images]})
        files["qupath_manifest.json"] = json.dumps({
            "images": images,
            "panels": [{"panel_id": "panel",
                        "channels": [{"marker": m, "channel_name": m} for m in PANCREAS_MARKERS]}]})
        files.update({k.replace("__", "."): v for k, v in replace.items()})
        for name, text in files.items():
            (source / name).write_text(text)
        project = {"schema": "multiplex5000-project/1", "id": "synthetic", "release": release,
                   "cohort": "io60_phenocycler", "stages": list(stages),
                   "automatic_gates": list(automatic),
                   "settings": {"xenium_root": str(tmp_path / "xenium")},
                   "sources": {**{NAMES[n]: str(source / n) for n in files},
                               "pancreas_qptiffs": str(raw), "pln_qptiffs": str(pln),
                               "xenium_records": str(copies)},
                   "releases_root": str(tmp_path / "releases")}
        path = tmp_path / "synthetic.yaml"
        path.write_text(yaml.safe_dump(project, sort_keys=False))
        return path
    return build
