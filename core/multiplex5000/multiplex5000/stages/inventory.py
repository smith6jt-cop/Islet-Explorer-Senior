"""Stage 00 -- the source inventory: every acquisition, every recorded relationship, what is missing.

What the release is built from, before anything is computed from it:

* **acquisitions** -- the PhenoCycler pancreas images (the QuPath manifest, checked against the
  cohort spec), the dedicated pLN slides (``pln_node_annotations.csv``), the Xenium region
  outputs (``xenium_paths.csv``) and the post-Xenium IF images (the restain tables). Each gets an
  ``acquisition_id`` made of its platform and the name that platform recorded for it --
  ``phenocycler:<QuPath image or pLN image file>``, ``xenium:<slide>:<region as experiment.xenium
  records it>``, ``if:<image file>`` -- never a name worked out from a filename. Raw files are
  checked against their records by size and modification time; nothing large is read.
* **panels and runs** -- what each instrument recorded (``multiplex5000.acquisition_records``):
  every qptiff's first page in both PhenoCycler image folders, and each Xenium region output's
  ``experiment.xenium``, ``gene_panel.json`` and ``metrics_summary.csv``. Every acquisition carries
  its panel, named the maintainer's way (``multiplex5000.panel_names``), its run, the run's date
  and the slide or slides run with it; ``runs.csv`` and the ``runs.png`` figure list every run. A
  Xenium output this machine holds only in part is read from the copy of its records taken from
  the complete output (the site's ``xenium_records``, ``MANIFEST.tsv`` beside them).
* **disease status and batches** -- each donor's disease status (ND, Aab+, T1D), from the
  project's donor metadata, else nPOD's case export, checked against each other where both have
  the donor; and the batch number each Xenium run name carries. Donors with no PhenoCycler
  acquisition and no notes line are Xenium-only, their rows kept apart in the coverage.
* **known damage** -- the damage the maintainer recorded, by acquisition and tissue, in the
  tracked ``tissue_condition.csv``; marked on exactly the damaged sections in the coverage, in the
  acquisitions and runs, and listed on the gate.
* **relationships** -- what ``multiplex5000.relationships`` reads from the donor notes,
  ``pln_blocks.csv`` and ``section_modes.csv``, with every disagreement between them.
* **blocks** -- the block each section was cut from, as the maintainer stated it (2026-10-21):
  pancreas sections from the pancreas block; a dedicated PhenoCycler pLN slide from a pLN block;
  a Xenium pLN section from a different block, except on a sequential section whose donor-notes
  line gives same-block pLN, where it is that same-block pLN.
* **coverage** -- by physical slide, one row per donor: the PhenoCycler pancreas slide (its
  pancreas; the pLN/lymphoid the donor notes put on it, same-block or different-block, checked
  against ``pln_blocks.csv``; its PhenoCycler scan; and, where Xenium ran on that same slide after
  PhenoCycler, that Xenium), the Xenium slide (its regions, its run and the IF scanned after it)
  and the dedicated pLN slide. Each acquisition says its panel, its run's date and the donor of
  the slide run with it; each donor's datasets and slides are counted -- one acquisition is one
  dataset, so a Xenium slide counts once however many regions it has. The rows of the sequential
  and the same-slide donors are kept apart. Each slide group opens with the slide's number and a
  thumbnail of the slide: a qptiff's own thumbnail page, a Xenium slide's overview scan (copied
  beside its records where this machine lacks it). As a table and the gate figure; the notes' cautions are
  quoted word for word on the gate, and nothing in them is counted or designated.

The gate blocks on: a repeated acquisition id, a cohort donor the notes do not describe, a cohort
image the QuPath manifest does not have, a raw image whose size is not its record's, any
disagreement between the sources that is not one already known and explained, an instrument
record that is missing, unreadable or at odds with itself, a pancreas image whose markers are not
its manifest panel's, a panel the naming rule does not cover, a run of more than two slides, a
recorded tissue condition that names no acquisition or an unknown tissue, a donor with no
recorded disease status or two sources that disagree on it, a Xenium pLN section
whose block the notes leave unsettled, and a Xenium output folder that names a region its record
does not, until the reason is known.

MultiPlex5000 D3 (2026-10-20).
"""

from __future__ import annotations

import csv
import json
import re
from collections import Counter
from pathlib import Path

from .. import acquisition_records, gates, panel_names, profiles, relationships, roots
from ..cohort import DONOR_EXCLUSIONS, excluded_on
from ..integration import tissues
from ..integration.manifest import parse_bundle_path
from ..integration.tissues import LYMPH_NODE, PANCREAS
from ..tasks import Input, Task

STAGE = "00_inventory"

#: The coverage is laid out by physical slide (maintainer, 2026-10-21): which sections were on one
#: slide, the block each was cut from, every platform acquisition of that slide with its run, and
#: totals that count datasets and slides, never sections. Left to right: the PhenoCycler pancreas
#: slide -- its sections as the donor notes give them, its PhenoCycler scan and, where Xenium ran on
#: that same slide after PhenoCycler (section_modes.csv's same_slide rows), that Xenium; the Xenium
#: slide -- its regions as xenium_paths.csv lists them, its run, and the IF scanned after it; and
#: the dedicated PhenoCycler pLN slide. One acquisition is one dataset, however many sections or
#: regions it covers.
PC_SLIDE = "PhenoCycler pancreas slide"
XENIUM_SLIDE = "Xenium slide"
PLN_SLIDE = "PhenoCycler pLN slide"
PHENOCYCLER = "PhenoCycler"
XENIUM = "Xenium"
POST_XENIUM_IF = "post-Xenium IF"
SLIDE_NUMBER = "slide"
PANCREAS_SECTION = "pancreas"
PLN_SECTION = "pLN/lymphoid"         # a dedicated pLN slide's sections, counted
DATASETS = "datasets"
SLIDES = "slides"
ACQUISITION_COLUMNS = (PHENOCYCLER, XENIUM, POST_XENIUM_IF)
#: The platform of each acquisition column, as the cohort's exclusions name them, and as the
#: coverage writes them.
PLATFORM_OF = {PHENOCYCLER: "phenocycler", XENIUM: "xenium", POST_XENIUM_IF: "if"}
PLATFORM_NAMES = {"phenocycler": "PhenoCycler", "xenium": "Xenium", "if": "IF"}
#: Columns that are not a section of the slide, so recorded damage never falls on them.
NOT_SECTIONS = (*ACQUISITION_COLUMNS, SLIDE_NUMBER)
#: The block a section was cut from: the pancreas block (a pLN in it is same-block) or another --
#: the "different block" of a different-block pLN, and the pLN block of a dedicated pLN slide.
PANCREAS_BLOCK = "pancreas block"
OTHER_BLOCK = "different block"
#: A section's fill says its tissue -- the tissue profiles' PANCREAS or LYMPH_NODE -- and nothing
#: else (maintainer, 2026-10-23: "use only color for tissue"); its block is the column it is
#: counted in, same-block or different-block pLN.
#: What a cell says for an acquisition a table records with no output listed for it, and for a
#: section whose block the rule does not settle.
NOT_LISTED = "no output listed"
UNSETTLED = "block not settled"
#: The Xenium design each donor's notes line gives; the coverage separates its rows by it. A
#: donor with no PhenoCycler acquisition and no notes line has Xenium only (6512 and 6618, added
#: by the maintainer 2026-10-22).
DESIGNS = {relationships.SEQUENTIAL: "sequential",
           relationships.SAME_SECTION_AFTER_PHENOCYCLER: "same slide"}
XENIUM_ONLY = "Xenium only"
#: Where each donor's disease status is read, in order -- the project's donor metadata, then
#: nPOD's case export for a donor the metadata does not have -- and its id and status columns.
#: The spellings each uses are the cohort profile's disease_status.
STATUS_SOURCES = (("donor_metadata", "donor_id", "disease.status"),
                  ("npod_cases", "case_id", "donor_type"))

#: Where stage 00 reads what the instruments recorded: the qptiff folders of both PhenoCycler
#: projects (every slide there can be another's run partner), and the copies of the records of
#: Xenium outputs this machine holds only in part.
QPTIFF_FOLDERS = ("pancreas_qptiffs", "pln_qptiffs")
#: The platform profile that names each platform's panels.
PANEL_PROFILES = {"phenocycler": "phenocycler", "xenium": "xenium", "if": "post_xenium_if"}
XENIUM_RECORDS = ("experiment.xenium", "gene_panel.json", "metrics_summary.csv")
#: The Xenium instrument's scan of the whole slide, every section on it -- the coverage's thumbnail
#: of a Xenium slide; read from the output's aux_outputs, else from the copy beside its records.
OVERVIEW = "overview_scan.png"
#: How tall a thumbnail is kept, in pixels; the figure draws it under an inch.
THUMBNAIL_PIXELS = 240

#: Xenium output folders that name a region other than the one their experiment.xenium records,
#: with the reason, once someone has looked. The gate blocks on any other.
_PANT_IS_PLN = ("run MM044: the instrument named this region PanT; on HiPerGator the complete output "
                "keeps that name, and this workstation's copy was renamed pLN. Its "
                "analysis_summary.html shows lymph node, as the folder name says (maintainer, "
                "2026-10-21)")
_MM046 = ("run MM046: the instrument named the regions Region_1 and Region_2; the folders, on "
          "HiPerGator and on this workstation, were renamed by tissue after export. Region_1 is the "
          "pancreas and Region_2 the lymph node, as the folder names say (maintainer, 2026-10-21)")
KNOWN_RENAMED_FOLDERS: dict[str, str] = {
    "output-XETG00298__0059865__pLN__20250328__202338": _PANT_IS_PLN,      # 6539
    "output-XETG00298__0059907__pLN__20250328__202338": _PANT_IS_PLN,      # 6436
    "output-XETG00298__0080824__6450_Panc__20251030__175543": _MM046,     # 6450, Region_1
    "output-XETG00298__0080824__6450_pLN__20251030__175544": _MM046,      # 6450, Region_2
    "output-XETG00298__0080836__6505_Panc__20251030__175544": _MM046,     # 6505, Region_1
    "output-XETG00298__0080836__6505_pLN__20251030__175544": _MM046,      # 6505, Region_2
}

#: Labels typed at an instrument that name another donor, and why the listed donor stands: a
#: Xenium cassette name by its slide's serial, a PhenoCycler SlideID by its qptiff. A label that
#: names another donor blocks the gate until it is listed here.
KNOWN_LABEL_ERRORS: dict[str, str] = {
    "0059945": "typed 6436 at the instrument; the slide is 6593's, as listed (maintainer, "
               "2026-10-23)",
    "0041393": "typed nPOD6553 at the instrument; the Xenium run sheet the maintainer supplied "
               "(Xenium_Visium request.xlsx, PanTissue5k, 2026-10-23) lists slide 0041393 as "
               "6533PB03/PLN02, as listed",
    "6523panc_Scan1.er.qptiff": "6523 was entered as 6623 at acquisition; the merged qptiff and "
                                "the QuPath project's images were renamed 6523, and only the raw "
                                "and intermediate cycle files on the share keep 6623 "
                                "(maintainer, 2026-10-23)",
}

#: The pLN/lymphoid the donor notes put on each PhenoCycler pancreas slide, by block: (the
#: column on that slide, the block as relationships names it, the DonorNote field that holds the
#: notes' count).
SAME_BLOCK = "same-block pLN/lymphoid"
DIFFERENT_BLOCK = "different-block pLN/lymphoid"
DESIGNATIONS = ((SAME_BLOCK, relationships.SAME_BLOCK, "same_block_pln"),
                (DIFFERENT_BLOCK, relationships.DIFFERENT_BLOCK, "different_block_pln"))
ANNOTATED = " in pln_blocks.csv"     # suffix of the column that holds pln_blocks.csv's pieces

_TRACKED = {  # name -> file name under core/phenocycler/data/integration
    "xenium_paths": "xenium_paths.csv", "section_modes": "section_modes.csv",
    "restain_sections": "restain_sections.csv",
    "restain_registration_inputs": "restain_registration_inputs.csv",
    "pln_node_annotations": "pln_node_annotations.csv",
    "tissue_condition": "tissue_condition.csv",
}


def _path(project, name: str, default):
    """The project's path for a source, else the site's, else the checkout's tracked copy."""
    return project.source(name) or default


def inputs(project) -> list[Input]:
    cohort = profiles.load_all()[("cohort", project.cohort)]
    found = [
        Input("cohort_spec", _path(project, "cohort_spec", profiles.REPO / cohort.membership)),
        Input("qupath_manifest", _path(project, "qupath_manifest",
                                       roots.INTEGRATION.parent / "qupath_manifest.json")),
    ]
    found += [Input(n, _path(project, n, roots.INTEGRATION / f)) for n, f in _TRACKED.items()]
    for name, required in (("donor_notes", True), ("donor_notes_copy", False),
                           ("pln_blocks", True), ("donor_metadata", True)):
        found.append(Input(name, _needed(project, name), required=required))
    if project.source("npod_cases"):
        found.append(Input("npod_cases", project.source("npod_cases"), required=False))
    found += _instrument_inputs(project, found)
    return found


def _needed(project, name: str) -> Path:
    path = project.source(name)
    if path is None:
        raise relationships.RelationshipError(
            f"{STAGE} needs the {name} source; name it in the site file or the project")
    return Path(path)


def _xenium_outputs(project, table: Path) -> list[tuple[str, Path]]:
    """(output folder name, where this machine holds it) for every recorded region output."""
    root = project.settings.get("xenium_root")
    out: dict[str, Path] = {}
    for row in _rows(table):
        recorded = row["source_path"].strip().rstrip("/")
        if not recorded or recorded.upper() == "NA":
            continue
        local = Path(recorded) if Path(recorded).is_absolute() or not root else Path(root) / recorded
        out.setdefault(local.name, local)
    return sorted(out.items())


def _instrument_inputs(project, declared: list[Input]) -> list[Input]:
    """Every qptiff of both PhenoCycler projects, read for its first page only (bound by size,
    modification time and the digest the QuPath manifest records), and the three records of
    every Xenium region output: from the output itself where this machine holds it complete,
    else from the copy of its records."""
    paths = {i.name: Path(i.path) for i in declared}
    manifest = json.loads(paths["qupath_manifest"].read_text())
    recorded = {Path(i["qptiff"]["path"]).name: i["qptiff"].get("sha256")
                for i in manifest.get("images", []) if i.get("qptiff", {}).get("path")}
    out = []
    for folder in QPTIFF_FOLDERS:
        for path in sorted(_needed(project, folder).glob("*.qptiff")):
            out.append(Input(f"qptiff:{path.name}", path, hash=False,
                             recorded_sha256=recorded.get(path.name)))
    copies = _needed(project, "xenium_records")
    out.append(Input("xenium_records_manifest", copies / "MANIFEST.tsv", required=False))
    overviews: dict[str, list[Path]] = {}
    for name, local in _xenium_outputs(project, paths["xenium_paths"]):
        held = local if (local / XENIUM_RECORDS[0]).is_file() else copies / name
        out += [Input(f"xenium:{name}:{f}", held / f, required=False) for f in XENIUM_RECORDS]
        serial, _ = parse_bundle_path(name)
        if serial:
            overviews.setdefault(serial, []).extend(
                [local / "aux_outputs" / OVERVIEW, copies / name / OVERVIEW])
    for serial, candidates in sorted(overviews.items()):
        found = next((c for c in candidates if c.is_file()), candidates[-1])
        out.append(Input(f"overview:{serial}", found, required=False))
    return out


def settings(project) -> dict:
    return {"xenium_root": project.settings.get("xenium_root")}


def _rows(path: Path, delimiter: str = ",") -> list[dict]:
    with Path(path).open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader((line for line in handle if not line.startswith("#")),
                                   delimiter=delimiter))


def _stat(path) -> tuple:
    if path is None:
        return None, None, None
    path = Path(path)
    if not path.exists():
        return False, None, None
    stat = path.stat()
    return True, stat.st_size, stat.st_mtime_ns


def _instrument_records(task: Task) -> dict:
    """Read every declared qptiff's first page and every Xenium output's three records; what
    cannot be read is kept with the reason, never skipped."""
    headers, unreadable = {}, {}
    for name, item in task.inputs.items():
        if name.startswith("qptiff:"):
            try:
                headers[Path(item.path).name] = acquisition_records.read_qptiff_header(item.path)
            except (acquisition_records.AcquisitionRecordError, OSError) as err:
                unreadable[Path(item.path).name] = str(err)
    xenium, missing = {}, {}
    for name in sorted(n for n in task.inputs if n.startswith("xenium:")):
        _, output, _ = name.split(":", 2)
        if output in xenium or output in missing:
            continue
        files = {f: Path(task.inputs[f"xenium:{output}:{f}"].path) for f in XENIUM_RECORDS}
        absent = [f for f, path in files.items() if not path.is_file()]
        if absent:
            missing[output] = absent
            continue
        try:
            xenium[output] = acquisition_records.read_xenium_records(*files.values())
        except (acquisition_records.AcquisitionRecordError, OSError, KeyError, ValueError) as err:
            unreadable[output] = str(err)
    copies = {}
    manifest = task.input("xenium_records_manifest")
    if manifest.is_file():
        for row in _rows(manifest, "\t"):
            if row["file"] in XENIUM_RECORDS:          # an overview scan sits under aux_outputs
                copies[row["output"]] = row["source"].rsplit("/", 1)[0]
    return {"headers": headers, "xenium": xenium, "missing": missing, "unreadable": unreadable,
            "copies": copies}


def _acquisitions(task: Task, spec: dict, manifest: dict, aliases: dict,
                  records: dict) -> list[dict]:
    loaded = profiles.load_all()
    rules = {p: loaded[("platform", profile)].panel_names for p, profile in PANEL_PROFILES.items()}
    out = []
    spec_images = {i["image_id"] for i in spec.get("images", [])}
    panels = {p["panel_id"]: p for p in manifest.get("panels", [])}

    def pheno(file: str, panel_id=None) -> dict:
        header = records["headers"].get(file)
        if header is None:
            return dict(run_note=records["unreadable"].get(file, "no qptiff of that name"))
        fields = dict(scan_time=header.scan_time, experiment=header.experiment,
                      experiment_date=header.experiment_date, run_position=header.position,
                      run=f"phenocycler:{header.experiment}:{header.scan_date}",
                      run_date=header.scan_date, software=header.software, qptiff=file,
                      slide_id=header.slide_id)
        try:
            fields["panel"] = panel_names.phenocycler(header.markers, rules["phenocycler"])
        except panel_names.PanelNameError as err:
            fields["panel_error"] = str(err)
        if panel_id is not None:
            listed = {panel_names.marker_name(c.get("marker") or c.get("channel_name"))
                      for c in panels.get(panel_id, {}).get("channels", [])}
            mine = {panel_names.marker_name(m) for m in header.markers}
            fields["markers_match_manifest_panel"] = listed == mine
        return fields

    for index, image in enumerate(manifest.get("images", [])):
        record = image.get("qptiff") or {}
        exists, size, mtime = _stat(record.get("path"))
        out.append(dict(
            acquisition_id=f"phenocycler:{image['image_id']}", platform="phenocycler",
            tissue=PANCREAS, roi="panc", donor_id=aliases.get(image["donor_id"], image["donor_id"]),
            donor_label=image["donor_id"], source="qupath_manifest", source_row=index + 1,
            path=record.get("path"), exists=exists, bytes=size,
            recorded_bytes=record.get("size_bytes"),
            size_matches_record=(None if size is None or record.get("size_bytes") is None
                                 else size == record["size_bytes"]),
            mtime_matches_record=(None if mtime is None or record.get("mtime_ns") is None
                                  else mtime == record["mtime_ns"]),
            recorded_digest=record.get("sha256"), recorded_digest_mode=record.get("hash_mode"),
            panel_id=image.get("panel_id"), segmentation_version=image.get("segmentation_version"),
            in_cohort_spec=image["image_id"] in spec_images,
            **pheno(Path(record.get("path", "")).name, image.get("panel_id"))))
    for index, row in enumerate(_rows(task.input("pln_node_annotations")), start=1):
        out.append(dict(
            acquisition_id=f"phenocycler:{row['image_file']}", platform="phenocycler",
            tissue=LYMPH_NODE, roi=row["roi"],
            donor_id=aliases.get(row["donor_id"], row["donor_id"]), donor_label=row["donor_id"],
            source="pln_node_annotations", source_row=index, path=None, exists=None,
            n_nodes=row.get("n_nodes", ""),
            note=f"IO60pLN slide; excluded={row.get('excluded', '')}", **pheno(row["image_file"])))
    root = task.settings.get("xenium_root")
    for index, row in enumerate(_rows(task.input("xenium_paths")), start=1):
        serial, folder_region = parse_bundle_path(row["source_path"])
        recorded = row["source_path"].strip().rstrip("/")
        if recorded and recorded.upper() != "NA" and not Path(recorded).is_absolute() and root:
            path = str(Path(root) / recorded)
        else:
            path = recorded if recorded and recorded.upper() != "NA" else None
        exists, size, _ = _stat(path)
        tissue = LYMPH_NODE if row["roi"].startswith("pln") else PANCREAS
        output = Path(path).name if path else None
        record = records["xenium"].get(output)
        fields: dict = {}
        if record is not None:
            fields = dict(run=f"xenium:{record.run_name}", run_date=record.run_date,
                          run_start=record.run_start, recorded_slide=record.slide_id,
                          cassette_name=record.experiment.get("cassette_name", ""),
                          batch_number=_batch_number(record.run_name),
                          recorded_region=record.region, folder_region=folder_region,
                          local_copy=("complete" if Path(path, XENIUM_RECORDS[0]).is_file()
                                      else "partial" if exists else "not on this machine"),
                          records_read_from=("this machine's complete copy"
                                             if Path(path, XENIUM_RECORDS[0]).is_file() else
                                             f"copy of {records['copies'].get(output, 'a source MANIFEST.tsv does not name')}"))
            try:
                fields["panel"] = panel_names.xenium(record.experiment, rules["xenium"])
            except panel_names.PanelNameError as err:
                fields["panel_error"] = str(err)
        region = record.region if record is not None else folder_region
        out.append(dict(
            acquisition_id=(f"xenium:{serial}:{region}" if serial
                            else f"xenium:unrecorded:{row['donor_id']}:{row['roi']}:{index}"),
            platform="xenium", tissue=tissue, roi=row["roi"],
            donor_id=aliases.get(row["donor_id"], row["donor_id"]), donor_label=row["donor_id"],
            source="xenium_paths", source_row=index, path=path, exists=exists, output=output,
            note=f"batch {row.get('batch', '')}" + ("" if serial else "; no bundle recorded"),
            **fields))
    sections = {r["image_file"]: r for r in _rows(task.input("restain_sections"))}
    for index, row in enumerate(_rows(task.input("restain_registration_inputs")), start=1):
        roles = sections.get(row["image_file"])
        recorded = [roles[c] for c in roles if c.startswith("channel_") and roles[c]] if roles else []
        fields = {}
        try:
            fields["panel"] = panel_names.post_xenium_if(rules["if"], recorded or None)
        except panel_names.PanelNameError as err:
            fields["panel_error"] = str(err)
        out.append(dict(
            acquisition_id=f"if:{row['image_file']}", platform="if", tissue=PANCREAS, **fields,
            roi=row["roi"], donor_id=aliases.get(row["donor_id"], row["donor_id"]),
            donor_label=row["donor_id"], source="restain_registration_inputs", source_row=index,
            path=None, exists=None, xenium_serial=row.get("xenium_serial"),
            note=("channel roles recorded in restain_sections.csv" if recorded
                  else "no channel roles recorded; the panel is the one every image carries")))
    return out


def _batch_number(run_name: str) -> str:
    """The batch number a Xenium run name carries -- MM035_Xenium-npodBatch1 is batch 1,
    MM48_batch14 batch 14 -- or "" for a run whose name carries none (the PhenoXeno runs)."""
    m = re.search(r"[Bb]atch[_-]?(\d+)", run_name or "")
    return str(int(m.group(1))) if m else ""


def _statuses(task: Task, cohort, donors: set) -> tuple[dict, list[str], list[str]]:
    """Each donor's disease status, as a level of the cohort profile's disease_status: from the
    project's donor metadata, else from nPOD's case export. Returns ({donor: (status, source)},
    the donors whose two sources disagree, the recorded spellings no level names)."""
    import pandas as pd

    if not cohort.disease_status:
        raise profiles.ProfileError(f"cohort profile {cohort.id} declares no disease_status; "
                                    f"{STAGE} reads each donor's status in its spellings")
    spellings = cohort.disease_status["spellings"]
    found: dict[str, tuple[str, str]] = {}
    disagree, unknown = [], []
    for source, id_column, status_column in STATUS_SOURCES:
        if source not in task.inputs:
            continue
        frame = pd.read_excel(task.input(source), dtype=str)
        missing = [c for c in (id_column, status_column) if c not in frame.columns]
        if missing:
            raise ValueError(f"{task.input(source)} has no {' or '.join(missing)} column")
        for label, value in zip(frame[id_column], frame[status_column], strict=True):
            if pd.isna(label):
                continue
            donor = cohort.donor_aliases.get(str(label).strip(), str(label).strip())
            if donor not in donors or pd.isna(value):
                continue
            level = spellings.get(str(value).strip())
            if level is None:
                unknown.append(f"{source} {donor}: {str(value).strip()!r}")
            elif donor not in found:
                found[donor] = (level, source)
            elif found[donor][0] != level:
                disagree.append(f"{donor}: {found[donor][1]} {found[donor][0]}, {source} {level}")
    return found, disagree, unknown


def _thumbnail(source: tuple, task: Task):
    """A small RGB array of one slide: a qptiff's own thumbnail page, or a Xenium slide's
    overview scan; None where the slide's image is not held or holds none."""
    import numpy as np
    from PIL import Image

    kind, name = source
    try:
        if kind == "qptiff":
            import tifffile

            with tifffile.TiffFile(task.input(f"qptiff:{name}")) as tf:
                page = next((pg for pg in tf.pages
                             if "<ImageType>Thumbnail</ImageType>" in (pg.description or "")), None)
                if page is None:
                    return None
                image = Image.fromarray(page.asarray())
        else:
            path = task.input(f"overview:{name}")
            if not Path(path).is_file():
                return None
            image = Image.open(path)
            image.load()
    except (KeyError, OSError, ValueError):
        return None
    image = image.convert("RGB")
    image.thumbnail((THUMBNAIL_PIXELS, THUMBNAIL_PIXELS))
    return np.asarray(image)


def _thumbnails(task: Task, cells: dict) -> tuple[dict, list[str]]:
    """{(donor, slide group): [thumbnail, ...]} for every slide cell with an image source -- a
    PhenoCycler slide that Xenium also ran on has both, PhenoCycler first -- and the images that
    are not held."""
    thumbs, missing = {}, []
    for (donor, group, column), cell in sorted(cells.items()):
        if column != SLIDE_NUMBER:
            continue
        for source in cell.get("images", []):
            image = _thumbnail(source, task)
            if image is None:
                missing.append(f"{donor} {group} ({source[1]})")
            else:
                thumbs.setdefault((donor, group), []).append(image)
        found = thumbs.get((donor, group), [])
        if len(found) == 2:      # a PhenoCycler slide's own image, then its Xenium overview
            found[1] = _oriented_like(found[0], found[1])
    return thumbs, missing


def _tissue_outline(image, size: int = 48):
    """An image's tissue, cropped to it and resampled to a square: bright against a dark slide."""
    import numpy as np
    from PIL import Image

    grey = np.asarray(Image.fromarray(image).convert("L"), dtype=float)
    mask = grey > grey.mean()
    rows, cols = np.flatnonzero(mask.any(axis=1)), np.flatnonzero(mask.any(axis=0))
    if not len(rows) or not len(cols):
        return None
    crop = mask[rows[0]:rows[-1] + 1, cols[0]:cols[-1] + 1].astype(np.uint8) * 255
    return np.asarray(Image.fromarray(crop).resize((size, size)), dtype=float) / 255


def _oriented_like(reference, image):
    """``image`` turned the way ``reference`` is. The two instruments hold a slide opposite ways
    round, so a same-slide Xenium overview can come out turned 180 degrees from the PhenoCycler's
    image of that slide; it is turned when that matches the reference's tissue better."""
    import numpy as np

    want, have = _tissue_outline(reference), _tissue_outline(image)
    if want is None or have is None:
        return image

    def agree(a, b):
        return float(np.corrcoef(a.ravel(), b.ravel())[0, 1])

    return image[::-1, ::-1] if agree(want, have[::-1, ::-1]) > agree(want, have) else image


def _label_carries_donor(label: str, a: dict) -> bool:
    """Whether a label typed at an instrument -- a Xenium cassette name, a PhenoCycler SlideID --
    carries the listed donor's number: HDL112, nPOD6591 and 0080579_6347_pLNpanc carry 112, 6591
    and 6347. It is a label, so it is compared, never used to name the donor."""
    numbers = {int(n) for n in re.findall(r"\d+", label or "")}
    return any(int(n) in numbers for donor in (a["donor_id"], a.get("donor_label", ""))
               for n in re.findall(r"\d+", donor))


def _slide_entry(s: dict) -> str:
    """One slide of a run as the runs sheet writes it: whose, and its slide number -- a Xenium
    slide's serial, which a PhenoCycler slide has where Xenium also ran on it. A PhenoCycler
    experiment's slide positions (Slide1, Slide2) are not slide numbers (maintainer, 2026-10-23)."""
    return s["label"] + (f" (slide {s['slide']})" if s["slide"] else "")


def _conditions(acquisitions: list[dict], rows: list[dict]) -> list[str]:
    """Mark each acquisition the maintainer recorded damage on (tissue_condition.csv), with the
    damaged tissues in the tissue profiles' names; return the rows that name no acquisition or a
    tissue no profile knows."""
    by_id = {a["acquisition_id"]: a for a in acquisitions}
    unmatched = []
    for row in rows:
        a = by_id.get(row["acquisition_id"])
        try:
            damaged = [tissues.normalize(t) for t in row["tissue"].split(";")]
        except tissues.TissueError as err:
            unmatched.append(f"{row['acquisition_id']} ({row.get('donor_id', '')}): {err}")
            continue
        if a is None:
            unmatched.append(f"{row['acquisition_id']} ({row.get('donor_id', '')})")
            continue
        a["condition"] = row["condition"]
        a["condition_when"] = row["when"]
        a["condition_tissues"] = "; ".join(damaged)
        a["known_condition"] = (f"{' and '.join(t.strip() for t in row['tissue'].split(';'))} "
                                f"{row['condition']} {row['when']}")
        a["condition_detail"] = row.get("detail", "")
    return unmatched


def _slide_label(acquisitions: list[dict]) -> str:
    """'6482 pancreas', or '6539 pancreas + pLN' for a Xenium slide with two tissues, with any
    recorded damage after it."""
    tissues = []
    for a in acquisitions:
        name = "pancreas" if a["tissue"] == PANCREAS else "pLN"
        if name not in tissues:
            tissues.append(name)
    damage = sorted({a["known_condition"] for a in acquisitions if a.get("known_condition")})
    return (f"{acquisitions[0]['donor_id']} {' + '.join(tissues)}"
            + (f" [{'; '.join(damage)}]" if damage else ""))


def _runs(acquisitions: list[dict], records: dict) -> list[dict]:
    """Every run of either platform, the slides run together in it and its panels. A PhenoCycler
    partner that is no acquisition here (a slide of a donor outside the cohort) is named by its
    qptiff and kept."""
    by_file = {a["qptiff"]: a for a in acquisitions if a.get("qptiff")}
    out = []
    groups = acquisition_records.phenocycler_runs(list(records["headers"].values()))
    for (experiment, _, day), headers in groups.items():
        slides = []
        for h in headers:
            name = Path(h.path).name
            acq = by_file.get(name)
            slides.append(dict(slide=acq.get("xenium_slide", "") if acq else "",
                               position=h.position, scan_time=h.scan_time,
                               donor=acq["donor_id"] if acq else h.slide_id,
                               label=(_slide_label([acq]) if acq else
                                      f"{h.slide_id} ({name}, no acquisition here)"),
                               acquisitions=[acq["acquisition_id"]] if acq else [],
                               panel=acq.get("panel", "") if acq else ""))
        out.append(dict(platform="phenocycler", run=f"phenocycler:{experiment}:{day}",
                        run_name=experiment, batch="", date=day, slides=slides))
    xenium: dict[str, list[dict]] = {}
    for a in acquisitions:
        if a["platform"] == "xenium" and a.get("run"):
            xenium.setdefault(a["run"], []).append(a)
    for run, members in xenium.items():
        by_slide: dict[str, list[dict]] = {}
        for a in members:
            by_slide.setdefault(a["recorded_slide"], []).append(a)
        slides = [dict(slide=s, position="", scan_time="", label=_slide_label(group),
                       donor=group[0]["donor_id"],
                       acquisitions=[a["acquisition_id"] for a in group],
                       panel="; ".join(sorted({a.get("panel", "") for a in group})))
                  for s, group in sorted(by_slide.items())]
        out.append(dict(platform="xenium", run=run, run_name=run.split(":", 1)[1],
                        batch=_batch_number(run.split(":", 1)[1]),
                        date=min(a["run_date"] for a in members), slides=slides))
    out.sort(key=lambda r: (r["platform"] != "phenocycler", r["date"], r["run"]))
    partners, run_with = {}, {}
    for r in out:
        for s in r["slides"]:
            others = [o for o in r["slides"] if o is not s]
            for acq in s["acquisitions"]:
                partners[acq] = "; ".join(_slide_entry(o) for o in others) or "run alone"
                run_with[acq] = ", ".join(o["donor"] for o in others)
    for a in acquisitions:
        if a["acquisition_id"] in partners:
            a["run_partners"] = partners[a["acquisition_id"]]
            a["run_with"] = run_with[a["acquisition_id"]]
    return out


def _blocks(acquisitions: list[dict], reading) -> None:
    """The block each acquisition's section was cut from, as the maintainer stated it (2026-10-21):
    the pancreas sections from the pancreas block; "All the phenocycler dedicated pLN slides are
    from a block" (a pLN block); "Xenium dedicated slides mostly put a pancreas and pLN section from
    different blocks", the exception being read from the donor notes "for the phenocycler slides
    that were sequential section": where the PhenoCycler slide was a sequential section of the
    Xenium slide and its line gives same-block pLN, the Xenium pLN section is that same-block pLN.
    A sequential line that gives both kinds leaves the block unsettled (blank), and the gate
    blocks on it."""
    for a in acquisitions:
        if a["platform"] == "if":
            continue
        if a["tissue"] == PANCREAS:
            a["block"] = PANCREAS_BLOCK
        elif a["platform"] == "phenocycler":
            a["block"] = OTHER_BLOCK          # a pLN block
        else:
            note = reading.donors.get(a["donor_id"])
            if note is None or note.xenium_design != relationships.SEQUENTIAL:
                a["block"] = OTHER_BLOCK
                continue
            same, different = note.same_block_pln or 0, note.different_block_pln or 0
            a["block"] = "" if same and different else PANCREAS_BLOCK if same else OTHER_BLOCK


def _donor_order(donor_id: str) -> tuple:
    return (not donor_id.isdigit(), donor_id.zfill(8))


def _xenium_slides(acquisitions: list[dict], reading) -> None:
    """A PhenoCycler pancreas slide that is also a Xenium slide (section_modes.csv same_slide) takes
    that slide's serial as its number; a donor with one such slide and one PhenoCycler scan only,
    as the slide grouping reads it."""
    serials: dict[str, set] = {}
    for m in reading.modes:
        if m.mode == "same_slide" and m.xenium_serial:
            serials.setdefault(m.donor_id, set()).add(m.xenium_serial)
    scans = Counter(a["donor_id"] for a in acquisitions
                    if a["platform"] == "phenocycler" and a["tissue"] == PANCREAS)
    for a in acquisitions:
        found = serials.get(a["donor_id"], set())
        if (a["platform"] == "phenocycler" and a["tissue"] == PANCREAS and len(found) == 1
                and scans[a["donor_id"]] == 1):
            a["xenium_slide"] = next(iter(found))


def _xenium_section(a: dict) -> str:
    """The section column a Xenium region is counted in: its pancreas, or its pLN by block."""
    if a["tissue"] == PANCREAS:
        return PANCREAS_SECTION
    return SAME_BLOCK if a.get("block") == PANCREAS_BLOCK else DIFFERENT_BLOCK


def _section_tissue(column: str) -> str:
    return PANCREAS if column == PANCREAS_SECTION else LYMPH_NODE


def _run_text(a: dict) -> str:
    """'2025-09-13 with 6551': the run's date and the donor of each slide run with it; 'alone'
    where nothing else was."""
    if not a.get("run_date"):
        return ""
    return f"{a['run_date']} " + (f"with {a['run_with']}" if a.get("run_with") else "alone")


def _acquisition_text(a: dict) -> str:
    """What an acquisition's cell says: its panel and its run. A post-Xenium IF image has no
    recorded run; its panel's channels go two to a line."""
    if a["platform"] == "if":
        channels = a.get("panel", "").split(", ")
        return "\n".join(", ".join(channels[i:i + 2]) for i in range(0, len(channels), 2))
    batch = f"batch {a['batch_number']}" if a.get("batch_number") else ""
    return "\n".join(part for part in (a.get("panel", ""), batch, _run_text(a)) if part)


def _coverage(acquisitions: list[dict], reading, spec, statuses: dict) -> dict:
    """The coverage, by physical slide.

    ``cells`` maps (donor, slide, column) to what that cell holds: its ``text``, the ``kind`` that
    colours it (a platform for an acquisition, a block for a section), whether recorded damage
    falls on it, and whether it stands for something a table records with no output listed.
    ``columns`` runs left to right as (slide, column, "section" or "acquisition"); ``groups``
    holds the donors by the Xenium design their notes give; ``table`` is coverage.csv, one row
    per donor; ``totals`` counts datasets by where they came from, centrally excluded donors not
    counted. A Xenium slide is one dataset however many regions it has; the IF scanned after
    Xenium is on that Xenium slide (its ``xenium_serial``)."""
    same_slide = {(m.donor_id, m.xenium_serial): m for m in reading.modes
                  if m.mode == "same_slide" and m.xenium_serial}
    scans: dict[str, list[str]] = {}
    for a in acquisitions:
        if a["platform"] == "phenocycler" and a["tissue"] == PANCREAS:
            scans.setdefault(a["donor_id"], []).append(a["acquisition_id"])

    def slide(a: dict):
        """(the slide group, a key for the physical slide) an acquisition is on, or None."""
        if a["platform"] == "phenocycler":
            return (PC_SLIDE if a["tissue"] == PANCREAS else PLN_SLIDE), a["acquisition_id"]
        if a["platform"] == "if":
            serial = a.get("xenium_serial")
        elif a["acquisition_id"].startswith("xenium:unrecorded:"):
            return None
        else:
            serial = a.get("recorded_slide") or a["acquisition_id"].split(":")[1]
        if not serial:
            return None
        if (a["donor_id"], serial) in same_slide:
            mine = scans.get(a["donor_id"], [])
            return PC_SLIDE, mine[0] if len(mine) == 1 else serial
        return XENIUM_SLIDE, serial

    cells: dict[tuple, dict] = {}
    datasets: dict[str, set] = {}
    slides: dict[str, set] = {}
    on_slide: dict[tuple, set] = {}
    regions: dict[tuple, set] = {}       # (donor, block column) -> its Xenium slide's pLN regions
    unlisted: list[tuple] = []           # a region xenium_paths.csv lists with no output
    unsettled: set[str] = set()

    def put(donor, group, column, text, kind, missing=False):
        cell = cells.setdefault((donor, group, column),
                                dict(parts=[], kind=kind, damaged=False, missing=missing))
        if text and text not in cell["parts"]:     # a slide's regions share one panel and run
            cell["parts"].append(text)

    for a in acquisitions:
        donor, where = a["donor_id"], slide(a)
        if a["platform"] == "xenium" and where is None:
            unlisted.append((donor, _xenium_section(a)))
            continue
        group, key = where or (XENIUM_SLIDE, None)
        if where is not None:
            slides.setdefault(donor, set()).add(where)
            on_slide.setdefault((donor, where), set()).add(a["platform"])
        if a["platform"] == "xenium":
            put(donor, group, XENIUM, _acquisition_text(a), XENIUM)
            datasets.setdefault(donor, set()).add((group, XENIUM, key))
            serial = a.get("recorded_slide") or a["acquisition_id"].split(":")[1]
            numbered = cells.get((donor, group, SLIDE_NUMBER))
            if serial and not (numbered and any(serial in part for part in numbered["parts"])):
                put(donor, group, SLIDE_NUMBER, serial, None)
            if serial:      # a same-slide Xenium's overview joins its PhenoCycler slide's image
                images = cells[(donor, group, SLIDE_NUMBER)].setdefault("images", [])
                if ("overview", serial) not in images:
                    images.append(("overview", serial))
            if group == XENIUM_SLIDE:
                if a["tissue"] == PANCREAS:
                    put(donor, group, PANCREAS_SECTION, "", PANCREAS)
                elif not a.get("block"):
                    unsettled.add(donor)
                else:
                    regions.setdefault((donor, _xenium_section(a)), set()).add(a["acquisition_id"])
            continue
        column = PHENOCYCLER if a["platform"] == "phenocycler" else POST_XENIUM_IF
        put(donor, group, column, _acquisition_text(a), column)
        if a["platform"] == "phenocycler" and a.get("slide_id"):
            put(donor, group, SLIDE_NUMBER, a.get("xenium_slide") or a["slide_id"], None)
            cells[(donor, group, SLIDE_NUMBER)].setdefault("images", []).insert(
                0, ("qptiff", a["qptiff"]))
        datasets.setdefault(donor, set()).add((group, column, a["acquisition_id"]))
        if group == PLN_SLIDE:
            put(donor, group, PLN_SECTION, a.get("n_nodes", ""), LYMPH_NODE)

    donors = ({a["donor_id"] for a in acquisitions} | set(reading.donors)
              | {ann.donor_id for ann in reading.pln})
    for donor in donors:
        note = reading.donors.get(donor)
        if note is not None or donor in scans:
            put(donor, PC_SLIDE, PANCREAS_SECTION, "", PANCREAS)
        for column, _, field in DESIGNATIONS:
            count = getattr(note, field) if note is not None else None
            if count:
                put(donor, PC_SLIDE, column, str(count), LYMPH_NODE)
    for (donor, column), found in regions.items():
        put(donor, XENIUM_SLIDE, column, str(len(found)), LYMPH_NODE)
    for donor, column in unlisted:
        if (donor, XENIUM_SLIDE, column) not in cells:
            put(donor, XENIUM_SLIDE, column, NOT_LISTED, None, missing=True)
    for donor in unsettled:
        put(donor, XENIUM_SLIDE, SAME_BLOCK, UNSETTLED, None, missing=True)
    listed = {(a["donor_id"], a.get("recorded_slide")) for a in acquisitions
              if a["platform"] == "xenium"}
    not_listed = []
    for (donor, serial), mode in sorted(same_slide.items()):
        if (donor, serial) not in listed:
            put(donor, PC_SLIDE, XENIUM, NOT_LISTED, None, missing=True)
            not_listed.append(f"{donor}: slide {serial} (section_modes.csv: {mode.evidence_text})")
    for cell in cells.values():
        cell["text"] = "\n".join(cell.pop("parts"))

    # recorded damage falls on the sections of the damaged tissue on that acquisition's slide
    when = set()
    for a in acquisitions:
        if not a.get("condition_tissues"):
            continue
        when.add(a.get("condition_when", ""))
        group = (slide(a) or (XENIUM_SLIDE, None))[0]
        damaged = set(a["condition_tissues"].split("; "))
        for (donor, g, column), cell in cells.items():
            if (donor == a["donor_id"] and g == group and column not in NOT_SECTIONS
                    and _section_tissue(column) in damaged and not cell["missing"]):
                cell["damaged"] = True

    columns = [(PC_SLIDE, SLIDE_NUMBER, "slide"),
               (PC_SLIDE, PANCREAS_SECTION, "section"), (PC_SLIDE, SAME_BLOCK, "section"),
               (PC_SLIDE, DIFFERENT_BLOCK, "section"), (PC_SLIDE, PHENOCYCLER, "acquisition"),
               (PC_SLIDE, XENIUM, "acquisition")]
    if any(g == PC_SLIDE and c == POST_XENIUM_IF for (_, g, c) in cells):
        columns.append((PC_SLIDE, POST_XENIUM_IF, "acquisition"))
    columns += [(XENIUM_SLIDE, SLIDE_NUMBER, "slide"), (XENIUM_SLIDE, PANCREAS_SECTION, "section"),
                (XENIUM_SLIDE, SAME_BLOCK, "section"), (XENIUM_SLIDE, DIFFERENT_BLOCK, "section"),
                (XENIUM_SLIDE, XENIUM, "acquisition"), (XENIUM_SLIDE, POST_XENIUM_IF, "acquisition"),
                (PLN_SLIDE, SLIDE_NUMBER, "slide"), (PLN_SLIDE, PLN_SECTION, "section"),
                (PLN_SLIDE, PHENOCYCLER, "acquisition")]

    order = [*DESIGNS.values(), XENIUM_ONLY]
    pheno = {a["donor_id"] for a in acquisitions if a["platform"] == "phenocycler"}
    grouped: dict[str, list[str]] = {}
    for donor in sorted(donors, key=_donor_order):
        note = reading.donors.get(donor)
        label = (DESIGNS.get(note.xenium_design, "no Xenium design in the notes") if note
                 else XENIUM_ONLY if donor not in pheno else "no Xenium design in the notes")
        grouped.setdefault(label, []).append(donor)
    groups = sorted(grouped.items(),
                    key=lambda kv: order.index(kv[0]) if kv[0] in order else len(order))

    annotated: Counter = Counter()
    for ann in reading.pln:
        annotated[(ann.donor_id, ann.block)] += ann.pieces
    cohort = set(map(str, spec.get("expected_donors", [])))
    table = []
    for label, members in groups:
        for donor in members:
            note = reading.donors.get(donor)
            status, source = statuses.get(donor, ("", ""))
            row = dict(donor_id=donor, design=label, disease_status=status,
                       disease_status_source=source, in_cohort_spec=donor in cohort,
                       centrally_excluded=donor in DONOR_EXCLUSIONS,
                       excluded_platforms=", ".join(p for p in PLATFORM_NAMES
                                                    if excluded_on(donor, p)))
            for group, column, kind in columns:
                row[f"{group}: {column}"] = _cell_text(cells.get((donor, group, column)), column, kind)
            row[DATASETS] = len(datasets.get(donor, ()))
            row[SLIDES] = len(slides.get(donor, ()))
            for column, block, _ in DESIGNATIONS:
                row[column + ANNOTATED] = annotated[(donor, block)]
            row["designation read from free text"] = bool(note is not None and note.pln_from_caution)
            row["caution in the notes"] = bool(note is not None and note.cautions)
            table.append(row)

    # an exclusion takes out only the platform it covers (cohort.EXCLUDED_PLATFORMS): 6457's
    # PhenoCycler data, never its Xenium
    per_column, not_counted = Counter(), Counter()
    for d in donors:
        for g, c, _ in datasets.get(d, ()):
            if excluded_on(d, PLATFORM_OF[c]):
                not_counted[d] += 1
            else:
                per_column[(g, c)] += 1
    totals = dict(per_column={f"{g}: {c}": per_column[(g, c)] for g, c, kind in columns
                              if kind == "acquisition"},
                  datasets=sum(per_column.values()),
                  slides=sum(1 for (d, _), platforms in on_slide.items()
                             if any(not excluded_on(d, p) for p in platforms)),
                  not_counted=dict(sorted(not_counted.items())))
    return dict(cells=cells, columns=columns, groups=groups, table=table, totals=totals,
                table_by_donor={r["donor_id"]: r for r in table}, not_listed=not_listed,
                damage_when=sorted(w for w in when if w))


def _cell_text(cell: dict | None, column: str, kind: str) -> str:
    """A coverage.csv cell: an acquisition's panel and run; a section's block (the notes' count,
    for the pLN/lymphoid on the PhenoCycler slide); "(damaged)" where recorded damage falls."""
    if cell is None:
        return ""
    if cell["missing"]:
        return cell["text"]
    if kind == "slide":
        return cell["text"].replace("\n", "; ")
    if kind == "acquisition":       # an IF panel's channels wrap; a panel and its run are two parts
        text = cell["text"].replace("\n", ", " if column == POST_XENIUM_IF else "; ")
    elif column in (SAME_BLOCK, DIFFERENT_BLOCK, PLN_SECTION):
        text = f"{cell['text']} pLN" if cell["text"] else "pLN"
    else:
        text = cell["kind"]
    return text + (" (damaged)" if cell["damaged"] else "")


def _designations(reading) -> list[str]:
    """Each pancreas the notes give any pLN/lymphoid -- or leave a block unsaid -- as written."""
    out = []
    for note in sorted(reading.donors.values(), key=lambda n: _donor_order(n.donor_id)):
        said = []
        for _, block, field in DESIGNATIONS:
            name = block.replace("_", "-")
            count = getattr(note, field)
            if count is None:
                said.append(f"{name} not stated")
            elif count:
                said.append(f"{count} {name}")
        if said:
            out.append(f"{note.donor_id}: {', '.join(said)}"
                       + (" (read from the line's free text)" if note.pln_from_caution else ""))
    return out


def _write_csv(path: Path, rows: list[dict]) -> None:
    columns: list[str] = []
    for row in rows:
        columns += [k for k in row if k not in columns]
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns or ["empty"])
        writer.writeheader()
        writer.writerows(rows)


def run(task: Task) -> dict:
    project = task.project
    cohort = profiles.load_all()[("cohort", project.cohort)]
    spec = json.loads(task.input("cohort_spec").read_text())
    manifest = json.loads(task.input("qupath_manifest").read_text())
    copy = task.input("donor_notes_copy")
    reading = relationships.load(
        notes=task.input("donor_notes"), notes_copy=copy if copy.exists() else False,
        pln_blocks=task.input("pln_blocks"), section_modes=task.input("section_modes"),
        aliases=cohort.donor_aliases, cohort_donors=spec.get("expected_donors", []))
    records = _instrument_records(task)
    acquisitions = _acquisitions(task, spec, manifest, cohort.donor_aliases, records)
    unmatched = _conditions(acquisitions, _rows(task.input("tissue_condition")))
    _blocks(acquisitions, reading)
    _xenium_slides(acquisitions, reading)
    runs = _runs(acquisitions, records)
    donors = ({a["donor_id"] for a in acquisitions} | set(reading.donors)
              | {ann.donor_id for ann in reading.pln})
    statuses, status_disagree, status_unknown = _statuses(task, cohort, donors)
    coverage = _coverage(acquisitions, reading, spec, statuses)
    coverage["thumbs"], unthumbed = _thumbnails(task, coverage["cells"])
    _write_csv(task.output("acquisitions.csv"), acquisitions)
    _write_csv(task.output("coverage.csv"), coverage["table"])
    _write_csv(task.output("runs.csv"), [dict(
        platform=r["platform"], run=r["run"], run_name=r["run_name"], batch=r["batch"],
        date=r["date"],
        slides=len(r["slides"]), panels="; ".join(sorted({s["panel"] for s in r["slides"]} - {""})),
        slide_numbers="; ".join(s["slide"] for s in r["slides"] if s["slide"]),
        slides_run_together=" | ".join(_slide_entry(s) for s in r["slides"]),
        acquisitions="; ".join(a for s in r["slides"] for a in s["acquisitions"]))
        for r in runs])
    _write_csv(task.output("donors.csv"), [dict(
        donor_id=d.donor_id, label=d.label, pancreas_region=d.pancreas_region,
        same_block_pln=d.same_block_pln, different_block_pln=d.different_block_pln,
        pln_from_free_text=d.pln_from_caution, xenium_design=d.xenium_design,
        post_xenium_if=d.post_xenium_if, cautions="; ".join(d.cautions), source=d.evidence.where)
        for d in reading.donors.values()])
    _write_csv(task.output("relationships.csv"), [dict(
        kind=r.kind, donor_id=r.donor_id, a=r.a, b=r.b, reason=r.reason,
        evidence="; ".join(f"{e.source} {e.where}" for e in r.evidence))
        for r in reading.relationships])
    known = relationships.KNOWN_DISAGREEMENTS
    disagreements = [dict(donor_id=c.donor_id, what=c.what, detail=c.detail,
                          known_reason=known.get((c.donor_id, c.what), ""))
                     for c in reading.conflicts]
    _write_csv(task.output("disagreements.csv"), disagreements)
    _figure(task, coverage)
    _runs_figure(task, runs)
    platforms = Counter(a["platform"] for a in acquisitions)
    pheno = [a for a in acquisitions if a["platform"] == "phenocycler"]
    xenium = [a for a in acquisitions if a["platform"] == "xenium" and a.get("output")]
    cassettes_differ = [a for a in xenium if a.get("cassette_name")
                        and not _label_carries_donor(a["cassette_name"], a)]
    slide_ids_differ = [a for a in pheno if a.get("slide_id")
                        and not _label_carries_donor(a["slide_id"], a)]
    ids = Counter(a["acquisition_id"] for a in acquisitions)
    spec_images = [i["image_id"] for i in spec.get("images", [])]
    in_manifest = {i["image_id"] for i in manifest.get("images", [])}
    xenium_cell = f"{XENIUM_SLIDE}: {XENIUM}"
    return {
        "acquisitions": dict(platforms), "donors_in_notes": len(reading.donors),
        "relationships": dict(Counter(r.kind for r in reading.relationships)),
        "repeated_acquisition_ids": sorted(k for k, n in ids.items() if n > 1),
        "cohort_donors_without_notes": sorted(set(map(str, spec.get("expected_donors", [])))
                                              - set(reading.donors)),
        "notes_donors_not_in_cohort": sorted(set(reading.donors)
                                             - set(map(str, spec.get("expected_donors", [])))),
        "cohort_images_not_in_manifest": sorted(set(spec_images) - in_manifest),
        "unknown_disagreements": [f"{d['donor_id']}: {d['what']}" for d in disagreements
                                  if not d["known_reason"]],
        "known_disagreements": [f"{d['donor_id']}: {d['what']} ({d['known_reason']})"
                                for d in disagreements if d["known_reason"]],
        "raw_size_differs": sorted(a["acquisition_id"] for a in acquisitions
                                   if a.get("size_matches_record") is False),
        "raw_mtime_differs": sorted(a["acquisition_id"] for a in acquisitions
                                    if a.get("mtime_matches_record") is False),
        "raw_not_found": sorted(a["acquisition_id"] for a in acquisitions
                                if a.get("exists") is False),
        "pln_designations": _designations(reading),
        "cautions": [f"{n.donor_id}: {'; '.join(n.cautions)}"
                     for n in sorted(reading.donors.values(), key=lambda n: _donor_order(n.donor_id))
                     if n.cautions],
        "qptiffs_read": len(records["headers"]),
        "phenocycler_unread": sorted(f"{a['acquisition_id']}: {a.get('run_note', '')}"
                                     for a in pheno if not a.get("run")),
        "qptiffs_unreadable": sorted(f"{k}: {v}" for k, v in records["unreadable"].items()
                                     if k.endswith(".qptiff")),
        "markers_not_manifest_panel": sorted(a["acquisition_id"] for a in pheno
                                             if a.get("markers_match_manifest_panel") is False),
        "panels_unnamed": sorted(f"{a['acquisition_id']}: {a['panel_error']}"
                                 for a in acquisitions if a.get("panel_error")),
        "panels": {p: dict(sorted(Counter(a["panel"] for a in acquisitions
                                          if a["platform"] == p and a.get("panel")).items()))
                   for p in PANEL_PROFILES},
        "xenium_records_missing": sorted(f"{k}: {', '.join(v)}"
                                         for k, v in records["missing"].items()),
        "xenium_records_unreadable": sorted(f"{k}: {v}" for k, v in records["unreadable"].items()
                                            if not k.endswith(".qptiff")),
        "xenium_records_disagree": sorted(
            [f"{k}: {d}" for k, r in records["xenium"].items() for d in r.differences]
            + [f"{a['output']}: experiment.xenium slide {a['recorded_slide']} != the folder's "
               f"{a['acquisition_id'].split(':')[1]}" for a in xenium
               if a.get("recorded_slide") and a["recorded_slide"] != a["acquisition_id"].split(":")[1]]),
        "xenium_folders_renamed": sorted(
            f"{a['output']}: experiment.xenium names region {a['recorded_region']}"
            for a in xenium if a.get("recorded_region") is not None
            and a["recorded_region"] != a["folder_region"] and a["output"] not in KNOWN_RENAMED_FOLDERS),
        "xenium_folders_renamed_known": sorted(
            f"{a['output']} ({KNOWN_RENAMED_FOLDERS[a['output']]})" for a in xenium
            if a.get("output") in KNOWN_RENAMED_FOLDERS),
        "slides_without_thumbnail": unthumbed,
        "xenium_cassettes_differ": sorted(
            f"{a['output']}: cassette {a['cassette_name']}, listed as donor {a['donor_label']}"
            for a in cassettes_differ if a["recorded_slide"] not in KNOWN_LABEL_ERRORS),
        "phenocycler_slide_labels_differ": sorted(
            f"{a['qptiff']}: SlideID {a['slide_id']}, listed as donor {a['donor_label']}"
            for a in slide_ids_differ if a["qptiff"] not in KNOWN_LABEL_ERRORS),
        "labels_differ_known": sorted(
            {f"slide {a['recorded_slide']}: cassette {a['cassette_name']}, listed as donor "
             f"{a['donor_label']} ({KNOWN_LABEL_ERRORS[a['recorded_slide']]})"
             for a in cassettes_differ if a["recorded_slide"] in KNOWN_LABEL_ERRORS}
            | {f"{a['qptiff']}: SlideID {a['slide_id']}, listed as donor {a['donor_label']} "
               f"({KNOWN_LABEL_ERRORS[a['qptiff']]})"
               for a in slide_ids_differ if a["qptiff"] in KNOWN_LABEL_ERRORS}),
        "xenium_held_in_part": sorted(f"{a['output']} ({a['local_copy']})" for a in xenium
                                      if a.get("local_copy") and a["local_copy"] != "complete"),
        "runs": {p: sum(r["platform"] == p for r in runs) for p in ("phenocycler", "xenium")},
        "runs_over_two_slides": sorted(f"{r['run']}: {len(r['slides'])} slides" for r in runs
                                       if len(r["slides"]) > 2),
        "slides_run_alone": sorted(f"{r['slides'][0]['label']}, {r['date']}" for r in runs
                                   if r["platform"] == "phenocycler" and len(r["slides"]) == 1),
        "conditions_unmatched": unmatched,
        "donors_without_status": sorted(d for d in donors if d not in statuses) + status_unknown,
        "status_disagreements": status_disagree,
        "disease_status": {level: sorted((d for d, (s, _) in statuses.items() if s == level),
                                         key=_donor_order)
                           for level in cohort.disease_status["levels"]},
        "xenium_only_donors": [f"{r['donor_id']} ({r[xenium_cell]})" for r in coverage["table"]
                               if r["design"] == XENIUM_ONLY],
        "xenium_batches": [f"batch {r['batch']}: {r['run_name']}, "
                           + " and ".join(s["donor"] for s in r["slides"])
                           for r in sorted(runs, key=lambda r: int(r["batch"] or 0))
                           if r["platform"] == "xenium" and r["batch"]],
        "xenium_pln_blocks_unsettled": sorted(
            f"{a['acquisition_id']} ({a['donor_id']})" for a in acquisitions
            if a["platform"] == "xenium" and a["tissue"] == LYMPH_NODE and a.get("block") == ""),
        "same_slide_xenium_not_listed": coverage["not_listed"],
        "datasets": {**coverage["totals"]["per_column"], "all": coverage["totals"]["datasets"]},
        "slides": coverage["totals"]["slides"],
        "datasets_not_counted": coverage["totals"]["not_counted"],
        "conditions": sorted(f"{a['acquisition_id']} ({a['donor_id']}): {a['known_condition']}"
                             for a in acquisitions if a.get("known_condition")),
        "partners_outside_this_inventory": sorted(
            f"{s['label']}, run with {', '.join(o['label'] for o in r['slides'] if o is not s)}"
            for r in runs for s in r["slides"] if not s["acquisitions"]),
    }


#: Fills: a platform colours an acquisition, a block colours a section -- never a quantity.
FILLS = {PHENOCYCLER: "#cfe0f3", XENIUM: "#cdebe2", POST_XENIUM_IF: "#e4dcf3",
         PANCREAS: "#a9d18e", LYMPH_NODE: "#f4b183"}
DAMAGED = "#b3261e"
#: The coverage figure's layout, in inches.
_HEADERS = {SAME_BLOCK: "same-block\npLN/lymphoid", DIFFERENT_BLOCK: "different-block\npLN/lymphoid",
            POST_XENIUM_IF: "post-Xenium\nIF"}
_SECTION_WIDTH = {PANCREAS_SECTION: 0.66, SAME_BLOCK: 0.92, DIFFERENT_BLOCK: 0.92,
                  PLN_SECTION: 0.8}
_ACQUISITION_WIDTH = {POST_XENIUM_IF: 1.12}
_GAP, _GROUP_GAP, _ROW, _LINE, _LEFT, _STATUS = 0.3, 0.34, 0.4, 0.118, 1.85, 0.55
_THUMB_ROW = 1.05        # a row that shows a thumbnail


def _width(column: str, kind: str, group: str = "") -> float:
    if kind == "slide":     # a PhenoCycler pancreas slide may be a Xenium slide too: two images
        return 1.3 if group == PC_SLIDE else 0.92
    if kind == "section":
        return _SECTION_WIDTH.get(column, 0.56)
    return _ACQUISITION_WIDTH.get(column, 1.36)


def _figure(task: Task, layout: dict) -> None:
    """One row per donor, the donors of each Xenium design apart; each physical slide's columns
    together under its name, a space between slides. A section is filled by its tissue and counted
    in its block's column, an acquisition by its platform and says its panel and run; the last two columns count
    each donor's datasets and slides, and the bottom row every column's datasets."""
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.patches import Patch, Rectangle

    cells, columns, groups = layout["cells"], layout["columns"], layout["groups"]
    status_x = _LEFT + 0.1
    xs, x = [], status_x + _STATUS + _GAP / 2
    for i, (group, column, kind) in enumerate(columns):
        if i and group != columns[i - 1][0]:
            x += _GAP
        xs.append((x, _width(column, kind, group)))
        x += _width(column, kind, group)
    totals = [(x + _GAP, 0.7), (x + _GAP + 0.7, 0.6)]
    right = x + _GAP + 1.3
    ys, y, separators = {}, 1.45, []
    for g, (_, members) in enumerate(groups):
        if g:
            separators.append(y + _GROUP_GAP / 2)
            y += _GROUP_GAP
        for donor in members:
            lines = max([sum(1 + (len(part) > 12 and "_" in part) for part in c["text"].split("\n"))
                         for k, c in cells.items() if k[0] == donor and c["text"]] or [1])
            thumbed = any(d == donor for d, _ in layout.get("thumbs", {}))
            ys[donor] = (y, max(_THUMB_ROW if thumbed else _ROW, _LINE * lines + 0.16))
            y += ys[donor][1]
    totals_y = y + 0.16
    width, height = right + 0.2, totals_y + _ROW + 0.62

    fig = plt.figure(figsize=(width, height))
    ax = fig.add_axes((0, 0, 1, 1))
    ax.set_xlim(0, width)
    ax.set_ylim(height, 0)
    ax.set_autoscale_on(False)          # a thumbnail must never move the axes
    ax.set_axis_off()
    ax.text(width / 2, 0.28, "Datasets per donor and the slides they came from (stage 00)",
            ha="center", va="center", fontsize=10.5)
    spans: dict[str, list[float]] = {}
    for (group, column, _), (x0, w) in zip(columns, xs, strict=True):
        spans.setdefault(group, [x0, x0 + w])[1] = x0 + w
        ax.text(x0 + w / 2, 1.13, _HEADERS.get(column, column), ha="center", va="center",
                fontsize=7.3, linespacing=1.1)
    for group, (x0, x1) in spans.items():
        ax.text((x0 + x1) / 2, 0.68, group, ha="center", va="center", fontsize=8.6,
                fontweight="bold")
        ax.plot([x0 + 0.03, x0 + 0.03, x1 - 0.03, x1 - 0.03], [0.9, 0.82, 0.82, 0.9],
                color="0.35", linewidth=0.9)
    for (x0, w), name in zip(totals, (DATASETS, SLIDES), strict=True):
        ax.text(x0 + w / 2, 1.13, name, ha="center", va="center", fontsize=7.3)
    ax.text(status_x + _STATUS / 2, 1.13, "status", ha="center", va="center", fontsize=7.3)
    for label, members in groups:
        y0, y1 = ys[members[0]][0], sum(ys[members[-1]])
        ax.text(0.2, (y0 + y1) / 2, label, rotation=90, ha="center", va="center", fontsize=8.6)
        ax.plot([0.36, 0.36], [y0 + 0.05, y1 - 0.05], color="0.35", linewidth=0.9)
    for s in separators:
        ax.plot([0.1, right], [s, s], color="0.55", linewidth=0.8)

    for donor, (y0, h) in ys.items():
        excluded = [p for p in PLATFORM_NAMES if excluded_on(donor, p)]
        ink = "0.5" if len(excluded) == len(PLATFORM_NAMES) else "black"
        ax.text(_LEFT, y0 + h / 2, donor + (
            " (excluded)" if ink != "black" else
            f" ({', '.join(PLATFORM_NAMES[p] for p in excluded)} excluded)" if excluded else ""),
            ha="right", va="center", fontsize=8, color=ink)
        ax.text(status_x + _STATUS / 2, y0 + h / 2,
                layout["table_by_donor"][donor]["disease_status"], ha="center", va="center",
                fontsize=8, color=ink)
        for (group, column, kind), (x0, w) in zip(columns, xs, strict=True):
            cell = cells.get((donor, group, column))
            if cell is None:
                continue
            platform = PLATFORM_OF.get(column) or ("xenium" if group == XENIUM_SLIDE
                                                   else "phenocycler")
            faded = excluded_on(donor, platform)
            cell_ink = "0.5" if faded else "black"
            box = ((x0 + 0.025, y0 + 0.025), w - 0.05, h - 0.05)
            middle = (x0 + w / 2, y0 + h / 2)
            if cell["missing"]:
                ax.add_patch(Rectangle(*box, fill=False, edgecolor="0.5", linewidth=0.8,
                                       linestyle=(0, (3, 2))))
                text = "\n".join(cell["text"].rsplit(" ", 1)) if kind == "section" else cell["text"]
                ax.text(*middle, text, ha="center", va="center", fontsize=6.2, color="0.4",
                        linespacing=1.1)
                continue
            if kind == "slide":
                text = "\n".join(part.replace("_", "\n", 1) if len(part) > 12 else part
                                 for part in cell["text"].split("\n"))
                thumbs = layout.get("thumbs", {}).get((donor, group))
                if not thumbs:
                    ax.text(*middle, text, ha="center", va="center", color=cell_ink, fontsize=6.6,
                            linespacing=1.15)
                    continue
                text_h = _LINE * (text.count("\n") + 1) + 0.04
                gap = 0.06
                room_w = (w - 0.1 - gap * (len(thumbs) - 1)) / len(thumbs)
                room_h = h - 0.08 - text_h
                sizes = [(t.shape[1] * s, t.shape[0] * s) for t in thumbs
                         for s in (min(room_w / t.shape[1], room_h / t.shape[0]),)]
                ix = x0 + (w - sum(iw for iw, _ in sizes) - gap * (len(thumbs) - 1)) / 2
                iy = y0 + 0.04
                for thumb, (iw, ih) in zip(thumbs, sizes, strict=True):
                    ax.imshow(thumb, extent=(ix, ix + iw, iy + ih, iy), origin="upper",
                              aspect="auto", interpolation="antialiased", zorder=2,
                              alpha=0.5 if faded else 1.0)
                    ix += iw + gap
                tallest = max(ih for _, ih in sizes)
                ax.text(x0 + w / 2, iy + tallest + 0.02 + text_h / 2, text, ha="center",
                        va="center", color=cell_ink, fontsize=6.4, linespacing=1.1)
                continue
            ax.add_patch(Rectangle(*box, facecolor=FILLS[cell["kind"]], edgecolor="none",
                                   alpha=0.5 if faded else 1.0))
            if cell["damaged"]:
                ax.add_patch(Rectangle(*box, fill=False, hatch="////", edgecolor=DAMAGED,
                                       linewidth=0))
                ax.add_patch(Rectangle(*box, fill=False, edgecolor=DAMAGED, linewidth=1.3))
            if cell["text"]:
                ax.text(*middle, cell["text"], ha="center", va="center", color=cell_ink,
                        fontsize=6.6 if kind == "acquisition" else 7.6, linespacing=1.15,
                        bbox=(dict(boxstyle="round,pad=0.12", facecolor="white", edgecolor="none")
                              if cell["damaged"] else None))
        for (x0, w), value in zip(totals, (layout["table_by_donor"][donor][DATASETS],
                                           layout["table_by_donor"][donor][SLIDES]), strict=True):
            ax.text(x0 + w / 2, y0 + h / 2, str(value), ha="center", va="center", fontsize=8.4,
                    fontweight="bold", color=ink)

    ax.plot([0.1, right], [totals_y - 0.08] * 2, color="0.25", linewidth=0.9)
    middle = totals_y + _ROW / 2
    left_out = sorted(layout["totals"]["not_counted"])
    ax.text(_LEFT, middle, "total" + (f" (excluded data of {', '.join(left_out)} not counted)"
                                      if left_out else ""),
            ha="right", va="center", fontsize=8, fontweight="bold")
    for (group, column, kind), (x0, w) in zip(columns, xs, strict=True):
        if kind == "acquisition":
            ax.text(x0 + w / 2, middle, str(layout["totals"]["per_column"][f"{group}: {column}"]),
                    ha="center", va="center", fontsize=8.4, fontweight="bold")
    for (x0, w), value in zip(totals, (layout["totals"]["datasets"], layout["totals"]["slides"]),
                              strict=True):
        ax.text(x0 + w / 2, middle, str(value), ha="center", va="center", fontsize=8.4,
                fontweight="bold")

    handles = [Patch(facecolor=FILLS[k], label=k) for k in (PHENOCYCLER, XENIUM, POST_XENIUM_IF)]
    handles += [Patch(facecolor=FILLS[PANCREAS], label=tissues.DISPLAY_NAME[PANCREAS]),
                Patch(facecolor=FILLS[LYMPH_NODE], label=PLN_SECTION)]
    if layout["damage_when"]:
        handles.append(Patch(facecolor="white", edgecolor=DAMAGED, hatch="////",
                             label="damaged " + "; ".join(layout["damage_when"])))
    fig.legend(handles=handles, loc="lower center", ncol=len(handles), fontsize=7.6,
               frameon=False, bbox_to_anchor=(0.5, 0.004))
    fig.savefig(gates.figure_path(task, "coverage.png"), dpi=150)
    plt.close(fig)


def _runs_figure(task: Task, runs: list[dict]) -> None:
    """Every run, oldest first within each platform: its date, panel and the slides run
    together. A PhenoCycler slide alone in its run is set apart (an odd number in its staining
    batch, the maintainer's account)."""
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.patches import Rectangle

    rows = []
    for r in runs:
        rows.append((r["date"], "PhenoCycler" if r["platform"] == "phenocycler" else "Xenium",
                     r["run_name"], r["batch"],
                     "; ".join(sorted({s["panel"] for s in r["slides"]} - {""})),
                     "   |   ".join(_slide_entry(s) for s in r["slides"]),
                     len(r["slides"])))
    palette = plt.get_cmap("tab20")
    colour = {name: palette(i % 20) for i, name in enumerate(sorted({row[4] for row in rows}))}
    n = max(len(rows), 1)
    fig, ax = plt.subplots(figsize=(17.0, 0.25 * (n + 2) + 1.4))
    ax.set_axis_off()
    step = 1.0 / (n + 1)
    columns = (0.0, 0.075, 0.155, 0.335, 0.385, 0.50)
    for x, head in zip(columns, ("run date", "platform", "run", "batch", "panel",
                                 "slides run together"), strict=True):
        ax.text(x, 1 - step / 2, head, fontsize=8, fontweight="bold", va="center",
                transform=ax.transAxes)
    for i, row in enumerate(rows):
        y = 1 - step * (i + 1.5)
        if i % 2:
            ax.add_patch(Rectangle((0, y - step / 2), 1, step, transform=ax.transAxes,
                                   color="0.95", zorder=0))
        ax.add_patch(Rectangle((columns[4] - 0.004, y - 0.4 * step), 0.11, 0.8 * step,
                               transform=ax.transAxes, color=colour[row[4]], alpha=0.35, zorder=1))
        alone = row[6] == 1
        for x, value in zip(columns, row[:6], strict=True):
            ax.text(x, y, value, fontsize=7.2, va="center", transform=ax.transAxes, zorder=2,
                    color="#9a4a00" if alone and x == columns[5] else "black",
                    style="italic" if alone and x == columns[5] else "normal")
    ax.set_title("Runs (stage 00): the slides run together, the panel and the run's date", fontsize=10)
    fig.text(0.01, 0.005,
             "PhenoCycler: slides scanned under one experiment on one day were run together; the"
             " date is the scan day; slide is the Xenium slide's serial where Xenium also ran on"
             " it. A slide alone in its run (italic) means an odd number in its staining batch.\n"
             "Xenium: slides with one run name in experiment.xenium; the date is its run_start_time,"
             " slide the slide's serial, batch the number its run name carries."
             " Panels are named the maintainer's way (multiplex5000.panel_names).", fontsize=7)
    fig.tight_layout(rect=(0, 0.04, 1, 1))
    fig.savefig(gates.figure_path(task, "runs.png"), dpi=150)
    plt.close(fig)


def gate(task: Task, summary: dict) -> gates.Gate:
    C = gates.Check
    checks = [
        C("acquisition ids are unique", summary["repeated_acquisition_ids"], "none", not
          summary["repeated_acquisition_ids"]),
        C("every cohort donor has a line in the donor notes",
          summary["cohort_donors_without_notes"], "none",
          not summary["cohort_donors_without_notes"]),
        C("every image in the cohort spec is in the QuPath manifest",
          summary["cohort_images_not_in_manifest"], "none",
          not summary["cohort_images_not_in_manifest"]),
        C("every raw image is the size its record says", summary["raw_size_differs"], "none",
          not summary["raw_size_differs"]),
        C("every disagreement between the recorded sources is known",
          summary["unknown_disagreements"], "none",
          not summary["unknown_disagreements"],
          detail="known ones are listed, with their reason, in relationships.KNOWN_DISAGREEMENTS"),
        C("every PhenoCycler image's first page is read: its experiment, slide and scan time",
          summary["phenocycler_unread"] + summary["qptiffs_unreadable"], "none",
          not (summary["phenocycler_unread"] or summary["qptiffs_unreadable"]),
          detail=f"{summary['qptiffs_read']} qptiffs read in both PhenoCycler image folders"),
        C("every pancreas image records the markers of its QuPath-manifest panel",
          summary["markers_not_manifest_panel"], "none", not summary["markers_not_manifest_panel"]),
        C("every panel is named by its platform's naming rule", summary["panels_unnamed"], "none",
          not summary["panels_unnamed"], detail="the rules are each platform profile's panel_names"),
        C("every listed Xenium output has its three instrument records",
          summary["xenium_records_missing"] + summary["xenium_records_unreadable"], "none",
          not (summary["xenium_records_missing"] or summary["xenium_records_unreadable"]),
          detail="experiment.xenium, gene_panel.json and metrics_summary.csv, from the output or "
                 "from the copy of its records taken from the complete output"),
        C("each Xenium output's records agree with each other and with its slide",
          summary["xenium_records_disagree"], "none", not summary["xenium_records_disagree"]),
        C("no run holds more than two slides", summary["runs_over_two_slides"], "none",
          not summary["runs_over_two_slides"]),
        C("every recorded tissue condition names an acquisition and a known tissue",
          summary["conditions_unmatched"], "none", not summary["conditions_unmatched"],
          detail="tissue_condition.csv; tissues in any spelling the tissue profiles accept, "
                 "';'-separated"),
        C("every donor has a recorded disease status", summary["donors_without_status"], "none",
          not summary["donors_without_status"],
          detail="the project's donor metadata, else nPOD's case export, in the spellings the "
                 "cohort profile's disease_status names"),
        C("the donor metadata and nPOD's case export agree on each donor's disease status",
          summary["status_disagreements"], "none", not summary["status_disagreements"]),
        C("every Xenium pLN section's block is settled by the donor notes",
          summary["xenium_pln_blocks_unsettled"], "none",
          not summary["xenium_pln_blocks_unsettled"],
          detail="a sequential section's pLN is same-block where its notes line gives same-block "
                 "pLN and different-block where it gives none; a line that gives both leaves it "
                 "unsettled"),
        C("every Xenium output folder names the region its experiment.xenium records",
          summary["xenium_folders_renamed"], "none", not summary["xenium_folders_renamed"],
          detail="the acquisition id takes the recorded region; the tissue is xenium_paths.csv's. "
                 "Known ones are listed, with their reason, in KNOWN_RENAMED_FOLDERS"),
        C("donors the notes describe that the cohort spec does not list",
          summary["notes_donors_not_in_cohort"], "reported", True, blocking=False),
        C("known disagreements", summary["known_disagreements"], "reported", True,
          blocking=False),
        C("raw images whose modification time is not their record's",
          len(summary["raw_mtime_differs"]), "reported", True, blocking=False,
          detail=", ".join(summary["raw_mtime_differs"][:6])),
        C("raw files not found on this machine", len(summary["raw_not_found"]), "reported", True,
          blocking=False, detail=", ".join(summary["raw_not_found"][:6])),
        C("pLN/lymphoid the donor notes designate, by pancreas", summary["pln_designations"],
          "reported", True, blocking=False,
          detail="same-block or different-block as each line writes it; a pancreas not listed "
                 "has none in its line. Checked against pln_blocks.csv by the disagreement check "
                 "above"),
        C("cautions in the donor notes, quoted", summary["cautions"], "reported", True,
          blocking=False, detail="word for word; nothing in them is counted or designated"),
        C("disease status, by donor", summary["disease_status"], "reported", True,
          blocking=False, detail="each donor's source is in coverage.csv"),
        C("Xenium-only donors", summary["xenium_only_donors"], "reported", True, blocking=False,
          detail="no PhenoCycler acquisition and no line in the donor notes"),
        C("Xenium batches", summary["xenium_batches"], "reported", True, blocking=False,
          detail="the batch number each Xenium run name carries; the PhenoXeno runs carry none, "
                 "and no PhenoCycler record carries one"),
        C("datasets, by the slide and platform they came from", summary["datasets"], "reported",
          True, blocking=False,
          detail=f"{summary['slides']} slides; one Xenium slide is one dataset however many "
                 "regions it has. An exclusion takes out only the platform it covers "
                 "(cohort.EXCLUDED_PLATFORMS); datasets it takes out: "
                 + (", ".join(f"{d} ({n})" for d, n in
                              summary["datasets_not_counted"].items()) or "none")),
        C("Xenium run on a PhenoCycler slide with no output listed",
          summary["same_slide_xenium_not_listed"], "reported", True, blocking=False,
          detail="section_modes.csv declares the slide same_slide; xenium_paths.csv lists no "
                 "output for it"),
        C("acquisitions per platform", summary["acquisitions"], "reported", True,
          blocking=False, detail="each Xenium region output counted, as the instruments list "
                                 "them; the datasets above count each slide once"),
        C("panels, as recorded and named", summary["panels"], "reported", True, blocking=False),
        C("runs per platform", summary["runs"], "reported", True, blocking=False,
          detail="every run, its date, panel and slides: runs.csv and runs.png"),
        C("known damaged tissue", summary["conditions"], "reported", True, blocking=False,
          detail="as the maintainer recorded it in tissue_condition.csv"),
        C("PhenoCycler slides run alone", summary["slides_run_alone"], "reported", True,
          blocking=False, detail="an odd number in the staining batch (the maintainer)"),
        C("run partners that are no acquisition here", summary["partners_outside_this_inventory"],
          "reported", True, blocking=False),
        C("Xenium outputs this machine holds only in part", len(summary["xenium_held_in_part"]),
          "reported", True, blocking=False,
          detail="their records are read from the copy taken from the complete output: "
                 + ", ".join(summary["xenium_held_in_part"][:4])
                 + (" ..." if len(summary["xenium_held_in_part"]) > 4 else "")),
        C("slides with no thumbnail", summary["slides_without_thumbnail"], "reported", True,
          blocking=False, detail="a PhenoCycler slide's thumbnail is its qptiff's own; a Xenium "
                                 "slide's is its overview scan, copied beside its records"),
        C("every PhenoCycler slide label carries the listed donor, or its error is known",
          summary["phenocycler_slide_labels_differ"], "none",
          not summary["phenocycler_slide_labels_differ"],
          detail="the SlideID typed when the slide was scanned (the qptiff's first page); the "
                 "donor stays the one the QuPath manifest or pln_node_annotations.csv lists. "
                 "Known errors are listed, with their reason, in KNOWN_LABEL_ERRORS"),
        C("every Xenium cassette name carries the listed donor, or its error is known",
          summary["xenium_cassettes_differ"], "none", not summary["xenium_cassettes_differ"],
          detail="the cassette name is the label typed at the instrument (experiment.xenium); "
                 "the donor stays the one xenium_paths.csv lists. Known errors are listed, with "
                 "their reason, in KNOWN_LABEL_ERRORS"),
        C("labels that name another donor, with their known reason",
          summary["labels_differ_known"], "reported", True, blocking=False),
        C("Xenium output folders renamed, with their known reason",
          summary["xenium_folders_renamed_known"], "reported", True, blocking=False),
    ]
    return gates.Gate(STAGE, "What the release is built from: every acquisition, every recorded "
                      "relationship and what is missing, before anything is computed.",
                      checks, ["coverage.png", "runs.png"])
