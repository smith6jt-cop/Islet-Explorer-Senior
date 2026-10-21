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
* **relationships** -- what ``multiplex5000.relationships`` reads from the donor notes,
  ``pln_blocks.csv`` and ``section_modes.csv``, with every disagreement between them.
* **coverage** -- for each donor, how many sources of each kind there are and, beside its pancreas
  image, the pLN/lymphoid the donor notes designate for that pancreas, same-block or
  different-block, checked against ``pln_blocks.csv``; as a table and the gate figure. The notes'
  cautions are marked on the figure and quoted word for word on the gate; nothing in them is
  counted or designated.

The gate blocks on: a repeated acquisition id, a cohort donor the notes do not describe, a cohort
image the QuPath manifest does not have, a raw image whose size is not its record's, any
disagreement between the sources that is not one already known and explained, an instrument
record that is missing, unreadable or at odds with itself, a pancreas image whose markers are not
its manifest panel's, a panel the naming rule does not cover, a run of more than two slides, and a
Xenium output folder that names a region its record does not, until the reason is known.

MultiPlex5000 D3 (2026-10-20).
"""

from __future__ import annotations

import csv
import json
from collections import Counter
from pathlib import Path

from .. import acquisition_records, gates, panel_names, profiles, relationships, roots
from ..integration.tissues import LYMPH_NODE, PANCREAS
from ..cohort import DONOR_EXCLUSIONS
from ..integration.manifest import parse_bundle_path
from ..tasks import Input, Task

STAGE = "00_inventory"

#: The coverage columns: one per kind of acquisition a donor can have.
COLUMNS = ("PhenoCycler pancreas", "PhenoCycler pLN slide", "Xenium pancreas", "Xenium pLN",
           "post-Xenium IF")
PANELS = " panel"      # suffix of the coverage column that names those acquisitions' panels
_KINDS = {("phenocycler", PANCREAS): COLUMNS[0], ("phenocycler", LYMPH_NODE): COLUMNS[1],
          ("xenium", PANCREAS): COLUMNS[2], ("xenium", LYMPH_NODE): COLUMNS[3],
          ("if", PANCREAS): COLUMNS[4]}

#: Where stage 00 reads what the instruments recorded: the qptiff folders of both PhenoCycler
#: projects (every slide there can be another's run partner), and the copies of the records of
#: Xenium outputs this machine holds only in part.
QPTIFF_FOLDERS = ("pancreas_qptiffs", "pln_qptiffs")
XENIUM_RECORDS = ("experiment.xenium", "gene_panel.json", "metrics_summary.csv")

#: Xenium output folders that name a region other than the one their experiment.xenium records,
#: with the reason, once someone has looked. The gate blocks on any other.
KNOWN_RENAMED_FOLDERS: dict[str, str] = {}

#: The pLN/lymphoid the donor notes designate for each pancreas, by block: (coverage column,
#: the block as relationships names it, the DonorNote field that holds the notes' count).
SAME_BLOCK = "same-block pLN/lymphoid"
DIFFERENT_BLOCK = "different-block pLN/lymphoid"
DESIGNATIONS = ((SAME_BLOCK, relationships.SAME_BLOCK, "same_block_pln"),
                (DIFFERENT_BLOCK, relationships.DIFFERENT_BLOCK, "different_block_pln"))
ANNOTATED = " in pln_blocks.csv"     # suffix of the column that holds pln_blocks.csv's pieces

#: The figure's columns, left to right: the pancreas image, what the notes designate with it,
#: then the other acquisitions.
SHOWN = (COLUMNS[0], SAME_BLOCK, DIFFERENT_BLOCK, *COLUMNS[1:])

_TRACKED = {  # name -> file name under core/phenocycler/data/integration
    "xenium_paths": "xenium_paths.csv", "section_modes": "section_modes.csv",
    "restain_sections": "restain_sections.csv",
    "restain_registration_inputs": "restain_registration_inputs.csv",
    "pln_node_annotations": "pln_node_annotations.csv",
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
                           ("pln_blocks", True)):
        found.append(Input(name, _needed(project, name), required=required))
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
    for name, local in _xenium_outputs(project, paths["xenium_paths"]):
        held = local if (local / XENIUM_RECORDS[0]).is_file() else copies / name
        out += [Input(f"xenium:{name}:{f}", held / f, required=False) for f in XENIUM_RECORDS]
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
            copies[row["output"]] = row["source"].rsplit("/", 1)[0]
    return {"headers": headers, "xenium": xenium, "missing": missing, "unreadable": unreadable,
            "copies": copies}


def _acquisitions(task: Task, spec: dict, manifest: dict, aliases: dict,
                  records: dict) -> list[dict]:
    rules = {p: profiles.load_all()[("platform", p)].panel_names for p in ("phenocycler", "xenium")}
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
                      run_date=header.scan_date, software=header.software, qptiff=file)
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
                          recorded_region=record.region, folder_region=folder_region,
                          local_copy=("complete" if Path(path, XENIUM_RECORDS[0]).is_file()
                                      else "partial" if exists else "not on this machine"),
                          records_read_from=("this machine's complete copy"
                                             if Path(path, XENIUM_RECORDS[0]).is_file() else
                                             f"copy of {records['copies'].get(output, '?')}"))
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
        out.append(dict(
            acquisition_id=f"if:{row['image_file']}", platform="if", tissue=PANCREAS,
            roi=row["roi"], donor_id=aliases.get(row["donor_id"], row["donor_id"]),
            donor_label=row["donor_id"], source="restain_registration_inputs", source_row=index,
            path=None, exists=None, xenium_serial=row.get("xenium_serial"),
            note=("channel roles recorded in restain_sections.csv" if roles
                  else "no channel roles recorded")))
    return out


def _slide_label(acquisitions: list[dict]) -> str:
    """'6482 pancreas', or '6539 pancreas + pLN' for a Xenium slide with two tissues."""
    tissues = []
    for a in acquisitions:
        name = "pancreas" if a["tissue"] == PANCREAS else "pLN"
        if name not in tissues:
            tissues.append(name)
    return f"{acquisitions[0]['donor_id']} {' + '.join(tissues)}"


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
            slides.append(dict(slide=h.slide_id, position=h.position, scan_time=h.scan_time,
                               label=(_slide_label([acq]) if acq else
                                      f"{h.slide_id} ({name}, no acquisition here)"),
                               acquisitions=[acq["acquisition_id"]] if acq else [],
                               panel=acq.get("panel", "") if acq else ""))
        out.append(dict(platform="phenocycler", run=f"phenocycler:{experiment}:{day}",
                        run_name=experiment, date=day, slides=slides))
    xenium: dict[str, list[dict]] = {}
    for a in acquisitions:
        if a["platform"] == "xenium" and a.get("run"):
            xenium.setdefault(a["run"], []).append(a)
    for run, members in xenium.items():
        by_slide: dict[str, list[dict]] = {}
        for a in members:
            by_slide.setdefault(a["recorded_slide"], []).append(a)
        slides = [dict(slide=s, position="", scan_time="", label=_slide_label(group),
                       acquisitions=[a["acquisition_id"] for a in group],
                       panel="; ".join(sorted({a.get("panel", "") for a in group})))
                  for s, group in sorted(by_slide.items())]
        out.append(dict(platform="xenium", run=run, run_name=run.split(":", 1)[1],
                        date=min(a["run_date"] for a in members), slides=slides))
    out.sort(key=lambda r: (r["platform"] != "phenocycler", r["date"], r["run"]))
    partners = {}
    for r in out:
        for s in r["slides"]:
            others = [o["label"] + (f" ({o['position']})" if o["position"] else "")
                      for o in r["slides"] if o is not s]
            for acq in s["acquisitions"]:
                partners[acq] = "; ".join(others) or "run alone"
    for a in acquisitions:
        if a["acquisition_id"] in partners:
            a["run_partners"] = partners[a["acquisition_id"]]
    return out


def _donor_order(donor_id: str) -> tuple:
    return (not donor_id.isdigit(), donor_id.zfill(8))


def _coverage(acquisitions, reading, spec) -> list[dict]:
    """One row per donor: its acquisitions by kind with the panels they were acquired with and,
    beside the pancreas image, the notes' pLN/lymphoid designations -- blank where the line does
    not say -- with the pieces ``pln_blocks.csv`` annotates for the same block, so a difference
    shows. A panel no record names is ``?``."""
    counts: dict[str, Counter] = {}
    named: dict[str, dict[str, Counter]] = {}
    for a in acquisitions:
        column = _KINDS[(a["platform"], a["tissue"])]
        counts.setdefault(a["donor_id"], Counter())[column] += 1
        if a["platform"] != "if":
            unrecorded = a["acquisition_id"].startswith("xenium:unrecorded:")
            named.setdefault(a["donor_id"], {}).setdefault(column, Counter())[
                "NA" if unrecorded else (a.get("panel") or "?")] += 1
    for donor in reading.donors:
        counts.setdefault(donor, Counter())
    annotated: Counter = Counter()
    for ann in reading.pln:
        counts.setdefault(ann.donor_id, Counter())
        annotated[(ann.donor_id, ann.block)] += ann.pieces
    cohort = set(map(str, spec.get("expected_donors", [])))
    rows = []
    for d in sorted(counts, key=_donor_order):
        note = reading.donors.get(d)
        row = dict(donor_id=d, in_cohort_spec=d in cohort, centrally_excluded=d in DONOR_EXCLUSIONS)
        row[COLUMNS[0]] = counts[d].get(COLUMNS[0], 0)
        for column, _, field in DESIGNATIONS:
            stated = getattr(note, field) if note is not None else None
            row[column] = "" if stated is None else stated
        row.update({c: counts[d].get(c, 0) for c in COLUMNS[1:]})
        for column in COLUMNS[:4]:
            row[column + PANELS] = "; ".join(
                label if n == 1 else f"{label} ({n})"
                for label, n in sorted(named.get(d, {}).get(column, Counter()).items()))
        for column, block, _ in DESIGNATIONS:
            row[column + ANNOTATED] = annotated[(d, block)]
        row["designation read from free text"] = bool(note is not None and note.pln_from_caution)
        row["caution in the notes"] = bool(note is not None and note.cautions)
        rows.append(row)
    return rows


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
    runs = _runs(acquisitions, records)
    coverage = _coverage(acquisitions, reading, spec)
    _write_csv(task.output("acquisitions.csv"), acquisitions)
    _write_csv(task.output("coverage.csv"), coverage)
    _write_csv(task.output("runs.csv"), [dict(
        platform=r["platform"], run=r["run"], run_name=r["run_name"], date=r["date"],
        slides=len(r["slides"]), panels="; ".join(sorted({s["panel"] for s in r["slides"]} - {""})),
        slides_run_together=" | ".join(
            s["label"] + (f" ({s['position']})" if s["position"] else "") for s in r["slides"]),
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
    ids = Counter(a["acquisition_id"] for a in acquisitions)
    spec_images = [i["image_id"] for i in spec.get("images", [])]
    in_manifest = {i["image_id"] for i in manifest.get("images", [])}
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
                   for p in ("phenocycler", "xenium")},
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
        "xenium_held_in_part": sorted(f"{a['output']} ({a['local_copy']})" for a in xenium
                                      if a.get("local_copy") and a["local_copy"] != "complete"),
        "runs": {p: sum(r["platform"] == p for r in runs) for p in ("phenocycler", "xenium")},
        "runs_over_two_slides": sorted(f"{r['run']}: {len(r['slides'])} slides" for r in runs
                                       if len(r["slides"]) > 2),
        "slides_run_alone": sorted(f"{r['slides'][0]['label']}, {r['date']}" for r in runs
                                   if r["platform"] == "phenocycler" and len(r["slides"]) == 1),
        "partners_outside_this_inventory": sorted(
            f"{s['label']}, run with {', '.join(o['label'] for o in r['slides'] if o is not s)}"
            for r in runs for s in r["slides"] if not s["acquisitions"]),
    }


#: The figure's colour for each column: acquisitions blue, the notes' designations by block.
_MAPS = {SAME_BLOCK: "Greens", DIFFERENT_BLOCK: "Oranges"}


def _cell(row: dict, column: str) -> tuple:
    """(value for the colour, the text, whether the notes and pln_blocks.csv differ)."""
    if column not in _MAPS:
        value = row[column]
        if not value:
            return 0, "–", False
        named = row.get(column + PANELS, "")
        return value, (named.replace("; ", "\n") if named else f"{value}"), False
    stated, annotated = row[column], row[column + ANNOTATED]
    if stated == "":
        return 0, "?", False
    text = f"{stated}" if stated else "–"
    if stated != annotated:
        return stated, f"{stated} / {annotated}", True
    return stated, text, False


def _figure(task: Task, coverage: list[dict]) -> None:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import numpy as np
    from matplotlib.patches import Rectangle

    donors = [r["donor_id"] for r in coverage]
    cells = [[_cell(r, c) for c in SHOWN] for r in coverage]
    # one colour scale for the acquisition counts and one for the designations
    top = {False: 2.0, True: 2.0}
    for row in cells:
        for column, (value, _, _) in zip(SHOWN, row, strict=True):
            top[column in _MAPS] = max(top[column in _MAPS], float(value))
    image = np.ones((len(coverage), len(SHOWN), 4))
    for i, row in enumerate(cells):
        for j, (value, _, _) in enumerate(row):
            if value:
                image[i, j] = plt.get_cmap(_MAPS.get(SHOWN[j], "Blues"))(
                    value / top[SHOWN[j] in _MAPS])
    fig, ax = plt.subplots(figsize=(12.5, 0.38 * max(len(donors), 4) + 2.4))
    ax.imshow(image, aspect="auto")
    for i, row in enumerate(cells):
        for j, (value, text, differs) in enumerate(row):
            scale = top[SHOWN[j] in _MAPS]
            ink = "0.55" if not value else ("white" if value >= 0.6 * scale else "black")
            ax.text(j, i, text, ha="center", va="center", fontsize=7 if "\n" in text else 7.5,
                    color=ink, linespacing=1.0)
            if differs:
                ax.add_patch(Rectangle((j - 0.5, i - 0.5), 1, 1, fill=False,
                                       edgecolor="#b3261e", linewidth=1.8))
    labels = [f"{r['donor_id']}{' (excluded)' if r['centrally_excluded'] else ''}"
              f"{' †' if r['caution in the notes'] else ''}"
              f"{'' if r['in_cohort_spec'] else ' *'}" for r in coverage]
    ax.set_yticks(range(len(donors)), labels, fontsize=8)
    ax.set_xticks(range(len(SHOWN)), [c.replace(" pLN/", "\npLN/") for c in SHOWN],
                  rotation=25, ha="right", fontsize=8)
    ax.set_title("Recorded sources per donor (stage 00)", fontsize=10)
    ax.set_xlabel(
        "Blue: acquisitions of each kind, named by the panel their instrument recorded (count in"
        " brackets where several share one; ? = no record names it; NA = xenium_paths.csv records"
        " no bundle); the shade is the count.\n"
        "Green / orange: the pLN/lymphoid the donor notes designate for each pancreas, same-block /"
        " different-block, checked against pln_blocks.csv;\n"
        "'n / m' (red edge) = notes / pln_blocks.csv where they differ; ? = the notes do not say."
        " † = the notes carry a caution, quoted on the gate. * = not in the cohort spec.",
        fontsize=7.5)
    fig.tight_layout()
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
                     r["run_name"],
                     "; ".join(sorted({s["panel"] for s in r["slides"]} - {""})) or "?",
                     "   |   ".join(s["label"] + (f" ({s['position']})" if s["position"] else "")
                                    for s in r["slides"]),
                     len(r["slides"])))
    palette = plt.get_cmap("tab20")
    colour = {name: palette(i % 20) for i, name in enumerate(sorted({row[3] for row in rows}))}
    n = max(len(rows), 1)
    fig, ax = plt.subplots(figsize=(14.5, 0.25 * (n + 2) + 1.4))
    ax.set_axis_off()
    step = 1.0 / (n + 1)
    columns = (0.0, 0.075, 0.155, 0.36, 0.475)
    for x, head in zip(columns, ("run date", "platform", "run", "panel", "slides run together"),
                       strict=True):
        ax.text(x, 1 - step / 2, head, fontsize=8, fontweight="bold", va="center",
                transform=ax.transAxes)
    for i, row in enumerate(rows):
        y = 1 - step * (i + 1.5)
        if i % 2:
            ax.add_patch(Rectangle((0, y - step / 2), 1, step, transform=ax.transAxes,
                                   color="0.95", zorder=0))
        ax.add_patch(Rectangle((columns[3] - 0.004, y - 0.4 * step), 0.11, 0.8 * step,
                               transform=ax.transAxes, color=colour[row[3]], alpha=0.35, zorder=1))
        alone = row[5] == 1
        for x, value in zip(columns, row[:5], strict=True):
            ax.text(x, y, value, fontsize=7.2, va="center", transform=ax.transAxes, zorder=2,
                    color="#9a4a00" if alone and x == columns[4] else "black",
                    style="italic" if alone and x == columns[4] else "normal")
    ax.set_title("Runs (stage 00): the slides run together, the panel and the run's date", fontsize=10)
    fig.text(0.01, 0.005,
             "PhenoCycler: slides scanned under one experiment on one day were run together; the"
             " date is the scan day; (SlideN) is the slide's position in the experiment. A slide"
             " alone in its run (italic) means an odd number in its staining batch.\n"
             "Xenium: slides with one run name in experiment.xenium; the date is its run_start_time."
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
        C("acquisitions per platform", summary["acquisitions"], "reported", True,
          blocking=False),
        C("panels, as recorded and named", summary["panels"], "reported", True, blocking=False),
        C("runs per platform", summary["runs"], "reported", True, blocking=False,
          detail="every run, its date, panel and slides: runs.csv and runs.png"),
        C("PhenoCycler slides run alone", summary["slides_run_alone"], "reported", True,
          blocking=False, detail="an odd number in the staining batch (the maintainer)"),
        C("run partners that are no acquisition here", summary["partners_outside_this_inventory"],
          "reported", True, blocking=False),
        C("Xenium outputs this machine holds only in part", len(summary["xenium_held_in_part"]),
          "reported", True, blocking=False,
          detail="their records are read from the copy taken from the complete output: "
                 + ", ".join(summary["xenium_held_in_part"][:4])
                 + (" ..." if len(summary["xenium_held_in_part"]) > 4 else "")),
        C("Xenium output folders renamed, with their known reason",
          summary["xenium_folders_renamed_known"], "reported", True, blocking=False),
    ]
    return gates.Gate(STAGE, "What the release is built from: every acquisition, every recorded "
                      "relationship and what is missing, before anything is computed.",
                      checks, ["coverage.png", "runs.png"])
