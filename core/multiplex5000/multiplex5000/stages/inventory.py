"""Stage 00 -- the source inventory: every acquisition, every recorded relationship, what is missing.

What the release is built from, before anything is computed from it:

* **acquisitions** -- the PhenoCycler pancreas images (the QuPath manifest, checked against the
  cohort spec), the dedicated pLN slides (``pln_node_annotations.csv``), the Xenium region
  outputs (``xenium_paths.csv``) and the post-Xenium IF images (the restain tables). Each gets an
  ``acquisition_id`` made of its platform and the name that platform recorded for it --
  ``phenocycler:<QuPath image>``, ``xenium:<serial>:<region>``, ``if:<image file>`` -- never a
  name worked out from a filename. Raw files are checked against their records by size and
  modification time; nothing large is read.
* **relationships** -- what ``multiplex5000.relationships`` reads from the donor notes,
  ``pln_blocks.csv`` and ``section_modes.csv``, with every disagreement between them.
* **coverage** -- for each donor, how many sources of each kind there are and, beside its pancreas
  image, the pLN/lymphoid the donor notes designate for that pancreas, same-block or
  different-block, checked against ``pln_blocks.csv``; as a table and the gate figure. The notes'
  cautions are marked on the figure and quoted word for word on the gate; nothing in them is
  counted or designated.

The gate blocks on: a repeated acquisition id, a cohort donor the notes do not describe, a cohort
image the QuPath manifest does not have, a raw image whose size is not its record's, and any
disagreement between the sources that is not one already known and explained.

MultiPlex5000 D3 (2026-10-20).
"""

from __future__ import annotations

import csv
import json
from collections import Counter
from pathlib import Path

from .. import gates, profiles, relationships, roots
from ..integration.tissues import LYMPH_NODE, PANCREAS
from ..cohort import DONOR_EXCLUSIONS
from ..integration.manifest import parse_bundle_path
from ..tasks import Input, Task

STAGE = "00_inventory"

#: The coverage columns: one per kind of source a donor can have.
COLUMNS = ("PhenoCycler pancreas", "PhenoCycler pLN slide", "Xenium pancreas", "Xenium pLN",
           "post-Xenium IF", "donor notes")

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
        path = project.source(name)
        if path is None:
            raise relationships.RelationshipError(
                f"{STAGE} needs the {name} source; name it in the site file or the project")
        found.append(Input(name, path, required=required))
    return found


def settings(project) -> dict:
    return {"xenium_root": project.settings.get("xenium_root")}


def _rows(path: Path) -> list[dict]:
    with Path(path).open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(line for line in handle if not line.startswith("#")))


def _stat(path) -> tuple:
    if path is None:
        return None, None, None
    path = Path(path)
    if not path.exists():
        return False, None, None
    stat = path.stat()
    return True, stat.st_size, stat.st_mtime_ns


def _acquisitions(task: Task, spec: dict, manifest: dict, aliases: dict) -> list[dict]:
    out = []
    spec_images = {i["image_id"] for i in spec.get("images", [])}
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
            in_cohort_spec=image["image_id"] in spec_images))
    for index, row in enumerate(_rows(task.input("pln_node_annotations")), start=1):
        out.append(dict(
            acquisition_id=f"phenocycler:{row['image_file']}", platform="phenocycler",
            tissue=LYMPH_NODE, roi=row["roi"],
            donor_id=aliases.get(row["donor_id"], row["donor_id"]), donor_label=row["donor_id"],
            source="pln_node_annotations", source_row=index, path=None, exists=None,
            note=f"IO60pLN slide; excluded={row.get('excluded', '')}"))
    root = task.settings.get("xenium_root")
    for index, row in enumerate(_rows(task.input("xenium_paths")), start=1):
        serial, region = parse_bundle_path(row["source_path"])
        recorded = row["source_path"].strip()
        if recorded and recorded.upper() != "NA" and not Path(recorded).is_absolute() and root:
            path = str(Path(root) / recorded)
        else:
            path = recorded if recorded and recorded.upper() != "NA" else None
        exists, size, _ = _stat(path)
        tissue = LYMPH_NODE if row["roi"].startswith("pln") else PANCREAS
        out.append(dict(
            acquisition_id=(f"xenium:{serial}:{region}" if serial
                            else f"xenium:unrecorded:{row['donor_id']}:{row['roi']}:{index}"),
            platform="xenium", tissue=tissue, roi=row["roi"],
            donor_id=aliases.get(row["donor_id"], row["donor_id"]), donor_label=row["donor_id"],
            source="xenium_paths", source_row=index, path=path, exists=exists,
            note=f"batch {row.get('batch', '')}" + ("" if serial else "; no bundle recorded")))
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


def _donor_order(donor_id: str) -> tuple:
    return (not donor_id.isdigit(), donor_id.zfill(8))


def _coverage(acquisitions, reading, spec) -> list[dict]:
    """One row per donor: its acquisitions by kind, its donor-notes line and, beside the pancreas
    image, the notes' pLN/lymphoid designations -- blank where the line does not say -- with the
    pieces ``pln_blocks.csv`` annotates for the same block, so a difference shows."""
    counts: dict[str, Counter] = {}
    kind = {("phenocycler", PANCREAS): COLUMNS[0], ("phenocycler", LYMPH_NODE): COLUMNS[1],
            ("xenium", PANCREAS): COLUMNS[2], ("xenium", LYMPH_NODE): COLUMNS[3],
            ("if", PANCREAS): COLUMNS[4]}
    for a in acquisitions:
        counts.setdefault(a["donor_id"], Counter())[kind[(a["platform"], a["tissue"])]] += 1
    for donor in reading.donors:
        counts.setdefault(donor, Counter())[COLUMNS[5]] += 1
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
    acquisitions = _acquisitions(task, spec, manifest, cohort.donor_aliases)
    coverage = _coverage(acquisitions, reading, spec)
    _write_csv(task.output("acquisitions.csv"), acquisitions)
    _write_csv(task.output("coverage.csv"), coverage)
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
    platforms = Counter(a["platform"] for a in acquisitions)
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
    }


#: The figure's colour for each column: acquisitions blue, the notes' designations by block.
_MAPS = {SAME_BLOCK: "Greens", DIFFERENT_BLOCK: "Oranges"}


def _cell(row: dict, column: str) -> tuple:
    """(value for the colour, the text, whether the notes and pln_blocks.csv differ)."""
    if column not in _MAPS:
        value = row[column]
        return value, (f"{value}" if value else "–"), False
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
    fig, ax = plt.subplots(figsize=(10.5, 0.32 * max(len(donors), 4) + 2.4))
    ax.imshow(image, aspect="auto")
    for i, row in enumerate(cells):
        for j, (value, text, differs) in enumerate(row):
            scale = top[SHOWN[j] in _MAPS]
            ink = "0.55" if not value else ("white" if value >= 0.6 * scale else "black")
            ax.text(j, i, text, ha="center", va="center", fontsize=8, color=ink)
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
        "Blue: acquisitions of each kind; the last column is the donor-notes line.\n"
        "Green / orange: the pLN/lymphoid the donor notes designate for each pancreas, same-block /"
        " different-block, checked against pln_blocks.csv;\n"
        "'n / m' (red edge) = notes / pln_blocks.csv where they differ; ? = the notes do not say."
        " † = the notes carry a caution, quoted on the gate. * = not in the cohort spec.",
        fontsize=7.5)
    fig.tight_layout()
    fig.savefig(gates.figure_path(task, "coverage.png"), dpi=150)
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
    ]
    return gates.Gate(STAGE, "What the release is built from: every acquisition, every recorded "
                      "relationship and what is missing, before anything is computed.",
                      checks, ["coverage.png"])
