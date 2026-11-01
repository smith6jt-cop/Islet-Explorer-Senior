"""The Source Inventory's MALDI part: each MALDI tissue, its inspection, and its place in the
coverage (2026-10-31).

Read only for a project that asks for MALDI (``settings.maldi: true``), from the ``maldi_spots``
source (a folder holding the export files and the ``SHA256SUMS`` that binds them,
``multiplex5000.maldi``); for every other project the Source Inventory is exactly what it was.

* **acquisitions** -- one per Tissue_id whose donor the IO60 sources know (the donor notes or
  any other acquisition), ``maldi:<Tissue_id>``; the two feature files are one acquisition, one
  dataset. A Tissue_id whose donor they do not know (6426, the HPAP sections) is counted as out
  of cohort and listed on the gate, not as a row.
* **inspection** -- which tissue each MALDI section is (pancreas or pLN/lymphoid) and which
  acquisition's slide it was cut next to are read from nothing but the maintainer's look at the
  ion images beside the donor's slides (``maldi_sections.png``), recorded in the tracked
  ``maldi_tissue_inspection.csv``. The MALDI sections are sequential slides from the same series
  as the PhenoCycler and Xenium slides (maintainer, 2026-10-31): each inspected row is a
  ``serial_section`` relationship -- different cells, never the same cells.
* **coverage** -- a MALDI section group beside the PhenoCycler, Xenium and dedicated-pLN slides:
  the tissue ids with their summed-ion images, the inspected tissue, and the acquisition with its
  panel. The physical MALDI slide is not recorded, so MALDI sections add datasets and no slides.

The gate blocks on: two export files that do not hold the same spots, a MALDI tissue of an IO60
donor that has not been inspected, and an inspection row that names no such tissue, an unknown
tissue, an acquisition that is not the same donor's, or nobody.
"""

from __future__ import annotations

import csv
from collections import Counter
from pathlib import Path

import numpy as np

from .. import gates, maldi, panel_names, profiles, relationships, roots
from ..cohort import excluded_on
from ..integration import tissues
from ..integration.tissues import LYMPH_NODE, PANCREAS
from ..tasks import Input, Task

MALDI = "MALDI"
MALDI_SLIDE = "MALDI section"
INSPECTION = "maldi_tissue_inspection"
INSPECTION_FILE = "maldi_tissue_inspection.csv"
INSPECTION_COLUMNS = ("tissue_id", "tissue", "neighbour_acquisition", "inspected_by", "note")
FILL = "#f7e3b5"
#: The thumbnail of a MALDI section in the coverage: its summed-ion image, this many pixels wide.
THUMB = 120


def source(project) -> Path | None:
    """The MALDI export, only for a project that asks for MALDI (``settings.maldi: true``): a site
    naming the source must not change any other project's Source Inventory."""
    if not project.settings.get("maldi"):
        return None
    folder = project.source("maldi_spots")
    if folder is None:
        raise relationships.RelationshipError(
            "settings.maldi needs the maldi_spots source; name it in the site file or the project")
    return folder


def inputs(project) -> list[Input]:
    folder = source(project)
    if folder is None:
        return []
    found = [Input("maldi_checksums", Path(folder) / maldi.CHECKSUMS)]
    for f in maldi.source_files(folder):
        found.append(Input(f"maldi_{f.role}", f.path, hash=False, recorded_sha256=f.sha256))
    found.append(Input(INSPECTION, project.source(INSPECTION) or roots.INTEGRATION / INSPECTION_FILE))
    return found


def enabled(task: Task) -> bool:
    return "maldi_checksums" in task.inputs


def _files(task: Task) -> list[maldi.SourceFile]:
    recorded = {f.role: f for f in maldi.source_files(task.input("maldi_checksums").parent)}
    return [maldi.SourceFile(role, task.input(f"maldi_{role}"), recorded[role].sha256)
            for role, _ in maldi.ROLES]


def _inspection_rows(path: Path) -> list[dict]:
    with Path(path).open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(line for line in handle if not line.startswith("#")))


def collect(task: Task, aliases: dict, known_donors: set[str], acquisitions: list[dict]) -> dict:
    """Everything the Source Inventory needs from the MALDI export: its acquisitions (IO60
    donors), the out-of-cohort tissues, the pairing of the two files, the inspection's reading and
    the summed-ion images."""
    files = _files(task)
    con = maldi.connect()
    feats = maldi.features(con, files)
    counts = Counter(role for role, _ in feats)
    rule = profiles.load_all()[("platform", "maldi")].panel_names
    panel, panel_error = None, None
    try:
        panel = panel_names.maldi(counts, rule)
    except panel_names.PanelNameError as err:
        panel_error = str(err)
    pairing = maldi.pairing_problems(con, files)
    rows, out_of_cohort = [], []
    for index, t in enumerate(maldi.tissue_summary(con, files), start=1):
        label, donor, rest = maldi.tissue_donor(t["tissue_id"], aliases)
        if donor not in known_donors:
            out_of_cohort.append(dict(tissue_id=t["tissue_id"], donor_label=label,
                                      n_spots=t["n_spots"]))
            continue
        rows.append(dict(
            acquisition_id=f"maldi:{t['tissue_id']}", platform="maldi", tissue_id=t["tissue_id"],
            donor_id=donor, donor_label=label, id_block_token=rest, source="maldi_spots",
            source_row=index, n_spots=t["n_spots"], rasters=t["rasters"],
            x_min_um=t["x_min_um"], x_max_um=t["x_max_um"], y_min_um=t["y_min_um"],
            y_max_um=t["y_max_um"], panel=panel or "", panel_error=panel_error or "",
            source_files="; ".join(f"{f.path.name} sha256 {f.sha256}" for f in files)))
    images = maldi.summed_ions(con, files, feats, [r["tissue_id"] for r in rows])
    glycogen = [(r, n) for r, n in feats if r == "glycogen"]
    glycogen_images = maldi.summed_ions(con, files, glycogen, [r["tissue_id"] for r in rows])
    glycan = [(r, n) for r, n in feats if r == "glycan"]
    glycan_images = maldi.summed_ions(con, files, glycan, [r["tissue_id"] for r in rows])
    for r in rows:
        _x, _y, summed = images[r["tissue_id"]]
        r["spots_no_ion"] = int((summed <= 0).sum())
        r["spots_on_tissue"] = int(maldi.on_tissue(summed).sum())
    reading = _read_inspection(task, rows, acquisitions)
    for r in rows:
        seen = reading["by_tissue"].get(r["tissue_id"])
        r["tissue"] = seen["tissue"] if seen else ""
        r["neighbour_acquisition"] = seen["neighbour_acquisition"] if seen else ""
        r["inspected_by"] = seen["inspected_by"] if seen else ""
    return dict(rows=rows, out_of_cohort=out_of_cohort, pairing=pairing, panel=panel,
                panel_error=panel_error, features=feats, images=images,
                glycogen_images=glycogen_images, glycan_images=glycan_images, **reading)


def _read_inspection(task: Task, rows: list[dict], acquisitions: list[dict]) -> dict:
    """The maintainer's inspection rows: those that name a MALDI tissue of this inventory, a
    tissue the tissue profiles know, an acquisition of the same donor, and who looked."""
    path = task.input(INSPECTION)
    tissues_here = {r["tissue_id"]: r for r in rows}
    by_id = {a["acquisition_id"]: a for a in acquisitions}
    by_tissue, unusable = {}, []
    if Path(path).is_file():
        for n, row in enumerate(_inspection_rows(path), start=2):
            tid = (row.get("tissue_id") or "").strip()
            filled = {k: (row.get(k) or "").strip() for k in INSPECTION_COLUMNS}
            if not (filled["tissue"] or filled["neighbour_acquisition"] or filled["inspected_by"]):
                continue                     # a template row the maintainer has not filled
            problems = []
            if tid not in tissues_here:
                problems.append("no IO60 MALDI tissue of that id")
            try:
                tissue = tissues.normalize(filled["tissue"])
            except tissues.TissueError:
                tissue = None
                problems.append(f"unknown tissue {filled['tissue']!r}")
            neighbour = by_id.get(filled["neighbour_acquisition"])
            if neighbour is None:
                problems.append(f"no acquisition {filled['neighbour_acquisition']!r}")
            elif tid in tissues_here and neighbour["donor_id"] != tissues_here[tid]["donor_id"]:
                problems.append(f"{filled['neighbour_acquisition']} is donor "
                                f"{neighbour['donor_id']}, not {tissues_here[tid]['donor_id']}")
            if not filled["inspected_by"]:
                problems.append("nobody named as inspecting it")
            if tid in by_tissue:
                problems.append("a second row for this tissue")
            if problems:
                unusable.append(f"row {n} ({tid or 'no tissue_id'}): {'; '.join(problems)}")
                continue
            by_tissue[tid] = dict(tissue=tissue, neighbour_acquisition=filled["neighbour_acquisition"],
                                  inspected_by=filled["inspected_by"], note=filled["note"], row=n)
    not_inspected = sorted(t for t in tissues_here if t not in by_tissue)
    return dict(by_tissue=by_tissue, inspection_unusable=unusable, not_inspected=not_inspected,
                inspection_path=str(path))


def relationships_for(task: Task, collected: dict) -> list[relationships.Relationship]:
    """Each inspected MALDI section as a serial section of the acquisition it neighbours."""
    path = Path(collected["inspection_path"])
    if not path.is_file():
        return []
    from ..artifacts import sha256_file

    digest = sha256_file(path)
    out = []
    for r in collected["rows"]:
        seen = collected["by_tissue"].get(r["tissue_id"])
        if seen is None:
            continue
        evidence = relationships.Evidence(INSPECTION, str(path), digest, f"row {seen['row']}")
        out.append(relationships.Relationship(
            relationships.SERIAL_SECTION, r["donor_id"],
            f"MALDI section {r['tissue_id']} ({tissues.DISPLAY_NAME.get(seen['tissue'], seen['tissue'])})",
            seen["neighbour_acquisition"], (evidence,),
            "sequential slides from the same series (maintainer, 2026-10-31); different cells"))
    return out


def _thumb(image: np.ndarray, top: float) -> np.ndarray:
    """A summed-ion image as a small RGB picture on the common scale ``top`` (``maldi.display``),
    no-spot pixels white."""
    import matplotlib

    matplotlib.use("Agg")
    from matplotlib import colormaps
    from PIL import Image

    shown = maldi.display(image, top)
    rgb = colormaps["magma"](np.nan_to_num(shown, nan=0))[..., :3]
    rgb[~np.isfinite(shown)] = 1.0
    picture = Image.fromarray((rgb * 255).astype(np.uint8))
    picture.thumbnail((THUMB, THUMB))
    return np.asarray(picture)


def _image(collected: dict, tissue_id: str, which: str = "images") -> np.ndarray:
    x, y, summed = collected[which][tissue_id]
    grid = maldi.spot_grid(x, y)
    return maldi.spot_image(x, y, summed, grid)


def extend_coverage(coverage: dict, collected: dict) -> None:
    """Add the MALDI section group to the coverage: cells, columns, the per-donor table, the
    totals and the thumbnails. MALDI sections add datasets, never slides: the physical MALDI
    slide is not recorded."""
    from .inventory import (DATASETS, NOT_INSPECTED, PANCREAS_SECTION, PLN_SECTION, SLIDE_NUMBER)

    cells, table = coverage["cells"], coverage["table_by_donor"]
    top = maldi.display_top([_image(collected, r["tissue_id"]) for r in collected["rows"]])
    by_donor: dict[str, list[dict]] = {}
    for r in collected["rows"]:
        by_donor.setdefault(r["donor_id"], []).append(r)
    for donor, mine in by_donor.items():
        mine = sorted(mine, key=lambda r: r["tissue_id"])
        cells[(donor, MALDI_SLIDE, SLIDE_NUMBER)] = dict(
            text="\n".join(r["tissue_id"] for r in mine), kind=None, damaged=False, missing=False)
        coverage.setdefault("thumbs", {})[(donor, MALDI_SLIDE)] = [
            _thumb(_image(collected, r["tissue_id"]), top) for r in mine]
        seen = Counter(r["tissue"] for r in mine if r["tissue"])
        for column, tissue in ((PANCREAS_SECTION, PANCREAS), (PLN_SECTION, LYMPH_NODE)):
            if seen.get(tissue):
                cells[(donor, MALDI_SLIDE, column)] = dict(text=str(seen[tissue]), kind=tissue,
                                                           damaged=False, missing=False)
        waiting = sum(not r["tissue"] for r in mine)
        if waiting:     # in the first tissue column with nothing inspected in it
            free = next((c for c in (PANCREAS_SECTION, PLN_SECTION)
                         if (donor, MALDI_SLIDE, c) not in cells), None)
            if free is not None:
                cells[(donor, MALDI_SLIDE, free)] = dict(
                    text=f"{waiting} {NOT_INSPECTED}" if waiting > 1 else NOT_INSPECTED,
                    kind=None, damaged=False, missing=True)
        cells[(donor, MALDI_SLIDE, MALDI)] = dict(
            text=(mine[0]["panel"].replace(" + ", "\n+ ")
                  + f"\n{len(mine)} section{'s' * (len(mine) > 1)}"), kind=MALDI,
            damaged=False, missing=False)
        if donor in table:
            table[donor][DATASETS] += len(mine)
    columns = [(MALDI_SLIDE, SLIDE_NUMBER, "slide"), (MALDI_SLIDE, PANCREAS_SECTION, "section"),
               (MALDI_SLIDE, PLN_SECTION, "section"), (MALDI_SLIDE, MALDI, "acquisition")]
    coverage["columns"] += columns
    for row in coverage["table"]:
        mine = by_donor.get(row["donor_id"], [])
        row[f"{MALDI_SLIDE}: tissue ids"] = "; ".join(r["tissue_id"] for r in mine)
        row[f"{MALDI_SLIDE}: {MALDI}"] = len(mine)
    counted = [r for r in collected["rows"] if not excluded_on(r["donor_id"], "maldi")]
    coverage["totals"]["per_column"][f"{MALDI_SLIDE}: {MALDI}"] = len(counted)
    coverage["totals"]["datasets"] += len(counted)


def figure(task: Task, collected: dict, acquisitions: list[dict]) -> None:
    """Each IO60 MALDI tissue beside its donor's slides, for the inspection
    (maldi_tissue_inspection.csv): per donor, each MALDI section's summed-ion and glycogen images,
    then the donor's PhenoCycler pancreas scan with its Tissue and pLN annotations, its dedicated
    pLN slide and its Xenium overview, each titled with the acquisition id an inspection row
    names. Every image is drawn as scanned."""
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    from . import inventory as inv

    by_donor: dict[str, list[dict]] = {}
    for r in collected["rows"]:
        by_donor.setdefault(r["donor_id"], []).append(r)
    neighbours: dict[str, list[dict]] = {}
    for a in acquisitions:
        if a["platform"] in ("phenocycler", "xenium") and a["donor_id"] in by_donor:
            neighbours.setdefault(a["donor_id"], []).append(a)
    donors = sorted(by_donor, key=inv._donor_order)
    panels = (("glycan_images", "glycans"), ("glycogen_images", "glycogen"))
    tops = {which: maldi.display_top([_image(collected, r["tissue_id"], which)
                                      for r in collected["rows"]]) for which, _ in panels}
    width = max([2 * len(by_donor[d]) + len(neighbours.get(d, [])) for d in donors] or [1])
    fig, axes = plt.subplots(max(len(donors), 1), width, figsize=(2.2 * width, 2.6 * max(len(donors), 1)),
                             squeeze=False)
    for ax in axes.ravel():
        ax.set_axis_off()
    for row, donor in enumerate(donors):
        k = 0
        for r in sorted(by_donor[donor], key=lambda r: r["tissue_id"]):
            for which, title in panels:
                image = _image(collected, r["tissue_id"], which)
                axes[row, k].imshow(maldi.display(image, tops[which]), cmap="magma", vmin=0,
                                    vmax=1, interpolation="nearest")
                axes[row, k].set_title(f"{r['tissue_id']}\nMALDI {title}, summed", fontsize=7)
                k += 1
        for a in sorted(neighbours.get(donor, []), key=lambda a: (a["platform"], a["acquisition_id"])):
            ax = axes[row, k]
            k += 1
            picture = None
            if a["platform"] == "phenocycler" and a.get("qptiff"):
                if a["tissue"] == PANCREAS:
                    try:
                        dapi, full = inv._phenocycler_image(task.input(f"qptiff:{a['qptiff']}"))
                    except (KeyError, OSError, ValueError):
                        dapi = None
                    if dapi is not None:
                        ax.imshow(dapi, cmap="gray")
                        name = f"annotations:{a['acquisition_id'].split(':', 1)[1]}"
                        if name in task.inputs and task.input(name).is_file():
                            sy, sx = dapi.shape[0] / full[0], dapi.shape[1] / full[1]
                            for tissue, ring in inv._annotation_rings(task.input(name), PANCREAS):
                                ax.plot(ring[:, 0] * sx, ring[:, 1] * sy, linewidth=0.6,
                                        color=inv.FILLS.get(tissue, "red"))
                        picture = True
                else:
                    picture = inv._thumbnail(("qptiff", a["qptiff"]), task)
            elif a["platform"] == "xenium":
                serial = a.get("recorded_slide") or a["acquisition_id"].split(":")[1]
                picture = inv._thumbnail(("overview", serial), task)
            if isinstance(picture, np.ndarray):
                ax.imshow(picture)
            # the exact id an inspection row names, wrapped where QuPath's image name is long
            ax.set_title(a["acquisition_id"].replace(" - ", "\n - ", 1), fontsize=6)
        axes[row, 0].text(-0.15, 0.5, donor, transform=axes[row, 0].transAxes, rotation=90,
                          ha="right", va="center", fontsize=9)
    fig.suptitle("Each IO60 MALDI section beside its donor's slides: record its tissue and the "
                 "acquisition whose slide it neighbours in maldi_tissue_inspection.csv\n"
                 "MALDI images on one common linear scale per feature set, as measured (0 = not "
                 "detected), with no per-section compensation", fontsize=10,
                 y=1 - 0.15 / fig.get_size_inches()[1], va="top")
    height = fig.get_size_inches()[1]
    fig.tight_layout(rect=(0, 0, 1, 1 - 1.5 / height))
    for i, (which, title) in enumerate(panels):
        maldi.colorbar(fig, (0.12 + 0.42 * i, 1 - 1.0 / height, 0.32, 0.12 / height), tops[which],
                       f"MALDI {title}, summed per spot")
    fig.savefig(gates.figure_path(task, "maldi_sections.png"), dpi=110)
    plt.close(fig)


def summary(collected: dict) -> dict:
    rows = collected["rows"]
    return {
        "maldi_acquisitions": len(rows),
        "maldi_donors": sorted({r["donor_id"] for r in rows}),
        "maldi_out_of_cohort": sorted(f"{o['tissue_id']} ({o['n_spots']} spots)"
                                      for o in collected["out_of_cohort"]),
        "maldi_pairing": collected["pairing"],
        "maldi_pairing_problems": {k: v for k, v in collected["pairing"].items()
                                   if not k.endswith("_rows") and v},
        "maldi_panel": collected["panel"],
        "maldi_panel_error": collected["panel_error"],
        "maldi_not_inspected": collected["not_inspected"],
        "maldi_inspection_unusable": collected["inspection_unusable"],
        "maldi_inspected": sorted(f"{t}: {v['tissue']}, next to {v['neighbour_acquisition']} "
                                  f"({v['inspected_by']})" for t, v in collected["by_tissue"].items()),
        "maldi_spots_no_ion": {r["tissue_id"]: r["spots_no_ion"] for r in rows if r["spots_no_ion"]},
    }


def checks(summary: dict) -> list[gates.Check]:
    C = gates.Check
    return [
        C("the two MALDI export files hold the same spots", summary["maldi_pairing_problems"],
          "none", not summary["maldi_pairing_problems"],
          detail="spot ids, positions, rasters and Tissue_ids, file against file: "
                 + ", ".join(f"{k} {v}" for k, v in summary["maldi_pairing"].items())),
        C("the MALDI panel is named by its platform's naming rule", summary["maldi_panel_error"]
          or summary["maldi_panel"], "a name", bool(summary["maldi_panel"]),
          detail="the platform profile's panel_names (profiles/platforms/maldi.yaml)"),
        C("every IO60 MALDI section has been inspected: its tissue and the slide it neighbours",
          summary["maldi_not_inspected"], "none", not summary["maldi_not_inspected"],
          detail="by eye, on maldi_sections.png, recorded in maldi_tissue_inspection.csv: the only "
                 "authority for a MALDI section's tissue and neighbour; the Tissue_id is never "
                 "parsed for either (maintainer, 2026-10-31)"),
        C("every MALDI inspection row names an IO60 MALDI section, a known tissue, an acquisition "
          "of the same donor and who inspected it", summary["maldi_inspection_unusable"], "none",
          not summary["maldi_inspection_unusable"]),
        C("MALDI sections inspected", summary["maldi_inspected"], "reported", True, blocking=False,
          detail="each a serial section of the acquisition it neighbours: different cells"),
        C("MALDI tissues outside the cohort, not ingested", summary["maldi_out_of_cohort"],
          "reported", True, blocking=False,
          detail="a Tissue_id whose donor no IO60 source knows (maintainer, 2026-10-31: IO60 "
                 "donors only)"),
        C("MALDI sections, by donor", f"{summary['maldi_acquisitions']} sections of "
          f"{len(summary['maldi_donors'])} donors", "reported", True, blocking=False,
          detail="the physical MALDI slide is not recorded: MALDI sections add datasets, not slides"),
        C("MALDI spots with no ion detected", summary["maldi_spots_no_ion"], "reported", True,
          blocking=False, detail=f"off the tissue by {maldi.ON_TISSUE_METHOD}; never a cutoff on "
                                 "a feature's value"),
    ]


def write(task: Task, collected: dict) -> None:
    from .inventory import _write_csv

    _write_csv(task.output("maldi_acquisitions.csv"), collected["rows"])
    _write_csv(task.output("maldi_features.csv"),
               [dict(role=r, feature=n) for r, n in collected["features"]])
    _write_csv(task.output("maldi_out_of_cohort.csv"), collected["out_of_cohort"])
