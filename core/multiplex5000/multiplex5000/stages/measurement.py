"""Stage 03 -- protein measurement: from each image's QuPath export to one value per cell and
marker, in the compartment the calibration review ruled for that marker.

This is the October own-compartment chain, moved (MultiPlex5000 D4). The maintainer, 2026-10-20:
"the October 14 own-compartment and calibration chain is the canonical protein path". Its
measurement half was the package's own pipeline -- ingest, geometry QC, REDSEA, expression --
with two changes, and this stage runs exactly that, in that order:

* **ingest** -- each cohort image's v3 export, a (cells, context) Parquet pair named by the QuPath
  manifest, into one partition per donor (``data/cells``).
* **geometry QC** -- which cells are measurable (``data/geometry_qc``). A cell with no nucleus is
  never admitted (``missing_geometry_policy = exclude``).
* **REDSEA** -- spillover compensation from each image's pixels and cell outlines, every channel in
  four compartments (``data/redsea``).
* **expression** -- one value per cell and marker (``data/expression``) and each marker's
  availability (``data/availability``). The two changes are both made here, as the October chain
  made them (``reprocess_protein_own_compartment_candidate.py``'s ``prepare``):
  - each marker comes from the compartment ``marker_compartment_ruling.csv`` names -- the
    calibration review's ruling of 2026-10-03, copied byte for byte (sha256 126ceba3...): 33
    Membrane, 18 Cytoplasm and 7 Nucleus; every other marker keeps the registry's compartment;
  - the analysis exclusions recorded on 2026-10-05 (the site's ``analysis_exclusions``) are applied
    to the cells the expression admits for analysis and estimation. Beyond geometry QC's own
    reasons, they are the 11,484 reviewed cells (maintainer: "exclude ALL 11,484 reviewed cells
    from ALL analyses"). They are applied after REDSEA, so they never change its spillover context.

Normalization (UniFORM shifts and scale) is stage 05's, not this one's.

The gate shows each marker per donor before and after REDSEA in its measured compartment, and each
donor's nucleus admission. It blocks on:
- a cohort image the manifest does not hold, or one not measured;
- object ids that are not unique or not the same across the stage's tables;
- an admitted cell with no nucleus;
- an exclusion naming a cell the export does not have;
- a ruling naming a marker or a compartment the registry does not know.
"""

from __future__ import annotations

import configparser
import json
import os
import shutil
from collections import Counter
from pathlib import Path
from types import SimpleNamespace

from .. import gates, profiles, roots
from ..cohort import excluded_on
from ..tasks import Input, Task

STAGE = "03_protein_measurement"

#: The package's own knowledge this stage reads, versioned with its code.
PACKAGE = Path(__file__).resolve().parents[1]
REGISTRY = PACKAGE / "marker_registry.json"
RULING = PACKAGE / "marker_compartment_ruling.csv"
#: The compartments a ruling may set; any other value ("not applicable") leaves the registry's.
RULED = ("Membrane", "Cytoplasm", "Nucleus")
#: The pipeline's own settings file, the only config.ini (CLAUDE.md).
CONFIG_INI = roots.REPO / "core" / "multiplex5000" / "config.ini"
#: The site source naming the recorded analysis exclusions: a folder with manifest.json and one
#: Parquet per donor, keyed by object_id.
EXCLUSIONS = "analysis_exclusions"
#: The config.ini sections whose values decide what this stage measures.
SETTINGS_SECTIONS = ("cells", "geometry_qc", "redsea")
#: Where the stage's tables go, under its folder.
DATA = "data"
TABLES = ("cells", "geometry_qc", "redsea", "expression", "availability")


class MeasurementError(RuntimeError):
    """Sources this stage cannot measure as they stand."""


# --------------------------------------------------------------------------------- inputs ----
def _cohort_spec(project) -> Path:
    cohort = profiles.load_all()[("cohort", project.cohort)]
    return project.source("cohort_spec") or profiles.REPO / cohort.membership


def _manifest_path(project) -> Path:
    return project.source("qupath_manifest") or roots.INTEGRATION.parent / "qupath_manifest.json"


def _measured(manifest, spec: dict) -> list:
    """The manifest's images this stage measures: the cohort spec's, less any donor whose
    PhenoCycler data the cohort excludes."""
    wanted = {i["image_id"] for i in spec.get("images", [])}
    return [i for i in manifest.images
            if i.image_id in wanted and not excluded_on(i.donor_id, "phenocycler")]


def _exclusion_files(folder: Path) -> dict[str, tuple[Path, str]]:
    """{donor as the cohort names it: (its exclusion Parquet, recorded sha256)}, from the folder's
    manifest.json -- which records each partition's native donor id, so none is worked out."""
    record = json.loads((folder / "manifest.json").read_text())
    return {str(d["donor_id_native"]): (Path(d["path"]), d.get("sha256"))
            for d in record.get("donors", [])}


def inputs(project) -> list[Input]:
    from ..artifacts import QuPathCohortManifest

    spec_path, manifest_path = _cohort_spec(project), _manifest_path(project)
    found = [Input("cohort_spec", spec_path), Input("qupath_manifest", manifest_path),
             Input("config_ini", CONFIG_INI), Input("marker_registry", REGISTRY),
             Input("compartment_ruling", RULING)]
    manifest = QuPathCohortManifest.read_json(manifest_path)
    images = _measured(manifest, json.loads(Path(spec_path).read_text()))
    for image in images:
        d = image.donor_id
        # The exports and images are large: bound by size and time, with the digest the manifest
        # records for each.
        files = {"cells": image.measurement_csv, "context": image.measurement_context,
                 "cell_geojson": image.cell_geojson, "qptiff": image.qptiff}
        if image.nucleus_geojson is not None and (
                image.nucleus_geojson.path != image.cell_geojson.path):
            files["nucleus_geojson"] = image.nucleus_geojson
        for kind, fp in files.items():
            if fp is not None:
                found.append(Input(f"{kind}:{d}", Path(fp.path), hash=False,
                                   recorded_sha256=fp.sha256))
    folder = project.source(EXCLUSIONS)
    if folder is None:
        raise MeasurementError(f"{STAGE} needs the {EXCLUSIONS} source; name it in the site file "
                               "or the project")
    found.append(Input("exclusions_manifest", Path(folder) / "manifest.json"))
    files = _exclusion_files(Path(folder))
    for image in images:
        if image.donor_id in files:
            path, sha = files[image.donor_id]
            found.append(Input(f"exclusions:{image.donor_id}", path, recorded_sha256=sha))
    return found


def settings(project) -> dict:
    from ..config import _ENV_OVERRIDES

    parser = configparser.ConfigParser(inline_comment_prefixes=("#", ";"))
    parser.read(CONFIG_INI)
    return {
        "config": {s: dict(parser.items(s)) for s in SETTINGS_SECTIONS if parser.has_section(s)},
        "jobs": int(project.settings.get("measurement_jobs", 1)),
        # A setting an environment variable would change silently is recorded, never hidden.
        "environment_overrides": {env: os.environ[env] for env, _ in _ENV_OVERRIDES.values()
                                  if os.environ.get(env)},
    }


# ---------------------------------------------------------------------- the compartments ----
def _ruling(path: Path) -> dict[str, str]:
    import pandas as pd

    frame = pd.read_csv(path, dtype=str)
    return dict(zip(frame["marker"], frame["measurement_compartment"], strict=True))


def own_compartment_registry(payload: dict, ruling: dict[str, str]) -> dict:
    """The registry with each ruled marker's compartment set to the ruling's, as the October chain
    set it; a marker the ruling leaves unset ("not applicable", or not named) keeps its own."""
    from copy import deepcopy

    own = deepcopy(payload)
    for marker in own["markers"]:
        if ruling.get(marker["name"]) in RULED:
            marker["compartment"] = ruling[marker["name"]]
    return own


def ruling_problems(payload: dict, ruling: dict[str, str]) -> list[str]:
    known = {m["name"] for m in payload["markers"]}
    allowed = set(payload.get("allowed_compartments", RULED))
    out = [f"{m}: not in the registry" for m in sorted(set(ruling) - known)]
    out += [f"{m}: {c} is not a compartment the registry allows" for m, c in sorted(ruling.items())
            if c in RULED and c not in allowed]
    return out


# ---------------------------------------------------------------------------------- run ----
def _apply_exclusions(qc, excluded: set):
    """The geometry QC table with the recorded exclusions taken out of analysis and estimation,
    as the October chain took them out; the table REDSEA used is left as it was."""
    out = qc.copy()
    hit = out["object_id"].astype(str).isin(excluded).to_numpy()
    out.loc[hit, "analysis_eligible"] = False
    out.loc[hit, "estimation_eligible"] = False
    return out, hit


def _before(cells_path: Path, image, registry, object_ids):
    """Each marker's value before REDSEA, in the compartment it is measured in: QuPath's own mean
    from the export, by the channel name the image records."""
    import pyarrow.parquet as pq

    names = {c.marker: c.channel_name for c in image.channel_map}
    columns = {}
    for spec in registry.markers:
        if spec.name in names:
            columns[spec.name] = f"{spec.compartment}: {names[spec.name]}: Mean"
    present = set(pq.ParquetFile(cells_path).schema_arrow.names)
    wanted = ["Object ID"] + sorted({c for c in columns.values() if c in present})
    table = pq.read_table(cells_path, columns=wanted).to_pandas().set_index("Object ID")
    table.index = table.index.astype(str)
    table = table.reindex(object_ids)
    return {m: table[c] for m, c in columns.items() if c in table}


def run(task: Task) -> dict:
    import numpy as np
    import pandas as pd

    from ..artifacts import QuPathCohortManifest
    from ..cells_parquet import build_cells_parquet
    from ..config import load_config
    from ..donor_pipeline import _geometry_config
    from ..expression import build_selected_expression, read_single_partition, write_donor_partition
    from ..geometry_stage import run_geometry_qc_for_image
    from ..marker_registry import registry_from_dict
    from ..redsea import RedseaParams, run_redsea

    spec = json.loads(task.input("cohort_spec").read_text())
    manifest = QuPathCohortManifest.read_json(task.input("qupath_manifest"))
    images = _measured(manifest, spec)
    wanted = {i["image_id"] for i in spec.get("images", [])}
    missing = sorted(wanted - {i.image_id for i in manifest.images})
    extra = sorted(i.image_id for i in manifest.images if i not in images)
    if extra:
        # build_cells_parquet ingests every image the manifest names.
        raise MeasurementError(f"the QuPath manifest holds images this stage does not measure: "
                               f"{extra[:5]}")

    payload = json.loads(task.input("marker_registry").read_text())
    ruling = _ruling(task.input("compartment_ruling"))
    own = own_compartment_registry(payload, ruling)
    (task.output(f"{DATA}/registry_own_compartment.json")).write_text(
        json.dumps(own, indent=2) + "\n")
    registry = registry_from_dict(own)

    cfg = load_config(task.input("config_ini"), data_dir=task.staging / DATA,
                      qupath_manifest=task.input("qupath_manifest"),
                      marker_registry=task.input("marker_registry"), n_jobs=task.settings["jobs"])
    build_cells_parquet(cfg)
    context = SimpleNamespace(config=cfg)
    for image in images:
        result = run_geometry_qc_for_image(read_single_partition(cfg.cells_dir, image.donor_id),
                                           image, config=_geometry_config(context, image))
        write_donor_partition(result.to_frame(), cfg.geometry_qc_dir, donor_id=image.donor_id)
    donors = [i.donor_id for i in images]
    run_redsea(cfg, donors, RedseaParams.from_config(cfg), n_jobs=cfg.n_jobs)
    if cfg.redsea_scratch.exists():           # its masks and intermediates are not outputs
        shutil.rmtree(cfg.redsea_scratch)
    scratch = cfg.data_dir / "scratch"
    if scratch.exists() and not any(scratch.iterdir()):
        scratch.rmdir()

    admission, markers, problems, unmatched, reasons = [], [], [], [], Counter()
    for image in images:
        d = image.donor_id
        cells = read_single_partition(cfg.cells_dir, d)
        redsea = read_single_partition(cfg.redsea_dir, d)
        qc = read_single_partition(cfg.geometry_qc_dir, d)
        excluded, table = set(), None
        if f"exclusions:{d}" in task.inputs:
            table = pd.read_parquet(task.input(f"exclusions:{d}"),
                                    columns=["object_id", "analysis_eligible",
                                             "analysis_exclusion_reason"])
            out = table.loc[~table["analysis_eligible"].astype(bool)]
            excluded = set(out["object_id"].astype(str))
            known = set(cells["object_id"].astype(str))
            stray = sorted(excluded - known)
            if stray:
                unmatched.append(f"{d}: {len(stray)} excluded cell(s) the export does not have, "
                                 f"e.g. {stray[0]}")
            reasons.update(f"{r}" for r in out["analysis_exclusion_reason"].fillna("(none)"))
        qc_x, hit = _apply_exclusions(qc, excluded)
        expression, availability = build_selected_expression(
            cells, redsea, qc_x, donor_id=d,
            acquired_markers={c.marker for c in image.channel_map}, registry=registry)
        write_donor_partition(expression, task.staging / DATA / "expression", donor_id=d)
        write_donor_partition(availability, task.staging / DATA / "availability", donor_id=d)

        ids = {name: frame["object_id"].astype(str) for name, frame in
               (("cells", cells), ("geometry_qc", qc), ("redsea", redsea),
                ("expression", expression))}
        for name, series in ids.items():
            if series.duplicated().any():
                problems.append(f"{d}: {name} repeats {int(series.duplicated().sum())} object id(s)")
        reference = set(ids["cells"])
        for name, series in ids.items():
            if set(series) != reference:
                problems.append(f"{d}: {name}'s object ids are not the cells'")

        eligible = expression["qc_analysis_eligible"].astype(bool).to_numpy()
        ratio = qc.set_index(qc["object_id"].astype(str))["nucleus_cell_area_ratio"].reindex(
            expression["object_id"].astype(str)).to_numpy(float)
        geometry_eligible = qc["analysis_eligible"].astype(bool).to_numpy()
        admission.append(dict(
            donor_id=d, cells=len(cells), geometry_admitted=int(geometry_eligible.sum()),
            excluded_by_ruling=int((geometry_eligible & hit).sum()),
            admitted=int(eligible.sum()),
            admitted_without_nucleus=int((eligible & ~(np.isfinite(ratio) & (ratio > 0))).sum()),
            **{f"not_admitted: {k}": v for k, v in Counter(
                qc.loc[~qc["analysis_eligible"].astype(bool), "analysis_reason"].fillna("(none)")
            ).items()}))

        before = _before(task.input(f"cells:{d}"), image, registry,
                         expression["object_id"].astype(str).to_numpy())
        available = dict(zip(availability["marker"], availability["available"].astype(bool),
                             strict=True))
        for spec_ in registry.markers:
            after = pd.to_numeric(expression.get(spec_.name), errors="coerce")
            row = dict(donor_id=d, marker=spec_.name, compartment=spec_.compartment,
                       spillover_policy=spec_.spillover_policy,
                       available=bool(available.get(spec_.name, False)), admitted=int(eligible.sum()))
            if after is not None and row["available"]:
                a = after.to_numpy(float)[eligible]
                b = (before[spec_.name].to_numpy(float)[eligible] if spec_.name in before
                     else np.full(a.shape, np.nan))
                row.update(finite=int(np.isfinite(a).sum()),
                           mean_before=float(np.nanmean(b)) if np.isfinite(b).any() else np.nan,
                           mean_after=float(np.nanmean(a)) if np.isfinite(a).any() else np.nan,
                           median_before=float(np.nanmedian(b)) if np.isfinite(b).any() else np.nan,
                           median_after=float(np.nanmedian(a)) if np.isfinite(a).any() else np.nan)
            markers.append(row)

    admission_frame = pd.DataFrame(admission).fillna(0)
    marker_frame = pd.DataFrame(markers)
    admission_frame.to_csv(task.output("nucleus_admission.csv"), index=False)
    marker_frame.to_csv(task.output("markers_before_after.csv"), index=False)
    _admission_figure(task, admission_frame)
    _marker_figure(task, marker_frame)

    measured = {name: sorted(p.name.split("=", 1)[1] for p in (task.staging / DATA / name).glob(
        "donor_id=*")) for name in TABLES}
    compartments = Counter(m["compartment"] for m in own["markers"])
    unavailable = (marker_frame.loc[~marker_frame["available"]].groupby("donor_id")["marker"]
                   .apply(lambda s: ", ".join(sorted(s))).to_dict())
    return {
        "cohort_images_not_in_manifest": missing,
        "donors": donors,
        "donors_not_measured": {name: sorted(set(donors) - set(found))
                                for name, found in measured.items() if set(donors) - set(found)},
        "object_id_problems": problems,
        "admitted_without_nucleus": {r["donor_id"]: r["admitted_without_nucleus"] for r in admission
                                     if r["admitted_without_nucleus"]},
        "missing_geometry_policy": cfg.cell_qc_missing_geometry_policy,
        "exclusions_unmatched": unmatched,
        "exclusions_by_reason": dict(reasons.most_common()),
        "excluded_by_ruling": {r["donor_id"]: r["excluded_by_ruling"] for r in admission},
        "ruling_problems": ruling_problems(payload, ruling),
        "compartments": dict(sorted(compartments.items())),
        "ruled_markers": {c: sorted(m for m, v in ruling.items() if v == c) for c in RULED},
        "admission": {r["donor_id"]: {k: r[k] for k in ("cells", "geometry_admitted",
                                                        "excluded_by_ruling", "admitted")}
                      for r in admission},
        "markers_unavailable": unavailable,
        "redsea": dict(task.settings["config"].get("redsea", {})),
        "environment_overrides": task.settings["environment_overrides"],
    }


# ------------------------------------------------------------------------------ figures ----
#: The fills of the admission figure: what happened to a cell, never how many.
ADMISSION_FILLS = {"admitted": "#7fb27f", "excluded_by_ruling": "#d99a4e",
                   "no nucleus": "#9aa6b2", "other geometry QC": "#cfd6dd"}


def _admission_figure(task: Task, frame) -> None:
    """Each donor's cells: admitted, excluded by the recorded ruling, or not admitted by geometry
    QC (no nucleus, or another reason)."""
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    frame = frame.sort_values("donor_id")
    reasons = [c for c in frame.columns if c.startswith("not_admitted: ")]
    no_nucleus = [c for c in reasons if "nucleus" in c]
    parts = {"admitted": frame["admitted"],
             "excluded_by_ruling": frame["excluded_by_ruling"],
             "no nucleus": frame[no_nucleus].sum(axis=1) if no_nucleus else 0 * frame["cells"],
             "other geometry QC": frame[[c for c in reasons if c not in no_nucleus]].sum(axis=1)
             if len(reasons) > len(no_nucleus) else 0 * frame["cells"]}
    fig, ax = plt.subplots(figsize=(10, 0.32 * len(frame) + 1.6))
    left = 0 * frame["cells"]
    for kind, values in parts.items():
        ax.barh(frame["donor_id"], values / 1e6, left=left / 1e6, color=ADMISSION_FILLS[kind],
                label=kind)
        left = left + values
    ax.invert_yaxis()
    ax.set_xlabel("cells (millions)")
    ax.set_title("Nucleus admission by donor (stage 03): every cell, admitted or why not",
                 fontsize=10)
    ax.legend(loc="lower right", fontsize=8, frameon=False)
    fig.tight_layout()
    fig.savefig(gates.figure_path(task, "nucleus_admission.png"), dpi=150)
    plt.close(fig)


def _marker_figure(task: Task, frame) -> None:
    """For each marker and donor, how REDSEA changed its mean over the admitted cells in the
    compartment it is measured in: log2 of the mean after over the mean before. Grey where the
    marker was not acquired, or no value could be compared."""
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import numpy as np

    available = frame.loc[frame["available"]].copy()
    with np.errstate(divide="ignore", invalid="ignore"):
        available["log2_change"] = np.log2(available["mean_after"] / available["mean_before"])
    grid = available.pivot(index="marker", columns="donor_id", values="log2_change")
    order = (frame.drop_duplicates("marker").set_index("marker")["compartment"])
    grid = grid.reindex(index=sorted(order.index, key=lambda m: (order[m], m.casefold())))
    labels = [f"{m} ({order[m]})" for m in grid.index]
    finite = np.abs(grid.to_numpy(float)[np.isfinite(grid.to_numpy(float))])
    span = float(np.nanpercentile(finite, 99)) if finite.size else 1.0
    span = max(span, 0.1)
    fig, ax = plt.subplots(figsize=(0.32 * max(len(grid.columns), 4) + 3, 0.16 * len(grid) + 1.8))
    cmap = plt.get_cmap("RdBu_r").with_extremes(bad="#d9d9d9")
    image = ax.imshow(np.ma.masked_invalid(grid.to_numpy(float)), aspect="auto", cmap=cmap,
                      vmin=-span, vmax=span)
    ax.set_xticks(range(len(grid.columns)), grid.columns, rotation=90, fontsize=7)
    ax.set_yticks(range(len(grid.index)), labels, fontsize=6)
    ax.set_title("REDSEA's change per marker and donor (stage 03): log2 mean after / mean "
                 "before, admitted cells,\nthe compartment the marker is measured in; grey = not "
                 "acquired or not comparable", fontsize=9)
    fig.colorbar(image, ax=ax, fraction=0.03, pad=0.01, label="log2 after / before")
    fig.tight_layout()
    fig.savefig(gates.figure_path(task, "markers_before_after.png"), dpi=150)
    plt.close(fig)


# --------------------------------------------------------------------------------- gate ----
def gate(task: Task, summary: dict) -> gates.Gate:
    C = gates.Check
    checks = [
        C("every cohort image is in the QuPath manifest", summary["cohort_images_not_in_manifest"],
          "none", not summary["cohort_images_not_in_manifest"]),
        C("every measured donor has all five tables", summary["donors_not_measured"], "none",
          not summary["donors_not_measured"],
          detail="cells, geometry QC, REDSEA, expression and availability, one partition each"),
        C("object ids are unique and the same across the stage's tables",
          summary["object_id_problems"], "none", not summary["object_id_problems"]),
        C("no admitted cell lacks a nucleus", summary["admitted_without_nucleus"], "none",
          not summary["admitted_without_nucleus"]
          and summary["missing_geometry_policy"] == "exclude",
          detail=f"missing_geometry_policy = {summary['missing_geometry_policy']}; a cell with no "
                 "nucleus is excluded from every calculation"),
        C("every recorded exclusion names a cell the export has", summary["exclusions_unmatched"],
          "none", not summary["exclusions_unmatched"]),
        C("the compartment ruling names only registry markers and compartments",
          summary["ruling_problems"], "none", not summary["ruling_problems"]),
        C("compartments measured", summary["compartments"], "reported", True, blocking=False,
          detail="marker_compartment_ruling.csv sets each ruled marker's; the rest keep the "
                 "registry's"),
        C("markers the ruling moves, by compartment", summary["ruled_markers"], "reported", True,
          blocking=False),
        C("nucleus admission by donor", summary["admission"], "reported", True, blocking=False,
          detail="cells; admitted by geometry QC; taken out by the recorded exclusions; admitted"),
        C("recorded exclusions, by reason", summary["exclusions_by_reason"], "reported", True,
          blocking=False, detail="every cell the 2026-10-05 table marks not analysis-eligible; "
                                 "beyond geometry QC's own reasons, the 11,484 reviewed cells"),
        C("markers not acquired, by donor", summary["markers_unavailable"], "reported", True,
          blocking=False, detail="unavailable, never zero"),
        C("REDSEA settings", summary["redsea"], "reported", True, blocking=False),
        C("settings an environment variable changed", summary["environment_overrides"],
          "reported", True, blocking=False),
    ]
    return gates.Gate(STAGE, "Protein measurement: one value per cell and marker, in the "
                      "compartment the calibration review ruled, after REDSEA.",
                      checks, ["nucleus_admission.png", "markers_before_after.png"])
