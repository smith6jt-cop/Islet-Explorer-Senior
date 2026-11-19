"""Cell Typing for the spleen primary cohort (``typing_route: spleen_protein``).

Maintainer, 2026-10-30 and 2026-11-01:

- **Sections:** every PC67 and PC33 section of the accepted Protein Measurement. Each panel uses its own
  registry and marker availability; level selection never pools sections from different panels.
- **Positivity** is binary, from each section's own cells: a marker's level is its known-negative percentile in the
  spleen lineage that cannot carry it (``lineage_levels``); interior coverage (``coverage``) serves the tiebreak.
- **Names** come from the 6551 spleen hierarchy (``tissues.spleen.load_rule_source``) through
  ``positivity.hierarchy``, never from the CODEX identity step. Endothelial is structure-only and waits for a vessel
  object; a cell matching two branches is settled by the CODEX canary rule or stays ambiguous.
- **Normalization** is a prerequisite (its acceptance is read), not an input: typing reads each marker's REDSEA value
  against that section's own level.
- **The gate is the first manual stop.**

Steps (maintainer 2026-11-01: the spleen registry starts over, measured on spleen, nothing carried from another
dataset):
1. levels per section from mutually exclusive lineages (``positivity.lineage_levels``): B and T cells from the data,
   each marker's known-negative percentile in the population that cannot carry it;
2. canary levels (in spleen the canary null is that same population, so the canary level is the level);
3. area coverage per section, on the section's window of its qptiff (the cell geometry the manifest names);
4. calls, the hierarchy and the tiebreak.

A required marker without a passing level on every section stops the run before naming: the method is repaired for
every section, not excused on one.

Outputs: ``data/cells/donor_id=<section>``, ``data/coverage/<section>.parquet``, ``tables/*.csv`` (options,
selections, levels, canary levels) and the gate figures.
"""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd

from .. import gates, roots
from ..tasks import Input, Task

STAGE = "06_typing"
ROUTE = "spleen_protein"
REGISTRY = ("phenocycler-pc67", "spleen")  # Default for existing PC67 review adapters.
TYPED_PANEL = "spleen_pc67"
REGISTRIES = {TYPED_PANEL: REGISTRY, "spleen_pc33": ("phenocycler-pc33", "spleen")}
PANEL_DIR = roots.REPO / "core/panels/data/tissues"
PANEL_FILES = {"panel_manifest": PANEL_DIR / "spleen/manifest.yml",
               "panel_header": PANEL_DIR / "spleen/pc67_export_header.tsv",
               "panel_header_pc33": PANEL_DIR / "spleen/pc33_export_header.tsv",
               "pln_manifest": PANEL_DIR / "pancreatic_lymph_node/manifest.yml",
               "rule_sheet": PANEL_DIR / "pancreatic_lymph_node/phenotype_rules.csv"}
CONTEXT = "typing_context"       # site source: <dir>/<section>.parquet (object_id, context, start_object)
WINDOW_MARGIN_PX = 16
TABLES = ("options", "selection", "differentiation_selection", "levels", "levels_differentiation")


class SpleenTypingError(ValueError):
    """The accepted inputs cannot be typed as declared."""


# --------------------------------------------------------------------------------- sections ----
def _spec(project) -> Path:
    from .panel_normalization import _spec as spec

    return spec(project)


def _panels(project) -> dict[str, list[str]]:
    from .panel_normalization import _panels as panels

    return panels(project)


def _manifest_path(project) -> Path:
    named = project.source("qupath_manifest")
    if not named:
        raise SpleenTypingError(f"{STAGE} needs the qupath_manifest source")
    return Path(named)


def _manifest(project):
    from ..artifacts import QuPathCohortManifest

    return QuPathCohortManifest.read_json(_manifest_path(project))


def _typed_panels(project) -> dict[str, list[str]]:
    panels = _panels(project)
    if not panels or set(panels) - REGISTRIES.keys():
        raise SpleenTypingError(f"the cohort spec names unsupported spleen panels: {sorted(set(panels) - REGISTRIES.keys())}")
    return {p: panels[p] for p in sorted(panels, key=lambda p: p != "spleen_pc33")}


def _typed(project) -> list[str]:
    return [tag for tags in _typed_panels(project).values() for tag in tags]


def _slides(project) -> dict[str, list]:
    """{qptiff path: manifest images on it} for every slide holding a typed section."""
    typed = set(_typed(project))
    slides: dict[str, list] = {}
    for image in _manifest(project).images:
        slides.setdefault(str(Path(image.qptiff.path).resolve()), []).append(image)
    return {q: imgs for q, imgs in slides.items() if typed & {i.donor_id for i in imgs}}


def _source_dir(project, name: str) -> Path:
    named = project.source(name)
    if not named:
        raise SpleenTypingError(f"{STAGE} needs the {name} source; name it in the site file")
    return Path(named)


def _partition(release: Path, table: str, donor: str) -> Path:
    files = sorted((release / "03_protein_measurement/data" / table / f"donor_id={donor}").glob("*.parquet"))
    if len(files) != 1:
        raise SpleenTypingError(f"{donor}: expected one Protein Measurement {table} partition, found {len(files)}")
    return files[0]


def settings(project) -> dict:
    from ..tissues import spleen

    panels = _typed_panels(project)
    return {"registries": {p: list(REGISTRIES[p]) for p in panels}, "panels": panels, "typed": _typed(project),
            "not_typed": {}, "workers": int(project.settings.get("typing_workers", 4)),
            "naming_policy": {p: spleen.naming_policy(PANEL_FILES["panel_manifest"], panel_id=p) for p in panels},
            "method": "per-panel, per-image levels, canary, coverage (positivity) -> 6551 spleen hierarchy (hierarchy)"}


def inputs(project) -> list[Input]:
    release = project.release_dir
    found = [Input("protein_receipt", release / "03_protein_measurement/receipt.json"),
             Input("protein_acceptance", release / "03_protein_measurement/gate/accepted.json"),
             Input("normalization_receipt", release / "05_normalization/receipt.json"),
             Input("normalization_acceptance", release / "05_normalization/gate/accepted.json"),
             Input("cohort_spec", _spec(project)), Input("qupath_manifest", _manifest_path(project))]
    found += [Input(name, path) for name, path in PANEL_FILES.items()]
    from ..positivity.registry import REGISTRY_DIR

    for panel_id in _typed_panels(project):
        registry_name = "__".join(REGISTRIES[panel_id])
        for suffix in ("", "__constants", "__proposals"):
            path = REGISTRY_DIR / f"{registry_name}{suffix}.csv"
            if suffix != "__proposals" or path.exists():
                found.append(Input(f"registry:{panel_id}{suffix}", path))
    provenance = roots.REPO / "_build_provenance.json"
    if provenance.is_file():
        found.append(Input("build_provenance", provenance))
    context = _source_dir(project, CONTEXT)
    for d in _typed(project):
        found += [Input(f"redsea:{d}", _partition(release, "redsea", d)),
                  Input(f"expression:{d}", _partition(release, "expression", d)),
                  Input(f"context:{d}", context / f"{d}.parquet")]
    for qptiff, images in _slides(project).items():
        stem = Path(qptiff).name.split(".")[0]
        found.append(Input(f"qptiff:{stem}", Path(qptiff), hash=False, recorded_sha256=images[0].qptiff.sha256))
        for image in images:
            if image.donor_id in _typed(project):
                found.append(Input(f"cell_geojson:{image.donor_id}", Path(image.cell_geojson.path), hash=False,
                                   recorded_sha256=image.cell_geojson.sha256))
    return found


# --------------------------------------------------------------------------------- steps ----
def _check_prerequisites(task: Task) -> None:
    spec = json.loads(task.input("cohort_spec").read_text())
    donors = sorted(str(i["donor_id"]) for i in spec["images"])
    receipt = json.loads(task.input("protein_receipt").read_text())
    if (receipt.get("stage") != "03_protein_measurement" or receipt.get("status") != "complete"
            or sorted(receipt.get("summary", {}).get("donors", [])) != donors):
        raise SpleenTypingError("the accepted Protein Measurement roster differs from the cohort spec")
    norm = json.loads(task.input("normalization_receipt").read_text())
    if norm.get("stage") != "05_normalization" or norm.get("status") != "complete":
        raise SpleenTypingError("Normalization is not complete")
    for name in ("protein_acceptance", "normalization_acceptance"):
        if not json.loads(task.input(name).read_text()).get("accepted_by"):
            raise SpleenTypingError(f"{name}: the gate is not accepted")


def _window(features: list[dict], shape: tuple[int, int]) -> tuple[int, int, int, int]:
    """(x0, y0, x1, y1): the cells' bounding box plus WINDOW_MARGIN_PX, inside the image."""
    xs, ys = [], []
    for f in features:
        g = f["geometry"]
        polys = g["coordinates"] if g["type"] == "MultiPolygon" else [g["coordinates"]]
        for poly in polys:
            ring = np.asarray(poly[0], float)
            xs += [ring[:, 0].min(), ring[:, 0].max()]
            ys += [ring[:, 1].min(), ring[:, 1].max()]
    H, W = shape
    m = WINDOW_MARGIN_PX
    return (max(0, int(np.floor(min(xs))) - m), max(0, int(np.floor(min(ys))) - m),
            min(W, int(np.ceil(max(xs))) + m + 1), min(H, int(np.ceil(max(ys))) + m + 1))


def _shifted(features: list[dict], x0: int, y0: int) -> list[dict]:
    def move(ring):
        return [[x - x0, y - y0] for x, y in ring]

    out = []
    for f in features:
        g = f["geometry"]
        if g["type"] == "Polygon":
            geom = {"type": "Polygon", "coordinates": [move(r) for r in g["coordinates"]]}
        else:
            geom = {"type": "MultiPolygon", "coordinates": [[move(r) for r in p] for p in g["coordinates"]]}
        out.append({"id": f.get("id"), "geometry": geom})
    return out


def _coverage_channel_index(image, array_shape: tuple[int, ...], axes: str) -> dict[str, int]:
    """Reject a channel map that could silently read the wrong OME-TIFF/qptiff plane."""
    if len(array_shape) != 3 or axes[-2:] != "YX" or any(n <= 0 for n in array_shape):
        raise SpleenTypingError(f"{image.donor_id}: coverage requires channel/Y/X image data, found {axes} {array_shape}")
    names = [c.channel_name for c in image.channel_map]
    indices = [c.channel_index for c in image.channel_map]
    if (not names or len(set(names)) != len(names) or len(set(indices)) != len(indices)
            or any(i < 0 or i >= array_shape[0] for i in indices)):
        raise SpleenTypingError(f"{image.donor_id}: channel map does not fit the image shape {array_shape}")
    return dict(zip(names, indices))


def _coverage(task: Task, image, tag: str, levels: pd.DataFrame, reg, redsea_ids) -> pd.DataFrame:
    """coverage.compute on the section's window of its qptiff, from the cell geometry the manifest names."""
    import tifffile
    import zarr

    from ..artifacts import iter_qupath_geojson_features
    from ..positivity import coverage as AC

    feats = [{"id": f.get("id"), "geometry": f["geometry"]}
             for f in iter_qupath_geojson_features(task.input(f"cell_geojson:{tag}"))]
    qptiff = Path(image.qptiff.path)
    with tifffile.TiffFile(task.input(f"qptiff:{qptiff.name.split('.')[0]}")) as tf:
        z = zarr.open(tf.series[0].aszarr(level=0), mode="r")
        index = _coverage_channel_index(image, z.shape, tf.series[0].axes)
        H, W = z.shape[-2:]
        x0, y0, x1, y1 = _window(feats, (H, W))

        def read(name):
            return np.asarray(z[index[name], y0:y1, x0:x1])

        frame = AC.compute(tag, _shifted(feats, x0, y0), (y1 - y0, x1 - x0), read, list(index), levels, reg)
    # every geometry object is drawn (each cell's pixels are its own), but the table keeps the canonical cells, as
    # REDSEA does: the GeoJSON can hold objects the measurement table does not
    frame = frame[frame["object_id"].isin(set(redsea_ids))].reset_index(drop=True)
    frame.to_parquet(task.output(f"data/coverage/{tag}.parquet"), index=False)      # the layout the source reads
    return frame


def _panel_tables(tables: dict, panel_id: str) -> dict:
    """Allow review adapters to read either aggregate or unchanged panel-local tables."""
    return {name: value.loc[value["panel_id"] == panel_id].drop(columns="panel_id").copy()
            if isinstance(value, pd.DataFrame) and "panel_id" in value else value
            for name, value in tables.items()}


def _named(task: Task, source, tag: str, image, reg, tables: dict, canary: pd.DataFrame, rules, branches,
           *, panel_id: str = TYPED_PANEL) -> pd.DataFrame:
    from ..identity import cell_uid
    from ..positivity import coverage as AC
    from ..positivity import hierarchy as H
    from ..positivity import identity as ID
    from ..tissues import spleen

    tables = _panel_tables(tables, panel_id)
    if "panel_id" in canary:
        canary = canary.loc[canary["panel_id"] == panel_id].drop(columns="panel_id")
    sel = tables["selection"].set_index("marker")
    used = [m for m in reg.table.index if m in sel.index and sel.loc[m, "decision"] == "use"]
    lv = tables["levels"][tables["levels"]["image"] == tag].set_index("marker")
    ld = tables["levels_differentiation"]
    ld = ld[ld["image"] == tag].set_index("marker")
    cl = canary[canary["image"] == tag].set_index("marker")
    typing = [m for m in ID.used_markers(tables["selection"], reg) if m in cl.index]
    names = H.sheet_names(rules, used)
    cols = sorted({reg.column(m, lv.loc[m, "measurement"]) for m in used}
                  | {reg.column(m, ld.loc[m, "measurement"]) for m in typing if m in ld.index})
    d = source.identity_frame(tag, cols, [])
    calls = H.marker_calls(d, reg, lv, used)
    policy = spleen.naming_policy(task.input("panel_manifest"), panel_id=panel_id)
    ev = H.evaluate_profiles(rules, {names[m]: calls[m] for m in used if m in names}, branches,
                             defer_exclusions=policy["defer_exclusions"])

    ev_markers = [m for m in typing if m in names]
    value = lambda m, meas: np.log1p(np.clip(np.nan_to_num(d[reg.column(m, meas)].to_numpy(float), nan=0.0), 0, None))
    V = np.column_stack([value(m, lv.loc[m, "measurement"]) for m in ev_markers])
    L = lv.loc[ev_markers, "level_log"].to_numpy(float)
    Lc = np.maximum(cl.loc[ev_markers, "canary_log"].to_numpy(float), L)
    F = V - L
    for j, m in enumerate(ev_markers):          # brightness on the differentiation reading, as identity reads it
        if m in ld.index and reg.differentiation(m):
            F[:, j] = value(m, ld.loc[m, "measurement"]) - float(ld.loc[m, "level_log"])
    cov = AC.read_validated(source.coverage(tag), source.redsea_ids(tag), tag, d["object_id"], ev_markers)
    evidence = {"columns": [names[m] for m in ev_markers], "P": np.column_stack([calls[m] == 1 for m in ev_markers]),
                "C": V > Lc, "F": F, "COV": cov[[AC.INTERIOR_PREFIX + m for m in ev_markers]].to_numpy(float),
                "interior_valid": cov[AC.INTERIOR_COUNT].to_numpy(np.int64) > 0, "typing": [names[m] for m in ev_markers]}
    out = H.name_cells(rules, ev, evidence, reg.const("CANARY_MARGIN"), reg.const("CANARY_MIN_DEPTH"),
                       reg.const("CANARY_TIE_INTERIOR_LEAD"), require_specific=policy["require_specific"],
                       inherit_depth=policy["inherit_depth"])
    out.insert(0, "cell_uid", [cell_uid(image.image_id, image.segmentation_version, str(o)) for o in d["object_id"]])
    out.insert(1, "object_id", d["object_id"].to_numpy())
    out.insert(2, "start_class", spleen.start_class(d))
    out["violations"] = _below_level(rules, out["refined_class"].to_numpy(object), calls, names)
    out["unsupported_t_cell"] = _unsupported_t_cells(rules, out["refined_class"].to_numpy(object), calls)
    for m in used:
        out[f"call__{m}"] = calls[m]
    for j, m in enumerate(ev_markers):
        out[f"canary__{m}"] = evidence["C"][:, j]
    out.insert(0, "donor_id", tag)
    out.to_parquet(task.output(f"data/cells/donor_id={tag}/data.parquet"), index=False)
    return out


def _below_level(rules, refined: np.ndarray, calls: dict, names: dict) -> np.ndarray:
    """Per cell, how many of its class's positive clauses (its rule and every ancestor) are not met by a marker above
    its level: 0 for every named cell, or the naming is wrong."""
    by_sheet = {names[m]: calls[m] for m in calls if m in names}
    by = rules.by_phenotype
    out = np.zeros(len(refined), np.int16)
    for cls in set(refined) - {""}:
        rows = refined == cls
        chain, current = [], cls
        while current in by:
            chain.append(by[current])
            current = by[current].parent
        for rule in chain:
            anyp = [m for m, t in rule.constraints if t == "anypos"]
            for m, t in rule.constraints:
                if t in ("pos", "allpos"):
                    out[rows] += by_sheet.get(m, np.full(len(refined), -1))[rows] != 1
            if anyp:
                hit = np.logical_or.reduce([by_sheet.get(m, np.zeros(len(refined)))[rows] == 1 for m in anyp])
                out[rows] += ~hit
    return out


def _unsupported_t_cells(rules, refined: np.ndarray, calls: dict) -> np.ndarray:
    """Independently check marker support for every final T descendant, whatever its label."""
    by_name, t_classes = rules.by_phenotype, set()
    for name in set(refined) - {""}:
        current = name
        while current in by_name:
            if current == "T_cell":
                t_classes.add(name)
                break
            current = by_name[current].parent
    absent = np.full(len(refined), -1, dtype=np.int8)
    supported = ((calls.get("CD3e", absent) == 1)
                 & ((calls.get("CD4", absent) == 1) | (calls.get("CD8", absent) == 1)))
    return np.isin(refined, list(t_classes)) & ~supported


def _write_panel_table(task: Task, name: str, frames: dict[str, pd.DataFrame]) -> None:
    """Panel-local CSVs retain their original schema; aggregates label their panel explicitly."""
    for panel_id, frame in frames.items():
        frame.to_csv(task.output(f"tables/by_panel/{panel_id}/{name}.csv"), index=False)
    aggregate = (next(iter(frames.values())) if len(frames) == 1 else
                 pd.concat([frame.assign(panel_id=panel_id) for panel_id, frame in frames.items()], ignore_index=True))
    aggregate.to_csv(task.output(f"tables/{name}.csv"), index=False)


def _file_record(path: Path) -> dict:
    from ..artifacts import sha256_file

    return {"path": str(path.resolve()), "bytes": path.stat().st_size, "sha256": sha256_file(path)}


def _review_code_identity() -> dict:
    from ..tasks import code_identity
    from ..positivity import coverage, hierarchy, identity, lineage_levels, registry, sources
    from ..tissues import lymph_node, spleen

    code = code_identity()
    code["files"] = {module.__name__: _file_record(Path(module.__file__))
                     for module in (coverage, hierarchy, identity, lineage_levels, registry, sources, lymph_node, spleen)}
    code["files"][__name__] = _file_record(Path(__file__))
    return code


def _write_section_review_ready(task: Task, tag: str, panel_id: str, cells: pd.DataFrame,
                                panel_summary: dict, code: dict) -> Path | None:
    """An atomic, explicitly provisional handoff; never a complete-stage receipt or acceptance."""
    import pyarrow.parquet as pq

    if not code.get("commit") or code.get("dirty") is not False:
        return None
    protein_receipt = task.input("protein_receipt")
    receipt = json.loads(protein_receipt.read_text())

    def accepted_measurement_record(name):
        path = task.input(f"{name}:{tag}")
        accepted = next((r for r in receipt.get("outputs", [])
                         if (protein_receipt.parent / r["path"]).resolve() == path.resolve()), None)
        actual = _file_record(path)
        if not accepted or any(actual[k] != accepted.get(k) for k in ("sha256", "bytes")):
            raise SpleenTypingError(f"{tag}: {name} input differs from the accepted Protein Measurement receipt")
        return actual

    expression = task.input(f"expression:{tag}")
    expression_record = accepted_measurement_record("expression")
    redsea_record = accepted_measurement_record("redsea")
    eligible = pq.ParquetFile(expression).read(columns=["object_id", "qc_analysis_eligible"]).to_pandas()
    if eligible["qc_analysis_eligible"].isna().any():
        raise SpleenTypingError(f"{tag}: missing analysis eligibility in accepted Protein Measurement")
    expected = pd.Index(eligible.loc[eligible["qc_analysis_eligible"], "object_id"])
    actual = pd.Index(cells["object_id"])
    uuid_match = (not expected.hasnans and not actual.hasnans and expected.is_unique and actual.is_unique
                  and len(expected) == len(actual) and expected.difference(actual).empty)
    per = _section_summary(cells)
    required_stopped = panel_summary["markers"]["required_stopped"]
    checks = [
        {"name": "every required marker has a passing level on every section", "value": required_stopped,
         "passed": not required_stopped},
        {"name": "cells typed", "value": len(cells), "passed": bool(uuid_match)},
        {"name": "cells named by a marker below its level", "value": per["named_below_level"],
         "passed": per["named_below_level"] == 0},
        {"name": "forbidden broad or double-negative T final calls", "value": per["forbidden_final_classes"],
         "passed": per["forbidden_final_classes"] == 0},
        {"name": "T-cell calls without CD3e and CD4-or-CD8 support", "value": per["unsupported_t_cells"],
         "passed": per["unsupported_t_cells"] == 0},
    ]
    if not all(check["passed"] for check in checks) or not code.get("commit"):
        return None
    panel_dir = task.staging / "tables/by_panel" / panel_id
    outputs = {"cells": _file_record(task.staging / f"data/cells/donor_id={tag}/data.parquet"),
               "coverage": _file_record(task.staging / f"data/coverage/{tag}.parquet")}
    outputs.update({f"tables/{name}": _file_record(panel_dir / f"{name}.csv")
                    for name in (*TABLES, "canary_levels")})
    inputs = set(PANEL_FILES) | {"protein_receipt", "protein_acceptance", "normalization_receipt",
                               "normalization_acceptance", "cohort_spec", "qupath_manifest", f"context:{tag}"}
    inputs.update(name for name in task.inputs if name.startswith(f"registry:{panel_id}") or name == "build_provenance")
    record = {
        "schema": "spleen-section-review-ready/1", "status": "provisional", "section": tag, "donor_id": tag,
        "panel": panel_id, "cells": len(cells), "code_commit": code["commit"],
        "code": {k: v for k, v in code.items() if k != "files"}, "code_files": code["files"],
        "attempt_id": task.staging.name, "attempt_dir": str(task.staging.resolve()),
        "outputs": outputs, "inputs": {name: _file_record(task.input(name)) for name in sorted(inputs)},
        "panel_adaptations": _file_record(panel_dir / "rule_adaptations.json"),
        "eligibility": {"expression_input": str(expression.resolve()), "expression_output": expression_record,
                        "redsea_output": redsea_record, "source_hashes_verified": True,
                        "contract": "qc_analysis_eligible == true", "expected_cells": len(expected),
                        "uuid_set_matches": bool(uuid_match)},
        "checks": checks,
        "note": "One finished donor for immediate review; full-stage publication and gate acceptance remain separate.",
    }
    path = task.output(f"review/ready/{tag}.json")
    temporary = path.with_suffix(".json.tmp")
    temporary.write_text(json.dumps(record, indent=2, sort_keys=True) + "\n")
    temporary.replace(path)
    print(f"Provisional review ready: {tag} ({len(cells):,} cells), {path}", flush=True)
    return path


def run(task: Task) -> dict:
    import matplotlib

    matplotlib.use("Agg")
    from ..positivity import lineage_levels as LL
    from ..positivity import registry as RG
    from ..positivity import sources as SRC
    from ..tissues import spleen

    _check_prerequisites(task)
    review_code = _review_code_identity()
    typed = task.settings["typed"]
    panels = {p: task.settings["panels"][p]
              for p in sorted(task.settings["panels"], key=lambda p: p != "spleen_pc33")}
    workers = task.settings["workers"]
    images = {i.donor_id: i for i in _manifest_from(task).images}
    if (set(typed) != {t for tags in panels.values() for t in tags}
            or set(typed) - images.keys()
            or any(images[tag].panel_id != panel_id for panel_id, tags in panels.items() for tag in tags)):
        raise SpleenTypingError("typed section/panel assignments disagree with the QuPath manifest")
    stage_dir = task.input("protein_receipt").parent
    source = SRC.ProteinMeasurementSource(stage_dir=stage_dir, context_dir=task.input(f"context:{typed[0]}").parent)
    summary = {"sections": typed, "not_typed": {}, "panels": {}, "stops": [], "lineages": {}}
    fitted = {}
    for panel_id, tags in panels.items():
        reg = RG.load(*REGISTRIES[panel_id])
        RG.require_reviewed(reg)
        rules = spleen.load_rule_source(task.input("panel_manifest"), panel_id=panel_id)
        branches = spleen.structure_only_branches(task.input("panel_manifest"), panel_id=panel_id)
        observed = spleen.observed_markers(task.input("panel_manifest"), panel_id=panel_id)
        tables = LL.choose(source, reg.spec, tags, workers=workers)
        fitted[panel_id] = (reg, rules, branches, tables)
        marker_summary = _marker_summary(tables["selection"], reg)
        summary["panels"][panel_id] = {
            "sections": tags, "registry": list(REGISTRIES[panel_id]), "markers": marker_summary,
            "observed_markers": sorted(observed),
            "unavailable_rules": sorted(name for name, gap in rules.missing_by_phenotype(observed).items() if gap),
            "structure_only": sorted(branches),
            "rule_adaptations": spleen.panel_adaptation_audit(task.input("panel_manifest"), panel_id=panel_id),
        }
        task.output(f"tables/by_panel/{panel_id}/rule_adaptations.json").write_text(
            json.dumps(summary["panels"][panel_id]["rule_adaptations"], indent=2, sort_keys=True) + "\n")
        summary["lineages"].update(tables["options"].groupby("image")[["b_cells", "t_cells"]].first().to_dict("index"))
        summary["stops"].extend(f"{panel_id}:{marker}" for marker in tables["stops"])
        _level_figures(task, source, reg, tables, tags)
    for name in TABLES:
        _write_panel_table(task, name, {p: values[3][name] for p, values in fitted.items()})
    summary["markers"] = {key: sorted({m for p in summary["panels"].values() for m in p["markers"][key]})
                          for key in ("use", "stop", "required_stopped")}
    summary["unavailable_rules"] = {p: s["unavailable_rules"] for p, s in summary["panels"].items()}
    summary["structure_only"] = {p: s["structure_only"] for p, s in summary["panels"].items()}
    if summary["markers"]["required_stopped"]:
        summary["typed"] = False
        return summary

    canaries = {p: LL.canary_levels(values[3], values[0].spec) for p, values in fitted.items()}
    _write_panel_table(task, "canary_levels", canaries)
    cells = {}
    covered = SRC.ProteinMeasurementSource(stage_dir=stage_dir, context_dir=source.context_dir,
                                           coverage_dir=task.staging / "data/coverage")
    for panel_id, tags in panels.items():
        reg, rules, branches, tables = fitted[panel_id]
        for tag in tags:
            _coverage(task, images[tag], tag, tables["levels"], reg, source.redsea_ids(tag))
            cells[tag] = _named(task, covered, tag, images[tag], reg, tables, canaries[panel_id], rules, branches,
                                panel_id=panel_id)
            _write_section_review_ready(task, tag, panel_id, cells[tag], summary["panels"][panel_id], review_code)
    summary["typed"] = True
    summary["per_section"] = {tag: _section_summary(c) for tag, c in cells.items()}
    _class_figure(task, cells)
    return summary


def _manifest_from(task: Task):
    from ..artifacts import QuPathCohortManifest

    return QuPathCohortManifest.read_json(task.input("qupath_manifest"))


def _marker_summary(selection: pd.DataFrame, reg) -> dict:
    sel = selection.set_index("marker")
    return {"use": sorted(m for m in sel.index if sel.loc[m, "decision"] == "use"),
            "stop": sorted(m for m in sel.index if sel.loc[m, "decision"] == "stop"),
            "required_stopped": sorted(m for m in sel.index if sel.loc[m, "decision"] == "stop"
                                       and m in reg.table.index and reg.table.loc[m, "required"] == "yes")}


def _section_summary(cells: pd.DataFrame) -> dict:
    return {"cells": int(len(cells)),
            "forbidden_final_classes": int(cells["refined_class"].isin(["Immune", "T_cell", "DN_T_cell"]).sum()),
            "unsupported_t_cells": int(cells["unsupported_t_cell"].sum()),

            "status": {k: int(v) for k, v in cells["status"].value_counts().items()},
            "cell_class": {k or "(none)": int(v) for k, v in cells["cell_class"].value_counts().items()},
            "refined_class": {k or "(none)": int(v) for k, v in cells["refined_class"].value_counts().items()},
            "tie_settled": int((cells["status"] == "tie_settled").sum()),
            "named_below_level": int((cells["violations"] > 0).sum()),
            "outside_regions": int((cells["start_class"] == "Tissue").sum())}


def _level_figures(task: Task, source, reg, tables: dict, typed: list[str]) -> None:
    """Per section, every typing marker's histogram (log1p of its positivity reading) with its chosen level."""
    import matplotlib.pyplot as plt

    from ..positivity import identity as ID

    markers = ID.used_markers(tables["selection"], reg) or list(tables["selection"]["marker"])
    lv = tables["levels"]
    for tag in typed:
        rows = lv[(lv["image"] == tag) & lv["marker"].isin(markers)].set_index("marker")
        cols = sorted({reg.column(m, rows.loc[m, "measurement"]) for m in rows.index})
        d = source.levels_frame(tag, cols, [])
        n = len(rows)
        ncol = 6
        fig, axes = plt.subplots(int(np.ceil(n / ncol)) or 1, ncol, figsize=(3 * ncol, 2.2 * max(1, np.ceil(n / ncol))),
                                 squeeze=False, layout="constrained")
        for ax, m in zip(axes.ravel(), rows.index):
            v = np.log1p(np.clip(d[reg.column(m, rows.loc[m, "measurement"])].dropna().to_numpy(float), 0, None))
            ax.hist(v[v > 0], bins=120, color="0.4")
            if np.isfinite(rows.loc[m, "level_log"]):
                ax.axvline(rows.loc[m, "level_log"], color="crimson", lw=1)
            ax.set_title(f"{m} ({rows.loc[m, 'method']}) {rows.loc[m, 'positive_fraction']:.1%}+", fontsize=7)
            ax.tick_params(labelsize=6)
        for ax in axes.ravel()[n:]:
            ax.axis("off")
        fig.suptitle(f"{tag}: typing-marker levels (red), log1p of the positivity reading", fontsize=9)
        path = gates.figure_path(task, f"levels/{tag}.png")
        path.parent.mkdir(parents=True, exist_ok=True)          # figure_path makes the gate folder, not subfolders
        fig.savefig(path, dpi=100)
        plt.close(fig)


def _class_figure(task: Task, cells: dict[str, pd.DataFrame]) -> None:
    import matplotlib.pyplot as plt

    frame = pd.DataFrame({tag: c["refined_class"].replace("", "(none)").value_counts(normalize=True)
                          for tag, c in cells.items()}).fillna(0.0)
    frame = frame.loc[frame.max(axis=1).sort_values(ascending=False).index]
    fig, ax = plt.subplots(figsize=(10, 5), layout="constrained")
    bottom = np.zeros(frame.shape[1])
    for cls, row in frame.iterrows():
        ax.bar(frame.columns, row.to_numpy(), bottom=bottom, label=cls)
        bottom += row.to_numpy()
    ax.set_ylabel("share of eligible cells")
    ax.legend(fontsize=7, ncol=2, bbox_to_anchor=(1.01, 1), loc="upper left")
    ax.tick_params(axis="x", rotation=45, labelsize=8)
    fig.savefig(gates.figure_path(task, "cell_classes.png"), dpi=120)
    plt.close(fig)


# --------------------------------------------------------------------------------- gate ----
def gate(task: Task, summary: dict) -> gates.Gate:
    per = summary.get("per_section", {})
    panels = summary.get("panels", {})
    required_stops = ({p: s["markers"]["required_stopped"] for p, s in panels.items()}
                      if panels else summary["markers"]["required_stopped"])
    checks = [
        gates.Check("every required marker has a passing level on every section",
                    required_stops, "none stopped within either panel", not summary["markers"]["required_stopped"]),
        gates.Check("cells typed", {t: v["cells"] for t, v in per.items()}, "every typed section",
                    bool(summary.get("typed")) and set(per) == set(summary["sections"])),
        gates.Check("cells named by a marker below its level", {t: v["named_below_level"] for t, v in per.items()},
                    "0", bool(per) and all(v["named_below_level"] == 0 for v in per.values())),
        gates.Check("forbidden broad or double-negative T final calls",
                    {t: v.get("forbidden_final_classes", 0) for t, v in per.items()}, "0",
                    bool(per) and all(v.get("forbidden_final_classes", 0) == 0 for v in per.values())),
        gates.Check("T-cell calls without CD3e and CD4-or-CD8 support",
                    {t: v.get("unsupported_t_cells") for t, v in per.items()}, "0",
                    bool(per) and all(v.get("unsupported_t_cells") == 0 for v in per.values())),
        gates.Check("final class by section", {t: v["refined_class"] for t, v in per.items()},
                    "reported", True, blocking=False),
        gates.Check("sections not typed (panel lacks the lineage roots)", summary["not_typed"], "listed", True,
                    blocking=False),
        gates.Check("markers whose level method stops (left out of typing)",
                    {p: s["markers"]["stop"] for p, s in panels.items()} if panels else summary["markers"]["stop"],
                    "reported; a stop on a required marker blocks", True, blocking=False),
        gates.Check("status by section", {t: v["status"] for t, v in per.items()}, "reported", True, blocking=False),
        gates.Check("broad class by section", {t: v["cell_class"] for t, v in per.items()}, "reported", True,
                    blocking=False),
        gates.Check("markers available and selected by panel",
                    {p: {"observed": s["observed_markers"], **s["markers"]} for p, s in panels.items()},
                    "reported separately; no marker availability pooled across panels", True, blocking=False),
        gates.Check("panel-specific rule adaptations",
                    {p: s["rule_adaptations"] for p, s in panels.items()}, "explicit unavailable clauses only", True,
                    blocking=False),
        gates.Check("rules unavailable by panel (marker not acquired)", summary.get("unavailable_rules", {}),
                    "unavailable, never negative", True, blocking=False),
        gates.Check("structure-only branches (no vessel object yet)", summary.get("structure_only", []),
                    "unavailable until the panel's vessel objects are accepted", True, blocking=False),
        gates.Check("cells outside every cell-derived region (start class Tissue)",
                    {t: v["outside_regions"] for t, v in per.items()}, "reported", True, blocking=False),
    ]
    figures = sorted(f"levels/{p.name}" for p in (task.staging / gates.GATE_DIR / "levels").glob("*.png"))
    if summary.get("typed"):
        figures.append("cell_classes.png")
    return gates.Gate(STAGE, "Spleen PC67 and PC33 Cell Typing: panel-specific per-image levels (positivity), named by the 6551 spleen "
                      "hierarchy; the first manual stop.", checks, figures)
