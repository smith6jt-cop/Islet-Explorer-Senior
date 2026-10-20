"""Pins and paths for the broad-typing bake-off lane.

Everything this lane reads is pinned by content here, once.  The run is the immutable
``30f6689e4d27d295f890``; the release pins the expression manifest's TOP-LEVEL ``content_id``
(the value ``phenospace.load_release`` compares), and records ``output.content_sha256`` beside
it.  Nothing here is hashed into ``run_id`` (``pipeline.py`` hashes ``phenocycler/`` and
``pyproject.toml`` only), and nothing in this lane writes ``data/runs/``.
"""
from __future__ import annotations

import json
import pathlib
import sys

HERE = pathlib.Path(__file__).resolve().parent
SCRIPTS = HERE.parent
PKG = SCRIPTS.parent                         # core/multiplex5000
REPO = PKG.parents[1]                        # repo root

for _sub in ("common", "celltype", "dct", "kintsugi", "xenium_celltype"):
    _p = str(SCRIPTS / _sub)
    if _p not in sys.path:
        sys.path.insert(0, _p)
if str(PKG) not in sys.path:
    sys.path.insert(0, str(PKG))

import panel  # noqa: E402  (scripts/celltype/panel.py -- the run pin and the panel authority)

RUN_ID = panel.RUN_ID
EXPRESSION = panel.EXPRESSION
EXPRESSION_TOTAL_ROWS = panel.EXPRESSION_TOTAL_ROWS
EXPRESSION_OUTPUT_SHA256 = panel.EXPRESSION_CONTENT_SHA256     # manifests/expression.json output.content_sha256
EXPRESSION_MANIFEST = PKG / "data" / "runs" / RUN_ID / "manifests" / "expression.json"
#: The value phenospace validates a release against (manifest top-level ``content_id``).
EXPRESSION_CONTENT_ID = "73cd7ee2084a937ce330241e3e6f3c037a946fa7cca709eb766c87630c704bf0"

DONORS = ("6436", "6442")                    # T1D, ND; both PhenoCycler batch 1
PIXEL_SIZE_UM = 0.5077663810243286

OUT = REPO / "data" / "development" / "broadtype"
UNIFORM = REPO / "data" / "development" / "celltype" / "uniform"        # Step-2 aligned values
USABILITY = REPO / "data" / "development" / "celltype" / "marker_donor_usability.parquet"
SHIFTS = REPO / "data" / "development" / "celltype" / "shifts.json"
REGISTRY = PKG / "multiplex5000" / "marker_registry.json"
RULES_37 = REPO / "phenotype_workflow.csv"
QUPATH_MANIFEST = PKG / "data" / "qupath_manifest.json"
ARCHIVE_MASKS = pathlib.Path("/data/codex_phenocycler/archive_2026-07-25/masks")
GEOMETRY_QC = PKG / "data" / "runs" / RUN_ID / "geometry_qc"
CELLS = PKG / "data" / "runs" / RUN_ID / "cells"

#: Registry ``kind`` values.  Identity inputs may use ``type`` only (the tiering rule).
KIND_TYPE, KIND_STATE, KIND_EXCLUDED, KIND_STRUCTURAL = "type", "state", "excluded", "structural"


def assert_pins() -> dict:
    """Fail loud if the run or the manifest moved.  Returns the manifest summary."""
    m = json.loads(EXPRESSION_MANIFEST.read_text())
    got = str(m.get("content_id", ""))
    if got != EXPRESSION_CONTENT_ID:
        raise RuntimeError(f"expression manifest content_id {got!r} != pinned {EXPRESSION_CONTENT_ID!r}")
    out_sha = str(m.get("output", {}).get("content_sha256", ""))
    if out_sha != EXPRESSION_OUTPUT_SHA256:
        raise RuntimeError(f"expression output sha {out_sha!r} != celltype pin {EXPRESSION_OUTPUT_SHA256!r}")
    latest = (PKG / "data" / "runs" / "LATEST").read_text().strip()
    if latest != RUN_ID:
        raise RuntimeError(f"runs/LATEST {latest!r} != RUN_ID {RUN_ID!r}; repoint every pin together")
    return {"run_id": RUN_ID, "content_id": got, "output_content_sha256": out_sha}


def registry() -> dict:
    return json.loads(REGISTRY.read_text())


def markers_by_kind() -> dict:
    """kind -> sorted marker names, straight from the registry."""
    out: dict = {}
    for m in registry()["markers"]:
        out.setdefault(m["kind"], []).append(m["name"])
    return {k: sorted(v) for k, v in out.items()}


def manifest_image(donor: str) -> dict:
    rows = [r for r in json.loads(QUPATH_MANIFEST.read_text())["images"] if str(r["donor_id"]) == str(donor)]
    if len(rows) != 1:
        raise KeyError(f"donor {donor}: {len(rows)} manifest rows")
    return rows[0]


def geojson_path(donor: str) -> pathlib.Path:
    return pathlib.Path(manifest_image(donor)["cell_geojson"]["path"])


def qptiff_path(donor: str) -> pathlib.Path:
    return pathlib.Path(manifest_image(donor)["qptiff"]["path"])


def expression_parquet(donor: str) -> pathlib.Path:
    return EXPRESSION / f"donor_id={donor}" / "data.parquet"


def uniform_parquet(donor: str) -> pathlib.Path:
    return UNIFORM / f"donor_id={donor}" / "data.parquet"


def archive_mask(donor: str) -> pathlib.Path:
    return ARCHIVE_MASKS / f"{donor}.tif"


def git_head() -> dict:
    import subprocess
    try:
        sha = subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=REPO, text=True).strip()
        dirty = bool(subprocess.check_output(["git", "status", "--porcelain"], cwd=REPO, text=True).strip())
    except Exception:  # noqa: BLE001
        sha, dirty = "", True
    return {"commit": sha, "dirty": dirty}
