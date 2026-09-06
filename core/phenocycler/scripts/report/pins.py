"""Pins and paths for the two-donor report lane (plan Phase 4).

Everything the report reads is named here once, by content where content exists.  The lane is
development only: it writes ``data/development/report/``, never ``data/runs/`` and never a stage
manifest, and nothing under ``scripts/**`` is hashed into ``run_id``.

Two rules shape this file.  First, the run and the Step-2 alignment are pinned by the same
constants the celltype and broadtype lanes already hold (``scripts/celltype/panel.py``,
``scripts/broadtype/pins.py``) -- imported, not copied, so a repoint of the run happens in one
place and every lane moves together.  Second, the Xenium inputs are resolved from
``data/integration/xenium_paths.csv`` (the bundle path authority) and the processed base a
build stamped with its own ``uns['prep']``; nothing is looked up by donor-name parsing of a
directory listing.

The figures this lane makes need NO identity release: P01-P05 and X01-X03 are QC.  The
registry latch ``NO_VERIFIED_TYPING`` stays down and every title says so (``DEV_LINE``).
"""
from __future__ import annotations

import hashlib
import importlib.util
import json
import pathlib
import subprocess
import sys
import time

HERE = pathlib.Path(__file__).resolve().parent
SCRIPTS = HERE.parent
PKG = SCRIPTS.parent                         # core/phenocycler
REPO = PKG.parents[1]                        # repo root

for _sub in ("common", "celltype", "restain", "xenium_celltype"):
    _p = str(SCRIPTS / _sub)
    if _p not in sys.path:
        sys.path.insert(0, _p)
if str(PKG) not in sys.path:
    sys.path.insert(0, str(PKG))

import panel  # noqa: E402  (scripts/celltype/panel.py -- the run pin and the panel authority)


def _load_broadtype_pins():
    """The broadtype spine, imported by file so its module name ``pins`` cannot shadow ours."""
    spec = importlib.util.spec_from_file_location("broadtype_pins", SCRIPTS / "broadtype" / "pins.py")
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


_bt = _load_broadtype_pins()

# --------------------------------------------------------------------------- the run
RUN_ID = panel.RUN_ID
RUN_ROOT = PKG / "data" / "runs" / RUN_ID
RUN_JSON = RUN_ROOT / "run.json"
MANIFESTS = {s: RUN_ROOT / "manifests" / f"{s}.json" for s in ("ingest", "geometry", "redsea", "expression")}
#: run.json's own content id (the run manifest), recorded so a rebuilt run.json is caught.
RUN_CONTENT_ID = "43ffce45dbcbd8b579b55df0ba61e0fbd444a31518d5741d09a93d5e323e9408"
#: Expression manifest top-level ``content_id`` (what phenospace validates a release against).
EXPRESSION_CONTENT_ID = _bt.EXPRESSION_CONTENT_ID
#: Expression manifest ``output.content_sha256`` (the celltype lane's pin).
EXPRESSION_OUTPUT_SHA256 = panel.EXPRESSION_CONTENT_SHA256

CELLS = RUN_ROOT / "cells"
REDSEA = RUN_ROOT / "redsea"
EXPRESSION = RUN_ROOT / "expression"
GEOMETRY_QC = RUN_ROOT / "geometry_qc"
INGEST_REJECTS = RUN_ROOT / "ingest_rejects.parquet"
QUPATH_MANIFEST = PKG / "data" / "qupath_manifest.json"
REGISTRY = PKG / "phenocycler" / "marker_registry.json"

# --------------------------------------------------------------------------- Step 2 and the frame
UNIFORM = _bt.UNIFORM                        # data/development/celltype/uniform/donor_id=<d>/data.parquet
SHIFTS = _bt.SHIFTS                          # data/development/celltype/shifts.json
FRAME = _bt.OUT / "frame"                    # data/development/broadtype/frame/
UNIFORM_NUM_BINS = 1024                      # uniform_adapter.NUM_BINS; asserted against shifts.json

# --------------------------------------------------------------------------- donors
DONORS = tuple(_bt.DONORS)                   # ("6436", "6442"): T1D, ND; both PhenoCycler batch 1
PIXEL_SIZE_UM = _bt.PIXEL_SIZE_UM
#: Xenium slide serial per donor.  Asserted against xenium_paths.csv and the base's uns['prep'].
XENIUM_SERIAL = {"6436": "0059907", "6442": "0059947"}

# --------------------------------------------------------------------------- Xenium
XENIUM_PATHS_CSV = PKG / "data" / "integration" / "xenium_paths.csv"
XENIUM_ROOT = pathlib.Path("/data/xenium")
XENIUM_PROCESSED = XENIUM_ROOT / "processed"
#: cell_feature_matrix.h5 must carry exactly these feature types with exactly these counts.
XENIUM_FEATURE_TYPES = {
    "Gene Expression": 5101,
    "Deprecated Codeword": 3291,
    "Negative Control Codeword": 609,
    "Unassigned Codeword": 582,
    "Negative Control Probe": 40,
    "Genomic Control": 21,
}
#: cells.parquet per-cell count columns and the feature type each one totals.
XENIUM_COUNT_COLUMNS = {
    "transcript_counts": "Gene Expression",
    "deprecated_codeword_counts": "Deprecated Codeword",
    "control_codeword_counts": "Negative Control Codeword",
    "unassigned_codeword_counts": "Unassigned Codeword",
    "control_probe_counts": "Negative Control Probe",
    "genomic_control_counts": "Genomic Control",
}

# --------------------------------------------------------------------------- output and text
OUT = REPO / "data" / "development" / "report"
PAIR_DIR = "_pair"
#: The restain Explorer bundle is under a DO_NOT_SEND hold; this lane never reads from it.
FORBIDDEN_SOURCE = REPO / "data" / "development" / "restain" / "explorer_bundle"
ISLET_YAML = REPO / "apps" / "islet" / "config" / "islet_default.yml"

DEV_LINE = ("development evidence — registry latch NO_VERIFIED_TYPING down; "
            "this method at this configuration on these donors")
#: The two markers of the mutable diagnostic roster (raw CD3e / Ker8_18), drawn in P05.
REDSEA_PAIR = ("CD3e", "Ker8_18")
NUCLEUS_MEASURED = tuple(panel.NUCLEUS_MEASURED)
METHOD_VERSION = "report-qc-figures-2026-09-05"


# --------------------------------------------------------------------------- lookups
def run_json() -> dict:
    return json.loads(RUN_JSON.read_text())


def manifest(stage: str) -> dict:
    return json.loads(MANIFESTS[stage].read_text())


def geometry_thresholds() -> dict:
    """``run.json`` -> ``config.canonical`` (a JSON string) -> ``geometry_qc``."""
    canonical = json.loads(run_json()["config"]["canonical"])
    return dict(canonical["geometry_qc"])


def assert_pins() -> dict:
    """Fail loud if the run, its manifests, or the Step-2 alignment moved."""
    r = run_json()
    if str(r.get("content_id")) != RUN_CONTENT_ID:
        raise RuntimeError(f"run.json content_id {r.get('content_id')!r} != pinned {RUN_CONTENT_ID!r}")
    m = manifest("expression")
    if str(m.get("content_id")) != EXPRESSION_CONTENT_ID:
        raise RuntimeError(f"expression manifest content_id {m.get('content_id')!r} != pinned {EXPRESSION_CONTENT_ID!r}")
    out_sha = str(m.get("output", {}).get("content_sha256", ""))
    if out_sha != EXPRESSION_OUTPUT_SHA256:
        raise RuntimeError(f"expression output sha {out_sha!r} != pinned {EXPRESSION_OUTPUT_SHA256!r}")
    latest = (RUN_ROOT.parent / "LATEST").read_text().strip()
    if latest != RUN_ID:
        raise RuntimeError(f"runs/LATEST {latest!r} != RUN_ID {RUN_ID!r}; repoint every pin together")
    sh = json.loads(SHIFTS.read_text())
    if int(sh.get("num_bins", -1)) != UNIFORM_NUM_BINS:
        raise RuntimeError(f"shifts.json num_bins {sh.get('num_bins')} != {UNIFORM_NUM_BINS}")
    return {"run_id": RUN_ID, "run_content_id": RUN_CONTENT_ID, "expression_content_id": EXPRESSION_CONTENT_ID,
            "expression_output_sha256": out_sha, "uniform_commit": sh.get("uniform_commit", "")}


def stage_content_ids() -> dict:
    """Top-level ``content_id`` of every stage manifest, for the footer."""
    return {s: str(json.loads(p.read_text()).get("content_id", "")) for s, p in MANIFESTS.items() if p.exists()}


def cohort_row(donor: str) -> dict:
    coh = panel.cohort()
    row = coh.loc[coh["donor_id"].astype(str) == str(donor)]
    if len(row) != 1:
        raise KeyError(f"donor {donor}: {len(row)} cohort rows")
    return row.iloc[0].to_dict()


def donor_status(donor: str) -> str:
    """The stored status token (``ND`` / ``AAB`` / ``T1D``); display goes through style.status_display."""
    return str(cohort_row(donor)["status"])


def manifest_image(donor: str) -> dict:
    return _bt.manifest_image(donor)


def qptiff_path(donor: str) -> pathlib.Path:
    return _bt.qptiff_path(donor)


def geojson_path(donor: str) -> pathlib.Path:
    return _bt.geojson_path(donor)


def channel_index(donor: str, marker: str) -> int:
    """qptiff plane index of a registry marker, from the manifest's channel map (never by name parsing)."""
    for c in manifest_image(donor)["channel_map"]:
        if c["marker"] == marker:
            return int(c["channel_index"])
    raise KeyError(f"donor {donor}: marker {marker!r} is not in the qptiff channel map")


def acquired_markers(donor: str) -> list:
    """Every marker the slide carries, in channel order (59 on both donors, DAPI and IAPP included)."""
    return [c["marker"] for c in sorted(manifest_image(donor)["channel_map"], key=lambda c: c["channel_index"])]


def expression_parquet(donor: str) -> pathlib.Path:
    return EXPRESSION / f"donor_id={donor}" / "data.parquet"


def geometry_parquet(donor: str) -> pathlib.Path:
    return GEOMETRY_QC / f"donor_id={donor}" / "data.parquet"


def cells_parquet(donor: str) -> pathlib.Path:
    return CELLS / f"donor_id={donor}" / "data_0.parquet"


def redsea_parquet(donor: str) -> pathlib.Path:
    return REDSEA / f"donor_id={donor}" / "data_0.parquet"


def uniform_parquet(donor: str) -> pathlib.Path:
    return UNIFORM / f"donor_id={donor}" / "data.parquet"


def xenium_bundle(donor: str, roi: str = "panc") -> pathlib.Path:
    import pandas as pd
    xen = pd.read_csv(XENIUM_PATHS_CSV, dtype=str)
    hit = xen[(xen["donor_id"] == str(donor)) & (xen["roi"] == roi)]
    hit = hit[hit["source_path"].notna() & (hit["source_path"] != "NA")]
    if len(hit) != 1:
        raise KeyError(f"donor {donor} roi={roi}: {len(hit)} bundle rows in {XENIUM_PATHS_CSV}")
    src = str(hit.iloc[0]["source_path"]).strip().rstrip("/")
    serial = XENIUM_SERIAL[str(donor)]
    if f"__{serial}__" not in src:
        raise RuntimeError(f"donor {donor}: bundle path {src!r} does not carry pinned serial {serial}")
    return XENIUM_ROOT / src


def xenium_base(donor: str, roi: str = "panc") -> pathlib.Path:
    serial = XENIUM_SERIAL[str(donor)]
    return XENIUM_PROCESSED / serial / roi / f"{serial}_phenotyped.h5ad"


def metrics_summary_path(donor: str) -> pathlib.Path | None:
    """``metrics_summary.csv`` is being transferred (plan 0.4); None means 'not transferred'."""
    p = xenium_bundle(donor) / "metrics_summary.csv"
    return p if p.exists() else None


def git_head() -> dict:
    try:
        sha = subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=REPO, text=True).strip()
        dirty = bool(subprocess.check_output(["git", "status", "--porcelain"], cwd=REPO, text=True).strip())
    except Exception:  # noqa: BLE001
        sha, dirty = "", True
    return {"commit": sha, "dirty": dirty}


def sha256_file(path: pathlib.Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 22), b""):
            h.update(chunk)
    return h.hexdigest()


_BASE_SHA: dict = {}


def base_sha256(donor: str) -> str:
    """Full sha256 of the Xenium base (1.2-1.6 GB; ~5 s), cached for the process."""
    p = xenium_base(donor)
    key = str(p)
    if key not in _BASE_SHA:
        _BASE_SHA[key] = sha256_file(p)
    return _BASE_SHA[key]


def provenance(donors=None, *, xenium: bool = True) -> dict:
    """Everything the footer prints: run id, content ids, base sha, code commit, UTC."""
    donors = list(donors or DONORS)
    prov = {
        "method_version": METHOD_VERSION,
        "run_id": RUN_ID,
        "run_content_id": RUN_CONTENT_ID,
        "expression_content_id": EXPRESSION_CONTENT_ID,
        "expression_output_sha256": EXPRESSION_OUTPUT_SHA256,
        "stage_content_ids": stage_content_ids(),
        "uniform_commit": json.loads(SHIFTS.read_text()).get("uniform_commit", "") if SHIFTS.exists() else "",
        "identity_release": None,          # QC figures need none; the latch stays down
        "code": git_head(),
        "written_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
        "donors": {},
    }
    for d in donors:
        row = {"status": donor_status(d)}
        if xenium:
            row.update(serial=XENIUM_SERIAL[str(d)], base=str(xenium_base(d)),
                       base_sha256=base_sha256(d) if xenium_base(d).exists() else "")
        prov["donors"][str(d)] = row
    return prov
