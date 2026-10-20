"""Pins and paths for the Xenium typing lane (plan 2026-09-05, Phase 2).

Everything this lane reads is pinned by content here, once.  The inputs are the per-donor
Xenium bases (``/data/xenium/processed/<serial>/panc/<serial>_phenotyped.h5ad``, built by
``scripts/xenium_celltype/build_phenotyped_base.py`` and stamped by
``analysis/xenium/scripts/finalize_panel.py``), the raw bundles named in
``data/integration/xenium_paths.csv`` (Xenium's own ``cells.parquet`` is the complete cell
set every release must cover), and the panel snapshot generated from ``core/panels`` -- the
ONE taxonomy authority (root ``CLAUDE.md``) -- whose ``SOURCE.txt`` tree id is pinned below.

Nothing here is hashed into a PhenoCycler ``run_id`` (``pipeline.py`` hashes ``phenocycler/``
and ``pyproject.toml`` only), nothing writes ``data/runs/``, and nothing writes a stage
manifest.  Output is development evidence under ``data/development/xenium_typing/``.

The donor is the join unit and ``cell_id`` (Xenium's own barcode) is the only cell key.  The
Xenium *run* (``MM0xx``, the third path segment of every bundle path) is the technical unit:
one instrument run carries two slides, and a slide serial encodes no donor and no region.
"""
from __future__ import annotations

import hashlib
import json
import pathlib
import re
import subprocess
import sys
import time

HERE = pathlib.Path(__file__).resolve().parent            # core/multiplex5000/scripts/xenium_typing
SCRIPTS = HERE.parent                                     # core/multiplex5000/scripts
PKG = SCRIPTS.parent                                      # core/multiplex5000
REPO = PKG.parents[1]                                     # repo root
XENIUM_ANALYSIS = REPO / "analysis" / "xenium"            # the Xenium lane (decisions #1-#17)

# Sibling lanes this one reuses, by directory (the repo's scripts/ convention).  ``PKG`` for
# the ``phenocycler`` package; ``analysis/xenium`` at the FRONT so ``scripts.tissue_markers``
# and ``scripts.workflow_paths`` resolve to that lane's namespace package.
for _sub in ("xenium_celltype", "xenium_celltype/hpg", "celltype", "common"):
    _p = str(SCRIPTS / _sub)
    if _p not in sys.path:
        sys.path.insert(0, _p)
if str(PKG) not in sys.path:
    sys.path.insert(0, str(PKG))
if str(XENIUM_ANALYSIS) not in sys.path:
    sys.path.insert(0, str(XENIUM_ANALYSIS))

#: The eight Xenium broad lineages -- the label space of every arm.  Order fixed for tables.
LINEAGES = ("Endocrine", "Exocrine", "Immune_B", "Immune_T", "Myeloid", "Neural", "Stromal", "Vascular")

#: Arm names.  A configuration change after a result is seen is a NEW arm name (protocol §0).
ARM_SCANVI = "scanvi-per-donor"           # X1: the incumbent per-donor identity-core scANVI sidecars
ARM_JOINT = "joint-run-aware"             # X2: one scVI over a hash subsample of every donor, batch = run
ARM_BASELINE = "xenium-analyzer-baseline"  # X3: the instrument's onboard graphclust / kmeans clusters
ARMS = (ARM_SCANVI, ARM_JOINT, ARM_BASELINE)

#: The interpreter for the scVI arm.  Every other module runs in the repo ``.venv``.
SCVI_ENV_PYTHON = pathlib.Path("/home/smith6jt/miniconda3/envs/scvi-env/bin/python")

# --------------------------------------------------------------------------- paths
OUT = REPO / "data" / "development" / "xenium_typing"
COHORT_MANIFEST = PKG / "data" / "qupath_manifest.json"
XENIUM_PATHS = PKG / "data" / "integration" / "xenium_paths.csv"
CONFIG_INI = PKG / "config.ini"                            # the ONLY config.ini (root CLAUDE.md)
PROCESSED_ROOT = pathlib.Path("/data/xenium/processed")
XENIUM_ROOT = pathlib.Path("/data/xenium")
SUBPANEL_DIR = REPO / "data" / "development" / "xenium_panels" / "subpanels"
PANEL_SOURCE_TXT = SUBPANEL_DIR / "SOURCE.txt"
RESTAIN_PROTEIN_RNA = REPO / "data" / "development" / "restain" / "protein_rna"
HIERARCHICAL_PHENOTYPE_PY = XENIUM_ANALYSIS / "scripts" / "hierarchical_phenotype.py"
STATE_PANELS_PY = XENIUM_ANALYSIS / "scripts" / "state_panels.py"

#: Panel snapshot pins (analysis/xenium standing note #16/#17): the git tree id of
#: ``core/panels/data`` the snapshot was generated from, and the upstream commit it equals.
PANEL_TREE_ID = "0348d149c4e4386488e06807f367cf8122c293c7"
PANEL_UPSTREAM_PIN = "XeniumPanelExplorer@076d1192a25946ac3f9cf941c98202d1f224579a"

#: The bundle directory name carries the slide serial (build_phenotyped_base._SERIAL).
_SERIAL = re.compile(r"output-XETG\d+__(\d+)__")
#: The instrument run is the third path segment: ``HuPrime5kPanTissue/<batch>/MM0xx_.../``.
_RUN = re.compile(r"^(MM\d{3})")

ROI = "panc"


class PinError(RuntimeError):
    """A pinned input moved, is absent, or disagrees with its pin -- never silently continued."""


# --------------------------------------------------------------------------- donors
def donors() -> tuple[str, ...]:
    """The acquired cohort, normalized for the Xenium donor -> bundle join.

    The live QuPath manifest includes acquisitions that the build specification can
    predate. Its zero-padded donor IDs (``0112``) refer to the same donor as ``112`` in
    ``xenium_paths.csv``. Central eligibility is checked after normalization.
    """
    from phenocycler.cohort import ensure_eligible_donors
    if not COHORT_MANIFEST.is_file():
        raise PinError(f"{COHORT_MANIFEST} is absent; the acquired cohort is unavailable")
    manifest = json.loads(COHORT_MANIFEST.read_text())
    cohort = tuple(dict.fromkeys(str(d).strip().lstrip("0")
                                 for d in manifest["expected_donors"]))
    if not cohort or any(not d for d in cohort):
        raise PinError(f"{COHORT_MANIFEST} has no valid acquired cohort")
    return ensure_eligible_donors(cohort, context="xenium typing lane")


# --------------------------------------------------------------------------- bundle paths
def paths_table(csv: pathlib.Path = XENIUM_PATHS, roi: str = ROI):
    """Cohort bundle rows for one ROI with a recorded path (``NA`` rows dropped).

    ``PhenoXeno`` rows describe separate same-slide acquisitions, not the MM0xx
    cohort bundles this lane types. A donor may have both in ``xenium_paths.csv``;
    those rows must never replace the cohort serial or introduce a second match.
    Whitespace and trailing slashes are stripped: the hand-maintained file carries both.
    """
    import pandas as pd
    if not pathlib.Path(csv).exists():
        raise PinError(f"{csv} is absent; it is the donor -> bundle manifest this lane keys on")
    t = pd.read_csv(csv, dtype=str)
    t = t[(t["roi"] == roi) & (t["batch"].str.strip() != "PhenoXeno")
          & t["source_path"].notna() & (t["source_path"].str.strip() != "NA")].copy()
    t["source_path"] = t["source_path"].str.strip().str.rstrip("/")
    t["donor_id"] = t["donor_id"].str.strip()
    duplicates = t.loc[t["donor_id"].duplicated(), "donor_id"].tolist()
    if duplicates:
        raise PinError(f"donor(s) {duplicates} have duplicate {roi} cohort bundle rows in {csv}")
    return t.reset_index(drop=True)


def serial_of(source_path: str) -> str:
    m = _SERIAL.search(str(source_path))
    if not m:
        raise PinError(f"no XETG serial in bundle path {source_path!r}")
    return m.group(1)


def xenium_run_of(source_path: str) -> str:
    """``MM0xx`` from the third path segment.  This is the scVI batch key (see joint_model)."""
    parts = [p for p in str(source_path).strip("/").split("/") if p]
    if len(parts) < 3:
        raise PinError(f"bundle path {source_path!r} has no run segment")
    m = _RUN.match(parts[2])
    if not m:
        raise PinError(f"third segment {parts[2]!r} of {source_path!r} is not an MM0xx run")
    return m.group(1)


def serial_map(csv: pathlib.Path = XENIUM_PATHS, roi: str = ROI) -> dict:
    t = paths_table(csv, roi)
    out = {}
    for d, p in zip(t["donor_id"], t["source_path"]):
        out[d] = serial_of(p)
    return out


def run_map(csv: pathlib.Path = XENIUM_PATHS, roi: str = ROI) -> dict:
    t = paths_table(csv, roi)
    return {d: xenium_run_of(p) for d, p in zip(t["donor_id"], t["source_path"])}


def source_path(donor: str, csv: pathlib.Path = XENIUM_PATHS, roi: str = ROI) -> str:
    t = paths_table(csv, roi)
    hit = t.loc[t["donor_id"] == str(donor), "source_path"]
    if len(hit) != 1:
        raise PinError(f"donor {donor} roi={roi}: {len(hit)} bundle rows in {csv}")
    return str(hit.iloc[0])


def serial(donor: str) -> str:
    return serial_of(source_path(donor))


def xenium_run(donor: str) -> str:
    return xenium_run_of(source_path(donor))


def bundle_dir(donor: str) -> pathlib.Path:
    return XENIUM_ROOT / source_path(donor)


def analysis_dir(donor: str) -> pathlib.Path:
    """The Xenium Analyzer output the Globus transfer is bringing in (may be absent)."""
    return bundle_dir(donor) / "analysis"


def cells_parquet(donor: str) -> pathlib.Path:
    return bundle_dir(donor) / "cells.parquet"


def base_path(donor: str) -> pathlib.Path:
    s = serial(donor)
    return PROCESSED_ROOT / s / ROI / f"{s}_phenotyped.h5ad"


def sidecar_path(donor: str) -> pathlib.Path:
    s = serial(donor)
    return PROCESSED_ROOT / s / ROI / f"{s}_hierarchical.h5ad"


def if_values_path(donor: str, marker: str, compartment: str = "cell") -> pathlib.Path:
    """Restain lane per-cell IF values (``scripts/restain/protein_rna.py``), aligned to the
    bundle's ``cells.parquet`` order.  Protein evidence; no positivity is implied by them."""
    return RESTAIN_PROTEIN_RNA / f"{donor}_{serial(donor)}_{compartment}_{marker}.npy"


# --------------------------------------------------------------------------- lane output
def arm_dir(arm: str, donor: str) -> pathlib.Path:
    return OUT / "arms" / arm / f"donor_id={donor}"


def gate_dir(arm: str) -> pathlib.Path:
    return OUT / "gate" / arm


def gate_json(arm: str, donor: str) -> pathlib.Path:
    return gate_dir(arm) / f"donor_id={donor}.json"


# --------------------------------------------------------------------------- panel snapshot
def panel_snapshot(source_txt: pathlib.Path = PANEL_SOURCE_TXT) -> dict:
    """Parse ``SOURCE.txt`` (``key: value`` lines) written by ``fetch_subpanels.py``."""
    if not pathlib.Path(source_txt).exists():
        raise PinError(f"{source_txt} is absent; regenerate the snapshot with "
                       f"analysis/xenium/scripts/fetch_subpanels.py")
    out = {}
    for line in pathlib.Path(source_txt).read_text().splitlines():
        if ":" in line and not line.startswith(" "):
            k, v = line.split(":", 1)
            out[k.strip()] = v.strip()
    return out


def identity_cores() -> dict:
    """lineage -> sorted identity-core genes, from the panel taxonomy the integration layer
    already validates (``vocab.load_panel_taxonomy``; ``panel_roles.csv`` broad_mask rows).

    These are IDENTITY genes only.  State subpanels never enter here (decision #17); the
    bases' ``identity_panel & is_state_only`` is asserted 0 by the arms that read them."""
    from phenocycler.config import load_config
    from phenocycler.integration.vocab import load_panel_taxonomy
    tax = load_panel_taxonomy(load_config(CONFIG_INI), "pancreas")
    got = set(tax.panel_to_lineage.values())
    if got != set(LINEAGES):
        raise PinError(f"panel taxonomy lineages {sorted(got)} != pinned {sorted(LINEAGES)}")
    return {L: list(tax.lineage_genes(L)) for L in LINEAGES}


def assert_pins(*, check_git_tree: bool = True) -> dict:
    """Fail loud if the snapshot stamp or the taxonomy tree moved.  Returns the stamp."""
    snap = panel_snapshot()
    if snap.get("tree_id") != PANEL_TREE_ID:
        raise PinError(f"panel snapshot tree_id {snap.get('tree_id')!r} != pinned {PANEL_TREE_ID!r}")
    if snap.get("upstream_pin") != PANEL_UPSTREAM_PIN:
        raise PinError(f"panel snapshot upstream_pin {snap.get('upstream_pin')!r} != pinned {PANEL_UPSTREAM_PIN!r}")
    if check_git_tree:
        try:
            live = subprocess.check_output(["git", "rev-parse", "HEAD:core/panels/data"], cwd=REPO,
                                           text=True, stderr=subprocess.DEVNULL).strip()
        except Exception:  # noqa: BLE001 -- no git: the stamp alone is the pin
            live = ""
        if live and live != PANEL_TREE_ID:
            raise PinError(f"live core/panels/data tree {live} != pinned {PANEL_TREE_ID}; "
                           f"regenerate the snapshot and repoint every pin together")
    if not XENIUM_PATHS.exists():
        raise PinError(f"{XENIUM_PATHS} is absent")
    return snap


# --------------------------------------------------------------------------- provenance
def sha256_file(path: pathlib.Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 22), b""):
            h.update(chunk)
    return h.hexdigest()


def git_head() -> dict:
    try:
        sha = subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=REPO, text=True).strip()
        dirty = bool(subprocess.check_output(["git", "status", "--porcelain"], cwd=REPO, text=True).strip())
    except Exception:  # noqa: BLE001
        sha, dirty = "", True
    return {"commit": sha, "dirty": dirty}


def utc_now() -> str:
    return time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())
