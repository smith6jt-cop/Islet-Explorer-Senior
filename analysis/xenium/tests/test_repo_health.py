"""Repo-health tests for the absorbed Xenium lane.

Cheap: no heavy imports, no real data. They catch the things that are annoying to debug
later — a notebook that will not open, a shell script SLURM will reject, a `utils` module
that does not import, and a panel snapshot that has drifted from the taxonomy it was
generated from.

**What this file no longer asserts, and why.** Three of its original checks tested things
that do not exist in the monorepo:

* `environment.yml` — dropped on absorption. It named the `xenium_analysis` conda env, which
  is **broken on this machine**: `import scanpy` and `import scvi` both raise
  `ImportError: cannot import name '_spropack' from 'scipy.sparse.linalg._propack'`. The lane
  runs in `scvi-env`. A test asserting a broken env is worse than no test.
* `config.ini` — dropped on absorption. Its `[phenotyping]` block mirrored *retired* marker
  panels (`Epsilon = GHRL`, the single-gene argmax that mislabelled ~150K cells and that
  decision #12 exists to prevent; `Acinar = AMY1A`). The root `CLAUDE.md` says there is
  exactly one `config.ini` and it is `core/phenocycler/config.ini`.
* the SLURM `xenium_analysis` env assertion — the batch scripts moved to `hpc/` and are
  HiPerGator-only. Their shell hygiene is still checked; which env they name is a cluster
  question this repo cannot answer.
"""

from __future__ import annotations

import json
import subprocess
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]          # analysis/xenium
REPO = Path(__file__).resolve().parents[3]          # monorepo root

SLURM_DIRS = [ROOT / "hpc" / "slurm", ROOT / "hpc" / "benchmark_slurm"]
NOTEBOOKS_DIR = ROOT / "notebooks"
SNAPSHOT = REPO / "data" / "development" / "xenium_panels" / "subpanels"


# --------------------------------------------------------------------------- shell hygiene
SLURM_SCRIPTS = sorted(p for d in SLURM_DIRS if d.exists() for p in d.glob("*.sh"))
BATCH_JOBS = [s for s in SLURM_SCRIPTS if "#SBATCH" in s.read_text(errors="replace")]


@pytest.mark.parametrize("script", SLURM_SCRIPTS, ids=lambda p: p.name)
class TestShellHygiene:
    def test_shebang(self, script: Path):
        first = script.read_bytes().splitlines()[0]
        assert first.startswith(b"#!") and b"bash" in first, f"{script.name} shebang"

    def test_lf_line_endings(self, script: Path):
        """CRLF is silently rejected by SLURM."""
        assert b"\r\n" not in script.read_bytes(), f"{script.name} has CRLF"


@pytest.mark.parametrize("script", BATCH_JOBS, ids=lambda p: p.name)
class TestSlurmBatchJobs:
    def test_has_job_name(self, script: Path):
        assert "#SBATCH --job-name=" in script.read_text()

    def test_has_resource_directives(self, script: Path):
        text = script.read_text()
        for d in ("--time=", "--mem", "--cpus-per-task"):
            assert d in text, f"{script.name} missing {d}"

    def test_outputs_to_logs_dir(self, script: Path):
        assert "--output=logs/" in script.read_text()


def test_the_hpc_lane_is_present_and_is_the_only_place_slurm_lives():
    assert SLURM_SCRIPTS, "expected SLURM scripts under hpc/"
    assert BATCH_JOBS, "expected at least one #SBATCH script"
    stray = sorted((ROOT / "scripts").rglob("*.sh"))
    stray_batch = [p for p in stray if "#SBATCH" in p.read_text(errors="replace")]
    assert not stray_batch, (
        f"SLURM jobs under scripts/: {[str(p) for p in stray_batch]}. They belong in hpc/, "
        "which is labelled HiPerGator-only — scripts/ is what runs here."
    )


# --------------------------------------------------------------------------- notebooks
NOTEBOOKS = sorted(NOTEBOOKS_DIR.glob("*.ipynb")) if NOTEBOOKS_DIR.exists() else []


@pytest.mark.parametrize("nb_path", NOTEBOOKS, ids=lambda p: p.name)
class TestNotebookValidity:
    def test_valid_json(self, nb_path: Path):
        nb = json.loads(nb_path.read_text())
        assert "cells" in nb and "metadata" in nb

    def test_has_python_kernelspec(self, nb_path: Path):
        ks = json.loads(nb_path.read_text()).get("metadata", {}).get("kernelspec", {})
        assert ks and "python" in ks.get("language", "").lower()

    def test_nbformat_version(self, nb_path: Path):
        assert json.loads(nb_path.read_text()).get("nbformat") == 4


def test_the_pipeline_notebooks_are_present_and_06_is_not():
    names = {p.name for p in NOTEBOOKS}
    required = {"01_preprocessing_v2.ipynb", "02_phenotyping.ipynb",
                "03_spatial_analysis.ipynb", "04_group_comparisons.ipynb",
                "05_tissue_comparisons.ipynb"}
    assert not required - names, f"missing: {sorted(required - names)}"
    # 06_xenium_phenocycler_integration.ipynb was an UNRUN skeleton -- placeholder
    # SAMPLE_NAME, a segmentation mask stamped from discs at the cell centroids, and one
    # global phase_cross_correlation shift. It is superseded by
    # core/phenocycler/phenocycler/integration/register.py (two registration tracks, a
    # fitted px->um scale, vertex-wise polygon warping). Dropped on absorption so there is
    # not a second, worse, runnable answer to the same question.
    assert "06_xenium_phenocycler_integration.ipynb" not in names


# --------------------------------------------------------------------------- no stale paths
def test_no_hipergator_paths_survive_in_executable_code():
    """The absorption's acceptance criterion, asserted rather than remembered.

    Documented HiPerGator *defaults* (`${VAR:-/blue/...}`) are allowed and are how the hpc/
    lane still works there; a bare hardcoded path is not. `docs/` is excluded: those are the
    historical handoffs, and rewriting them would falsify the record.
    """
    out = subprocess.run(
        ["grep", "-rn", "blue/maigan", "--include=*.py", "--include=*.ipynb", "--include=*.R",
         str(ROOT)],
        capture_output=True, text=True).stdout.splitlines()
    me = str(Path(__file__).resolve())
    bad = [ln for ln in out
           if not ln.startswith(me)                       # this file names the pattern
           and ":-/blue" not in ln                        # documented HPG default
           and 'Sys.getenv("R_LIBS_USER"' not in ln]     # ditto, R side
    assert not bad, "hardcoded HiPerGator paths:\n" + "\n".join(bad)


# --------------------------------------------------------------------------- panel snapshot
def _snapshot_tree_id() -> str | None:
    src = SNAPSHOT / "SOURCE.txt"
    if not src.exists():
        return None
    for line in src.read_text().splitlines():
        if line.startswith("tree_id:"):
            return line.split(":", 1)[1].strip()
    return None


@pytest.mark.skipif(not (SNAPSHOT / "SOURCE.txt").exists(),
                    reason="no panel snapshot; run analysis/xenium/scripts/fetch_subpanels.py")
def test_the_panel_snapshot_matches_the_core_panels_it_was_generated_from():
    """The de-duplication's guard.

    `core/panels` is the one authority for the panel taxonomy and the snapshot under
    `data/development/xenium_panels/` is a derived, gitignored cache. A path is not a pin, so
    the stamp records `git rev-parse HEAD:core/panels/data` -- the tree object id of the
    taxonomy itself, which moves if and only if a byte of it moves. Comparing them means a
    snapshot generated before an edit to `core/panels` fails here rather than silently
    disagreeing with its own source, which is exactly the drift the upstream script's pinned
    SHA existed to prevent.
    """
    live = subprocess.run(["git", "-C", str(REPO), "rev-parse", "HEAD:core/panels/data"],
                          capture_output=True, text=True, check=True).stdout.strip()
    stamped = _snapshot_tree_id()
    assert stamped == live, (
        f"panel snapshot is stale: stamped {stamped}, core/panels/data is now {live}. "
        "Re-run analysis/xenium/scripts/fetch_subpanels.py."
    )


@pytest.mark.skipif(not (SNAPSHOT / "panel_roles.csv").exists(),
                    reason="no panel snapshot; run analysis/xenium/scripts/fetch_subpanels.py")
def test_the_derived_allowlist_agrees_with_each_tissue_module():
    """`_panel_roles.assert_consistent` — the upstream drift guard, run in CI.

    It fails loud when the vendored manifest disagrees with a tissue module's
    SUBPANEL_ALLOWLIST. Running it here means a bad materialisation is caught by the test
    suite instead of by a phenotyping run three hours in.
    """
    import importlib
    import sys
    sys.path.insert(0, str(ROOT / "scripts"))
    from tissue_markers import _panel_roles as PR
    for tissue in ("pancreas", "thymus"):
        PR.assert_consistent(importlib.import_module(f"tissue_markers.{tissue}"), tissue)


# --------------------------------------------------------------------------- utils imports
def test_xenium_explorer_export_imports_without_scanpy():
    import importlib
    import sys
    sys.path.insert(0, str(ROOT))
    mod = importlib.import_module("utils.xenium_explorer_export")
    for name in ("export_for_xenium_explorer", "export_groups_to_csv",
                 "export_groups_to_zarr", "generate_color_palette"):
        assert hasattr(mod, name), f"utils.xenium_explorer_export.{name} missing"
