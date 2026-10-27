"""Build the tracked panel and workflow knowledge into an installable package.

The checkout keeps a single source under core/panels and core/multiplex5000.
Wheel and sdist builds take a snapshot of those tracked files; no second
maintained copy is kept in the repository.
"""

import json
import subprocess
from pathlib import Path
from shutil import copy2

from setuptools import setup
from setuptools.command.build_py import build_py as _build_py
from setuptools.command.sdist import sdist as _sdist

PROJECT = Path(__file__).resolve().parent
REPOSITORY = PROJECT.parents[1]
SNAPSHOT = PROJECT / "_source_data"
# One tracked copy of the selected pancreas hierarchy remains in scripts while
# its data inputs are moved to Stage 06. A wheel carries a build snapshot of
# precisely that code so a site never imports another checkout's scripts.
PANCREAS_TYPING_SCRIPTS = (
    "core/multiplex5000/scripts/typing/run_calibration_arm.py",
    "core/multiplex5000/scripts/typing/subpanel.py",
    "core/multiplex5000/scripts/typing/build_cells.py",
    "core/multiplex5000/scripts/typing/refine.py",
    "core/multiplex5000/scripts/typing/immune_subtype.py",
    "core/multiplex5000/scripts/typing/subpanel_io.py",
    "core/multiplex5000/scripts/typing/superseded.py",
    "core/multiplex5000/scripts/typing/context.py",
    "core/multiplex5000/scripts/typing/islet_identity.py",
    "core/multiplex5000/scripts/typing/islet_membership.py",
    "core/multiplex5000/scripts/typing/lineages.py",
    "core/multiplex5000/scripts/typing/cohort.py",
    "core/multiplex5000/scripts/typing/build_islet_geometry.py",
    "core/multiplex5000/scripts/celltype/panel.py",
    "core/multiplex5000/scripts/celltype/rules.py",
    "core/multiplex5000/scripts/celltype/score.py",
    "core/multiplex5000/scripts/celltype/run_step3.py",
    "core/multiplex5000/scripts/celltype/alignment_block.py",
    "core/multiplex5000/scripts/celltype/marker_name_crosswalk.json",
    "phenotype_workflow.csv",
)
DIRECTORIES = (
    "core/panels/data",
    "core/phenocycler/data/integration",
    "core/multiplex5000/profiles",
    "core/multiplex5000/projects",
    "core/multiplex5000/sites",
    "core/multiplex5000/notebooks",
)
FILES = (
    "core/multiplex5000/config.ini",
    "core/multiplex5000/cohort-spec.json",
    "core/multiplex5000/cohort-spec-v3.json",
    "core/multiplex5000/multiplex5000/marker_registry.json",
    *PANCREAS_TYPING_SCRIPTS,
)


def _source_root() -> Path:
    if (REPOSITORY / "core/panels/data").is_dir():
        return REPOSITORY
    if (SNAPSHOT / "core/panels/data").is_dir():
        return SNAPSHOT
    raise RuntimeError("MultiPlex5000 needs its tracked panel data to build")


def _knowledge_files() -> list[tuple[Path, Path]]:
    root = _source_root()
    if root == REPOSITORY:
        # A workstation can hold scratch files beside a panel or notebook.
        # Only the committed source tree may enter a release artifact.
        result = subprocess.run(
            ["git", "-C", str(root), "ls-files", "-z", "--", *DIRECTORIES, *FILES],
            check=True, capture_output=True)
        paths = [root / name.decode() for name in result.stdout.split(b"\0") if name]
    else:
        paths = [root / name for name in FILES]
        for name in DIRECTORIES:
            paths.extend(path for path in (root / name).rglob("*") if path.is_file())
    missing = [str(path.relative_to(root)) for path in paths if not path.is_file()]
    missing.extend(name for name in FILES if root / name not in paths)
    if missing:
        raise RuntimeError(f"MultiPlex5000 build lacks tracked resources: {missing}")
    return [(path.relative_to(root), path) for path in sorted(set(paths))]


def _provenance() -> dict:
    root = _source_root()
    if root == SNAPSHOT:
        return json.loads((SNAPSHOT / "_build_provenance.json").read_text())
    commit = subprocess.run(["git", "-C", str(root), "rev-parse", "HEAD"],
                            check=True, capture_output=True, text=True).stdout.strip()
    tracked = subprocess.run(
        ["git", "-C", str(root), "status", "--porcelain", "--",
         "core/multiplex5000", "core/panels/data", "core/phenocycler/data/integration"],
        check=True, capture_output=True, text=True).stdout
    return {"commit": commit, "dirty": bool(tracked.strip())}


def _write_provenance(target: Path) -> None:
    target.mkdir(parents=True, exist_ok=True)
    (target / "_build_provenance.json").write_text(json.dumps(_provenance(), sort_keys=True) + "\n")


class build_py(_build_py):
    def run(self):
        super().run()
        target = Path(self.build_lib) / "multiplex5000" / "_source"
        for relative, source in _knowledge_files():
            destination = target / relative
            destination.parent.mkdir(parents=True, exist_ok=True)
            copy2(source, destination)
        _write_provenance(target)

    def get_outputs(self, include_bytecode=1):
        output = list(super().get_outputs(include_bytecode))
        target = Path(self.build_lib) / "multiplex5000" / "_source"
        return output + [str(target / relative) for relative, _ in _knowledge_files()] + [
            str(target / "_build_provenance.json")]


class sdist(_sdist):
    def make_release_tree(self, base_dir, files):
        super().make_release_tree(base_dir, files)
        target = Path(base_dir) / "_source_data"
        for relative, source in _knowledge_files():
            destination = target / relative
            destination.parent.mkdir(parents=True, exist_ok=True)
            copy2(source, destination)
        _write_provenance(target)


setup(cmdclass={"build_py": build_py, "sdist": sdist})
