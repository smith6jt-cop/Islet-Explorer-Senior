"""Projects: the release being built, its cohort, the stages that run, and which gates are automatic.

A project file is one small YAML file, ``core/multiplex5000/projects/<id>.yaml``. Both routes
read it -- the Snakemake graph and the notebooks -- so they build the same release the same way.
Releases are built under the site's ``releases_root`` (``/data/io60-current/releases`` unless the
site says otherwise), one folder per release, one folder per stage inside it.

MultiPlex5000 D3 (2026-10-20).
"""

from __future__ import annotations

import re
from dataclasses import dataclass, field
from pathlib import Path
from typing import Optional

import yaml

from . import roots

PROJECT_SCHEMA = "multiplex5000-project/1"
PROJECTS_DIR = Path(__file__).resolve().parents[1] / "projects"
_RELEASE_ID = re.compile(r"^[a-z0-9][a-z0-9._-]{0,79}$")
_REQUIRED = ("id", "release", "cohort", "stages", "automatic_gates")
_OPTIONAL = ("settings", "sources", "releases_root")


class ProjectError(ValueError):
    """A project file this package cannot honour."""


@dataclass(frozen=True)
class Project:
    path: Path
    id: str
    release: str
    cohort: str
    stages: tuple[str, ...]
    automatic_gates: tuple[str, ...]
    releases_root: Path
    settings: dict = field(default_factory=dict)
    sources: dict = field(default_factory=dict)

    @property
    def release_dir(self) -> Path:
        return self.releases_root / self.release

    def previous(self, stage: str) -> Optional[str]:
        """The stage before ``stage`` in this project, or None for the first."""
        if stage not in self.stages:
            raise ProjectError(f"{stage} is not one of this project's stages {list(self.stages)}")
        index = self.stages.index(stage)
        return self.stages[index - 1] if index else None

    def source(self, name: str) -> Optional[Path]:
        """A recorded source outside the repository: the project's own path, else the site's."""
        if name in self.sources:
            return Path(self.sources[name])
        return roots.source(name)


def resolve(name_or_path) -> Path:
    candidate = Path(name_or_path)
    if candidate.suffix in (".yaml", ".yml") or candidate.is_file():
        return candidate.resolve()
    return (PROJECTS_DIR / f"{name_or_path}.yaml").resolve()


def load(name_or_path="io60", *, releases_root=None) -> Project:
    """Read and check a project file, by id (``projects/<id>.yaml``) or by path."""
    from . import stages as _stages

    path = resolve(name_or_path)
    try:
        data = yaml.safe_load(path.read_text())
    except (OSError, yaml.YAMLError) as err:
        raise ProjectError(f"{path}: cannot be read ({err})") from err
    if not isinstance(data, dict) or data.get("schema") != PROJECT_SCHEMA:
        raise ProjectError(f"{path}: expected schema: {PROJECT_SCHEMA}")
    body = {k: v for k, v in data.items() if k != "schema"}
    unknown = sorted(set(body) - set(_REQUIRED) - set(_OPTIONAL))
    missing = [k for k in _REQUIRED if k not in body]
    if unknown or missing:
        raise ProjectError(f"{path}: unknown key(s) {unknown}, missing key(s) {missing}")
    if body["id"] != path.stem:
        raise ProjectError(f"{path}: id {body['id']!r} must equal the file name {path.stem!r}")
    if not isinstance(body["release"], str) or not _RELEASE_ID.match(body["release"]):
        raise ProjectError(f"{path}: release must be lower-case letters, digits, '.', '_' or '-'")
    for key in ("stages", "automatic_gates"):
        if not isinstance(body[key], list) or not all(isinstance(s, str) for s in body[key]):
            raise ProjectError(f"{path}: {key} must be a list of stage ids")
    if not body["stages"]:
        raise ProjectError(f"{path}: a project runs at least one stage")
    order = [s for s, _ in _stages.PLAN]
    for stage in body["stages"]:
        _stages.get(stage)                       # unknown or not yet built: refused, with why
    if [order.index(s) for s in body["stages"]] != sorted(order.index(s) for s in body["stages"]):
        raise ProjectError(f"{path}: stages must keep the plan's order {order}")
    extra = sorted(set(body["automatic_gates"]) - set(body["stages"]))
    if extra:
        raise ProjectError(f"{path}: automatic_gates names stage(s) the project does not run: "
                           f"{extra}")
    from . import profiles

    if ("cohort", body["cohort"]) not in profiles.load_all():
        raise ProjectError(f"{path}: no cohort profile {body['cohort']!r}")
    settings = body.get("settings") or {}
    sources = body.get("sources") or {}
    if not isinstance(settings, dict) or not isinstance(sources, dict):
        raise ProjectError(f"{path}: settings and sources are mappings")
    for name, value in sources.items():
        if not (isinstance(value, str) and Path(value).is_absolute()):
            raise ProjectError(f"{path}: source {name} must be an absolute path")
    root = releases_root or body.get("releases_root")
    if root is not None and not Path(root).is_absolute():
        raise ProjectError(f"{path}: releases_root must be an absolute path")
    return Project(path, body["id"], body["release"], body["cohort"], tuple(body["stages"]),
                   tuple(body["automatic_gates"]), Path(root) if root else roots.releases_root(),
                   dict(settings), {k: Path(v) for k, v in sources.items()})
