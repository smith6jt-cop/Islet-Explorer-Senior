"""Where the data lives, spelled once -- and set by the site, not by where the code is installed.

Two kinds of root live here, and they are found differently on purpose:

* **Tracked inventories** (``INTEGRATION``) are versioned with the code, so they follow the checkout
  this package is installed from: the editable install of ``core/multiplex5000``, whose ``parents[3]``
  is the repository root.
* **Data roots** (``DATA_ROOT``, ``DEVELOPMENT``, ``TYPING``, the layers directory) come from a
  **site file**: the JSON file named by ``$MPLX5K_SITE_FILE``, else
  ``~/.config/multiplex5000/site.json``. Deriving them from the install location was a hazard: on
  2026-10-03 the shared ``.venv`` was re-installed editable from the v3 worktree, and from then on
  these roots pointed into the v3 worktree's ``data/`` for every checkout using that ``.venv``. A
  site file names the data by absolute path, so which checkout the code comes from no longer
  decides which data it reads (MultiPlex5000 plan, phase D1).

With no site file the data roots fall back to this checkout's ``data/`` -- the right answer for CI
and a fresh clone -- and a ``SiteFallbackWarning`` says so, because a fallback must never be
silent. A site file that is named but missing, malformed, relative, or pointing at a folder that
does not exist is refused: an explicit setting that cannot be honoured is an error, not a fallback.
``scripts/env_check.py`` reports which source is in effect.

These are DATA roots, not code locations. A module that imports its sibling scripts keeps finding
them from its own ``__file__``; routing code paths through here would make a worktree run the main
checkout's code.

Measured 2026-09-21: 197 files compute their own data roots from ``Path(__file__)`` in seven
spellings. They move to these roots as their stage moves into MultiPlex5000 (phases D2-D11).
"""

from __future__ import annotations

import json
import os
import warnings
from pathlib import Path
from typing import Optional

#: The repository root of the checkout this package is installed from (a CODE location).
REPO = Path(__file__).resolve().parents[3]

#: Hand-recorded inventories and the generated vocabulary crosswalk -- tracked, so they follow the
#: checkout.
INTEGRATION = REPO / "core" / "phenocycler" / "data" / "integration"

#: The environment variable naming a site file, and the per-user file used when it is unset.
SITE_FILE_ENV = "MPLX5K_SITE_FILE"
DEFAULT_SITE_FILE = Path.home() / ".config" / "multiplex5000" / "site.json"

#: The one site-file format this module reads.
SITE_SCHEMA = "multiplex5000-site/1"

#: The layer artifacts the app's cards read, unless ``$ISLET_LAYERS_DIR`` or the site names another.
LAYERS_DEFAULT = Path("/data/islet_layers")


class SiteError(RuntimeError):
    """A site file that is named but cannot be honoured."""


class SiteFallbackWarning(UserWarning):
    """No site file: the data roots fall back to this checkout's ``data/``."""


def site_file() -> Optional[Path]:
    """The site file in effect: ``$MPLX5K_SITE_FILE`` if set (it must exist), else the per-user file
    if present, else ``None``."""
    named = os.environ.get(SITE_FILE_ENV, "").strip()
    if named:
        path = Path(named).expanduser()
        if not path.is_file():
            raise SiteError(f"${SITE_FILE_ENV} names {path}, which is not a file")
        return path
    return DEFAULT_SITE_FILE if DEFAULT_SITE_FILE.is_file() else None


def read_site(path: Path) -> dict:
    """Parse and check one site file; refuse anything it cannot honour."""
    try:
        data = json.loads(Path(path).read_text())
    except (OSError, ValueError) as err:
        raise SiteError(f"{path}: cannot be read as JSON ({err})") from err
    if not isinstance(data, dict) or data.get("schema") != SITE_SCHEMA:
        raise SiteError(f"{path}: expected \"schema\": \"{SITE_SCHEMA}\"")
    for key in ("site", "data_root"):
        if not isinstance(data.get(key), str) or not data[key].strip():
            raise SiteError(f"{path}: \"{key}\" is required and must be a non-empty string")
    data_root = Path(data["data_root"])
    if not data_root.is_absolute():
        raise SiteError(f"{path}: data_root {data_root} must be an absolute path")
    if not data_root.is_dir():
        raise SiteError(f"{path}: data_root {data_root} is not a directory")
    layers = data.get("layers_dir")
    if layers is not None and not (isinstance(layers, str) and Path(layers).is_absolute()):
        raise SiteError(f"{path}: layers_dir, when given, must be an absolute path")
    releases = data.get("releases_root")
    if releases is not None and not (isinstance(releases, str) and Path(releases).is_absolute()):
        raise SiteError(f"{path}: releases_root, when given, must be an absolute path")
    snakemake = data.get("snakemake")
    if snakemake is not None and not (isinstance(snakemake, str) and Path(snakemake).is_absolute()):
        raise SiteError(f"{path}: snakemake, when given, must be an absolute path")
    sources = data.get("sources", {})
    if not isinstance(sources, dict) or not all(
            isinstance(k, str) and isinstance(v, str) and Path(v).is_absolute()
            for k, v in sources.items()):
        raise SiteError(f"{path}: sources, when given, must map names to absolute paths")
    return {"site": data["site"], "data_root": data_root, "path": Path(path),
            "layers_dir": Path(layers) if layers else None,
            "releases_root": Path(releases) if releases else None,
            "snakemake": Path(snakemake) if snakemake else None,
            "sources": {k: Path(v) for k, v in sources.items()}}


_SITE_PATH = site_file()

#: The parsed site file, or ``None`` when the checkout fallback is in effect.
SITE: Optional[dict] = read_site(_SITE_PATH) if _SITE_PATH else None

if SITE is None:
    #: The repository-level ``data/`` directory: the site's, else this checkout's.
    DATA_ROOT = REPO / "data"
    #: Where ``DATA_ROOT`` came from, in words, for reports and the environment check.
    DATA_ROOT_SOURCE = (f"checkout fallback ({REPO}): no site file -- set ${SITE_FILE_ENV} or "
                        f"write {DEFAULT_SITE_FILE}")
    warnings.warn(f"phenocycler.roots: {DATA_ROOT_SOURCE}", SiteFallbackWarning, stacklevel=2)
else:
    DATA_ROOT = SITE["data_root"]
    DATA_ROOT_SOURCE = f"site {SITE['site']} ({SITE['path']})"

#: Mutable development outputs, ``data/development`` under the data root.
DEVELOPMENT = DATA_ROOT / "development"

#: The typing lane's outputs: ``cells``, ``app``, ``immune``, ``geometry`` and the rest.
TYPING = DEVELOPMENT / "typing"


#: Where releases are built and kept when the site does not say (the MultiPlex5000 plan's root).
RELEASES_DEFAULT = Path("/data/io60-current/releases")
RELEASES_ROOT_ENV = "MPLX5K_RELEASES_ROOT"


def releases_root() -> Path:
    """``$MPLX5K_RELEASES_ROOT``, else the site's ``releases_root``, else
    ``/data/io60-current/releases``; read when called, not at import."""
    named = os.environ.get(RELEASES_ROOT_ENV)
    if named:
        return Path(named)
    if SITE is not None and SITE.get("releases_root") is not None:
        return SITE["releases_root"]
    return RELEASES_DEFAULT


def source(name: str) -> Optional[Path]:
    """Where this site keeps a recorded source that lives outside the repository (the site file's
    ``sources``), or ``None`` when the site does not say -- which is unknown, not absent."""
    if SITE is None:
        return None
    return SITE["sources"].get(name)


def layers_dir() -> Path:
    """``$ISLET_LAYERS_DIR``, else the site's ``layers_dir``, else ``/data/islet_layers``; read when
    called, not at import."""
    named = os.environ.get("ISLET_LAYERS_DIR")
    if named:
        return Path(named)
    if SITE is not None and SITE["layers_dir"] is not None:
        return SITE["layers_dir"]
    return LAYERS_DEFAULT
