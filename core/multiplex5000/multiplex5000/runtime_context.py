"""Read-only runtime and location reporting, without executing scientific work."""
from __future__ import annotations

import hashlib
import json
import platform
import stat
import subprocess
import sys
from pathlib import Path

from . import __version__, roots


def path_status(path: Path) -> dict:
    """Describe availability on this host; absence here is not global data loss."""
    path = path.absolute()
    result = {"path": str(path)}
    try:
        info = path.stat()
    except FileNotFoundError:
        return {**result, "status": "not_found_on_this_host"}
    except OSError as error:
        return {**result, "status": "inaccessible", "error": str(error)}
    return {**result, "status": "available", "bytes": info.st_size,
            "mtime_ns": info.st_mtime_ns, "regular_file": stat.S_ISREG(info.st_mode)}


def runtime_context() -> dict:
    git = {"available": False, "revision": None, "dirty": None, "changes": []}
    # A wheel's bundled source snapshot is not its containing directory's Git tree.
    if (roots.REPO / ".git").exists():
        try:
            prefix = ["git", "--no-optional-locks", "-C", str(roots.REPO)]
            revision = subprocess.run([*prefix, "rev-parse", "HEAD"], check=True,
                                      capture_output=True, text=True, timeout=5).stdout.strip()
            status = subprocess.run([*prefix, "status", "--porcelain=v1", "--untracked-files=normal"],
                                    check=True, capture_output=True, text=True, timeout=5).stdout
            git.update(available=True, revision=revision, dirty=bool(status),
                       changes=status.splitlines())
        except (OSError, subprocess.SubprocessError) as error:
            git["error"] = str(error)
    provenance_path = roots.REPO / "_build_provenance.json"
    build_provenance = None
    if provenance_path.is_file():
        try:
            contents = provenance_path.read_bytes()
            build_provenance = {"path": str(provenance_path),
                                "sha256": hashlib.sha256(contents).hexdigest(),
                                "record": json.loads(contents)}
        except (OSError, ValueError) as error:
            build_provenance = {"path": str(provenance_path), "error": str(error)}
    return {"executable": sys.executable, "python": platform.python_version(),
            "package_version": __version__, "package_file": str(Path(__file__).with_name("__init__.py")),
            "repository": str(roots.REPO), "git": git, "build_provenance": build_provenance,
            "site_file": str(roots.SITE["path"].resolve()) if roots.SITE else None,
            "site": roots.SITE["site"] if roots.SITE else None,
            "data_root": str(roots.DATA_ROOT), "data_root_source": roots.DATA_ROOT_SOURCE}
