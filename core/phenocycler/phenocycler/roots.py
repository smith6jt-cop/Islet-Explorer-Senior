"""Where the data lives, spelled once.

Every root here is derived from the checkout this package is installed from: the editable install
of ``core/phenocycler``, whose ``parents[3]`` is the repository root.  Data lives only in the main
checkout -- every ``data/`` tree is gitignored -- so a script run from a worktree against the repo
``.venv`` reads the same data as one run from the main checkout.

These are DATA roots, not code locations.  A module that imports its sibling scripts keeps finding
them from its own ``__file__``; routing code paths through here would make a worktree run the main
checkout's code.

Measured 2026-09-21: 197 files compute their own data roots from ``Path(__file__)`` in seven
spellings.  This module is used only by the files the 2026-09-21 plan's Step 6 touched; the rest
are not swept.
"""

from __future__ import annotations

import os
from pathlib import Path

#: The repository root of the checkout this package is installed from.
REPO = Path(__file__).resolve().parents[3]

#: Hand-recorded inventories and the generated vocabulary crosswalk.
INTEGRATION = REPO / "core" / "phenocycler" / "data" / "integration"

#: Mutable development outputs, a level above ``core/``.
DEVELOPMENT = REPO / "data" / "development"

#: The typing lane's outputs: ``cells``, ``app``, ``immune``, ``geometry`` and the rest.
TYPING = DEVELOPMENT / "typing"

#: The layer artifacts the app's cards read, unless ``$ISLET_LAYERS_DIR`` names another.
LAYERS_DEFAULT = Path("/data/islet_layers")


def layers_dir() -> Path:
    """``$ISLET_LAYERS_DIR``, else ``/data/islet_layers``; read when called, not at import."""
    return Path(os.environ.get("ISLET_LAYERS_DIR", str(LAYERS_DEFAULT)))
