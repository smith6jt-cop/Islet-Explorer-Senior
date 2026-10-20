"""``phenocycler`` is a temporary import alias for ``multiplex5000``.

The package was renamed with its directory on 2026-10-20 (MultiPlex5000 D1):
``core/phenocycler/phenocycler`` became ``core/multiplex5000/multiplex5000``. Until every importer
uses the new name, ``import phenocycler`` returns the ``multiplex5000`` package itself and
``import phenocycler.<module>`` returns the module ``multiplex5000.<module>`` -- the same objects,
not copies -- so module state, classes, monkeypatches and ``importlib.reload`` behave as they did.
``python -m phenocycler.<module>`` runs ``multiplex5000.<module>``.

Nothing inside ``multiplex5000`` imports this name, and ``tests/current/test_phenocycler_alias.py``
holds it to that. New code imports ``multiplex5000``. This package goes when ``git grep`` finds no
importer of ``phenocycler``.
"""

from __future__ import annotations

import importlib
import importlib.abc
import importlib.machinery
import importlib.util
import sys

import multiplex5000

_OLD = __name__
_NEW = multiplex5000.__name__


def _real_name(fullname: str) -> str:
    return _NEW + fullname[len(_OLD):]


class _Alias(importlib.abc.MetaPathFinder, importlib.abc.Loader):
    """Find ``phenocycler.<x>`` and load it as the already-importable ``multiplex5000.<x>``."""

    is_phenocycler_alias = True

    def find_spec(self, fullname, path=None, target=None):
        if not fullname.startswith(_OLD + "."):
            return None
        real_spec = importlib.util.find_spec(_real_name(fullname))
        if real_spec is None:
            return None
        spec = importlib.machinery.ModuleSpec(
            fullname, self, origin=real_spec.origin,
            is_package=real_spec.submodule_search_locations is not None)
        spec.loader_state = _real_name(fullname)
        return spec

    def create_module(self, spec):
        module = importlib.import_module(spec.loader_state)
        # Keep the module's own spec: the import system is about to point __spec__ at this alias.
        spec.loader_state = module.__spec__
        return module

    def exec_module(self, module):
        module.__spec__ = module.__spec__.loader_state

    # ``python -m phenocycler.<module>``: runpy asks the loader for the module's code.
    def get_code(self, fullname):
        real = _real_name(fullname)
        return importlib.util.find_spec(real).loader.get_code(real)

    def get_source(self, fullname):
        real = _real_name(fullname)
        return importlib.util.find_spec(real).loader.get_source(real)


# First on the path, so a submodule is never found a second time under the old name.
if not any(getattr(finder, "is_phenocycler_alias", False) for finder in sys.meta_path):
    sys.meta_path.insert(0, _Alias())

sys.modules[__name__] = multiplex5000
