"""The tissue hooks of the per-image typing method, and the check that a registry only names what its tissue has.

Every typing tissue module in ``multiplex5000.tissues`` supplies:
- ``start_class(cells) -> array``: each cell's start class from its context (pancreas: its endocrine object's class;
  spleen: its region object), Tissue where it has none;
- ``reference_populations(registry) -> frozenset``: the named populations (``levels.populations``) the tissue's
  registry may use as a level or check reference;
- ``selection_sections(roster) -> list``: the sections options are chosen on (pancreas: the ND images; spleen: all);
- ``structure_objects() -> tuple``: the structure classes whose objects can place cells.
"""
from __future__ import annotations

import importlib
from types import ModuleType

from . import registry as RG

HOOKS = ("start_class", "reference_populations", "selection_sections", "structure_objects")


def hooks(tissue: str) -> ModuleType:
    module = importlib.import_module(f"multiplex5000.tissues.{tissue}")
    missing = [h for h in HOOKS if not callable(getattr(module, h, None))]
    if missing:
        raise RG.RegistryError(f"tissue {tissue!r} has no typing hook(s) {missing}")
    return module


def check_registry(reg: RG.Registry, module: ModuleType) -> None:
    """Refuse a registry that names a reference population its tissue does not have."""
    allowed = module.reference_populations(reg)
    named = set()
    for column in ("level_reference", "check_reference"):
        if column in reg.table.columns:
            for value in reg.table[column]:
                named.update(RG.split(str(value).replace("+", ";")))
    unknown = sorted(named - set(allowed))
    if unknown:
        raise RG.RegistryError(f"{reg.name}: reference populations its tissue does not have: {unknown}")
