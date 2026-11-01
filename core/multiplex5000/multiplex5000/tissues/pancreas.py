"""Pancreas hooks of the per-image typing method (``multiplex5000.positivity``).

The pancreas-only parts of the CODEX method, behind the four hooks every typing tissue supplies:
- ``start_class``: the class of the cell's endocrine object (Islet, EndocrineCluster, EndocrineSingle), else Tissue;
- ``reference_populations``: the named populations a pancreas registry may use as references, including the
  endocrine and exocrine ones (hormone-positive endocrine cells, exocrine and acinar epithelium);
- ``selection_sections``: the ND images (owner, Senior 10-01: test classifiers on the non-diabetic controls, because
  disease changes the reference populations);
- ``structure_objects``: the structure classes whose objects place cells (Blood vessel, Nerve, Lymphatic).
"""
from __future__ import annotations

from collections.abc import Mapping

import numpy as np
import pandas as pd

POPULATIONS = frozenset({"hormone", "leukocyte", "exocrine", "epithelial_away", "structure_cells", "t_cells",
                         "acinar_away", "exocrine_clean", "leukocyte_clean"})
STRUCTURES = ("Blood vessel", "Nerve", "Lymphatic")


def start_class(cells: pd.DataFrame) -> np.ndarray:
    return cells["endocrine_object_class"].fillna("Tissue").to_numpy(dtype=object)


def reference_populations(registry) -> frozenset[str]:
    return POPULATIONS


def selection_sections(roster: Mapping[str, Mapping[str, str]]) -> list[str]:
    """`roster`: {section: {"donor_status": ...}}; the ND sections, in roster order."""
    return [s for s, r in roster.items() if r.get("donor_status") == "ND"]


def structure_objects() -> tuple[str, ...]:
    return STRUCTURES
