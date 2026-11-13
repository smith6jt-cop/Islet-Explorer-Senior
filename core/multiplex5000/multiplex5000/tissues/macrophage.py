"""Pancreas macrophage marker profiles in the shared typing hierarchy."""

from __future__ import annotations

from itertools import product

import numpy as np

# The PhenoCycler M2Gal3 antibody reads galectin-3 (LGALS3 in the CODEX panel).
PROFILE_MARKERS = ("CD68", "CD163", "M2Gal3")


def label(signs: tuple[bool, bool, bool]) -> str:
    """Use the owner's CODEX label spelling with Senior's marker alias resolved."""
    a, b, c = ("+" if x else "-" for x in signs)
    return f"CD68{a}CD163{b}Macrophage_LGALS3{c}"


def labels() -> tuple[str, ...]:
    """The six measured profiles with CD68 or CD163 positive."""
    return tuple(label(v) for v in product((True, False), repeat=3) if v[0] or v[1])


def current_profile_names(series):
    """Pool archived CD4-split names into their measured three-marker profiles."""
    legacy = {name + suffix: name for name in labels() for suffix in ("CD4+", "CD4-")}
    return series.astype(object).replace(legacy)


def classify_profile(values, threshold: float = 0.10):
    """Return a profile only when all three markers were measured on the cell."""
    block = values.reindex(columns=PROFILE_MARKERS).to_numpy(dtype=float)
    measured = np.isfinite(block).all(axis=1)
    positive = block >= threshold
    eligible = measured & (positive[:, 0] | positive[:, 1])
    result = np.full(len(values), None, dtype=object)
    result[eligible] = [label(tuple(row)) for row in positive[eligible]]
    return result
