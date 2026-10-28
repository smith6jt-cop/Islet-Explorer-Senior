"""Pancreas macrophage marker profiles in the shared typing hierarchy."""

from __future__ import annotations

from collections.abc import Mapping
from itertools import product

import numpy as np

# The PhenoCycler M2Gal3 antibody reads galectin-3 (LGALS3 in the CODEX panel).
PROFILE_MARKERS = ("CD68", "CD163", "M2Gal3", "CD4")


def label(signs: tuple[bool, bool, bool, bool]) -> str:
    """Use the owner's CODEX label spelling with Senior's marker alias resolved."""
    a, b, c, d = ("+" if x else "-" for x in signs)
    return f"CD68{a}CD163{b}Macrophage_LGALS3{c}CD4{d}"


def labels() -> tuple[str, ...]:
    """The 12 measured profiles with CD68 or CD163 positive."""
    return tuple(label(v) for v in product((True, False), repeat=4) if v[0] or v[1])


def classify_profile(values, threshold: float | Mapping[str, float] = 0.10):
    """Classify measured profiles using levels in the same units as ``values``.

    A marker-to-level mapping supports image-specific positivity levels. The
    scalar form preserves the committed Senior method until its producer is
    replaced; callers must not infer missing levels from that default.
    """
    block = values.reindex(columns=PROFILE_MARKERS).to_numpy(dtype=float)
    if isinstance(threshold, Mapping):
        missing = sorted(set(PROFILE_MARKERS) - threshold.keys())
        if missing:
            raise ValueError(f"macrophage positivity levels missing markers: {missing}")
        levels = np.array([threshold[marker] for marker in PROFILE_MARKERS], dtype=float)
    else:
        levels = np.full(len(PROFILE_MARKERS), threshold, dtype=float)
    if not np.isfinite(levels).all():
        raise ValueError("macrophage positivity levels must be finite")
    measured = np.isfinite(block).all(axis=1)
    positive = block >= levels
    eligible = measured & (positive[:, 0] | positive[:, 1])
    result = np.full(len(values), None, dtype=object)
    result[eligible] = [label(tuple(row)) for row in positive[eligible]]
    return result
