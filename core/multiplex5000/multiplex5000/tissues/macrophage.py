"""Pancreas macrophage identity and Gal3 subtypes in the protein hierarchy."""

from __future__ import annotations

from collections.abc import Mapping

import numpy as np

PROFILE_MARKERS = ("CD68", "CD163", "M2Gal3")


def labels() -> tuple[str, str]:
    return ("Macrophage_Gal3+", "Macrophage_Gal3-")


def archived_profile_groups(name):
    """The two former schemas both reported three mutually exclusive anchor profiles."""
    sign = {labels()[0]: "+", labels()[1]: "-"}[name]
    six = [anchor + "Macrophage_LGALS3" + sign
           for anchor in ("CD68+CD163+", "CD68+CD163-", "CD68-CD163+")]
    twelve = [profile + suffix for profile in six for suffix in ("CD4+", "CD4-")]
    return six, twelve


def current_profile_names(series):
    """Pool either archived schema; retain missing and unrecognized source values."""
    old = {prior: name for name in labels() for group in archived_profile_groups(name) for prior in group}
    return series.astype(object).replace(old)


def classify_profile(values, threshold: float | Mapping[str, float] = 0.10):
    """A measured CD68+ or CD163+ establishes the parent. M2Gal3 determines its subtype."""
    block = values.reindex(columns=PROFILE_MARKERS).to_numpy(dtype=float)
    if isinstance(threshold, Mapping):
        levels = np.array([threshold.get(marker, np.nan) for marker in PROFILE_MARKERS], dtype=float)
    else:
        levels = np.full(len(PROFILE_MARKERS), threshold, dtype=float)
    # A missing channel or channel-specific level cannot serve as a witness.
    # Another calibrated positive anchor independently establishes the parent.
    measured = np.isfinite(block) & np.isfinite(levels)
    positive = measured & (block >= levels)
    parent = positive[:, 0] | positive[:, 1]
    result = np.full(len(values), None, dtype=object)
    result[parent & measured[:, 2]] = np.where(positive[parent & measured[:, 2], 2], *labels())
    return result
