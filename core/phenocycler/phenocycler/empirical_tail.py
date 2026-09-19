"""Finite-sample upper-tail primitives retained for RESTORE research.

RETIRED -- NOT A FALLBACK.  Recorded here 2026-09-18 because the docstring below reads as a live
specification and is not one.  RESTORE was retired on 2026-08-14 and is on the never-revive list
(root ``CLAUDE.md``): the RESTORE pair webs, the Gate 2-5 apparatus and the GMM/Otsu gate searches
do not come back as fallbacks.  Nothing in this module produced the typing on disk -- that is the
subpanel classification under ``core/phenocycler/scripts/typing/`` -- and the pipeline stages that
would call it (``controls``, ``calibrate``, ``type``) are blocked by ``NO_VERIFIED_TYPING``, a
latch that has not moved.  Kept as a record of what was implemented and measured.

These functions are independent of any marker roster, pair declaration, or
production workflow.  They remain available to the explicit synthetic
RESTORE evaluator while pair design is under development.
"""

from __future__ import annotations

from collections.abc import Sequence

import numpy as np


def _empirical_order_index(n: int, alpha: float) -> int | None:
    """Return the conformal upper-tail order index when ``n`` resolves alpha."""

    if not (0.0 < alpha < 1.0):
        raise ValueError("marker_alpha must be between 0 and 1")
    if alpha < 1.0 / (n + 1):
        return None
    index = int(np.ceil((1.0 - alpha) * (n + 1))) - 1
    return int(np.clip(index, 0, n - 1))


def empirical_upper_threshold(log_values: np.ndarray, alpha: float) -> float:
    """Return the finite-sample empirical upper-tail order statistic."""

    values = np.asarray(log_values, dtype=float)
    if values.ndim != 1 or len(values) < 2:
        raise ValueError("log_values must contain at least two values")
    if not np.isfinite(values).all():
        raise ValueError("log_values must be finite")
    index = _empirical_order_index(len(values), alpha)
    if index is None:
        raise ValueError(f"{len(values)} controls cannot resolve marker alpha {alpha:g}")
    return float(np.partition(values, index)[index])


def empirical_tail_probability(
    log_values: Sequence[float] | np.ndarray,
    sorted_control_log_values: Sequence[float] | np.ndarray,
) -> np.ndarray:
    """Return finite-sample empirical upper-tail probabilities with +1 correction."""

    values = np.asarray(log_values, dtype=float)
    controls = np.asarray(sorted_control_log_values, dtype=float)
    if controls.ndim != 1 or len(controls) < 2:
        raise ValueError("sorted controls must contain at least two values")
    if not np.isfinite(controls).all() or np.any(np.diff(controls) < 0):
        raise ValueError("sorted controls must be finite and nondecreasing")
    left = np.searchsorted(controls, values, side="left")
    return (1.0 + len(controls) - left) / (len(controls) + 1.0)


__all__ = ["empirical_tail_probability", "empirical_upper_threshold"]
