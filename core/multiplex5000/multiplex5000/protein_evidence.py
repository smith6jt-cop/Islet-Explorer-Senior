"""Review-only protein positivity evidence within one physical section.

The existing IO60 typing lane requires a marker to exceed 0.10 of its span
from the measured background mode to a robust bright end. Here that candidate
is fitted separately for each marker and section, as the maintainer selected
for the pLN/spleen pilot on 2026-10-26. It consumes Stage 05's UniFORM shift;
there is no pooled cutoff, positive-share target, or final cell-type call.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from .rna_evidence import DETECTED, UNAVAILABLE, UNDETECTED

REQUIRE_MIN = 0.10
ROBUST_TOP_Q = 99.9
MODE_BINS = 512
MIN_BACKGROUND_SUPPORT = 1_000


class ProteinEvidenceError(ValueError):
    """A section's marker values cannot be interpreted as declared."""


@dataclass(frozen=True)
class ProteinEvidence:
    section: str
    marker: str
    states: np.ndarray
    floor: float | None
    top: float | None
    cutoff: float | None
    shift_scale: float | None
    n_background_support: int
    unavailable_reason: str | None


def candidate_calls(section: str, marker: str, values: np.ndarray | None,
                    nucleus_eligible: np.ndarray, *,
                    shift_scale: float | None) -> ProteinEvidence:
    """Fit one marker's candidate cutoff using this section's eligible cells only.

    ``None`` values or shift mean unavailable, not negative. Individual null
    or negative measurements stay unavailable too. At least 1,000 positive
    finite values are needed to locate a measured background mode. The 0.10
    fraction is a review candidate; Stage 06 controls must be reviewed before
    any cell classes can be accepted.
    """
    if not section or not marker:
        raise ProteinEvidenceError("physical section and marker must be named")
    eligible = np.asarray(nucleus_eligible)
    if eligible.ndim != 1 or eligible.dtype != np.bool_:
        raise ProteinEvidenceError("nucleus eligibility must be one-dimensional Boolean")
    states = np.full(len(eligible), UNAVAILABLE, dtype=np.int8)

    def unavailable(reason: str, support: int = 0) -> ProteinEvidence:
        return ProteinEvidence(section, marker, states, None, None, None,
                               shift_scale, support, reason)

    if values is None:
        return unavailable("marker_not_measured")
    raw = np.asarray(values, dtype=np.float64)
    if raw.shape != eligible.shape:
        raise ProteinEvidenceError("marker values and nucleus admission differ in cell order")
    if shift_scale is None:
        return unavailable("uniform_shift_unavailable")
    if not np.isfinite(shift_scale) or shift_scale <= 0:
        raise ProteinEvidenceError("UniFORM shift must be positive and finite")
    aligned = raw * shift_scale
    valid = eligible & np.isfinite(aligned) & (aligned >= 0)
    positive = aligned[valid & (aligned > 0)]
    support = len(positive)
    if support < MIN_BACKGROUND_SUPPORT:
        return unavailable("insufficient_background_support", support)
    log = np.log1p(positive)
    if float(log.min()) == float(log.max()):
        return unavailable("no_intensity_range", support)
    counts, edges = np.histogram(log, bins=MODE_BINS)
    # Zero intensities were excluded before histogramming, so the first bin
    # may contain the actual background mode and must remain eligible.
    # No donor contributes values to another section's floor or bright end.
    floor = float(np.expm1(edges[int(np.argmax(counts)) + 1]))
    top = float(np.percentile(aligned[valid], ROBUST_TOP_Q))
    if not (np.isfinite(floor) and np.isfinite(top) and top > floor):
        return unavailable("bright_end_not_above_background", support)
    cutoff = floor + REQUIRE_MIN * (top - floor)
    states[valid] = np.where(aligned[valid] >= cutoff, DETECTED, UNDETECTED)
    return ProteinEvidence(section, marker, states, floor, top, cutoff,
                           shift_scale, support, None)
