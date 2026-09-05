"""The pre-registered bar.  Written before any arm ran; never edited after a result is seen.

No key encodes a share of cells, a prior, a floor or cap on a label's frequency, or any donor
or disease ordering.  ``assert_frequency_free`` is the executable form of that sentence and
the test suite calls it; a new key must pass it or the test fails.
"""
from __future__ import annotations

CRITERIA = {
    # Step 1's tail definition, verbatim (scripts/celltype/usability.py tail_overlap: k = max(20, q*n)).
    "tail_q": 0.001,
    "tail_min_k": 20,
    # Evidence floors: on a statistic's n, never on a share.  Unjudgeable is never a pass.
    "min_calls": 50,
    "min_tail_in_pool": 20,
    # Nulls.
    "null_draws": 200,
    "null_quantile": 0.95,
    # Spatial coherence graph (usability.BARS morans_k).
    "morans_k": 12,
    # Stability.
    "n_replicates": 5,
    "primary_seed": 0,
    "subsample_frac": 0.9,
    # Frame: 8 blind + 8 islet-anchored 1024-px cores per donor.
    "n_blind": 8,
    "n_islet": 8,
    "core_px": 1024,
    "window_px": 2048,
    # Verification and survival fractions (fractions of JUDGED items, never of cells).
    "heldout_frac": 2 / 3,
    "survive_frac": 2 / 3,
}

FORBIDDEN_KEY_TOKENS = ("share", "prevalence", "frequency", "fraction_of_cells", "prior",
                        "expected", "min_share", "max_share", "collapse", "order", "status", "disease")


def assert_frequency_free(criteria: dict = CRITERIA) -> None:
    for k in criteria:
        low = k.lower()
        if any(tok in low for tok in FORBIDDEN_KEY_TOKENS):
            raise AssertionError(f"criterion key {k!r} looks like an expected-frequency term")
    for k in ("heldout_frac", "survive_frac", "subsample_frac"):
        v = criteria[k]
        if not (0 < v <= 1):
            raise AssertionError(f"{k} must be a fraction of judged items in (0, 1], got {v}")


assert_frequency_free()
