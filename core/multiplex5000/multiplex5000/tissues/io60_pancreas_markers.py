"""Marker eligibility for IO60 pancreas protein identity and its derived representations."""

# The channel remains in measurement records but must not name, veto or refine a cell.
EXCLUDED_FROM_TYPING = frozenset({"EpCAM"})
