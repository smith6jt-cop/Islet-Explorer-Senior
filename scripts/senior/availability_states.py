#!/usr/bin/env python3
"""The availability-state vocabulary, in one place, on the Python side.

WHY THIS FILE EXISTS. `apps/R/ieCore/R/ie_availability.R` declares six states in `IE_STATES`
(line 31). Four Python producers each wrote their own five-element copy of that list:

    scripts/senior/build_availability.py       STATES = (... five ...)
    scripts/senior/build_restain_layer.py      STATES = (... five ...)
    scripts/senior/build_trajectory_layer.py   STATES = {... five ...}
    scripts/senior/build_grid_composition.py   STATES = (... five ...)

Each was missing `uncalled`. A producer cannot write a state its own guard rejects, so a cell
the typer looked at and resolved nothing for had nowhere to go and could only be carried as one
of the other five -- a hole reported as an answer, or an answer reported as a hole. The list is
now derived here and nowhere else, so the two sides cannot drift again.

ORDER IS LOAD-BEARING AND IS NOT ALPHABETICAL. These are factor levels on the R side and a card
legend enumerates them in order. `ie_availability.R:28` records the rule in its own words:
`uncalled` is APPENDED, not inserted, so no existing level moves. This tuple is that order, and
a new state is appended here too.

WHAT THIS FILE DOES NOT DO. Widening a producer's vocabulary lets `uncalled` be carried; it does
not make any producer emit it. Which rows become `uncalled` -- and in which capabilities -- is a
question about what each builder measures, not about this list, and it is deliberately left to
the maintainer. No row's state changes because of this module.
"""

from __future__ import annotations

#: The availability states, in the order a legend walks them.
#:
#: Mirrors `IE_STATES` in `apps/R/ieCore/R/ie_availability.R` exactly, including order. The R
#: side is the declaring authority; this is the Python reading of it.
#:
#:   measured     the value exists
#:   blank        never typed
#:   unavailable  measured, then dropped at QC -- something WAS measured; it did not survive
#:   withheld     measurable and deliberately not released
#:   failed       attempted and failed
#:   uncalled     typed, and the typing resolved no identity. NOT an absence of measurement and
#:                NOT an absence of typing: the producer looked, and recorded that it found
#:                nothing. Appended 2026-09-18.
AVAILABILITY_STATES: tuple[str, ...] = (
    "measured",
    "blank",
    "unavailable",
    "withheld",
    "failed",
    "uncalled",
)

#: The one state that carries a count. Every other state carries a reason and a null count --
#: a missing measurement is unavailable, never zero.
STATE_MEASURED = "measured"
STATE_BLANK = "blank"
STATE_UNAVAILABLE = "unavailable"
STATE_WITHHELD = "withheld"
STATE_FAILED = "failed"
STATE_UNCALLED = "uncalled"

__all__ = [
    "AVAILABILITY_STATES",
    "STATE_MEASURED", "STATE_BLANK", "STATE_UNAVAILABLE",
    "STATE_WITHHELD", "STATE_FAILED", "STATE_UNCALLED",
]
