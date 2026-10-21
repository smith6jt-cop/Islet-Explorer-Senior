"""How a recorded panel is written in every output (maintainer, 2026-10-21).

* **PhenoCycler** -- the number of markers the image records, DAPI not counted; then ``+`` or
  ``-`` and the one marker that tells two panels of that size apart, as the image spells it; or
  ``MetSI`` where the metabolic spike-in is in the panel; a number with neither is followed by the
  word ``markers``. So the IO60 panels are ``57 markers``, ``58+ CD38``, ``58- b-Catenin1`` and
  ``69MetSI``, and the dedicated pLN slides' panel is ``55 markers``.
* **Xenium** -- ``5kOnly`` for the predesigned 5K panel alone, or ``5k+<custom genes>(<custom
  panel code>)`` for an add-on panel, e.g. ``5k+100(6KTRAF)``.

The facts each rule needs -- which markers are not counted, which single markers carry the sign,
which markers are the metabolic spike-in, what a predesigned panel is called -- are data in the
platform profiles' ``panel_names``. This module only applies them to what an image or a Xenium
record says, and refuses a panel the rule does not cover rather than guess its name.
"""

from __future__ import annotations

import re
from collections.abc import Iterable


class PanelNameError(ValueError):
    """A recorded panel that the naming rule does not cover."""


def marker_name(name: str) -> str:
    """One spelling of a marker name: an image writes ``Bcl-2``, the QuPath manifest ``Bcl_2``."""
    return re.sub(r"[\s\-]+", "_", name.strip())


def phenocycler(markers: Iterable[str], rule: dict) -> str:
    recorded = {marker_name(m): m for m in markers}       # one spelling -> as the image writes it
    names = set(recorded)
    counted = len(names - {marker_name(m) for m in rule["not_counted"]})
    spike = {marker_name(m) for m in rule["metabolic_spike_in"]}
    if spike & names:
        if not spike <= names:
            raise PanelNameError(f"part of the metabolic spike-in: {sorted(spike & names)}")
        return f"{counted}{rule['metabolic_suffix']}"
    signs = [(sign, recorded[marker_name(marker)]) for sign, marker in rule["single_marker"].items()
             if marker_name(marker) in names]
    if len(signs) > 1:
        raise PanelNameError(f"both single markers without the metabolic spike-in: {signs}")
    if signs:
        sign, marker = signs[0]
        return f"{counted}{sign} {marker}"
    return f"{counted} {rule['count_word']}"


def xenium(record: dict, rule: dict) -> str:
    """A Xenium panel from its ``experiment.xenium`` fields."""
    predesigned = record.get("panel_num_targets_predesigned")
    base = rule["predesigned"].get(str(predesigned))
    if base is None:
        raise PanelNameError(f"no name for a predesigned panel of {predesigned} targets")
    custom = record.get("panel_num_targets_custom") or 0
    if not custom:
        return f"{base}{rule['alone_suffix']}"
    design = record.get("panel_design_id")
    if not design:
        raise PanelNameError("an add-on panel without a design id")
    return f"{base}+{custom}({design})"
