"""How a recorded panel is written in every output (maintainer, 2026-10-21).

* **PhenoCycler** -- the number of markers the image records, DAPI not counted; then ``MetSI``
  where the whole metabolic spike-in is in the panel; then every other marker the image adds to
  its tissue's reference panel with ``+`` and every reference marker it lacks with ``-`` (maintainer,
  2026-10-24: "each marker needs a + or -"), added before lacking, each alphabetical, as the image
  or the reference spells it. A panel that is its reference is the number and the word ``markers``.
  So the IO60 pancreas panels are ``57 markers``, ``58 +CD38``, ``58 +b-Catenin1`` and ``69MetSI
  +b-Catenin1 +CD141 +CD21 +CD38 +Histone H3 (p Ser28) -IAPP -Ker8-18``, and the dedicated pLN
  slides' panel is ``55 markers``.
* **Xenium** -- ``5kOnly`` for the predesigned 5K panel alone, or ``5k+<custom genes>(<custom
  panel code>)`` for an add-on panel, e.g. ``5k+100(6KTRAF)``.
* **post-Xenium IF** -- one panel for every image (maintainer, 2026-10-21): ``INS, GCG, failed
  marker, AF``. Where an image's channel roles are recorded, they must be that panel.

The facts each rule needs -- which markers are not counted, each tissue's reference panel,
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


def phenocycler(markers: Iterable[str], rule: dict, tissue: str) -> str:
    """A PhenoCycler panel from the markers an image records, against ``tissue``'s reference."""
    not_counted = {marker_name(m) for m in rule["not_counted"]}
    recorded = {marker_name(m): m for m in markers if marker_name(m) not in not_counted}
    if tissue not in rule["reference"]:
        raise PanelNameError(f"no reference panel for {tissue}")
    reference = {marker_name(m): m for m in rule["reference"][tissue]
                 if marker_name(m) not in not_counted}
    head = f"{len(recorded)}"
    spike = {marker_name(m) for m in rule["metabolic_spike_in"]}
    if spike & set(recorded):
        if not spike <= set(recorded):
            raise PanelNameError(f"part of the metabolic spike-in: {sorted(spike & set(recorded))}")
        head += rule["metabolic_suffix"]
    else:
        spike = set()

    def alphabetical(spellings: dict, keys: set) -> list[str]:
        return sorted((spellings[k] for k in keys), key=str.casefold)

    signed = ([f"+{m}" for m in alphabetical(recorded, set(recorded) - set(reference) - spike)]
              + [f"-{m}" for m in alphabetical(reference, set(reference) - set(recorded))])
    if not signed and not spike:
        return f"{head} {rule['count_word']}"
    return " ".join([head, *signed])


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


def post_xenium_if(rule: dict, roles: Iterable[str] | None = None) -> str:
    """The post-Xenium IF panel's name. ``roles`` are the channel roles recorded for one image,
    where any are; a recorded panel that is not the declared one is refused, never renamed."""
    channels = list(rule["channels"])
    if roles is not None:
        recorded = sorted(set(roles) - set(rule["not_named"]))
        if recorded != sorted(channels):
            raise PanelNameError(f"the channel roles recorded are {recorded}, not the panel every "
                                 f"post-Xenium IF image carries, {sorted(channels)}")
    return ", ".join(rule["written"].get(c, c) for c in channels)
