"""A donor has one pancreas and 1-4 lymph nodes; every one is a separate section.

Three things have to hold for that cohort, and each of them used to fail quietly rather than
loudly, which is why they are pinned here:

* the ROI enumeration reaches the 4th node (it was hardcoded to three);
* a section outside the enumeration is *reported*, not skipped;
* the qptiff a section registers against is that section's, not its donor's first.
"""

from __future__ import annotations

import pytest
from phenocycler.integration.tissues import DEFAULT_MAX_PLN, rois

from phenocycler import load_config

pytest.importorskip("pandas")


# --------------------------------------------------------------------------- #
# ROI enumeration
# --------------------------------------------------------------------------- #

def test_the_enumeration_reaches_the_fourth_node():
    cfg = load_config()
    pln = cfg.rois_for_tissue("pancreatic_lymph_node")
    assert pln == ["pln_1", "pln_2", "pln_3", "pln_4"]
    assert "pln_4" in cfg.resolve_rois()          # a bare run must cover it


def test_the_node_count_is_configuration_not_code():
    """It describes the cohort — how many nodes were cut — so it must be settable."""
    assert load_config(max_pln=2).rois_for_tissue("pancreatic_lymph_node") == ["pln_1", "pln_2"]
    assert len(load_config(max_pln=6).rois_for_tissue("pancreatic_lymph_node")) == 6
    assert load_config(max_pln=0).rois_for_tissue("pancreatic_lymph_node") == ["pln_1"]


def test_one_definition_not_two():
    """`config.rois_for_tissue` and `tissues.rois` disagreeing is how pln_4 got lost."""
    cfg = load_config()
    assert cfg.rois_for_tissue("pancreatic_lymph_node") == list(
        rois("pancreatic_lymph_node", max_pln=cfg.max_pln))
    assert cfg.rois_for_tissue("pancreas") == list(rois("pancreas"))
    assert DEFAULT_MAX_PLN == 4


def test_pancreas_is_unaffected_and_spleen_uses_its_profile():
    cfg = load_config()
    assert cfg.rois_for_tissue("pancreas") == ["panc"]
    assert cfg.rois_for_tissue("spleen") == ["spl"]


def test_a_node_beyond_the_enumeration_is_still_parseable():
    """The cap bounds iteration, not the vocabulary — so raising it is all that is needed."""
    from phenocycler.sections import key_for, parse

    assert parse("6539pln5").roi == "pln_5"
    assert key_for("6539", "pln_5") == "6539pln5"


# --------------------------------------------------------------------------- #
# The section's own qptiff
# --------------------------------------------------------------------------- #

def _images(tmp_path):
    for name in ("6539_Scan1.er.qptiff", "6539pLN_Scan1.er.qptiff",
                 "6539pLN2_Scan1.er.qptiff"):
        (tmp_path / name).write_bytes(b"")
    return tmp_path


def test_each_section_registers_against_its_own_image(tmp_path):
    """The bug: `{donor}*` glob + `hits[0]`, and `_` sorts before `p`.

    Every lymph node was handed `6539_Scan1` — the pancreas. For a lymph node that is the
    only registration track there is, because the point-set track needs islets to seed on.
    """
    from phenocycler.integration.rasterize import find_qptiff

    d = _images(tmp_path)
    assert find_qptiff(d, "6539", roi="panc").name == "6539_Scan1.er.qptiff"
    assert find_qptiff(d, "6539", roi="pln_1").name == "6539pLN_Scan1.er.qptiff"
    assert find_qptiff(d, "6539", roi="pln_2").name == "6539pLN2_Scan1.er.qptiff"


def test_a_missing_section_image_returns_nothing_rather_than_the_wrong_one(tmp_path):
    """Silence beats a confident misregistration against a different piece of tissue."""
    from phenocycler.integration.rasterize import find_qptiff

    d = _images(tmp_path)
    assert find_qptiff(d, "6539", roi="pln_3") is None      # no 6539pLN3 image exists


def test_the_donor_only_form_still_works_for_callers_without_an_roi(tmp_path):
    from phenocycler.integration.rasterize import find_qptiff

    d = _images(tmp_path)
    assert find_qptiff(d, "6539") is not None
    assert find_qptiff(d, "9999", roi="panc") is None


def test_an_roi_that_names_no_section_does_not_fall_back_to_the_donor_glob(tmp_path):
    """Once an ROI is given it is honoured strictly.

    An earlier version caught the `SectionError`, blanked the wanted key and fell through to
    the `{donor}*` glob — which is the very "take the donor's first image" behaviour the roi
    parameter exists to prevent. An unresolvable ROI must yield nothing, not the pancreas.
    """
    from phenocycler.integration.rasterize import find_qptiff

    d = _images(tmp_path)
    assert find_qptiff(d, "6539", roi="not_a_region") is None
    assert find_qptiff(d, "6539", roi="spleen_2") is None
