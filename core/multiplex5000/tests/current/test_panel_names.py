"""Panel names, the maintainer's way (2026-10-21), from the rules in the platform profiles.

    .venv/bin/python -m pytest core/multiplex5000/tests/current/test_panel_names.py -q

The marker sets are the five the IO60 qptiffs record, written as the images write them; the
expected names are the ones the maintainer chose: 57 markers, 58+ CD38, 58- b-Catenin1, 69MetSI
and 55 markers for the dedicated pLN slides; and 5kOnly or 5k+100(<custom panel code>) for Xenium.
"""

from __future__ import annotations

import pytest

from multiplex5000 import acquisition_records, panel_names, profiles

RULES = {p: profiles.load_all()[("platform", p)].panel_names
         for p in ("phenocycler", "xenium", "post_xenium_if")}

#: The 57 markers every IO60 pancreas panel shares, and DAPI.
BASE = ["DAPI"] + [f"M{i}" for i in range(55)] + ["IAPP", "Ker8-18"]
METABOLIC = ["ASCT2", "ATP5A", "CS", "G6PD", "GLUT1", "HK1", "IDH2", "LDHA", "SDHA"]


@pytest.mark.parametrize("markers, name", [
    (BASE, "57 markers"),
    (BASE + ["CD38"], "58+ CD38"),
    (BASE + ["b-Catenin1"], "58- b-Catenin1"),          # as the image spells it
    ([m for m in BASE if m not in ("IAPP", "Ker8-18")] + METABOLIC
     + ["CD38", "b_Catenin1", "CD141", "CD21", "Histone H3 (p Ser28)"], "69MetSI"),
    (["DAPI"] + [f"N{i}" for i in range(55)], "55 markers"),
])
def test_phenocycler_panels_are_named_the_maintainers_way(markers, name):
    assert panel_names.phenocycler(markers, RULES["phenocycler"]) == name


@pytest.mark.parametrize("markers, message", [
    (BASE + METABOLIC[:4], "part of the metabolic spike-in"),
    (BASE + ["CD38", "b_Catenin1"], "both single markers"),
])
def test_a_panel_the_rule_does_not_cover_is_refused_not_guessed(markers, message):
    with pytest.raises(panel_names.PanelNameError, match=message):
        panel_names.phenocycler(markers, RULES["phenocycler"])


@pytest.mark.parametrize("record, name", [
    ({"panel_num_targets_predesigned": 5001, "panel_num_targets_custom": 100,
      "panel_design_id": "6KTRAF"}, "5k+100(6KTRAF)"),
    ({"panel_num_targets_predesigned": 5001, "panel_num_targets_custom": 100,
      "panel_design_id": "Y9QRJ6"}, "5k+100(Y9QRJ6)"),
    ({"panel_num_targets_predesigned": 5001, "panel_num_targets_custom": 0,
      "panel_design_id": "hAtlas_v1.1"}, "5kOnly"),
])
def test_xenium_panels_are_named_the_maintainers_way(record, name):
    assert panel_names.xenium(record, RULES["xenium"]) == name


def test_an_unknown_predesigned_xenium_panel_is_refused():
    with pytest.raises(panel_names.PanelNameError, match="480"):
        panel_names.xenium({"panel_num_targets_predesigned": 480}, RULES["xenium"])


def _header(path, experiment, date, scan, position):
    return acquisition_records.QptiffHeader(path, "s", experiment, position, date, scan, "f",
                                            ("DAPI",))


def test_slides_scanned_under_one_experiment_on_one_day_are_one_run():
    headers = [_header("a", "Exp", "d1", "2026-04-01 14:37:54", "Slide1"),
               _header("b", "Exp", "d1", "2026-04-01 14:57:03", "Slide2"),
               _header("c", "Exp", "d1", "2026-04-04 15:04:12", "Slide2"),     # alone that day
               _header("d", "Other", "d2", "2026-04-01 15:10:00", "Slide1")]   # another experiment
    runs = acquisition_records.phenocycler_runs(headers)
    assert [[h.path for h in group] for group in runs.values()] == [["a", "b"], ["c"], ["d"]]


def test_every_post_xenium_if_image_has_the_one_panel_and_another_is_refused():
    rule = RULES["post_xenium_if"]
    assert panel_names.post_xenium_if(rule) == "INS, GCG, failed marker, AF"
    # as restain_sections.csv records the channel roles: DAPI's is not part of the name
    assert panel_names.post_xenium_if(rule, ["NUCLEUS", "GCG", "FAILED", "INS", "AF"]) == (
        "INS, GCG, failed marker, AF")
    with pytest.raises(panel_names.PanelNameError, match="SST"):
        panel_names.post_xenium_if(rule, ["NUCLEUS", "GCG", "SST", "INS", "AF"])
