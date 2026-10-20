"""``section_relationship`` -- a filtering and sorting attribute, derived and never declared.

Some donors carry pancreas and pancreatic lymph node in one QuPath image; most do not.  That is
worth filtering on, and it is read off the export's own ``Parent`` column rather than asserted:
the per-image annotation census ``rebuild_donor_from_revised_export.py`` writes is the evidence,
and where there is none the column says ``unavailable(...)`` and names which evidence was absent.

The attribute licenses no distance.  Islet-to-lymphoid distance is not wanted, and two tissues in
one image share a pixel grid and no anatomy; ``refuse_cross_tissue_distance()`` is untouched by
anything here, and one test holds that builder's guard in place so this column cannot become an
excuse to weaken it.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parents[4]
sys.path.insert(0, str(REPO))

from scripts.senior import build_sections_layer as sections  # noqa: E402

BUNDLE = "HuPrime5kPanTissue/6KTRAF/MM035/output-XETG00298__{serial}__{label}__20240919__203848"


@pytest.fixture()
def inventory(tmp_path):
    """Four donors: one same-block, one same-slide, one both, one with neither."""
    directory = tmp_path / "integration"
    directory.mkdir()
    (directory / "xenium_paths.csv").write_text(
        "donor_id,roi,source_path\n"
        f"115,panc,{BUNDLE.format(serial='0041117', label='PanT')}\n"
        f"115,pln_1,{BUNDLE.format(serial='0041117', label='pLN')}\n"
        f"6533,panc,{BUNDLE.format(serial='0041393', label='Pancreas')}\n"
        f"6533,pln_1,{BUNDLE.format(serial='0041393', label='pLN')}\n"
        f"6539,panc,{BUNDLE.format(serial='0059865', label='Panc')}\n"
        f"6374,panc,{BUNDLE.format(serial='0041032', label='Pan')}\n"
    )
    (directory / "restain_sections.csv").write_text(
        "image_file,donor_id,roi,xenium_serial,pixel_size_um\n"
        "HDL115_Scan1.qptiff,115,panc,0041117,0.2505\n"
        "0059865_Scan1.qptiff,6539,panc,0059865,0.2505\n"
    )
    return directory


@pytest.fixture()
def evidence(tmp_path):
    path = tmp_path / "section_parents.json"
    path.write_text(json.dumps({"donors": {
        "6533": {
            "image": "0080870_6533_pLNpanc_Scan3.er.qptiff - resolution #1",
            "annotation_parent_counts": {"Tissue": 931910, "pLN": 34841},
            "pln_annotation_classes": ["pLN"],
        },
        "115": {
            "image": "HDL115pancLN_Scan1.er.qptiff - resolution #1",
            "annotation_parent_counts": {"Tissue": 656835, "pLN": 261761},
            "pln_annotation_classes": ["pLN"],
        },
        "6374": {
            "image": "6374_Scan1.er.qptiff - resolution #1",
            "annotation_parent_counts": {"Tissue": 100},
            "pln_annotation_classes": [],
        },
    }}))
    return path


def rows_of(table):
    return {(r["donor_id"], r["roi"], r["modality"]): r for r in table.to_pylist()}


def test_the_schema_carries_the_attribute_and_its_basis(inventory, evidence):
    table = sections.build_sections(inventory, evidence)
    assert table.schema == sections.SECTIONS_SCHEMA
    assert "section_relationship" in table.schema.names
    assert "section_relationship_basis" in table.schema.names


def test_same_block_comes_from_the_annotation_census(inventory, evidence):
    """6533's pancreas image carries a pLN annotation beside the tissue annotation."""
    rows = rows_of(sections.build_sections(inventory, evidence))
    row = rows[("6533", "panc", "xenium")]
    assert row["section_relationship"] == "same_block"
    assert "pln_annotation_in_pancreas_image" in row["section_relationship_basis"]
    assert "pLN:34841" in row["section_relationship_basis"]


def test_same_slide_comes_from_the_recorded_restain_serial(inventory, evidence):
    """The post-Xenium IF and the Xenium row are one physical section, and the serial says so."""
    rows = rows_of(sections.build_sections(inventory, evidence))
    for key in (("6539", "panc", "restain_if"), ("6539", "panc", "xenium")):
        assert rows[key]["section_relationship"] == "same_slide"
        assert "0059865" in rows[key]["section_relationship_basis"]


def test_a_row_with_both_keeps_both_in_the_basis(inventory, evidence):
    """Precedence picks one value for the filter; nothing measured is dropped."""
    rows = rows_of(sections.build_sections(inventory, evidence))
    row = rows[("115", "panc", "restain_if")]
    assert row["section_relationship"] == "same_block"
    assert "same_block:" in row["section_relationship_basis"]
    assert "same_slide:" in row["section_relationship_basis"]


def test_no_evidence_is_a_named_unavailable_not_a_different_block(inventory, evidence):
    """Nothing in this checkout records block identity, so nothing may claim it."""
    rows = rows_of(sections.build_sections(inventory, evidence))
    row = rows[("6374", "panc", "xenium")]
    assert row["section_relationship"] is None
    assert row["section_relationship_basis"] == "unavailable(no_pln_annotation_in_this_donors_image)"

    values = set(sections.build_sections(inventory, evidence)["section_relationship"].to_pylist())
    assert "different_block" not in values


def test_a_donor_absent_from_the_census_says_so(inventory, evidence):
    payload = json.loads(evidence.read_text())
    payload["donors"].pop("6374")
    evidence.write_text(json.dumps(payload))
    rows = rows_of(sections.build_sections(inventory, evidence))
    assert rows[("6374", "panc", "xenium")]["section_relationship_basis"] == (
        "unavailable(donor_absent_from_annotation_census)")


def test_an_absent_census_still_builds_the_layer(inventory, tmp_path):
    """A layer whose relationship is unknown still has to build, and has to say why."""
    table = sections.build_sections(inventory, tmp_path / "nothing_here.json")
    rows = rows_of(table)
    assert rows[("6533", "panc", "xenium")]["section_relationship"] is None
    assert rows[("6533", "panc", "xenium")]["section_relationship_basis"] == (
        "unavailable(annotation_census_not_built)")
    assert rows[("6539", "panc", "restain_if")]["section_relationship"] == "same_slide"


def test_a_malformed_census_is_not_read_as_evidence(inventory, tmp_path):
    broken = tmp_path / "broken.json"
    broken.write_text("{not json")
    assert sections.read_parent_evidence(broken) == {}


def test_the_census_donor_spelling_is_normalised(inventory, tmp_path):
    """The census records the manifest's '0112'; the layer's donor column says '112'."""
    path = tmp_path / "census.json"
    path.write_text(json.dumps({"donors": {"0112": {
        "image": "0081366_HDL112_pLNpanc_Scan3.er.qptiff - resolution #1",
        "annotation_parent_counts": {"Tissue": 641639, "pLN": 617396},
        "pln_annotation_classes": ["pLN"],
    }}}))
    assert set(sections.read_parent_evidence(path)) == {"112"}


def test_the_attribute_does_not_change_which_sections_exist(inventory, evidence):
    """It is an attribute of the rows, not a filter on them."""
    with_evidence = sections.build_sections(inventory, evidence)
    without = sections.build_sections(inventory, None)
    assert with_evidence.num_rows == without.num_rows
    keys = ("donor_id", "roi", "modality", "image_file", "xenium_serial", "pixel_size_um",
            "frame")
    assert with_evidence.select(keys) == without.select(keys)


def test_the_cross_tissue_distance_refusal_is_untouched():
    """This column is for filtering and sorting. It is not a licence to measure across tissues."""
    source = (REPO / "scripts/senior/build_distance_profiles.py").read_text()
    assert "def refuse_cross_tissue_distance" in source
    assert "def _assert_no_cross_tissue_distance" in source
