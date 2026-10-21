"""Stage 00 -- the source inventory -- on a synthetic project: what it lists, and what blocks its
gate (MultiPlex5000 D3).

    .venv/bin/python -m pytest core/multiplex5000/tests/current/test_stage_inventory.py -q
"""

from __future__ import annotations

import csv
import json

import pytest

from multiplex5000 import gates, stages, tasks
from multiplex5000 import project as P
from multiplex5000.relationships import RelationshipError

STAGE = stages.get("00_inventory")


def _run(path):
    project = P.load(path)
    result = tasks.run_stage(STAGE, project)
    folder = result.path
    checks = {c["name"]: c for c in gates.load(folder)["checks"]}
    return project, folder, checks


def _rows(path):
    with open(path, newline="") as handle:
        return list(csv.DictReader(handle))


def test_the_inventory_lists_every_acquisition_by_its_recorded_name(synthetic_project):
    project, folder, checks = _run(synthetic_project())
    acquisitions = {r["acquisition_id"]: r for r in _rows(folder / "acquisitions.csv")}
    assert set(acquisitions) == {
        "phenocycler:6374_Scan1 - resolution #1", "phenocycler:HDL112_Scan3 - resolution #1",
        "phenocycler:6374pLN_Scan1.qptiff", "xenium:0041000:Pan", "xenium:0041001:Pancreas",
        "xenium:0041001:pLN", "if:if_6374.qptiff"}
    assert acquisitions["xenium:0041001:pLN"]["donor_id"] == "0112", "112 is resolved by alias"
    assert acquisitions["xenium:0041001:pLN"]["tissue"] == "pancreatic_lymph_node"
    assert acquisitions["phenocycler:6374_Scan1 - resolution #1"]["size_matches_record"] == "True"
    coverage = {r["donor_id"]: r for r in _rows(folder / "coverage.csv")}
    assert coverage["0112"]["Xenium pLN"] == "1" and coverage["6374"]["post-Xenium IF"] == "1"
    kinds = {r["kind"] for r in _rows(folder / "relationships.csv")}
    assert {"same_section", "serial_section", "different_block", "shared_slide"} <= kinds
    assert all(c["passed"] for c in checks.values() if c["blocking"])
    assert (folder / "gate" / "coverage.png").stat().st_size > 1000


DESIGNATED = ("6374: pancTail and 1 pLN same-block; bulk T and B cell aggregates kept in Tissue, "
              "sequential Xenium\n"
              "HDL112: pancTail and 1 pLN different-block, area of pLN is about the same as panc, "
              "same-slide post-phenocycler Xenium\n")
BLOCKS = ("image,entry_id,donor,pln_annotation,object_id,pieces,area_mm2,overlap_with_tissue_um2,"
          "gap_to_tissue_um,measured_relation,block_from_notes,notes\n"
          "img112,1,HDL112,1,uuid-1,1,1.0,0.0,664.7,separate,different-block,x\n"
          "img6374,2,6374,1,uuid-2,1,1.0,0.0,0.0,abuts Tissue (gap 0),same-block,x\n")


def test_coverage_carries_the_notes_designations_beside_each_pancreas(synthetic_project):
    # the notes as the maintainer keeps them since 2026-10-21: no header, a blank first line
    _, folder, checks = _run(synthetic_project(
        notes__txt="\n" + DESIGNATED, notes_copy__md="IO60panc2nd Notes\n" + DESIGNATED,
        pln_blocks__csv=BLOCKS))
    coverage = {r["donor_id"]: r for r in _rows(folder / "coverage.csv")}
    assert list(coverage["6374"])[3:6] == ["PhenoCycler pancreas", STAGE.SAME_BLOCK,
                                           STAGE.DIFFERENT_BLOCK]
    assert (coverage["6374"][STAGE.SAME_BLOCK], coverage["6374"][STAGE.DIFFERENT_BLOCK]) == ("1", "0")
    assert (coverage["0112"][STAGE.SAME_BLOCK], coverage["0112"][STAGE.DIFFERENT_BLOCK]) == ("0", "1")
    assert coverage["6374"][STAGE.SAME_BLOCK + STAGE.ANNOTATED] == "1"
    assert coverage["6374"]["caution in the notes"] == "True"
    designated = checks["pLN/lymphoid the donor notes designate, by pancreas"]
    assert designated["value"] == ["0112: 1 different-block", "6374: 1 same-block"]
    quoted = checks["cautions in the donor notes, quoted"]
    assert quoted["value"] == ["0112: area of pLN is about the same as panc",
                               "6374: bulk T and B cell aggregates kept in Tissue"]
    assert not designated["blocking"] and not quoted["blocking"]
    assert all(c["passed"] for c in checks.values() if c["blocking"])


def test_the_figure_shows_a_difference_and_what_the_notes_leave_unsaid():
    row = {STAGE.SAME_BLOCK: 3, STAGE.SAME_BLOCK + STAGE.ANNOTATED: 0,
           STAGE.DIFFERENT_BLOCK: "", STAGE.DIFFERENT_BLOCK + STAGE.ANNOTATED: 0,
           "Xenium pLN": 2}
    assert STAGE._cell(row, STAGE.SAME_BLOCK) == (3, "3 / 0", True)       # notes / pln_blocks
    assert STAGE._cell(row, STAGE.DIFFERENT_BLOCK) == (0, "?", False)     # the line does not say
    assert STAGE._cell(row, "Xenium pLN") == (2, "2", False)


def test_the_receipt_binds_every_source_it_read(synthetic_project):
    _, folder, _ = _run(synthetic_project())
    receipt = json.loads((folder / "receipt.json").read_text())
    names = {i["name"] for i in receipt["inputs"]}
    assert {"cohort_spec", "qupath_manifest", "xenium_paths", "section_modes", "donor_notes",
            "pln_blocks", "pln_node_annotations", "restain_sections",
            "restain_registration_inputs", "donor_notes_copy"} == names
    assert all(i.get("sha256") for i in receipt["inputs"])


@pytest.mark.parametrize("replace, check", [
    ({"xenium_paths__csv": "donor_id,roi,batch,source_path\n"
      "6374,panc,B1,bundles/output-XETG00298__0041000__Pan__20240101__000000\n"
      "6374,panc,B1,bundles/output-XETG00298__0041000__Pan__20240101__000000\n"},
     "acquisition ids are unique"),
    ({"notes__txt": "IO60panc2nd Notes\n6374: pancTail and 1 pLN same-block, sequential Xenium\n"
      "HDL112: pancTail and 1 pLN different-block, same-slide post-phenocycler Xenium\n",
      "notes_copy__md": "IO60panc2nd Notes\n6374: pancTail and 1 pLN same-block, sequential Xenium\n"
      "HDL112: pancTail and 1 pLN different-block, same-slide post-phenocycler Xenium\n"},
     "every disagreement between the recorded sources is known"),
    ({"notes__txt": "IO60panc2nd Notes\n6374: pancTail no pLN, sequential Xenium\n",
      "notes_copy__md": "IO60panc2nd Notes\n6374: pancTail no pLN, sequential Xenium\n"},
     "every cohort donor has a line in the donor notes"),
])
def test_what_blocks_the_gate(synthetic_project, replace, check):
    project, folder, checks = _run(synthetic_project(**replace))
    assert not checks[check]["passed"] and checks[check]["blocking"]
    with pytest.raises(gates.GateError, match="blocking"):
        gates.accept(folder, by="tester")


def test_a_raw_image_that_is_not_the_size_its_record_says_blocks(synthetic_project, tmp_path):
    path = synthetic_project()
    (tmp_path / "raw" / "6374.qptiff").write_bytes(b"x" * 99)
    _, _, checks = _run(path)
    check = checks["every raw image is the size its record says"]
    assert not check["passed"] and check["value"] == ["phenocycler:6374_Scan1 - resolution #1"]


def test_the_stage_needs_its_sources(synthetic_project, monkeypatch):
    project = P.load(synthetic_project())
    stripped = P.Project(project.path, project.id, project.release, project.cohort,
                         project.stages, project.automatic_gates, project.releases_root,
                         project.settings, {k: v for k, v in project.sources.items()
                                            if k != "pln_blocks"})
    monkeypatch.setattr(P.roots, "source", lambda name: None)
    with pytest.raises(RelationshipError, match="pln_blocks"):
        STAGE.inputs(stripped)
