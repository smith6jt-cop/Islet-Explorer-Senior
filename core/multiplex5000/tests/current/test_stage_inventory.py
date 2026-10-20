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
