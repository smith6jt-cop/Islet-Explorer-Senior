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
    _, folder, checks = _run(synthetic_project(
        notes__txt="IO60panc2nd Notes\n" + DESIGNATED,
        notes_copy__md="IO60panc2nd Notes\n" + DESIGNATED, pln_blocks__csv=BLOCKS))
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
    inputs = {i["name"]: i for i in receipt["inputs"]}
    tables = {"cohort_spec", "qupath_manifest", "xenium_paths", "section_modes", "donor_notes",
              "pln_blocks", "pln_node_annotations", "restain_sections",
              "restain_registration_inputs", "donor_notes_copy", "xenium_records_manifest"}
    qptiffs = {"qptiff:6374.qptiff", "qptiff:0112.qptiff", "qptiff:6374pLN_Scan1.qptiff"}
    outputs = ("output-XETG00298__0041000__Pan__20240101__000000",
               "output-XETG00298__0041001__Pancreas__20240101__000000",
               "output-XETG00298__0041001__pLN__20240101__000000")
    records = {f"xenium:{o}:{f}" for o in outputs for f in STAGE.XENIUM_RECORDS}
    assert set(inputs) == tables | qptiffs | records
    assert all(inputs[n].get("sha256") for n in tables | records), "tables and records are hashed"
    # a qptiff is bound by size and time, with the digest the QuPath manifest records for it
    assert all(inputs[n]["bytes"] and not inputs[n].get("sha256") for n in qptiffs)
    assert inputs["qptiff:6374.qptiff"]["recorded_sha256"] == "0" * 64
    # the output held here only in part is read from the copy of the complete output's records
    held = inputs[f"xenium:{outputs[2]}:experiment.xenium"]["path"]
    assert "xenium_records" in held and "bundles" in inputs[f"xenium:{outputs[0]}:experiment.xenium"]["path"]


def test_panels_runs_and_partners_come_from_the_instrument_records(synthetic_project):
    _, folder, checks = _run(synthetic_project())
    acquisitions = {r["acquisition_id"]: r for r in _rows(folder / "acquisitions.csv")}
    pancreas = acquisitions["phenocycler:6374_Scan1 - resolution #1"]
    assert (pancreas["panel"], pancreas["run_date"], pancreas["run_position"]) == ("2+", "2026-01-02",
                                                                                   "Slide1")
    assert pancreas["run_partners"] == "0112 pancreas (Slide2)"
    assert pancreas["markers_match_manifest_panel"] == "True"
    pln = acquisitions["phenocycler:6374pLN_Scan1.qptiff"]
    assert (pln["panel"], pln["run_partners"]) == ("1", "run alone")
    xenium = acquisitions["xenium:0041001:pLN"]
    assert (xenium["panel"], xenium["run"], xenium["run_date"]) == ("5k+100(6KTRAF)",
                                                                    "xenium:R1_synthetic",
                                                                    "2024-01-01")
    assert xenium["local_copy"] == "partial" and xenium["records_read_from"].startswith("copy of hpg:")
    assert acquisitions["xenium:0041000:Pan"]["run_partners"] == "0112 pancreas + pLN"
    runs = {r["run"]: r for r in _rows(folder / "runs.csv")}
    assert set(runs) == {"phenocycler:Synthetic_Pancreas:2026-01-02",
                         "phenocycler:Synthetic_pLN:2026-01-09", "xenium:R1_synthetic"}
    assert runs["phenocycler:Synthetic_Pancreas:2026-01-02"]["slides"] == "2"
    coverage = {r["donor_id"]: r for r in _rows(folder / "coverage.csv")}
    assert "donor notes" not in coverage["6374"]
    assert (coverage["0112"]["PhenoCycler pancreas panel"], coverage["0112"]["Xenium pLN panel"]) == (
        "2+", "5k+100(6KTRAF)")
    assert checks["PhenoCycler slides run alone"]["value"] == ["6374 pLN, 2026-01-09"]
    assert checks["panels, as recorded and named"]["value"] == {
        "phenocycler": {"1": 1, "2+": 2}, "xenium": {"5k+100(6KTRAF)": 3}}
    assert all(c["passed"] for c in checks.values() if c["blocking"])
    assert (folder / "gate" / "runs.png").stat().st_size > 1000


def _edit_json(path, **changes):
    data = json.loads(path.read_text())
    data.update(changes)
    path.write_text(json.dumps(data))


def test_a_run_of_more_than_two_slides_blocks_and_names_a_partner_outside_the_inventory(
        synthetic_project, tmp_path):
    import datetime

    from conftest import PANCREAS_MARKERS, write_qptiff

    path = synthetic_project()
    write_qptiff(tmp_path / "raw" / "9999_Scan1.qptiff", slide="9999",
                 experiment="Synthetic_Pancreas", position="Slide3",
                 experiment_date="2025-12-01T09:00:00Z",
                 scanned=datetime.datetime(2026, 1, 2, 14, 40), markers=PANCREAS_MARKERS)
    _, _, checks = _run(path)
    assert checks["no run holds more than two slides"]["value"] == [
        "phenocycler:Synthetic_Pancreas:2026-01-02: 3 slides"]
    assert not checks["no run holds more than two slides"]["passed"]
    assert checks["run partners that are no acquisition here"]["value"][0].startswith(
        "9999 (9999_Scan1.qptiff, no acquisition here)")


@pytest.mark.parametrize("edit, check", [
    (lambda root: (root / "xenium_records" / "output-XETG00298__0041001__pLN__20240101__000000"
                   / "experiment.xenium").unlink(),
     "every listed Xenium output has its three instrument records"),
    (lambda root: (root / "xenium" / "bundles" / "output-XETG00298__0041000__Pan__20240101__000000"
                   / "metrics_summary.csv").write_text(
        "run_name,region_name,panel_design_id\nanother_run,Pan,6KTRAF\n"),
     "each Xenium output's records agree with each other and with its slide"),
    (lambda root: _edit_json(root / "sources" / "qupath_manifest.json", panels=[
        {"panel_id": "panel", "channels": [{"marker": "DAPI"}, {"marker": "CD3e"}]}]),
     "every pancreas image records the markers of its QuPath-manifest panel"),
])
def test_what_the_instrument_records_block(synthetic_project, tmp_path, edit, check):
    path = synthetic_project()
    edit(tmp_path)
    _, _, checks = _run(path)
    assert not checks[check]["passed"] and checks[check]["blocking"]


def test_a_renamed_xenium_folder_blocks_until_its_reason_is_known(synthetic_project, tmp_path,
                                                                   monkeypatch):
    path = synthetic_project()
    output = "output-XETG00298__0041000__Pan__20240101__000000"
    _edit_json(tmp_path / "xenium" / "bundles" / output / "experiment.xenium", region_name="Region_1")
    (tmp_path / "xenium" / "bundles" / output / "metrics_summary.csv").write_text(
        "run_name,region_name,panel_design_id\nR1_synthetic,Region_1,6KTRAF\n")
    _, folder, checks = _run(path)
    check = checks["every Xenium output folder names the region its experiment.xenium records"]
    assert not check["passed"] and check["value"] == [
        f"{output}: experiment.xenium names region Region_1"]
    ids = {r["acquisition_id"] for r in _rows(folder / "acquisitions.csv")}
    assert "xenium:0041000:Region_1" in ids, "the id takes the region the instrument recorded"
    monkeypatch.setitem(STAGE.KNOWN_RENAMED_FOLDERS, output, "renamed after export")
    tasks.run_stage(STAGE, P.load(path), force=True)
    checks = {c["name"]: c for c in gates.load(folder)["checks"]}
    assert checks["every Xenium output folder names the region its experiment.xenium records"]["passed"]
    assert checks["Xenium output folders renamed, with their known reason"]["value"] == [
        f"{output} (renamed after export)"]


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
