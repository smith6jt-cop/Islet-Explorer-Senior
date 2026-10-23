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


def _at(slide: str, column: str) -> str:
    """The coverage.csv column of one slide's section or acquisition."""
    return f"{slide}: {column}"


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
    assert coverage["0112"][_at(STAGE.PC_SLIDE, STAGE.XENIUM)] == "5k+100(6KTRAF); 2024-01-01 with 6374"
    assert coverage["6374"][_at(STAGE.XENIUM_SLIDE, STAGE.POST_XENIUM_IF)] == (
        "INS, GCG, failed marker, AF"), "no channel roles recorded for it: the one IF panel"
    # disease status from the project's metadata (112 as 0112, by its declared alias)
    assert [(coverage[d]["disease_status"], coverage[d]["disease_status_source"])
            for d in ("6374", "0112")] == [("T1D", "donor_metadata"), ("ND", "donor_metadata")]
    kinds = {r["kind"] for r in _rows(folder / "relationships.csv")}
    assert {"same_section", "serial_section", "different_block", "shared_slide"} <= kinds
    assert all(c["passed"] for c in checks.values() if c["blocking"])
    assert (folder / "gate" / "coverage.png").stat().st_size > 1000
    assert (folder / "gate" / "xenium_slides.png").stat().st_size > 1000


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
    same, different = (_at(STAGE.PC_SLIDE, c) for c in (STAGE.SAME_BLOCK, STAGE.DIFFERENT_BLOCK))
    columns = list(coverage["6374"])
    first = columns.index(_at(STAGE.PC_SLIDE, "pancreas"))
    assert columns[first:first + 3] == [_at(STAGE.PC_SLIDE, "pancreas"), same, different]
    assert (coverage["6374"][same], coverage["6374"][different]) == ("1 pLN", "")
    assert (coverage["0112"][same], coverage["0112"][different]) == ("", "1 pLN")
    assert (coverage["6374"]["design"], coverage["0112"]["design"]) == ("sequential", "same slide")
    assert coverage["6374"][STAGE.SAME_BLOCK + STAGE.ANNOTATED] == "1"
    assert coverage["6374"]["caution in the notes"] == "True"
    designated = checks["pLN/lymphoid the donor notes designate, by pancreas"]
    assert designated["value"] == ["0112: 1 different-block", "6374: 1 same-block"]
    quoted = checks["cautions in the donor notes, quoted"]
    assert quoted["value"] == ["0112: area of pLN is about the same as panc",
                               "6374: bulk T and B cell aggregates kept in Tissue"]
    assert not designated["blocking"] and not quoted["blocking"]
    assert all(c["passed"] for c in checks.values() if c["blocking"])


def _cell(text="", kind=None, damaged=False, missing=False):
    return dict(text=text, kind=kind, damaged=damaged, missing=missing)


def test_a_coverage_cell_says_what_it_holds_never_a_count_and_marks_only_its_damage():
    text = STAGE._cell_text
    assert text(_cell("3", STAGE.LYMPH_NODE), STAGE.DIFFERENT_BLOCK, "section") == "3 pLN"
    assert text(_cell("2", STAGE.LYMPH_NODE), STAGE.PLN_SECTION, "section") == "2 pLN"
    assert text(None, STAGE.SAME_BLOCK, "section") == ""                    # the line does not say
    assert text(_cell("", STAGE.PANCREAS, damaged=True), "pancreas", "section") == (
        "pancreas (damaged)")
    assert text(_cell("5kOnly\n2026-09-10 with 0112", STAGE.XENIUM), STAGE.XENIUM,
                "acquisition") == "5kOnly; 2026-09-10 with 0112"
    assert text(_cell(STAGE.NOT_LISTED, missing=True), STAGE.XENIUM, "acquisition") == STAGE.NOT_LISTED
    assert STAGE._run_text({"run_date": "2025-09-10", "run_with": ""}) == "2025-09-10 alone"
    assert STAGE._run_text({"run_date": "2025-09-11", "run_with": "6579"}) == "2025-09-11 with 6579"


def test_the_receipt_binds_every_source_it_read(synthetic_project):
    _, folder, _ = _run(synthetic_project())
    receipt = json.loads((folder / "receipt.json").read_text())
    inputs = {i["name"]: i for i in receipt["inputs"]}
    tables = {"cohort_spec", "qupath_manifest", "xenium_paths", "section_modes", "donor_notes",
              "pln_blocks", "pln_node_annotations", "restain_sections",
              "restain_registration_inputs", "donor_notes_copy", "xenium_records_manifest",
              "tissue_condition", "xenium_lymphoid_inspection", "donor_metadata", "npod_cases"}
    qptiffs = {"qptiff:6374.qptiff", "qptiff:0112.qptiff", "qptiff:6374pLN_Scan1.qptiff"}
    outputs = ("output-XETG00298__0041000__Pan__20240101__000000",
               "output-XETG00298__0041001__Pancreas__20240101__000000",
               "output-XETG00298__0041001__pLN__20240101__000000")
    records = {f"xenium:{o}:{f}" for o in outputs for f in STAGE.XENIUM_RECORDS}
    overviews = {"overview:0041000", "overview:0041001"}       # one per Xenium slide
    fovs = {f"fovs:{o}" for o in outputs}                      # one per region output
    annotations = {"annotations:6374_Scan1 - resolution #1",  # one per QuPath pancreas image
                   "annotations:HDL112_Scan3 - resolution #1"}
    assert set(inputs) == tables | qptiffs | records | overviews | fovs | annotations
    assert all(inputs[n].get("sha256") for n in tables | records | overviews | fovs | annotations), (
        "tables, records, overview scans, field-of-view locations and annotations are hashed")
    assert inputs["annotations:6374_Scan1 - resolution #1"]["path"].endswith(
        "6374_Scan1_-_resolution_1.geojson"), "the file name QuPath's export gives the image"
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
    assert (pancreas["panel"], pancreas["run_date"], pancreas["run_position"]) == (
        "2+ CD38", "2026-01-02", "Slide1")
    # a PhenoCycler experiment's slide position is no slide number; 0112's slide is Xenium 0041001
    assert pancreas["run_partners"] == "0112 pancreas (slide 0041001)"
    assert pancreas["markers_match_manifest_panel"] == "True"
    pln = acquisitions["phenocycler:6374pLN_Scan1.qptiff"]
    assert (pln["panel"], pln["run_partners"]) == ("2 markers", "run alone")
    xenium = acquisitions["xenium:0041001:pLN"]
    assert (xenium["panel"], xenium["run"], xenium["run_date"]) == ("5k+100(6KTRAF)",
                                                                    "xenium:R1_synthetic",
                                                                    "2024-01-01")
    assert xenium["local_copy"] == "partial" and xenium["records_read_from"].startswith("copy of hpg:")
    assert acquisitions["xenium:0041000:Pan"]["run_partners"] == (
        "0112 pancreas + pLN [pancreas damaged before Xenium] (slide 0041001)")
    runs = {r["run"]: r for r in _rows(folder / "runs.csv")}
    assert set(runs) == {"phenocycler:Synthetic_Pancreas:2026-01-02",
                         "phenocycler:Synthetic_pLN:2026-01-09", "xenium:R1_synthetic"}
    assert runs["phenocycler:Synthetic_Pancreas:2026-01-02"]["slides"] == "2"
    coverage = {r["donor_id"]: r for r in _rows(folder / "coverage.csv")}
    assert "donor notes" not in coverage["6374"]
    assert coverage["0112"][_at(STAGE.PC_SLIDE, STAGE.PHENOCYCLER)] == "2+ CD38; 2026-01-02 with 6374"
    # 0112's Xenium ran on its PhenoCycler slide (section_modes.csv); 6374's on a slide of its own
    assert coverage["0112"][_at(STAGE.PC_SLIDE, STAGE.XENIUM)] == "5k+100(6KTRAF); 2024-01-01 with 6374"
    assert not coverage["0112"][_at(STAGE.XENIUM_SLIDE, STAGE.XENIUM)]
    assert coverage["6374"][_at(STAGE.XENIUM_SLIDE, STAGE.XENIUM)] == (
        "5k+100(6KTRAF); 2024-01-01 with 0112")
    assert coverage["6374"][_at(STAGE.PLN_SLIDE, STAGE.PHENOCYCLER)] == "2 markers; 2026-01-09 alone"
    # each slide's number as its instrument recorded it: a qptiff's SlideID, a Xenium serial --
    # and 0112's PhenoCycler slide is also Xenium slide 0041001
    slide = STAGE.SLIDE_NUMBER
    assert [coverage["6374"][_at(group, slide)] for group in (STAGE.PC_SLIDE, STAGE.XENIUM_SLIDE,
                                                                 STAGE.PLN_SLIDE)] == [
        "6374", "0041000", "6374pLN"]
    assert coverage["0112"][_at(STAGE.PC_SLIDE, slide)] == "0041001"
    assert runs["phenocycler:Synthetic_Pancreas:2026-01-02"]["slides_run_together"] == (
        "6374 pancreas | 0112 pancreas (slide 0041001)")
    assert runs["xenium:R1_synthetic"]["slide_numbers"] == "0041000; 0041001"
    # a dataset is an acquisition, never a section: 0112's one Xenium run covers its pancreas and
    # its pLN, on the slide its PhenoCycler scan is of
    assert [(coverage[d][STAGE.DATASETS], coverage[d][STAGE.SLIDES]) for d in ("6374", "0112")] == [
        ("4", "3"), ("2", "1")]
    assert checks["datasets, by the slide and platform they came from"]["value"] == {
        _at(STAGE.PC_SLIDE, STAGE.PHENOCYCLER): 2, _at(STAGE.PC_SLIDE, STAGE.XENIUM): 1,
        _at(STAGE.XENIUM_SLIDE, STAGE.XENIUM): 1, _at(STAGE.XENIUM_SLIDE, STAGE.POST_XENIUM_IF): 1,
        _at(STAGE.PLN_SLIDE, STAGE.PHENOCYCLER): 1, "all": 6}
    assert not any(v in ("?", "*") or v.endswith((" ?", " *")) for r in coverage.values()
                   for v in r.values()), "nothing is left to a question mark or a footnote star"
    assert checks["PhenoCycler slides run alone"]["value"] == ["6374 pLN, 2026-01-09"]
    assert checks["panels, as recorded and named"]["value"] == {
        "phenocycler": {"2 markers": 1, "2+ CD38": 2}, "xenium": {"5k+100(6KTRAF)": 3},
        "if": {"INS, GCG, failed marker, AF": 1}}
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


def test_recorded_damage_is_marked_and_a_row_naming_no_acquisition_blocks(synthetic_project):
    _, folder, checks = _run(synthetic_project())
    acquisitions = {r["acquisition_id"]: r for r in _rows(folder / "acquisitions.csv")}
    assert acquisitions["xenium:0041001:Pancreas"]["known_condition"] == "pancreas damaged before Xenium"
    assert not acquisitions["xenium:0041001:pLN"].get("known_condition")
    coverage = {r["donor_id"]: r for r in _rows(folder / "coverage.csv")}
    pancreas, pln = (_at(STAGE.PC_SLIDE, c) for c in ("pancreas", STAGE.DIFFERENT_BLOCK))
    # the damaged Xenium ran on 0112's PhenoCycler slide: the mark is on that slide's pancreas
    assert coverage["0112"][pancreas] == "pancreas (damaged)"
    assert coverage["0112"][pln] == "1 pLN", "only the damaged tissue's section is marked"
    assert "damaged" not in coverage["0112"][_at(STAGE.PC_SLIDE, STAGE.XENIUM)]
    runs = {r["run"]: r for r in _rows(folder / "runs.csv")}
    assert "0112 pancreas + pLN [pancreas damaged before Xenium]" in runs["xenium:R1_synthetic"][
        "slides_run_together"]
    assert checks["known damaged tissue"]["value"] == [
        "xenium:0041001:Pancreas (0112): pancreas damaged before Xenium"]
    header = "acquisition_id,donor_id,tissue,condition,when,detail,recorded_by,recorded_on\n"
    _, folder, _ = _run(synthetic_project(release="r2", tissue_condition__csv=header + (
        "xenium:0041001:Pancreas,112,pancreas; pLN,damaged,before Xenium,x,t,2026-01-01\n")))
    coverage = {r["donor_id"]: r for r in _rows(folder / "coverage.csv")}
    assert (coverage["0112"][pancreas], coverage["0112"][pln]) == ("pancreas (damaged)",
                                                                   "1 pLN (damaged)")
    _, _, checks = _run(synthetic_project(release="r3", tissue_condition__csv=header + (
        "xenium:9999999:Panc,7000,pancreas,damaged,before Xenium,x,t,2026-01-01\n"
        "xenium:0041001:pLN,112,pancreas and pLN,damaged,before Xenium,x,t,2026-01-01\n")))
    check = checks["every recorded tissue condition names an acquisition and a known tissue"]
    assert not check["passed"] and check["value"][0] == "xenium:9999999:Panc (7000)"
    assert check["value"][1].startswith("xenium:0041001:pLN (112): unknown tissue 'pancreas and pLN'")


def test_a_xenium_slides_pln_lymphoid_is_read_from_its_inspection_only(synthetic_project,
                                                                       tmp_path):
    from conftest import INSPECTION_HEADER, SOURCES, write_fovs, write_xenium_records

    output = "output-XETG00298__0041000__pLN__20240101__000000"
    # the notes give 6374's sequential PhenoCycler slide a same-block pLN; that does not carry
    # over to its Xenium slide (maintainer, 2026-10-23)
    notes = ("IO60panc2nd Notes\n6374: pancTail and 1 pLN same-block, sequential Xenium\n"
             "HDL112: pancTail and 1 pLN different-block, same-slide post-phenocycler Xenium\n")

    def run(inspection, release):
        path = synthetic_project(release=release, notes__txt=notes, notes_copy__md=notes,
                                 pln_blocks__csv=BLOCKS, xenium_lymphoid_inspection__csv=inspection,
                                 xenium_paths__csv=(SOURCES["xenium_paths.csv"]
                                                    + f"6374,pln_1,B1,bundles/{output}\n"))
        write_xenium_records(tmp_path / "xenium" / "bundles" / output, run="R1_synthetic",
                             start="2024-01-01T10:00:00Z", slide="0041000", region="pLN")
        write_fovs(tmp_path / "xenium" / "bundles" / output / "aux_outputs"
                   / "overview_scan_fov_locations.json", (1, 8, 4, 3))
        _, folder, checks = _run(path)
        return ({r["acquisition_id"]: r for r in _rows(folder / "acquisitions.csv")},
                {r["donor_id"]: r for r in _rows(folder / "coverage.csv")}, checks)

    xenium = [_at(STAGE.XENIUM_SLIDE, c) for c in (STAGE.SAME_BLOCK, STAGE.DIFFERENT_BLOCK)]
    # inspected: nothing lymphoid in the pancreas area, so the pLN region is a section of its own
    acquisitions, coverage, checks = run(SOURCES["xenium_lymphoid_inspection.csv"], "r1")
    assert acquisitions["xenium:0041000:pLN"]["block"] == STAGE.OTHER_BLOCK
    assert acquisitions["xenium:0041000:Pan"]["lymphoid_in_pancreas"] == "0"
    assert [coverage["6374"][c] for c in xenium] == ["", "1 pLN"]
    # 0112's Xenium ran on its PhenoCycler slide, whose sections the notes give; a dedicated
    # PhenoCycler pLN slide is from a pLN block
    assert acquisitions["xenium:0041001:pLN"]["block"] == STAGE.OTHER_BLOCK
    assert acquisitions["phenocycler:6374pLN_Scan1.qptiff"]["block"] == STAGE.OTHER_BLOCK
    assert coverage["6374"][_at(STAGE.PLN_SLIDE, STAGE.PLN_SECTION)] == "1 pLN", (
        "pln_node_annotations.csv records 1 node on 6374's pLN slide")
    # the Xenium slide's two regions are one run of one slide: one dataset
    assert (coverage["6374"][STAGE.DATASETS], coverage["6374"][STAGE.SLIDES]) == ("4", "3")
    assert all(c["passed"] for c in checks.values() if c["blocking"])
    assert checks["pLN/lymphoid seen in each Xenium slide's pancreas area"]["value"] == [
        "0041000 (6374): 0 in the pancreas area (t 2026-01-01)"]

    # inspected: two lymphoid areas in the pancreas area, the pLN region one of them
    acquisitions, coverage, checks = run(
        INSPECTION_HEADER + "0041000,6374,2,pLN 1,t,2026-01-02,embedded node\n", "r2")
    assert acquisitions["xenium:0041000:pLN"]["block"] == STAGE.PANCREAS_BLOCK
    assert [coverage["6374"][c] for c in xenium] == ["2 pLN", ""]
    assert checks["pLN/lymphoid seen in each Xenium slide's pancreas area"]["value"] == [
        "0041000 (6374): 2 in the pancreas area, pLN 1 part of that section (t 2026-01-02)"]
    assert all(c["passed"] for c in checks.values() if c["blocking"])


def test_a_xenium_slide_not_inspected_blocks_and_says_so(synthetic_project):
    from conftest import INSPECTION_HEADER

    _, folder, checks = _run(synthetic_project(xenium_lymphoid_inspection__csv=INSPECTION_HEADER))
    check = checks["every Xenium slide's pancreas area has been inspected for pLN/lymphoid"]
    # 0112's slide 0041001 is its PhenoCycler slide (section_modes.csv), so only 6374's is asked
    assert not check["passed"] and check["blocking"] and check["value"] == ["0041000 (6374)"]
    coverage = {r["donor_id"]: r for r in _rows(folder / "coverage.csv")}
    assert coverage["6374"][_at(STAGE.XENIUM_SLIDE, STAGE.SAME_BLOCK)] == STAGE.NOT_INSPECTED
    assert (folder / "gate" / "xenium_slides.png").stat().st_size > 1000


@pytest.mark.parametrize("row, problem", [
    ("0049999,7000,0,,t,2026-01-01,", "slide 0049999: no Xenium slide of its own has that serial"),
    ("0041001,112,0,,t,2026-01-01,", "slide 0041001: no Xenium slide of its own has that serial"),
    ("0041000,6374,some,,t,2026-01-01,", "slide 0041000: 'some' is not a count"),
    ("0041000,6374,1,pLN9,t,2026-01-01,", "slide 0041000: pLN9 is not one of its pLN regions"),
    ("0041000,6374,1,,,,", "slide 0041000: names no inspector and date"),
])
def test_an_inspection_row_that_cannot_be_applied_blocks(synthetic_project, row, problem):
    from conftest import INSPECTION_HEADER

    _, _, checks = _run(synthetic_project(
        xenium_lymphoid_inspection__csv=INSPECTION_HEADER + row + "\n"))
    check = checks["every inspection row names a Xenium slide, its own pLN regions, a count and "
                   "who inspected it"]
    assert not check["passed"] and check["blocking"] and check["value"] == [problem]


def test_a_xenium_run_on_a_phenocycler_slide_with_no_output_listed_is_shown_and_reported(
        synthetic_project):
    _, folder, checks = _run(synthetic_project(section_modes__csv=(
        "donor_id,roi,xenium_serial,mode,evidence,note\n"
        "112,panc,0041999,same_slide,the maintainer,\n")))
    coverage = {r["donor_id"]: r for r in _rows(folder / "coverage.csv")}
    assert coverage["0112"][_at(STAGE.PC_SLIDE, STAGE.XENIUM)] == STAGE.NOT_LISTED
    # 0041001 is no longer 0112's PhenoCycler slide, so its Xenium is a slide of its own
    assert coverage["0112"][_at(STAGE.XENIUM_SLIDE, STAGE.XENIUM)] == (
        "5k+100(6KTRAF); 2024-01-01 with 6374")
    assert (coverage["0112"][STAGE.DATASETS], coverage["0112"][STAGE.SLIDES]) == ("2", "2")
    check = checks["Xenium run on a PhenoCycler slide with no output listed"]
    assert check["value"] == ["0112: slide 0041999 (section_modes.csv: the maintainer)"]
    assert not check["blocking"]


ROLES = ("image_file,donor_id,roi,xenium_serial,channel_dapi,channel_opal520,channel_opal570,"
         "channel_opal690,channel_sampleaf\n")


@pytest.mark.parametrize("roles, passed", [
    ("NUCLEUS,GCG,FAILED,INS,AF", True),          # the declared panel, as restain_sections.csv records it
    ("NUCLEUS,GCG,SST,INS,AF", False),            # another panel: refused, never renamed
])
def test_every_post_xenium_if_image_carries_the_one_declared_panel(synthetic_project, roles,
                                                                    passed):
    _, folder, checks = _run(synthetic_project(
        restain_sections__csv=ROLES + f"if_6374.qptiff,6374,panc,0041000,{roles}\n"))
    acquisition = {r["acquisition_id"]: r for r in _rows(folder / "acquisitions.csv")}[
        "if:if_6374.qptiff"]
    check = checks["every panel is named by its platform's naming rule"]
    assert check["passed"] is passed
    if passed:
        assert acquisition["panel"] == "INS, GCG, failed marker, AF"
    else:
        assert "SST" in acquisition["panel_error"] and check["value"][0].startswith("if:if_6374")


def test_a_cassette_name_that_names_another_donor_blocks_until_its_error_is_known(
        synthetic_project, tmp_path, monkeypatch):
    from conftest import write_xenium_records

    path = synthetic_project()
    bundles = tmp_path / "xenium" / "bundles"
    for output, slide, region, cassette in (
            ("output-XETG00298__0041000__Pan__20240101__000000", "0041000", "Pan", "nPOD6374"),
            ("output-XETG00298__0041001__Pancreas__20240101__000000", "0041001", "Pancreas",
             "nPOD6553")):
        write_xenium_records(bundles / output, run="R1_synthetic", start="2024-01-01T10:00:00Z",
                             slide=slide, region=region, cassette=cassette)
    _, folder, checks = _run(path)
    acquisitions = {r["acquisition_id"]: r for r in _rows(folder / "acquisitions.csv")}
    assert acquisitions["xenium:0041000:Pan"]["cassette_name"] == "nPOD6374"
    name = "every Xenium cassette name carries the listed donor, or its error is known"
    # nPOD6374 carries 6374; nPOD6553 does not carry 112, the donor the list gives that slide
    assert not checks[name]["passed"] and checks[name]["blocking"]
    assert checks[name]["value"] == [
        "output-XETG00298__0041001__Pancreas__20240101__000000: cassette nPOD6553, listed as "
        "donor 112"]
    monkeypatch.setitem(STAGE.KNOWN_LABEL_ERRORS, "0041001", "typed in error")
    tasks.run_stage(STAGE, P.load(path), force=True)
    checks = {c["name"]: c for c in gates.load(folder)["checks"]}
    assert checks[name]["passed"]
    assert checks["labels that name another donor, with their known reason"]["value"] == [
        "slide 0041001: cassette nPOD6553, listed as donor 112 (typed in error)"]


@pytest.mark.parametrize("run, batch", [
    ("MM035_Xenium-npodBatch1", "1"), ("MM036_XeniumBatch2_ReRun", "2"), ("MM44_batch10", "10"),
    ("MM48_batch14", "14"), ("6533_HDL112_pLNpanc_PhenoXeno", ""),
    ("PhenoXeno_pLNpanc_6580_6347", ""),
])
def test_a_xenium_run_names_its_batch_number_or_none(run, batch):
    assert STAGE._batch_number(run) == batch


def test_two_sources_that_disagree_on_a_disease_status_block(synthetic_project):
    _, _, checks = _run(synthetic_project(npod=[{"case_id": "6374", "donor_type": "No Diabetes"}]))
    check = checks["the donor metadata and nPOD's case export agree on each donor's disease status"]
    assert not check["passed"] and check["value"] == ["6374: donor_metadata T1D, npod_cases ND"]


def test_a_donor_with_no_recorded_disease_status_blocks(synthetic_project):
    _, _, checks = _run(synthetic_project(metadata=[{"donor_id": "6374",
                                                     "disease.status": "T1D"}]))
    check = checks["every donor has a recorded disease status"]
    assert not check["passed"] and check["value"] == ["0112"]


def test_a_xenium_only_donor_has_its_own_rows_status_and_batch(synthetic_project, tmp_path):
    from conftest import NPOD_ROWS, SOURCES, write_xenium_records

    output = "output-XETG00298__0049999__Region_1__20260416__180001"
    path = synthetic_project(
        npod=[*NPOD_ROWS, {"case_id": "7000", "donor_type": "Autoantibody Positive"}],
        xenium_paths__csv=SOURCES["xenium_paths.csv"] + f"7000,panc,XeniumOnly,bundles/{output}\n",
        xenium_lymphoid_inspection__csv=(SOURCES["xenium_lymphoid_inspection.csv"]
                                         + "0049999,7000,0,,t,2026-01-01,\n"))
    write_xenium_records(tmp_path / "xenium" / "bundles" / output, run="MM48_batch14",
                         start="2026-04-16T17:59:37Z", slide="0049999", region="Region_1",
                         design="Y9QRJ6", cassette="7000")
    _, folder, checks = _run(path)
    coverage = {r["donor_id"]: r for r in _rows(folder / "coverage.csv")}
    row = coverage["7000"]
    assert (row["design"], row["disease_status"], row["disease_status_source"]) == (
        "Xenium only", "Aab+", "npod_cases")
    assert row[_at(STAGE.XENIUM_SLIDE, STAGE.XENIUM)] == "5k+100(Y9QRJ6); batch 14; 2026-04-16 alone"
    assert (row[STAGE.DATASETS], row[STAGE.SLIDES]) == ("1", "1")
    assert checks["Xenium-only donors"]["value"] == [
        "7000 (5k+100(Y9QRJ6); batch 14; 2026-04-16 alone)"]
    assert checks["Xenium batches"]["value"] == ["batch 14: MM48_batch14, 7000"]
    assert checks["disease status, by donor"]["value"] == {"ND": ["0112"], "Aab+": ["7000"],
                                                           "T1D": ["6374"]}
    assert all(c["passed"] for c in checks.values() if c["blocking"])


def test_a_phenocycler_slide_label_that_names_another_donor_blocks_until_its_error_is_known(
        synthetic_project, tmp_path, monkeypatch):
    import datetime

    from conftest import PANCREAS_MARKERS, write_qptiff

    path = synthetic_project()
    write_qptiff(tmp_path / "raw" / "6374.qptiff", slide="6623panc",
                 experiment="Synthetic_Pancreas", position="Slide1",
                 experiment_date="2025-12-01T09:00:00Z",
                 scanned=datetime.datetime(2026, 1, 2, 14, 0), markers=PANCREAS_MARKERS)
    _, folder, checks = _run(path)
    name = "every PhenoCycler slide label carries the listed donor, or its error is known"
    assert not checks[name]["passed"] and checks[name]["blocking"]
    assert checks[name]["value"] == ["6374.qptiff: SlideID 6623panc, listed as donor 6374"]
    monkeypatch.setitem(STAGE.KNOWN_LABEL_ERRORS, "6374.qptiff", "entered in error")
    tasks.run_stage(STAGE, P.load(path), force=True)
    checks = {c["name"]: c for c in gates.load(folder)["checks"]}
    assert checks[name]["passed"]
    assert checks["labels that name another donor, with their known reason"]["value"] == [
        "6374.qptiff: SlideID 6623panc, listed as donor 6374 (entered in error)"]


def test_each_slide_shows_a_thumbnail_from_its_own_image(synthetic_project):
    _, folder, checks = _run(synthetic_project())
    # the qptiffs carry a Thumbnail page and both Xenium slides an overview scan; the pLN
    # slide's qptiff has none, and is named rather than drawn blank
    assert checks["slides with no thumbnail"]["value"] == [
        "6374 PhenoCycler pLN slide (6374pLN_Scan1.qptiff)"]
    assert not checks["slides with no thumbnail"]["blocking"]
    assert (folder / "gate" / "coverage.png").stat().st_size > 1000


def test_a_phenocycler_slide_that_xenium_also_ran_on_shows_both_images(synthetic_project,
                                                                       monkeypatch):
    drawn = {}
    figure = STAGE._figure

    def keep(task, layout):
        drawn.update(layout["thumbs"])
        figure(task, layout)

    monkeypatch.setattr(STAGE, "_figure", keep)
    _run(synthetic_project())
    # 0112's PhenoCycler slide is Xenium slide 0041001 (section_modes.csv): its qptiff's
    # thumbnail, then the Xenium overview scan; 6374's PhenoCycler slide has only its qptiff's
    assert len(drawn[("0112", STAGE.PC_SLIDE)]) == 2
    assert drawn[("0112", STAGE.PC_SLIDE)][0].shape[:2] == (8, 6), "the qptiff's comes first"
    assert len(drawn[("6374", STAGE.PC_SLIDE)]) == 1
    assert len(drawn[("6374", STAGE.XENIUM_SLIDE)]) == 1


def test_an_exclusion_takes_out_only_the_platform_it_covers(synthetic_project, monkeypatch):
    """As 6457's does (the maintainer, 2026-10-22: "the xenium for 6457 is not excluded")."""
    from multiplex5000 import cohort

    monkeypatch.setitem(cohort.DONOR_EXCLUSIONS, "6374", "a PhenoCycler reason")
    monkeypatch.setitem(cohort.EXCLUDED_PLATFORMS, "6374", ("phenocycler",))
    _, folder, checks = _run(synthetic_project())
    row = {r["donor_id"]: r for r in _rows(folder / "coverage.csv")}["6374"]
    assert row["excluded_platforms"] == "phenocycler"
    datasets = checks["datasets, by the slide and platform they came from"]
    # 6374's PhenoCycler pancreas and pLN scans are not counted; its Xenium slide and IF are
    assert datasets["value"] == {
        _at(STAGE.PC_SLIDE, STAGE.PHENOCYCLER): 1, _at(STAGE.PC_SLIDE, STAGE.XENIUM): 1,
        _at(STAGE.XENIUM_SLIDE, STAGE.XENIUM): 1, _at(STAGE.XENIUM_SLIDE, STAGE.POST_XENIUM_IF): 1,
        _at(STAGE.PLN_SLIDE, STAGE.PHENOCYCLER): 0, "all": 4}
    assert "6374 (2)" in datasets["detail"]


def test_a_same_slide_xenium_overview_is_turned_the_way_the_phenocycler_image_is():
    import numpy as np

    reference = np.zeros((40, 20, 3), np.uint8)
    reference[2:12, 2:8] = 200            # a small piece at the top left
    reference[16:38, 4:18] = 200          # the large section below it
    turned = reference[::-1, ::-1].copy()
    assert np.array_equal(STAGE._oriented_like(reference, turned), reference)
    assert np.array_equal(STAGE._oriented_like(reference, reference), reference)


def test_the_phenocycler_scan_is_read_at_its_smallest_level_tall_enough_to_inspect(tmp_path):
    import numpy as np
    import tifffile

    path = tmp_path / "scan.qptiff"
    reduced = "<PerkinElmer-QPI-ImageDescription><ImageType>ReducedResolution</ImageType>"
    tifffile.imwrite(path, np.zeros((2000, 10), np.uint8))
    for height, value in ((1000, 7), (500, 9)):        # each level's first page, then a second
        for channel in range(2):
            tifffile.imwrite(path, np.full((height, 5), value + channel, np.uint8), append=True,
                             description=reduced + "</PerkinElmer-QPI-ImageDescription>")
    image, full = STAGE._phenocycler_image(path)
    # 1000 rows is the smallest level at least INSPECTION_PIXELS (900) tall; its first page
    assert image.shape == (1000, 5) and full == (2000, 10)


def test_annotation_outlines_are_read_by_tissue(tmp_path):
    from conftest import write_annotations

    square = [[[0, 0], [4, 0], [4, 4], [0, 4], [0, 0]]]
    path = write_annotations(tmp_path / "a.geojson", ("Tissue", square), ("pLN", square),
                             ("Other", square))
    rings = STAGE._annotation_rings(path, "pancreas")
    # Tissue is the image's own tissue, pLN the tissue profiles' lymph node, Other not drawn
    assert [t for t, _ in rings] == ["pancreas", "pancreatic_lymph_node"]
    assert rings[0][1].shape == (5, 2)
    assert STAGE._export_tag("6374_Scan1.er.qptiff - resolution #1") == (
        "6374_Scan1.er.qptiff_-_resolution_1")


@pytest.mark.parametrize("donor, scans, sequential, said", [
    ("0112", [{"qptiff": "0112.qptiff"}], False,
     "no PhenoCycler section\nrecorded as sequential\nto this slide"),
    ("7000", [], True, "no PhenoCycler slide"),
    ("6457", [], True, "PhenoCycler data excluded"),
])
def test_a_phenocycler_slide_is_paired_only_as_a_recorded_sequential_section(donor, scans,
                                                                           sequential, said):
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots()
    # a same-slide donor's PhenoCycler slide is the one Xenium ran on after it, not a section
    # next to its Xenium slide of its own (maintainer, 2026-10-23: "some of them look very
    # different e.g. 6347")
    STAGE._draw_phenocycler(None, ax, donor, scans, sequential)
    assert [t.get_text() for t in ax.texts] == [said]
    plt.close(fig)


@pytest.mark.parametrize("acquisition, label", [
    ({"tissue": "pancreas", "roi": "panc"}, "pancreas"),
    ({"tissue": "pancreatic_lymph_node", "roi": "pln_1"}, "pLN 1"),
    ({"tissue": "pancreatic_lymph_node", "roi": "pln_3"}, "pLN 3"),
])
def test_a_xenium_region_is_labelled_by_its_tissue_not_its_typed_name(acquisition, label):
    # run MM044 typed PanT for its pLN regions; the label says what xenium_paths.csv says
    assert STAGE._region_label({**acquisition, "acquisition_id": "xenium:0059907:PanT"}) == label
