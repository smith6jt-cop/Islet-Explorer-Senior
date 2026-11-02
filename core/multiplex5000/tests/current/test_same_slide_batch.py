"""Exact pair projects progress independently without invented scientific approval."""

import csv
import json
from dataclasses import asdict
from pathlib import Path

import numpy as np
import pandas as pd
import pyarrow as pa
import pyarrow.parquet as pq
import pytest

from multiplex5000 import gates, project, tasks
from multiplex5000 import same_slide_batch as batch
from multiplex5000.artifacts import sha256_file
from multiplex5000.integration.transform import Transform
from multiplex5000.same_slide_sources import load_pairs
from multiplex5000.stages import same_slide_registration as registration


def _write(path, value):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value))


def _table(path, values):
    path.parent.mkdir(parents=True, exist_ok=True)
    pq.write_table(pa.Table.from_pandas(pd.DataFrame(values), preserve_index=False), path)


def _fixture(tmp_path, *, count=3, accepted=True, last_transform=False, gene_only=False):
    workspace = tmp_path / "pair_releases"
    configs, rows, projects = [], [], []
    for i in range(count):
        donor = f"HDL{101 + i}"
        region = str(101 + i)
        entry = f"{donor}_PC67_0110618-S{i + 1}"
        acquisition = f"xenium:0110618:{region}"
        output = f"/synthetic/output-XETG00298__0110618__{region}__20260826__224731"
        source = tmp_path / f"source_{i}.yaml"
        _write(source, {"schema": "multiplex5000-same-slide-sources/1", "tissue": "spleen",
                       "relationship": "same_section_after_phenocycler", "protein_frame": "crop_local_um",
                       "xenium_frame": "region_local_um", "expected_pairs": 1,
                       "section_table": {"path": "/synthetic/sections.csv", "sha256": "0" * 64},
                       "pairs": [{"donor_id": donor, "protein_entry": entry,
                                  "protein_export": f"/synthetic/{entry}.csv.gz",
                                  "protein_image": "/synthetic/0110618.qptiff",
                                  "xenium_acquisition_id": acquisition, "xenium_output": output}]})
        pair = load_pairs(source)[0]
        row = json.loads(json.dumps(asdict(pair), default=str))
        row["pair_id"] = f"synthetic-pair-{i}"
        rows.append(row)
        release = f"candidate-{i}"
        folder = workspace / release
        sample = f"sample-{i}"  # Explicit region table key; never parsed from donor.
        cells_p = folder / f"03_protein_measurement/cells/{entry}.parquet"
        cells_r = folder / f"04_rna_ingest/cells/{sample}.parquet"
        p_uid = [f"phenocycler:0110618:{entry}:seg:{j}" for j in [1, 2]]
        x_uid = [f"{acquisition}:cells:{j}" for j in [1, 2]]
        _table(cells_p, {"cell_uid": p_uid, "object_id": ["1", "2"], "donor_id": [donor] * 2,
                         "tissue": ["spleen"] * 2, "acquisition_id": ["phenocycler:0110618"] * 2,
                         "x_um": [100., 110.], "y_um": [0., 0.], "protein_eligible": [True, True]})
        _table(cells_r, {"cell_uid": x_uid, "donor_id": [donor] * 2, "tissue": ["spleen"] * 2,
                         "acquisition_id": [acquisition] * 2, "x_um": [.5, 10.5], "y_um": [0., 0.],
                         "typing_eligible": [True, True]})
        with (folder / "04_rna_ingest/regions.csv").open("w", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=["acquisition_id", "sample_key", "donor_id", "tissue", "source_path"])
            writer.writeheader()
            writer.writerow({"acquisition_id": acquisition, "sample_key": sample,
                             "donor_id": donor, "tissue": "spleen", "source_path": output})
        calls_p = folder / f"06_typing/protein/{entry}.parquet"
        calls_r = folder / f"06_typing/rna/{sample}.parquet"
        _table(calls_p, {"cell_uid": p_uid, "object_id": ["1", "2"],
                         "candidate_phenotype": ["T cell", None],
                         "candidate_state": ["unique_candidate", "unavailable"]})
        rna = {"cell_uid": x_uid, "typing_eligible": [True, True]}
        rna.update({"gene:CD3E": [1, -1]} if gene_only else {"cell_class": ["T cell", None]})
        _table(calls_r, rna)
        _write(folder / "06_typing/receipt.json", {
            "stage": "06_typing", "status": "complete", "release": release,
            "inputs": [{"name": name, "sha256": sha256_file(path)} for name, path in [
                ("same_slide_source_manifest", source), (f"protein_cells:{entry}", cells_p),
                (f"rna_cells:{sample}", cells_r)]],
            "outputs": [{"path": str(path.relative_to(folder / "06_typing")), "sha256": sha256_file(path)}
                        for path in [calls_p, calls_r]]})
        _write(folder / "06_typing/gate/gate.json", {"stage": "06_typing", "release": release,
                                                    "checks": [], "passed": True})
        gates.accept(folder / "06_typing", by="synthetic-test")
        transform_root = tmp_path / "transforms"
        if i != count - 1 or last_transform:
            Transform(affine=np.array([[1., 0., 100.], [0., 1., 0.]]), meta={
                "protein_entry": entry, "xenium_acquisition_id": acquisition,
                "moving_frame": "region_local_um", "fixed_frame": "crop_local_um",
                "reviewed": True}).save(transform_root / entry)
        project_path = tmp_path / f"pair_{i}.yaml"
        _write(project_path, {"schema": "multiplex5000-project/1", "id": f"pair_{i}",
                              "release": release, "cohort": "spleen_same_slide_pilot",
                              "stages": ["06_typing", "07_registration", "08_integration"],
                              "automatic_gates": [], "releases_root": str(workspace),
                              "sources": {"same_slide_sources": str(source), "same_slide_registrations": str(transform_root)},
                              "settings": {"rna_acquisitions": [acquisition], "same_slide_hierarchy_map": {
                                  "tissue": "spleen", "level": "lineage", "protein_column": "candidate_phenotype",
                                  "rna_column": "cell_class", "protein": {"T cell": "T"}, "rna": {"T cell": "T"}}}})
        selected = project.load(project_path)
        projects.append(selected)
        if accepted and (i != count - 1 or last_transform):
            run = tasks.run_stage(registration, selected)
            gates.accept(run.path, by="synthetic-test")
        configs.append({"pair_id": row["pair_id"], "project": str(project_path), "required_evidence": "independent_calls"})
    roster = tmp_path / "roster.json"
    _write(roster, {"explicit_pairs": rows, "unresolved": []})
    config = tmp_path / "batch.json"
    _write(config, {"schema": batch.SCHEMA, "roster": {"path": str(roster), "sha256": sha256_file(roster)},
                    "workspace_root": str(workspace), "pairs": configs})
    return config, projects


def test_two_ready_pairs_complete_while_missing_geometry_remains_accounted(tmp_path):
    config, projects = _fixture(tmp_path)
    result = batch.coordinate(config, tmp_path / "ledger", action="integration")
    assert result["counts"] == {"completed": 2, "blocked": 1}
    assert result["declared_pairs"] == result["accounted_pairs"] == 3
    assert result["all_pairs_accounted"] and not result["full_cohort_candidate_completion"]
    assert result["pairs"][0]["summary"]["candidate_links"] == 2
    assert result["pairs"][0]["identity_comparison"] == "declared"
    assert any(x["code"] == "registration_source_not_ready" for x in result["pairs"][2]["requirements"])
    assert all(not (p.release_dir / "08_integration/gate/accepted.json").exists() for p in projects)
    assert (projects[0].release_dir / "08_integration/receipt.json").is_file()
    assert (projects[1].release_dir / "08_integration/receipt.json").is_file()
    assert not (projects[2].release_dir / "08_integration").exists()


def test_preflight_does_not_write_project_files_and_reports_ready_pairs(tmp_path):
    config, projects = _fixture(tmp_path)
    before = {str(p): sha256_file(p) for project in projects for p in project.release_dir.rglob("*") if p.is_file()}
    result = batch.coordinate(config, tmp_path / "ledger")
    after = {str(p): sha256_file(p) for project in projects for p in project.release_dir.rglob("*") if p.is_file()}
    assert before == after
    assert result["pairs"][0]["integration"] == "ready_to_run"
    assert result["counts"] == {"blocked": 3}


def test_registration_runs_independently_and_stops_at_its_own_review_gate(tmp_path):
    config, projects = _fixture(tmp_path, accepted=False)
    result = batch.coordinate(config, tmp_path / "ledger", action="registration")
    assert result["pairs"][0]["registration"] == "published_awaiting_review"
    assert result["pairs"][1]["registration"] == "published_awaiting_review"
    assert all(not (p.release_dir / "07_registration/gate/accepted.json").exists() for p in projects)
    assert all(not (p.release_dir / "08_integration").exists() for p in projects)


def test_stale_pair_inputs_do_not_invalidate_other_pairs(tmp_path):
    config, projects = _fixture(tmp_path)
    calls = next((projects[0].release_dir / "06_typing/rna").glob("*.parquet"))
    table = pq.read_table(calls).to_pandas()
    table["cell_class"] = "changed"
    pq.write_table(pa.Table.from_pandas(table), calls)
    result = batch.coordinate(config, tmp_path / "ledger", action="integration")
    assert [row["status"] for row in result["pairs"]] == ["blocked", "completed", "blocked"]
    assert result["pairs"][0]["requirements"][0]["code"] == "independent_typing_not_ready"


def test_reuse_rechecks_settings_and_bound_outputs(tmp_path):
    config, projects = _fixture(tmp_path, count=2, last_transform=True)
    first = batch.coordinate(config, tmp_path / "first", action="integration")
    assert first["full_cohort_candidate_completion"]
    first_receipt = sha256_file(projects[0].release_dir / "08_integration/receipt.json")
    second = batch.coordinate(config, tmp_path / "second", action="integration")
    assert second["counts"] == {"completed": 2}
    assert sha256_file(projects[0].release_dir / "08_integration/receipt.json") == first_receipt
    spec = json.loads(projects[0].path.read_text())
    spec["settings"]["same_slide_hierarchy_map"]["level"] = "another-level"
    _write(projects[0].path, spec)
    path = projects[1].release_dir / "08_integration/pairing_counts.json"
    path.write_text(path.read_text() + "\n")
    third = batch.coordinate(config, tmp_path / "third", action="integration")
    assert third["counts"] == {"blocked": 2}
    assert "settings differ" in third["pairs"][0]["requirements"][0]["detail"]
    assert "output pairing_counts.json changed" in third["pairs"][1]["requirements"][0]["detail"]


def test_gene_evidence_cannot_satisfy_independent_call_mode(tmp_path):
    config, projects = _fixture(tmp_path, count=1, last_transform=True, gene_only=True)
    result = batch.coordinate(config, tmp_path / "ledger", action="integration")
    assert result["pairs"][0]["status"] == "unavailable"
    assert result["pairs"][0]["registration"] == "accepted"
    assert result["pairs"][0]["identity_comparison"] == "unavailable"
    assert not (projects[0].release_dir / "08_integration").exists()


def test_explicit_spatial_gene_mode_does_not_claim_an_identity_comparison(tmp_path):
    config, projects = _fixture(tmp_path, count=1, last_transform=True, gene_only=True)
    spec = json.loads(config.read_text())
    spec["pairs"][0]["required_evidence"] = "spatial_gene_evidence"
    _write(config, spec)
    project_spec = json.loads(projects[0].path.read_text())
    project_spec["settings"].pop("same_slide_hierarchy_map")
    _write(projects[0].path, project_spec)
    result = batch.coordinate(config, tmp_path / "ledger", action="integration")
    assert result["counts"] == {"completed": 1}
    assert result["pairs"][0]["identity_comparison"] == "unavailable"
    assert not result["scientific_approval"]


@pytest.mark.parametrize("change", ["duplicate", "omit", "roster_hash", "shared_root"])
def test_invalid_batch_never_silently_reduces_or_overlaps_the_cohort(tmp_path, change):
    config, projects = _fixture(tmp_path, accepted=False)
    spec = json.loads(config.read_text())
    if change == "duplicate":
        spec["pairs"].append(spec["pairs"][0])
    elif change == "omit":
        spec["pairs"].pop()
    elif change == "roster_hash":
        spec["roster"]["sha256"] = "0" * 64
    else:
        spec["pairs"][1]["project"] = spec["pairs"][0]["project"]
    _write(config, spec)
    with pytest.raises(batch.BatchError):
        batch.coordinate(config, tmp_path / "ledger")
    assert not (tmp_path / "ledger").exists()


@pytest.mark.parametrize("change", ["frame", "automatic_gate", "outside_workspace", "wrong_pair"])
def test_unsafe_or_mismatched_pair_configuration_fails_closed(tmp_path, change):
    config, projects = _fixture(tmp_path, accepted=False)
    spec = json.loads(projects[0].path.read_text())
    if change == "frame":
        folder = projects[0].source("same_slide_registrations") / "HDL101_PC67_0110618-S1"
        transform = Transform.load(folder)
        transform.meta["fixed_frame"] = "full_slide_um"
        transform.save(folder)
    elif change == "automatic_gate":
        spec["automatic_gates"] = ["07_registration"]
    elif change == "outside_workspace":
        spec["releases_root"] = str(tmp_path / "production")
    else:
        spec["settings"]["rna_acquisitions"] = ["xenium:0110618:102"]
    _write(projects[0].path, spec)
    result = batch.coordinate(config, tmp_path / "ledger", action="registration")
    assert result["pairs"][0]["status"] == "blocked"
    assert not (projects[0].release_dir / "07_registration").exists()
    assert result["pairs"][1]["registration"] == "published_awaiting_review"


def test_missing_adapter_and_unresolved_relationships_are_preserved(tmp_path):
    config, _ = _fixture(tmp_path, count=1, accepted=False)
    spec = json.loads(config.read_text())
    roster_path = Path(spec["roster"]["path"])
    roster = json.loads(roster_path.read_text())
    roster["explicit_pairs"].append({"pair_id": "io60-recorded", "tissue": "pancreas",
                                     "relationship": "same_section", "rna_acquisition_ids": ["xenium:exact:record"]})
    roster["unresolved"] = [{"tissue": "pancreatic_lymph_node", "reason": "physical relationship unavailable"}]
    _write(roster_path, roster)
    spec["roster"]["sha256"] = sha256_file(roster_path)
    spec["pairs"].append({"pair_id": "io60-recorded", "project": None, "required_evidence": "independent_calls"})
    _write(config, spec)
    result = batch.coordinate(config, tmp_path / "ledger")
    assert result["counts"] == {"blocked": 1, "unavailable": 1}
    assert result["pairs"][1]["requirements"][0]["code"] == "exact_section_adapter_unavailable"
    assert result["unresolved_relationships"] == roster["unresolved"]
    assert result["accounted_pairs"] == 2


def test_cli_uses_nonzero_status_for_incomplete_cohort_and_fresh_ledgers(tmp_path, capsys):
    config, _ = _fixture(tmp_path, count=1, accepted=False)
    out = tmp_path / "ledger"
    assert batch.main(["preflight", "--config", str(config), "--out", str(out)]) == 1
    assert json.loads(capsys.readouterr().out)["accounted_pairs"] == 1
    with pytest.raises(FileExistsError):
        batch.coordinate(config, out)


def test_ledger_cannot_occupy_a_pair_stage_output(tmp_path):
    config, projects = _fixture(tmp_path, count=1, accepted=False)
    with pytest.raises(batch.BatchError, match="outside per-pair"):
        batch.coordinate(config, projects[0].release_dir / "07_registration")
    assert not (projects[0].release_dir / "07_registration").exists()


def test_mixed_cell_frames_are_rejected_during_preflight(tmp_path):
    config, projects = _fixture(tmp_path, count=1, accepted=False, last_transform=True)
    project = projects[0]
    path = next((project.release_dir / "03_protein_measurement/cells").glob("*.parquet"))
    cells = pq.read_table(path).to_pandas()
    cells["coordinate_frame"] = ["crop_local_um", "full_slide_um"]
    pq.write_table(pa.Table.from_pandas(cells), path)
    receipt_path = project.release_dir / "06_typing/receipt.json"
    receipt = json.loads(receipt_path.read_text())
    next(row for row in receipt["inputs"] if row["name"].startswith("protein_cells:"))["sha256"] = sha256_file(path)
    _write(receipt_path, receipt)
    gates.accept(project.release_dir / "06_typing", by="synthetic-test")
    result = batch.coordinate(config, tmp_path / "ledger")
    assert result["pairs"][0]["status"] == "blocked"
    assert "coordinate frame" in result["pairs"][0]["requirements"][0]["detail"]


def test_source_roster_cannot_duplicate_one_physical_pair_under_new_ids(tmp_path):
    config, _ = _fixture(tmp_path, count=1, accepted=False)
    spec = json.loads(config.read_text())
    roster_path = Path(spec["roster"]["path"])
    roster = json.loads(roster_path.read_text())
    duplicate = dict(roster["explicit_pairs"][0], pair_id="duplicate-physical-pair")
    roster["explicit_pairs"].append(duplicate)
    _write(roster_path, roster)
    spec["roster"]["sha256"] = sha256_file(roster_path)
    spec["pairs"].append({"pair_id": duplicate["pair_id"], "project": None, "required_evidence": "independent_calls"})
    _write(config, spec)
    with pytest.raises(batch.BatchError, match="physical pair"):
        batch.coordinate(config, tmp_path / "ledger")
