"""Tiny synthetic contracts; these tests never open current biological sources."""
from copy import deepcopy
import json
from pathlib import Path
from types import SimpleNamespace

import pyarrow as pa
import pyarrow.parquet as pq
import pytest

from analysis.shared import protein_semantic_view as view


def write_json(path, value):
    path.write_text(json.dumps(value, indent=2) + "\n")


def bind(path):
    return {"path": str(path.resolve()), "sha256": view.base_view.sha256(path)}


def claim(name, level="qualified_broad_lineage_programme", programme=True,
          broad_lineage=True, structure=False):
    return {"annotation_name": name, "evidence_level": level,
            "reason": "Explicit synthetic reviewer decision",
            "eligibility": dict(programme=programme, broad_lineage=broad_lineage,
                                structure=structure)}


def refresh(fixture, *, accept=True):
    """Rebind deliberately edited synthetic compact inputs."""
    for key, value in fixture.values.items():
        write_json(fixture.paths[key], value)
        fixture.contract[key] = bind(fixture.paths[key])
    if accept:
        acceptance = {"status": "ROOT_ACCEPTED_PROTEIN_SEMANTIC_DECISIONS",
                      "bindings": {k: fixture.contract[k] for k in fixture.values},
                      "source_hold": deepcopy(fixture.contract["source_hold"])}
        write_json(fixture.acceptance, acceptance)
    fixture.contract["root_acceptance"] = bind(fixture.acceptance)
    write_json(fixture.path, fixture.contract)


@pytest.fixture
def fixture(tmp_path):
    scopes, tables = [], {}
    for lane, n in view.PROTEIN_SCOPES.items():
        for i in range(n):
            donor = "6523" if i == 0 else f"synthetic-{i}"
            region = "ordinary_pancreas" if lane == "pancreas_protein" else "synthetic-region"
            sid = f"{lane}|{donor}|{region}"
            path = tmp_path / f"{len(scopes)}.parquet"
            scope = {"scope_id": sid, "lane": lane, "donor_id": donor,
                     "path": str(path), "sha256": "", "cell_key_column": "object_id",
                     "claim_column": "candidate_cell_class" if lane == "pancreas_protein" else "candidate_lineage",
                     "QC_column": "qc_analysis_eligible", "QC_kind": "boolean", "filters": {},
                     "claim_status": "synthetic_unchanged_native_claim"}
            if lane == "pancreas_protein":
                tables[sid] = pa.table({"object_id": [f"same-{j}" for j in range(7)],
                    "candidate_cell_class": ["Blood vessel", "Immune", "Unknown", "Blood vessel",
                                             "Blood vessel", "Blood vessel", "Unspecified"],
                    "qc_analysis_eligible": [True, True, True, True, False, None, True],
                    "analysis_eligible": [True, True, True, False, False, None, True]})
                scope.update(additional_admission_column="analysis_eligible",
                             expected_additional_admission_rows=4)
            scopes.append(scope)
    pln_scopes = [s for s in scopes if s["lane"] != "pancreas_protein"]
    inventory, decisions, unit_keys = [], [], []
    rows = {s["scope_id"]: [] for s in pln_scopes}
    for i in range(609):
        s = pln_scopes[i % 31]
        prior = "Structural" if i < 17 else "T_cell" if i == 17 else "B_Plasma"
        key = dict(scope_id=s["scope_id"], cluster_id=f"c{i:03}", prior_candidate_lineage=prior)
        unit_keys.append(key)
        parts = [key]
        if i < 40:
            parts.append(dict(key, prior_candidate_lineage="T/NK_unresolved"))
        for part in parts:
            inventory.append(dict(part, current_native_QC_claim_rows=1))
            rows[s["scope_id"]].append({"object_id": "repeated-across-scopes",
                "candidate_lineage": part["prior_candidate_lineage"],
                "cluster_id": part["cluster_id"], "qc_analysis_eligible": True})
            if part["prior_candidate_lineage"] == "Structural":
                d = claim("Structure-associated pattern", "structural_or_contextual_pattern",
                          broad_lineage=False, structure=True)
            elif part["prior_candidate_lineage"] == "T/NK_unresolved":
                d = claim("T/NK unresolved", "identity_unresolved", programme=False, broad_lineage=False)
            else:
                d = claim("T-associated annotation" if prior == "T_cell" else "B-associated annotation")
            if i == 100:
                d = None
            if i != 101:  # Both explicit null and absent root decisions must remain unavailable.
                decisions.append(dict(part, decision=d))
    for s in pln_scopes:
        source = rows[s["scope_id"]] + [{"object_id": "failed", "candidate_lineage": "Excluded",
                                       "cluster_id": None, "qc_analysis_eligible": False}]
        tables[s["scope_id"]] = pa.Table.from_pylist(source)
    for s in scopes:
        t = tables[s["scope_id"]]
        pq.write_table(t, s["path"])
        s.update(sha256=view.base_view.sha256(s["path"]), expected_source_rows=len(t),
                 expected_QC_rows=sum(x is True for x in t["qc_analysis_eligible"].to_pylist()))
    registry = {"schema_version": 1, "scopes": scopes, "manifest_bindings": [],
                "selection_status": "synthetic", "across_scope_counts_are_additive": False}
    values = {"registry": registry, "PAN_route_dictionary": [
        {"prior_candidate_cell_class": "Blood vessel", "decision": claim("Vascular structure-associated",
            "structure_associated_allocation", broad_lineage=False, structure=True)},
        {"prior_candidate_cell_class": "Immune", "decision": claim("Broad immune annotation",
            "affirmative_broad_immune_allocation")},
        {"prior_candidate_cell_class": "Unknown", "decision": claim("Identity unresolved",
            "identity_unresolved", programme=False, broad_lineage=False)}],
        "pLN_unit_claim_inventory": inventory, "pLN_root_decisions": decisions}
    contract = {"schema_version": 1, "kind": "protein_semantic_overlay",
        "across_scope_counts_are_additive": False,
        "base_files": [{k: s[k] for k in ("scope_id", "path", "sha256")} for s in scopes],
        "source_hold": {"scope_id": view.HELD_SCOPE, "reason": "Synthetic unresolved source identity"}}
    f = SimpleNamespace(path=tmp_path/"overlay.json", acceptance=tmp_path/"acceptance.json",
        paths={k: tmp_path/(k+".json") for k in values}, values=values,
        contract=contract, tables=tables, scopes=scopes, units=unit_keys)
    refresh(f)
    return f


def test_compact_validation_does_not_read_base_files(fixture, monkeypatch):
    for s in fixture.scopes:
        Path(s["path"]).unlink()
    monkeypatch.setattr(view.base_view, "iter_scope_batches", lambda *a, **k: pytest.fail("base read"))
    summary = view.load_overlay(fixture.path).summary()
    assert summary["n_base_scopes"] == 56 and summary["n_pLN_inventory_parts"] == 649
    assert summary["n_pLN_root_decisions_unavailable"] == 2
    assert summary["base_cell_files_read_or_hashed"] is False


def test_exact_prior_parts_preserve_unresolved_and_structural_meanings(fixture):
    overlay = view.load_overlay(fixture.path)
    key = fixture.units[17]
    table = fixture.tables[key["scope_id"]]
    out = view.resolve_batch(overlay, key["scope_id"], table)
    assert out.select(table.column_names).equals(table)
    same_cluster = [r for r in out.to_pylist() if r["cluster_id"] == key["cluster_id"]]
    assert [(r["candidate_lineage"], r["eligible_for_broad_lineage"]) for r in same_cluster] == [
        ("T_cell", True), ("T/NK_unresolved", False)]
    structural = view.resolve_batch(overlay, fixture.units[0]["scope_id"],
                                   fixture.tables[fixture.units[0]["scope_id"]]).to_pylist()[0]
    assert structural["eligible_for_programme"] is True
    assert structural["eligible_for_structure"] is True
    assert structural["eligible_for_broad_lineage"] is False
    # Every inherited part survives independently, including in qualified units.
    retained = []
    for sid, scope in overlay.scopes.items():
        if scope["lane"] != "pancreas_protein":
            retained.extend(view.resolve_batch(overlay, sid, fixture.tables[sid]).to_pylist())
    tnks = [r for r in retained if r["candidate_lineage"] == "T/NK_unresolved"]
    structures = [r for r in retained if r["candidate_lineage"] == "Structural"]
    assert len(tnks) == 40 and len(structures) == 17
    assert all(r["eligible_for_broad_lineage"] is False for r in tnks + structures)


def test_missing_decisions_remain_nullable_not_negative_or_qc_identity(fixture):
    overlay = view.load_overlay(fixture.path)
    for i in (100, 101):
        k = fixture.units[i]
        out = view.resolve_batch(overlay, k["scope_id"], fixture.tables[k["scope_id"]])
        row = next(r for r in out.to_pylist() if r["cluster_id"] == k["cluster_id"])
        assert row["qc_analysis_eligible"] is True
        assert row["annotation_name"] is None
        assert row["eligible_for_broad_lineage"] is None
        assert row["semantic_decision_status"] == "root_decision_unavailable"
        assert row["eligible_for_native_provenance"] is True


def test_joint_purpose_counts_and_native_admissions_are_distinct(fixture):
    overlay = view.load_overlay(fixture.path)
    sid = fixture.scopes[1]["scope_id"]
    c = view.scope_counts(overlay, sid, purpose="structure", batch_size=2)
    assert (c["source_rows"], c["native_QC_rows"], c["native_additional_admission_rows"]) == (7, 5, 4)
    assert (c["purpose_eligible_rows"], c["purpose_ineligible_rows"], c["purpose_unavailable_rows"]) == (1, 4, 2)
    assert c["purpose_annotation_counts"] == [{"annotation_name": "Vascular structure-associated",
        "evidence_level": "structure_associated_allocation", "n": 1}]
    assert sum(r["n"] for r in c["joint_status_counts"]) == 7
    assert view.scope_counts(overlay, sid, purpose="broad_lineage")["purpose_eligible_rows"] == 1
    assert view.scope_counts(overlay, sid, purpose="native_provenance")["purpose_eligible_rows"] == 7


def test_hold_is_only_exact_ordinary_pan6523_not_dedicated_pln6523(fixture):
    overlay = view.load_overlay(fixture.path)
    held = view.scope_counts(overlay, view.HELD_SCOPE, purpose="programme")
    assert held["held_source_rows"] == 7 and held["purpose_eligible_rows"] == 0
    assert held["native_QC_rows"] == 5 and held["native_additional_admission_rows"] == 4
    pln_sid = fixture.units[0]["scope_id"]
    assert "|6523|" in pln_sid
    c = view.scope_counts(overlay, pln_sid, purpose="programme")
    assert c["held_source_rows"] == 0 and c["purpose_eligible_rows"] > 0


def test_same_uuid_and_cluster_in_another_scope_does_not_inherit_decision(fixture):
    overlay = view.load_overlay(fixture.path)
    k = fixture.units[17]
    other_sid = fixture.units[18]["scope_id"]
    table = pa.table({"object_id": ["repeated-across-scopes"], "candidate_lineage": ["T_cell"],
                      "cluster_id": [k["cluster_id"]], "qc_analysis_eligible": [True]})
    assert view.resolve_batch(overlay, other_sid, table)["annotation_name"].to_pylist() == [None]


def test_changed_decisions_rehashed_in_overlay_need_fresh_exact_acceptance(fixture):
    fixture.values["pLN_root_decisions"][0]["decision"]["annotation_name"] = "Changed semantic claim"
    refresh(fixture, accept=False)
    with pytest.raises(ValueError, match="does not endorse the exact"):
        view.load_overlay(fixture.path)


def test_explicit_null_purpose_and_all_null_pending_contract_remain_unavailable(fixture):
    fixture.values["PAN_route_dictionary"][1]["decision"]["eligibility"]["broad_lineage"] = None
    refresh(fixture)
    overlay = view.load_overlay(fixture.path)
    sid = fixture.scopes[1]["scope_id"]
    row = view.resolve_batch(overlay, sid, fixture.tables[sid]).to_pylist()[1]
    assert row["semantic_decision_status"] == "root_decision_present"
    assert row["eligible_for_programme"] is True and row["eligible_for_broad_lineage"] is None
    for key in ("PAN_route_dictionary", "pLN_root_decisions"):
        for record in fixture.values[key]:
            record["decision"] = None
    refresh(fixture)
    del fixture.contract["root_acceptance"]
    write_json(fixture.path, fixture.contract)
    pending = view.load_overlay(fixture.path)
    assert pending.summary()["n_pLN_root_decisions_unavailable"] == 649
    assert view.scope_counts(pending, sid, purpose="programme")["purpose_unavailable_rows"] == 7


def test_present_decision_cannot_omit_acceptance(fixture):
    del fixture.contract["root_acceptance"]
    write_json(fixture.path, fixture.contract)
    with pytest.raises(ValueError, match="file binding must be an object"):
        view.load_overlay(fixture.path)


@pytest.mark.parametrize("change,match", [
    ("unit", "exactly609"), ("count", "per-scope denominators"),
    ("key", "foreign pLN root"), ("pin", "Base pin differs"),
    ("missing_base", "All56"), ("hold", "Only the exact"),
    ("pattern_identity", "cannot admit a named lineage"),
    ("numeric_boolean", "true, false or null"),
])
def test_contract_rejects_wrong_scope_counts_pins_and_semantic_types(fixture, change, match):
    if change == "unit":
        # Preserve649parts but change one shared-unit part into an extra physical unit.
        fixture.values["pLN_unit_claim_inventory"][1]["cluster_id"] = "extra-unit"
    elif change == "count":
        fixture.values["pLN_unit_claim_inventory"][0]["current_native_QC_claim_rows"] = 2
    elif change == "key":
        fixture.values["pLN_root_decisions"][0]["cluster_id"] = "foreign-unit"
    elif change == "pin":
        fixture.contract["base_files"][0]["sha256"] = "0"*64
    elif change == "missing_base":
        fixture.contract["base_files"].pop()
    elif change == "hold":
        fixture.contract["source_hold"]["scope_id"] = fixture.units[0]["scope_id"]
    elif change == "pattern_identity":
        fixture.values["PAN_route_dictionary"][0]["decision"]["eligibility"]["broad_lineage"] = True
    elif change == "numeric_boolean":
        fixture.values["PAN_route_dictionary"][1]["decision"]["eligibility"]["broad_lineage"] = 1
    refresh(fixture)
    with pytest.raises(ValueError, match=match):
        view.load_overlay(fixture.path)


def test_acceptance_must_accept_exact_hold_and_decision_status(fixture):
    acceptance = json.loads(fixture.acceptance.read_text())
    acceptance["source_hold"]["reason"] = "Different hold"
    write_json(fixture.acceptance, acceptance)
    fixture.contract["root_acceptance"] = bind(fixture.acceptance)
    write_json(fixture.path, fixture.contract)
    with pytest.raises(ValueError, match="source hold differs"):
        view.load_overlay(fixture.path)
    acceptance["status"] = "PENDING_REVIEW"
    write_json(fixture.acceptance, acceptance)
    fixture.contract["root_acceptance"] = bind(fixture.acceptance)
    write_json(fixture.path, fixture.contract)
    with pytest.raises(ValueError, match="explicitly accept"):
        view.load_overlay(fixture.path)


def test_native_bool_types_and_admission_constraints_survive_overlay(fixture):
    overlay = view.load_overlay(fixture.path)
    sid = fixture.scopes[1]["scope_id"]
    table = fixture.tables[sid]
    bad = table.set_column(2, "qc_analysis_eligible", pa.array([1]*len(table)))
    with pytest.raises(ValueError, match="QC must remain boolean"):
        view.resolve_batch(overlay, sid, bad)
    bad = table.set_column(3, "analysis_eligible", pa.array([True]*len(table)))
    with pytest.raises(ValueError, match="cannot bypass native QC"):
        view.resolve_batch(overlay, sid, bad)
    derived = view.resolve_batch(overlay, sid, table)
    with pytest.raises(ValueError, match="overwrite existing semantic"):
        view.resolve_batch(overlay, sid, derived)


def test_selected_base_hash_change_rejected_only_on_explicit_count(fixture):
    overlay = view.load_overlay(fixture.path)
    sid = fixture.scopes[1]["scope_id"]
    pq.write_table(fixture.tables[sid].slice(0, 1), overlay.scopes[sid]["path"])
    with pytest.raises(ValueError, match="Selected cell table changed"):
        view.scope_counts(overlay, sid, purpose="programme")


def test_cli_requires_purpose_and_will_not_overwrite_receipt(fixture, monkeypatch, capsys, tmp_path):
    monkeypatch.setattr("sys.argv", ["semantic", "--overlay", str(fixture.path), "validate"])
    view.main()
    assert json.loads(capsys.readouterr().out)["selection_activated"] is False
    argv = ["semantic", "--overlay", str(fixture.path), "counts", "--scope-id", fixture.scopes[1]["scope_id"]]
    monkeypatch.setattr("sys.argv", argv)
    with pytest.raises(SystemExit):
        view.main()
    output = tmp_path/"new-counts.json"
    monkeypatch.setattr("sys.argv", argv + ["--purpose", "structure", "--output", str(output)])
    view.main()
    assert json.loads(output.read_text())["purpose_eligible_rows"] == 1
    with pytest.raises(FileExistsError):
        view.main()
