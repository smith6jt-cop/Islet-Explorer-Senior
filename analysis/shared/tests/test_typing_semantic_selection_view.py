"""Synthetic-only schema2 entries and independently accepted proof fixtures."""
from copy import deepcopy
import json
from pathlib import Path

import pyarrow as pa
import pyarrow.parquet as pq
import pytest

from analysis.shared import typing_semantic_selection_view as view
from analysis.shared.tests.test_protein_semantic_view import (
    fixture as protein_fixture, bind, claim, refresh, write_json,
)


def freeze_entry(f, *, accept=True):
    write_json(f.readback_path, f.readback)
    f.entry["protein_semantics"]["independent_full_row_readback"] = bind(f.readback_path)
    if accept:
        f.readback_acceptance = {"status": view.ACCEPTANCE_STATUS,
            "full56_semantics_accepted": True,
            "bindings": {k: {name: f.readback["bindings"][k][name] for name in ("path", "sha256")} for k in (
                "registry", "overlay_contract", "actual_consumer", "actual_consumer_base_helper", "root_acceptance")}}
        f.readback_acceptance["bindings"]["independent_full_row_readback"] = bind(f.readback_path)
        write_json(f.readback_acceptance_path, f.readback_acceptance)
    f.entry["root_readback_acceptance"] = bind(f.readback_acceptance_path)
    write_json(f.entry_path, f.entry)


@pytest.fixture
def semantic_fixture(protein_fixture, tmp_path):
    f = protein_fixture
    registry = f.values["registry"]
    for lane, n in (("pancreas_RNA", 25), ("pLN_RNA", 27)):
        for i in range(n):
            section = f"synthetic-{i}"
            sid = f"{lane}|serial-{i}|{section}"
            path = tmp_path/(sid.replace("|", "-")+".parquet")
            t = pa.table({"serial_qualified_cell_key": [f"{sid}|same", f"{sid}|other"],
                "section_id": [section, section], "candidate_lineage": ["B", "Unresolved"],
                "subcluster_id": ["c0", None], "candidate_lineage_analysis_eligible": [True, False]})
            pq.write_table(t, path)
            s = {"scope_id": sid, "lane": lane, "donor_id": str(i), "section_id": section,
                 "path": str(path), "sha256": bind(path)["sha256"],
                 "cell_key_column": "serial_qualified_cell_key", "claim_column": "candidate_lineage",
                 "QC_column": "subcluster_id", "QC_kind": "nonnull", "expected_source_rows": 2,
                 "expected_QC_rows": 1, "additional_admission_column": "candidate_lineage_analysis_eligible",
                 "expected_additional_admission_rows": 1, "filters": {},
                 "claim_status": "synthetic accepted RNA broad tier"}
            if lane == "pLN_RNA":
                s.update(storage_scope="section_specific_candidate", source_serial=f"serial-{i}",
                         filters={"section_id": section}, native_source_filters={"section_key": section})
            registry["scopes"].append(s)
            f.tables[sid] = t
    registry["current_RNA_lane_generations"] = {"pancreas": "r13", "pLN": "r10"}
    # The final schema2 gate needs an explicit decision for all649 parts.
    by_key = {tuple(r[k] for k in view.protein.KEY_FIELDS): r for r in f.values["pLN_root_decisions"]}
    for part in f.values["pLN_unit_claim_inventory"]:
        key = tuple(part[k] for k in view.protein.KEY_FIELDS)
        if key not in by_key:
            f.values["pLN_root_decisions"].append({**dict(zip(view.protein.KEY_FIELDS, key)),
                                                 "decision": claim("B-associated annotation")})
        elif by_key[key]["decision"] is None:
            by_key[key]["decision"] = claim("B-associated annotation")
    refresh(f)
    f.previous_path = tmp_path/"previous-selection-inputs.json"
    write_json(f.previous_path, {"registry": f.contract["registry"],
        "current_RNA_lane_generations": registry["current_RNA_lane_generations"]})
    f.readback_path, f.readback_acceptance_path = tmp_path/"readback.json", tmp_path/"root-readback-acceptance.json"
    f.entry_path = tmp_path/"schema2.json"
    f.entry = {"schema_version": 2, "kind": "typing_semantic_selection",
        "required_capabilities": sorted(view.CAPABILITIES),
        "base_registry": f.contract["registry"], "previous_selection_inputs": bind(f.previous_path),
        "consumer": bind(Path(view.__file__)), "scopes": deepcopy(registry["scopes"]),
        "current_RNA_lane_generations": deepcopy(registry["current_RNA_lane_generations"]),
        "across_scope_counts_are_additive": False, "biological_identity_accuracy_validated": False,
        "protein_semantics": {"overlay": bind(f.path), "semantic_api": bind(Path(view.protein.__file__)),
            "root_acceptance": bind(f.acceptance),
            "required_scientific_purposes": list(view.protein.SCIENTIFIC_PURPOSES),
            "explicit_native_provenance_purpose": "native_provenance"},
        "RNA_semantics": {"policy": "unchanged_selected_RNA_scope_contracts", "base_api": bind(Path(view.base.__file__))}}
    overlay = view.protein.load_overlay(f.path)
    f.readback = {"schema_version": 1, "status": view.READBACK_STATUS, "n_scopes": 56,
        "selection_gate": {k: True for k in view.READBACK_GATES},
        "n_derived_columns": 14, "expected_semantics_from_consumer_helpers": False,
        "bindings": {"registry": f.contract["registry"], "overlay_contract": bind(f.path),
            "actual_consumer": bind(Path(view.protein.__file__)),
            "actual_consumer_base_helper": bind(Path(view.base.__file__)), "root_acceptance": bind(f.acceptance),
            "extra_frozen_audit_source": {"path": "synthetic unused provenance", "sha256": "1"*64}},
        "coverage": {"n_pLN_units": 609, "n_pLN_prior_claim_parts": 649,
            "n_pLN_ordered_UUID_membership_hashes": 609,
            "ordered_UUID_membership_hashes_verified": True, "exact_inventory_parts_verified": True},
        "across_scope_counts_are_additive": False, "biological_identity_accuracy_validated": False,
        "scopes": []}
    # Actual reader binding() carries sizes; the entry/acceptance pins normalize
    # to path/SHA. Additional receipt metadata must not break the interface.
    for name in ("registry", "overlay_contract", "actual_consumer", "actual_consumer_base_helper", "root_acceptance"):
        f.readback["bindings"][name] = dict(f.readback["bindings"][name],
            size_bytes=Path(f.readback["bindings"][name]["path"]).stat().st_size)
    for sid, s in overlay.scopes.items():
        purposes = {p: view.protein.scope_counts(overlay, sid, purpose=p) for p in view.protein.PURPOSES}
        c = purposes["programme"]
        f.readback["scopes"].append({"scope_id": sid, "base_file": {"path": s["path"], "sha256": s["sha256"]},
            "source_rows": c["source_rows"], "native_QC_rows": c["native_QC_rows"],
            "native_additional_admission_rows": c["native_additional_admission_rows"],
            "original_row_order_preserved": True, "original_columns_values_preserved": True,
            "Arrow_schema_metadata_preserved": True,
            "purpose_totals": {p: {k: c["purpose_"+k] for k in
                                  ("eligible_rows", "ineligible_rows", "unavailable_rows")}
                               for p, c in purposes.items()}})
    f.readback["n_source_rows_nonadditive"] = sum(s["source_rows"] for s in f.readback["scopes"])
    f.readback["n_native_QC_observations_nonadditive"] = sum(s["native_QC_rows"] for s in f.readback["scopes"])
    freeze_entry(f)
    return f


def test_legacy_loader_rejects2_but_ignores_overlay_fields_in1(semantic_fixture):
    f = semantic_fixture
    with pytest.raises(ValueError, match="Unsupported typing registry schema"):
        view.base.load_registry(f.entry_path)
    legacy = deepcopy(f.values["registry"])
    legacy["protein_semantics"] = {"mandatory": "silently ignored by legacy loader"}
    path = f.entry_path.parent/"unsafe-v1-example.json"
    write_json(path, legacy)
    assert len(view.base.load_registry(path)["scopes"]) == 108


def test_compact_validate_reads_no_base_cells_and_metadata_is_not_legacy_scope(semantic_fixture, monkeypatch):
    f = semantic_fixture
    for scope in f.entry["scopes"]:
        Path(scope["path"]).unlink()
    monkeypatch.setattr(view.base, "iter_scope_batches", lambda *a, **k: pytest.fail("base cells opened"))
    selection = view.load_semantic_selection(f.entry_path)
    assert selection.summary()["base_cell_files_read_or_hashed"] is False
    metadata = selection.scope_metadata(f.scopes[1]["scope_id"])
    assert not isinstance(metadata, dict)
    with pytest.raises(TypeError):
        view.base.scope_counts(metadata)
    with pytest.raises(TypeError):
        metadata.fields["path"] = "mutated"


def test_protein_counts_require_purpose_and_never_call_legacy_counts(semantic_fixture, monkeypatch):
    f = semantic_fixture
    selection = view.load_semantic_selection(f.entry_path)
    monkeypatch.setattr(view.base, "scope_counts", lambda *a, **k: pytest.fail("legacy protein claim count"))
    sid = f.scopes[1]["scope_id"]
    with pytest.raises(ValueError, match="explicit supported purpose"):
        selection.counts(sid)
    assert selection.counts(sid, purpose="structure")["purpose_eligible_rows"] == 1
    assert selection.counts(sid, purpose="native_provenance")["purpose_eligible_rows"] == 7
    held = selection.counts(view.protein.HELD_SCOPE, purpose="programme")
    assert held["purpose_eligible_rows"] == 0 and held["native_QC_rows"] == 5


def test_wrapper_passes_bounded8192_batch_default_without_changing_pinned_API(semantic_fixture, monkeypatch):
    selection = view.load_semantic_selection(semantic_fixture.entry_path)
    sid = semantic_fixture.scopes[1]["scope_id"]
    def counts(*args, **kwargs):
        assert kwargs["batch_size"] == 8192
        return {"bounded": True}
    monkeypatch.setattr(view.protein, "scope_counts", counts)
    assert selection.counts(sid, purpose="programme") == {"bounded": True}
    def batches(*args, **kwargs):
        assert kwargs["batch_size"] == 8192
        yield from ()
    monkeypatch.setattr(view.protein, "iter_resolved_batches", batches)
    assert list(selection.iter_annotations(sid)) == []
    assert list(selection.iter_for_purpose(sid, purpose="programme")) == []


def test_all52_RNA_scope_dictionaries_and_counts_unchanged(semantic_fixture):
    f = semantic_fixture
    selection = view.load_semantic_selection(f.entry_path)
    RNA = [s for s in f.entry["scopes"] if "RNA" in s["lane"]]
    assert len(RNA) == 52
    for s in RNA:
        assert dict(selection.scope_metadata(s["scope_id"]).fields) == s
        assert selection.counts(s["scope_id"]) == view.base.scope_counts(s)
    with pytest.raises(ValueError, match="omit protein purpose"):
        selection.counts(RNA[0]["scope_id"], purpose="broad_lineage")


def test_annotation_stream_keeps_all_rows_and_explicit_filter_chooses_true_only(semantic_fixture):
    f = semantic_fixture
    selection = view.load_semantic_selection(f.entry_path)
    sid = f.scopes[1]["scope_id"]
    original = f.tables[sid]
    all_rows = pa.concat_tables(list(selection.iter_annotations(sid, batch_size=2)))
    assert all_rows.select(original.column_names).equals(original)
    assert all_rows["eligible_for_structure"].to_pylist() == [True, False, False, False, False, None, None]
    selected = pa.concat_tables(list(selection.iter_for_purpose(sid, purpose="structure", batch_size=2)))
    assert len(selected) == 1 and selected["eligible_for_structure"].to_pylist() == [True]
    assert len(pa.concat_tables(list(selection.iter_for_purpose(sid, purpose="native_provenance")))) == 7


@pytest.mark.parametrize("change,match", [
    ("schema1", "schema2"), ("missing_capability", "capability"),
    ("unknown_capability", "capability"), ("missing_overlay", "file binding"),
    ("consumer_hash", "Changed semantic input"), ("RNA_dict", "scope dictionaries"),
    ("protein_dict", "scope dictionaries"), ("scope_order", "scope dictionaries"),
    ("purpose_default", "purposes"), ("RNA_policy", "RNA semantics"),
])
def test_entry_cannot_silently_drop_semantics_or_change_native_scopes(semantic_fixture, change, match):
    f = semantic_fixture
    if change == "schema1": f.entry["schema_version"] = 1
    elif change == "missing_capability": f.entry["required_capabilities"].pop()
    elif change == "unknown_capability": f.entry["required_capabilities"][0] = "raw_QC_equals_identity"
    elif change == "missing_overlay": f.entry["protein_semantics"].pop("overlay")
    elif change == "consumer_hash": f.entry["consumer"]["sha256"] = "0"*64
    elif change == "RNA_dict": f.entry["scopes"][-1]["expected_QC_rows"] += 1
    elif change == "protein_dict": f.entry["scopes"][0]["additional_admission_semantics"] = "altered"
    elif change == "scope_order": f.entry["scopes"].reverse()
    elif change == "purpose_default": f.entry["protein_semantics"]["required_scientific_purposes"] = ["known_lineage"]
    elif change == "RNA_policy": f.entry["RNA_semantics"]["policy"] = "reinterpret"
    write_json(f.entry_path, f.entry)
    with pytest.raises(ValueError, match=match): view.load_semantic_selection(f.entry_path)


@pytest.mark.parametrize("change,match", [
    ("completion", "successful independent"), ("missing_scope", "all56"),
    ("source_pin", "source file differs"), ("row_order", "preserve original"),
    ("column_values", "preserve original"), ("Arrow_metadata", "preserve original"),
    ("ordered_UUID", "UUID-membership"), ("parts", "UUID-membership"),
    ("purpose_sum", "partition retained"), ("native_provenance", "Native provenance"),
    ("native_QC", "denominator mismatch"), ("native_admission", "admission denominator"),
    ("source_total", "source row total"), ("QC_total", "QC observation total"),
    ("semantic_pin", "bindings differ"),
    ("false_gate", "operational selection gate"),
    ("missing_derived_gate", "operational selection gate"),
    ("consumer_expected", "all14 derived"),
    ("missing_derived_count", "all14 derived"),
])
def test_readback_requires_full_original_preservation_coverage_and_joint_counts(semantic_fixture, change, match):
    f = semantic_fixture
    row = f.readback["scopes"][0]
    if change == "completion": f.readback["status"] = "COMPLETED"
    elif change == "missing_scope": f.readback["scopes"].pop()
    elif change == "source_pin": row["base_file"]["sha256"] = "0"*64
    elif change == "row_order": row["original_row_order_preserved"] = False
    elif change == "column_values": row["original_columns_values_preserved"] = False
    elif change == "Arrow_metadata": row["Arrow_schema_metadata_preserved"] = False
    elif change == "ordered_UUID": f.readback["coverage"]["ordered_UUID_membership_hashes_verified"] = False
    elif change == "parts": f.readback["coverage"]["n_pLN_prior_claim_parts"] = 648
    elif change == "purpose_sum": row["purpose_totals"]["programme"]["unavailable_rows"] += 1
    elif change == "native_provenance": row["purpose_totals"]["native_provenance"].update(eligible_rows=6,ineligible_rows=1)
    elif change == "native_QC": row["native_QC_rows"] += 1
    elif change == "native_admission": row["native_additional_admission_rows"] += 1
    elif change == "source_total": f.readback["n_source_rows_nonadditive"] += 1
    elif change == "QC_total": f.readback["n_native_QC_observations_nonadditive"] += 1
    elif change == "semantic_pin": f.readback["bindings"]["actual_consumer"]["sha256"] = "0"*64
    elif change == "false_gate": f.readback["selection_gate"]["all_native_QC_admissions_and_exclusions_preserved"] = False
    elif change == "missing_derived_gate": f.readback["selection_gate"].pop("all_derived_values_independently_reconstructed_and_exact")
    elif change == "consumer_expected": f.readback["expected_semantics_from_consumer_helpers"] = True
    elif change == "missing_derived_count": f.readback.pop("n_derived_columns")
    freeze_entry(f)
    with pytest.raises(ValueError, match=match): view.load_semantic_selection(f.entry_path)


def test_rehashing_readback_does_not_reuse_root_acceptance(semantic_fixture):
    f = semantic_fixture
    f.readback["review_note"] = "Changed frozen receipt"
    freeze_entry(f, accept=False)
    with pytest.raises(ValueError, match="does not endorse these exact"):
        view.load_semantic_selection(f.entry_path)


def test_receipt_all5_binding_sizes_are_allowed_but_root_pins_stay_normalized(semantic_fixture):
    f = semantic_fixture
    names = ("registry", "overlay_contract", "actual_consumer", "actual_consumer_base_helper", "root_acceptance")
    assert all("size_bytes" in f.readback["bindings"][name] for name in names)
    assert all(set(f.readback_acceptance["bindings"][name]) == {"path", "sha256"} for name in names)
    assert view.load_semantic_selection(f.entry_path).summary()["n_scopes"] == 108


@pytest.mark.parametrize("field,value", [("path", "/wrong/reader/source.py"), ("sha256", "0"*64)])
def test_receipt_binding_size_does_not_mask_wrong_path_or_SHA(semantic_fixture, field, value):
    f = semantic_fixture
    f.readback["bindings"]["actual_consumer"][field] = value
    freeze_entry(f)
    with pytest.raises(ValueError, match="bindings differ"):
        view.load_semantic_selection(f.entry_path)


def test_root_acceptance_requires_full56_semantics_not_completion(semantic_fixture):
    f = semantic_fixture
    f.readback_acceptance["full56_semantics_accepted"] = False
    write_json(f.readback_acceptance_path, f.readback_acceptance)
    f.entry["root_readback_acceptance"] = bind(f.readback_acceptance_path)
    write_json(f.entry_path, f.entry)
    with pytest.raises(ValueError, match="explicitly accept independent full56 semantics"):
        view.load_semantic_selection(f.entry_path)


def test_final_selection_rejects_pending649_root_decisions(semantic_fixture):
    f = semantic_fixture
    f.values["pLN_root_decisions"][0]["decision"] = None
    refresh(f)
    f.entry["base_registry"] = f.contract["registry"]
    f.entry["protein_semantics"]["overlay"] = bind(f.path)
    f.entry["protein_semantics"]["root_acceptance"] = bind(f.acceptance)
    write_json(f.entry_path, f.entry)
    with pytest.raises(ValueError, match="all649"):
        view.load_semantic_selection(f.entry_path)


def test_cli_validation_and_no_default_protein_or_overwrite(semantic_fixture, monkeypatch, capsys, tmp_path):
    f = semantic_fixture
    argv = ["semantic-selection", "--selection", str(f.entry_path)]
    monkeypatch.setattr("sys.argv", argv + ["validate"])
    view.main()
    assert json.loads(capsys.readouterr().out)["selection_activated"] is False
    counts = argv + ["counts", "--scope-id", f.scopes[1]["scope_id"]]
    monkeypatch.setattr("sys.argv", counts)
    with pytest.raises(ValueError, match="explicit supported purpose"): view.main()
    output = tmp_path/"new-counts.json"
    monkeypatch.setattr("sys.argv", counts + ["--purpose", "structure", "--output", str(output)])
    view.main()
    assert json.loads(output.read_text())["purpose_eligible_rows"] == 1
    with pytest.raises(FileExistsError): view.main()
