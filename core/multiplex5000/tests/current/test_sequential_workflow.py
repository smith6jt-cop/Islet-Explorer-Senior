import json
from pathlib import Path

import pandas as pd
import pytest
from shapely.geometry import box

from multiplex5000.integration.sequential_workflow import (
    POLICY, SequentialError, digest, inspect_bundle, verify_bundle, replay_assignment,
)
from multiplex5000.integration.serial_context import geometry_descriptor, assign_multicell_geometry


def write(path, value):
    path.write_text(json.dumps(value, indent=2) + "\n")


def manifest(root):
    write(root / "artifact_manifest.json", [dict(path=str(p.relative_to(root)), bytes=p.stat().st_size,
                                               sha256=digest(p)) for p in sorted(root.rglob("*"))
                                            if p.is_file() and p.name != "artifact_manifest.json"])


@pytest.fixture
def bundle(tmp_path):
    root = tmp_path / "bundle with spaces"; root.mkdir()
    results = root / "results" / "D1"; results.mkdir(parents=True)
    native = tmp_path / "native"; native.mkdir()
    specs = []; refs = []; tables = {}
    for role in ("RNA", "protein"):
        p = native / role; p.mkdir()
        rows = []
        for suffix, cls, n, geometry in (("multi", "Islet", 2, box(0, 0, 40, 40)),
                                       ("single", "EndocrineSingle", 1, box(100, 100, 110, 110))):
            rows.append(dict(structure_id=role + suffix, object_class=cls, n_whole_cells=n,
                             geometry_wkb=geometry.wkb, role=role, donor_id="D1", proposal_evidence_available=True,
                             coordinate_frame=role + ":um", **geometry_descriptor(geometry)))
        table = pd.DataFrame(rows); tables[role] = table
        table.to_parquet(p / "objects.parquet", index=False)
        pd.DataFrame([dict(cell_id=role + str(i), structure_id=role + ("multi" if i < 2 else "single"))
                      for i in range(3)]).to_parquet(p / "object_members.parquet", index=False)
        qc = dict(n_objects=2, n_multiply_owned_cells=0, n_cut_cells=0, n_construction_count_disagreements=0,
                  n_single_geometry_disagreements=0, n_single_constructor_geometry_disagreements=0)
        write(p / "verified.json", qc); write(p / "image_qc.json", {"complete": True})
        spec = dict(role=role, donor_id="D1", source_serial="00001", acquisition_id=role + ":join", output=str(p))
        if role == "protein":
            transform = native / "transform.json"
            write(transform, {"affine": [[1, 0, 0], [0, 1, 0], [0, 0, 1]], "has_displacement": False})
            spec.update(transform_path=str(transform), transform_sha256=digest(transform), historical_verdict="PASS")
        specs.append(spec)
        refs.append(dict(role=role, donor_id="D1", source_serial="00001", output=str(p), n_objects=2,
                         object_sha256=digest(p / "objects.parquet"), membership_sha256=digest(p / "object_members.parquet"),
                         native_verification_sha256=digest(p / "verified.json"), image_review_sha256=digest(p / "image_qc.json")))
    for weight in (.25, .5, .75):
        assignment = assign_multicell_geometry(tables["RNA"], tables["protein"], position_weight=weight)
        if weight == .25:
            assignment["historical_anatomy_annotation"] = "not recomputed in assignment replay"
        assignment.to_parquet(results / f"assignment_position_weight_{weight}.parquet", index=False)
    tables["RNA"].iloc[:1].to_parquet(results / "RNA_registered_objects.parquet", index=False)
    pairs = assign_multicell_geometry(tables["RNA"], tables["protein"], position_weight=.25)
    pairs = pairs.assign(donor_id="D1", RNA_object_class="Islet", protein_object_class="Islet",
                         individual_correspondence_validated=False, same_cell_identity_claim=False,
                         correspondence_policy=POLICY, eligible_for_exploratory_regional_review=True,
                         anatomy_alternative_review_flag=False, ambiguous_geometry=False)
    for name in ("all_pairs.parquet", "regional_review_candidates.parquet"):
        pairs.to_parquet(root / name, index=False)
    pairs.to_parquet(results / "pairs.parquet", index=False)
    pd.concat([tables[role].iloc[:1] for role in tables]).to_parquet(root / "review_candidate_objects.parquet", index=False)
    write(root / "analysis_plan.json", {"correspondence_policy": POLICY})
    write(root / "summary.json", dict(correspondence_policy=POLICY, n_multicell_candidates=1,
                                     n_final_review_candidates=1, n_single_involving_pairs=0,
                                     n_native_objects=4, individual_pairs_validated=False, production_activated=False))
    write(root / "run_specs.json", specs); write(root / "native_output_reuse.json", refs)
    write(root / "final_verification.json", {"status": "passed"})
    write(root / "historical_preservation_reference.json", [])
    write(root / "release_decisions.json", [dict(donor_id="D1", n_pairs=1, n_final_review_candidates=1,
          correspondence_policy=POLICY, accepted_RNA_categories_for_regional_review=["Islet"],
          release_accepted_for_regional_review=True, production_accepted=False, individual_pairs_validated=False,
          pair_sha256=digest(results / "pairs.parquet"))])
    (root / "REPORT.html").write_text("original report")
    manifest(root)
    return root


def test_inspection_and_verification_never_change_sources(bundle):
    before = {str(p): (digest(p), p.stat().st_mtime_ns) for p in bundle.rglob("*") if p.is_file()}
    result = inspect_bundle(bundle)
    assert result["native_constructions"] == 2 and result["hashes_verified"] is False
    verified = verify_bundle(bundle)
    assert verified["counts"]["n_candidates"] == 1
    assert verified["external_files_hash_verified"] == 9
    assert verified["individual_pairs_validated"] is False
    assert {str(p): (digest(p), p.stat().st_mtime_ns) for p in bundle.rglob("*") if p.is_file()} == before


def test_undocumented_presentation_change_is_rejected(bundle):
    (bundle / "REPORT.html").write_text("unrecorded change")
    with pytest.raises(SequentialError, match="Undocumented artifact change"):
        verify_bundle(bundle)


def test_hash_bound_notice_is_reported_and_original_preserved(bundle):
    original = digest(bundle / "REPORT.html")
    (bundle / "REPORT.before_final_WARN_adjudication.html").write_text("original report")
    (bundle / "REPORT.html").write_text("original report with dated notice")
    write(bundle / "FINAL_WARN_RELEASE.json", dict(old_report_original_sha256=original,
          old_report_with_notice_sha256=digest(bundle / "REPORT.html"),
          old_all_pairs_sha256=digest(bundle / "all_pairs.parquet")))
    result = verify_bundle(bundle)
    assert len(result["documented_presentation_amendments"]) == 1
    (bundle / "REPORT.before_final_WARN_adjudication.html").write_text("changed backup")
    with pytest.raises(SequentialError, match="Pinned file changed"):
        verify_bundle(bundle)


def test_single_rows_cannot_be_hidden_by_a_new_artifact_manifest(bundle):
    q = pd.read_parquet(bundle / "all_pairs.parquet")
    q["RNA_object_class"] = "EndocrineSingle"
    q.to_parquet(bundle / "all_pairs.parquet", index=False); manifest(bundle)
    with pytest.raises(SequentialError, match="Single or unknown category"):
        verify_bundle(bundle)


def test_review_membership_is_checked_beyond_file_hashes(bundle):
    q = pd.read_parquet(bundle / "regional_review_candidates.parquet")
    q["moving_id"] = "different_native_object"
    q.to_parquet(bundle / "regional_review_candidates.parquet", index=False); manifest(bundle)
    with pytest.raises(SequentialError, match="exact eligible subset"):
        verify_bundle(bundle)


def test_manifest_cannot_escape_bundle(bundle):
    rows = json.loads((bundle / "artifact_manifest.json").read_text())
    rows.append(dict(path="../outside", bytes=1, sha256="bad")); write(bundle / "artifact_manifest.json", rows)
    with pytest.raises(SequentialError, match="Unsafe bundle manifest path"):
        inspect_bundle(bundle)


def test_bounded_replay_reuses_native_singles_without_pairing_them(bundle, tmp_path):
    output = tmp_path / "replay with spaces"
    result = replay_assignment(bundle, output, donor="D1")
    assert result["native_objects_reused"] == {"RNA": 2, "protein": 2}
    assert result["RNA_multicell_geometries_verified"] == 1
    assert all(r["n_pairs"] == 1 and r["exact_membership_and_values"] for r in result["assignment_checks"])
    assert result["assignment_checks"][0]["downstream_annotation_columns_not_recomputed"] == ["historical_anatomy_annotation"]
    assert result["native_construction_executed"] is False
    with pytest.raises(SequentialError, match="Output already exists"):
        replay_assignment(bundle, output, donor="D1")


def test_native_change_is_rejected_before_replay_output_creation(bundle, tmp_path):
    specs = json.loads((bundle / "run_specs.json").read_text())
    Path(specs[0]["output"], "object_members.parquet").write_bytes(b"changed")
    output = tmp_path / "must not exist"
    with pytest.raises(SequentialError, match="Pinned file changed"):
        replay_assignment(bundle, output, donor="D1")
    assert not output.exists()
