"""Read, verify and replay a completed sequential multicell correspondence bundle.

This interface reuses native constructions. It never executes QuPath, changes
registration, accepts scientific results, or modifies an existing result bundle.
The legacy ``acquisition_id`` key is preserved as a source join identifier; it is
not interpreted as a count of runs, slides, sections or regions.
"""
from __future__ import annotations

import hashlib
import json
from datetime import datetime, timezone
from pathlib import Path

POLICY = "multicell_only_no_single_identity_v1"
CLASSES = frozenset(("Islet", "EndocrineCluster"))
WEIGHTS = (0.25, 0.5, 0.75)


class SequentialError(ValueError):
    """A bundle cannot safely be reused under its recorded contract."""


def digest(path):
    h = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def _json(path):
    try:
        return json.loads(Path(path).read_text())
    except (OSError, ValueError) as error:
        raise SequentialError(f"Cannot read JSON {path}: {error}") from error


def _require(condition, message):
    if not condition:
        raise SequentialError(message)


def _check_hash(path, expected):
    try:
        actual = digest(path)
    except OSError as error:
        raise SequentialError(f"Cannot read pinned file {path}: {error}") from error
    _require(actual == expected, f"Pinned file changed: {path} (expected {expected}, found {actual})")
    return actual


def _bundle_path(root, relative):
    p = Path(relative)
    _require(not p.is_absolute() and ".." not in p.parts and p.parts,
             f"Unsafe bundle manifest path: {relative}")
    result = (root / p).resolve()
    _require(result.is_relative_to(root), f"Manifest path escapes bundle: {relative}")
    return result


def _key(spec):
    return (str(spec["role"]), str(spec["donor_id"]), str(spec["source_serial"]))


def _read(root):
    root = Path(root).expanduser().resolve()
    plan = _json(root / "analysis_plan.json")
    _require(plan.get("correspondence_policy") == POLICY,
             "Select the corrected multicell bundle; its explicit correspondence policy is missing or incompatible")
    summary = _json(root / "summary.json")
    _require(summary.get("correspondence_policy") == POLICY, "Summary policy differs from analysis plan")
    _require(not summary.get("production_activated") and not summary.get("individual_pairs_validated"),
             "This reuse workflow expects provisional correspondence, not an accepted production release")
    specs = _json(root / "run_specs.json")
    keys = [_key(s) for s in specs]
    _require(len(keys) == len(set(keys)), "Duplicate native section/modality keys")
    _require(all(s["role"] in ("RNA", "protein") for s in specs), "Unknown native modality")
    manifest = _json(root / "artifact_manifest.json")
    names = [row["path"] for row in manifest]
    _require(len(names) == len(set(names)), "Duplicate artifact manifest entries")
    resolved = [_bundle_path(root, name) for name in names]
    _require(len(resolved) == len(set(resolved)), "Artifact manifest paths alias one another")
    required = {"analysis_plan.json", "summary.json", "run_specs.json", "native_output_reuse.json",
                "final_verification.json", "all_pairs.parquet", "regional_review_candidates.parquet",
                "release_decisions.json", "review_candidate_objects.parquet"}
    _require(required <= set(names), f"Manifest lacks required artifacts: {sorted(required - set(names))}")
    return root, plan, summary, specs, manifest


def inspect_bundle(bundle):
    """Describe bindings and availability without hashing files or writing output."""
    root, plan, summary, specs, manifest = _read(bundle)
    return {
        "schema": "sequential_bundle_inspection_v1", "bundle": str(root),
        "correspondence_policy": POLICY, "summary": summary,
        "native_constructions": len(specs),
        "RNA_sections": sum(s["role"] == "RNA" for s in specs),
        "protein_sections": sum(s["role"] == "protein" for s in specs),
        "artifact_count": len(manifest), "artifact_bytes": sum(r["bytes"] for r in manifest),
        "missing_artifacts": [r["path"] for r in manifest if not _bundle_path(root, r["path"]).is_file()],
        "native_sources": [{"role": s["role"], "donor_id": s["donor_id"],
                            "source_serial": s["source_serial"], "region": s.get("roi"),
                            "legacy_source_join_id": s.get("acquisition_id"),
                            "path": s["output"], "available": Path(s["output"]).is_dir()} for s in specs],
        "hashes_verified": False, "individual_pairs_validated": False,
        "production_activated": False,
    }


def _report_amendment(root, record, actual):
    """Recognize the explicitly hash-bound historical WARN notice only."""
    if record["path"] != "REPORT.html" or not (root / "FINAL_WARN_RELEASE.json").is_file():
        return None
    path = root / "FINAL_WARN_RELEASE.json"
    receipt = _json(path)
    if (receipt.get("old_report_original_sha256") != record["sha256"] or
            receipt.get("old_report_with_notice_sha256") != actual):
        return None
    backup = root / "REPORT.before_final_WARN_adjudication.html"
    _check_hash(backup, record["sha256"])
    _check_hash(root / "all_pairs.parquet", receipt["old_all_pairs_sha256"])
    return {"path": record["path"], "manifest_sha256": record["sha256"], "current_sha256": actual,
            "amendment_receipt": str(path), "amendment_receipt_sha256": digest(path),
            "original_backup": str(backup),
            "scope": "Historical presentation notice only; does not identify the latest downstream scientific adjudication"}


def _pairs_contract(root, summary):
    import pandas as pd
    pairs = pd.read_parquet(root / "all_pairs.parquet")
    review = pd.read_parquet(root / "regional_review_candidates.parquet")
    decisions = _json(root / "release_decisions.json")
    keys = ["donor_id", "moving_id", "fixed_id"]
    for name, table in (("all pairs", pairs), ("review candidates", review)):
        _require(not table[keys].isna().any().any(), f"Missing native pair identifiers in {name}")
        _require(not table.duplicated(keys).any(), f"Duplicate pair identities in {name}")
        _require(not table.duplicated(["donor_id", "moving_id"]).any() and
                 not table.duplicated(["donor_id", "fixed_id"]).any(), f"Not one-to-one within donor: {name}")
        _require(table.RNA_object_class.isin(CLASSES).all() and table.protein_object_class.isin(CLASSES).all(),
                 f"Single or unknown category in {name}")
        _require(not table.individual_correspondence_validated.any() and not table.same_cell_identity_claim.any(),
                 f"Unexpected individual identity claim in {name}")
        _require(table.correspondence_policy.eq(POLICY).all(), f"Incompatible pair policy in {name}")
    _require(len(pairs) == summary["n_multicell_candidates"], "Candidate count differs from summary")
    _require(len(review) == summary["n_final_review_candidates"], "Review count differs from summary")
    _require(summary["n_single_involving_pairs"] == 0, "Summary contains Single correspondence")
    expected = pairs[pairs.eligible_for_exploratory_regional_review].reset_index(drop=True)
    try:
        pd.testing.assert_frame_equal(expected, review.reset_index(drop=True), check_exact=True)
    except AssertionError as error:
        raise SequentialError(f"Review table differs from exact eligible subset: {error}") from error
    _require(not review.anatomy_alternative_review_flag.any() and not review.ambiguous_geometry.any(),
             "Flagged ambiguity entered final exploratory review subset")
    _require(all(d.get("correspondence_policy") == POLICY and not d.get("production_accepted") and
                 not d.get("individual_pairs_validated") for d in decisions), "Decision policy or acceptance state is incompatible")
    _require(len({d["donor_id"] for d in decisions}) == len(decisions), "Duplicate donor decision")
    for d in decisions:
        donor = str(d["donor_id"])
        rows = pairs[pairs.donor_id.astype(str).eq(donor)]
        chosen = review[review.donor_id.astype(str).eq(donor)]
        _require(len(rows) == d["n_pairs"] and len(chosen) == d["n_final_review_candidates"],
                 f"Donor counts differ from decision: {donor}")
        _require(chosen.RNA_object_class.isin(d["accepted_RNA_categories_for_regional_review"]).all(),
                 f"Review category lacks regional support: {donor}")
        _require(d["release_accepted_for_regional_review"] or chosen.empty,
                 f"Diagnostic-only donor entered regional review: {donor}")
        _check_hash(root / "results" / donor / "pairs.parquet", d["pair_sha256"])
    objects = pd.read_parquet(root / "review_candidate_objects.parquet")
    _require(len(objects) == 2 * len(review) and objects.object_class.isin(CLASSES).all(),
             "Review object count or class differs from paired membership")
    for role, column in (("RNA", "moving_id"), ("protein", "fixed_id")):
        actual = objects[objects.role.eq(role)]
        _require(set(zip(actual.donor_id.astype(str), actual.structure_id)) ==
                 set(zip(review.donor_id.astype(str), review[column])), f"Review {role} geometry membership differs")
    return {"n_candidates": len(pairs), "n_final_review_candidates": len(review),
            "n_review_geometries": len(objects), "n_single_involving_pairs": 0, "n_donor_decisions": len(decisions)}


def verify_bundle(bundle, *, verify_native=True):
    """Hash all bundled artifacts and validate reusable scientific contracts.

    Native verification hashes all 39 recorded constructions for the current
    bundle. It does not reread raw assay matrices or repeat construction. Exact
    membership and transformed geometry comparisons are performed by the bounded
    ``replay_assignment`` action for the explicitly selected donor.
    """
    root, plan, summary, specs, manifest = _read(bundle)
    amendments = []
    for row in manifest:
        path = _bundle_path(root, row["path"])
        try:
            actual = digest(path)
        except OSError as error:
            raise SequentialError(f"Cannot read bundle artifact {path}: {error}") from error
        if actual != row["sha256"]:
            amendment = _report_amendment(root, row, actual)
            _require(amendment is not None, f"Undocumented artifact change: {path}")
            amendments.append(amendment)
        else:
            _require(path.stat().st_size == row["bytes"], f"Artifact size differs: {path}")
    counts = _pairs_contract(root, summary)
    refs = _json(root / "native_output_reuse.json")
    _require({_key(x) for x in refs} == {_key(s) for s in specs} and len(refs) == len(specs),
             "Native reuse roster differs from run specifications")
    by_key = {_key(s): s for s in specs}
    native_hashes = 0
    if verify_native:
        for reference in refs:
            spec = by_key[_key(reference)]
            _require(Path(reference["output"]).resolve() == Path(spec["output"]).resolve(), "Native output binding differs")
            p = Path(reference["output"])
            for filename, field in (("objects.parquet", "object_sha256"),
                                    ("object_members.parquet", "membership_sha256"),
                                    ("verified.json", "native_verification_sha256"),
                                    ("image_qc.json", "image_review_sha256")):
                _check_hash(p / filename, reference[field]); native_hashes += 1
            qc = _json(p / "verified.json")
            _require(qc["n_objects"] == reference["n_objects"], "Native object counts differ")
            for key in ("n_multiply_owned_cells", "n_cut_cells", "n_construction_count_disagreements",
                        "n_single_geometry_disagreements", "n_single_constructor_geometry_disagreements"):
                _require(qc[key] == 0, f"Native construction failed {key}: {p}")
        for s in specs:
            if s["role"] == "protein":
                _require(s["historical_verdict"] in ("PASS", "WARN"), "Failed or unreviewed transform entered correspondence")
                _check_hash(s["transform_path"], s["transform_sha256"]); native_hashes += 1
        historical = _json(root / "historical_preservation_reference.json")
        for row in historical:
            _check_hash(row["path"], row["sha256"]); native_hashes += 1
    _require(sum(r["n_objects"] for r in refs) == summary["n_native_objects"], "Native total differs from summary")
    from multiplex5000.runtime_context import runtime_context
    return {"schema": "sequential_bundle_verification_v1", "status": "passed",
            "verified_utc": datetime.now(timezone.utc).isoformat(), "bundle": str(root),
            "artifact_manifest_sha256": digest(root / "artifact_manifest.json"),
            "correspondence_policy": POLICY, "artifact_count": len(manifest),
            "documented_presentation_amendments": amendments, "counts": counts,
            "native_sources_hash_verified": bool(verify_native), "external_files_hash_verified": native_hashes,
            "native_constructions": len(refs), "native_objects": sum(r["n_objects"] for r in refs),
            "workflow_source_file": str(Path(__file__).resolve()), "workflow_source_sha256": digest(__file__),
            "runtime": runtime_context(), "individual_pairs_validated": False, "production_activated": False}


def replay_assignment(bundle, output, *, donor, position_weights=WEIGHTS):
    """Recompute a bounded donor's assignments without changing source results."""
    import numpy as np
    import pandas as pd
    import shapely
    from shapely.affinity import affine_transform
    from .serial_context import geometry_descriptor, assign_multicell_geometry, multicell_objects

    root, plan, summary, specs, manifest = _read(bundle)
    output = Path(output).expanduser().resolve()
    _require(not output.exists(), f"Output already exists: {output}")
    _require(not output.is_relative_to(root), "Replay output must be outside the preserved bundle")
    _require(tuple(position_weights) == WEIGHTS, "Replay uses the three recorded position weights: 0.25, 0.5, 0.75")
    donor = str(donor)
    selected = [s for s in specs if str(s["donor_id"]) == donor]
    rna = [s for s in selected if s["role"] == "RNA"]
    protein = [s for s in selected if s["role"] == "protein"]
    _require(len(rna) == len(protein) == 1, f"Donor {donor} needs one unambiguous RNA/protein section pair")
    verification = verify_bundle(root, verify_native=False)
    refs = {_key(r): r for r in _json(root / "native_output_reuse.json")}
    native = {}
    for spec in (rna[0], protein[0]):
        p = Path(spec["output"]); ref = refs[_key(spec)]
        for filename, field in (("objects.parquet", "object_sha256"), ("object_members.parquet", "membership_sha256"),
                                ("verified.json", "native_verification_sha256"), ("image_qc.json", "image_review_sha256")):
            _check_hash(p / filename, ref[field])
        objects = pd.read_parquet(p / "objects.parquet")
        members = pd.read_parquet(p / "object_members.parquet")
        _require(not objects.structure_id.duplicated().any(), "Duplicate native object identifiers")
        _require(not members.cell_id.duplicated().any(), "A native cell belongs to more than one object")
        _require(members.structure_id.isin(objects.structure_id).all(), "Unknown object in native membership")
        counts = members.groupby("structure_id").size().reindex(objects.structure_id, fill_value=0).to_numpy()
        _require(np.array_equal(counts, objects.n_whole_cells.to_numpy()), "Native whole-cell memberships disagree with object counts")
        native[spec["role"]] = objects
    ps = protein[0]
    _require(ps["historical_verdict"] in ("PASS", "WARN"), "Replay requires a retained PASS/WARN transform")
    _check_hash(ps["transform_path"], ps["transform_sha256"])
    transform = _json(ps["transform_path"])
    _require(not transform.get("has_displacement"), "Replay requires the recorded affine transform")
    matrix = np.asarray(transform["affine"])
    params = [matrix[0, 0], matrix[0, 1], matrix[1, 0], matrix[1, 1], matrix[0, 2], matrix[1, 2]]
    a = native["RNA"]
    a = multicell_objects(a[a.proposal_evidence_available]).copy()
    b = multicell_objects(native["protein"]).copy()
    a["geometry"] = [affine_transform(g, params) for g in shapely.from_wkb(a.geometry_wkb)]
    descriptor = pd.DataFrame([geometry_descriptor(g) for g in a.geometry])
    for field in descriptor:
        a[field] = descriptor[field].to_numpy()
    a = a.sort_values("structure_id").reset_index(drop=True)
    b = b.sort_values("structure_id").reset_index(drop=True)
    frozen = pd.read_parquet(root / "results" / donor / "RNA_registered_objects.parquet").set_index("structure_id")
    _require(set(frozen.index) == set(a.structure_id), "Transformed RNA membership differs")
    _require(np.all(shapely.equals(a.geometry.to_numpy(), shapely.from_wkb(frozen.loc[a.structure_id].geometry_wkb))),
             "Transformed RNA geometries differ from the frozen result")
    computed = {}; checks = []
    for weight in WEIGHTS:
        result = assign_multicell_geometry(a, b, position_weight=weight, max_dist_um=200, area_ratio=2.5)
        previous = pd.read_parquet(root / "results" / donor / f"assignment_position_weight_{weight}.parquet")
        keys = ["moving_id", "fixed_id"]
        # The original audit adds pair annotations to the primary assignment
        # DataFrame in place. Replay validates the algorithm's complete output;
        # anatomy/control annotations are explicitly outside this bounded run.
        _require(set(result.columns) <= set(previous.columns), "Frozen assignment lacks algorithm output columns")
        annotations = sorted(set(previous.columns) - set(result.columns))
        try:
            pd.testing.assert_frame_equal(result.sort_values(keys).reset_index(drop=True),
                                          previous[result.columns].sort_values(keys).reset_index(drop=True), check_exact=True)
        except AssertionError as error:
            raise SequentialError(f"Bounded assignment replay differs at position weight {weight}: {error}") from error
        computed[weight] = result
        checks.append({"position_weight": weight, "n_pairs": len(result), "exact_membership_and_values": True,
                       "assignment_columns_compared": list(result.columns),
                       "downstream_annotation_columns_not_recomputed": annotations})
    output.mkdir(parents=True, exist_ok=False)
    for weight, table in computed.items():
        table.to_parquet(output / f"assignment_position_weight_{weight}.parquet", index=False)
    receipt = {"schema": "sequential_assignment_replay_v1", "status": "passed", "bundle": str(root),
               "output": str(output), "donor_id": donor, "correspondence_policy": POLICY,
               "artifact_manifest_sha256": verification["artifact_manifest_sha256"],
               "workflow_source_file": verification["workflow_source_file"],
               "workflow_source_sha256": verification["workflow_source_sha256"],
               "native_section_bindings": [{k: s.get(k) for k in ("role", "donor_id", "source_serial", "roi", "acquisition_id", "output")} for s in (rna[0], protein[0])],
               "native_objects_reused": {role: len(table) for role, table in native.items()},
               "RNA_multicell_geometries_verified": len(a), "assignment_checks": checks,
               "native_construction_executed": False, "transforms_refit": False,
               "controls_recomputed": False, "scientific_review_decisions_changed": False,
               "individual_pairs_validated": False, "production_activated": False,
               "runtime": verification["runtime"],
               "completed_utc": datetime.now(timezone.utc).isoformat()}
    (output / "replay_receipt.json").write_text(json.dumps(receipt, indent=2) + "\n")
    return receipt
