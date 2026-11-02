#!/usr/bin/env python3
"""Reconcile a complete PROT-1 candidate against bound inputs and unchanged primary calls."""
import argparse
import ast
import hashlib
import json
import subprocess
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pyarrow.parquet as pq


def sha(path):
    with path.open("rb") as handle:
        return hashlib.file_digest(handle, "sha256").hexdigest()


def check_optional_declaration(source):
    matches = [node for node in ast.walk(ast.parse(source))
               if isinstance(node, ast.Call) and isinstance(node.func, ast.Name)
               and node.func.id == "Input" and len(node.args) >= 2
               and ast.unparse(node.args[0]) == "f'pln_blocks:{donor}'"]
    if (len(matches) != 1
            or ast.unparse(matches[0].args[1]) != "upstream / f'05_normalization/alignment_block/donor_id={donor}.parquet'"
            or not any(keyword.arg == "required" and isinstance(keyword.value, ast.Constant)
                       and keyword.value.value is False for keyword in matches[0].keywords)):
        raise ValueError("producer does not explicitly declare the known optional pLN input contract")


def verify_optional_declaration(snapshot_path):
    producer = Path(__file__).resolve().parents[2] / "core/multiplex5000/multiplex5000/stages/io60_mixed_protein_pilot.py"
    snapshot = json.loads(snapshot_path.read_text())
    key = str(producer.relative_to(Path(__file__).resolve().parents[2]))
    matches = [item for item in snapshot["sources"] if item["path"] == key]
    if len(matches) != 1 or matches[0]["sha256"] != sha(producer):
        raise ValueError("optional input declaration changed since the producer snapshot")
    check_optional_declaration(producer.read_text())
    return {"path": str(producer), "sha256": sha(producer),
            "snapshot_path": str(snapshot_path), "snapshot_sha256": sha(snapshot_path),
            "input_name_template": "pln_blocks:{donor}", "required": False,
            "path_template": "{upstream_release}/05_normalization/alignment_block/donor_id={donor}.parquet"}


def verify_launch_bindings(path, receipt, snapshot_path):
    launch = json.loads(path.read_text())
    before, consumed = launch["inputs"], receipt["inputs"]
    if (len({row["name"] for row in before}) != len(before)
            or len({row["name"] for row in consumed}) != len(consumed)
            or {row["name"]: row for row in before} != {row["name"]: row for row in consumed}):
        raise ValueError("stage receipt inputs differ from the full pre-launch binding roster")
    if launch["settings"] != receipt["settings"]:
        raise ValueError("stage settings differ from the pre-launch contract")
    if (Path(launch["project"]["path"]) != Path(receipt["project"])
            or launch["project"]["sha256"] != sha(Path(receipt["project"]))
            or Path(launch["producer_snapshot"]["path"]) != snapshot_path
            or launch["producer_snapshot"]["sha256"] != sha(snapshot_path)):
        raise ValueError("launch project or producer snapshot binding changed")
    return {"path": str(path), "sha256": sha(path), "inputs": len(before)}


def verify_audit_results(root, expected_donors, receipt):
    """Require explicit evidence; the generic validator's exit code is insufficient."""
    summary = json.loads((root / "summary.json").read_text())
    rows = json.loads((root / "per_acquisition.json").read_text())
    if summary["hash_failures"]:
        raise ValueError("receipt hash verification failed")
    if summary["identity_failures"]:
        raise ValueError("mixed acquisition identity verification failed")
    if any(row.get("identity_check") != "pass" for row in rows):
        raise ValueError("every acquisition needs an explicit passing identity check")
    if (len(rows) != len(expected_donors)
            or {row.get("donor_id") for row in rows} != set(expected_donors)
            or summary["acquisitions"] != len(rows)):
        raise ValueError("acquisition identity roster differs from the declared cohort")
    by_donor = {row["donor_id"]: row for row in rows}
    optional_absences = []
    for check in summary["unavailable_hashes"]:
        bindings = [item for item in receipt["inputs"] if item["path"] == check["path"]]
        item = bindings[0] if len(bindings) == 1 else {}
        donor = item.get("name", "").removeprefix("pln_blocks:")
        expected_path = Path(receipt["settings"]["upstream_release"]) / f"05_normalization/alignment_block/donor_id={donor}.parquet"
        # Stage 06 explicitly declares this one input family required=False.
        # Its absence is immaterial only for sections with no pLN cells at all.
        if (check.get("role") != "input" or check.get("expected") is not None
                or check.get("actual") is not None or item.get("exists") is not False
                or item.get("sha256") is not None or donor not in by_donor
                or item.get("name") != f"pln_blocks:{donor}"
                or sum(row.get("name") == item.get("name") for row in receipt["inputs"]) != 1
                or Path(check["path"]) != expected_path
                or by_donor[donor].get("embedded_pln") != 0
                or Path(check["path"]).exists() or Path(check["path"]).is_symlink()):
            raise ValueError(f"required or unclassified unavailable receipt binding: {check['path']}")
        optional_absences.append({"input_name": item["name"], "path": item["path"], "exists": False,
                                  "donor_id": donor, "embedded_pln_cells": 0,
                                  "status": "optional_absent_no_pln_cells",
                                  "receipt_input": item,
                                  "identity_audit": "release_audit/per_acquisition.json",
                                  "identity_audit_sha256": sha(root / "per_acquisition.json"),
                                  "evidence": "exact receipt input remains absent; full aligned-table identity audit has zero pLN cells; stage inputs declare pln_blocks required=False"})
    if len({row["path"] for row in optional_absences}) != len(optional_absences):
        raise ValueError("duplicate optional absence in receipt audit")
    return summary, rows, optional_absences


def main():
    from validate_typing_release import spatial

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--release", required=True, type=Path)
    parser.add_argument("--out", required=True, type=Path)
    parser.add_argument("--baseline", required=True, type=Path)
    parser.add_argument("--producer-snapshot", type=Path,
                        help="exact source snapshot for this attempt; default is the run's original snapshot")
    parser.add_argument("--input-bindings", required=True, type=Path,
                        help="full input bindings saved immediately before this attempt")
    args = parser.parse_args()
    stage = args.release / "06_typing"
    assert (stage / "receipt.json").is_file(), "a failed or partial stage is not a candidate"
    args.out.mkdir(parents=True, exist_ok=True)
    for folder in ("samples", "overlays"):
        (args.out / "refinement_states" / folder).mkdir(parents=True, exist_ok=True)
    baseline = json.loads(args.baseline.read_text())
    receipt = json.loads((stage / "receipt.json").read_text())
    snapshot_path = (args.producer_snapshot or args.release.parents[1] / "producer_snapshot.json").resolve()
    launch_binding = verify_launch_bindings(args.input_bindings.resolve(), receipt, snapshot_path)
    optional_declaration = verify_optional_declaration(snapshot_path)
    subprocess.run([sys.executable, str(Path(__file__).with_name("validate_typing_release.py")),
                    "--root", str(args.release), "--out", str(args.out / "release_audit"), "--kind", "mixed"], check=True)
    audit_summary, audit_rows, optional_absences = verify_audit_results(
        args.out / "release_audit", receipt["settings"]["donors"], receipt)
    logs = {r["donor_id"]: r for r in json.loads((stage / "pancreas/refinement/refine_log.json").read_text())}
    assert set(logs) == set(receipt["settings"]["donors"]) == {r["donor_id"] for r in baseline["rows"]}
    assert len(audit_rows) == len(logs) and {row["donor_id"] for row in audit_rows} == set(logs)
    records = []
    for before in baseline["rows"]:
        donor = before["donor_id"]
        original = Path(before["primary_output"]["path"])
        assert sha(original) == before["primary_output"]["sha256"], f"{donor}: original primary changed"
        old = pq.ParquetFile(original).read(columns=["object_id", "cell_class"]).to_pandas()
        path = stage / f"pancreas/cells/donor_id={donor}/data.parquet"
        cells = pq.ParquetFile(path).read(columns=["object_id", "cell_class", "refined_class", "refine_state",
                                                 "refine_pass", "refine_note", "demoted", "qc_analysis_eligible"]).to_pandas()
        joined = cells.merge(old, on="object_id", how="outer", validate="one_to_one", indicator=True, suffixes=("", "_original"))
        assert joined._merge.eq("both").all() and joined.cell_class.equals(joined.cell_class_original)
        assert joined.qc_analysis_eligible.all()
        unavailable = joined.refine_state.str.startswith("unavailable_")
        assert joined.loc[unavailable, "refined_class"].equals(joined.loc[unavailable, "cell_class"])
        assert not joined.loc[unavailable, "demoted"].any()
        assert joined.loc[unavailable, "refine_note"].str.len().gt(0).all()
        assert joined.refine_state.notna().all()
        log = logs[donor]
        assert joined.refine_state.value_counts().to_dict() == log["refinement_states"]
        decisions = log["marker_availability"]["decisions"]
        assert log["marker_availability"]["contract"] == "stage03_and_stage05"
        assert len({row["marker"] for row in decisions}) == len(decisions)
        carried = {row["marker"] for row in decisions if row["carried"]}
        requested = {row["marker"] for row in decisions}
        assert carried == set(log["markers_used"])
        for decision in decisions:
            assert decision["carried"] == (decision["column_present"] and decision["panel_available"]
                                            and decision["normalization_status"] in {"aligned", "single_sample"})
        for source in log["marker_availability"]["sources"]:
            assert sha(Path(source["path"])) == source["sha256"]
        incomplete = {cls: sorted((set(markers) & requested) - carried)
                      for cls, markers in log["own_markers"].items()
                      if (set(markers) & requested) - carried}
        for cls in incomplete:
            assert not joined.loc[joined.cell_class.eq(cls), "demoted"].any()
            assert not joined.loc[joined.cell_class.eq(cls), "refine_state"].eq("evaluated").any()
            assert not ((joined.refined_class == cls) & (joined.cell_class != cls)).any()
            assert not any(group["cell_class"] == cls for group in log["groups"])
        side = pq.ParquetFile(stage / f"pancreas/refinement/donor_id={donor}/data.parquet").read(columns=["object_id", "refine_state", "concordance"]).to_pandas()
        assert len(side) == len(joined) and side.object_id.is_unique
        graph_unavailable = side.refine_state.isin(["unavailable_neighbours", "unavailable_markers", "unavailable_measurement", "unavailable_fragment_measurement", "not_applicable"])
        assert np.isfinite(side.loc[~graph_unavailable, "concordance"]).all()
        mixed = pq.ParquetFile(stage / f"mixed/donor_id={donor}/data.parquet").read(columns=["object_id", "cell_uid", "typing_state", "refine_state", "refine_note", "cell_class", "refined_class", "tissue"]).to_pandas()
        retained = mixed.loc[mixed.object_id.isin(cells.object_id)].merge(
            cells[["object_id", "refine_state", "refine_note"]], on="object_id", validate="one_to_one", suffixes=("", "_side"))
        assert retained.refine_state.equals(retained.refine_state_side) and retained.refine_note.equals(retained.refine_note_side)
        assert mixed.loc[mixed.tissue.eq("pancreatic_lymph_node"), ["cell_class", "refined_class", "refine_state"]].isna().all().all()
        sample = pd.read_csv(args.out / f"release_audit/samples/{donor}.csv.gz")
        sample = sample.merge(mixed[["cell_uid", "refine_state"]], on="cell_uid", validate="one_to_one")
        sample["review_state"] = sample.refine_state.fillna(sample.typing_state)
        spatial(args.out / "refinement_states", donor, sample, "review_state")
        records.append({"donor_id": donor, "primary_cells": len(cells), "identity_check": "pass",
                        "primary_calls_unchanged": True, "refinement_states": log["refinement_states"],
                        "measured_cells": log["n_measured"], "demoted_cells": log["n_demoted"],
                        "classes_with_unavailable_defining_markers": incomplete,
                        "carried_markers": log["markers_used"], "marker_availability": log["marker_availability"],
                        "primary_output": {"path": str(path), "sha256": sha(path)}})
        print(donor, "pass", log["refinement_states"], flush=True)
    provenance = json.loads((stage / "pancreas/refinement/representation_provenance.json").read_text())
    assert provenance["fit_donors"] == receipt["settings"]["representation_donors"]
    assert Path(provenance["parameter_path"]).is_file() and Path(provenance["scale_path"]).is_file()
    assert sha(Path(provenance["parameter_path"])) == provenance["parameter_sha256"]
    assert sha(Path(provenance["scale_path"])) == provenance["scale_sha256"]
    result = {"status": "pass", "scope": "computational candidate; no biological approval or release activation",
              "acquisitions": len(records), "primary_cells": sum(r["primary_cells"] for r in records),
              "receipt_hash_failures": audit_summary["hash_failures"],
              "receipt_unavailable_hashes": audit_summary["unavailable_hashes"],
              "receipt_optional_absences": optional_absences,
              "optional_input_declaration": optional_declaration,
              "input_launch_bindings": launch_binding,
              "receipt_required_unavailable_hashes": [],
              "mixed_identity_failures": audit_summary["identity_failures"],
              "verifier": {"path": str(Path(__file__).resolve()), "sha256": sha(Path(__file__))},
              "generic_validator": {"path": str(Path(__file__).with_name("validate_typing_release.py").resolve()),
                                    "sha256": sha(Path(__file__).with_name("validate_typing_release.py"))},
              "input_baseline": {"path": str(args.baseline), "sha256": sha(args.baseline)},
              "receipt": {"path": str(stage / "receipt.json"), "sha256": sha(stage / "receipt.json")}, "rows": records}
    (args.out / "remediation_checks.json").write_text(json.dumps(result, indent=2, allow_nan=False) + "\n")


if __name__ == "__main__":
    main()
