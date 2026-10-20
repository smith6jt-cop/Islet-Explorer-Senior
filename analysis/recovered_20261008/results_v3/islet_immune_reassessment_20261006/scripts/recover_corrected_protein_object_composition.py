#!/usr/bin/env python3
"""Recover only object composition while preserving the failed execution ledger."""
import os
for _name in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMBA_NUM_THREADS"):
    os.environ[_name] = "6"
from pathlib import Path
import argparse
import importlib.util
import json
import shutil
import sys

import pandas as pd
import pyarrow.parquet as pq

HERE = Path(__file__).resolve().parent
spec = importlib.util.spec_from_file_location("protein_recovery_controller", HERE / "run_protein_own_frozen_reference_repair.py")
C = importlib.util.module_from_spec(spec); spec.loader.exec_module(C)
H, P = C.H, C.P
REPO = HERE.parents[3]
BUILDER = P.TYPING / "build_endocrine_object_composition.py"
INVENTORY = Path("/data/xenium/io60_cohort_20261003/input_inventory.json")
SUCCESS = ("full_row_release", "current_generation_assignments", "admitted_marker_display")

def read(path, columns):
    return pq.ParquetFile(path).read(columns=columns).to_pandas()

def audit_composition(root, composition, typing, assignments, donors):
    """Compare every object member with independent terminal assignment admission."""
    records = []
    for donor in donors:
        native = typing / f"cells/donor_id={donor}/data.parquet"
        target = composition / f"donor_id={donor}/data.parquet"
        assignment = assignments / f"donor_id={str(int(donor))}/data.parquet"
        cells = read(native, ["object_id", "endocrine_object_id", "endocrine_object_class", "cell_class", "refined_class", "is_endocrine"])
        cells = cells.loc[cells.endocrine_object_id.notna()].copy()
        if cells.object_id.isna().any() or cells.object_id.duplicated().any():
            raise ValueError(f"Invalid original composition member keys: {donor}")
        a = read(assignment, ["object_id", "analysis_eligible"])
        if a.object_id.duplicated().any(): raise ValueError(f"Repeated assignment keys: {donor}")
        cells["admitted"] = cells.object_id.map(a.set_index("object_id").analysis_eligible)
        if cells.admitted.isna().any(): raise ValueError(f"Missing current assignment for original object member: {donor}")
        cells["admitted"] = cells.admitted.astype(bool)
        objects = read(target, None)
        if objects.endocrine_object_id.duplicated().any() or set(objects.endocrine_object_id) != set(cells.endocrine_object_id):
            raise ValueError(f"Object universe changed: {donor}")
        groups = {str(k): g for k, g in cells.groupby("endocrine_object_id", sort=False)}
        n_multi = 0; n_multi_with_exclusion = 0; n_ordinary_uncertain = 0
        uncertainty = {"Endocrine (unresolved subtype)", "Endocrine (ambiguous subtype)"}
        for row in objects.to_dict("records"):
            raw = groups[str(row["endocrine_object_id"])]
            kept = raw.loc[raw.admitted]; omitted = raw.loc[~raw.admitted]
            classes = raw.endocrine_object_class.dropna().unique()
            if len(classes) != 1 or row["original_object_class"] != classes[0]:
                raise ValueError(f"Original object class changed: {donor}/{row['endocrine_object_id']}")
            expected = {"n_cells_in_object": len(raw), "n_original_cells_in_object": len(raw),
                "n_analysis_admitted_cells_in_object": len(kept), "n_analysis_excluded_cells_in_object": len(omitted),
                "endocrine_denominator": int(kept.is_endocrine.fillna(False).sum())}
            for field, value in expected.items():
                if int(row[field]) != value: raise ValueError(f"Exact member census changed: {donor} {field}")
            if bool(row["composition_analysis_eligible"]) != bool(len(kept)):
                raise ValueError(f"Composition admission changed: {donor}")
            excluded = set(str(row["excluded_cell_ids"]).split("|")) if row["excluded_cell_ids"] else set()
            if excluded != set(omitted.object_id): raise ValueError(f"Excluded object member keys changed: {donor}")
            for field, label in (("primary_count__", "cell_class"), ("refined_count__", "refined_class")):
                counts = {str(k).lower().replace(" ", "_"): int(n) for k, n in kept[label].astype(str).value_counts().items()}
                for column in objects.columns:
                    if column.startswith(field):
                        if int(row[column]) != counts.get(column.removeprefix(field), 0):
                            raise ValueError(f"Admitted class count changed: {donor} {column}")
            unresolved = int((kept.is_endocrine.fillna(False) & ~kept.cell_class.isin({"Alpha", "Beta", "Delta"})).sum())
            if int(row["count_endocrine__unresolved_subtype"]) != unresolved:
                raise ValueError(f"Ordinary endocrine uncertainty lost: {donor}")
            identity_unresolved = int(kept.cell_class.isin({"Ambiguous", "Unresolved", "Unresolved (no lineage evidence)"} | uncertainty).sum())
            if int(row["identity_unresolved_count"]) != identity_unresolved:
                raise ValueError(f"Admitted identity uncertainty count changed: {donor}")
            for subtype in ("alpha", "beta", "delta"):
                available = bool(row[f"available__endocrine__{subtype}_fraction"])
                value = row[f"feature__endocrine__{subtype}_fraction"]
                if (unresolved or not expected["endocrine_denominator"]) and (available or pd.notna(value)):
                    raise ValueError(f"Unsupported complete hormone composition manufactured: {donor}")
            n_ordinary_uncertain += int(kept.cell_class.isin(uncertainty).sum())
            if classes[0] == "EndocrineSingle" and len(raw) > 1:
                n_multi += 1; n_multi_with_exclusion += int(bool(len(omitted)))
        records.append(dict(donor_id=donor, n_objects=len(objects), n_exact_original_member_keys=len(cells),
            n_analysis_admitted_members=int(cells.admitted.sum()), n_analysis_excluded_members=int((~cells.admitted).sum()),
            n_multicell_single_annotations=n_multi, n_multicell_single_annotations_with_excluded_members=n_multi_with_exclusion,
            n_admitted_ordinary_endocrine_uncertainty=n_ordinary_uncertain, source=str(native), source_sha256=H.sha(native),
            assignment=str(assignment), assignment_sha256=H.sha(assignment), output=str(target), output_sha256=H.sha(target)))
        print("exact object membership/admission parity", donor, len(objects), len(cells), flush=True)
    receipt = dict(status="passed", donors=records, donor_count=len(records),
        n_objects=sum(r["n_objects"] for r in records), n_exact_original_member_keys=sum(r["n_exact_original_member_keys"] for r in records),
        identity_or_thresholds_changed=False, accuracy_validated=False,
        interpretation="Exact source object membership and terminal current-generation assignment contribution parity; no independent biological validation")
    H.save(root / "recovered_endocrine_object_composition_integrity.json", receipt)
    return receipt

def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--root", type=Path, default=C.DEFAULT_OUT)
    ap.add_argument("--review-receipt", type=Path, required=True)
    ap.add_argument("--preflight-only", action="store_true")
    args = ap.parse_args(); root = args.root.resolve()
    config = C.freeze(root); C.check_sources(config)
    review = json.loads(args.review_receipt.read_text())
    if review.get("status") != "passed" or review.get("approved_controller_sha256") != H.sha(__file__) or review.get("approved_composition_builder_sha256") != H.sha(BUILDER):
        raise ValueError("Independent review does not approve these exact recovery/controller sources")
    prior_path = root / "protein_release_followup_run.json"
    prior = json.loads(prior_path.read_text()); prior_sha = H.sha(prior_path)
    if prior["status"] != "failed" or [s["key"] for s in prior["stages"][:3]] != list(SUCCESS) or any(s["status"] != "complete" for s in prior["stages"][:3]) or prior["stages"][3]["key"] != "endocrine_object_composition" or prior["stages"][3]["status"] != "failed":
        raise ValueError("Recovery requires exactly three preserved successful stages and the failed object stage")
    typing = root / "own_compartment/typing"; run = json.loads((typing / "run.json").read_text())
    if run["status"] != "complete" or len(run["donors"]) != 25 or set(run["donors"]) != set(config["donors"]):
        raise ValueError("Corrected all25 source typing is not terminal")
    for stage in prior["stages"][:3]:
        for output in stage["outputs"]:
            if H.sha(output["path"]) != output["sha256"]: raise ValueError("Preserved successful output changed: " + output["path"])
    for source in prior["source_scripts"]:
        if H.sha(source["snapshot"]) != source["sha256"]: raise ValueError("Prior executed source snapshot changed")
    for item in prior["immutable_inputs"]:
        if H.sha(item["path"]) != item["sha256"]: raise ValueError("Prior immutable input changed")
    ledger = root / "protein_release_followup_recovery_run.json"
    out = root / "endocrine_object_composition_current_generation_recovery_20261008"
    inputs = [dict(path=str(p), sha256=H.sha(p)) for p in (prior_path, args.review_receipt.resolve(),
        root / "pancreas_cell_release_current_generation/manifest.json", root / "full_row_protein_release_integrity.json",
        root / "pancreas_analysis_assignments_current_generation/manifest.json", root / "native_marker_display_current_generation/native_display_sampling_manifest.json",
        root / "current_generation_endocrine_single_exclusions.parquet", root / "current_generation_endocrine_single_exclusions_manifest.json")]
    paths = [Path(__file__).resolve(), BUILDER, P.TYPING / "protein_analysis_admission.py", P.TYPING / "build_islet_geometry.py",
        P.TYPING / "context.py", P.TYPING / "islet_membership.py", P.TYPING / "subpanel_io.py"]
    sources = [dict(path=str(p), sha256=H.sha(p)) for p in paths]
    print("OBJECT-ONLY RECOVERY PREFLIGHT PASSED; THREE SUCCESSFUL STAGES PRESERVED", flush=True)
    if args.preflight_only: return
    if ledger.exists() or out.exists(): raise FileExistsError("Preserve any prior recovery evidence before retry")
    snapshot = root / "composition_recovery_source_snapshot"; snapshot.mkdir()
    for i, source in enumerate(sources):
        target = snapshot / f"{i:02d}_{Path(source['path']).name}"; shutil.copy2(source["path"], target); source["snapshot"] = str(target)
    record = dict(status="running", controller_pid=os.getpid(), started_at_UTC=C.now(), current_stage="endocrine_object_composition",
        recovery_of=dict(path=str(prior_path), sha256=prior_sha, failed_stage="endocrine_object_composition",
            preserved_successful_stage_keys=list(SUCCESS), preserved_partial_composition=prior["output_roots"]["composition"]),
        source_scripts=sources, immutable_inputs=inputs, donors=run["donors"],
        output_roots={**prior["output_roots"], "composition": str(out)}, stages=prior["stages"][:3],
        threads=6, RSS_limit_GiB=32, release_promoted=False, accuracy_validated=False)
    H.save(ledger, record)
    env = os.environ.copy(); env.update(IO60_GEOMETRY_ROOT=run["geometry_root"], IO60_TYPING_ROOT=str(typing),
        IO60_EXPRESSION_ROOT=run["expression_root"], IO60_ALIGNED_ROOT=run["aligned_root"], IO60_SCALE_PATH=run["scale"])
    cmd = [sys.executable, str(BUILDER), "--typing-root", str(typing), "--aligned-root", run["aligned_root"], "--inventory", str(INVENTORY),
        "--out", str(out), "--analysis-status", "current_generation_source_corrected_biological_review_pending",
        "--cell-exclusions", str(root / "current_generation_endocrine_single_exclusions.parquet")]
    stage = dict(key="endocrine_object_composition", started_at_UTC=C.now())
    try:
        stage.update(H.run(cmd, env, root / "endocrine_object_composition_recovery_20261008.log", 32))
        if stage["returncode"]: raise RuntimeError("Recovered object stage exited " + str(stage["returncode"]))
        for item in [*sources, *inputs]:
            if H.sha(item["path"]) != item["sha256"]: raise ValueError("Frozen recovery source/input changed")
        C.check_sources(config)
        manifest = json.loads((out / "manifest.json").read_text())
        if manifest["typing_run_sha256"] != H.sha(typing / "run.json") or manifest["composition_builder_sha256"] != H.sha(BUILDER) or manifest["cell_exclusion_overlay_sha256"] != H.sha(root / "current_generation_endocrine_single_exclusions.parquet"):
            raise ValueError("Recovered composition manifest source binding changed")
        if {q["donor_id"] for q in manifest["per_donor"]} != set(run["donors"]) or len(manifest["per_donor"]) != 25:
            raise ValueError("Recovered composition donor scope changed")
        stage["object_integrity"] = audit_composition(root, out, typing, Path(prior["output_roots"]["assignments"]), run["donors"])
        stage.update(status="complete", finished_at_UTC=C.now(), outputs=H.artifacts([out]),
            manifest=str(out / "manifest.json"), manifest_sha256=H.sha(out / "manifest.json"))
    except Exception as exc:
        stage.update(status="failed", error=str(exc), finished_at_UTC=C.now()); record["stages"].append(stage)
        record.update(status="failed", finished_at_UTC=C.now()); H.save(ledger, record); raise
    record["stages"].append(stage); record.update(status="complete", current_stage=None, finished_at_UTC=C.now()); H.save(ledger, record)
    print("OBJECT COMPOSITION RECOVERY COMPLETE; FAILED HISTORY PRESERVED", flush=True)

if __name__ == "__main__":
    main()
