#!/usr/bin/env python3
"""Assemble and audit fresh full-row protein views, display and object composition."""
import os
for _name in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMBA_NUM_THREADS"):
    os.environ[_name] = "6"
from pathlib import Path
import argparse
import importlib.util
import json
import shutil
import sys

import numpy as np
import pandas as pd
import pyarrow.parquet as pq

HERE = Path(__file__).resolve().parent
spec = importlib.util.spec_from_file_location("protein_all25_controller", HERE / "run_protein_own_frozen_reference_repair.py")
C = importlib.util.module_from_spec(spec); spec.loader.exec_module(C)
H, P = C.H, C.P
REPO = HERE.parents[3]
BUILDER = REPO / "core/phenocycler/scripts/typing/build_pancreas_cell_release.py"
MATERIALIZER = REPO / "analysis/xenium/scripts/materialize_protein_analysis_assignments.py"
DISPLAY = REPO / "analysis/xenium/scripts/prepare_native_typing_display.py"
COMPOSITION = REPO / "core/phenocycler/scripts/typing/build_endocrine_object_composition.py"
INPUT_MANIFEST = REPO / "results/v3/release_20261003/current_input_manifest.json"
INVENTORY = Path("/data/xenium/io60_cohort_20261003/input_inventory.json")
SINGLE_REASONS = {"endocrine_single_no_usable_hormone", "endocrine_single_tied_highest_hormones"}
RAW_META = {"donor_id", "object_id", "image", "parent_raw", "qupath_classification", "X_centroid", "Y_centroid",
    "nucleus_area_um2", "cell_area", "cell_solidity", "nucleus_cell_area_ratio", "signed_dist_islet_um", "islet_num",
    "cell_region", "qc_analysis_eligible", "qc_estimation_eligible", "analysis_reason", "estimation_reason", "in_alignment_block"}

def read(path, columns=None):
    return pq.ParquetFile(path).read(columns=columns).to_pandas()

def exact_sample(path, columns, keys):
    pieces = []
    for batch in pq.ParquetFile(path).iter_batches(columns=columns, batch_size=100000):
        frame = batch.to_pandas()
        keep = frame.object_id.isin(keys)
        if keep.any(): pieces.append(frame.loc[keep])
    result = pd.concat(pieces, ignore_index=True) if pieces else pd.DataFrame(columns=columns)
    if result.object_id.duplicated().any() or set(result.object_id) != set(keys):
        raise ValueError(f"Exact source sample did not match once: {path}")
    return result.set_index("object_id").loc[keys]

def audit_release(root, release, run):
    rows = []
    for donor in run["donors"]:
        source = Path(run["geometry_root"]) / f"revised_cells/donor_id={donor}/data.parquet"
        target = release / f"donor_id={donor}/data.parquet"
        original = read(source, ["object_id", "qc_analysis_eligible"])
        qc_fields = ["qc_analysis_eligible", "qc_estimation_eligible", "analysis_reason", "estimation_reason"]
        released = read(target, ["object_id", "qc_analysis_eligible", "analysis_reason", "nucleus_area_um2", "cell_class", "excluded_reason",
            *["original_geometry_" + c for c in qc_fields]])
        for label, frame in (("original", original), ("release", released)):
            if frame.object_id.isna().any() or frame.object_id.duplicated().any(): raise ValueError(f"Invalid {label} UUIDs: {donor}")
        if len(original) != len(released) or set(original.object_id) != set(released.object_id):
            raise ValueError(f"Full revised-geometry universe changed: {donor}")
        original_qc = original.set_index("object_id").qc_analysis_eligible.fillna(False).astype(bool)
        qc = released.object_id.map(original_qc)
        if not qc.eq(released.qc_analysis_eligible.fillna(False).astype(bool)).all(): raise ValueError(f"Release QC mask changed: {donor}")
        if (qc & (~np.isfinite(released.nucleus_area_um2) | released.nucleus_area_um2.le(0))).any():
            raise ValueError(f"No-nucleus object became calculation-eligible: {donor}")
        if not qc.eq(released.cell_class.notna()).all(): raise ValueError(f"Primary availability differs from new geometry QC: {donor}")
        authoritative_path = Path(run["original_geometry_root"]) / f"revised_cells/donor_id={donor}/data.parquet"
        authoritative = read(authoritative_path, ["object_id", *qc_fields])
        if authoritative.object_id.isna().any() or authoritative.object_id.duplicated().any() or set(authoritative.object_id) != set(released.object_id):
            raise ValueError(f"Authoritative original QC provenance universe changed: {donor}")
        keyed = authoritative.set_index("object_id")
        for field in qc_fields:
            a = released.object_id.map(keyed[field]); b = released["original_geometry_" + field]
            if not np.array_equal(a.isna(), b.isna()) or not a.astype("string").fillna("").eq(b.astype("string").fillna("")).all():
                raise ValueError(f"Original geometry QC provenance field changed: {donor} {field}")
        keys = sorted(original.object_id)[:2048]
        source_columns = pq.ParquetFile(source).schema.names
        mapping = {c: ("raw_source_parent_islet_num" if c == "islet_num" else c if c in RAW_META else "raw_" + c)
                   for c in source_columns if c not in {"object_id", "donor_id"}}
        a = exact_sample(source, ["object_id", *mapping], keys)
        b = exact_sample(target, ["object_id", *mapping.values()], keys)
        for field, output in mapping.items():
            left, right = a[field], b[output]
            if pd.api.types.is_numeric_dtype(left) or pd.api.types.is_bool_dtype(left):
                equal = np.array_equal(pd.to_numeric(left).to_numpy(float), pd.to_numeric(right).to_numpy(float), equal_nan=True)
            else:
                equal = left.fillna("").astype(str).equals(right.fillna("").astype(str).rename(field))
            if not equal: raise ValueError(f"Original observation changed in release: {donor} {field}")
        rows.append(dict(donor_id=donor, source=str(source), source_sha256=H.sha(source), release=str(target),
            release_sha256=H.sha(target), n_full_UUIDs=len(original), n_QC_eligible=int(qc.sum()),
            n_primary_excluded=int(released.cell_class.eq("Excluded").sum()), n_exact_observation_sample=len(keys),
            n_source_columns_verified=len(mapping), n_excluded_single=int(released.excluded_reason.isin(SINGLE_REASONS).sum())))
        rows[-1].update(original_geometry_QC_source=str(authoritative_path), original_geometry_QC_sha256=H.sha(authoritative_path),
            n_original_QC_exact_keys=len(authoritative), original_geometry_QC_fields_preserved=qc_fields,
            n_original_QC_eligible=int(authoritative.qc_analysis_eligible.fillna(False).sum()),
            n_already_ineligible_current_reason_changed=int((~released.original_geometry_qc_analysis_eligible.fillna(False).astype(bool) &
                released.original_geometry_analysis_reason.fillna("").ne(released.analysis_reason.fillna(""))).sum()))
        print("full-row release universe/raw observation parity", donor, len(original), flush=True)
    receipt = dict(status="passed", donor_count=len(rows), n_full_UUIDs=sum(r["n_full_UUIDs"] for r in rows),
        n_exact_observation_sample=sum(r["n_exact_observation_sample"] for r in rows), donors=rows,
        all_raw_observations_retained=True, exact_new_geometry_QC_mask_preserved=True, no_nucleus_admission=False,
        original_geometry_QC_flags_and_reasons_preserved_on_all_keys=True,
        accuracy_validated=False, interpretation="Full UUID/QC universe and exact raw-source sampled parity; no biological validation")
    H.save(root / "full_row_protein_release_integrity.json", receipt)
    return receipt

def single_overlay(root, assignments, typing, donors):
    records = []; parts = []
    for donor in donors:
        assignment = assignments / f"donor_id={str(int(donor))}/data.parquet"
        fields = ["donor_id", "object_id", "endocrine_object_id", "endocrine_object_class", "qc_analysis_eligible", "cell_class", "analysis_eligible", "excluded_reason", "analysis_exclusion_reason"]
        frame = read(assignment, fields)
        selected = frame.loc[frame.qc_analysis_eligible.fillna(False) & frame.endocrine_object_class.eq("EndocrineSingle") & ~frame.analysis_eligible.fillna(False)].copy()
        if not selected.cell_class.eq("Excluded").all() or not selected.excluded_reason.isin(SINGLE_REASONS).all():
            raise ValueError(f"Non-native or unsupported exclusion entered the fresh Single overlay: {donor}")
        calls = read(typing / f"cells/donor_id={donor}/data.parquet", ["object_id", "endocrine_object_class", "excluded_reason"])
        expected = calls.loc[calls.endocrine_object_class.eq("EndocrineSingle") & calls.excluded_reason.isin(SINGLE_REASONS)]
        if set(selected.object_id) != set(expected.object_id): raise ValueError(f"Assignment/native Single exclusion keys differ: {donor}")
        selected = selected[["donor_id", "object_id", "endocrine_object_id", "excluded_reason"]].rename(columns={"excluded_reason": "exclusion_reason"})
        if selected.endocrine_object_id.isna().any() or selected.object_id.duplicated().any(): raise ValueError(f"Invalid current excluded Single membership: {donor}")
        parts.append(selected); records.append(dict(donor_id=donor, n_exact_single_keys=len(selected), assignment_sha256=H.sha(assignment)))
    result = pd.concat(parts, ignore_index=True)
    if result.duplicated(["donor_id", "object_id"]).any() or result.duplicated(["donor_id", "endocrine_object_id"]).any():
        raise ValueError("Fresh Single overlay repeats a cell/object key")
    path = root / "current_generation_endocrine_single_exclusions.parquet"
    result.to_parquet(path, index=False)
    manifest = json.loads((assignments / "manifest.json").read_text())
    receipt = dict(status="passed", path=str(path), sha256=H.sha(path), n_exact_single_keys=len(result), donors=records,
        source_assignments_manifest=str(assignments / "manifest.json"), source_assignments_manifest_sha256=H.sha(assignments / "manifest.json"),
        mandatory_registry_census=manifest["census"], scope="New terminal typed/QC-eligible unsupported Singles only; mandatory keys already absent upstream remain documented independently",
        no_all_QC_failures_reclassified_as_Singles=True)
    H.save(root / "current_generation_endocrine_single_exclusions_manifest.json", receipt)
    return path

def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--root", type=Path, default=C.DEFAULT_OUT)
    ap.add_argument("--preflight-only", action="store_true")
    args = ap.parse_args(); root = args.root.resolve()
    config = C.freeze(root); C.check_sources(config)
    typing = root / "own_compartment/typing"; typing_run = typing / "run.json"
    run = json.loads(typing_run.read_text())
    if run["status"] != "complete" or set(run["donors"]) != set(config["donors"]) or len(run["donors"]) != 25:
        raise ValueError("Fresh corrected all25 typing is not terminal")
    if run["generation_config_sha256"] != H.sha(root / "generation_config.json"): raise ValueError("Typing generation binding differs")
    release = root / "pancreas_cell_release_current_generation"
    assignments = root / "pancreas_analysis_assignments_current_generation"
    display = root / "native_marker_display_current_generation"
    composition = root / "endocrine_object_composition_current_generation"
    ledger = root / "protein_release_followup_run.json"
    paths = [Path(__file__), HERE / "run_protein_own_frozen_reference_repair.py", HERE / "resume_protein_compartment_candidate.py",
        BUILDER, MATERIALIZER, DISPLAY, COMPOSITION, P.TYPING / "protein_analysis_admission.py", P.TYPING / "build_islet_geometry.py",
        P.TYPING / "context.py", P.TYPING / "islet_membership.py", P.TYPING / "subpanel_io.py"]
    sources = [dict(path=str(p), sha256=H.sha(p)) for p in paths]
    if not run.get("original_geometry_root"): raise ValueError("Corrected release requires explicit authoritative original QC provenance")
    inputs = [dict(path=str(p), sha256=H.sha(p)) for p in (typing_run, INPUT_MANIFEST, INVENTORY, root / "generation_config.json", P.ASSIGN / "manifest.json",
        *[Path(run["original_geometry_root"]) / f"revised_cells/donor_id={d}/data.parquet" for d in run["donors"]])]
    print("FRESH ALL25 RELEASE/FOLLOW-UP PREFLIGHT PASSED", flush=True)
    if args.preflight_only: return
    if ledger.exists() or any(p.exists() for p in (release, assignments, display, composition)):
        raise FileExistsError("Preserve prior partial outputs and execution receipts before any retry")
    snapshot = root / "release_followup_source_snapshot"; snapshot.mkdir()
    for i, source in enumerate(sources):
        target = snapshot / f"{i:02d}_{Path(source['path']).name}"; shutil.copy2(source["path"], target); source["snapshot"] = str(target)
    record = dict(status="running", controller_pid=os.getpid(), started_at_UTC=C.now(), source_scripts=sources, immutable_inputs=inputs,
        output_roots=dict(release=str(release), assignments=str(assignments), display=str(display), composition=str(composition)),
        donors=run["donors"], stages=[], threads=6, RSS_limit_GiB=180, release_promoted=False, accuracy_validated=False)
    H.save(ledger, record)
    env = os.environ.copy(); env.update(IO60_GEOMETRY_ROOT=run["geometry_root"], IO60_TYPING_ROOT=str(typing),
        IO60_EXPRESSION_ROOT=run["expression_root"], IO60_ALIGNED_ROOT=run["aligned_root"], IO60_SCALE_PATH=run["scale"])
    steps = [
        ("full_row_release", [sys.executable, str(BUILDER), "--arm-root", str(typing), "--out", str(release), "--input-manifest", str(INPUT_MANIFEST)], [release]),
        ("current_generation_assignments", [sys.executable, str(MATERIALIZER), "--source", str(release), "--out", str(assignments), "--policy", "current_generation", "--mandatory-exclusions", str(P.ASSIGN)], [assignments]),
        ("admitted_marker_display", [sys.executable, str(DISPLAY), "--source", str(release), "--out", str(display), "--exclusions", str(assignments)], [display]),
        ("endocrine_object_composition", [sys.executable, str(COMPOSITION), "--typing-root", str(typing), "--aligned-root", run["aligned_root"], "--inventory", str(INVENTORY), "--out", str(composition), "--analysis-status", "current_generation_source_corrected_biological_review_pending", "--cell-exclusions", str(root / "current_generation_endocrine_single_exclusions.parquet")], [composition])]
    for name, command, outputs in steps:
        C.check_sources(config)
        for item in [*sources, *inputs]:
            if H.sha(item["path"]) != item["sha256"]: raise ValueError(f"Frozen follow-up source/input changed: {item['path']}")
        stage = dict(key=name, started_at_UTC=C.now()); record["current_stage"] = name; H.save(ledger, record)
        print("START", name, flush=True)
        try:
            stage.update(H.run(command, env, root / (name + "_current_generation.log"), 180))
            if stage["returncode"]: raise RuntimeError(f"{name} exited {stage['returncode']}")
            if name == "full_row_release":
                stage["full_row_integrity"] = audit_release(root, release, run)
            elif name == "current_generation_assignments":
                summary = json.loads((assignments / "manifest.json").read_text())
                if summary["policy"] != "current_generation" or len(summary["donors"]) != 25:
                    raise ValueError("Fresh assignment view was not materialized with corrected all25 current policy")
                approved = sum(q["typing_excluded_cells"] for q in json.loads((P.ASSIGN / "manifest.json").read_text())["donors"])
                if summary["census"]["mandatory_exclusion_keys"] != approved: raise ValueError("Mandatory exact-key registry size changed")
                stage["computed_census"] = summary["census"]
                single_overlay(root, assignments, typing, run["donors"])
            stage.update(status="complete", finished_at_UTC=C.now(), outputs=H.artifacts(outputs))
        except Exception as exc:
            stage.update(status="failed", error=str(exc), finished_at_UTC=C.now()); record["stages"].append(stage)
            record.update(status="failed", finished_at_UTC=C.now()); H.save(ledger, record); raise
        record["stages"].append(stage); H.save(ledger, record)
        print("COMPLETE", name, "peak GiB", stage["peak_process_group_RSS_GiB"], flush=True)
    record.update(status="complete", current_stage=None, finished_at_UTC=C.now()); H.save(ledger, record)
    print("FULL-ROW RELEASE, FRESH ASSIGNMENTS, ADMITTED DISPLAY AND OBJECT COMPOSITION COMPLETE", flush=True)

if __name__ == "__main__":
    main()
