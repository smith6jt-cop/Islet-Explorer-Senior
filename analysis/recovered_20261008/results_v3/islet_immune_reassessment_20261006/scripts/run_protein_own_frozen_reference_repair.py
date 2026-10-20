#!/usr/bin/env python3
"""Isolated all-current-donor own-source repair with original reference IDs.

Preserves the historical generation. New own-source multipliers, floors and spans
are fitted on the declared same-block population. Biological accuracy is not
established by completing these processing and integrity stages.
"""
import os
for name in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMBA_NUM_THREADS"):
    os.environ[name] = "6"
from pathlib import Path
import argparse
import importlib.util
import json
import shutil
import sys
from datetime import datetime, timezone
import numpy as np
import pandas as pd
import pyarrow.parquet as pq

HERE = Path(__file__).resolve().parent
spec = importlib.util.spec_from_file_location("repair_helpers", HERE / "resume_protein_compartment_candidate.py")
H = importlib.util.module_from_spec(spec); spec.loader.exec_module(H)
P = H.P
METHOD = "own_compartment"
OLD = P.R / "cell_typing_validation/candidate_compartment_reprocessing"
DEFAULT_OUT = Path("/data/io60-calibration/own_compartment_frozen_reference_repair_20261008")
EXPECTED = {
    "subpanel.py": "f6e49e749600627d9af6abf6371621b2368ab7ed0d975b481a0a76dcfd3ccde9",
    "refine.py": "0ca08447ca5af60e538ce04ac9bc0354e012a99aa81f23ad7a4f0bb10000d7cf",
}

def now():
    return datetime.now(timezone.utc).isoformat()

def current_donors():
    sys.path.insert(0, str(P.ROOT / "core/phenocycler"))
    from phenocycler.cohort import DONOR_EXCLUSIONS
    manifest = json.loads((P.DATA / "qupath_manifest.json").read_text())
    denied = {P.canonical(d) for d in DONOR_EXCLUSIONS}
    donors = [str(i["donor_id"]) for i in manifest["images"] if P.canonical(i["donor_id"]) not in denied]
    if len(donors) != 25 or len(set(donors)) != 25 or set(donors) != set(manifest["expected_donors"]):
        raise ValueError("Current manifest is not the reviewed exact 25-donor cohort")
    return donors

def source_files():
    return [Path(__file__), HERE / "reprocess_protein_own_compartment_candidate.py",
            HERE / "resume_protein_compartment_candidate.py",
            HERE / "verify_protein_candidate_expression_sources.py",
            HERE / "verify_protein_candidate_primary_integrity.py",
            HERE / "audit_protein_repaired_primary_evidence.py",
            HERE / "audit_protein_candidate_missing_negatives.py",
            *[P.TYPING / n for n in ("subpanel.py", "build_cells.py", "refine.py", "immune_subtype.py", "subpanel_io.py", "islet_membership.py")],
            *[P.CELLTYPE / n for n in ("panel.py", "run_step2.py", "run_step3.py", "alignment_block.py", "score.py", "rules.py", "usability.py", "uniform_adapter.py")],
            *[P.ROOT / "core/phenocycler/phenocycler" / n for n in ("expression.py", "marker_registry.py", "cohort.py")]]

def freeze(out):
    donors = current_donors()
    for name, expected in EXPECTED.items():
        if H.sha(P.TYPING / name) != expected:
            raise ValueError(f"Corrected producer is not the reviewed generation: {name}")
    config = out / "generation_config.json"
    if config.exists():
        result = json.loads(config.read_text())
        if result["donors"] != donors:
            raise ValueError("Generation cohort changed")
        check_sources(result)
        return result
    out.mkdir(parents=True, exist_ok=True)
    reference = P.ROOT / "data/development/celltype/shifts.json"
    prior = json.loads(reference.read_text())
    if set(prior["donors"]) != set(donors):
        raise ValueError("Original reference source does not cover current exact cohort")
    snapshot = out / "producer_source_snapshot"
    snapshot.mkdir()
    scripts = []
    for i, path in enumerate(source_files()):
        archived = snapshot / f"{i:02d}_{path.name}"
        shutil.copy2(path, archived)
        scripts.append(dict(path=str(path), sha256=H.sha(path), snapshot=str(archived)))
    refcopy = out / "reference_donor_identity_source.json"
    shutil.copy2(reference, refcopy)
    inputs = [reference, P.DATA / "qupath_manifest.json", P.REGISTRY, P.POLICY,
              P.BLOCK / "manifest.json", P.GEOMETRY / "subpanel/subpanel_log.json"]
    result = dict(donors=donors, method=METHOD, created_at_UTC=now(),
        source_scripts=scripts, immutable_inputs=[dict(path=str(p), sha256=H.sha(p)) for p in inputs],
        original_reference_path=str(reference), original_reference_sha256=H.sha(reference),
        pinned_reference_copy=str(refcopy), pinned_reference_copy_sha256=H.sha(refcopy),
        original_reference_donor_ids={m: q["reference_donor"] for m, q in prior["shifts"].items()},
        geometry_root=str(out / "geometry"), current_block_sidecars=str(P.BLOCK),
        mandatory_assignments=str(P.ASSIGN), historical_diagnostic_root=str(OLD),
        threads=6, RSS_limit_GiB=180, accuracy_validated=False,
        calibration="Original reference donor IDs pinned; new own-source shifts and new own-source floor/top/span; no inherited Cell floor/span",
        block_scope="Same-block pancreas and attached pLN included; separate-block tissue on same slide excluded from this calibration",
        missingness="Unavailable measurements stay unavailable, with per-marker finite alignment populations",
        primary_CD8="Matching new own-source shifted values", refinement_CD8="Matching inverse multiplier recovers unshifted own-source CD8; inverse-converted floor/top/span; fresh fitted representation",
        producer_repairs=["finite observed declared negatives", "ordinary endocrine subtype abstention for incomplete/no-positive/tied hormone evidence", "direct Single policy preserved"],
        baseline_modified=False, release_promoted=False)
    H.save(config, result)
    return result

def check_sources(config):
    for row in [*config["source_scripts"], *config["immutable_inputs"]]:
        if H.sha(row["path"]) != row["sha256"]:
            raise ValueError(f"Frozen producer/input changed: {row['path']}")
    if H.sha(config["pinned_reference_copy"]) != config["pinned_reference_copy_sha256"]:
        raise ValueError("Pinned reference copy changed")
    if "production_reference_binding" in config:
        binding = config["production_reference_binding"]
        if H.sha(binding["path"]) != binding["sha256"]:
            raise ValueError("Selected production reference binding changed")
        receipt = json.loads(Path(binding["path"]).read_text())
        for row in receipt["immutable_binding_inputs"]:
            if H.sha(row["path"]) != row["sha256"]:
                raise ValueError(f"Selected reference/block/exclusion input changed: {row['path']}")

def prepare(out, donors):
    prior_path = OLD / "source_manifest.json"
    prior = json.loads(prior_path.read_text())
    selected = [r for r in prior["records"] if r["arm"] == METHOD]
    reuse = [r["donor_id_native"] for r in selected]
    if set(reuse) != set(P.ND) or len(reuse) != 9:
        raise ValueError("Unexpected diagnostic source reuse scope")
    remaining = [d for d in donors if d not in reuse]
    if len(remaining) != 16:
        raise ValueError("Expected sixteen newly prepared donors")
    # Prepare new sources with the same current selected-expression implementation.
    P.prepare(out, remaining, methods=[METHOD])
    newly = json.loads((out / "source_manifest.json").read_text())
    reuse_receipt = []
    for record in selected:
        donor = record["donor_id_native"]
        for row in [*record["source_files"].values(), record["assignment"]]:
            if H.sha(row["path"]) != row["sha256"]:
                raise ValueError(f"Reused original measurement/exclusion source changed: {row['path']}")
        row = dict(record)
        for name in ("expression", "availability"):
            source = Path(record[name])
            if H.sha(source) != record[name + "_sha256"]:
                raise ValueError(f"Reused own source changed: {source}")
            target = out / METHOD / name / f"donor_id={donor}/data.parquet"
            target.parent.mkdir(parents=True, exist_ok=True); target.symlink_to(source.resolve())
            row[name] = str(target); row[name + "_sha256"] = H.sha(target)
        original_geo = P.GEOMETRY / f"revised_cells/donor_id={donor}/data.parquet"
        prepared_geo = OLD / f"geometry/revised_cells/donor_id={donor}/data.parquet"
        cols = ["object_id", "qc_analysis_eligible", "qc_estimation_eligible"]
        baseline = P.read(original_geo, cols).set_index("object_id")
        prepared = P.read(prepared_geo, cols).set_index("object_id")
        denied = P.read(record["assignment"]["path"], ["object_id", "analysis_eligible"])
        bad = set(denied.loc[~denied.analysis_eligible.astype(bool), "object_id"])
        baseline.loc[baseline.index.isin(bad), cols[1:]] = False
        if not baseline.equals(prepared.reindex(baseline.index)):
            raise ValueError(f"Reused geometry/QC universe differs from current source: {donor}")
        target = out / f"geometry/revised_cells/donor_id={donor}/data.parquet"
        target.parent.mkdir(parents=True, exist_ok=True); target.symlink_to(prepared_geo.resolve())
        context = OLD / METHOD / f"typing/scratch/subpanel_context/{donor}__geometry.parquet"
        ctx = P.read(context, ["object_id"])
        eligible = set(prepared.index[prepared.qc_analysis_eligible.astype(bool)])
        if ctx.object_id.duplicated().any() or set(ctx.object_id) != eligible:
            raise ValueError(f"Reusable geometry context UUID universe differs: {donor}")
        ctx_target = out / METHOD / f"typing/scratch/subpanel_context/{donor}__geometry.parquet"
        ctx_target.parent.mkdir(parents=True, exist_ok=True); shutil.copy2(context, ctx_target)
        row["source_reuse"] = dict(original_manifest=str(prior_path), original_manifest_sha256=H.sha(prior_path),
            original_geometry=str(original_geo), original_geometry_sha256=H.sha(original_geo),
            reused_geometry=str(prepared_geo), geometry_sha256=H.sha(prepared_geo),
            geometry_context_source=str(context), geometry_context_sha256=H.sha(context),
            reused_marker_results=False, levels_recomputed=True)
        newly["records"].append(row); reuse_receipt.append(row["source_reuse"])
        print("VERIFIED OWN SOURCE/GEOMETRY REUSE", donor, flush=True)
    newly.update(donors=donors, arms=[METHOD], status="expression_complete",
                 source_reuse= reuse_receipt, generation="Corrected all25 own-source generation")
    newly["records"].sort(key=lambda r: donors.index(r["donor_id_native"]))
    old_qc = pd.read_csv(OLD / "expression_measurement_QC.csv", dtype={"donor_id": str})
    new_qc = pd.read_csv(out / "expression_measurement_QC.csv", dtype={"donor_id": str})
    pd.concat([new_qc, old_qc.loc[old_qc.arm.eq(METHOD)]], ignore_index=True).to_csv(out / "expression_measurement_QC.csv", index=False)
    H.save(out / "source_manifest.json", newly)

def alignment_preflight(out, donors):
    shifts_path = out / METHOD / "shifts.json"
    shifts = json.loads(shifts_path.read_text())
    pinned = json.loads((out / "reference_donor_identity_source.json").read_text())["shifts"]
    if shifts["donors"] != donors:
        raise ValueError("New shifts cohort differs")
    rows = []
    for donor in donors:
        source = out / METHOD / f"expression/donor_id={donor}/data.parquet"
        aligned = out / METHOD / f"uniform/donor_id={donor}/data.parquet"
        names = pq.ParquetFile(aligned).schema.names
        markers = [m for m in shifts["shifts"] if m in names]
        left = H.sample(source, ["object_id", *markers])
        right = H.sample(aligned, ["object_id", "in_alignment_block", *markers])
        pair = left.merge(right, on="object_id", validate="one_to_one", suffixes=("_source", "_aligned"))
        if len(pair) != len(left):
            raise ValueError(f"Aligned exact UUID sample changed: {donor}")
        inside = pair.in_alignment_block.to_numpy(bool)
        for marker in markers:
            entry = shifts["shifts"][marker]
            if entry["reference_donor"] != pinned[marker]["reference_donor"]:
                raise ValueError(f"Original reference donor changed: {marker}")
            multiplier = entry["scale"][donor]
            expected = pair[marker + "_source"].to_numpy(float).copy()
            if donor in entry.get("unavailable_donors", {}):
                expected[:] = np.nan
            elif multiplier != "unavailable":
                expected *= multiplier
            expected[~inside] = np.nan
            if not np.allclose(expected, pair[marker + "_aligned"].to_numpy(float), rtol=1e-7, atol=1e-8, equal_nan=True):
                raise ValueError(f"Aligned source/multiplier/missingness mismatch: {donor} {marker}")
        multiplier = float(shifts["shifts"]["CD8"]["scale"][donor])
        restored = pair.CD8_aligned.to_numpy(float) / multiplier
        expected = pair.CD8_source.to_numpy(float).copy(); expected[~inside] = np.nan
        if multiplier <= 0 or not np.isfinite(multiplier) or not np.allclose(restored, expected, rtol=1e-7, atol=1e-8, equal_nan=True):
            raise ValueError(f"Matching inverse CD8 does not recover own-source input: {donor}")
        rows.append(dict(donor_id=donor, sampled_exact_UUIDs=len(pair), n_markers=len(markers),
            original_reference_IDs_preserved=True, source_multiplier_and_block_mask_parity=True,
            inverse_CD8_matches_own_source=True, expression_sha256=H.sha(source), aligned_sha256=H.sha(aligned)))
    H.save(out / "alignment_and_CD8_preflight.json", dict(status="passed", checked_at_UTC=now(), rows=rows,
        shifts_sha256=H.sha(shifts_path), scale_sha256=H.sha(out / METHOD / "scale.parquet"),
        accuracy_validated=False, interpretation="Bounded measurement-unit, reference-identity and mask integrity"))

def prepared_input_preflight(out, donors):
    manifest = json.loads((out / "source_manifest.json").read_text())
    if manifest["donors"] != donors or manifest["arms"] != [METHOD] or len(manifest["records"]) != 25:
        raise ValueError("Prepared source manifest is not exact all25 own-only scope")
    checked = {}
    for record in manifest["records"]:
        entries = [*record["source_files"].values(), record["assignment"],
            *[dict(path=record[k], sha256=record[k + "_sha256"]) for k in ("expression", "availability")]]
        for entry in entries:
            path = entry["path"]
            if path in checked:
                if checked[path] != entry["sha256"]: raise ValueError(f"Conflicting prepared source hash: {path}")
                continue
            if H.sha(path) != entry["sha256"]: raise ValueError(f"Prepared/source/exclusion input changed: {path}")
            checked[path] = entry["sha256"]
    geometries = []
    for donor in donors:
        path = out / f"geometry/revised_cells/donor_id={donor}/data.parquet"
        geometries.append(dict(donor_id=donor, path=str(path), sha256=H.sha(path)))
    H.save(out / "prepared_input_integrity.json", dict(status="passed", checked_at_UTC=now(),
        donors=donors, exact_prepared_own_only_cohort=True,
        source_manifest_sha256=H.sha(out / "source_manifest.json"),
        original_and_prepared_hash_checks=[dict(path=p, sha256=h) for p, h in checked.items()],
        geometry_inputs=geometries, accuracy_validated=False))

def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--out", type=Path, default=DEFAULT_OUT)
    ap.add_argument("--through", choices=("prepare", "calibration", "primary"), default="prepare")
    ap.add_argument("--worker", choices=("prepare", "prepared_input_preflight", "alignment_preflight"))
    ap.add_argument("--initialize-only", action="store_true")
    ap.add_argument("--review-receipt", type=Path)
    args = ap.parse_args(); out = args.out.resolve()
    if not out.is_relative_to(Path("/data")):
        raise ValueError("Heavy new outputs must reside under /data")
    config = freeze(out); donors = config["donors"]
    if args.initialize_only:
        print(out / "generation_config.json"); return
    if args.worker:
        if args.worker == "prepare": prepare(out, donors)
        elif args.worker == "prepared_input_preflight": prepared_input_preflight(out, donors)
        else: alignment_preflight(out, donors)
        return
    env = os.environ.copy(); typing = out / METHOD / "typing"; typing.mkdir(parents=True, exist_ok=True)
    env.update(IO60_EXPRESSION_ROOT=str(out / METHOD / "expression"), IO60_ALIGNED_ROOT=str(out / METHOD / "uniform"),
        IO60_SCALE_PATH=str(out / METHOD / "scale.parquet"), IO60_TYPING_ROOT=str(typing),
        IO60_GEOMETRY_ROOT=str(out / "geometry"), IO60_SUBPANEL_LOG=str(P.GEOMETRY / "subpanel/subpanel_log.json"),
        IO60_SHIFTS_PATH=str(out / METHOD / "shifts.json"), IO60_REFINE_CD8_ALIGNMENT="unshifted")
    ledger = out / "run.json"
    record = json.loads(ledger.read_text()) if ledger.exists() else dict(donors=donors, stages=[], config=str(out / "generation_config.json"), config_sha256=H.sha(out / "generation_config.json"), started_at_UTC=now(), accuracy_validated=False)
    if record["config_sha256"] != H.sha(out / "generation_config.json"):
        raise ValueError("Generation configuration changed")
    stages = [
        ("prepare_sources", [sys.executable, str(Path(__file__)), "--out", str(out), "--worker", "prepare"], [out / "source_manifest.json", out / "expression_measurement_QC.csv"]),
        ("source_parity", [sys.executable, str(HERE / "verify_protein_candidate_expression_sources.py"), "--candidate-root", str(out)], [out / "resume_expression_source_parity_QC.json", out / "resume_expression_source_parity_QC.csv"])]
    da = ["--donors", *donors]
    if args.through in ("calibration", "primary"):
        stages += [
            ("prepared_input_preflight", [sys.executable, str(Path(__file__)), "--out", str(out), "--worker", "prepared_input_preflight"], [out / "prepared_input_integrity.json"]),
            ("uniform_and_new_scale", [sys.executable, str(HERE / "reprocess_protein_own_compartment_candidate.py"), "--worker", "uniform", "--out", str(out), "--arm", METHOD, "--current-cohort", "--references-from", config["pinned_reference_copy"], *da], [out / METHOD / "shifts.json", out / METHOD / "scale.parquet", out / METHOD / "scale_provenance.json", out / METHOD / "unavailable_alignment_histograms.json", out / METHOD / "uniform"]),
            ("alignment_preflight", [sys.executable, str(Path(__file__)), "--out", str(out), "--worker", "alignment_preflight"], [out / "alignment_and_CD8_preflight.json"])]
    if args.through == "primary":
        if not args.review_receipt:
            raise ValueError("Independent generation review receipt required before corrected primary launch")
        review = json.loads(args.review_receipt.read_text())
        if review.get("status") != "passed" or review.get("generation_config_sha256") != H.sha(out / "generation_config.json"):
            raise ValueError("Independent review does not approve this exact generation configuration")
        record["independent_review"] = dict(path=str(args.review_receipt), sha256=H.sha(args.review_receipt))
        stages += [
            ("subpanel", [sys.executable, str(P.TYPING / "subpanel.py"), *da, "--membership", "geometry", "--emit-detail", "--work", str(typing / "scratch/subpanel_context"), "--audit-out", str(typing / "subpanel_audit.json"), "--out", str(typing / "subpanel_rebuilt")], [typing / "subpanel_rebuilt", typing / "subpanel_audit.json"]),
            ("build_cells_primary", [sys.executable, str(P.TYPING / "build_cells.py"), *da, "--calls", str(typing / "subpanel_rebuilt"), "--out", str(typing / "cells")], [typing / "cells"]),
            ("primary_key_integrity", [sys.executable, str(HERE / "verify_protein_candidate_primary_integrity.py"), "--candidate-root", str(out)], [out / "resume_primary_integrity_QC.json", out / "resume_primary_integrity_QC.csv"]),
            ("primary_evidence_accounting", [sys.executable, str(HERE / "audit_protein_repaired_primary_evidence.py"), "--root", str(out)], [out / "repaired_primary_evidence_QC.json", out / "repaired_primary_donor_evidence_QC.csv", out / "repaired_primary_stratum_evidence_QC.csv", out / "repaired_primary_routes_statuses_QC.csv"]),
            ("primary_negative_evidence_integrity", [sys.executable, str(HERE / "audit_protein_candidate_missing_negatives.py"), "--candidate-root", str(out), "--assert-no-unavailable-negatives"], [out / "resume_missing_negative_primary_definition_QC.json", out / "resume_missing_negative_primary_definition_counts.csv", out / "resume_missing_negative_primary_definition_exact_keys.parquet"])]
    record.update(status="running", through=args.through, controller_pid=os.getpid(), RSS_limit_GiB=180, threads=6)
    H.save(ledger, record)
    for name, command, outputs in stages:
        check_sources(config)
        prior = next((s for s in record["stages"] if s["key"] == name and s.get("status") == "complete"), None)
        if prior:
            for item in prior["outputs"]:
                if H.sha(item["path"]) != item["sha256"]: raise ValueError(f"Completed artifact changed: {item['path']}")
            print("VERIFIED COMPLETE", name, flush=True); continue
        if any(s["key"] == name for s in record["stages"]):
            raise ValueError(f"A failed stage exists: {name}; preserve and inspect before a new generation")
        stage = dict(key=name, started_at_UTC=now()); record["current_stage"] = name; H.save(ledger, record)
        print("START", name, flush=True)
        try:
            stage.update(H.run(command, env, out / (name + ".log"), 180))
            if stage["returncode"]: raise RuntimeError(f"{name} exited {stage['returncode']}")
            stage.update(status="complete", finished_at_UTC=now(), outputs=H.artifacts(outputs))
        except Exception as exc:
            stage.update(status="failed", error=str(exc), finished_at_UTC=now()); record["stages"].append(stage)
            record.update(status="failed"); H.save(ledger, record); raise
        record["stages"].append(stage); H.save(ledger, record)
        if name == "build_cells_primary":
            H.save(typing / "run.json", dict(status="primary_complete", donors=donors,
                source_scripts=config["source_scripts"], stages=record["stages"],
                accuracy_validated=False, baseline_unchanged=True, integrity_QA_pending=True))
        print("COMPLETE", name, "peak GiB", stage["peak_process_group_RSS_GiB"], flush=True)
    record.update(status=args.through + "_complete", current_stage=None, finished_at_UTC=now()); H.save(ledger, record)
    H.save(typing / "run.json", dict(status=args.through + "_complete", donors=donors, source_scripts=config["source_scripts"], stages=record["stages"], accuracy_validated=False, baseline_unchanged=True))
    print("GENERATION COMPLETE THROUGH", args.through, flush=True)

if __name__ == "__main__":
    main()
