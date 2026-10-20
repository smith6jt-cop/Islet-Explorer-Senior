#!/usr/bin/env python3
"""Finish corrected all25 own-source typing with fresh unshifted-CD8 refinement."""
import os
for name in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMBA_NUM_THREADS"):
    os.environ[name] = "6"
from pathlib import Path
import argparse
import importlib.util
import json
import shutil
import sys

HERE = Path(__file__).resolve().parent
spec = importlib.util.spec_from_file_location("all25_primary", HERE / "run_protein_own_frozen_reference_repair.py")
C = importlib.util.module_from_spec(spec); spec.loader.exec_module(C)
spec = importlib.util.spec_from_file_location("diagnostic_refine_helpers", HERE / "complete_protein_compartment_refinement.py")
D = importlib.util.module_from_spec(spec); spec.loader.exec_module(D)
H = C.H; P = C.P

def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--out", type=Path, default=C.DEFAULT_OUT)
    ap.add_argument("--preflight-only", action="store_true")
    args = ap.parse_args(); out = args.out.resolve()
    config = C.freeze(out); donors = config["donors"]; C.check_sources(config)
    typing = out / C.METHOD / "typing"
    if json.loads((out / "run.json").read_text())["status"] != "primary_complete":
        raise ValueError("All25 corrected primary stages and integrity QA are not complete")
    for name in ("resume_primary_integrity_QC.json", "repaired_primary_evidence_QC.json"):
        if json.loads((out / name).read_text())["status"] != "passed":
            raise ValueError(f"Primary integrity receipt did not pass: {name}")
    negative = json.loads((out / "resume_missing_negative_primary_definition_QC.json").read_text())
    if negative["n_flagged_exact_keys"] != 0:
        raise ValueError("Primary definition calls retain unavailable declared negative evidence")
    # Verify the immutable calibration receipt instead of rewriting its timestamp
    # and invalidating the primary controller's recorded artifact checksum.
    alignment = json.loads((out / "alignment_and_CD8_preflight.json").read_text())
    if alignment["status"] != "passed" or alignment["shifts_sha256"] != H.sha(out / C.METHOD / "shifts.json") or alignment["scale_sha256"] != H.sha(out / C.METHOD / "scale.parquet"):
        raise ValueError("Recorded own-source alignment/unit receipt is absent or mismatched")
    if {r["donor_id"] for r in alignment["rows"]} != set(donors):
        raise ValueError("Recorded alignment receipt omits a current donor")
    for row in alignment["rows"]:
        donor = row["donor_id"]
        for directory, key in (("expression", "expression_sha256"), ("uniform", "aligned_sha256")):
            if H.sha(out / C.METHOD / directory / f"donor_id={donor}/data.parquet") != row[key]:
                raise ValueError(f"Verified own-source refinement input changed: {donor} {directory}")
    finite_units = out / "CD8_finite_own_source_unit_parity.json"
    receipt = json.loads(finite_units.read_text())
    if receipt["status"] != "passed" or receipt["n_exact_finite_UUIDs"] != 25 * 2048 or receipt["shifts_sha256"] != H.sha(out / C.METHOD / "shifts.json"):
        raise ValueError("Independent finite in-block own-source CD8 parity is absent or mismatched")
    print("ALL25 PRIMARY INTEGRITY AND OWN-SOURCE CD8 PREFLIGHT PASSED", flush=True)
    if args.preflight_only: return
    ledger = out / "refinement_run.json"
    if ledger.exists():
        raise FileExistsError("Refinement generation already exists; inspect actual processes and preserve partial artifacts before retry")
    for name in ("refinement", "cells_with_refinement", "immune"):
        if (typing / name).exists(): raise FileExistsError(typing / name)
    initial = typing / "cells"
    preserved = typing / "cells_primary_preserved"
    shutil.copytree(initial, preserved)
    for donor in donors:
        if H.sha(initial / f"donor_id={donor}/data.parquet") != H.sha(preserved / f"donor_id={donor}/data.parquet"):
            raise ValueError(f"Preserved primary file differs: {donor}")
    shutil.copy2(typing / "run.json", typing / "run_primary_preserved.json")
    snapshot = out / "refinement_controller_source_snapshot"; snapshot.mkdir()
    sources = []
    for i, source in enumerate((Path(__file__), HERE / "complete_protein_compartment_refinement.py", HERE / "verify_protein_own_CD8_finite_units.py")):
        target = snapshot / f"{i:02d}_{source.name}"; shutil.copy2(source, target)
        sources.append(dict(path=str(source), sha256=H.sha(source), snapshot=str(target)))
    record = dict(status="running", donors=donors, controller_pid=os.getpid(), started_at_UTC=C.now(),
        source_scripts=sources, generation_config_sha256=H.sha(out / "generation_config.json"),
        primary_preserved=str(preserved), stages=[], threads=6, RSS_limit_GiB=180,
        fresh_representation=True, primary_CD8="matching shifted own-source units",
        refinement_CD8="unshifted own-source units with matching inverse-converted CD8 floor/top/span",
        baseline_modified=False, release_promoted=False, accuracy_validated=False)
    H.save(ledger, record)
    shifts = out / C.METHOD / "shifts.json"
    env = os.environ.copy(); env.update(IO60_EXPRESSION_ROOT=str(out / C.METHOD / "expression"),
        IO60_ALIGNED_ROOT=str(out / C.METHOD / "uniform"), IO60_SCALE_PATH=str(out / C.METHOD / "scale.parquet"),
        IO60_TYPING_ROOT=str(typing), IO60_GEOMETRY_ROOT=str(out / "geometry"),
        IO60_SUBPANEL_LOG=str(P.GEOMETRY / "subpanel/subpanel_log.json"),
        IO60_SHIFTS_PATH=str(shifts), IO60_REFINE_CD8_ALIGNMENT="unshifted")
    da = ["--donors", *donors]
    stages = [
        ("refine", [sys.executable, str(P.TYPING / "refine.py"), *da, "--out", str(typing / "refinement"), "--cd8-alignment", "unshifted", "--shifts", str(shifts)], [typing / "refinement"]),
        ("build_cells_with_refinement", [sys.executable, str(P.TYPING / "build_cells.py"), *da, "--calls", str(typing / "subpanel_rebuilt"), "--refine", str(typing / "refinement"), "--audit", str(typing / "subpanel_audit.json"), "--out", str(typing / "cells_with_refinement")], [typing / "cells_with_refinement"]),
        ("immune_subtype", [sys.executable, str(P.TYPING / "immune_subtype.py"), *da], [typing / "immune"])]
    for name, command, outputs in stages:
        C.check_sources(config)
        for source in sources:
            if H.sha(source["path"]) != source["sha256"]: raise ValueError(f"Finalization controller changed: {source['path']}")
        stage = dict(key=name, started_at_UTC=C.now()); record["current_stage"] = name; H.save(ledger, record)
        H.save(typing / "run.json", dict(status="running", donors=donors, controller_pid=os.getpid(),
            current_stage=name, primary_preserved=str(preserved), accuracy_validated=False))
        print("START", name, flush=True)
        try:
            stage.update(H.run(command, env, typing / (name + "_all25_unshifted_CD8.log"), 180))
            if stage["returncode"]: raise RuntimeError(f"{name} exited {stage['returncode']}")
            stage.update(status="complete", finished_at_UTC=C.now(), outputs=H.artifacts(outputs))
            if name == "refine":
                receipt = json.loads((typing / "refinement/representation_provenance.json").read_text())
                if receipt["fit_state"] != "fresh" or receipt["fit_donors"] != sorted(donors) or receipt["cd8_alignment"] != "unshifted":
                    raise ValueError("Representation was not fitted freshly on exact all25 unshifted-CD8 inputs")
                stage["representation_provenance"] = receipt
            elif name == "build_cells_with_refinement":
                initial.rename(typing / "cells_primary_original")
                initial.symlink_to(typing / "cells_with_refinement", target_is_directory=True)
                D.check_join(typing, donors)
        except Exception as exc:
            stage.update(status="failed", error=str(exc), finished_at_UTC=C.now()); record["stages"].append(stage)
            record["status"] = "failed"; H.save(ledger, record)
            H.save(typing / "run.json", dict(status="failed", donors=donors, stage=stage, accuracy_validated=False)); raise
        record["stages"].append(stage); H.save(ledger, record)
        print("COMPLETE", name, "peak GiB", stage["peak_process_group_RSS_GiB"], flush=True)
    record.update(status="complete", current_stage=None, finished_at_UTC=C.now()); H.save(ledger, record)
    H.save(typing / "run.json", dict(status="complete", donors=donors, stages=record["stages"],
        primary_preserved=str(preserved), generation_config_sha256=record["generation_config_sha256"],
        source_scripts={q["path"]: q["sha256"] for q in [*config["source_scripts"], *sources]},
        primary_generation_run=str(out / "run.json"), primary_generation_run_sha256=H.sha(out / "run.json"),
        refinement_generation_run=str(ledger), refinement_generation_run_sha256=H.sha(ledger),
        geometry_root=str(out / "geometry"), original_geometry_root=str(P.GEOMETRY),
        aligned_root=str(out / C.METHOD / "uniform"), scale=str(out / C.METHOD / "scale.parquet"),
        aligned_files={d: H.sha(out / C.METHOD / f"uniform/donor_id={d}/data.parquet") for d in donors},
        scale_sha256=H.sha(out / C.METHOD / "scale.parquet"), expression_root=str(out / C.METHOD / "expression"),
        registry_criteria=str(P.GEOMETRY / "subpanel/subpanel_log.json"),
        registry_criteria_sha256=H.sha(P.GEOMETRY / "subpanel/subpanel_log.json"),
        source_manifest=str(out / "source_manifest.json"), source_manifest_sha256=H.sha(out / "source_manifest.json"),
        mandatory_assignment_manifest=str(P.ASSIGN / "manifest.json"),
        mandatory_assignment_manifest_sha256=H.sha(P.ASSIGN / "manifest.json"),
        baseline_unchanged=True, accuracy_validated=False, release_promoted=False))
    print("ALL25 CORRECTED PROTEIN PROCESSING COMPLETE; BIOLOGICAL VALIDATION REMAINS SEPARATE", flush=True)

if __name__ == "__main__":
    main()
