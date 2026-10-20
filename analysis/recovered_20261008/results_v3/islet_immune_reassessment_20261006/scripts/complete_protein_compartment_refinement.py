#!/usr/bin/env python3
"""Complete diagnostic refinement with matching unshifted-CD8 measurement units.

Preserves primary calls, fits fresh representation parameters and retains the
candidate/diagnostic status. It cannot promote these nine-ND comparisons.
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
import pyarrow.parquet as pq

HERE = Path(__file__).resolve().parent
spec = importlib.util.spec_from_file_location("resume_compartment", HERE / "resume_protein_compartment_candidate.py")
H = importlib.util.module_from_spec(spec); spec.loader.exec_module(H)
P = H.P

def now():
    return datetime.now(timezone.utc).isoformat()

def cd8_preflight(candidate, donors):
    rows = []
    for method in H.ARMS:
        shifts_path = candidate / method / "shifts.json"
        shifts = json.loads(shifts_path.read_text())
        for donor in donors:
            source = candidate / method / f"expression/donor_id={donor}/data.parquet"
            aligned = candidate / method / f"uniform/donor_id={donor}/data.parquet"
            raw = H.sample(source, ["object_id", "CD8"])
            values = H.sample(aligned, ["object_id", "CD8", "in_alignment_block"])
            pair = raw.merge(values, on="object_id", validate="one_to_one", suffixes=("_source", "_aligned"))
            if len(pair) != len(raw):
                raise ValueError(f"CD8 preflight UUID mismatch: {method} {donor}")
            multiplier = float(shifts["shifts"]["CD8"]["scale"][donor])
            if not np.isfinite(multiplier) or multiplier <= 0:
                raise ValueError(f"CD8 matching multiplier invalid: {method} {donor}")
            expected = pair.CD8_source.to_numpy(float).copy()
            expected[~pair.in_alignment_block.to_numpy(bool)] = np.nan
            restored = pair.CD8_aligned.to_numpy(float) / multiplier
            if not np.allclose(restored, expected, rtol=1e-7, atol=1e-8, equal_nan=True):
                raise ValueError(f"CD8 inverse does not recover own-method expression: {method} {donor}")
            rows.append(dict(method=method, donor_id=donor, sampled_exact_UUIDs=len(pair),
                inverse_multiplier=1 / multiplier, preserved_missingness=True, matches_own_method_expression=True,
                expression_sha256=H.sha(source), aligned_sha256=H.sha(aligned), shifts_sha256=H.sha(shifts_path)))
    H.save(candidate / "resume_CD8_own_method_refinement_preflight.json", dict(status="passed", checked_at_UTC=now(),
        rows=rows, accuracy_validated=False, interpretation="Matching measurement-unit and block-mask verification; no biological validation"))

def primary_copy(path, donors):
    source = path / "cells"; target = path / "cells_primary_preserved"
    if not target.exists():
        shutil.copytree(source, target)
    for donor in donors:
        left = source / f"donor_id={donor}/data.parquet"
        right = target / f"donor_id={donor}/data.parquet"
        if H.sha(left) != H.sha(right):
            raise ValueError(f"Primary-preservation mismatch: {donor}")
    return target

def check_join(path, donors):
    rows = []
    for donor in donors:
        primary = path / f"cells_primary_preserved/donor_id={donor}/data.parquet"
        joined = path / f"cells/donor_id={donor}/data.parquet"
        columns = ["object_id", "cell_class", "qc_analysis_eligible", "endocrine_route", "is_endocrine"]
        left = pq.ParquetFile(primary).read(columns=columns).to_pandas()
        right = pq.ParquetFile(joined).read(columns=columns + ["refined_class"]).to_pandas()
        pair = left.merge(right, on="object_id", how="outer", validate="one_to_one", suffixes=("_primary", "_joined"), indicator=True)
        if not pair._merge.eq("both").all():
            raise ValueError(f"Refinement join changed UUID universe: {donor}")
        for column in columns[1:]:
            if not pair[column + "_primary"].fillna("").equals(pair[column + "_joined"].fillna("")):
                raise ValueError(f"Refinement join changed primary {column}: {donor}")
        if right.refined_class.isna().any():
            raise ValueError(f"Null refined_class: {donor}")
        rows.append(dict(donor_id=donor, n_exact_UUIDs=len(left), primary_fields_unchanged=True,
                         primary_sha256=H.sha(primary), joined_sha256=H.sha(joined)))
    H.save(path / "primary_refinement_join_integrity.json", dict(status="passed", rows=rows,
        interpretation="Primary identity, route and QC preservation across keyed refinement join; biological accuracy unknown"))

def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--candidate-root", type=Path, default=P.R / "cell_typing_validation/candidate_compartment_reprocessing")
    ap.add_argument("--preflight-only", action="store_true")
    ap.add_argument("--rss-limit-gib", type=float, default=180)
    args = ap.parse_args(); candidate = args.candidate_root.resolve()
    manifest = json.loads((candidate / "source_manifest.json").read_text()); donors = manifest["donors"]
    if set(donors) != set(P.ND) or len(donors) != len(P.ND):
        raise ValueError("Refinement continuation is restricted to the exact prepared nine ND donors")
    for method in H.ARMS:
        path = candidate / method / "typing"
        if not path.resolve().is_relative_to(Path("/data")):
            raise ValueError("Large refinement outputs must reside on /data")
        if json.loads((path / "run.json").read_text())["status"] != "primary_complete":
            raise ValueError(f"Primary processing is not complete or already refined: {method}")
    cd8_preflight(candidate, donors)
    print("MATCHING OWN-METHOD UNSHIFTED CD8 PREFLIGHT PASSED", flush=True)
    if args.preflight_only:
        return
    ledger = candidate / "resume_refinement_run_20261008.json"
    if ledger.exists():
        raise FileExistsError("A refinement continuation already exists; inspect its actual process/output state before resuming")
    scripts = [Path(__file__), P.TYPING / "refine.py", P.TYPING / "build_cells.py", P.TYPING / "immune_subtype.py",
               P.TYPING / "subpanel_io.py", P.CELLTYPE / "score.py", P.CELLTYPE / "rules.py"]
    hashes = {str(path): H.sha(path) for path in scripts}
    record = dict(status="running", controller_pid=os.getpid(), donors=donors, started_at_UTC=now(),
        source_scripts=hashes, stages=[], threads=6, RSS_limit_GiB=args.rss_limit_gib,
        primary_CD8="matching shifted own-method values", refinement_CD8="unshifted own-method values",
        refinement_immune_CD8_scale="matching inverse donor multiplier applied to frozen floor/top/span",
        fresh_representation=True, baseline_modified=False, accuracy_validated=False,
        scope="Nine ND diagnostic comparison; all25 frozen-reference release remains separate")
    H.save(ledger, record)
    for method in H.ARMS:
        path = candidate / method / "typing"; primary = primary_copy(path, donors)
        shutil.copy2(path / "run.json", path / "run_primary_preserved.json")
        for name in ("refinement", "cells_with_refinement", "immune"):
            if (path / name).exists():
                raise FileExistsError(f"Pending refinement stage already has output: {path / name}")
        env = os.environ.copy(); shifts = candidate / method / "shifts.json"
        env.update(IO60_EXPRESSION_ROOT=str(candidate / method / "expression"), IO60_ALIGNED_ROOT=str(candidate / method / "uniform"),
            IO60_SCALE_PATH=str(candidate / method / "scale.parquet"), IO60_TYPING_ROOT=str(path),
            IO60_GEOMETRY_ROOT=str(candidate / "geometry"), IO60_SUBPANEL_LOG=str(P.GEOMETRY / "subpanel/subpanel_log.json"),
            IO60_SHIFTS_PATH=str(shifts), IO60_REFINE_CD8_ALIGNMENT="unshifted")
        da = ["--donors", *donors]
        stages = [
            ("refine", [sys.executable, str(P.TYPING / "refine.py"), *da, "--out", str(path / "refinement"),
                        "--cd8-alignment", "unshifted", "--shifts", str(shifts)], [path / "refinement"]),
            ("build_cells_with_refinement", [sys.executable, str(P.TYPING / "build_cells.py"), *da,
                        "--calls", str(path / "subpanel_rebuilt"), "--refine", str(path / "refinement"),
                        "--audit", str(path / "subpanel_audit.json"), "--out", str(path / "cells_with_refinement")], [path / "cells_with_refinement"]),
            ("immune_subtype", [sys.executable, str(P.TYPING / "immune_subtype.py"), *da], [path / "immune"])]
        for name, command, outputs in stages:
            for source_path, expected in hashes.items():
                if H.sha(source_path) != expected:
                    raise ValueError(f"Refinement producer changed during run: {source_path}")
            key = method + ":" + name
            record["current_stage"] = key; H.save(ledger, record)
            H.save(path / "run.json", dict(status="running", donors=donors, current_stage=key,
                controller_pid=os.getpid(), source_scripts=hashes, primary_preserved=str(primary), accuracy_validated=False))
            print("START", key, flush=True)
            stage = dict(key=key, started_at_UTC=now())
            try:
                stage.update(H.run(command, env, path / (name + "_unshifted_CD8_20261008.log"), args.rss_limit_gib))
                if stage["returncode"]:
                    raise RuntimeError(f"Stage failed: {key} returncode={stage['returncode']}")
                stage.update(status="complete", finished_at_UTC=now(), outputs=H.artifacts(outputs))
                if name == "refine":
                    receipt = json.loads((path / "refinement/representation_provenance.json").read_text())
                    if receipt["fit_state"] != "fresh" or receipt["fit_donors"] != sorted(donors) or receipt["cd8_alignment"] != "unshifted":
                        raise ValueError(f"Refinement representation does not match fresh prepared cohort: {method}")
                    stage["representation_provenance"] = receipt
                elif name == "build_cells_with_refinement":
                    # Keep the original primary directory and make its active logical
                    # path refer to the final joined output before immune subtyping.
                    (path / "cells").rename(path / "cells_primary_original")
                    (path / "cells").symlink_to(path / "cells_with_refinement", target_is_directory=True)
                    check_join(path, donors)
            except Exception as exc:
                stage.update(status="failed", error=str(exc), finished_at_UTC=now())
                record["stages"].append(stage); record.update(status="failed")
                H.save(ledger, record); H.save(path / "run.json", dict(status="failed", donors=donors, stage=stage))
                raise
            record["stages"].append(stage); H.save(ledger, record)
            print("COMPLETE", key, "peak GiB", stage["peak_process_group_RSS_GiB"], flush=True)
        H.save(path / "run.json", dict(status="complete", donors=donors, source_scripts=hashes,
            stages=[s for s in record["stages"] if s["key"].startswith(method + ":")],
            primary_source_snapshot=str(candidate / "diagnostic_primary_source_snapshot.json"),
            primary_preserved=str(primary), baseline_unchanged=True, accuracy_validated=False, scope=record["scope"]))
    record.update(status="complete", finished_at_UTC=now(), current_stage=None); H.save(ledger, record)
    H.save(candidate / "comparison_status.json", dict(status="processing_complete", donors=donors,
        primary_producer="preserved pre-repair missing-negative diagnostic", refinement_CD8="unshifted with matching units",
        candidate_release_selected=False, baseline_unchanged=True, accuracy_validated=False,
        interpretation="Diagnostic comparison only; declared negative-gap audit and all25 frozen-reference source repair remain required"))

if __name__ == "__main__":
    main()
