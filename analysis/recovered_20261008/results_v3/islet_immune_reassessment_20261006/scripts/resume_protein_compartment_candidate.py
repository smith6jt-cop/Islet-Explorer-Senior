#!/usr/bin/env python3
"""Resume the isolated nine-ND diagnostic comparison with explicit verified stages.

Preserves interrupted outputs and primary calls. This does not select a release or
validate biological accuracy. Newly written large outputs reside on /data.
"""
import os
for name in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMBA_NUM_THREADS"):
    os.environ[name] = "6"
from pathlib import Path
import argparse
import hashlib
import importlib.util
import json
import shutil
import signal
import subprocess
import sys
import time
from datetime import datetime, timezone
import numpy as np
import pyarrow.parquet as pq

HERE = Path(__file__).resolve().parent
spec = importlib.util.spec_from_file_location("compartment_candidate", HERE / "reprocess_protein_own_compartment_candidate.py")
P = importlib.util.module_from_spec(spec)
spec.loader.exec_module(P)
ARMS = ("Cell_control", "own_compartment")

def now():
    return datetime.now(timezone.utc).isoformat()

def sha(path):
    return P.sha(path)

def save(path, value):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_name(path.name + ".tmp")
    tmp.write_text(json.dumps(value, indent=2, allow_nan=False) + "\n")
    tmp.replace(path)

def sample(path, columns, n=2048):
    return next(pq.ParquetFile(path).iter_batches(batch_size=n, columns=columns)).to_pandas()

def preflight(candidate):
    manifest = json.loads((candidate / "source_manifest.json").read_text())
    donors = manifest["donors"]
    if set(donors) != set(P.ND) or len(donors) != len(P.ND):
        raise ValueError("This diagnostic comparison requires exactly the prepared nine ND donors")
    checks = []
    for path, expected in [(P.POLICY, manifest["policy_sha256"]),
                           (P.REGISTRY, manifest["baseline_registry_sha256"]),
                           (P.BLOCK / "manifest.json", manifest["current_block_manifest_sha256"])]:
        actual = sha(path)
        checks.append(dict(path=str(path), expected_sha256=expected, actual_sha256=actual))
        if actual != expected:
            raise ValueError(f"Prepared input changed: {path}")
    # The original manifest recorded hashes before the interrupted processing.
    for record in manifest["records"]:
        for key in ("expression", "availability"):
            path = Path(record[key]); actual = sha(path)
            expected = record[key + "_sha256"]
            checks.append(dict(path=str(path), expected_sha256=expected, actual_sha256=actual))
            if actual != expected:
                raise ValueError(f"Prepared {key} changed: {path}")
        # Both copies of each source record must still reference the exact immutable
        # measurement and mandatory-exclusion inputs from which they were built.
        if record["arm"] == ARMS[0]:
            for entry in [*record["source_files"].values(), record["assignment"]]:
                path = Path(entry["path"]); actual = sha(path)
                checks.append(dict(path=str(path), expected_sha256=entry["sha256"], actual_sha256=actual))
                if actual != entry["sha256"]:
                    raise ValueError(f"Original input changed: {path}")
    # Check that the completed Cell calibration actually applies its recorded
    # multiplier and masks other-block cells; file existence alone is insufficient.
    arm = candidate / ARMS[0]
    receipt = json.loads((arm / "scale_provenance.json").read_text())
    if receipt["donors"] != donors:
        raise ValueError("Completed Cell calibration uses a different donor cohort")
    for key in ("shifts", "scale"):
        path = arm / (key + (".json" if key == "shifts" else ".parquet"))
        if sha(path) != receipt[key + "_sha256"]:
            raise ValueError(f"Completed Cell calibration changed: {path}")
    shifts = json.loads((arm / "shifts.json").read_text())["shifts"]
    parity = []
    for donor in donors:
        source = arm / f"expression/donor_id={donor}/data.parquet"
        aligned = arm / f"uniform/donor_id={donor}/data.parquet"
        names = pq.ParquetFile(aligned).schema.names
        markers = [m for m in shifts if m in names]
        left = sample(source, ["object_id"] + markers)
        right = sample(aligned, ["object_id", "in_alignment_block"] + markers)
        joined = left.merge(right, on="object_id", validate="one_to_one", suffixes=("_source", "_aligned"))
        if len(joined) != len(left):
            raise ValueError(f"Completed Cell alignment UUID sample changed: {donor}")
        inside = joined.in_alignment_block.to_numpy(bool)
        for marker in markers:
            multiplier = shifts[marker]["scale"][donor]
            expected = joined[marker + "_source"].to_numpy(float).copy()
            if multiplier != "unavailable":
                expected *= multiplier
            expected[~inside] = np.nan
            observed = joined[marker + "_aligned"].to_numpy(float)
            if not np.allclose(expected, observed, rtol=1e-7, atol=1e-8, equal_nan=True):
                raise ValueError(f"Completed Cell alignment multiplier/missingness mismatch: {donor} {marker}")
        parity.append(dict(donor_id=donor, sampled_UUIDs=len(joined), markers_checked=len(markers),
                           aligned_sha256=sha(aligned), exact_recorded_multiplier_and_missingness=True))
    result = dict(status="passed", checked_at_UTC=now(), donors=donors,
                  hash_checks=checks, sampled_alignment_parity=parity,
                  accuracy_validated=False, interpretation="Input and numerical integrity only")
    save(candidate / "resume_preflight_20261008.json", result)
    return donors

def redirect(candidate, storage, arm, name):
    target = candidate / arm / name
    destination = storage / arm / name
    if target.is_symlink():
        if target.resolve() != destination.resolve():
            raise ValueError(f"Unexpected candidate storage link: {target}")
        return
    destination.parent.mkdir(parents=True, exist_ok=True)
    if target.exists():
        archive = storage / "interrupted_20261006" / arm / name
        if archive.exists():
            raise FileExistsError(archive)
        archive.parent.mkdir(parents=True, exist_ok=True)
        shutil.move(str(target), str(archive))
    if destination.exists():
        raise FileExistsError(destination)
    destination.mkdir(parents=True)
    target.symlink_to(destination, target_is_directory=True)

def rss_group(pgid):
    total = 0
    for proc in Path("/proc").iterdir():
        if not proc.name.isdecimal():
            continue
        try:
            tail = (proc / "stat").read_text().rsplit(")", 1)[1].split()
            if int(tail[2]) != pgid:
                continue
            total += int(tail[21]) * os.sysconf("SC_PAGE_SIZE")
        except (OSError, ValueError, IndexError):
            continue
    return total

def run(command, env, log, rss_limit_gib):
    peak = 0
    with log.open("w") as stream:
        process = subprocess.Popen(command, cwd=P.ROOT, env=env, stdout=stream,
                                   stderr=subprocess.STDOUT, start_new_session=True)
        while process.poll() is None:
            peak = max(peak, rss_group(process.pid))
            if peak > rss_limit_gib * (1 << 30):
                os.killpg(process.pid, signal.SIGTERM)
                try:
                    process.wait(timeout=15)
                except subprocess.TimeoutExpired:
                    os.killpg(process.pid, signal.SIGKILL)
                    process.wait()
                raise RuntimeError(f"Stage RSS exceeded {rss_limit_gib} GiB; preserved log {log}")
            time.sleep(2)
    return dict(command=command, log=str(log), returncode=process.returncode,
                peak_process_group_RSS_GiB=round(peak / (1 << 30), 3))

def artifacts(paths):
    found = []
    for root in paths:
        files = sorted(root.rglob("*.parquet")) if root.is_dir() else [root]
        for path in files:
            if not path.is_file():
                raise FileNotFoundError(path)
            found.append(dict(path=str(path), bytes=path.stat().st_size, sha256=sha(path)))
    if not found:
        raise ValueError("Stage produced no recorded artifacts")
    return found

def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--candidate-root", type=Path, default=P.R / "cell_typing_validation/candidate_compartment_reprocessing")
    ap.add_argument("--storage-root", type=Path, default=Path("/data/io60-calibration/candidate_compartment_reprocessing_resume_20261008"))
    ap.add_argument("--through", choices=("calibration", "primary", "immune"), default="primary")
    ap.add_argument("--preflight-only", action="store_true")
    ap.add_argument("--rss-limit-gib", type=float, default=180)
    args = ap.parse_args(); candidate = args.candidate_root.resolve(); storage = args.storage_root.resolve()
    if not storage.is_relative_to(Path("/data")):
        raise ValueError("Large resumed outputs must be under /data")
    donors = preflight(candidate)
    print("PREPARED INPUT AND CELL ALIGNMENT PREFLIGHT PASSED", flush=True)
    if args.preflight_only:
        return
    redirect(candidate, storage, ARMS[1], "uniform")
    for arm in ARMS:
        redirect(candidate, storage, arm, "typing")
    source_files = [HERE / "reprocess_protein_own_compartment_candidate.py", Path(__file__),
                    *[P.TYPING / name for name in ("subpanel.py", "build_cells.py", "refine.py", "immune_subtype.py", "subpanel_io.py")],
                    *[P.CELLTYPE / name for name in ("panel.py", "run_step2.py", "run_step3.py", "alignment_block.py")]]
    source_hashes = {str(path): sha(path) for path in source_files}
    ledger = candidate / "resume_run_20261008.json"
    record = json.loads(ledger.read_text()) if ledger.exists() else dict(donors=donors, stages=[],
                source_scripts=source_hashes, started_at_UTC=now(), baseline_modified=False,
                accuracy_validated=False, scope="Nine ND diagnostic comparison; not an all25 frozen-reference release")
    if record["source_scripts"] != source_hashes:
        raise ValueError("Resume script/producer changed; inspect and start a new recorded resume")
    record.update(status="running", controller_pid=os.getpid(), through=args.through,
                  RSS_limit_GiB=args.rss_limit_gib, threads=6)
    save(ledger, record)
    for arm in ARMS:
        env = os.environ.copy()
        env.update(IO60_EXPRESSION_ROOT=str(candidate / arm / "expression"),
                   IO60_ALIGNED_ROOT=str(candidate / arm / "uniform"), IO60_SCALE_PATH=str(candidate / arm / "scale.parquet"),
                   IO60_TYPING_ROOT=str(candidate / arm / "typing"), IO60_GEOMETRY_ROOT=str(candidate / "geometry"),
                   IO60_SUBPANEL_LOG=str(P.GEOMETRY / "subpanel/subpanel_log.json"))
        path = candidate / arm / "typing"; da = ["--donors", *donors]
        stages = []
        if arm == ARMS[1]:
            stages.append(("uniform", [sys.executable, str(HERE / "reprocess_protein_own_compartment_candidate.py"), "--worker", "uniform", "--out", str(candidate), "--arm", arm, *da],
                           [candidate / arm / "uniform", candidate / arm / "scale.parquet", candidate / arm / "shifts.json"]))
        if args.through != "calibration":
            stages += [
                ("subpanel", [sys.executable, str(P.TYPING / "subpanel.py"), *da, "--membership", "geometry", "--emit-detail", "--refit", "--work", str(path / "scratch/subpanel_context"), "--audit-out", str(path / "subpanel_audit.json"), "--out", str(path / "subpanel_rebuilt")], [path / "subpanel_rebuilt", path / "subpanel_audit.json"]),
                ("build_cells_primary", [sys.executable, str(P.TYPING / "build_cells.py"), *da, "--calls", str(path / "subpanel_rebuilt"), "--out", str(path / "cells")], [path / "cells"])]
        # Refinement remains an explicit later stage so a pending producer contract
        # (e.g. CD8 representation policy) cannot be silently bypassed.
        if args.through == "immune":
            raise ValueError("Refinement continuation must be added only after the producer representation contract is resolved")
        for name, command, outputs in stages:
            key = arm + ":" + name
            prior = next((s for s in record["stages"] if s["key"] == key and s["returncode"] == 0), None)
            if prior:
                for output in prior["outputs"]:
                    if sha(output["path"]) != output["sha256"]:
                        raise ValueError(f"Completed resumed artifact changed: {output['path']}")
                print("VERIFIED COMPLETED STAGE", key, flush=True)
                continue
            stage = dict(key=key, started_at_UTC=now(), status="running")
            save(path / "run.json", dict(status="running", donors=donors, current_stage=key,
                controller_pid=os.getpid(), source_scripts=source_hashes, baseline_unchanged=True,
                accuracy_validated=False, scope=record["scope"]))
            record["current_stage"] = key; save(ledger, record)
            print("START", key, flush=True)
            try:
                result = run(command, env, path / (name + "_resume_20261008.log"), args.rss_limit_gib)
                stage.update(result, finished_at_UTC=now(), status="complete" if result["returncode"] == 0 else "failed")
                if result["returncode"]:
                    raise RuntimeError(f"{key} exited {result['returncode']}")
                stage["outputs"] = artifacts(outputs)
            except Exception as exc:
                stage.update(status="failed", error=str(exc), finished_at_UTC=now())
                record["stages"].append(stage); record.update(status="failed", current_stage=key)
                save(ledger, record); save(path / "run.json", dict(status="failed", donors=donors, stage=stage))
                raise
            record["stages"].append(stage); save(ledger, record)
            print("COMPLETE", key, "peak GiB", stage["peak_process_group_RSS_GiB"], flush=True)
        save(path / "run.json", dict(status="primary_complete" if args.through == "primary" else "calibration_complete",
            donors=donors, source_scripts=source_hashes, stages=[s for s in record["stages"] if s["key"].startswith(arm + ":")],
            baseline_unchanged=True, accuracy_validated=False, scope=record["scope"]))
    record.update(status=args.through + "_complete", finished_at_UTC=now(), current_stage=None)
    save(ledger, record)
    print("RESUME COMPLETE THROUGH", args.through, flush=True)

if __name__ == "__main__":
    main()
