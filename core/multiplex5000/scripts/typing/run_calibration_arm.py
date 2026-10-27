#!/usr/bin/env python3
"""Run one isolated arm of the REDSEA/raw × CD8 × scale typing comparison.

Each arm gets independent calls, refinement, cell tables, immune subtyping, and a resumable
stage log. The registry classification criteria are read from their pinned source file, while
alignment, scale, and geometry are declared separately in ``run.json``.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import os
import subprocess
import sys
from datetime import datetime, timezone
from pathlib import Path


HERE = Path(__file__).resolve().parent
REPO = HERE.parents[3]
SCRIPTS = HERE
DEFAULT_GEOMETRY = REPO / "data/development/typing"
DEFAULT_REGISTRY_LOG = DEFAULT_GEOMETRY / "subpanel/subpanel_log.json"


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for block in iter(lambda: f.read(8 * 1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def atomic_json(path: Path, value: dict) -> None:
    temp = path.with_suffix(path.suffix + ".tmp")
    temp.write_text(json.dumps(value, indent=2, sort_keys=True) + "\n")
    temp.replace(path)


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--aligned", required=True, type=Path)
    ap.add_argument("--scale", required=True, type=Path)
    ap.add_argument("--out", required=True, type=Path)
    ap.add_argument("--geometry-root", type=Path, default=DEFAULT_GEOMETRY)
    ap.add_argument("--registry-log", type=Path, default=DEFAULT_REGISTRY_LOG)
    ap.add_argument("--refinement-cd8-alignment", choices=["aligned", "unshifted"],
                    default="aligned",
                    help="CD8 handling for refinement only; primary calls use --aligned")
    ap.add_argument("--shifts", type=Path,
                    help="matching UniFORM shifts for the declared aligned root")
    ap.add_argument("--donors", nargs="*", default=None)
    ap.add_argument("--representation-donors", nargs="+", default=None,
                    help="fit refinement on these recorded donors; default is the full cohort")
    ap.add_argument("--threads", type=int, default=6,
                    help="CPU threads for numerical libraries in this arm")
    ap.add_argument("--force", action="store_true")
    args = ap.parse_args()

    aligned = args.aligned.resolve()
    scale = args.scale.resolve()
    geometry = args.geometry_root.resolve()
    registry = args.registry_log.resolve()
    out = args.out.resolve()
    shifts = args.shifts.resolve() if args.shifts else None
    if not aligned.is_dir() or not scale.is_file() or not registry.is_file():
        raise FileNotFoundError("aligned root, scale table, and registry log must all exist")
    if args.refinement_cd8_alignment == "unshifted" and (shifts is None or not shifts.is_file()):
        raise FileNotFoundError("unshifted refinement CD8 requires --shifts for this aligned input")
    donors_on_disk = sorted(p.name.split("=", 1)[1] for p in aligned.glob("donor_id=*")
                            if (p / "data.parquet").is_file())
    donors = list(args.donors) if args.donors else donors_on_disk
    if not donors or len(donors) != len(set(donors)):
        raise ValueError("donor list must be nonempty and unique")
    missing = sorted(set(donors) - set(donors_on_disk))
    if missing or "6457" in donors:
        raise ValueError(f"missing aligned donors={missing}; 6457 is explicitly excluded")
    representation_donors = (list(args.representation_donors)
                             if args.representation_donors else donors)
    if (len(representation_donors) != len(set(representation_donors))
            or not set(representation_donors) <= set(donors)):
        raise ValueError("representation donors must be unique members of the run's donor scope")
    missing_geometry = sorted(d for d in donors
                              if not (geometry / f"revised_cells/donor_id={d}/data.parquet").is_file()
                              and not (geometry / f"expression/donor_id={d}/data.parquet").is_file()
                              and not (REPO / f"core/phenocycler/data/runs/30f6689e4d27d295f890/"
                                       f"expression/donor_id={d}/data.parquet").is_file())
    if missing_geometry:
        raise FileNotFoundError(f"geometry unavailable for {missing_geometry}")

    out.mkdir(parents=True, exist_ok=True)
    existing_cell_donors = {p.name.split("=", 1)[1] for p in (out / "cells").glob("donor_id=*")
                            if (p / "data.parquet").is_file()}
    if existing_cell_donors - set(donors):
        raise ValueError("output contains cell partitions outside the declared donor cohort; "
                         "use a fresh output directory to avoid stale representation inputs")
    manifest_path = out / "run.json"
    if manifest_path.exists() and not args.force:
        raise FileExistsError(f"{manifest_path} exists; use --force to replace the run")
    script_paths = [Path(__file__), SCRIPTS / "subpanel.py", SCRIPTS / "build_cells.py",
                    SCRIPTS / "refine.py", SCRIPTS / "immune_subtype.py",
                    REPO / "core/multiplex5000/multiplex5000/tissues/t_cell.py"]
    run = {
        "schema_version": 1,
        "status": "running",
        "started_utc": datetime.now(timezone.utc).isoformat(),
        "donors": donors,
        "representation_donors": representation_donors,
        "aligned_root": str(aligned),
        "aligned_files": {d: sha256(aligned / f"donor_id={d}/data.parquet") for d in donors},
        "scale": str(scale),
        "scale_sha256": sha256(scale),
        "refinement_cd8_alignment": args.refinement_cd8_alignment,
        "shifts_path": str(shifts) if shifts else None,
        "shifts_sha256": sha256(shifts) if shifts else None,
        "geometry_root": str(geometry),
        "registry_criteria": str(registry),
        "registry_criteria_sha256": sha256(registry),
        "cpu_threads": args.threads,
        "source_scripts": {str(p): sha256(p) for p in script_paths},
        "stages": [],
        "typing_root": str(out),
    }
    atomic_json(manifest_path, run)

    env = os.environ.copy()
    env.update({
        "IO60_ALIGNED_ROOT": str(aligned),
        "IO60_SCALE_PATH": str(scale),
        "IO60_TYPING_ROOT": str(out),
        "IO60_GEOMETRY_ROOT": str(geometry),
        "IO60_SUBPANEL_LOG": str(registry),
        "NUMBA_NUM_THREADS": str(args.threads),
        "OMP_NUM_THREADS": str(args.threads),
        "OPENBLAS_NUM_THREADS": str(args.threads),
        "MKL_NUM_THREADS": str(args.threads),
    })
    if shifts:
        env["IO60_SHIFTS_PATH"] = str(shifts)
    donor_args = ["--donors", *donors]
    stages = [
        ("subpanel", [sys.executable, str(SCRIPTS / "subpanel.py"), *donor_args,
                       "--membership", "geometry", "--emit-detail", "--refit",
                       "--work", str(out / "scratch/subpanel_context"),
                       "--audit-out", str(out / "subpanel_audit.json"),
                       "--out", str(out / "subpanel_rebuilt")]),
        ("build_cells_primary", [sys.executable, str(SCRIPTS / "build_cells.py"), *donor_args,
                                  "--calls", str(out / "subpanel_rebuilt"),
                                  "--out", str(out / "cells")]),
        ("refine", [sys.executable, str(SCRIPTS / "refine.py"), *donor_args,
                    "--out", str(out / "refinement"),
                    "--cd8-alignment", args.refinement_cd8_alignment,
                    *(["--representation-donors", *args.representation_donors]
                      if args.representation_donors else []),
                    *(["--shifts", str(shifts)] if shifts else [])]),
        ("build_cells_with_refinement", [sys.executable, str(SCRIPTS / "build_cells.py"),
                                          *donor_args,
                                          "--calls", str(out / "subpanel_rebuilt"),
                                          "--refine", str(out / "refinement"),
                                          "--audit", str(out / "subpanel_audit.json"),
                                          "--out", str(out / "cells")]),
        ("immune_subtype", [sys.executable, str(SCRIPTS / "immune_subtype.py"), *donor_args]),
    ]
    try:
        for name, command in stages:
            log = out / f"{name}.log"
            with log.open("w") as fh:
                proc = subprocess.run(command, cwd=REPO, env=env, stdout=fh,
                                      stderr=subprocess.STDOUT, text=True)
            record = {"name": name, "command": command, "log": str(log),
                      "returncode": proc.returncode,
                      "finished_utc": datetime.now(timezone.utc).isoformat()}
            run["stages"].append(record)
            atomic_json(manifest_path, run)
            if proc.returncode:
                run["status"] = "failed"
                run["failed_stage"] = name
                atomic_json(manifest_path, run)
                print(f"FAILED {name}; see {log}", file=sys.stderr)
                return proc.returncode
            print(f"OK {name} -> {log}", flush=True)
        run["status"] = "complete"
        run["finished_utc"] = datetime.now(timezone.utc).isoformat()
        atomic_json(manifest_path, run)
        print(f"complete: {out}")
        return 0
    except BaseException:
        run["status"] = "interrupted"
        run["interrupted_utc"] = datetime.now(timezone.utc).isoformat()
        atomic_json(manifest_path, run)
        raise


if __name__ == "__main__":
    raise SystemExit(main())
