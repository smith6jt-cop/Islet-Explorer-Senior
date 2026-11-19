#!/usr/bin/env python3
"""Apply the pinned port's existing exact parity tests to every declared CODEX image."""
import argparse
import csv
import json
import os
import subprocess
import sys
from pathlib import Path


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--snapshot", required=True, type=Path)
    parser.add_argument("--run", required=True, type=Path)
    parser.add_argument("--out", required=True, type=Path)
    args = parser.parse_args()
    args.out.mkdir(parents=True, exist_ok=False)
    with (args.run / "typing/levels.csv").open() as handle:
        tags = sorted({r["image"] for r in csv.DictReader(handle)})
    result = []
    for tag in tags:
        env = dict(os.environ, PYTHONPATH=str(args.snapshot.resolve() / "core/multiplex5000"),
                   MPLX5K_CODEX_PARITY_RUN=str(args.run), MPLX5K_CODEX_PARITY_IMAGE=tag,
                   OPENBLAS_NUM_THREADS="4", OMP_NUM_THREADS="4")
        with (args.out / f"{tag}.log").open("w") as log:
            run = subprocess.run([sys.executable, "-m", "pytest", "-q",
                str(args.snapshot / "core/multiplex5000/tests/current/test_positivity_parity.py"),
                f"--junitxml={args.out / (tag + '.xml')}"], env=env, stdout=log, stderr=subprocess.STDOUT, check=False)
        result.append({"image": tag, "returncode": run.returncode,
                       "status": "pass" if run.returncode == 0 else "failed_or_blocked",
                       "log": str(args.out / f"{tag}.log")})
        (args.out / "summary.json").write_text(json.dumps(result, indent=2) + "\n")
        print(tag, run.returncode, flush=True)


if __name__ == "__main__":
    main()
