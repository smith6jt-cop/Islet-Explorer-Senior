#!/usr/bin/env python3
"""Build the pancreas protein typing generation the analysis lane reads, then pin it.

The chain is the October 8 one, run with the producers recovered byte for byte into this branch
(`recovered_v3_producers_20261018.json`) and the 2026-10-08 generation's aligned values, scale, shifts
and geometry. Two rulings change it:

  2026-10-18  Ker8_18 is a veto in no immune definition (`subpanel_io.UNSHARED_NEGATIVES`), so the
              four same-slide donors are typed like every donor.
  2026-10-19  CD8 carries UniFORM's shift on all 25 donors in the refinement as well as the primary
              calls. The 2026-10-08 chain divided it back out for the refinement only.

The primary typing (subpanel, primary build) does not read the refinement, so a generation that only
reruns the refinement reuses it from `--primary-from` as hard links. Steps, in the order the original
controllers ran them: refinement, refinement join, immune subtypes, cell release, analysis admission,
conservative candidate. The last step pins every donor's release, admission and candidate file by
SHA-256 in `<root>/contract.json`, with the row counts qualify must reproduce.

    python -m analysis.cohesive.typing_generation        # from the repository root
"""
from __future__ import annotations

import argparse
import importlib.util
import json
import os
import subprocess
import sys
import time
from concurrent.futures import ThreadPoolExecutor
from datetime import datetime, timezone
from pathlib import Path

import pyarrow.compute as pc
import pyarrow.parquet as pq

from analysis.cohesive.common import binding, require

GEN = Path("/data/io60-calibration/own_compartment_frozen_reference_repair_20261008")
V3 = Path("/home/smith6jt/IO60panc2nd/Islet-Explorer-Senior-v3")
BRANCH = Path(__file__).resolve().parents[2]
TYPING_SRC = BRANCH / "core/multiplex5000/scripts/typing"
PY = "/home/smith6jt/IO60panc2nd/Islet-Explorer-Senior/.venv/bin/python"
ROOT = Path("/data/io60-analysis/typing-20261019/run")
PRIMARY_FROM = Path("/data/io60-analysis/typing-20261018/run/typing")
DONORS = ["0112", "115", "6347", "6374", "6380", "6436", "6442", "6450", "6461", "6476", "6482", "6493", "6505",
          "6521", "6523", "6533", "6534", "6539", "6551", "6566", "6579", "6580", "6591", "6593", "6623"]
RULINGS = ("2026-10-18: Ker8_18 is a veto in no immune definition (shared epithelial negatives, all 25)",
           "2026-10-19: CD8 carries UniFORM's shift on all 25 donors, in the refinement as in the primary calls")
CANDIDATE_BUILDER = (V3 / "results/v3/islet_immune_reassessment_20261006/cell_typing_validation/"
                     "peer_review_typing_20261008/protein/build_conservative_protein_candidate.py")
#: What the primary typing wrote; reused unchanged when only the refinement reruns.
PRIMARY = ("subpanel_rebuilt", "subpanel_audit.json")
PRIMARY_CELLS = "cells_primary_original"


def environment(typing: Path) -> dict:
    return {"IO60_EXPRESSION_ROOT": str(GEN / "own_compartment/expression"),
            "IO60_ALIGNED_ROOT": str(GEN / "own_compartment/uniform"),
            "IO60_SCALE_PATH": str(GEN / "own_compartment/scale.parquet"), "IO60_GEOMETRY_ROOT": str(GEN / "geometry"),
            "IO60_SUBPANEL_LOG": str(V3 / "data/development/typing/subpanel/subpanel_log.json"),
            "IO60_SHIFTS_PATH": str(GEN / "own_compartment/shifts.json"), "IO60_TYPING_ROOT": str(typing),
            "OMP_NUM_THREADS": "6", "OPENBLAS_NUM_THREADS": "6", "MKL_NUM_THREADS": "6", "NUMBA_NUM_THREADS": "6"}


class Generation:
    def __init__(self, root: Path):
        self.root, self.typing = root, root / "typing"
        self.env = environment(self.typing)
        self.ledger = {"created_UTC": datetime.now(timezone.utc).isoformat(), "ruling": "; ".join(RULINGS),
                       "producers": str(BRANCH / "analysis/cohesive/recovered_v3_producers_20261018.json"),
                       "base_generation": str(GEN), "environment": self.env, "steps": []}

    def save(self) -> None:
        (self.root / "generation.json").write_text(json.dumps(self.ledger, indent=2) + "\n")

    def run(self, name: str, command: list[str]) -> None:
        started = time.time()
        log = self.root / f"{name}.log"
        print(f"START {name}", flush=True)
        with log.open("w") as handle:
            result = subprocess.run(command, env={**os.environ, **self.env}, stdout=handle,
                                    stderr=subprocess.STDOUT, cwd="/tmp")
        self.ledger["steps"].append({"name": name, "command": command, "returncode": result.returncode,
                                     "seconds": round(time.time() - started, 1), "log": str(log)})
        self.save()
        require(result.returncode == 0, f"{name} failed ({result.returncode}); see {log}")
        print(f"DONE  {name} in {time.time() - started:.0f}s", flush=True)

    def reuse_primary(self, source: Path) -> None:
        """Hard links, so the reused files stay whole if the earlier generation is deleted."""
        for name in (*PRIMARY, PRIMARY_CELLS):
            require((source / name).exists(), f"Primary typing missing: {source / name}")
        self.typing.mkdir(parents=True)
        for name in PRIMARY:
            subprocess.run(["cp", "-al", str(source / name), str(self.typing / name)], check=True)
        subprocess.run(["cp", "-al", str(source / PRIMARY_CELLS), str(self.typing / "cells")], check=True)
        self.ledger["steps"].append({"name": "reuse_primary", "source": str(source),
                                     "linked": [*PRIMARY, f"{PRIMARY_CELLS} -> cells"]})
        self.save()

    def candidate(self) -> None:
        """The original conservative builder, executed from its own file, with only its two folders moved."""
        spec = importlib.util.spec_from_file_location("conservative_candidate", CANDIDATE_BUILDER)
        module = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(module)
        module.SOURCE = self.root / "pancreas_analysis_assignments_current_generation"
        module.OUT = self.root / "conservative_assignment_candidate"
        started = time.time()
        module.main()
        self.ledger["steps"].append({"name": "conservative_candidate", "builder": str(CANDIDATE_BUILDER),
                                     "builder_sha256": binding(CANDIDATE_BUILDER)["sha256"],
                                     "seconds": round(time.time() - started, 1)})
        self.save()


def counts(path: Path) -> dict:
    table = pq.read_table(path, columns=["qc_analysis_eligible", "analysis_eligible"])
    return {"expected_source_rows": table.num_rows,
            "expected_QC_rows": int(pc.sum(pc.fill_null(table["qc_analysis_eligible"], False)).as_py()),
            "expected_additional_admission_rows": int(pc.sum(pc.fill_null(table["analysis_eligible"], False)).as_py())}


def donor_record(root: Path, donor: str) -> dict:
    padded = donor.zfill(4) if donor == "112" else donor
    release = root / "pancreas_cell_release_current_generation" / f"donor_id={padded}" / "data.parquet"
    admission = root / "pancreas_analysis_assignments_current_generation" / f"donor_id={donor}" / "data.parquet"
    candidate = root / "conservative_assignment_candidate" / f"donor_id={donor}" / "data.parquet"
    for path in (release, admission, candidate):
        require(path.is_file(), f"Missing generation file: {path}")
    record = {"release": binding(release), "admission": binding(admission), "candidate": binding(candidate),
              **counts(candidate)}
    require(record == {**record, **counts(admission)}, f"Candidate and admission counts differ: {donor}")
    return record


def pin(root: Path) -> Path:
    """Write `<root>/contract.json`: what qualify reads, each file by SHA-256."""
    generation = json.loads((root / "generation.json").read_text())
    require("finished_UTC" in generation, "The generation did not finish")
    donors = sorted(p.name.split("=", 1)[1] for p in (root / "pancreas_analysis_assignments_current_generation").glob("donor_id=*"))
    require(len(donors) == len(DONORS), f"Expected {len(DONORS)} donors, found {len(donors)}")
    with ThreadPoolExecutor(16) as pool:
        records = dict(zip(donors, pool.map(lambda d: donor_record(root, d), donors)))
    out = root / "contract.json"
    out.write_text(json.dumps({
        "generation_root": str(root), "generation_manifest": binding(root / "generation.json"),
        "ruling": generation["ruling"], "base_generation": generation["base_generation"],
        "donors": records}, indent=2, sort_keys=True) + "\n")
    return out


def main(argv=None) -> None:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--root", type=Path, default=ROOT, help="new generation folder; must not exist")
    ap.add_argument("--primary-from", type=Path, default=PRIMARY_FROM,
                    help="typing folder of the generation whose primary typing is reused")
    args = ap.parse_args(argv)
    require(not args.root.exists(), f"{args.root} exists; a generation is written once")
    args.root.mkdir(parents=True)
    g = Generation(args.root)
    T, donors = g.typing, ["--donors", *DONORS]
    g.reuse_primary(args.primary_from)
    g.run("refine", [PY, str(TYPING_SRC / "refine.py"), *donors, "--out", str(T / "refinement"),
                     "--cd8-alignment", "aligned", "--shifts", g.env["IO60_SHIFTS_PATH"]])
    g.run("build_cells_with_refinement", [PY, str(TYPING_SRC / "build_cells.py"), *donors,
          "--calls", str(T / "subpanel_rebuilt"), "--refine", str(T / "refinement"),
          "--audit", str(T / "subpanel_audit.json"), "--out", str(T / "cells_with_refinement")])
    (T / "cells").rename(T / PRIMARY_CELLS)
    (T / "cells").symlink_to("cells_with_refinement")
    g.run("immune_subtype", [PY, str(TYPING_SRC / "immune_subtype.py"), *donors])
    original = json.loads((GEN / "own_compartment/typing/run.json").read_text())
    (T / "run.json").write_text(json.dumps({
        "status": "complete", "donors": DONORS, "geometry_root": str(GEN / "geometry"),
        "original_geometry_root": str(V3 / "data/development/typing"),
        "aligned_root": str(GEN / "own_compartment/uniform"), "scale": str(GEN / "own_compartment/scale.parquet"),
        "expression_root": str(GEN / "own_compartment/expression"), "source_scripts": original["source_scripts"],
        "producers_recovered_into": str(BRANCH), "ruling": g.ledger["ruling"], "refinement_cd8_alignment": "aligned",
        "baseline_unchanged": True, "accuracy_validated": False, "release_promoted": False}, indent=2) + "\n")
    g.run("cell_release", [PY, str(TYPING_SRC / "build_pancreas_cell_release.py"), "--arm-root", str(T),
          "--out", str(args.root / "pancreas_cell_release_current_generation"),
          "--input-manifest", str(V3 / "results/v3/release_20261003/current_input_manifest.json")])
    g.run("analysis_admission", [PY, str(BRANCH / "analysis/xenium/scripts/materialize_protein_analysis_assignments.py"),
          "--source", str(args.root / "pancreas_cell_release_current_generation"),
          "--out", str(args.root / "pancreas_analysis_assignments_current_generation"), "--policy", "current_generation",
          "--mandatory-exclusions", "/data/io60-calibration/pancreas_analysis_assignments_single_exclusions_20261005"])
    g.candidate()
    g.ledger["finished_UTC"] = datetime.now(timezone.utc).isoformat()
    g.save()
    print(f"GENERATION COMPLETE; pinned in {pin(args.root)}", flush=True)


if __name__ == "__main__":
    sys.exit(main())
