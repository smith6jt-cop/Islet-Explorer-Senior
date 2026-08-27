#!/usr/bin/env python
"""Decide whether the run is a result, and say exactly what is missing when it is not.

`sbatch` exiting 0 means a shell ran. It does not mean 20 donors were embedded, or that the
shards of a donor-arm tiled 1..M, or that a merge read the shards it thinks it read. Array
jobs are especially good at hiding this: one task of 120 can fail, the array reports as
completed, and the merge for that donor-arm will happily rbind two shards out of three --
which is why `kodama_merge.R` refuses a partial tile, and why this checks the same thing
independently of the job that wrote it.

Nine checks, all fail-closed:
  1. every cohort donor in qupath_manifest.json appears in the work manifest
  2. every donor has prep.json marked complete, with the recorded cell count
  3. the PCA input files exist and their row counts agree with prep.json
  4. every (donor, arm) has all n_shards present, with seeds tiling 1..M exactly and no
     shard written against a different input digest
  5. every (donor, arm) has an embedding whose cell_ids match the prep cell_ids exactly
  6. every embedding is finite -- a NaN coordinate would cluster silently
  7. every (donor, arm) has group calls and a composition for every k
  8. composition fractions sum to 1 per (donor, arm, k) and no lineage was awarded to a
     group the calls file says did not clear
  9. the KODAMA build recorded in every shard matches kodama.lock

Exit status is the answer: 0 = complete and internally consistent, 2 = not.
"""

from __future__ import annotations

import argparse
import json
import sys
from collections import defaultdict
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent


class Report:
    def __init__(self) -> None:
        self.problems: list[str] = []
        self.notes: list[str] = []

    def fail(self, msg: str) -> None:
        self.problems.append(msg)

    def note(self, msg: str) -> None:
        self.notes.append(msg)

    def ok(self) -> bool:
        return not self.problems


def _json(p: Path):
    try:
        return json.loads(p.read_text())
    except Exception:                                  # noqa: BLE001
        return None


def verify(run: Path, repo: Path, r: Report) -> dict:
    man = _json(run / "work_manifest.json")
    if man is None:
        r.fail(f"no readable work_manifest.json under {run}")
        return {}

    cohort = set(json.loads((repo / "data" / "qupath_manifest.json").read_text())
                 ["expected_donors"])
    listed = {d["donor_id"] for d in man["donors"]}
    if cohort - listed:
        r.fail(f"(1) cohort donors absent from the manifest: {sorted(cohort - listed)}")
    if listed - cohort:
        r.fail(f"(1) manifest lists donors outside the cohort: {sorted(listed - cohort)}")

    lock = _json(HERE / "kodama.lock") or {}
    M = int(man["kodama"]["M"])
    n_shards = int(man["kodama"]["n_shards"])
    arms = man["arms"]
    ks_expected: set[int] | None = None
    summary = defaultdict(dict)

    for d in man["donors"]:
        donor, roi = d["donor_id"], d["roi"]
        dest = run / f"donor_id={donor}" / f"roi={roi}"
        prep = _json(dest / "input" / "prep.json")
        if not prep or not prep.get("complete"):
            r.fail(f"(2) {donor}: prep.json missing or not complete")
            continue
        n = int(prep["n_cells"])
        summary[donor]["n_cells"] = n
        if prep.get("subsampled"):
            r.fail(f"(2) {donor}: prep.json records a SUBSAMPLE ({prep.get('n_used')} of "
                   f"{prep.get('n_available')}); this run is defined on all cells")

        ids_p = dest / "input" / "cell_id.csv"
        if not ids_p.is_file():
            r.fail(f"(3) {donor}: input/cell_id.csv is absent")
            continue
        ids = pd.read_csv(ids_p)["cell_id"].astype(str)
        if len(ids) != n:
            r.fail(f"(3) {donor}: cell_id.csv has {len(ids):,} rows, prep.json says {n:,}")
        for f in ("pcs.mtx", "coords.csv"):
            if not (dest / "input" / f).is_file():
                r.fail(f"(3) {donor}: input/{f} is absent")

        for arm in arms:
            key = f"{donor}/{arm}"
            # --- 4. shards tile 1..M ---
            shards = sorted((dest).glob(f"res_{arm}_shard*.json"))
            metas = [_json(p) for p in shards]
            metas = [m for m in metas if m and m.get("complete")]
            if len(metas) != n_shards:
                r.fail(f"(4) {key}: {len(metas)} complete shards, manifest says {n_shards}")
            ks = sorted({k for m in metas
                         for k in range(m["k_global_first"], m["k_global_last"] + 1)})
            if metas and ks != list(range(1, M + 1)):
                r.fail(f"(4) {key}: shard seeds cover {len(ks)} iterations, not 1..{M}")
            digests = {m.get("input_digest") for m in metas}
            if len(digests) > 1:
                r.fail(f"(4) {key}: shards were computed against {len(digests)} different "
                       f"input digests; they are not shards of one run")
            shas = {m.get("kodama_version"): m.get("kodama_version") for m in metas}
            if lock.get("kodama_version") and set(shas) - {lock["kodama_version"]}:
                r.fail(f"(9) {key}: shards report KODAMA {sorted(shas)}, "
                       f"lock says {lock['kodama_version']}")

            # --- 5/6. embedding ---
            emb_p = dest / f"embedding_{arm}.csv"
            kj = _json(dest / f"kodama_{arm}.json")
            if not emb_p.is_file() or not (kj and kj.get("complete")):
                r.fail(f"(5) {key}: no completed embedding")
                continue
            emb = pd.read_csv(emb_p)
            if len(emb) != n:
                r.fail(f"(5) {key}: embedding has {len(emb):,} rows, prep says {n:,}")
            elif not emb["cell_id"].astype(str).equals(ids):
                r.fail(f"(5) {key}: embedding cell_ids are not the prep cell_ids in order")
            bad = int((~np.isfinite(emb[["dim1", "dim2"]].to_numpy())).sum())
            if bad:
                r.fail(f"(6) {key}: {bad} non-finite embedding coordinates")
            if kj.get("n_shards") != n_shards:
                r.fail(f"(4) {key}: kodama_{arm}.json merged {kj.get('n_shards')} shards, "
                       f"manifest says {n_shards}")

            # --- 7/8. calls and composition ---
            calls_p, comp_p = dest / f"group_calls_{arm}.csv", dest / f"composition_{arm}.csv"
            if not (calls_p.is_file() and comp_p.is_file()):
                r.fail(f"(7) {key}: group_calls/composition missing")
                continue
            calls, comp = pd.read_csv(calls_p), pd.read_csv(comp_p)
            ks_here = set(comp["k"].unique())
            ks_expected = ks_here if ks_expected is None else ks_expected
            if ks_here != ks_expected:
                r.fail(f"(7) {key}: k sweep {sorted(ks_here)} != {sorted(ks_expected)}")
            for k, g in comp.groupby("k"):
                s = float(g["frac"].sum())
                if abs(s - 1.0) > 1e-6:
                    r.fail(f"(8) {key} k={k}: composition fractions sum to {s:.6f}")
            awarded = set(calls.loc[calls.lineage.notna() & (calls.lineage != ""), "lineage"])
            claimed = set(comp.loc[comp.lineage != "Unavailable", "lineage"])
            if claimed - awarded:
                r.fail(f"(8) {key}: composition claims {sorted(claimed - awarded)} which no "
                       f"group was awarded")
            unavail = comp.loc[comp.lineage == "Unavailable"].groupby("k")["frac"].sum()
            summary[donor][arm] = {
                "unavailable_by_k": {int(k): round(float(v), 4) for k, v in unavail.items()},
                "kodama_seconds": (kj.get("dissimilarity_seconds", 0) or 0)
                                  + (kj.get("visualization_seconds", 0) or 0),
            }
            r.note(f"{key}: Unavailable " +
                   ", ".join(f"k={int(k)} {v:.1%}" for k, v in unavail.items()))
    return dict(summary)


def assemble(run: Path) -> None:
    """One cohort composition table and one arm-difference table. No verdict, just the numbers."""
    parts = []
    for p in sorted(run.glob("donor_id=*/roi=*/composition_*.csv")):
        parts.append(pd.read_csv(p))
    if not parts:
        return
    allc = pd.concat(parts, ignore_index=True)
    allc.to_csv(run / "cohort_composition.csv", index=False)

    wide = allc.pivot_table(index=["donor_id", "k", "lineage"], columns="arm",
                            values="frac").reset_index()
    if {"spatial", "nonspatial"} <= set(wide.columns):
        # The arm difference is the measurement of how much location leaks into identity.
        # It is reported per donor and per k, never averaged into a single reassuring number.
        wide["spatial_minus_nonspatial"] = wide["spatial"] - wide["nonspatial"]
    wide.to_csv(run / "cohort_arm_difference.csv", index=False)
    print(f"[assemble] wrote cohort_composition.csv ({len(allc):,} rows) and "
          f"cohort_arm_difference.csv")


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--run", type=Path, required=True)
    ap.add_argument("--repo", type=Path, default=HERE.parents[2])
    ap.add_argument("--assemble", action="store_true")
    a = ap.parse_args(argv)

    r = Report()
    summary = verify(a.run, a.repo, r)
    if a.assemble:
        assemble(a.run)

    (a.run / "verification.json").write_text(json.dumps(
        {"complete": r.ok(), "problems": r.problems, "donors": summary}, indent=2))

    for n in r.notes:
        print(f"  {n}")
    print()
    if r.ok():
        print(f"RUN VERIFIED — {len(summary)} donors, every check passed.")
        return 0
    print(f"RUN INCOMPLETE — {len(r.problems)} problem(s):")
    for p in r.problems:
        print(f"  {p}")
    print("\nNothing here is a result until this exits 0.")
    return 2


if __name__ == "__main__":
    raise SystemExit(main())
