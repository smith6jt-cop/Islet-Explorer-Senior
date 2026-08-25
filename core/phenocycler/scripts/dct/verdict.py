"""The release verdict, computed on the held-out donors and nowhere else.

Per region, pooled per donor, at the aggregation the protocol fixed in advance: a donor
passes only if the pooled statistic clears the bar **and** at least 6 of its 8 regions clear
it individually, so one large region cannot carry a donor.  Lane PASS requires every
evaluation donor.

`--selection` reruns the same thing on 6539 / 6450 / 6523 for comparison.  Those numbers are
NOT the verdict -- the configuration was chosen on them -- and are printed only so the gap
between fitted and held-out performance is visible rather than assumed to be zero.

    python scripts/dct/verdict.py                           # from core/phenocycler/
"""
from __future__ import annotations

import json
import pathlib
import sys

import numpy as np
import pandas as pd

HERE = pathlib.Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from arms import OUT                                        # noqa: E402
from concordance_gate import CEILING, GateInputs, run_gate   # noqa: E402

CLASSES = tuple(json.loads((HERE / "broad_class_crosswalk.json").read_text())["entries"])
EVALUATION = ["115", "6521", "6593", "6534", "6591"]
SELECTION = ["6539", "6450", "6523"]
BLOCKING = ["sign", "specificity", "detectability", "argmax_coherence",
            "impossible_mass", "null_separation"]
MIN_REGIONS = 6


def _inputs(frames, label):
    d = pd.concat(frames, ignore_index=True) if len(frames) > 1 else frames[0]
    P = d[list(CLASSES)].to_numpy(float)
    P = P / P.sum(1, keepdims=True)
    meas = {c: d[c].to_numpy(float) for c in d.columns
            if c not in CLASSES and pd.api.types.is_numeric_dtype(d[c])}
    return GateInputs(P, CLASSES, meas, label=label), len(d)


def judge(donors, log_name, checkpoint="2026-06-23-ptft", arm="A0-baseline"):
    recs = [json.loads(l) for l in (OUT / log_name).read_text().splitlines() if l.strip()]
    errs = [r for r in recs if "error" in r]
    if errs:
        print(f"!! {len(errs)} runs errored -- recorded, and each is a donor that cannot pass:")
        for e in errs:
            print(f"   {e['donor']}/r{e['region']}: {e['error'][:140]}")
    recs = [r for r in recs if "error" not in r and r["checkpoint"] == checkpoint and r["arm"] == arm]

    out = {}
    for donor in donors:
        rows = sorted([r for r in recs if r["donor"] == donor], key=lambda r: r["region"])
        if not rows:
            print(f"   {donor}: no runs -- a missing donor is a failed donor, not an absent one")
            out[donor] = {"passed": False, "reason": "no runs"}
            continue
        frames, per_region = [], {}
        for r in rows:
            f = OUT / donor / f"arm_{arm}_{checkpoint}_r{r['region']}.parquet"
            d = pd.read_parquet(f)
            frames.append(d)
            gi, n = _inputs([d], f"{donor}/r{r['region']}")
            g = run_gate(gi, restricted=True)
            per_region[r["region"]] = {k: g.criteria[k][0] for k in g.criteria} | {"n": n}
        gi, n = _inputs(frames, donor)
        pooled = run_gate(gi, restricted=True)
        n_clear = sum(1 for v in per_region.values() if all(v[k] for k in BLOCKING))
        passed = pooled.passed and n_clear >= MIN_REGIONS
        out[donor] = {"passed": passed, "pooled": pooled, "per_region": per_region,
                      "n_cells": n, "regions_clear": n_clear}
    return out


def report(res, title):
    print("\n" + "=" * 104)
    print(title)
    print(f"{'donor':<8}{'cells':>9}{'regions':>9}  " + "".join(f"{k[:11]:>13}" for k in BLOCKING) + "   verdict")
    for donor, v in res.items():
        if "pooled" not in v:
            print(f"{donor:<8}{'—':>9}{'—':>9}   {v['reason']}")
            continue
        flags = "".join(f"{('PASS' if v['pooled'].criteria[k][0] else 'FAIL'):>13}" for k in BLOCKING)
        print(f"{donor:<8}{v['n_cells']:>9,}{v['regions_clear']:>6}/8  {flags}"
              f"   {'PASS' if v['passed'] else 'FAIL'}")
    lane = bool(res) and all(v.get("passed") for v in res.values())
    print(f"\nLANE: {'PASS' if lane else 'FAIL'}  "
          f"({sum(1 for v in res.values() if v.get('passed'))}/{len(res)} donors)")
    return lane


def main():
    print(f"Release gate — restricted label space (26 of 51), 28 registry-derived anchors.")
    print(f"Permutation null and the {CEILING} plane-identity ceiling are printed per donor below.")
    ev = judge(EVALUATION, "evaluation_c2.jsonl")
    lane = report(ev, "HELD OUT — the verdict. Nothing was fitted to these donors.")

    print("\nPer-criterion detail, pooled per held-out donor:")
    for donor, v in ev.items():
        if "pooled" not in v:
            continue
        print(f"\n  {donor}  (null |rho|<={v['pooled'].null['rho_abs']:.3f}, "
              f"|AUROC-0.5|<={v['pooled'].null['auroc_dev']:.3f})")
        for k in BLOCKING + ["confidence_informativeness"]:
            ok, why = v["pooled"].criteria[k]
            print(f"    [{'PASS' if ok else 'FAIL'}] {k:<26} {why[:150]}")

    if "--selection" in sys.argv:
        sel = judge(SELECTION, "factorial_c1.jsonl")
        report(sel, "SELECTION donors — NOT the verdict; the configuration was chosen on these.")
    return 0 if lane else 1


if __name__ == "__main__":
    sys.exit(main())
