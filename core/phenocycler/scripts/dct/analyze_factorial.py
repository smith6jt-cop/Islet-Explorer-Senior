"""C1 readout: which confound, if any, moves the gate — and on how many donors.

A pass/fail flip is too coarse to screen with: every arm fails most criteria, so counting
flips would report "nothing moved" while an arm was carrying half the anchors across.  So
this rescores each arm's saved parquet into CONTINUOUS statistics and compares each arm to
`A0-baseline` on the same donor, region and checkpoint.

The stopping rule from the protocol: an arm counts as moving an axis only if it moves it
**on at least 2 of the 3 screening donors, in the same direction**.  A single-donor effect
on a 3-donor screen is noise until shown otherwise -- and 6539, which is in the screen, is
a cohort outlier on the one quantity this study measures.

    python scripts/dct/analyze_factorial.py                  # from core/phenocycler/
"""
from __future__ import annotations

import json
import pathlib
import sys
from collections import defaultdict

import numpy as np
import pandas as pd

HERE = pathlib.Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from arms import ARMS, CHECKPOINTS, OUT                        # noqa: E402
from concordance_gate import CRITERIA, GateInputs, load_anchors, run_gate  # noqa: E402

CLASSES = tuple(json.loads((HERE / "broad_class_crosswalk.json").read_text())["entries"])
ANCHORS = load_anchors()
BASELINE = "A0-baseline"

#: Higher is better for every one of these except impossible_mass.
STATS = ["sign_positive", "rank1", "auroc_pass", "auroc_mean", "coherent", "impossible_mass"]
LOWER_IS_BETTER = {"impossible_mass"}


def score_one(parquet):
    d = pd.read_parquet(parquet)
    P = d[list(CLASSES)].to_numpy(float)
    P = P / P.sum(1, keepdims=True)
    meas = {c: d[c].to_numpy(float) for c in d.columns
            if c not in CLASSES and pd.api.types.is_numeric_dtype(d[c])}
    r = run_gate(GateInputs(P, CLASSES, meas), restricted=True)
    au = [r.auroc.get(m, np.nan) for m, _ in ANCHORS]
    bars = [CRITERIA["auroc_immune"] if c == "Immune" else CRITERIA["auroc_other"]
            for _, c in ANCHORS]
    coh = [v for v in r.coherence.values() if all(np.isfinite(x) for x in v)]
    return {
        "n_cells": len(d),
        "sign_positive": sum(1 for m, c in ANCHORS
                             if np.isfinite(r.rho.get(m, {}).get(c, np.nan))
                             and r.rho[m][c] > 0),
        "rank1": sum(1 for m, _ in ANCHORS if r.rank.get(m, 99) == 1),
        "auroc_pass": sum(1 for a, b in zip(au, bars) if np.isfinite(a) and a >= b),
        "auroc_mean": float(np.nanmean(au)) if np.isfinite(au).any() else np.nan,
        "coherent": sum(1 for t, b in coh if t > b),
        "impossible_mass": r.impossible_mass_mean,
    }


def main():
    log = OUT / "factorial_c1.jsonl"
    recs = [json.loads(l) for l in log.read_text().splitlines() if l.strip()]
    errs = [r for r in recs if "error" in r]
    recs = [r for r in recs if "error" not in r]
    if errs:
        print(f"!! {len(errs)} arm/region combinations errored — recorded, not hidden:")
        for e in errs[:10]:
            print(f"   {e['donor']}/{e['arm']}/{e['checkpoint']}/r{e['region']}: {e['error'][:140]}")
        print()

    S = {}
    for r in recs:
        f = OUT / r["donor"] / f"arm_{r['arm']}_{r['checkpoint']}_r{r['region']}.parquet"
        if not f.exists():
            print(f"   missing parquet for {r['donor']}/{r['arm']}/{r['checkpoint']}/r{r['region']}")
            continue
        S[(r["donor"], r["arm"], r["checkpoint"], r["region"])] = score_one(f)

    donors = sorted({k[0] for k in S})
    ckpts = [c for c in CHECKPOINTS if any(k[2] == c for k in S)]
    print(f"{len(S)} scored runs over {len(donors)} donors, {len(ckpts)} checkpoints\n")

    # ---- baseline table -------------------------------------------------
    print("=" * 100)
    print(f"BASELINE ({BASELINE}), mean over regions — 28 anchors, restricted label space")
    print(f"{'donor':<8}{'checkpoint':<20}{'sign':>6}{'rank1':>7}{'AUROC>=bar':>12}"
          f"{'mean AUROC':>12}{'coherent':>10}{'impossible':>12}")
    for d in donors:
        for c in ckpts:
            v = [S[k] for k in S if k[0] == d and k[1] == BASELINE and k[2] == c]
            if not v:
                continue
            m = {s: np.mean([x[s] for x in v]) for s in STATS}
            print(f"{d:<8}{c:<20}{m['sign_positive']:>5.1f}/28{m['rank1']:>6.1f}/28"
                  f"{m['auroc_pass']:>11.1f}/28{m['auroc_mean']:>12.3f}"
                  f"{m['coherent']:>9.1f}{m['impossible_mass']:>12.3f}")

    # ---- arm deltas against the baseline --------------------------------
    for c in ckpts:
        print("\n" + "=" * 100)
        print(f"ARM DELTAS vs {BASELINE} — checkpoint {c}")
        print(f"{'arm':<22}{'donors moved':>14}  " + "".join(f"{s:>16}" for s in STATS))
        for arm in ARMS:
            if arm.name == BASELINE:
                continue
            per_donor = {}
            for d in donors:
                deltas = defaultdict(list)
                for k in S:
                    if k[0] != d or k[1] != arm.name or k[2] != c:
                        continue
                    b = S.get((d, BASELINE, c, k[3]))
                    if b is None:
                        continue
                    for s in STATS:
                        deltas[s].append(S[k][s] - b[s])
                if deltas:
                    per_donor[d] = {s: float(np.mean(v)) for s, v in deltas.items()}
            if not per_donor:
                continue
            avg = {s: float(np.mean([p[s] for p in per_donor.values()])) for s in STATS}
            # "Moved" is tied to the pre-registered bar, not to a made-up score: the arm
            # either carries anchors across their own thresholds, or brings the
            # impossible-class mass under the criterion-5 bar it was failing.
            def improved(p, base):
                crossed_mass = (base["impossible_mass"] > CRITERIA["impossible_mass_mean"]
                                and base["impossible_mass"] + p["impossible_mass"]
                                <= CRITERIA["impossible_mass_mean"])
                anchors_gained = sum(p[s] for s in ("sign_positive", "rank1", "auroc_pass", "coherent"))
                return crossed_mass or anchors_gained >= 4      # >=1 anchor on each count statistic
            base_by_donor = {d: {s: float(np.mean([S[k][s] for k in S
                                                   if k[0] == d and k[1] == BASELINE and k[2] == c]))
                                 for s in STATS} for d in per_donor}
            n_moved = sum(1 for d, p in per_donor.items() if improved(p, base_by_donor[d]))
            flag = "**" if n_moved >= 2 else "  "
            print(f"{flag}{arm.name:<20}{n_moved:>8}/{len(per_donor):<5}  "
                  + "".join(f"{avg[s]:>+16.3f}" for s in STATS))

    print("\n** = moved the gate on >= 2 of the screening donors in the same direction.")
    print("   Nothing flagged means domain/panel incompatibility is the surviving explanation,")
    print("   C2 is not run, and the verdict is recorded with the effect sizes that were not reached.")


if __name__ == "__main__":
    main()
