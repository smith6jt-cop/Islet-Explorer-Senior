"""The Xenium producer verdict: one table over the arms' gate records, the protocol's mechanical rule,
and the producer statement.  Nothing here consults the labels themselves.

Rule (protocol 2.6 = 1.6 with the cohort tie-break): a lineage is VERIFIED on a donor iff its core was
awarded on the arm's own cells, the label is spatially coherent against the within-section permutation
null, and the verdict is invariant across the arm's replicates (unjudged stability is never a pass).
An arm SURVIVES iff n_verified >= ceil(survive_frac x n_lineages) on BOTH donors.  Survivors rank by
(i) verified lineages summed over the donors, (ii) mean award margin (real - null p95) over the verified
lineages, (iii) lineages awarded on every cohort donor the arm ran on (0 when it ran on the pair only),
(iv) fewest external pins, then lowest primary wall time.  No survivor -> no producer; the latch stays down.

    cd core/multiplex5000 && ../../.venv/bin/python scripts/xenium_typing/verdict.py [--donors 6436 6442]
"""
from __future__ import annotations

import argparse
import json
import math
import pathlib
import sys
import time

import numpy as np

HERE = pathlib.Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import pins  # noqa: E402

SURVIVE_FRAC = 2.0 / 3.0
DONORS = ("6436", "6442")
SCOPE_SENTENCE = "this arm at this configuration on donors {donors}"


def load(arm: str, donor: str) -> dict | None:
    p = pins.gate_json(arm, donor)
    return json.loads(p.read_text()) if p.exists() else None


def margin(v: dict) -> float:
    ms = []
    for L, r in v.get("lineages", {}).items():
        if r.get("verified"):
            a = r.get("award", {})
            if a.get("real_own") is not None and a.get("null_p95_own") is not None:
                ms.append(float(a["real_own"]) - float(a["null_p95_own"]))
    return float(np.mean(ms)) if ms else float("nan")


def external_pins(v: dict) -> int:
    pinned = (v.get("inputs") or {})
    return len([k for k in pinned if k.endswith("_sha256") or k.endswith("_commit")])


def build(donors=DONORS, arms=pins.ARMS) -> dict:
    n_lin = len(pins.LINEAGES)
    need = math.ceil(SURVIVE_FRAC * n_lin)
    rows, per_label = {}, []
    for arm in arms:
        vs = {d: load(arm, d) for d in donors}
        if any(v is None for v in vs.values()):
            rows[arm] = {"judged": False, "missing": [d for d, v in vs.items() if v is None]}
            continue
        nver = {d: int(v["n_verified"]) for d, v in vs.items()}
        survives = all(nver[d] >= need for d in donors)
        shortfall = "; ".join(f"{d}: {nver[d]} verified < {need}" for d in donors if nver[d] < need) or "-"
        rows[arm] = {"judged": True, "verified": nver, "verified_lineages": {d: v["verified_lineages"] for d, v in vs.items()},
                     "needed": need, "survives": survives, "shortfall": shortfall,
                     "mean_award_margin": float(np.nanmean([margin(v) for v in vs.values()])),
                     "cohort_awarded_all": 0, "n_external_pins": max(external_pins(v) for v in vs.values()),
                     "n_cells": {d: int(v["n_cells"]) for d, v in vs.items()},
                     "caveat": next(iter(vs.values())).get("caveat", "")}
        for d, v in vs.items():
            for L, r in v["lineages"].items():
                a, s, st = r["award"], r["spatial"], r["stability"]
                why = []
                if not a.get("awarded"):
                    why.append("core not awarded" if a.get("judgeable", True) else f"unjudgeable ({a.get('reason', 'no calls')})")
                if not s.get("clear"):
                    why.append("not spatially coherent" if s.get("judgeable", True) else f"spatial unjudgeable ({s.get('reason', '')})")
                if st.get("invariant") is None:
                    why.append("stability unjudged (no replicates)")
                elif not st.get("invariant"):
                    why.append("not invariant across replicates")
                per_label.append({"arm": arm, "donor": d, "lineage": L, "n_calls": int(v["label_counts"].get(L, 0)),
                                  "verified": bool(r["verified"]), "blocking": "; ".join(why) or "-"})
    survivors = [a for a, r in rows.items() if r.get("survives")]
    key = lambda a: (-sum(rows[a]["verified"].values()), -np.nan_to_num(rows[a]["mean_award_margin"], nan=-1e9),  # noqa: E731
                     -rows[a]["cohort_awarded_all"], rows[a]["n_external_pins"])
    ranked = sorted(survivors, key=key)
    producer = ranked[0] if ranked else None
    statement = (f"producer: {producer} -- ranked survivors {ranked}" if producer
                 else f"no arm survives (needs >= {need} of {n_lin} lineages verified on BOTH donors); no producer is chosen and the registry latch stays down")
    return {"method_version": "xenium-typing-verdict-2026-09-07", "written_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
            "donors": list(donors), "n_lineages": n_lin, "survive_frac": SURVIVE_FRAC, "verified_needed_per_donor": need,
            "scope_sentence": SCOPE_SENTENCE.format(donors=" and ".join(donors)),
            "rule": ("verified = core awarded AND spatially coherent AND replicate-invariant; survive iff verified >= "
                     "ceil(survive_frac x n_lineages) on BOTH donors; rank by verified sum, mean award margin, cohort coverage, "
                     "external pins, wall time"),
            "arms": rows, "per_label": per_label, "survivors_ranked": ranked, "producer": producer, "producer_statement": statement}


def to_markdown(s: dict) -> str:
    out = [f"# Xenium producer verdict ({s['method_version']}, {s['written_utc']})", "",
           f"Scope: {s['scope_sentence']}. Survival needs >= {s['verified_needed_per_donor']} of {s['n_lineages']} lineages verified on both donors.", "",
           "| arm | " + " | ".join(f"{d} verified" for d in s["donors"]) + " | survives | shortfall | mean award margin | external pins |", "|---|" + "---|" * (len(s["donors"]) + 4)]
    for a, r in s["arms"].items():
        if not r.get("judged"):
            out.append(f"| {a} | " + " | ".join("not judged" for _ in s["donors"]) + f" | no | missing gate for {r['missing']} | - | - |"); continue
        out.append(f"| {a} | " + " | ".join(str(r["verified"][d]) for d in s["donors"]) + f" | {'yes' if r['survives'] else 'no'} | {r['shortfall']} | {r['mean_award_margin']:.3f} | {r['n_external_pins']} |")
    out += ["", "## Per-lineage (arm x donor x lineage)", "", "| arm | donor | lineage | n calls | verified | blocking |", "|---|---|---|---|---|---|"]
    out += [f"| {p['arm']} | {p['donor']} | {p['lineage']} | {p['n_calls']:,} | {'yes' if p['verified'] else 'no'} | {p['blocking']} |" for p in s["per_label"]]
    out += ["", f"**{s['producer_statement']}**", "", f"Rule: {s['rule']}.", ""]
    cav = next((r.get("caveat") for r in s["arms"].values() if r.get("caveat")), "")
    if cav:
        out += ["Family-wise caveat (verbatim from the gate): " + cav, ""]
    return "\n".join(out)


def main(argv=None) -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--donors", nargs="+", default=list(DONORS))
    a = ap.parse_args(argv)
    s = build(tuple(a.donors))
    d = pins.OUT / "gate"; d.mkdir(parents=True, exist_ok=True)
    (d / "summary.json").write_text(json.dumps(s, indent=1, default=str) + "\n")
    (d / "summary.md").write_text(to_markdown(s))
    print(f"[xenium verdict] {s['producer_statement']}")
    print(f"[xenium verdict] wrote {d / 'summary.json'} and {d / 'summary.md'}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
