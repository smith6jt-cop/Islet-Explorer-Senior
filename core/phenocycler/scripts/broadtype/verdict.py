"""The verdict: the survival rule and the A.3 producer rule, applied mechanically to what the
screen recorded, and the summary appended BELOW the protocol's line.

    cd core/phenocycler && ../../.venv/bin/python scripts/broadtype/verdict.py [--scope frame] [--no-append]

Why this is a separate module from ``screen.py``: the screen measures, the verdict decides, and
the decision must be re-derivable from the JSON alone -- anyone can re-run this on the recorded
statistics without re-running an arm.  Nothing here reads a call, a value, or a QuPath class; it
reads ``verdict/<arm>/donor_id=<d>.json`` and applies the numbers frozen in ``criteria.py``.

The rules, verbatim from the protocol (Section 3):
  * a label is verified when ``gate.judge_label`` said so (spatial coherence, held-out enrichment
    where gated, stability where replicates exist -- all blocking);
  * an arm survives iff verified >= ceil(survive_frac x n_declared) labels on BOTH donors
    (7 declared labels, 2/3 -> 5);
  * producer choice among survivors, in order: (i) verified-label count summed over both donors,
    (ii) mean held-out margin (h - null p95) over the markers that cleared, (iii) fewest external
    pins, then lowest primary wall time.  Agreement with another arm never enters.
  * no survivor -> "no producer chosen; latch stays down".

Every table carries the scope sentence and, where an arm's records carry one, the noise
false-award rate its tests measured.  ``append_to_protocol`` only ever appends after Section 4's
header; the test asserts the text above the line is byte-identical before and after.
"""
from __future__ import annotations

import argparse
import json
import math
import pathlib
import sys
import time

HERE = pathlib.Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import declarations  # noqa: E402
import pins  # noqa: E402
from criteria import CRITERIA  # noqa: E402

METHOD_VERSION = "broadtype-verdict-2026-09-05"
PROTOCOL = pins.REPO / "docs" / "BROADTYPE_PROTOCOL_2026-09.md"
SECTION4_HEADER = "# 4. Results"
#: pinned keys every record carries; anything beyond them is an external pin (ranking key iii)
BASE_PINS = ("run_id", "expression_content_id", "declarations")
NO_PRODUCER = "no producer chosen; latch stays down"

STATUS_VERIFIED, STATUS_UNJUDGEABLE, STATUS_FAILED = "verified", "unjudgeable", "failed"


def scope_sentence(donors=pins.DONORS) -> str:
    return "this arm at this configuration on donors " + " and ".join(str(d) for d in donors)


def _f(x, nd: int = 3) -> str:
    """Format a statistic that may have come back from JSON as None (NaN is written as null)."""
    if x is None:
        return "nan"
    try:
        v = float(x)
    except (TypeError, ValueError):
        return str(x)
    return "nan" if not math.isfinite(v) else f"{v:.{nd}f}"


def _ceil_frac(frac: float, n: int) -> int:
    """ceil(frac * n) without the float artefact that makes 2/3 * 6 land on 4.0000000000000004."""
    return int(math.ceil(round(frac * n, 9)))


# --------------------------------------------------------------------------- one label
def label_status(res: dict) -> tuple[str, list[str]]:
    """Classify one ``gate.judge_label`` result as verified / unjudgeable / failed, with the
    blocking reasons written out.

    ``failed`` means a blocking statistic was judged and did not clear (or the arm itself did not
    produce calls); ``unjudgeable`` means every shortfall is an evidence floor (too few calls or
    tail cells) -- which is never a pass, but is a different finding from a failure.  A label
    with both kinds of shortfall is reported ``failed``: an arm does not get credit for the
    statistic it starved.
    """
    if res.get("error"):
        return STATUS_FAILED, [str(res["error"])]
    if res.get("verified"):
        return STATUS_VERIFIED, []
    failed, unjudgeable = [], []
    sp = res.get("spatial") or {}
    if sp:
        if not sp.get("judgeable"):
            unjudgeable.append(f"spatial: {sp.get('reason') or 'not judgeable'}")
        elif not sp.get("clear"):
            failed.append(f"spatial: I={_f(sp.get('I'))} <= p95={_f(sp.get('null_p'))} (n_calls={sp.get('n_calls')})")
    if res.get("held_out_gated"):
        ho = res.get("heldout") or {}
        per = ho.get("per_marker") or []
        judged = [p for p in per if p.get("judgeable")]
        if len(judged) < 2:
            short = "; ".join(f"{p.get('marker')}: {p.get('reason')}" for p in per if not p.get("judgeable"))
            unjudgeable.append(f"heldout: {len(judged)} judgeable marker(s), need >= 2 ({short or 'no markers'})")
        elif not ho.get("verified"):
            miss = [p for p in judged if not p.get("clear")]
            detail = "; ".join(f"{p['marker']}: h={_f(p.get('h'))} <= p95={_f(p.get('null_p'))}" for p in miss)
            failed.append(f"heldout: {ho.get('n_clear')}/{ho.get('n_judged')} markers clear, need {ho.get('n_required')} ({detail})")
            starved = [p for p in per if not p.get("judgeable")]
            if starved:
                unjudgeable.append("heldout (not judged): " + "; ".join(f"{p.get('marker')}: {p.get('reason')}" for p in starved))
    st = res.get("stability") or {}
    reps = st.get("per_replicate") or []
    if reps and st.get("invariant") is False:
        flips = [r["replicate"] for r in reps if r.get("spatial_clear") != sp.get("clear")]
        failed.append(f"stability: spatial verdict flips on {flips}")
    if failed:
        return STATUS_FAILED, failed + unjudgeable
    if unjudgeable:
        return STATUS_UNJUDGEABLE, unjudgeable
    return STATUS_FAILED, ["not verified; no blocking statistic recorded"]


def heldout_margins(res: dict) -> list[float]:
    """h - null_p95 for every held-out marker that cleared (ranking key ii)."""
    out = []
    for p in (res.get("heldout") or {}).get("per_marker") or []:
        if p.get("judgeable") and p.get("clear") and p.get("h") is not None and p.get("null_p") is not None:
            out.append(float(p["h"]) - float(p["null_p"]))
    return out


# --------------------------------------------------------------------------- loading
def verdict_dir(out: pathlib.Path | None = None) -> pathlib.Path:
    return (out or pins.OUT) / "verdict"


def load_verdicts(scope: str = "frame", out: pathlib.Path | None = None) -> dict:
    """(arm, donor) -> the screen's JSON, for one scope.  Frame verdicts are ``donor_id=<d>.json``;
    any other scope carries a ``__<scope>`` suffix so the two never overwrite each other."""
    found = {}
    for p in sorted(verdict_dir(out).glob("*/donor_id=*.json")):
        try:
            j = json.loads(p.read_text())
        except json.JSONDecodeError as e:  # a half-written file from a screen still running
            print(f"[verdict] skip {p}: {e}", file=sys.stderr)
            continue
        if j.get("scope", "frame") != scope:
            continue
        found[(str(j["arm"]), str(j["donor_id"]))] = j
    return found


def _n_external_pins(j: dict) -> int:
    if j.get("n_external_pins") is not None:
        return int(j["n_external_pins"])
    n = 0
    for r in j.get("records_used") or []:
        n = max(n, len(set((r.get("pinned") or {})) - set(BASE_PINS)))
    return n


def _primary_wall(j: dict) -> float | None:
    cfg = (j.get("configs") or {}).get("primary") or {}
    w = cfg.get("wall_seconds")
    if w is None:
        w = (cfg.get("record") or {}).get("seconds")
    return None if w is None else float(w)


# --------------------------------------------------------------------------- the rules
def build_summary(verdicts: dict, *, donors=pins.DONORS, scope: str = "frame", dec: dict | None = None) -> dict:
    dec = dec or declarations.load()
    donors = [str(d) for d in donors]
    n_declared = len(dec["labels"])
    need = _ceil_frac(CRITERIA["survive_frac"], n_declared)
    arms = sorted({a for a, _ in verdicts})
    per_arm = {}
    for arm in arms:
        rows, n_ver, margins, walls, pins_n, far = {}, {}, [], [], 0, {}
        for d in donors:
            j = verdicts.get((arm, d))
            if j is None:
                rows[d] = None
                n_ver[d] = None
                continue
            labs = {}
            for L in dec["labels"]:
                res = (j.get("labels") or {}).get(L) or {"error": "label_not_judged"}
                status, reasons = (res.get("status"), res.get("blocking_reasons"))
                if status is None:
                    status, reasons = label_status(res)
                m = heldout_margins(res)
                margins.extend(m)
                labs[L] = {"status": status, "blocking_reasons": list(reasons or []),
                           "heldout_margins": m, "n_primary_calls": res.get("n_primary_calls")}
            rows[d] = labs
            n_ver[d] = sum(1 for v in labs.values() if v["status"] == STATUS_VERIFIED)
            pins_n = max(pins_n, _n_external_pins(j))
            w = _primary_wall(j)
            if w is not None:
                walls.append(w)
            if j.get("noise_false_award_rate"):
                far[d] = j["noise_false_award_rate"]
        screened = [d for d in donors if rows[d] is not None]
        survives = len(screened) == len(donors) and all(n_ver[d] >= need for d in donors)
        why = []
        for d in donors:
            if rows[d] is None:
                why.append(f"{d}: not screened")
            elif n_ver[d] < need:
                why.append(f"{d}: {n_ver[d]} verified < {need}")
        per_arm[arm] = {
            "per_donor": rows,
            "n_verified": n_ver,
            "verified_total": sum(v for v in n_ver.values() if v is not None),
            "mean_heldout_margin": (sum(margins) / len(margins)) if margins else None,
            "n_cleared_heldout_markers": len(margins),
            "n_external_pins": pins_n,
            "primary_wall_seconds": (sum(walls) if walls else None),
            "noise_false_award_rate": far or None,
            "survives": bool(survives),
            "survival_shortfall": why,
        }
    survivors = [a for a in arms if per_arm[a]["survives"]]

    def key(a):
        s = per_arm[a]
        margin = s["mean_heldout_margin"]
        wall = s["primary_wall_seconds"]
        return (-s["verified_total"],
                -(margin if margin is not None else float("-inf")),
                s["n_external_pins"],
                wall if wall is not None else float("inf"))

    ranked = sorted(survivors, key=key)
    producer = ranked[0] if ranked else None
    return {
        "method_version": METHOD_VERSION,
        "written_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
        "scope": scope,
        "scope_sentence": scope_sentence(donors),
        "donors": donors,
        "criteria": dict(CRITERIA),
        "declarations": dec.get("method_version"),
        "n_declared_labels": n_declared,
        "verified_needed_per_donor": need,
        "rule": ("survive iff verified >= ceil(survive_frac x n_declared) on BOTH donors; producer = first "
                 "survivor by (i) verified total, (ii) mean held-out margin over cleared markers, (iii) fewest "
                 "external pins then lowest primary wall time; agreement never chooses"),
        "arms": per_arm,
        "survivors_ranked": [{"arm": a, "rank_key": [float(k) if isinstance(k, float) else k for k in key(a)]} for a in ranked],
        "producer": producer,
        "producer_statement": (f"producer: {producer}" if producer else NO_PRODUCER),
    }


# --------------------------------------------------------------------------- markdown
def _md_table(header: list[str], rows: list[list[str]]) -> str:
    esc = lambda s: str(s).replace("|", "\\|").replace("\n", " ")  # noqa: E731
    out = ["| " + " | ".join(header) + " |", "|" + "|".join("---" for _ in header) + "|"]
    out += ["| " + " | ".join(esc(c) for c in r) + " |" for r in rows]
    return "\n".join(out)


def _far_line(summary: dict, arms: list[str]) -> str:
    parts = []
    for a in arms:
        far = summary["arms"][a].get("noise_false_award_rate")
        if far:
            parts.append(f"{a}: " + "; ".join(f"{d}: {json.dumps(v, default=str)}" for d, v in far.items()))
    if not parts:
        return "Noise false-award rate: not recorded in any arm's records (only null-gated cluster arms measure one)."
    return "Noise false-award rate measured on synthetic noise (from the arms' records): " + " | ".join(parts) + "."


def render_markdown(summary: dict) -> str:
    sent = summary["scope_sentence"]
    arms = sorted(summary["arms"])
    donors = summary["donors"]
    lines = [
        f"Scope: {sent}. Frame: `{summary['scope']}`. Declarations `{summary['declarations']}`; "
        f"{summary['n_declared_labels']} declared labels; survival needs >= {summary['verified_needed_per_donor']} "
        f"verified per donor on both donors (criteria `survive_frac` = {summary['criteria']['survive_frac']:.4f}). "
        f"Verdict `{summary['method_version']}`, written {summary['written_utc']}.",
        "",
        "### Per-label verdicts (arm x donor x label)",
        "",
    ]
    rows = []
    for a in arms:
        for d in donors:
            labs = summary["arms"][a]["per_donor"].get(d)
            if labs is None:
                rows.append([a, d, "(all)", "not screened", "no verdict JSON for this donor"])
                continue
            for L, v in labs.items():
                rows.append([a, d, L, v["status"], "; ".join(v["blocking_reasons"]) or "-"])
    lines.append(_md_table(["arm", "donor", "label", "status", "blocking reasons"], rows) if rows
                 else "_no verdict records found_")
    lines += ["", f"Scope of every row above: {sent}. {_far_line(summary, arms)}", "", "### Survival", ""]
    rows = []
    for a in arms:
        s = summary["arms"][a]
        rows.append([a, *[("-" if s["n_verified"].get(d) is None else str(s["n_verified"][d])) for d in donors],
                     str(summary["verified_needed_per_donor"]), "yes" if s["survives"] else "no",
                     "; ".join(s["survival_shortfall"]) or "-",
                     _f(s["mean_heldout_margin"]), str(s["n_external_pins"]), _f(s["primary_wall_seconds"], 1)])
    lines.append(_md_table(["arm", *[f"{d} verified" for d in donors], "needed", "survives", "shortfall",
                            "mean held-out margin", "external pins", "primary wall s"], rows) if rows
                 else "_no arms screened_")
    lines += ["", f"Scope of the survival row: {sent}. {_far_line(summary, arms)}", "", "### Producer (A.3)", ""]
    if summary["producer"]:
        lines.append(f"**{summary['producer_statement']}** -- ranked survivors: "
                     + ", ".join(f"{r['arm']} {r['rank_key']}" for r in summary["survivors_ranked"])
                     + f". Rule: {summary['rule']}. Scope: {sent}.")
    else:
        lines.append(f"**{NO_PRODUCER}** -- no arm verified >= {summary['verified_needed_per_donor']} labels on both donors. "
                     f"Rule: {summary['rule']}. Scope: {sent}.")
    return "\n".join(lines) + "\n"


def append_to_protocol(md: str, protocol: pathlib.Path = PROTOCOL, *, stamp: str | None = None) -> str:
    """Append the summary under Section 4.  Append only: the file's previous content is a byte
    prefix of the new content, and a protocol without the Section 4 header is refused."""
    protocol = pathlib.Path(protocol)
    before = protocol.read_text()
    if SECTION4_HEADER not in before:
        raise RuntimeError(f"{protocol}: '{SECTION4_HEADER}' not found; refusing to append anywhere else")
    stamp = stamp or time.strftime("%Y-%m-%d %H:%M UTC", time.gmtime())
    block = f"\n## Screen summary, appended {stamp}\n\n{md.rstrip()}\n"
    with protocol.open("a") as fh:
        fh.write(("" if before.endswith("\n") else "\n") + block)
    after = protocol.read_text()
    if not after.startswith(before):
        raise RuntimeError("append changed text above the line; this must never happen")
    return block


def write_summary(summary: dict, out: pathlib.Path | None = None) -> tuple[pathlib.Path, pathlib.Path]:
    d = verdict_dir(out)
    d.mkdir(parents=True, exist_ok=True)
    pj, pm = d / "summary.json", d / "summary.md"
    pj.write_text(json.dumps(summary, indent=1, default=str) + "\n")
    pm.write_text(render_markdown(summary))
    return pj, pm


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description="broadtype verdict: survival + producer rule from the screen's records")
    ap.add_argument("--scope", choices=("frame", "donor"), default="frame")
    ap.add_argument("--donors", nargs="+", default=list(pins.DONORS), help="both must be screened for survival")
    ap.add_argument("--protocol", default=str(PROTOCOL))
    ap.add_argument("--no-append", action="store_true", help="write summary.{json,md} only")
    a = ap.parse_args(argv)
    verdicts = load_verdicts(a.scope)
    summary = build_summary(verdicts, donors=a.donors, scope=a.scope)
    pj, pm = write_summary(summary)
    print(f"[verdict] {len(verdicts)} (arm, donor) verdicts at scope {a.scope}; "
          f"survivors: {[r['arm'] for r in summary['survivors_ranked']] or 'none'}; {summary['producer_statement']}")
    print(f"[verdict] wrote {pj} and {pm}")
    if not a.no_append:
        append_to_protocol(pm.read_text(), pathlib.Path(a.protocol))
        print(f"[verdict] appended under Section 4 of {a.protocol}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
