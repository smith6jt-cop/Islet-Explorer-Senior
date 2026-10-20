"""The driver and the verdict must apply the frozen rules mechanically, write the protocol only
below its line, and drive an arm end to end through the same subprocess path a real arm takes.

Nothing here touches cohort data.  The end-to-end test plants seven spatially clustered label
blobs with high values on each label's declared markers, points every pin at a tmp tree, and
runs ``screen.main`` against a fake arm module whose CLI is ``arm_base.cli`` -- so the record
schema, alignment, ablations, subsample replicates, the gate, the state tier and the agreement
CSV are all exercised exactly as they will be on donors 6436 and 6442.
"""
from __future__ import annotations

import json
import math
import pathlib
import shutil
import sys

import numpy as np
import pandas as pd
import pytest

HERE = pathlib.Path(__file__).resolve().parents[1]
sys.path.insert(0, str(HERE))
import arm_base  # noqa: E402
import declarations  # noqa: E402
import panel  # noqa: E402  (scripts/celltype/panel.py, on sys.path via pins)
import pins  # noqa: E402
import screen  # noqa: E402
import verdict  # noqa: E402
from criteria import CRITERIA  # noqa: E402

DEC = declarations.load()
LABELS = list(DEC["labels"])
DONORS = list(pins.DONORS)


# --------------------------------------------------------------------------- synthetic verdict JSONs
def _heldout(markers, kind, margin):
    per = []
    for m in markers:
        if kind == "verified":
            per.append({"marker": m, "judgeable": True, "clear": True, "h": 3.0 + margin, "null_p": 3.0, "reason": ""})
        elif kind == "unjudgeable":
            per.append({"marker": m, "judgeable": False, "clear": False, "h": None, "null_p": None,
                        "reason": f"{m}: tail cells in pool: only 5 usable observations, need >= 20"})
        else:
            per.append({"marker": m, "judgeable": True, "clear": False, "h": 1.0, "null_p": 1.5, "reason": ""})
    judged = [p for p in per if p["judgeable"]]
    clear = [p for p in judged if p["clear"]]
    need = int(math.ceil(CRITERIA["heldout_frac"] * len(judged))) if judged else None
    return {"per_marker": per, "n_judged": len(judged), "n_clear": len(clear), "n_required": need,
            "verified": bool(judged and len(judged) >= 2 and len(clear) >= need), "reason": ""}


def _label_result(L, kind, margin=0.5):
    spec = DEC["labels"][L]
    gated = spec["held_out_gated"]
    if kind == "verified":
        spatial = {"I": 0.4, "null_p": 0.02, "n_calls": 500, "judgeable": True, "clear": True, "reason": ""}
    elif kind == "unjudgeable":
        spatial = {"I": None, "null_p": None, "n_calls": 10, "judgeable": False, "clear": False,
                   "reason": f"{L}: calls: only 10 usable observations, need >= 50"}
    else:
        spatial = {"I": 0.01, "null_p": 0.02, "n_calls": 500, "judgeable": True, "clear": False, "reason": ""}
    res = {"label": L, "declared_markers": spec["markers"], "held_out_gated": gated,
           "heldout": _heldout(spec["markers"], kind, margin) if gated else {"per_marker": [], "verified": None},
           "spatial": spatial,
           "stability": {"per_replicate": [{"replicate": f"subsample:{i}", "spatial_clear": spatial["clear"], "judgeable": True}
                                           for i in range(1, 6)], "invariant": True},
           "n_primary_calls": spatial["n_calls"]}
    res["verified"] = bool(spatial["clear"] and (res["heldout"]["verified"] if gated else True))
    return res


def _verdict(arm, donor, kinds, *, margin=0.5, ext_pins=0, wall=10.0, far=None):
    """A verdict JSON in the shape screen.py writes; ``kinds`` maps label -> verified/failed/unjudgeable.
    No ``status`` key: the verdict must classify from the statistics itself."""
    return {"method_version": "test", "arm": arm, "donor_id": donor, "scope": "frame",
            "labels": {L: _label_result(L, kinds.get(L, "verified"), margin) for L in LABELS},
            "n_external_pins": ext_pins, "configs": {"primary": {"wall_seconds": wall}},
            "records_used": [], "noise_false_award_rate": far}


def _write(out, j):
    p = out / "verdict" / j["arm"] / f"donor_id={j['donor_id']}.json"
    p.parent.mkdir(parents=True, exist_ok=True)
    p.write_text(json.dumps(j))


def test_survival_rule_needs_both_donors_and_the_ceiling(tmp_path, monkeypatch):
    out = tmp_path / "out"
    monkeypatch.setattr(pins, "OUT", out)
    all_ok = {}
    _write(out, _verdict("alpha", DONORS[0], all_ok, ext_pins=2))
    _write(out, _verdict("alpha", DONORS[1], all_ok, ext_pins=2))
    _write(out, _verdict("beta", DONORS[0], all_ok, margin=1.0))
    _write(out, _verdict("beta", DONORS[1], {"Endocrine": "failed"}, margin=1.0))
    _write(out, _verdict("gamma", DONORS[0], all_ok))
    _write(out, _verdict("gamma", DONORS[1], {"T_NK": "unjudgeable", "B_Plasma": "unjudgeable", "Myeloid": "unjudgeable"}))
    _write(out, _verdict("delta", DONORS[0], all_ok))                       # one donor only
    _write(out, _verdict("eps", DONORS[0], {"T_NK": "failed", "B_Plasma": "failed"}, far={"false_award_rate_noise": 0.01}))
    _write(out, _verdict("eps", DONORS[1], {"T_NK": "failed", "B_Plasma": "unjudgeable"}, far={"false_award_rate_noise": 0.0}))
    # a verdict at another scope must not be counted at frame scope
    other = _verdict("zeta", DONORS[0], all_ok); other["scope"] = "donor"
    p = out / "verdict" / "zeta" / f"donor_id={DONORS[0]}__donor.json"; p.parent.mkdir(parents=True); p.write_text(json.dumps(other))

    vs = verdict.load_verdicts("frame")
    assert ("zeta", DONORS[0]) not in vs and len(vs) == 9
    s = verdict.build_summary(vs)
    assert s["n_declared_labels"] == 7 and s["verified_needed_per_donor"] == 5
    a = s["arms"]
    assert a["alpha"]["survives"] and a["beta"]["survives"] and a["eps"]["survives"]
    assert not a["gamma"]["survives"] and a["gamma"]["n_verified"][DONORS[1]] == 4
    assert not a["delta"]["survives"] and f"{DONORS[1]}: not screened" in a["delta"]["survival_shortfall"]
    # classification from the statistics, not from a status key
    g = a["gamma"]["per_donor"][DONORS[1]]
    assert g["T_NK"]["status"] == "unjudgeable" and "spatial:" in g["T_NK"]["blocking_reasons"][0]
    e = a["eps"]["per_donor"][DONORS[0]]
    assert e["T_NK"]["status"] == "failed" and any(r.startswith("spatial: I=") for r in e["T_NK"]["blocking_reasons"])
    assert a["beta"]["per_donor"][DONORS[1]]["Endocrine"]["status"] == "failed"
    # ranking (i): verified total 14 > 13 > 10; producer is mechanical
    assert [r["arm"] for r in s["survivors_ranked"]] == ["alpha", "beta", "eps"]
    assert s["producer"] == "alpha" and s["producer_statement"] == "producer: alpha"
    md = verdict.render_markdown(s)
    assert md.count(verdict.scope_sentence()) >= 3          # every table and the producer line
    assert "false_award_rate_noise" in md and "eps:" in md
    assert "| gamma | " in md and "| delta | " + DONORS[1] + " | (all) | not screened" in md
    pj, pm = verdict.write_summary(s)
    assert json.loads(pj.read_text())["producer"] == "alpha" and pm.read_text() == md


def test_ranking_tiebreakers_and_no_producer():
    def summary(*vs):
        return verdict.build_summary({(v["arm"], v["donor_id"]): v for v in vs})
    both = lambda arm, **kw: [_verdict(arm, d, {}, **kw) for d in DONORS]  # noqa: E731
    # (ii) equal verified totals -> larger mean held-out margin wins
    s = summary(*both("a", margin=0.5), *both("b", margin=1.0))
    assert [r["arm"] for r in s["survivors_ranked"]] == ["b", "a"]
    # (iii) equal margins -> fewer external pins
    s = summary(*both("a", ext_pins=2), *both("b", ext_pins=0))
    assert s["producer"] == "b"
    # (iii) equal pins -> lower primary wall time (summed over donors)
    s = summary(*both("a", wall=10.0), *both("b", wall=5.0))
    assert s["producer"] == "b" and s["arms"]["b"]["primary_wall_seconds"] == 10.0
    # no survivor -> the latch statement, verbatim
    s = summary(*[_verdict("a", d, {"T_NK": "failed", "B_Plasma": "failed", "Myeloid": "failed"}) for d in DONORS])
    assert s["producer"] is None and s["producer_statement"] == verdict.NO_PRODUCER
    assert verdict.NO_PRODUCER in verdict.render_markdown(s)
    # an empty screen still renders
    s = summary()
    assert s["producer"] is None and "no verdict records" in verdict.render_markdown(s)


def test_label_status_prefers_failed_over_unjudgeable_and_reports_stability():
    r = _label_result("T_NK", "failed")
    r["heldout"]["per_marker"][0].update(judgeable=False, clear=False, reason="CD3e: tail cells in pool: only 3")
    status, reasons = verdict.label_status(r)
    assert status == "failed" and any(x.startswith("heldout:") for x in reasons) and any("not judged" in x for x in reasons)
    r = _label_result("Endocrine", "verified")
    r["stability"]["per_replicate"][2]["spatial_clear"] = False
    r["stability"]["invariant"] = False
    r["verified"] = False
    status, reasons = verdict.label_status(r)
    assert status == "failed" and reasons == ["stability: spatial verdict flips on ['subsample:3']"]
    assert verdict.label_status({"error": "primary_not_available (rc=1)"}) == ("failed", ["primary_not_available (rc=1)"])


# --------------------------------------------------------------------------- the protocol append
def test_append_never_modifies_text_above_the_line(tmp_path):
    proto = tmp_path / "PROTOCOL.md"
    shutil.copy(verdict.PROTOCOL, proto)
    before = proto.read_text()
    frozen = before.split(verdict.SECTION4_HEADER)[0] + verdict.SECTION4_HEADER
    verdict.append_to_protocol("### a table\n\n| x |\n|---|\n| 1 |\n", proto, stamp="T1")
    mid = proto.read_text()
    assert mid.startswith(before) and mid.startswith(frozen)
    assert mid.count(verdict.SECTION4_HEADER) == 1 and "## Screen summary, appended T1" in mid
    verdict.append_to_protocol("second\n", proto, stamp="T2")
    after = proto.read_text()
    assert after.startswith(mid) and after.startswith(frozen)
    assert after.index("appended T1") < after.index("appended T2")
    bare = tmp_path / "no_header.md"
    bare.write_text("# something else\n")
    with pytest.raises(RuntimeError):
        verdict.append_to_protocol("x", bare, stamp="T3")
    assert bare.read_text() == "# something else\n"


# --------------------------------------------------------------------------- the plan
def test_config_plan_ablates_held_out_markers_only_and_replicates_by_kind(monkeypatch):
    present = ["GCG", "INS", "SST", "CD31", "Caveolin", "CD3e", "CD4", "SMA"]
    monkeypatch.setattr(arm_base, "identity_markers_for", lambda d: present)
    plan = dict(screen.config_plan("x", False, DEC))
    assert "primary" in plan and plan["primary"] == ["--seed", "0"]
    abl = sorted(c.split(":", 1)[1] for c in plan if c.startswith("ablate:"))
    assert abl == ["CD31", "CD3e", "CD4", "Caveolin", "SMA"]          # no hormone (Endocrine is not held-out gated)
    subs = [c for c in plan if c.startswith("subsample:")]
    assert subs == [f"subsample:{i}" for i in range(1, 6)] and "--subsample-frac" in plan["subsample:1"]
    seeds = [c for c in dict(screen.config_plan("x", True, DEC)) if c.startswith("seed:")]
    assert seeds == [f"seed:{i}" for i in range(1, 6)]
    assert screen.inspect_arm("not-here", HERE / "arms")["present"] is False


# --------------------------------------------------------------------------- end to end on a fake arm
IDENTITY = {"Endocrine": ["GCG", "INS", "SST"], "Epithelial_nonendocrine": ["E_cadherin", "EpCAM"],
            "Vascular": ["CD31", "Caveolin"], "Stromal": ["SMA", "Vimentin"], "T_NK": ["CD3e", "CD4"],
            "B_Plasma": ["CD20", "CD79a"], "Myeloid": ["CD68", "CD163"]}
STATE = ["Ki67", "HLA_DR"]
DONOR = "9001"
CENTRES = [(300, 300), (300, 900), (300, 1500), (900, 300), (900, 900), (900, 1500), (1500, 900)]

FAKE_ARM = '''
"""Fake arm for the screen test: mean within-donor rank over each label's acquired declared
markers, argmax with a margin.  Not a method -- it exists so the subprocess path, alignment,
ablations, replicates, gate, state tier and agreement CSV are exercised on planted structure."""
from __future__ import annotations
import os, pathlib, sys
import numpy as np, pandas as pd
sys.path.insert(0, __BROADTYPE__)
import pins, arm_base, declarations, panel
# The test's world: point the lane's pins at the synthetic tree (env set by the test).
pins.OUT = pathlib.Path(os.environ["BT_TEST_OUT"])
pins.uniform_parquet = lambda d: pathlib.Path(os.environ["BT_TEST_UNIFORM"]) / f"donor_id={d}" / "data.parquet"
pins.assert_pins = lambda: {"run_id": "synthetic"}
_markers = os.environ["BT_TEST_MARKERS"].split(",")
panel.donor_markers = lambda donor=None: list(_markers)

ARM = "fake-arm"
STOCHASTIC = False
MARGIN = 0.15

def run(ctx):
    dec = declarations.load()
    R = arm_base.within_donor_rank(ctx.X)
    S = {}
    for L, spec in dec["labels"].items():
        cols = [m for m in spec["markers"] if m in ctx.markers]
        if cols:
            S[L] = R[cols].mean(axis=1).to_numpy()
    S = pd.DataFrame(S)
    A = S.to_numpy(float)
    order = np.argsort(-np.where(np.isnan(A), -np.inf, A), axis=1)
    top = S.columns.to_numpy()[order[:, 0]]
    s1 = np.take_along_axis(A, order[:, [0]], 1)[:, 0]
    s2 = np.take_along_axis(A, order[:, [1]], 1)[:, 0]
    ok = np.isfinite(s1)
    clear = ok & (s1 - np.where(np.isfinite(s2), s2, -np.inf) >= MARGIN)
    call = np.where(clear, top, "Ambiguous").astype(object)
    reason = np.where(clear, "", "low_margin").astype(object)
    call[~ok] = "Unavailable"
    reason[~ok] = "no_identity_marker_measured"
    return pd.DataFrame({"object_id": ctx.ids, "label": call, "score": s1, "reason": reason})

if __name__ == "__main__":
    raise SystemExit(arm_base.cli(ARM, run, env=sys.executable,
                                  default_params={"margin": MARGIN, "features": "within_donor_rank mean over declared markers"}))
'''


def _synthetic_world(tmp_path, per_blob=120, seed=0):
    """Two regions x seven spatial blobs; each blob's cells carry +4 on their label's declared
    markers.  Endocrine blobs sit in QuPath 'core' (audit column only).  Three cells have no
    identity measurement at all and must come back Unavailable, never dropped."""
    rng = np.random.default_rng(seed)
    markers = [m for ms in IDENTITY.values() for m in ms]
    rows = []
    for region in range(2):
        for (L, ms), (cx, cy) in zip(IDENTITY.items(), CENTRES):
            X = rng.normal(0.0, 1.0, (per_blob, len(markers)))
            for m in ms:
                X[:, markers.index(m)] += 4.0
            for i in range(per_blob):
                rows.append({"object_id": f"{region}-{L}-{i:04d}", "region_idx": region,
                             "X_centroid": cx + 5000 * region + rng.normal(0, 45), "Y_centroid": cy + rng.normal(0, 45),
                             "cell_region": "core" if L == "Endocrine" else "tissue",
                             "islet_num": 1 if L == "Endocrine" else 0, **dict(zip(markers, X[i]))})
    for i in range(3):
        rows.append({"object_id": f"nan-{i}", "region_idx": 0, "X_centroid": 100.0 + i, "Y_centroid": 100.0,
                     "cell_region": "tissue", "islet_num": 0, **{m: np.nan for m in markers}})
    tab = pd.DataFrame(rows)
    # a sprinkle of single missing measurements (NaN, never 0)
    for _ in range(10):
        tab.loc[rng.integers(len(tab) - 3), markers[rng.integers(len(markers))]] = np.nan
    for m in STATE:
        tab[m] = rng.normal(0, 1, len(tab))
    tab.insert(0, "donor_id", DONOR)
    tab.insert(2, "qc_analysis_eligible", True)
    uniform = tmp_path / "uniform"
    (uniform / f"donor_id={DONOR}").mkdir(parents=True)
    tab.drop(columns=["region_idx"]).to_parquet(uniform / f"donor_id={DONOR}" / "data.parquet", index=False)
    out = tmp_path / "out"
    (out / "frame" / "cells").mkdir(parents=True)
    fc = tab[["object_id", "region_idx"]].copy()
    fc["draw"] = np.where(fc["region_idx"] == 0, "blind", "islet")
    fc["in_expression"] = True
    fc.to_parquet(out / "frame" / "cells" / f"donor_id={DONOR}.parquet", index=False)
    arms = tmp_path / "arms"
    arms.mkdir()
    (arms / "fake_arm.py").write_text(FAKE_ARM.replace("__BROADTYPE__", repr(str(HERE))))
    return out, uniform, arms, markers + STATE, tab


def test_fake_arm_end_to_end_through_screen(tmp_path, monkeypatch, capsys):
    out, uniform, arms, all_markers, tab = _synthetic_world(tmp_path)
    monkeypatch.setattr(pins, "OUT", out)
    monkeypatch.setattr(pins, "uniform_parquet", lambda d: uniform / f"donor_id={d}" / "data.parquet")
    monkeypatch.setattr(pins, "assert_pins", lambda: {"run_id": "synthetic"})
    monkeypatch.setattr(panel, "donor_markers", lambda donor=None: list(all_markers))
    monkeypatch.setenv("BT_TEST_OUT", str(out))
    monkeypatch.setenv("BT_TEST_UNIFORM", str(uniform))
    monkeypatch.setenv("BT_TEST_MARKERS", ",".join(all_markers))

    rc = screen.main(["--donors", DONOR, "--arms", "fake-arm", "not-written-yet", "--arms-dir", str(arms),
                      "--python", sys.executable])
    assert rc == 0
    text = capsys.readouterr().out
    assert "skip not-written-yet" in text and "is absent" in text

    pj = screen.verdict_json_path("fake-arm", DONOR, "frame")
    raw = pj.read_text()
    assert "NaN" not in raw                                       # NaN is written as null, never as a bare token
    j = json.loads(raw)
    assert j["n_cells"] == len(tab) and j["n_regions"] == 2 and j["scope"] == "frame"
    assert j["scope_sentence"] == verdict.scope_sentence()
    # the plan: primary + one ablation per held-out marker acquired + 5 subsample replicates
    held_out = sorted(m for L, ms in IDENTITY.items() if DEC["labels"][L]["held_out_gated"] for m in ms)
    assert sorted(j["ablation_markers"]) == held_out and len(held_out) == 12
    assert j["replicate_kind"] == "subsample"
    assert set(j["configs"]) == {"primary", *(f"ablate:{m}" for m in held_out), *(f"subsample:{i}" for i in range(1, 6))}
    for cid, c in j["configs"].items():
        assert c["ok"] and c["returncode"] == 0 and c["wall_seconds"] > 0 and c["arm_seconds"] is not None, cid
        assert pathlib.Path(c["calls"]).exists() and pathlib.Path(c["log"]).exists()
    # tiering evidence in the records: the ablated marker is gone from identity, no state marker ever entered
    recs = {r["config_id"]: r for r in j["records_used"]}
    assert "CD3e" not in recs["ablate:CD3e"]["identity_markers"] and "CD3e" in recs["primary"]["identity_markers"]
    assert not (set(recs["primary"]["identity_markers"]) & set(STATE)) and j["state_marker_in_identity_inputs"] == []
    # planted structure recovered and verified for every declared label
    assert set(j["labels"]) == set(LABELS)
    for L, r in j["labels"].items():
        assert r["status"] == "verified" and r["verified"], (L, r["status"], r["blocking_reasons"])
        assert r["spatial"]["clear"] and r["stability"]["invariant"] and len(r["stability"]["per_replicate"]) == 5
        assert all(x["n_in_replicate"] == int(round(0.9 * len(tab))) for x in r["stability"]["per_replicate"])
        if r["held_out_gated"]:
            assert r["heldout"]["n_judged"] == len(IDENTITY[L]) and r["heldout"]["n_clear"] == len(IDENTITY[L])
            assert all(p["h"] > p["null_p"] for p in r["heldout"]["per_marker"])
        assert j["primary_label_counts"][L] >= 200
    assert j["n_verified"] == 7
    assert j["labels"]["Endocrine"]["islet_core_audit"]["audit_only"] and j["labels"]["Endocrine"]["islet_core_audit"]["enrichment"] > 1
    assert j["primary_label_counts"].get("Unavailable", 0) >= 3       # the three unmeasured cells got a row
    # the state tier: appended after identity, one row per frozen cell, continuous ranks only
    cws = pd.read_parquet(j["calls_with_state"]["path"])
    assert len(cws) == len(tab) and set(cws["object_id"]) == set(tab["object_id"])
    assert {"state_Ki67_donor_rank", "state_HLA_DR_donor_rank"} <= set(cws.columns)
    assert cws["state_Ki67_donor_rank"].between(0, 1).all()
    assert j["calls_with_state"]["state_columns"] == ["state_Ki67_donor_rank", "state_HLA_DR_donor_rank"]
    # agreement CSV: one arm -> header only, descriptive
    agr = pd.read_csv(screen.agreement_path(DONOR, "frame"))
    assert list(agr.columns)[:4] == ["donor_id", "scope", "arm_a", "arm_b"] and len(agr) == 0
    # a second run with --skip-existing reuses every calls parquet and runs no subprocess
    rc = screen.main(["--donors", DONOR, "--arms", "fake-arm", "--arms-dir", str(arms), "--skip-existing", "--no-agreement"])
    assert rc == 0
    j2 = json.loads(pj.read_text())
    assert all(c["reused"] and c["wall_seconds"] is None for c in j2["configs"].values())
    assert j2["n_verified"] == 7

    # the verdict on this one synthetic donor: survives, is the producer, appends only below the line
    proto = tmp_path / "PROTOCOL.md"
    shutil.copy(verdict.PROTOCOL, proto)
    before = proto.read_text()
    rc = verdict.main(["--donors", DONOR, "--protocol", str(proto)])
    assert rc == 0
    s = json.loads((out / "verdict" / "summary.json").read_text())
    assert s["producer"] == "fake-arm" and s["arms"]["fake-arm"]["survives"]
    assert s["arms"]["fake-arm"]["n_verified"][DONOR] == 7 and s["arms"]["fake-arm"]["n_cleared_heldout_markers"] == 12
    assert s["arms"]["fake-arm"]["mean_heldout_margin"] > 0 and s["arms"]["fake-arm"]["n_external_pins"] == 0
    after = proto.read_text()
    assert after.startswith(before) and "producer: fake-arm" in after and f"donors {DONOR}" in after
    # and with the real donor list the synthetic arm cannot survive: neither donor was screened
    s2 = verdict.build_summary(verdict.load_verdicts("frame"))
    assert not s2["arms"]["fake-arm"]["survives"] and s2["producer_statement"] == verdict.NO_PRODUCER
