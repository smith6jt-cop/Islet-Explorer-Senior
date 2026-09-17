"""The screen driver: every arm is run the same way and judged by the same bar.

    cd core/phenocycler && ../../.venv/bin/python scripts/broadtype/screen.py \\
        --donors 6436 6442 --arms cosine-broad8 tasp ... [--scope frame] [--skip-existing]

Why one driver instead of each arm judging itself: the bar (``criteria.py`` + ``gate.py``) has to
be applied identically to nine methods that live in five environments.  An arm module only ever
sees ``ctx.X`` -- the aligned values of the declared identity (type) markers, minus the ablated
one.  The driver alone holds the ORIGINAL values (the withheld marker's tail is defined on them,
never on the ablated input), the frame geometry (``X_centroid``/``Y_centroid``, ``region_idx``),
the audit column ``cell_region`` (QuPath parentage; audit only, never an input), and it is the
driver -- not the arm -- that appends the state tier, AFTER identity, so no state marker can
enter a call.  Every arm runs through its own CLI as a subprocess so an R or foreign-venv arm
never loads into this process, and its wall time is measured from outside.

Per arm x donor at ``--scope frame`` (default; ``donor`` = whole donor):
  1. primary        seed 0                                       config ``primary``
  2. ablations      one per declared marker of every held-out-gated label acquired on the donor,
                    seed 0                                       config ``ablate:<marker>``
  3. replicates     seeds 1..5 when the arm module declares ``STOCHASTIC = True``; otherwise five
                    90 % hash subsamples                         config ``seed:<n>`` / ``subsample:<n>``
  4. ``gate.judge_label`` for every declared label, with ``other_tails`` = the union of every OTHER
     label's declared-marker tails (descriptive contamination) and the islet-core audit for
     Endocrine (audit only);
  5. the state tier (``arm_base.state_columns``) appended to the primary calls;
  6. cross-arm ``gate.agreement`` on the primary calls of every arm found on disk -- descriptive,
     never a producer.

Outputs, all under ``data/development/broadtype/`` (never ``data/runs/``, never a stage manifest):
  verdict/<arm>/donor_id=<d>.json                    every statistic, the records used, wall times
  verdict/<arm>/donor_id=<d>/calls_with_state.parquet primary calls + state_<m>_donor_rank columns
  verdict/agreement_<d>.csv                           ARI / AMI for every arm pair (descriptive)
  logs/screen/<arm>__<d>__<config>.log                each subprocess's stdout+stderr
(at ``--scope donor`` the JSON, directory and CSV carry a ``__donor`` suffix.)

An arm whose module is absent is skipped with a printed reason (arms are written concurrently);
one whose subprocess fails is recorded with its return code and judged on whatever it produced
(a missing primary makes every label ``failed: primary_not_available``).  A record whose
``identity_markers`` contain a state marker fails every label: the tiering rule is enforced here
as well as in the tests.

Recorded deviation from the protocol text: Section 2's whole-donor fallback for a marker whose
Step-1 tail has fewer than ``min_tail_in_pool`` cells inside the frame is NOT automated here.  Such
a marker is reported unjudgeable in the frame verdict, beside its ``tails_in_frame`` count, and
``--scope donor`` runs the whole-donor screen into a separate verdict file for ``verdict.py
--scope donor``.  Nothing is silently substituted.
"""
from __future__ import annotations

import argparse
import ast
import importlib.util
import json
import math
import pathlib
import subprocess
import sys
import time

import numpy as np
import pandas as pd

HERE = pathlib.Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import arm_base  # noqa: E402
import declarations  # noqa: E402
import gate  # noqa: E402
import pins  # noqa: E402
import records  # noqa: E402
import verdict as verdict_mod  # noqa: E402  (label_status: one classification for the JSON and the tables)
from criteria import CRITERIA  # noqa: E402

METHOD_VERSION = "broadtype-screen-2026-09-05"
ARMS_DIR = HERE / "arms"
BASE_PINS = verdict_mod.BASE_PINS


# --------------------------------------------------------------------------- paths
def _suffix(scope: str) -> str:
    return "" if scope == "frame" else f"__{scope}"


def verdict_json_path(arm: str, donor: str, scope: str) -> pathlib.Path:
    return pins.OUT / "verdict" / arm / f"donor_id={donor}{_suffix(scope)}.json"


def state_dir(arm: str, donor: str, scope: str) -> pathlib.Path:
    return pins.OUT / "verdict" / arm / f"donor_id={donor}{_suffix(scope)}"


def agreement_path(donor: str, scope: str) -> pathlib.Path:
    return pins.OUT / "verdict" / f"agreement_{donor}{_suffix(scope)}.csv"


def log_dir() -> pathlib.Path:
    return pins.OUT / "logs" / "screen"


def calls_path(arm: str, donor: str, scope: str, config_id: str) -> pathlib.Path:
    """The parquet ``records.write_calls`` produces for this configuration (``:`` -> ``_``)."""
    return records.out_dir(arm, donor, scope) / f"calls__{config_id.replace(':', '_')}.parquet"


# --------------------------------------------------------------------------- arm modules
def arm_module_path(arm: str, arms_dir: pathlib.Path | str = ARMS_DIR) -> pathlib.Path:
    return pathlib.Path(arms_dir) / f"{arm.replace('-', '_')}.py"


def _source_constants(path: pathlib.Path) -> dict:
    """``ARM`` / ``STOCHASTIC`` read from the module's source when importing it is not possible
    in this interpreter (an arm still being written, or one importing a foreign environment)."""
    out = {}
    for node in ast.parse(path.read_text()).body:
        if isinstance(node, ast.Assign) and len(node.targets) == 1 and isinstance(node.targets[0], ast.Name):
            if node.targets[0].id in ("ARM", "STOCHASTIC"):
                try:
                    out[node.targets[0].id] = ast.literal_eval(node.value)
                except (ValueError, SyntaxError):
                    pass
    return out


def inspect_arm(arm: str, arms_dir: pathlib.Path | str = ARMS_DIR) -> dict:
    """Import the arm module by file path and read ``ARM`` and ``STOCHASTIC``.

    The module's own ``ARM`` must equal the name asked for: a renamed configuration is a new arm,
    and running module X under name Y would file X's calls under Y's verdict.
    """
    path = arm_module_path(arm, arms_dir)
    info = {"arm": arm, "module": str(path), "present": path.exists(), "stochastic": False,
            "import_error": None, "declared_arm": None}
    if not path.exists():
        return info
    try:
        spec = importlib.util.spec_from_file_location(f"broadtype_arm__{path.stem}", path)
        mod = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(mod)
        info["stochastic"] = bool(getattr(mod, "STOCHASTIC", False))
        info["declared_arm"] = getattr(mod, "ARM", None)
    except BaseException as e:  # noqa: BLE001  (a SystemExit at import is as informative as an ImportError)
        info["import_error"] = f"{type(e).__name__}: {e}"
        consts = _source_constants(path)
        info["stochastic"] = bool(consts.get("STOCHASTIC", False))
        info["declared_arm"] = consts.get("ARM")
    if info["declared_arm"] is not None and str(info["declared_arm"]) != arm:
        raise ValueError(f"{path} declares ARM={info['declared_arm']!r} but was asked for as {arm!r}; "
                         "a renamed configuration is a new arm")
    return info


# --------------------------------------------------------------------------- the plan
def ablation_markers(donor: str, dec: dict) -> list:
    """One ablation per declared marker of every held-out-gated label acquired on the donor."""
    have = set(arm_base.identity_markers_for(donor))
    out = []
    for spec in dec["labels"].values():
        if not spec["held_out_gated"]:
            continue
        out.extend(m for m in spec["markers"] if m in have and m not in out)
    return out


def replicate_configs(stochastic: bool) -> list:
    """(config_id, extra argv) for the 5 replicates: seeds for a stochastic arm, 90 % hash
    subsamples (distinct by replicate name, which salts the hash) for a deterministic one."""
    n, seed0 = CRITERIA["n_replicates"], CRITERIA["primary_seed"]
    if stochastic:
        return [(f"seed:{i}", ["--seed", str(i), "--replicate", f"seed:{i}"]) for i in range(1, n + 1)]
    frac = CRITERIA["subsample_frac"]
    return [(f"subsample:{i}", ["--seed", str(seed0), "--subsample-frac", repr(frac), "--replicate", f"subsample:{i}"])
            for i in range(1, n + 1)]


def config_plan(donor: str, stochastic: bool, dec: dict) -> list:
    seed0 = str(CRITERIA["primary_seed"])
    plan = [("primary", ["--seed", seed0])]
    plan += [(f"ablate:{m}", ["--seed", seed0, "--ablate", m]) for m in ablation_markers(donor, dec)]
    plan += replicate_configs(stochastic)
    return plan


# --------------------------------------------------------------------------- running
def run_config(python: str, module: pathlib.Path, arm: str, donor: str, scope: str, config_id: str,
               extra: list, *, skip_existing: bool) -> dict:
    """One subprocess of the arm's CLI; wall time measured here.  A fresh run deletes the previous
    calls parquet first so a failed run can never be masked by a stale file."""
    path = calls_path(arm, donor, scope, config_id)
    if skip_existing and path.exists():
        return {"config_id": config_id, "calls": str(path), "reused": True, "wall_seconds": None,
                "returncode": None, "log": None, "ok": True}
    if path.exists():
        path.unlink()
    ld = log_dir()
    ld.mkdir(parents=True, exist_ok=True)
    log = ld / f"{arm}__{donor}__{config_id.replace(':', '_')}.log"
    cmd = [python, str(module), "--donor", str(donor), "--scope", scope, *extra]
    t0 = time.time()
    with log.open("w") as fh:
        fh.write("$ " + " ".join(cmd) + "\n")
        fh.flush()
        rc = subprocess.run(cmd, cwd=str(pins.PKG), stdout=fh, stderr=subprocess.STDOUT).returncode
    wall = time.time() - t0
    ok = rc == 0 and path.exists()
    out = {"config_id": config_id, "calls": str(path) if path.exists() else None, "reused": False,
           "wall_seconds": wall, "returncode": rc, "log": str(log), "ok": ok, "command": cmd}
    if rc == 0 and not path.exists():
        out["note"] = "exit 0 but no calls parquet at the expected path"
    return out


def load_records(arm: str, donor: str, scope: str) -> dict:
    """config_id -> the latest ArmRecord written for it (records.jsonl is append-only)."""
    p = records.out_dir(arm, donor, scope) / "records.jsonl"
    out = {}
    if p.exists():
        for line in p.read_text().splitlines():
            if line.strip():
                r = json.loads(line)
                out[str(r["config_id"])] = r
    return out


def false_award_rate(recs: dict) -> dict | None:
    """The noise false-award rate an arm's tests measured, if its records carry one (arms B, F)."""
    found = {}
    for r in recs.values():
        for src in ("params", "pinned"):
            for k, v in (r.get(src) or {}).items():
                if "false_award" in str(k).lower():
                    found[str(k)] = v
    return found or None


# --------------------------------------------------------------------------- alignment
def aligned_labels(path: pathlib.Path | str, ids: np.ndarray) -> tuple:
    """Labels aligned to the frozen cell order; a cell absent from the file (a subsample
    replicate) is Unavailable for that configuration, so the gate's ``judged`` mask excludes it."""
    df = pd.read_parquet(path, columns=["object_id", "label"])
    s = pd.Series(df["label"].astype(str).to_numpy(), index=pd.Index(df["object_id"].astype(str)))
    if s.index.has_duplicates:
        raise ValueError(f"{path}: duplicate object_id")
    lab = s.reindex(pd.Index(np.asarray(ids).astype(str)))
    n_present = int(lab.notna().sum())
    return lab.fillna(records.UNAVAILABLE).to_numpy().astype(str), n_present


def calls_with_state(primary_path: pathlib.Path | str, donor: str, ids: np.ndarray) -> pd.DataFrame:
    """The primary calls in frozen order with the state tier appended AFTER identity."""
    calls = pd.read_parquet(primary_path)
    calls["object_id"] = calls["object_id"].astype(str)
    calls = calls.set_index("object_id").reindex(pd.Index(np.asarray(ids).astype(str), name="object_id")).reset_index()
    miss = calls["label"].isna()
    calls.loc[miss, "label"] = records.UNAVAILABLE
    calls.loc[miss, "reason"] = "arm_no_output"
    calls["reason"] = calls["reason"].fillna("")
    state = arm_base.state_columns(donor, np.asarray(ids))
    return calls.merge(state, on="object_id", how="left")


def _jsonable(o):
    """numpy scalars -> Python; NaN/inf -> null (json would otherwise emit non-standard NaN)."""
    if isinstance(o, dict):
        return {str(k): _jsonable(v) for k, v in o.items()}
    if isinstance(o, (list, tuple)):
        return [_jsonable(v) for v in o]
    if isinstance(o, np.ndarray):
        return _jsonable(o.tolist())
    if isinstance(o, (bool, np.bool_)):
        return bool(o)
    if isinstance(o, (int, np.integer)):
        return int(o)
    if isinstance(o, (float, np.floating)):
        v = float(o)
        return v if math.isfinite(v) else None
    if isinstance(o, pathlib.Path):
        return str(o)
    return o


# --------------------------------------------------------------------------- judging
def tails_in_frame_for(donor: str, markers: list) -> dict | None:
    p = pins.OUT / "frame" / "tails_in_frame.csv"
    if not p.exists():
        return None
    t = pd.read_csv(p, dtype={"donor_id": str})
    t = t[(t["donor_id"] == str(donor)) & t["marker"].isin(markers)]
    return {r.marker: {"n_tail_in_frame": int(r.n_tail_in_frame), "judgeable_in_frame": bool(r.judgeable_in_frame)}
            for r in t.itertuples()}


def judge(arm: str, donor: str, scope: str, configs: dict, recs: dict, ctx, dec: dict) -> dict:
    """Every declared label through ``gate.judge_label`` on the arm's calls for this donor."""
    ids = ctx.ids
    values = ctx.X                                              # ORIGINAL aligned values, all identity markers
    xy = ctx.context[["X_centroid", "Y_centroid"]].to_numpy(float)
    groups = ctx.region_idx if scope == "frame" else None      # within-region nulls on the frame
    cell_region = ctx.context["cell_region"].to_numpy() if "cell_region" in ctx.context.columns else None
    have = list(values.columns)
    used = {cid: recs[cid] for cid in configs if cid in recs and configs[cid].get("calls")}
    leak = sorted({m for r in used.values() for m in r.get("identity_markers", [])} & set(dec["state_markers"]))
    out = {"labels": {}, "records_used": list(used.values()),
           "n_external_pins": max([len(set(r.get("pinned") or {}) - set(BASE_PINS)) for r in used.values()] or [0]),
           "noise_false_award_rate": false_award_rate(used),
           "state_marker_in_identity_inputs": leak}
    primary_cfg = configs.get("primary") or {}
    primary_path = primary_cfg.get("calls")
    if leak or not primary_path or not pathlib.Path(primary_path).exists():
        err = (f"state_marker_in_identity_inputs: {leak}" if leak
               else f"primary_not_available (rc={primary_cfg.get('returncode')})")
        for L, spec in dec["labels"].items():
            res = {"label": L, "common_lineage": spec["common_lineage"], "declared_markers": spec["markers"],
                   "held_out_gated": spec["held_out_gated"], "verified": False, "error": err}
            res["status"], res["blocking_reasons"] = verdict_mod.label_status(res)
            out["labels"][L] = res
        out["primary_label_counts"] = {}
        return out

    primary, _ = aligned_labels(primary_path, ids)
    ablations, replicates, replicate_n = {}, {}, {}
    for cid, cfg in configs.items():
        if cid == "primary" or not cfg.get("calls") or not pathlib.Path(cfg["calls"]).exists():
            continue
        lab, n = aligned_labels(cfg["calls"], ids)
        if cid.startswith("ablate:"):
            ablations[cid.split(":", 1)[1]] = lab
        else:
            replicates[cid], replicate_n[cid] = lab, n
    # Step-1 tails on the ORIGINAL values, once per acquired identity marker
    tails = {m: gate.tail_mask(values[m].to_numpy(float)) for m in have}
    for L, spec in dec["labels"].items():
        declared_here = [m for m in spec["markers"] if m in have]
        others = {m: tails[m] for L2, s2 in dec["labels"].items() if L2 != L for m in s2["markers"] if m in tails}
        res = gate.judge_label(L, primary, ablations, replicates, values=values, xy=xy, groups=groups, ids=ids,
                               declared_markers=declared_here, held_out_gated=bool(spec["held_out_gated"]),
                               other_tails=others, cell_region=cell_region, seed=CRITERIA["primary_seed"])
        res["common_lineage"] = spec["common_lineage"]
        res["markers_not_acquired"] = [m for m in spec["markers"] if m not in have]
        res["n_primary_calls"] = int((primary == L).sum())
        for pm in res["heldout"]["per_marker"]:
            cfg = configs.get(f"ablate:{pm['marker']}")
            if pm.get("reason") == "ablation_not_run" and cfg is not None and cfg.get("returncode") not in (None, 0):
                pm["reason"] = f"ablation_failed(rc={cfg['returncode']})"
        for r in res["stability"]["per_replicate"]:
            r["n_in_replicate"] = replicate_n.get(r["replicate"])
        res["status"], res["blocking_reasons"] = verdict_mod.label_status(res)
        out["labels"][L] = res
    counts = pd.Series(primary).value_counts()
    out["primary_label_counts"] = {str(k): int(v) for k, v in counts.items()}
    return out


# --------------------------------------------------------------------------- per arm x donor
def screen_arm_donor(info: dict, donor: str, scope: str, dec: dict, *, python: str, skip_existing: bool) -> dict:
    arm = info["arm"]
    module = pathlib.Path(info["module"])
    t_start = time.time()
    plan = config_plan(donor, info["stochastic"], dec)
    configs = {}
    for cid, extra in plan:
        r = run_config(python, module, arm, donor, scope, cid, extra, skip_existing=skip_existing)
        configs[cid] = r
        state = "reused" if r["reused"] else f"{r['wall_seconds']:.1f}s rc={r['returncode']}"
        print(f"[screen] {arm} {donor} {cid}: {state}" + ("" if r["ok"] else "  <-- FAILED"), flush=True)
    recs = load_records(arm, donor, scope)
    for cid, r in configs.items():
        rec = recs.get(cid)
        r["arm_seconds"] = rec.get("seconds") if rec else None
        r["record_written_utc"] = rec.get("written_utc") if rec else None
        r["record"] = rec
    ctx = arm_base.load_context(donor, scope)
    judged = judge(arm, donor, scope, configs, recs, ctx, dec)
    n_ver = sum(1 for v in judged["labels"].values() if v.get("verified"))
    # state tier, appended after identity
    sd = state_dir(arm, donor, scope)
    sd.mkdir(parents=True, exist_ok=True)
    cws = None
    if configs["primary"].get("calls") and pathlib.Path(configs["primary"]["calls"]).exists():
        df = calls_with_state(configs["primary"]["calls"], donor, ctx.ids)
        p = sd / "calls_with_state.parquet"
        df.to_parquet(p, index=False)
        cws = {"path": str(p), "state_columns": [c for c in df.columns if c.startswith("state_")],
               "note": "state tier appended AFTER identity; never an identity input; continuous ranks, nothing thresholded"}
    out = {
        "method_version": METHOD_VERSION,
        "arm": arm, "donor_id": str(donor), "scope": scope, "stochastic": bool(info["stochastic"]),
        "module": info["module"], "module_import_error": info["import_error"],
        "scope_sentence": verdict_mod.scope_sentence(),
        "criteria": dict(CRITERIA), "declarations": dec["method_version"],
        "pins": {"run_id": pins.RUN_ID, "expression_content_id": pins.EXPRESSION_CONTENT_ID},
        "n_cells": int(len(ctx.ids)),
        "n_regions": int(len(np.unique(ctx.region_idx))) if ctx.region_idx is not None else None,
        "identity_markers_on_donor": list(ctx.markers),
        "ablation_markers": [cid.split(":", 1)[1] for cid in configs if cid.startswith("ablate:")],
        "replicate_kind": "seed" if info["stochastic"] else "subsample",
        "configs": configs,
        **judged,
        "n_verified": n_ver,
        "tails_in_frame": tails_in_frame_for(donor, list(ctx.markers)) if scope == "frame" else None,
        "calls_with_state": cws,
        "wall_seconds_total": time.time() - t_start,
        "code": pins.git_head(),
        "written_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
    }
    pj = verdict_json_path(arm, donor, scope)
    pj.parent.mkdir(parents=True, exist_ok=True)
    pj.write_text(json.dumps(_jsonable(out), indent=1, default=str) + "\n")
    summary = ", ".join(f"{L}={v['status']}" for L, v in judged["labels"].items())
    print(f"[screen] {arm} {donor}: {n_ver}/{len(dec['labels'])} verified ({summary}); wrote {pj}", flush=True)
    out["_ids"] = ctx.ids
    return out


def cross_arm_agreement(donor: str, scope: str, ids: np.ndarray) -> pd.DataFrame:
    """ARI / AMI for every pair of arms with a primary on disk.  Descriptive only."""
    found = {}
    for p in sorted((pins.OUT / "arms").glob(f"*/donor_id={donor}/{scope}/calls__primary.parquet")):
        found[p.parents[2].name], _ = aligned_labels(p, ids)
    rows = []
    names = sorted(found)
    for i, a in enumerate(names):
        for b in names[i + 1:]:
            rows.append({"donor_id": donor, "scope": scope, "arm_a": a, "arm_b": b,
                         **gate.agreement(found[a], found[b]), "descriptive_only": True})
    df = pd.DataFrame(rows, columns=["donor_id", "scope", "arm_a", "arm_b", "ari", "ami", "n", "descriptive_only"])
    p = agreement_path(donor, scope)
    p.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(p, index=False)
    print(f"[screen] agreement {donor}: {len(names)} arms with a primary, {len(df)} pairs -> {p}", flush=True)
    return df


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description="broadtype screen: run every arm the same way, judge by the frozen bar")
    ap.add_argument("--donors", nargs="+", default=list(pins.DONORS))
    ap.add_argument("--arms", nargs="+", required=True, help="arm names (module arms/<name with - as _>.py)")
    ap.add_argument("--scope", choices=("frame", "donor"), default="frame")
    ap.add_argument("--skip-existing", action="store_true", help="reuse an existing calls parquet for a config")
    ap.add_argument("--arms-dir", default=str(ARMS_DIR))
    ap.add_argument("--python", default=sys.executable, help="interpreter for the arm CLIs (they spawn their own envs)")
    ap.add_argument("--no-agreement", action="store_true")
    a = ap.parse_args(argv)
    pins.assert_pins()
    dec = declarations.load()
    ids_by_donor = {}
    for arm in a.arms:
        info = inspect_arm(arm, a.arms_dir)
        if not info["present"]:
            print(f"[screen] skip {arm}: module {info['module']} is absent (not written yet?)", flush=True)
            continue
        if info["import_error"]:
            print(f"[screen] {arm}: module import failed in this interpreter ({info['import_error']}); "
                  f"ARM/STOCHASTIC read from source, running its CLI anyway", flush=True)
        for d in a.donors:
            res = screen_arm_donor(info, str(d), a.scope, dec, python=a.python, skip_existing=a.skip_existing)
            ids_by_donor[str(d)] = res["_ids"]
    if not a.no_agreement:
        for d in a.donors:
            ids = ids_by_donor.get(str(d))
            if ids is None:
                ids = arm_base.load_context(str(d), a.scope).ids
            cross_arm_agreement(str(d), a.scope, ids)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
