"""Arm C `tacit` -- TACIT's learned per-type cut, applied to every cell.

    cd core/multiplex5000 && ../../.venv/bin/python scripts/broadtype/arms/tacit.py --donor 6436 --scope frame

Method (frozen in README.md / the protocol, Section 3):

1. Signature = TYPE x MARKER matrix, first column literally ``cell_type``; entry 1 where a declared
   label's identity marker is present in ``ctx.markers`` (the ablated marker is already gone), 0
   otherwise.  A label left with no present marker is dropped from the signature and recorded -- it
   cannot be emitted.  Identity inputs are the declared type markers only.
2. Input = ``ctx.X`` with negative values left-censored at 0 (count recorded).  The package CLR-
   normalises the clustering path internally, and CLR needs non-negative input; the frame's aligned
   values are already non-negative, so on this data the censor is a recorded no-op.
3. Thresholds are learned on a 200 000-cell hash subsample -- ``gate.hash_subsample(ids, 200000,
   salt="tacit")`` over the cells with a COMPLETE identity measurement -- or on all complete cells
   when there are fewer, by ``TACIT(data_expression, r, p, Signature)`` (R 4.5.1, TACIT v1.0.0).
   ``p = 10`` UMAP dimensions; ``r`` is calibrated on the package's own pre-processing against the
   README's "microcluster sizes averaging 0.1-0.5 % of cells" (the exact r and the resulting
   cluster count are recorded in ``learned``; see ``tacit_run.R`` for the ladder rule).
4. The learned cut is then applied to ALL cells through TACIT's own per-type relevance score,
   ``ct = X_censored %*% t(Signature)``, reimplemented in R as the identical expression and compared
   cell for cell against the package's ``Group_threshold_data`` on the subsample (must match).
   A cell positive (strict ``>``) for exactly one type gets that label; several -> Ambiguous
   (multi_type); none -> Ambiguous(no_type_above_threshold).
5. TACIT returns its own per-cell labels for the subsample (microcluster cleaning + kNN resolution
   of mixed cells); those are used for the subsample cells, the thresholds for the rest, and both
   are written to ``tacit_detail__<config>.parquet`` beside the calls.  TACIT's "Others" becomes
   Ambiguous with the reason the score table gives (no type above threshold / multi_type).
6. ``score`` = the winning type's score minus its threshold (a raw-intensity margin, so it is
   comparable within a type, not across types; negative for a TACIT cluster-level label whose own
   cell sits below the cut).  Unavailable -> NaN.

Deviations from the task text, recorded here and in the report:

* The task described the score as computed "after the same CLR normalisation the package applied".
  The installed package does not CLR-normalise ``ct``: ``ct <- data_anb_matrix %*% t(Signature_matrix)``
  is evaluated on the raw ``data_expression`` (CLR feeds only the Seurat clustering object).  "Exactly
  as TACIT computes ct" wins, so the score is on the censored raw values, and the 0-mismatch check
  against the package's own binaries proves it.
* A cell with any NaN identity marker is Unavailable(missing_measurement), never 0: Seurat cannot
  take NaN, and a type whose marker is missing cannot be judged for that cell.  On the 6436 frame
  every such cell has ALL identity markers missing (no aligned measurement exists for it).
* ``p`` is ``min(10, n_markers)`` because uwot refuses more components than features; on the real
  panel (28 or 27 markers) this is always 10.  The clamp exists for the synthetic test.
* TACIT v1.0.0 errors when no microcluster (or no cell) is mixed (a ``1:nrow`` loop over an empty
  frame); that only happens on artificially clean data and is reported as tacit_failed, not patched.

Stochasticity: TACIT fixes ``RunUMAP seed.use = 42`` and ``FindClusters random.seed = 1`` internally;
the only unseeded step is a ``runif`` jitter in its kNN resolution, which ``set.seed(ctx.seed)`` pins.
A ``seed:<n>`` replicate therefore moves only that jitter.  The subsample is hash-drawn (content-keyed,
salt "tacit") and does not depend on the seed.

Wall clock: the R subprocess is capped (``wall_clock_cap_s``); on timeout or any package error every
cell is Unavailable(reason="tacit_failed:<why>") and the reason is printed.

Parameters learned at run time are appended to ``PARAMS["learned"]`` -- ``arm_base.cli`` hands the
same dict object to ``write`` after ``run`` returns, so the record carries them -- and written in full
to ``tacit_learned__<config>.json`` beside the calls.  ``ctx.extra["tacit"]`` may override run
parameters (``r``, ``p``, ``n_subsample``, ``band``, ``r_ladder``, ``wall_clock_cap_s``,
``scratch_root``); the synthetic test uses that to stay small.  A production run passes none.
"""
from __future__ import annotations

import json
import os
import pathlib
import shutil
import subprocess
import sys
import time

import numpy as np
import pandas as pd

HERE = pathlib.Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))
import arm_base  # noqa: E402
import declarations  # noqa: E402
import pins  # noqa: E402
import records  # noqa: E402
from gate import hash_subsample  # noqa: E402

ARM = "tacit"
RSCRIPT = "/usr/local/bin/Rscript"
R_RUNNER = HERE / "tacit_run.R"

#: External pins.  TACIT is installed from GitHub tag v1.0.0 (DESCRIPTION Version 0.0.0.1).
TACIT_PIN = {"package": "TACIT", "remote": "github:huynhkl953/TACIT", "ref": "v1.0.0",
             "sha": "206ae19c888a6fe59ddc069c16da8a96c4403933", "description_version": "0.0.0.1"}

PARAMS: dict = {
    "rscript": RSCRIPT,
    "runner": str(R_RUNNER.relative_to(pins.REPO)),
    "tacit": TACIT_PIN,
    "r_stack": {"R": "4.5.1", "Seurat": "5.5.0", "SeuratObject": "5.4.0", "uwot": "0.2.4",
                "segmented": "2.1-4", "caret": "7.0-1", "dplyr": "1.1.4", "arrow": "24.0.0"},
    "signature": "TYPE x MARKER 0/1; first column `cell_type`; 1 iff the marker is declared for the label "
                 "and present in ctx.markers; labels with no present marker dropped (recorded)",
    "censor_at": 0.0,
    "n_subsample": 200_000,
    "subsample_salt": "tacit",
    "subsample_pool": "cells with a complete identity measurement (any NaN -> Unavailable)",
    "p": 10,
    "r": None,                                   # None -> calibrated (see tacit_run.R); recorded in `learned`
    "r_ladder": [1, 2, 3, 5, 8, 10, 15, 20, 30, 50, 80, 120],
    "r_band_mean_cluster_size_frac": [0.001, 0.005],   # README: microclusters average 0.1-0.5 % of cells
    "r_rule": "first ladder value whose mean microcluster size <= 0.5 % of the subsample, on TACIT's own "
              "pre-processing (deterministic: RunUMAP seed.use=42, FindClusters random.seed=1)",
    "score_expression": "ct = X_censored %*% t(Signature)  (TACIT's `ct`, raw input, no CLR)",
    "positivity": "score > threshold (strict; threshold_function's final comparison)",
    "label_rule": "exactly one positive type -> label; >1 -> Ambiguous(multi_type); 0 -> "
                  "Ambiguous(no_type_above_threshold); subsample cells take TACIT's own label, "
                  "'Others' -> Ambiguous with the score-table reason",
    "score": "winning type's score - its threshold (raw-intensity margin; NaN for Unavailable)",
    "wall_clock_cap_s": 90 * 60,
    "seed_use": "set.seed(ctx.seed) before TACIT (pins its runif jitter only)",
}


# ------------------------------------------------------------------ pieces
def build_signature(markers: list, labels: dict | None = None) -> tuple[pd.DataFrame, dict]:
    """TYPE x MARKER 0/1 with the first column literally ``cell_type``.

    Only markers present in ``markers`` get a column, so an ablation is transparent; a label whose
    declared markers are all absent is dropped and returned in the second element with a reason.
    """
    labels = labels if labels is not None else declarations.load()["labels"]
    present = list(markers)
    rows, dropped = [], {}
    for L, spec in labels.items():
        mine = [m for m in spec["markers"] if m in present]
        if not mine:
            dropped[L] = "no declared identity marker present in the input (ablation or panel)"
            continue
        rows.append({"cell_type": L, **{m: (1 if m in mine else 0) for m in present}})
    sig = pd.DataFrame(rows, columns=["cell_type", *present])
    return sig, dropped


def censor(X: np.ndarray, at: float = 0.0) -> tuple[np.ndarray, int]:
    """Left-censor at ``at``; NaN stays NaN.  Returns the censored copy and the count changed."""
    neg = np.isfinite(X) & (X < at)
    out = np.where(neg, at, X)
    return out, int(neg.sum())


def labels_from_scores(scores: np.ndarray, thresholds: np.ndarray, types: list,
                       tacit_label: np.ndarray | None = None) -> pd.DataFrame:
    """The frozen label rule on a per-type score table (rows = cells, cols = ``types``).

    ``tacit_label`` (object array, None/NaN where absent) overrides a cell's label with TACIT's
    own call: a declared type is kept (score = its margin), "Others" or a "::"-joined mixture becomes
    Ambiguous with the reason the score table supports.  Returns label, score, reason, n_pos,
    label_threshold (the pure threshold rule, kept for the detail record) and tacit_used.
    """
    n = scores.shape[0]
    margin = scores - thresholds[None, :]
    pos = scores > thresholds[None, :]
    n_pos = pos.sum(axis=1)
    best = np.argmax(margin, axis=1) if scores.shape[1] else np.zeros(n, dtype=int)
    best_margin = margin[np.arange(n), best] if scores.shape[1] else np.full(n, np.nan)
    label_thr = np.where(n_pos == 1, np.asarray(types, dtype=object)[np.argmax(pos, axis=1)], records.AMBIGUOUS).astype(object)
    reason_thr = np.where(n_pos == 1, "", np.where(n_pos > 1, "multi_type", "no_type_above_threshold")).astype(object)
    label = label_thr.copy(); reason = reason_thr.copy(); score = best_margin.astype(float).copy()
    used = np.zeros(n, dtype=bool)
    if tacit_label is not None:
        tl = pd.Series(tacit_label, dtype=object)
        has = tl.notna().to_numpy() & (tl.astype(str) != "").to_numpy()
        idx_type = {t: i for i, t in enumerate(types)}
        for i in np.flatnonzero(has):
            t = str(tl.iloc[i])
            used[i] = True
            if t in idx_type:
                label[i] = t; reason[i] = ""; score[i] = margin[i, idx_type[t]]
            elif t == "Others" or "::" in t:
                label[i] = records.AMBIGUOUS
                reason[i] = "multi_type" if n_pos[i] > 1 else ("no_type_above_threshold" if n_pos[i] == 0 else "tacit_others")
            else:
                label[i] = records.AMBIGUOUS; reason[i] = f"tacit_label_undeclared:{t}"
    return pd.DataFrame({"label": label, "score": score, "reason": reason, "n_pos": n_pos,
                         "label_threshold": label_thr, "reason_threshold": reason_thr, "tacit_used": used})


def _unavailable(ids: np.ndarray, reason: str) -> pd.DataFrame:
    return pd.DataFrame({"object_id": pd.Series(ids).astype(str), "label": records.UNAVAILABLE,
                         "score": np.nan, "reason": reason})[list(records.CALL_COLUMNS)]


def _short(msg: str, n: int = 240) -> str:
    msg = " ".join(str(msg).split())
    return msg if len(msg) <= n else msg[:n - 3] + "..."


def r_env_string() -> str:
    try:
        v = subprocess.run([RSCRIPT, "--version"], capture_output=True, text=True, timeout=60)
        line = (v.stdout or v.stderr).strip().splitlines()[0]
    except Exception:  # noqa: BLE001
        line = "Rscript unavailable"
    return f"{RSCRIPT} ({line}); TACIT {TACIT_PIN['ref']}@{TACIT_PIN['sha'][:8]}"


# ------------------------------------------------------------------ the arm
def run(ctx: arm_base.ArmContext) -> pd.DataFrame:
    t_start = time.time()
    over = dict(ctx.extra.get("tacit", {})) if isinstance(ctx.extra, dict) else {}
    n_subsample = int(over.get("n_subsample", PARAMS["n_subsample"]))
    p = int(over.get("p", PARAMS["p"]))
    r = over.get("r", PARAMS["r"])
    band = list(over.get("band", PARAMS["r_band_mean_cluster_size_frac"]))
    ladder = list(over.get("r_ladder", PARAMS["r_ladder"]))
    cap_s = float(over.get("wall_clock_cap_s", PARAMS["wall_clock_cap_s"]))
    scratch_root = pathlib.Path(over.get("scratch_root", pins.OUT / "scratch" / ARM))
    keep_scratch = bool(over.get("keep_scratch", os.environ.get("BROADTYPE_KEEP_SCRATCH", "") == "1"))

    ids = np.asarray(ctx.ids).astype(str)
    learned: dict = {"n_cells": int(len(ids)), "config_id": ctx.config_id, "ablated": ctx.ablated}
    PARAMS["learned"] = learned

    sig, dropped = build_signature(ctx.markers)
    learned["signature_types"] = list(sig["cell_type"]) if len(sig) else []
    learned["labels_dropped_no_marker"] = dropped
    if len(sig) == 0:
        learned["outcome"] = "ablation_removes_required_input"
        return _unavailable(ids, "ablation_removes_required_input")
    types = list(sig["cell_type"])

    X = ctx.X[ctx.markers].to_numpy(dtype=float)
    complete = ~np.isnan(X).any(axis=1)
    Xc, n_censored = censor(X, PARAMS["censor_at"])
    learned["n_censored_values"] = n_censored
    learned["n_incomplete_cells"] = int((~complete).sum())
    if complete.sum() == 0:
        learned["outcome"] = "no_complete_cells"
        return _unavailable(ids, "missing_measurement")

    comp_idx = np.flatnonzero(complete)
    sub_local = hash_subsample(ids[comp_idx], n_subsample, salt=PARAMS["subsample_salt"])
    in_sub = np.zeros(len(ids), dtype=bool)
    in_sub[comp_idx[sub_local]] = True
    learned["n_subsample"] = int(in_sub.sum())
    learned["subsample_is_all_complete_cells"] = bool(in_sub.sum() == len(comp_idx))

    scratch = scratch_root / f"{ctx.donor_id}__{ctx.scope}__{ctx.config_id.replace(':', '_')}__{os.getpid()}"
    scratch.mkdir(parents=True, exist_ok=True)
    try:
        inp = pd.DataFrame(Xc[comp_idx], columns=ctx.markers)
        inp.insert(0, "in_subsample", in_sub[comp_idx])
        inp.insert(0, "object_id", ids[comp_idx])
        inp.to_parquet(scratch / "input.parquet", index=False)
        sig.to_csv(scratch / "signature.csv", index=False)
        (scratch / "params.json").write_text(json.dumps({
            "seed": int(ctx.seed), "p": p, "r": r, "r_ladder": ladder, "band_lo": band[0], "band_hi": band[1]}))

        t0 = time.time()
        try:
            proc = subprocess.run([RSCRIPT, str(R_RUNNER), str(scratch)], capture_output=True, text=True,
                                  timeout=cap_s, cwd=str(scratch))
        except subprocess.TimeoutExpired:
            learned["outcome"] = f"tacit_failed:timeout_after_{int(cap_s)}s"
            learned["r_seconds"] = time.time() - t0
            return _unavailable(ids, learned["outcome"])
        except FileNotFoundError as e:
            learned["outcome"] = f"tacit_failed:rscript_missing:{_short(e)}"
            return _unavailable(ids, learned["outcome"])
        learned["r_seconds"] = time.time() - t0
        (scratch / "driver_stderr.log").write_text(proc.stderr or "")
        if proc.returncode != 0:
            why = ""
            fj = scratch / "failure.json"
            if fj.exists():
                f = json.loads(fj.read_text())
                why = f"{f.get('phase', '?')}:{f.get('error', '')}"
            else:
                why = (proc.stderr or proc.stdout or f"exit {proc.returncode}").strip().splitlines()[-1:] or ["?"]
                why = why[0]
            learned["outcome"] = f"tacit_failed:{_short(why)}"
            return _unavailable(ids, learned["outcome"])

        sc = pd.read_parquet(scratch / "scores.parquet")
        lj = json.loads((scratch / "learned.json").read_text())
        learned.update({k: v for k, v in lj.items()})
        sc["object_id"] = sc["object_id"].astype(str)
        sc = sc.set_index("object_id").reindex(pd.Index(ids[comp_idx]))
        if sc[f"score__{types[0]}"].isna().any():
            learned["outcome"] = "tacit_failed:scores_missing_for_some_cells"
            return _unavailable(ids, learned["outcome"])
        scores = sc[[f"score__{t}" for t in types]].to_numpy(dtype=float)
        thr = np.array([float(lj["thresholds"][t]) for t in types])
        tl = sc["tacit_label"].to_numpy(dtype=object)
        lab = labels_from_scores(scores, thr, types, tacit_label=tl)

        calls = pd.DataFrame({"object_id": ids, "label": records.UNAVAILABLE, "score": np.nan, "reason": "missing_measurement"})
        calls["label"] = calls["label"].astype(object); calls["reason"] = calls["reason"].astype(object)
        calls.loc[comp_idx, "label"] = lab["label"].to_numpy()
        calls.loc[comp_idx, "score"] = lab["score"].to_numpy()
        calls.loc[comp_idx, "reason"] = lab["reason"].to_numpy()

        # record both answers on the subsample and the whole score table beside the calls
        both = lab["tacit_used"].to_numpy()
        agree = (lab["label"].to_numpy() == lab["label_threshold"].to_numpy()) & both
        learned["n_tacit_labelled_cells"] = int(both.sum())
        learned["n_tacit_agrees_with_threshold_rule"] = int(agree.sum())
        learned["n_tacit_label_below_own_threshold"] = int((both & lab["label"].isin(types).to_numpy() & (lab["score"].to_numpy() < 0)).sum())
        learned["label_counts"] = calls["label"].value_counts().to_dict()
        learned["reason_counts"] = calls.loc[calls["reason"] != "", "reason"].value_counts().to_dict()
        learned["seconds_total"] = time.time() - t_start
        learned["outcome"] = "ok"
        detail = pd.DataFrame({"object_id": ids[comp_idx], "in_subsample": in_sub[comp_idx],
                               "cluster": sc["cluster"].to_numpy(), "tacit_label": tl,
                               "label_threshold": lab["label_threshold"].to_numpy(), "reason_threshold": lab["reason_threshold"].to_numpy(),
                               "label": lab["label"].to_numpy(), "score": lab["score"].to_numpy(), "n_pos": lab["n_pos"].to_numpy()})
        for t in types:
            detail[f"score__{t}"] = scores[:, types.index(t)]
        out_dir = pathlib.Path(over["detail_dir"]) if "detail_dir" in over else records.out_dir(ARM, ctx.donor_id, ctx.scope)
        out_dir.mkdir(parents=True, exist_ok=True)
        stem = ctx.config_id.replace(":", "_")
        detail.to_parquet(out_dir / f"tacit_detail__{stem}.parquet", index=False)
        (out_dir / f"tacit_learned__{stem}.json").write_text(json.dumps(learned, indent=1, default=str) + "\n")
        if (scratch / "tacit_stdout.log").exists():
            shutil.copyfile(scratch / "tacit_stdout.log", out_dir / f"tacit_stdout__{stem}.log")
        return calls[list(records.CALL_COLUMNS)]
    finally:
        if not keep_scratch:
            shutil.rmtree(scratch, ignore_errors=True)


if __name__ == "__main__":
    raise SystemExit(arm_base.cli(ARM, run, env=r_env_string(), default_params=PARAMS))
