"""Arm C `tacit` on planted structure -- no cohort data.

The load-bearing test drives the real R pipeline (Rscript + TACIT v1.0.0) through ``tacit.run`` on
a 3 000-cell matrix built from six DECLARED identity markers -- Endocrine (GCG, INS), T_NK (CD3e,
CD4), Vascular (CD31, Caveolin) -- with planted doublets, planted no-signal cells, a few NaN rows and
a few negative values.  Doublets and no-signal cells are not decoration: TACIT v1.0.0 errors when
no microcluster is mixed (an empty-frame ``1:nrow`` loop), so a perfectly clean plant would test the
package's bug and nothing else.  The subsample (2 000 of 3 000) is smaller than the matrix on
purpose, so both label paths run: TACIT's own labels inside the subsample, the learned thresholds
outside it.

Everything that needs R skips with a stated reason when Rscript or TACIT is missing; the label rule,
the signature builder and the censor are tested without R.
"""
from __future__ import annotations

import json
import pathlib
import shutil
import subprocess
import sys

import numpy as np
import pandas as pd
import pytest

HERE = pathlib.Path(__file__).resolve().parents[1]
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(HERE / "arms"))

import arm_base  # noqa: E402
import declarations  # noqa: E402
import records  # noqa: E402
import tacit  # noqa: E402

DECLARED = declarations.load()["labels"]
MARKERS = ["GCG", "INS", "CD3e", "CD4", "CD31", "Caveolin"]
PLANT = {"Endocrine": ["GCG", "INS"], "T_NK": ["CD3e", "CD4"], "Vascular": ["CD31", "Caveolin"]}


def _r_available() -> str | None:
    """None when R + TACIT + Seurat + arrow are usable, else the reason to skip."""
    if not pathlib.Path(tacit.RSCRIPT).exists():
        return f"{tacit.RSCRIPT} not found"
    try:
        p = subprocess.run([tacit.RSCRIPT, "-e",
                            "suppressPackageStartupMessages({library(TACIT); library(Seurat); library(arrow); library(segmented)}); cat('ok')"],
                           capture_output=True, text=True, timeout=300)
    except Exception as e:  # noqa: BLE001
        return f"Rscript probe failed: {e}"
    if p.returncode != 0 or "ok" not in p.stdout:
        return f"TACIT/Seurat/arrow not loadable in R: {(p.stderr or '').strip().splitlines()[-1:]}"
    return None


def _planted(n: int = 3000, seed: int = 0):
    """Three pure types, 5 % Endocrine+T_NK doublets, 5 % no-signal, 10 NaN rows, 20 negatives."""
    rng = np.random.default_rng(seed)
    truth = rng.choice(["Endocrine", "T_NK", "Vascular", "doublet", "none"], size=n, p=[.3, .3, .3, .05, .05])
    X = np.abs(rng.normal(50, 15, size=(n, len(MARKERS))))
    col = {m: i for i, m in enumerate(MARKERS)}
    for L, ms in PLANT.items():
        sel = (truth == L) | ((truth == "doublet") & (L in ("Endocrine", "T_NK")))
        for m in ms:
            X[sel, col[m]] += 400
    nan_rows = rng.choice(n, 10, replace=False)
    X[nan_rows, :] = np.nan
    ok_rows = np.setdiff1d(np.arange(n), nan_rows)
    neg = (rng.choice(ok_rows, 20, replace=False), rng.integers(0, len(MARKERS), 20))
    X[neg] = -5.0
    ids = np.array([f"cell{i:05d}" for i in range(n)])
    ctx = arm_base.ArmContext(donor_id="synthetic", scope="frame", ids=ids,
                              X=pd.DataFrame(X, columns=MARKERS), context=pd.DataFrame({"object_id": ids}),
                              markers=list(MARKERS), seed=0)
    return ctx, truth, nan_rows


# ----------------------------------------------------------------- without R
def test_signature_shape_and_dropped_labels():
    sig, dropped = tacit.build_signature(MARKERS)
    assert list(sig.columns) == ["cell_type", *MARKERS]
    assert set(sig["cell_type"]) == set(PLANT)
    assert set(dropped) == set(DECLARED) - set(PLANT)
    vals = sig[MARKERS].to_numpy()
    assert set(np.unique(vals)) <= {0, 1}
    for _, row in sig.iterrows():
        assert [m for m in MARKERS if row[m] == 1] == PLANT[row["cell_type"]]


def test_signature_ablation_is_transparent():
    sig, dropped = tacit.build_signature([m for m in MARKERS if m != "INS"])
    assert "INS" not in sig.columns
    assert sig.set_index("cell_type").loc["Endocrine"].sum() == 1
    sig2, dropped2 = tacit.build_signature(["GCG", "INS"])
    assert list(sig2["cell_type"]) == ["Endocrine"] and set(dropped2) == set(DECLARED) - {"Endocrine"}


def test_censor_counts_only_negatives_and_keeps_nan():
    X = np.array([[1.0, -2.0], [np.nan, 0.0], [-0.5, 3.0]])
    Xc, n = tacit.censor(X)
    assert n == 2
    assert np.isnan(Xc[1, 0]) and Xc[0, 1] == 0.0 and Xc[2, 0] == 0.0 and Xc[0, 0] == 1.0


def test_label_rule_on_a_hand_table():
    types = ["A", "B"]
    thr = np.array([10.0, 20.0])
    scores = np.array([[15.0, 5.0],    # A only
                       [15.0, 25.0],   # both
                       [5.0, 5.0],     # none
                       [10.0, 21.0]])  # A exactly at threshold (not >), B above -> B
    lab = tacit.labels_from_scores(scores, thr, types)
    assert list(lab["label"]) == ["A", records.AMBIGUOUS, records.AMBIGUOUS, "B"]
    assert list(lab["reason"]) == ["", "multi_type", "no_type_above_threshold", ""]
    assert lab["score"].iloc[0] == 5.0 and lab["score"].iloc[3] == 1.0
    assert lab["score"].iloc[2] == -5.0            # the leading candidate's (negative) margin
    # TACIT's own labels override: a type is kept with its margin, Others -> Ambiguous with the table's reason
    tl = np.array(["B", "Others", "Others", None], dtype=object)
    lab2 = tacit.labels_from_scores(scores, thr, types, tacit_label=tl)
    assert list(lab2["label"]) == ["B", records.AMBIGUOUS, records.AMBIGUOUS, "B"]
    assert lab2["score"].iloc[0] == -15.0 and list(lab2["reason"][:3]) == ["", "multi_type", "no_type_above_threshold"]
    assert list(lab2["tacit_used"]) == [True, True, True, False]


# ----------------------------------------------------------------- with R
@pytest.mark.skipif(_r_available() is not None, reason=str(_r_available()))
def test_planted_types_recovered_through_the_r_pipeline(tmp_path):
    ctx, truth, nan_rows = _planted()
    ctx.extra["tacit"] = {"r": 2, "p": 5, "n_subsample": 2000, "wall_clock_cap_s": 900,
                          "scratch_root": tmp_path / "scratch", "detail_dir": tmp_path / "detail"}
    calls = tacit.run(ctx)
    learned = tacit.PARAMS["learned"]
    assert learned["outcome"] == "ok", learned

    # record shape: every cell, in order, with the common columns and a reason for every non-call
    assert list(calls.columns) == list(records.CALL_COLUMNS)
    assert list(calls["object_id"]) == list(ctx.ids)
    records.validate_calls(calls, expected_ids=ctx.ids, labels=set(DECLARED))
    assert set(calls["label"]) <= set(PLANT) | {records.AMBIGUOUS, records.UNAVAILABLE}

    # the arm's own bookkeeping
    assert learned["r_used"] == 2 and learned["p_used"] == 5 and learned["n_subsample"] == 2000
    assert learned["n_binary_mismatch_vs_tacit"] == 0            # our `ct` == the package's Group_threshold_data
    assert learned["n_censored_values"] == 20
    assert learned["n_incomplete_cells"] == 10
    assert set(learned["labels_dropped_no_marker"]) == set(DECLARED) - set(PLANT)
    assert set(learned["thresholds"]) == set(PLANT) and all(np.isfinite(v) for v in learned["thresholds"].values())

    # missing measurement is Unavailable, never 0
    nan_calls = calls.iloc[nan_rows]
    assert (nan_calls["label"] == records.UNAVAILABLE).all() and (nan_calls["reason"] == "missing_measurement").all()
    assert nan_calls["score"].isna().all()

    # planted pure types recovered on BOTH paths (TACIT's labels inside the subsample, thresholds outside)
    detail = pd.read_parquet(tmp_path / "detail" / "tacit_detail__primary.parquet").set_index("object_id")
    in_sub = detail["in_subsample"].reindex(ctx.ids).fillna(False).to_numpy(bool)
    assert in_sub.sum() == 2000 and (~in_sub).sum() == 1000
    lab = calls["label"].to_numpy()
    complete = ~np.isin(np.arange(len(ctx.ids)), nan_rows)
    for L in PLANT:
        for path, mask in (("tacit", in_sub), ("threshold", ~in_sub)):
            sel = (truth == L) & complete & mask
            frac = (lab[sel] == L).mean()
            # 0.85 on the threshold path: a cut learned on 2,000 cells and applied outside them
            # recovers a planted type at ~0.89 (measured 2026-09-05); the bar records that, it does not tune it
            assert frac >= (0.95 if path == "tacit" else 0.85), f"{L} via {path}: {frac:.3f}"
    called = calls["label"].isin(PLANT).to_numpy()
    assert np.isfinite(calls.loc[called, "score"]).all()

    # The native subsample path must preserve TACIT's labels, including cluster-level calls
    # below the cell's own threshold. The hand-table test above covers that override too.
    native = detail["tacit_label"].reindex(ctx.ids)
    native_type = in_sub & native.isin(PLANT).to_numpy()
    assert native.notna().to_numpy().tolist() == in_sub.tolist()
    assert np.array_equal(lab[native_type], native.to_numpy()[native_type])
    assert (lab[in_sub & ~native_type] == records.AMBIGUOUS).all()
    for L in PLANT:
        sel = native_type & (native.to_numpy() == L)
        margins = detail[f"score__{L}"].reindex(ctx.ids).to_numpy()[sel] - learned["thresholds"][L]
        np.testing.assert_allclose(calls["score"].to_numpy()[sel], margins)

    # The other path uses the learned thresholds for every complete cell outside the subsample.
    threshold_labels = detail["label_threshold"].reindex(ctx.ids).to_numpy()
    threshold_reasons = detail["reason_threshold"].reindex(ctx.ids).to_numpy()
    outside = complete & ~in_sub
    assert np.array_equal(lab[outside], threshold_labels[outside])
    assert np.array_equal(calls["reason"].to_numpy()[outside], threshold_reasons[outside])

    # doublets are Ambiguous(multi_type); no-signal cells are Ambiguous(no_type_above_threshold)
    # the package resolves a multi-positive cell by kNN inside the subsample (its own design); only the
    # threshold path outside the subsample can say Ambiguous(multi_type), so that is where it is asserted
    dbl = (truth == "doublet") & complete & ~in_sub
    assert ((lab[dbl] == records.AMBIGUOUS) & (calls["reason"].to_numpy()[dbl] == "multi_type")).mean() >= 0.75   # measured 0.78
    non = (truth == "none") & complete
    # This is a threshold-path property, not a promise about TACIT's cluster/kNN overrides:
    # the pinned package labels 16 planted no-signal cells below every threshold (2026-09-17).
    # Keep the 0.8 cutoff on the learned decision for ALL no-signal cells and on actual calls
    # outside the subsample; above we check native subsample labels cell for cell.
    assert ((threshold_labels[non] == records.AMBIGUOUS) &
            (threshold_reasons[non] == "no_type_above_threshold")).mean() >= 0.8
    non_outside = non & outside
    assert ((lab[non_outside] == records.AMBIGUOUS) &
            (calls["reason"].to_numpy()[non_outside] == "no_type_above_threshold")).mean() >= 0.8

    # scratch is gone; the detail sidecars exist where we asked
    assert not any((tmp_path / "scratch").glob("*")) or not (tmp_path / "scratch").exists()
    assert (tmp_path / "detail" / "tacit_learned__primary.json").exists()
    lj = json.loads((tmp_path / "detail" / "tacit_learned__primary.json").read_text())
    assert lj["versions"]["TACIT_sha"] == tacit.TACIT_PIN["sha"]


@pytest.mark.skipif(_r_available() is not None, reason=str(_r_available()))
def test_wall_clock_cap_yields_unavailable_for_every_cell(tmp_path):
    ctx, _truth, _nan = _planted(n=600)
    ctx.extra["tacit"] = {"r": 2, "p": 5, "n_subsample": 600, "wall_clock_cap_s": 0.01,
                          "scratch_root": tmp_path / "scratch", "detail_dir": tmp_path / "detail"}
    calls = tacit.run(ctx)
    assert (calls["label"] == records.UNAVAILABLE).all()
    assert calls["reason"].str.startswith("tacit_failed:timeout_after_").all()
    assert calls["score"].isna().all()
    assert tacit.PARAMS["learned"]["outcome"].startswith("tacit_failed:timeout")
    shutil.rmtree(tmp_path / "scratch", ignore_errors=True)
