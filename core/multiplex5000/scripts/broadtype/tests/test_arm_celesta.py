"""Arm D ``celesta`` on planted data: the record shape, the reasons, and recovery of the plant.

No cohort data.  The signature tests need only the declarations; the driver test needs R with
CELESTA installed and skips cleanly without it.
"""
from __future__ import annotations

import json
import pathlib
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
import celesta  # noqa: E402

DEC = declarations.load()
#: One declared marker per label -- small enough for a fast R run, wide enough that every
#: label's negatives list has members present (so the signature carries 0s as well as 1s).
ONE_PER_LABEL = {"T_NK": "CD3e", "B_Plasma": "CD20", "Myeloid": "CD68", "Endocrine": "INS",
                 "Epithelial_nonendocrine": "E_cadherin", "Vascular": "CD31", "Stromal": "SMA"}
MARKERS = list(ONE_PER_LABEL.values())


# --------------------------------------------------------------------------- signature
def test_prior_is_flat_round_one_with_1_0_na_codes():
    prior = celesta.build_prior(MARKERS, DEC)
    assert list(prior.columns[:2]) == ["cell_type", "lineage_level"]
    assert list(prior.columns[2:]) == MARKERS
    assert set(prior["cell_type"]) == set(DEC["labels"])
    assert list(prior["lineage_level"]) == [f"1_0_{k}" for k in range(1, len(DEC["labels"]) + 1)]
    for _, row in prior.iterrows():
        spec = DEC["labels"][row["cell_type"]]
        for m in MARKERS:
            v = row[m]
            if m in spec["markers"]:
                assert v == 1
            elif m in spec["negatives"]:
                assert v == 0
            else:
                assert pd.isna(v)
    # every row has a positive marker; no row is 0/NA only
    assert ((prior[MARKERS] == 1).sum(axis=1) >= 1).all()


def test_prior_uses_every_declared_marker_when_all_are_present():
    prior = celesta.build_prior(DEC["identity_markers"], DEC)
    for _, row in prior.iterrows():
        spec = DEC["labels"][row["cell_type"]]
        assert sorted(m for m in DEC["identity_markers"] if pd.notna(row[m]) and row[m] == 1) == sorted(spec["markers"])
        assert sorted(m for m in DEC["identity_markers"] if pd.notna(row[m]) and row[m] == 0) == sorted(spec["negatives"])


def test_ablation_removes_the_column_and_keeps_the_label():
    markers = [m for m in DEC["identity_markers"] if m != "SMA"]      # Stromal keeps Vimentin
    prior = celesta.build_prior(markers, DEC)
    assert "SMA" not in prior.columns
    assert "Stromal" in set(prior["cell_type"])
    assert celesta.dropped_labels(markers, DEC) == []


def test_label_with_no_surviving_marker_is_dropped_not_zero_only():
    markers = [m for m in MARKERS if m != "SMA"]
    prior = celesta.build_prior(markers, DEC)
    assert "Stromal" not in set(prior["cell_type"])
    assert celesta.dropped_labels(markers, DEC) == ["Stromal"]
    assert list(prior["lineage_level"]) == [f"1_0_{k}" for k in range(1, len(prior) + 1)]


def test_recorded_defaults_are_frequency_free_and_the_prior_is_flagged():
    from criteria import FORBIDDEN_KEY_TOKENS
    for k in celesta.CELESTA_DEFAULTS:
        assert not any(t in k.lower() for t in FORBIDDEN_KEY_TOKENS if t != "prior"), k
    assert "equal.proportions=TRUE" in celesta.DEFAULT_PARAMS["in_method_prior"]
    assert celesta.DEFAULT_PARAMS["filter_cells"] is False


# --------------------------------------------------------------------------- planted context
def _planted(seed: int = 0, per_label: int = 60, n_nan: int = 5, n_small: int = 6):
    """Region 00: seven planted groups, mildly clustered in space, one bright declared marker
    each; a few cells with a NaN measurement; three negative values to censor.  Region 01: too
    small.  Region 02: a second planted region, so the multi-file join is exercised."""
    rng = np.random.default_rng(seed)
    labels = list(ONE_PER_LABEL)
    rows, truth, region, xy = [], [], [], []
    for r, n_per in (("00", per_label), ("02", per_label // 2)):
        centres = rng.uniform(80, 420, size=(len(labels), 2))
        for li, L in enumerate(labels):
            hi = ONE_PER_LABEL[L]
            for _ in range(n_per):
                v = np.exp(rng.normal(np.log(5.0), 0.5, size=len(MARKERS)))
                v[MARKERS.index(hi)] = np.exp(rng.normal(np.log(300.0), 0.3))
                rows.append(v); truth.append(L); region.append(r)
                xy.append(centres[li] + rng.normal(0, 40, size=2))
    X = np.array(rows)
    X[10, 1] = -1.0; X[11, 2] = -2.0; X[12, 3] = -0.5                 # censored at 0, counted (finite rows)
    truth = np.array(truth, dtype=object); region = np.array(region, dtype=object)
    xy = np.array(xy)
    # NaN cells in region 00
    nan_idx = np.arange(n_nan)
    X[nan_idx, 0] = np.nan
    # a too-small region
    Xs = np.exp(rng.normal(np.log(5.0), 0.5, size=(n_small, len(MARKERS))))
    X = np.vstack([X, Xs]); truth = np.concatenate([truth, ["-"] * n_small])
    region = np.concatenate([region, ["01"] * n_small]); xy = np.vstack([xy, rng.uniform(0, 100, size=(n_small, 2)) + 5000])
    ids = np.array([f"cell-{i:05d}" for i in rng.permutation(len(X))], dtype=object)   # ids not in row order
    ctx = arm_base.ArmContext(
        donor_id="synthetic", scope="frame", ids=ids,
        X=pd.DataFrame(X, columns=MARKERS),
        context=pd.DataFrame({"object_id": ids, "X_centroid": xy[:, 0], "Y_centroid": xy[:, 1]}),
        markers=MARKERS, region_idx=np.array([int(r) for r in region]))
    return ctx, truth, region, nan_idx


needs_r = pytest.mark.skipif(not celesta.r_available(), reason="Rscript with CELESTA (+zeallot, spdep) not available")


@needs_r
def test_driver_recovers_the_plant_and_reasons_every_gap(tmp_path):
    ctx, truth, region, nan_idx = _planted()
    calls = celesta.run(ctx, scratch_root=tmp_path / "scratch", diagnostics_dir=tmp_path / "diag")
    # shape and coverage
    assert list(calls.columns) == list(records.CALL_COLUMNS)
    assert calls["object_id"].tolist() == list(ctx.ids)
    records.validate_calls(calls, expected_ids=ctx.ids, labels=set(DEC["labels"]))
    lab = calls["label"].to_numpy(); rsn = calls["reason"].to_numpy(); sc = calls["score"].to_numpy()
    # gaps carry the recorded reasons
    assert (lab[nan_idx] == records.UNAVAILABLE).all() and (rsn[nan_idx] == celesta.REASON_MISSING).all()
    small = region == "01"
    assert (lab[small] == records.UNAVAILABLE).all() and (rsn[small] == celesta.REASON_SMALL).all()
    # called cells have a probability in [0, 1]; Unavailable has NaN
    called = ~np.isin(lab, [records.AMBIGUOUS, records.UNAVAILABLE])
    assert called.any()
    assert np.all((sc[called] >= 0) & (sc[called] <= 1))
    assert np.isnan(sc[lab == records.UNAVAILABLE]).all()
    assert set(rsn[lab == records.AMBIGUOUS]) <= {celesta.REASON_UNASSIGNED}
    # the plant is recovered in both regions: majority per group is the planted label, few wrong
    planted = ~small
    planted[nan_idx] = False
    for r in ("00", "02"):
        for L in ONE_PER_LABEL:
            g = planted & (region == r) & (truth == L)
            got = lab[g]
            hit = (got == L).mean()
            wrong = (~np.isin(got, [L, records.AMBIGUOUS, records.UNAVAILABLE])).mean()
            assert hit >= 0.7, f"region {r} {L}: recovered {hit:.2f}"
            assert wrong <= 0.1, f"region {r} {L}: {wrong:.2f} called as another label"
    # diagnostics beside the calls: censor counts, per-region status, the package defaults it ran
    diag = json.loads((tmp_path / "diag" / "diagnostics__primary.json").read_text())
    assert diag["n_censored_negative_by_marker"] == {MARKERS[1]: 1, MARKERS[2]: 1, MARKERS[3]: 1}
    assert diag["n_missing_measurement"] == len(nan_idx)
    assert diag["regions"]["01"]["status"] == "too_small"
    assert diag["regions"]["00"]["status"] == "ok" and diag["regions"]["02"]["status"] == "ok"
    assert diag["env"]["celesta_sha"] == celesta.CELESTA_COMMIT
    assert diag["env"]["rmixmod_seed"] == ctx.seed == 0
    assert diag["env"]["defaults"]["high_expression_threshold_anchor"] == 0.7
    assert not (tmp_path / "scratch").exists() or not any((tmp_path / "scratch").iterdir()), "scratch not cleaned"


@needs_r
def test_result_is_invariant_to_row_order(tmp_path):
    """FitGmmModel thins zeros by position and Rmixmod's EM is seeded by the arm; sorting by
    object_id plus the supplied seed make one input give one answer."""
    ctx, *_ = _planted(seed=1, per_label=30, n_nan=0, n_small=0)
    # plant exact zeros so the order-dependent thinning path is exercised (> 10 % zeros)
    X = ctx.X.to_numpy().copy()
    rng = np.random.default_rng(7)
    z = rng.random(X.shape) < 0.15
    X[z] = 0.0
    ctx.X = pd.DataFrame(X, columns=MARKERS)
    a = celesta.run(ctx, scratch_root=tmp_path / "s1", diagnostics_dir=tmp_path / "d1")
    perm = rng.permutation(len(ctx.ids))
    ctx2 = arm_base.ArmContext(donor_id=ctx.donor_id, scope=ctx.scope, ids=ctx.ids[perm],
                               X=ctx.X.iloc[perm].reset_index(drop=True),
                               context=ctx.context.iloc[perm].reset_index(drop=True),
                               markers=MARKERS, region_idx=ctx.region_idx[perm])
    b = celesta.run(ctx2, scratch_root=tmp_path / "s2", diagnostics_dir=tmp_path / "d2")
    m = a.merge(b, on="object_id", suffixes=("_a", "_b"))
    assert len(m) == len(a)
    assert (m["label_a"] == m["label_b"]).all()
    assert np.allclose(m["score_a"].fillna(-1), m["score_b"].fillna(-1))


def test_no_marker_left_is_unavailable_for_every_cell(tmp_path):
    ctx, *_ = _planted(per_label=5, n_nan=0, n_small=0)
    ctx.X = ctx.X.iloc[:, :0]
    ctx.markers = []
    calls = celesta.run(ctx, scratch_root=tmp_path / "scratch", diagnostics_dir=tmp_path / "diag")
    assert (calls["label"] == records.UNAVAILABLE).all()
    assert (calls["reason"] == celesta.REASON_ABLATION).all()
    assert len(calls) == len(ctx.ids)
