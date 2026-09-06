"""Arm B on synthetic cells: planted disjoint marker blocks come back as their labels, noise
awards nothing beyond the arm's own recorded false-award rate, and the record has the contract's
shape.  No cohort data is read; the SOM runs for real in ``external/.venv-som``.
"""
import pathlib
import sys

import numpy as np
import pandas as pd
import pytest

HERE = pathlib.Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))
sys.path.insert(0, str(HERE.parent / "arms"))

import arm_base  # noqa: E402
import declarations  # noqa: E402
import som_kintsugi as arm  # noqa: E402
from records import AMBIGUOUS, CALL_COLUMNS, UNAVAILABLE  # noqa: E402

pytestmark = pytest.mark.skipif(not arm.SOM_PYTHON.exists(), reason="external/.venv-som not built")

DEC = declarations.load()
LABELS = list(DEC["labels"])
MARKERS = list(DEC["identity_markers"])
PER_LABEL = 600


def _ctx(X: pd.DataFrame, markers: list, seed: int = 0) -> arm_base.ArmContext:
    ids = np.array([f"cell-{i:06d}" for i in range(len(X))], dtype=object)
    return arm_base.ArmContext(donor_id="synthetic", scope="frame", ids=ids, X=X.reset_index(drop=True),
                               context=pd.DataFrame({"object_id": ids}), markers=list(markers), seed=seed)


def _planted(seed: int = 0) -> tuple[pd.DataFrame, np.ndarray]:
    """Seven populations, each bright in exactly its own label's declared markers."""
    rng = np.random.default_rng(seed)
    blocks, truth = [], []
    for L in LABELS:
        X = rng.normal(0.0, 1.0, size=(PER_LABEL, len(MARKERS)))
        own = [MARKERS.index(m) for m in DEC["labels"][L]["markers"]]
        X[:, own] += 3.0
        blocks.append(X)
        truth += [L] * PER_LABEL
    return pd.DataFrame(np.vstack(blocks), columns=MARKERS), np.array(truth, dtype=object)


def _noise(seed: int = 1) -> pd.DataFrame:
    rng = np.random.default_rng(seed)
    return pd.DataFrame(rng.normal(0.0, 1.0, size=(PER_LABEL * len(LABELS), len(MARKERS))), columns=MARKERS)


def _check_shape(calls: pd.DataFrame, ctx: arm_base.ArmContext) -> None:
    assert list(calls.columns) == list(CALL_COLUMNS)
    assert pd.Index(calls["object_id"].astype(str)).equals(pd.Index(ctx.ids.astype(str)))
    assert set(calls["label"]) <= set(LABELS) | {AMBIGUOUS, UNAVAILABLE}
    flagged = calls["label"].isin([AMBIGUOUS, UNAVAILABLE])
    assert (calls.loc[flagged, "reason"].fillna("") != "").all(), "every Ambiguous/Unavailable needs a reason"
    assert (calls.loc[~flagged, "reason"].fillna("") == "").all()
    assert np.isfinite(calls.loc[~flagged, "score"]).all(), "an awarded cell carries its contrast-minus-p95 score"


@pytest.fixture(scope="module")
def planted_run(tmp_path_factory):
    arm.SCRATCH_ROOT = tmp_path_factory.mktemp("scratch_planted")
    X, truth = _planted()
    # a handful of cells with one missing measurement: unavailable, never zero
    X.iloc[:25, MARKERS.index("CD3e")] = np.nan
    ctx = _ctx(X, MARKERS)
    calls = arm.run(ctx)
    return ctx, calls, truth, dict(ctx.extra["som_kintsugi"])


@pytest.fixture(scope="module")
def noise_run(tmp_path_factory):
    arm.SCRATCH_ROOT = tmp_path_factory.mktemp("scratch_noise")
    ctx = _ctx(_noise(), MARKERS)
    calls = arm.run(ctx)
    return ctx, calls, dict(ctx.extra["som_kintsugi"])


def test_record_shape_and_reasons(planted_run):
    ctx, calls, _truth, _res = planted_run
    _check_shape(calls, ctx)


def test_non_finite_cells_are_unavailable_not_imputed(planted_run):
    _ctx_, calls, _truth, res = planted_run
    head = calls.iloc[:25]
    assert (head["label"] == UNAVAILABLE).all()
    assert (head["reason"] == arm.NON_FINITE_REASON).all()
    assert res["n_non_finite"] == 25 and res["non_finite_by_marker"] == {"CD3e": 25}


def test_planted_blocks_recover_their_labels(planted_run):
    _ctx_, calls, truth, res = planted_run
    ok = calls["reason"] != arm.NON_FINITE_REASON
    for L in LABELS:
        sel = ok.to_numpy() & (truth == L)
        got = calls.loc[sel, "label"]
        frac = float((got == L).mean())
        assert frac >= 0.9, f"{L}: {frac:.2f} recovered; labels seen {got.value_counts().to_dict()}"
    # nothing planted was called as another lineage
    wrong = calls.loc[ok.to_numpy() & ~calls["label"].isin([AMBIGUOUS, UNAVAILABLE]), "label"].to_numpy() != truth[ok.to_numpy() & ~calls["label"].isin([AMBIGUOUS, UNAVAILABLE]).to_numpy()]
    assert wrong.mean() < 0.02, f"{wrong.mean():.3f} of awarded cells carry a different planted label"


def test_node_choice_is_by_seed_stability_and_recorded(planted_run):
    _ctx_, _calls, _truth, res = planted_run
    assert res["chosen_nodes"] in arm.NODES_GRID
    table = {int(k): v["mean_pairwise_ari"] for k, v in res["ari_by_nodes"].items()}
    assert set(table) == set(arm.NODES_GRID)
    assert all(len(v["pairs"]) == 10 for v in res["ari_by_nodes"].values())     # C(5, 2)
    assert res["chosen_nodes"] == arm.choose_nodes(table)
    assert len(res["fits"]) == len(arm.NODES_GRID) * arm.N_SEEDS
    assert arm.DEFAULT_PARAMS["result"] is not None and "false_award" in res


def test_choose_nodes_breaks_ties_toward_fewest():
    assert arm.choose_nodes({9: 0.5, 16: 0.9, 25: 0.9, 36: 0.1, 49: 0.2}) == 16
    assert arm.choose_nodes({9: 1.0, 16: 1.0}) == 9


def test_ablating_a_label_s_last_marker_excludes_it_from_declared():
    markers = [m for m in MARKERS if m not in DEC["labels"]["Stromal"]["markers"]]
    declared, excluded = arm.declared_columns(markers)
    assert "Stromal" not in declared and "Stromal" in excluded
    assert set(declared) == set(LABELS) - {"Stromal"}
    for L, cols in declared.items():
        assert [markers[j] for j in cols] == [m for m in DEC["labels"][L]["markers"] if m in markers]


def _pooled_award(fits: list, nodes: int) -> tuple[float, int, int]:
    """awarded / eligible pooled over the real fits at one node count (5 seeds)."""
    sel = [f["award"] for f in fits if f["nodes"] == nodes]
    n_aw, n_el = sum(a["n_awarded"] for a in sel), sum(a["n_eligible"] for a in sel)
    return n_aw / n_el, n_aw, n_el


def test_noise_awards_nothing_above_the_recorded_false_award_rate(noise_run):
    """The headline rate is the ARM's own false-award rate: the SOM re-fitted on a per-column
    shuffled copy of the block over the same nodes x seeds grid and annotated with the same gate,
    pooled over the 5 seeds at the chosen node count.  Noise-only input, clustered and annotated
    for real over its own 5 seeds at that node count, must not award more of its clusters than
    that (an allowance of 0.15 for two pooled finite samples of the same rate; a single 16-node
    fit estimates the rate to only +/- 0.12, which is why both sides are pooled).  The
    annotation-only rate (shuffled features, primary clusters) is recorded beside it and sits
    BELOW the arm's -- a SOM cluster is selected to be co-extreme in the features the gate then
    tests, which the gate's random-set null does not model -- so it is asserted to exist, not to
    bound anything."""
    ctx, calls, res = noise_run
    _check_shape(calls, ctx)
    fa = res["false_award"]
    recorded = fa["false_award_rate"]
    assert np.isfinite(recorded) and np.isfinite(fa["annotation_only"]["rate"])
    assert set(fa["arm_by_nodes"]) == {str(n) for n in arm.NODES_GRID}
    assert all(len(v["per_seed"]) == arm.N_SEEDS for v in fa["arm_by_nodes"].values())
    assert recorded == fa["arm_by_nodes"][str(res["chosen_nodes"])]["rate"]
    rate_noise, n_aw, n_el = _pooled_award(res["fits"], res["chosen_nodes"])
    assert rate_noise <= recorded + 0.15, (
        f"noise-only input awarded {n_aw}/{n_el} clusters ({rate_noise:.2f}) over 5 seeds at "
        f"{res['chosen_nodes']} nodes against a recorded arm false-award rate of {recorded:.2f}")


def test_planted_structure_awards_far_more_than_its_structureless_twin(planted_run):
    """On planted data the shuffled twin must lose the structure: the arm's recorded false-award
    rate stays below the primary award rate and far fewer cells receive a lineage."""
    _ctx_, _calls, _truth, res = planted_run
    prim = res["primary_award"]
    fa = res["false_award"]["arm_by_nodes"][str(res["chosen_nodes"])]
    assert prim["rate"] >= 0.9
    assert fa["rate"] < prim["rate"]
    assert fa["frac_cells_awarded"] < prim["frac_cells_awarded"]
