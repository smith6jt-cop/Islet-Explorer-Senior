"""Seven silent failures, each a null that used to be written as a number or a verdict.

  * build_galectin_layer -- a blank QuPath count became the integer 0 (``fillna(0)``), so a cell
    count that was never exported read the same as a measured zero;
  * build_galectin_layer -- a negative difference between two object-classifier passes was
    clipped to 0 (``clip(lower=0)``), so a disagreement read as "no metabolically-active cells";
  * build_galectin_layer -- a blank ``Num Detections`` left every proportion NaN on a row whose
    state still said ``measured``, with nothing naming the missing denominator;
  * islet_identity_mixture -- a null draw whose mixture fit raised was dropped, and the p-value
    was taken over the survivors with no record of how many were lost;
  * islet_identity_mixture -- a null that lost EVERY draw gave four NaN p-values, and the verdict
    chain fell through them to ``one population``;
  * build_galectin_layer -- a negative, non-finite or fractional count or denominator was not a
    count and was cast, divided by, or raised on, rather than refused where the column is read;
  * islet_identity_mixture -- a control or target whose own fit returned a non-finite statistic
    had NaN p-values, and the verdict chain fell through them to ``one population``.

Every frame here is two synthetic rows: one with the value present, one without.  A row that
HAD a measurement must be numerically identical before and after, and that is asserted against
the old arithmetic written out inline, not assumed.  No biological number appears anywhere;
nothing here is a verification target.

Run with
    cd core/multiplex5000 && PYTHONDONTWRITEBYTECODE=1 .venv/bin/python -m pytest \\
      -p no:cacheprovider -q tests/current/test_silent_failures_galectin_mixture.py
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd

REPO = Path(__file__).resolve().parents[4]
sys.path.insert(0, str(REPO / "core/multiplex5000/scripts/typing"))

import build_galectin_layer as G       # noqa: E402
import islet_identity_mixture as M     # noqa: E402

COUNT_NAMES = ("count_Galectin_positive", "count_Galectin_macrophage",
               "count_Macrophage", "count_Metabolically_active")


def _annotation_rows(**columns) -> pd.DataFrame:
    """Two annotation rows with every count present; ``columns`` overrides one of them."""
    frame = pd.DataFrame({
        "case_id": ["115", "115"],
        "islet_key": ["Islet_1", "Islet_2"],
        "region": ["core", "core"],
        "Num Detections": [10, 10],
        "Num M2Galectin3": [5.0, 5.0],
        "Num Iba1_CD68: M2Galectin3": [2.0, 2.0],
        "Num Iba1_CD68": [3.0, 3.0],
    })
    for name, values in columns.items():
        frame[name] = values
    return frame


def _old_counts(frame: pd.DataFrame) -> pd.DataFrame:
    """The arithmetic as it stood before 2026-09-21, for the rows it was right about."""
    out = frame[["case_id", "islet_key", "region"]].copy()
    out["n_cells_qupath"] = frame["Num Detections"]
    for column, label in G.COUNT_COLUMNS.items():
        out[f"count_{label}"] = frame[column].fillna(0).astype(int)
    out["count_Metabolically_active"] = (out["count_Galectin_positive"]
                                         - out["count_Galectin_macrophage"]).clip(lower=0)
    for label in list(G.COUNT_COLUMNS.values()) + ["Metabolically_active"]:
        out[f"prop_{label}"] = np.where(out["n_cells_qupath"] > 0,
                                        out[f"count_{label}"] / out["n_cells_qupath"], np.nan)
    return out


# --------------------------------------------------------------------------- site 1: blank count


def test_a_blank_qupath_count_is_unavailable_not_zero():
    frame = _annotation_rows(**{"Num M2Galectin3": [5.0, np.nan]})
    out = G.derive_counts(frame)

    # The row that had every measurement is what it always was.
    assert [int(out[c].iloc[0]) for c in COUNT_NAMES] == [5, 2, 3, 3]
    assert out["prop_Galectin_positive"].iloc[0] == 5 / 10
    assert out["prop_Metabolically_active"].iloc[0] == 3 / 10
    assert out[G.COUNT_STATE].iloc[0] == G.MEASURED
    assert pd.isna(out[G.COUNT_REASON].iloc[0])

    # The blank one is null -- not 0 -- and the row says which count and why.
    assert pd.isna(out["count_Galectin_positive"].iloc[1])
    assert pd.isna(out["count_Metabolically_active"].iloc[1]), "derived from the blank"
    assert pd.isna(out["prop_Galectin_positive"].iloc[1])
    assert pd.isna(out["prop_Metabolically_active"].iloc[1])
    # The two passes that WERE exported on that row are still counts.
    assert int(out["count_Galectin_macrophage"].iloc[1]) == 2
    assert int(out["count_Macrophage"].iloc[1]) == 3
    assert out[G.COUNT_STATE].iloc[1] == G.UNAVAILABLE
    reason = out[G.COUNT_REASON].iloc[1]
    assert "Num M2Galectin3" in reason and "count_Galectin_positive" in reason

    # No column was renamed or removed: the R cards read these by name.  The state and reason
    # columns share NEITHER prefix, because `apps/islet/R/00_globals.R` (startsWith "count_")
    # and `scripts/typing/trajectory_islet.py` (c[len("count_"):]) both glob by it.
    assert {c for c in out.columns if c.startswith("count_")} == set(COUNT_NAMES)
    assert not G.COUNT_STATE.startswith(("count_", "prop_"))
    assert not G.COUNT_REASON.startswith(("count_", "prop_"))
    assert {c for c in out.columns if c.startswith("prop_")} == {
        c.replace("count_", "prop_") for c in COUNT_NAMES}


def test_every_blank_count_is_named_and_the_others_are_kept():
    """A row can lose one count without losing the rest; each null names its own input."""
    frame = _annotation_rows(**{"Num Iba1_CD68": [3.0, np.nan]})
    out = G.derive_counts(frame)
    assert pd.isna(out["count_Macrophage"].iloc[1])
    # The difference does not depend on the macrophage denominator, so it is still a count.
    assert int(out["count_Metabolically_active"].iloc[1]) == 3
    assert int(out["count_Galectin_positive"].iloc[1]) == 5
    assert out[G.COUNT_STATE].iloc[1] == G.UNAVAILABLE
    assert "Num Iba1_CD68" in out[G.COUNT_REASON].iloc[1]
    assert "count_Macrophage" in out[G.COUNT_REASON].iloc[1]


def test_a_blank_qupath_denominator_is_unavailable_and_a_zero_one_stays_measured():
    """`Num Detections` is the denominator of every proportion.  Blank, it is UNAVAILABLE, and
    the row must say so and name it: until 2026-09-21 the four counts were read, the state was
    `measured`, and every proportion was silently NaN.  Exactly 0 with a present count is a
    different thing -- a measured region with no detections -- and stays as it is today: NaN
    proportion, `measured`, no reason.  Whether a count agrees with a zero denominator is not
    something this build checks, and that is unchanged."""
    blank = G.derive_counts(_annotation_rows(**{"Num Detections": [10, np.nan]}))

    # The row with a denominator is what it always was.
    assert blank[G.COUNT_STATE].iloc[0] == G.MEASURED
    assert pd.isna(blank[G.COUNT_REASON].iloc[0])
    assert blank["prop_Galectin_positive"].iloc[0] == 5 / 10

    # Blank: the four counts were read and are kept; no proportion; the state names why.
    assert pd.isna(blank["n_cells_qupath"].iloc[1])
    assert [int(blank[c].iloc[1]) for c in COUNT_NAMES] == [5, 2, 3, 3]
    for c in COUNT_NAMES:
        assert pd.isna(blank[c.replace("count_", "prop_")].iloc[1]), c
    assert blank[G.COUNT_STATE].iloc[1] == G.UNAVAILABLE
    reason = blank[G.COUNT_REASON].iloc[1]
    assert "Num Detections" in reason and "n_cells_qupath" in reason

    # Zero: present, and 0 -- the proportion is undefined but nothing is missing.
    zero = G.derive_counts(_annotation_rows(**{"Num Detections": [10, 0]}))
    assert int(zero["n_cells_qupath"].iloc[1]) == 0
    assert [int(zero[c].iloc[1]) for c in COUNT_NAMES] == [5, 2, 3, 3]
    for c in COUNT_NAMES:
        assert np.isnan(zero[c.replace("count_", "prop_")].iloc[1]), c
    assert (zero[G.COUNT_STATE] == G.MEASURED).all()
    assert zero[G.COUNT_REASON].isna().all()


# ----------------------------------------------------------------- site 2: negative difference


def test_a_classifier_disagreement_is_unavailable_not_clipped_to_zero():
    frame = _annotation_rows(**{"Num M2Galectin3": [5.0, 2.0],
                                "Num Iba1_CD68: M2Galectin3": [2.0, 5.0]})
    out = G.derive_counts(frame)

    # Agreement: the difference is a count, and the row is measured.
    assert int(out["count_Metabolically_active"].iloc[0]) == 3
    assert out["prop_Metabolically_active"].iloc[0] == 3 / 10
    assert out[G.COUNT_STATE].iloc[0] == G.MEASURED

    # Disagreement: a negative difference is not a count and must not become 0.
    assert pd.isna(out["count_Metabolically_active"].iloc[1])
    assert pd.isna(out["prop_Metabolically_active"].iloc[1])
    # The two inputs were read and are kept as read; only the derived count is withheld.
    assert int(out["count_Galectin_positive"].iloc[1]) == 2
    assert int(out["count_Galectin_macrophage"].iloc[1]) == 5
    assert out[G.COUNT_STATE].iloc[1] == G.UNAVAILABLE
    reason = out[G.COUNT_REASON].iloc[1]
    assert "count_Metabolically_active" in reason
    assert "Num M2Galectin3" in reason and "Num Iba1_CD68: M2Galectin3" in reason
    assert "(5)" in reason and "(2)" in reason, "the reason carries both numbers"


# ------------------------------------------------------ both sites: measured rows are unchanged


def test_rows_with_every_measurement_are_numerically_identical_to_the_old_arithmetic():
    frame = pd.DataFrame({
        "case_id": ["115", "115", "6442"],
        "islet_key": ["Islet_1", "Islet_2", "Islet_7"],
        "region": ["core", "peri", "core"],
        "Num Detections": [10, 7, 0],
        "Num M2Galectin3": [5.0, 0.0, 4.0],
        "Num Iba1_CD68: M2Galectin3": [2.0, 0.0, 4.0],
        "Num Iba1_CD68": [3.0, 1.0, 6.0],
    })
    new, old = G.derive_counts(frame), _old_counts(frame)
    for c in COUNT_NAMES:
        assert new[c].astype("int64").tolist() == old[c].astype("int64").tolist(), c
        p = c.replace("count_", "prop_")
        assert new[p].dtype == old[p].dtype == np.dtype("float64"), p
        assert np.array_equal(new[p].to_numpy(), old[p].to_numpy(), equal_nan=True), p
    assert (new[G.COUNT_STATE] == G.MEASURED).all()
    assert new[G.COUNT_REASON].isna().all()
    # A zero denominator still gives no proportion, exactly as before.
    assert np.isnan(new["prop_Galectin_positive"].iloc[2])


# --------------------------------------------------------------- site 3: dropped null draws


def _counting_stub(raise_on_call: int | None):
    """A stand-in for `fit_statistics` that returns the seed it was given, so each draw's
    statistic is known in advance, and raises on exactly one call when asked."""
    calls = {"n": 0}

    def fit(x: np.ndarray, seed: int) -> tuple[float, float]:
        calls["n"] += 1
        if calls["n"] == raise_on_call:
            raise ValueError("the mixture fit did not converge")
        return float(seed), float(seed) / 10
    return fit


def test_a_dropped_null_draw_is_counted_and_the_pvalue_is_still_over_the_survivors():
    pool = np.arange(40, dtype=float)[:, None]
    # Draw b is fitted with seed + b + 1, so five draws from seed 0 carry LR 1, 2, 3, 4, 5.
    lost_one = M.reference_null(pool, n=5, draws=5, seed=0, fit=_counting_stub(raise_on_call=2))
    intact = M.reference_null(pool, n=5, draws=5, seed=0, fit=_counting_stub(raise_on_call=None))

    assert lost_one.n_draws_requested == 5
    assert lost_one.n_dropped_draws == 1
    assert intact.n_draws_requested == 5 and intact.n_dropped_draws == 0

    # The survivors are exactly the draws whose fit returned, unchanged in value ...
    assert lost_one.values.shape == (4, 2)
    assert lost_one.values[:, 0].tolist() == [1.0, 3.0, 4.0, 5.0]
    assert np.array_equal(lost_one.values, intact.values[[0, 2, 3, 4]])

    # ... and the p-value is the old computation over those survivors, untouched.
    observed = 3.5
    survivors = lost_one.values[:, 0]
    expected = (np.sum(survivors >= observed) + 1) / (len(survivors) + 1)
    assert M.pvalue(survivors, observed) == expected == 3 / 5


def test_the_draw_counts_sit_beside_every_pvalue_the_module_writes(monkeypatch):
    """Through `test_population`, the function that writes the rows: the joint null loses one
    draw, the per-feature nulls lose none, and each row says so beside its p-values."""
    rng = np.random.default_rng(0)
    space = ["CD99", "nucleus_elongation"]
    ref = pd.DataFrame(rng.normal(size=(600, 2)), columns=space)
    und = pd.DataFrame(rng.normal(size=(300, 2)), columns=space)

    # Calls 1..4 are the joint null's draws (the third raises); 5 and 6 are control and target.
    # Returning the call number makes control and target sit ABOVE every null draw, so the
    # control clears and the verdict is the plain one -- the point here is the bookkeeping.
    calls = {"n": 0}

    def fit(x: np.ndarray, seed: int) -> tuple[float, float]:
        calls["n"] += 1
        if calls["n"] == 3:
            raise ValueError("the mixture fit did not converge")
        return 100.0 + calls["n"], 100.0 + calls["n"]
    monkeypatch.setattr(M, "fit_statistics", fit)

    row, feats, _, _ = M.test_population(ref, und, space, "115", "large",
                                         draws=4, draws_1d=4, seed=0)

    assert row["n_draws_requested"] == 4
    assert row["n_dropped_draws"] == 1
    # Three survivors (calls 1, 2, 4), control above all of them: (0 + 1) / (3 + 1).
    assert row["p_lr_control"] == row["p_sep_control"] == 1 / 4
    assert row["p_lr_target"] == row["p_sep_target"] == 1 / 4
    assert row["verdict"] == "one population"

    assert len(feats) == len(space)
    for feat in feats:
        assert feat["n_draws_requested"] == 4
        assert feat["n_dropped_draws"] == 0
        assert feat["p_lr_control"] == 1 / 5, "four survivors, nothing dropped"



# ------------------------------------------------------------------- site 4: an empty null


def _raising_stub():
    """A stand-in for `fit_statistics` that raises on every call: a null with no survivors."""
    def fit(x: np.ndarray, seed: int) -> tuple[float, float]:
        raise ValueError("the mixture fit did not converge")
    return fit


def test_a_null_that_lost_every_draw_is_empty_and_its_pvalue_is_unavailable():
    pool = np.arange(40, dtype=float)[:, None]
    empty = M.reference_null(pool, n=5, draws=5, seed=0, fit=_raising_stub())
    assert empty.n_draws_requested == empty.n_dropped_draws == 5
    assert empty.values.shape == (0, 2)
    # No null, no p-value: unavailable, not 1.0 and not 0.0.
    assert np.isnan(M.pvalue(empty.values[:, 0], 3.5))


def test_an_empty_null_is_unfit_and_never_reads_as_a_population_verdict(monkeypatch):
    """Through `test_population`: every joint-null draw raises, control and target fit.  Four
    NaN p-values compare False against 0.05 at every step of the verdict chain, and until
    2026-09-21 that fell through to `one population` -- a missing measurement read as a result.
    The verdict must be an UNFIT, spelled the way the module spells its other UNFIT."""
    rng = np.random.default_rng(0)
    space = ["CD99", "nucleus_elongation"]
    ref = pd.DataFrame(rng.normal(size=(600, 2)), columns=space)
    und = pd.DataFrame(rng.normal(size=(300, 2)), columns=space)

    draws = 4
    calls = {"n": 0}

    def fit(x: np.ndarray, seed: int) -> tuple[float, float]:
        calls["n"] += 1
        if calls["n"] <= draws:            # the joint null's draws, every one of them
            raise ValueError("the mixture fit did not converge")
        return 100.0 + calls["n"], 100.0 + calls["n"]
    monkeypatch.setattr(M, "fit_statistics", fit)

    row, feats, _, g2 = M.test_population(ref, und, space, "115", "large",
                                          draws=draws, draws_1d=4, seed=0)

    assert row["n_draws_requested"] == row["n_dropped_draws"] == draws
    for key in ("p_lr_control", "p_sep_control", "p_lr_target", "p_sep_target",
                "null_lr_median", "null_sep_median"):
        assert np.isnan(row[key]), key
    assert row["verdict"] == "UNFIT: null empty"
    assert row["verdict"].startswith("UNFIT: "), "the module's other UNFIT is `UNFIT: control rejects`"
    assert row["verdict"] not in ("one population", "displaced subpopulation", "two populations")
    assert g2 is None, "no posterior is fitted on an UNFIT stratum"
    # The control and target statistics themselves were computed and are reported as they are.
    assert row["lr_control"] == 105.0 and row["lr_target"] == 106.0
    # The per-feature nulls were intact; each ledger says so.
    assert len(feats) == len(space)
    for feat in feats:
        assert feat["n_draws_requested"] == 4 and feat["n_dropped_draws"] == 0


# ------------------------------------------- site 5: a value that cannot be a count, at the read


def _refused(out: pd.DataFrame, i: int, column: str, output: str, value: str) -> str:
    """The row is unavailable and its reason names the export column, the output it feeds, and
    the value as read -- the module's idiom for a blank count, with the value in place of `blank`."""
    assert out[G.COUNT_STATE].iloc[i] == G.UNAVAILABLE
    reason = out[G.COUNT_REASON].iloc[i]
    assert f"`{column}`" in reason and output in reason and value in reason, reason
    assert "blank" not in reason, "it was not blank; it was refused"
    return reason


def test_a_denominator_that_cannot_be_a_count_is_refused_at_the_read():
    """`Num Detections` is a count.  Until 2026-09-22 a negative, non-finite or fractional one
    passed straight into the proportions: measured on this build, -5 gave `measured` with every
    proportion NaN, inf gave `measured` with every proportion 0.0, and 2.5 gave `measured` with
    prop_Galectin_positive 2.0 -- a proportion of two.  Each is refused where the column is read:
    the denominator is null, the proportions have none, the state is `unavailable`, and the
    reason names the column, the output and the value.  The four counts on the row were read
    and are kept."""
    for bad, shown in ((-5, "-5"), (np.inf, "inf"), (2.5, "2.5")):
        out = G.derive_counts(_annotation_rows(**{"Num Detections": [10, bad]}))

        # The row with a denominator is what it always was.
        assert out[G.COUNT_STATE].iloc[0] == G.MEASURED, bad
        assert pd.isna(out[G.COUNT_REASON].iloc[0])
        assert out["prop_Galectin_positive"].iloc[0] == 5 / 10

        assert pd.isna(out["n_cells_qupath"].iloc[1]), bad
        assert [int(out[c].iloc[1]) for c in COUNT_NAMES] == [5, 2, 3, 3]
        for c in COUNT_NAMES:
            assert pd.isna(out[c.replace("count_", "prop_")].iloc[1]), (bad, c)
        reason = _refused(out, 1, "Num Detections", "n_cells_qupath", shown)
        assert "no denominator" in reason
        assert str(out["n_cells_qupath"].dtype) == "Int64", "only admitted values reach the cast"


def test_a_count_that_cannot_be_a_count_is_refused_at_the_read_and_not_derived_from():
    """`Num M2Galectin3` of -1: until 2026-09-22 it was cast to the integer -1, carried a
    proportion of -0.1, and the row's reason was the DISAGREEMENT one -- `count_Metabolically_
    active` reported the two passes at odds, when the input was never a count.  `inf` and 2.5
    did not get that far: the `Int64` cast raised (OverflowError, TypeError) on the whole
    frame.  All three are refused at the read: the count is null, so is everything derived from
    it, and the reason names the input, not the derivation."""
    for bad, shown in ((-1.0, "-1"), (np.inf, "inf"), (2.5, "2.5")):
        out = G.derive_counts(_annotation_rows(**{"Num M2Galectin3": [5.0, bad]}))

        assert out[G.COUNT_STATE].iloc[0] == G.MEASURED, bad
        assert int(out["count_Galectin_positive"].iloc[0]) == 5

        assert pd.isna(out["count_Galectin_positive"].iloc[1]), bad
        assert pd.isna(out["prop_Galectin_positive"].iloc[1]), bad
        assert pd.isna(out["count_Metabolically_active"].iloc[1]), "derived from the refused count"
        assert pd.isna(out["prop_Metabolically_active"].iloc[1])
        # The two passes that WERE counts on that row are still counts.
        assert int(out["count_Galectin_macrophage"].iloc[1]) == 2
        assert int(out["count_Macrophage"].iloc[1]) == 3
        reason = _refused(out, 1, "Num M2Galectin3", "count_Galectin_positive", shown)
        assert "exceeds" not in reason, "a refused input is not a classifier disagreement"
        assert str(out["count_Galectin_positive"].dtype) == "Int64"


def test_an_integral_float_and_a_numeric_string_are_the_integer_they_spell():
    """A CSV column with one blank reads as float64 and one with one stray character reads as
    object, so 10.0 and "10" are how a count of 10 arrives.  They are admitted as 10; the stray
    character on the same column is refused by name; a boolean column is not a count at all."""
    out = G.derive_counts(_annotation_rows(**{"Num Detections": [10.0, "10"]}))
    assert out["n_cells_qupath"].tolist() == [10, 10]
    assert (out[G.COUNT_STATE] == G.MEASURED).all()

    out = G.derive_counts(_annotation_rows(**{"Num Iba1_CD68": ["3", "three"]}))
    assert int(out["count_Macrophage"].iloc[0]) == 3
    assert pd.isna(out["count_Macrophage"].iloc[1])
    _refused(out, 1, "Num Iba1_CD68", "count_Macrophage", "'three'")

    out = G.derive_counts(_annotation_rows(**{"Num Iba1_CD68": [True, False]}))
    assert out["count_Macrophage"].isna().all(), "True is not a count of 1"
    _refused(out, 0, "Num Iba1_CD68", "count_Macrophage", "True")
    _refused(out, 1, "Num Iba1_CD68", "count_Macrophage", "False")


# ------------------------------------------------- site 6: a statistic that is not finite


def test_a_control_or_target_whose_statistic_is_not_finite_is_unfit(monkeypatch):
    """Through `test_population`: the joint null is intact (four draws, none dropped), but the
    CONTROL's fit (call 5) or the TARGET's (call 6) returns a non-finite statistic -- a zero
    pooled SD gives a NaN separation.  `pvalue` returns NaN for it, NaN compares False at every
    step of the verdict chain, and until 2026-09-22 the chain fell through to `one population`
    on a stratum whose control or target was never scored.  The verdict must be an UNFIT,
    spelled the way the module spells its other UNFITs, so `main()` counts it with them and not
    among the testable donors.  Nothing about the statistic or p-value changes: they are NaN
    and reported as NaN."""
    rng = np.random.default_rng(0)
    space = ["CD99", "nucleus_elongation"]
    ref = pd.DataFrame(rng.normal(size=(600, 2)), columns=space)
    und = pd.DataFrame(rng.normal(size=(300, 2)), columns=space)
    draws = 4

    for nan_call, side, other in ((5, "control", "target"), (6, "target", "control")):
        calls = {"n": 0}

        def fit(x: np.ndarray, seed: int) -> tuple[float, float]:
            calls["n"] += 1
            if calls["n"] == nan_call:
                return np.nan, np.nan
            return 100.0 + calls["n"], 100.0 + calls["n"]
        monkeypatch.setattr(M, "fit_statistics", fit)

        row, feats, _, g2 = M.test_population(ref, und, space, "115", "large",
                                              draws=draws, draws_1d=4, seed=0)

        assert row["n_draws_requested"] == draws and row["n_dropped_draws"] == 0, side
        assert np.isnan(row[f"lr_{side}"]) and np.isnan(row[f"sep_{side}"]), side
        assert np.isnan(row[f"p_lr_{side}"]) and np.isnan(row[f"p_sep_{side}"]), side
        # The other side was scored, above every null draw: (0 + 1) / (4 + 1).
        assert row[f"p_lr_{other}"] == row[f"p_sep_{other}"] == 1 / 5, other
        assert row["verdict"] == "UNFIT: statistic not finite", side
        assert row["verdict"].startswith("UNFIT: "), "spelled like `UNFIT: control rejects`"
        assert row["verdict"] not in ("one population", "displaced subpopulation",
                                      "two populations"), side
        assert g2 is None, "no posterior is fitted on an UNFIT stratum"
        assert len(feats) == len(space)
        for feat in feats:
            assert feat["n_draws_requested"] == 4 and feat["n_dropped_draws"] == 0
