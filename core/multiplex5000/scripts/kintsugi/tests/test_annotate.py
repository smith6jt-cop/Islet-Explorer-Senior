"""Annotation: the hierarchy resolves nested candidates, the gate rests on the background."""
import pathlib
import sys

import numpy as np
import pandas as pd
import pytest

HERE = pathlib.Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent.parent / "celltype"))
sys.path.insert(0, str(HERE.parent))

import annotate                                     # noqa: E402
import rules                                        # noqa: E402


@pytest.fixture(scope="module")
def anc():
    return annotate.ancestry(rules.load_rules())


def test_nested_candidates_resolve_to_the_most_specific(anc):
    """flatten() composes with ancestors, so a leaf always drags its parents in."""
    call, why = annotate.resolve(["Imm-macrophage", "Immune", "Mac-M2-like"], anc)
    assert call == "Mac-M2-like" and why == ""


def test_sibling_candidates_resolve_to_their_common_ancestor(anc):
    """Four T subsets the evidence cannot separate still support 'this is a T cell'."""
    call, why = annotate.resolve(
        ["Imm-T-cell", "Immune", "T-CD8", "T-CD8-cytotoxic", "T-exhausted",
         "T-progenitor-memory"], anc)
    assert call == "Imm-T-cell"
    assert why == "resolved_to_common_ancestor"


def test_candidates_in_different_branches_are_ambiguous(anc):
    """Endothelial and Immune share no ancestor below `all`; picking one would be a guess."""
    call, why = annotate.resolve(["Endo-vascular", "Endothelial", "Immune"], anc)
    assert call == ""
    assert why == "candidates_span_multiple_branches"


def test_no_candidates_is_ambiguous_with_a_reason(anc):
    call, why = annotate.resolve([], anc)
    assert call == "" and why == "no_reference_requirements_met"


def test_ancestry_refuses_a_cycle():
    table = pd.DataFrame({"parent": ["B", "A"], "phenotype": ["A", "B"], "CD3e": [None, None]})
    with pytest.raises(ValueError, match="cycle"):
        annotate.ancestry(table)


def test_background_mode_finds_the_mode_not_the_minimum():
    """The floor is the whole gate: min() is the REDSEA clamp, not the background."""
    rng = np.random.default_rng(0)
    background = rng.normal(100.0, 5.0, 200_000)
    positive = rng.normal(2000.0, 200.0, 8_000)
    values = np.concatenate([background, positive, np.zeros(5_000)])
    mode = annotate.background_mode(values)
    assert 80.0 < mode < 130.0, mode
    assert mode > values.min() + 50.0


def test_gate_is_measured_from_the_background_mode():
    """With the right floor, background cells are not 'elevated' and positives are."""
    rng = np.random.default_rng(1)
    markers = ["M0", "M1"]
    values = np.column_stack([
        np.concatenate([rng.normal(100, 5, 50_000), rng.normal(3000, 200, 2_000)]),
        np.concatenate([rng.normal(100, 5, 51_000), rng.normal(2500, 200, 1_000)]),
    ])
    scale = annotate.marker_scale(values, markers)
    scaled = (values - scale.loc[markers, "floor"].to_numpy()) / \
        scale.loc[markers, "span"].to_numpy()
    hot = scaled >= annotate.MIN_ELEVATED
    assert hot[:50_000, 0].mean() < 0.02          # background is not called positive
    assert hot[50_000:, 0].mean() > 0.90          # the positive population is
    assert hot[:51_000, 1].mean() < 0.02


def test_a_marker_with_no_positive_population_is_not_gateable():
    """A known and inherited property, recorded rather than papered over.

    The scale is ``[background mode, p99.9]``.  For a marker carrying no positive population
    the robust top IS the background's own upper tail, so the span collapses and a third of
    pure background clears ``MIN_ELEVATED``.  ``celltype/score.marker_scale`` has the same
    form, so this is inherited rather than introduced -- and it is precisely what Step 1's
    usability verdict exists to catch: a marker with no spatial structure above its own
    permutation null and no partner tail is ``unusable`` and never reaches this gate.

    Pinned so that if the scale is ever changed, the change is deliberate.
    """
    rng = np.random.default_rng(2)
    flat = rng.normal(100.0, 5.0, 52_000)[:, None]
    scale = annotate.marker_scale(flat, ["FLAT"])
    scaled = (flat - scale.loc["FLAT", "floor"]) / scale.loc["FLAT", "span"]
    assert (scaled >= annotate.MIN_ELEVATED).mean() > 0.20, (
        "if this no longer holds the scale changed; check it was on purpose")
