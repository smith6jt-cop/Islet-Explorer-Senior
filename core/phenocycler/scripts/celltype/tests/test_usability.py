"""The statistics must behave on data whose answer is known, and fail closed on data that
has no answer.  These run without touching the cohort."""
from __future__ import annotations

import pathlib
import sys

import json

import numpy as np
import pandas as pd
import pytest

HERE = pathlib.Path(__file__).resolve().parents[1]
sys.path.insert(0, str(HERE))
import panel                                   # noqa: E402
import usability as U                          # noqa: E402


@pytest.fixture(scope="module")
def graph():
    rng = np.random.default_rng(0)
    xy = rng.random((20_000, 2)) * 1000.0
    return xy, U.neighbour_graph(xy, 12)


# ------------------------------------------------------------------ spatial coherence
def test_morans_i_is_zero_on_noise_and_high_on_a_blob(graph):
    xy, nbr = graph
    rng = np.random.default_rng(1)
    assert abs(U.morans_i(rng.lognormal(size=len(xy)), nbr)) < 0.01
    blob = np.exp(-((xy[:, 0] - 500) ** 2 + (xy[:, 1] - 500) ** 2) / 2e4)
    assert U.morans_i(blob, nbr) > 0.9


def test_permutation_null_brackets_the_noise_case(graph):
    xy, nbr = graph
    rng = np.random.default_rng(2)
    v = rng.lognormal(size=len(xy))
    null = U.morans_null([v], nbr, draws=10)
    assert 0 < null < 0.02, "the null must be small but not zero"
    assert abs(U.morans_i(v, nbr)) <= null * 3, "unpermuted noise sits inside its own null"


# ---------------------------------------------------------------------- tail overlap
def test_tail_overlap_is_one_for_independent_and_maximal_for_identical():
    rng = np.random.default_rng(3)
    a = rng.random(100_000)
    assert U.tail_overlap(a, a, 0.001) == pytest.approx(1 / 0.001, rel=0.01)
    lo = U.tail_overlap(a, rng.random(100_000), 0.001)
    assert lo < 5, "independent markers must not look enriched"


def test_tail_overlap_is_invariant_to_any_monotone_rescaling():
    """The whole point: it presupposes no scale, so a per-donor shift cannot change it."""
    rng = np.random.default_rng(4)
    a, b = rng.lognormal(size=50_000), rng.lognormal(size=50_000)
    base = U.tail_overlap(a, b, 0.001)
    assert U.tail_overlap(a * 137.0, b, 0.001) == base
    assert U.tail_overlap(np.log1p(a), b * 0.003, 0.001) == base


def test_tail_overlap_fails_closed_when_underpowered():
    rng = np.random.default_rng(5)
    n = U.MIN_CELLS - 1
    assert np.isnan(U.tail_overlap(rng.random(n), rng.random(n), 0.001))


# ------------------------------------------------------- acquisition-group structure
def test_acquisition_group_ratio_separates_group_steps_from_donor_spread():
    donors = [f"d{i}" for i in range(8)]
    groups = pd.Series(["A"] * 4 + ["B"] * 4, index=donors)
    stepped = pd.Series([1.0, 1.1, 0.9, 1.05, 500.0, 520.0, 480.0, 510.0], index=donors)
    scattered = pd.Series([1.0, 9.0, 0.5, 4.0, 2.0, 7.0, 0.8, 5.0], index=donors)
    assert U.acquisition_group_ratio(stepped, groups)["acq_ratio"] > U.BARS["acq_ratio_flag"]
    assert U.acquisition_group_ratio(scattered, groups)["acq_ratio"] < U.BARS["acq_ratio_flag"]


# ---------------------------------------------------------------------- the verdict
def test_partners_label_the_same_cell_and_controls_do_not():
    """The definition, asserted: hormones mark different endocrine cells, so they are
    controls; keratins mark the same epithelium, so they are partners."""
    partners = {(a, b) for a, b, _ in panel.LINEAGE_PARTNERS}
    partners |= {(b, a) for a, b, _ in panel.LINEAGE_PARTNERS}
    controls = {(a, b) for a, b in panel.DECLARED_CONTROLS}
    controls |= {(b, a) for a, b in panel.DECLARED_CONTROLS}
    for pair in [("INS", "GCG"), ("GCG", "SST"), ("INS", "SST")]:
        assert pair not in partners and pair in controls, pair
    for pair in [("Ker8_18", "Pan_Cytokeratin"), ("CD3e", "Granzyme_B"), ("Ker8_18", "EpCAM")]:
        assert pair in partners, pair


def _row(**kw):
    base = {"n": U.MIN_CELLS * 10, "zero_frac": 0.05, "morans_null": 0.001,
            "morans_rank": 0.2, "morans_top1": 0.2, "partner_best": 100.0,
            "control_worst": 0.5, "acq_ratio": 1.1}
    return base | kw


def test_verdict_passes_only_on_evidence_it_has():
    assert U.judge_marker(_row()).passed


@pytest.mark.parametrize("bad", [
    {"n": 10},                                          # under-powered
    {"zero_frac": 0.9},                                 # mostly exact zeros
    {"morans_rank": np.nan, "morans_top1": np.nan},     # no spatial evidence at all
    {"morans_rank": 0.0005, "morans_top1": 0.0005},     # inside its own null
])
def test_verdict_fails_closed_on_absent_or_underpowered_evidence(bad):
    assert not U.judge_marker(_row(**bad)).passed


def test_verdict_judges_the_measurement_and_nothing_else():
    """A pair that does not resolve, or a marker that is not comparable across acquisition
    groups, is REPORTED -- neither makes the measurement invalid.  CD11c is the real case:
    spatially coherent, cleanly aligned, and its declared partner CD209 marks other cells."""
    row = _row(partner_best=1.2, control_worst=0.0, acq_ratio=35.0)
    assert U.judge_marker(row).passed
    f = U.flag_marker(row)
    assert f["flag_pair"] and f["flag_acquisition_group"]
    assert f["flag_pair_why"] and f["flag_acquisition_group_why"]


def test_a_marker_with_no_declared_partner_is_not_flagged_for_it():
    f = U.flag_marker(_row(partner_best=np.nan, control_worst=np.nan))
    assert not f["flag_pair"]


def test_a_marker_on_everything_is_not_condemned_by_its_scattered_bright_tail():
    """E-cadherin's case: low morans_top1, real morans_rank, strong partner tail."""
    v = U.judge_marker(_row(morans_top1=0.0009, morans_rank=0.15))
    assert v.passed, v.reason


def test_partner_contradicted_by_its_control_is_flagged_not_condemned():
    row = _row(partner_best=0.5, control_worst=30.0)
    assert U.judge_marker(row).passed
    assert "more than its partner" in U.flag_marker(row)["flag_pair_why"]


# ---------------------------------------------------------------------------- panel
def test_every_declared_marker_exists_in_the_registry():
    known = set(panel.markers(include_variant=True))
    named = {m for a, b, _ in panel.LINEAGE_PARTNERS for m in (a, b)}
    named |= {m for a, b in panel.DECLARED_CONTROLS for m in (a, b)}
    named |= set(panel.CONTROL_REPRESENTATIVES.values())
    assert not named - known, f"declared but not in the registry: {sorted(named - known)}"


def test_every_marker_with_a_partner_has_two_controls_and_they_are_not_partners():
    partners, controls = panel.partner_map(), panel.control_map()
    for m in partners:
        got = controls.get(m, [])
        assert len(got) >= 2, f"{m} has no negative control; an enrichment alone is unreadable"
        assert not set(got) & {p for p, _ in partners[m]}, f"{m}: a partner is used as a control"
        assert m not in got


def test_derived_controls_come_from_a_different_compartment():
    """Cells of one compartment sit together and spillover crosses touching boundaries, so a
    same-compartment marker is a poor negative control.  Measured: non-partner pairs where
    both are immune overlap at median 1.61x against 0.00x otherwise."""
    comp, controls = panel.compartments(), panel.control_map()
    declared = {(a, b) for a, b in panel.DECLARED_CONTROLS}
    declared |= {(b, a) for a, b in panel.DECLARED_CONTROLS}
    for m, got in controls.items():
        for c in got:
            if (m, c) in declared:
                continue          # explicitly declared pairs may share a compartment
            assert not (comp.get(m) and comp.get(m) == comp.get(c)), f"{m}/{c} share {comp[m]}"


def test_compartment_is_inherited_from_partners_when_the_registry_declines():
    """Granzyme_B is a state marker with no registry compartment, but it partners CD8 and
    CD3e, so it must not be handed CD3e as a negative control."""
    comp = panel.compartments()
    assert comp["Granzyme_B"] == "Immune"
    assert "CD3e" not in panel.control_map()["Granzyme_B"]


def test_cohort_is_the_twenty_eligible_donors_in_four_acquisition_groups():
    c = panel.cohort()
    assert len(c) == 20 and c["donor_id"].is_unique
    assert not {"6457", "6579"} & set(c["donor_id"]), "centrally excluded donors must not appear"
    assert set(c["acq_group"]) == {"2.3.1/CD38", "2.4.0/CD38",
                                   "2.4.1/b_Catenin1", "2.4.1/neither"}


def test_panel_variant_markers_are_out_by_default():
    assert not set(panel.PANEL_VARIANT_MARKERS) & set(panel.markers())
    assert set(panel.PANEL_VARIANT_MARKERS) <= set(panel.markers(include_variant=True))
    assert not set(panel.REGISTRY_EXCLUDED) & set(panel.markers(include_variant=True))


# ------------------------------------------------------------------------ Step 2 / UniFORM
import uniform_adapter as UA                        # noqa: E402


def test_uniform_checkout_is_the_pinned_commit():
    """A drifted checkout can change the shift with nothing in the output able to show it."""
    reg, norm, pre = UA._import_uniform()
    assert all(hasattr(reg, f) for f in ("automatic_registration", "compute_correlation_shifts"))
    assert hasattr(norm, "calculate_shift_in_log_pixels")
    assert hasattr(pre, "fit_two_component_gmm")


def test_registration_recovers_a_known_scale_factor():
    """Three donors of the same distribution, scaled by 1x / 4x / 0.25x.

    This is the whole method on data whose answer is known: if it cannot recover a scale it
    was handed, no number it produces on the cohort is worth reading.
    """
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    reg, norm, _ = UA._import_uniform()
    rng = np.random.default_rng(0)
    base = rng.lognormal(3, 1, 50_000)
    truth = {"A": 1.0, "B": 4.0, "C": 0.25}
    vals = {d: base * f for d, f in truth.items()}
    donors = list(vals)

    counts, ranges = UA.marker_histograms(vals)
    hist = {"M": {"counts": [counts[d] for d in donors], "bin_edges": [None] * len(donors)}}
    shift_map, refs, _ = reg.automatic_registration(hist, ["M"], ["M"], donors,
                                                    None, None, UA.NUM_BINS)
    plt.close("all")
    log_shift = norm.calculate_shift_in_log_pixels({"M": ranges}, ["M"], UA.NUM_BINS, shift_map)
    scale = dict(zip(donors, (float(np.exp(-x)) for x in log_shift["M"])))

    ref_f = truth[refs[0]]
    for d in donors:
        assert scale[d] == pytest.approx(ref_f / truth[d], rel=0.10), (d, scale[d])

    before = UA.background_mode_spread(counts, ranges)
    after = UA.background_mode_spread(counts, ranges, dict(zip(donors, shift_map["M"])))
    assert after["sd_log"] < before["sd_log"] / 5


def test_the_shift_is_a_pure_multiplicative_rescaling():
    """Monotone, so every within-donor rank statistic is invariant under it.

    That is the property the Step 1 tail statistics rely on to check Step 2 without
    circularity, so it is asserted rather than assumed.
    """
    rng = np.random.default_rng(1)
    raw = rng.lognormal(2, 1.5, 20_000)
    partner = rng.lognormal(2, 1.5, 20_000)
    for s in (0.013, 1.0, 87.0):
        assert np.array_equal(np.argsort(raw), np.argsort(raw * s))
        assert U.tail_overlap(raw * s, partner, 0.001) == U.tail_overlap(raw, partner, 0.001)


def test_background_mode_spread_is_zero_for_identical_donors():
    rng = np.random.default_rng(2)
    v = rng.lognormal(3, 1, 30_000)
    counts, ranges = UA.marker_histograms({"A": v.copy(), "B": v.copy(), "C": v.copy()})
    s = UA.background_mode_spread(counts, ranges)
    assert s["sd_log"] == 0.0 and s["fold"] == pytest.approx(1.0)


# ------------------------------------------------------------------- the run this reads
def test_the_expression_stage_is_the_one_this_lane_was_built_against():
    """Pin the data, not just the run id.

    The run choice is immaterial -- both runs on disk hardlink the same expression files --
    but that is only true while the content hash holds.  Pinning it means a swapped run or a
    mutated stage fails here instead of silently changing every number downstream.
    """
    m = json.loads((panel.PKG / "data" / "runs" / panel.RUN_ID /
                    "manifests" / "expression.json").read_text())
    assert m["output"]["content_sha256"] == panel.EXPRESSION_CONTENT_SHA256
    assert m["output"]["total_rows"] == panel.EXPRESSION_TOTAL_ROWS
    assert m["method_version"] == "selected-expression-v1"


def test_the_other_run_carries_the_same_expression_bytes():
    """The stated justification, asserted rather than remembered."""
    other = "c8412757d12de36ca056"
    root = panel.PKG / "data" / "runs"
    if not (root / other).exists():
        pytest.skip(f"{other} not on this machine")
    a = json.loads((root / panel.RUN_ID / "manifests" / "expression.json").read_text())
    b = json.loads((root / other / "manifests" / "expression.json").read_text())
    assert a["output"]["content_sha256"] == b["output"]["content_sha256"]
    assert a["output"]["schema_sha256"] == b["output"]["schema_sha256"]
