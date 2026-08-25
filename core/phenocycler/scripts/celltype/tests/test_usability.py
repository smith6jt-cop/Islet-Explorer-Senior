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


# ------------------------------------------------------------ the rule table (Step 3)
import rules as R                                    # noqa: E402

#: The rule table lives outside version control and is revised by the maintainer, so its
#: absence is a normal state rather than a failure.  These tests check the machinery against
#: whichever table is present; with none, there is nothing to check.
needs_rules = pytest.mark.skipif(not R.RULES_CSV.exists(),
                                 reason=f"no rule table at {R.RULES_CSV}")


@needs_rules
def test_crosswalk_is_frozen_complete_injective_and_agrees_with_the_regex():
    """A regex that silently misses one name turns that rule's constraint into no
    constraint, which fails open.  The frozen map is the guard; this asserts it is still
    right, so a v5 rename fails here instead of dropping evidence."""
    cw = R.load_crosswalk()
    cols = R._marker_columns()
    assert set(cols) <= set(cw), f"unmapped marker columns: {set(cols) - set(cw)}"
    assert len(set(cw.values())) == len(cw), "crosswalk is not injective"
    known = set(panel.markers(include_variant=True)) | set(panel.REGISTRY_EXCLUDED)
    assert not set(cw.values()) - known
    assert all(cw[c] == t for c, t in R.build_crosswalk().items())


@needs_rules
def test_flattening_composes_a_phenotype_with_its_ancestors():
    """A row is a partial statement, whichever layout the table uses.

    In the 59-plex table a child drops between 2 and 8 of its parent's constraints and
    overrides up to 4, so reading a row alone loses real evidence.
    """
    rules = R.load_rules()
    flat = R.flatten(rules)
    row = {r["phenotype"]: {c: r[c] for c in rules.columns
                            if c not in ("parent", "phenotype") and pd.notna(r[c])}
           for _, r in rules.iterrows()}
    par = dict(zip(rules["phenotype"], rules["parent"]))
    gained = [len(flat[ph]) - len(row[ph]) for ph in row if par.get(ph) not in (None, "all")]
    assert gained and max(gained) > 0, "no phenotype gained a constraint from its ancestors"
    for ph, p in par.items():
        if p in flat:
            assert set(flat[p]) <= set(flat[ph]), f"{ph} lost a constraint of {p}"


@needs_rules
def test_the_hierarchy_is_read_from_whichever_layout_the_table_uses():
    rules = R.load_rules()
    assert list(rules.columns[:2]) == ["parent", "phenotype"]
    assert rules["phenotype"].is_unique
    known = set(rules["phenotype"]) | {"all"}
    assert set(rules["parent"]) <= known, set(rules["parent"]) - known


def test_flattening_refuses_a_cycle_rather_than_looping():
    bad = pd.DataFrame({"parent": ["B", "A"], "phenotype": ["A", "B"], "CD3e": ["pos", "neg"]})
    with pytest.raises(ValueError, match="cycle"):
        R.flatten(bad)


def test_a_variant_marker_is_encoded_and_costs_donors_not_the_rule():
    """The cohort is three panel variants, so availability is a per-donor property.

    A rule naming CD38 is evaluable on the 13 donors that carry it, not unencodable.  Built
    over the cohort intersection instead, its constraint vanishes and B_cell becomes
    identical to CD38_high_B with nothing to show it -- which is the fail-open this guards.
    """
    t = pd.DataFrame({
        "parent": ["all", "B_cell"], "phenotype": ["B_cell", "CD38_high_B"],
        "CD20": ["pos", None], "CD38": [None, "pos"], "CD3e": ["neg", None]})
    flat = R.flatten(t)
    M, dropped = R.reference_matrix(flat, panel.markers(include_variant=True), "zero")
    assert not dropped
    assert not (M["B_cell"] == M["CD38_high_B"]).all()

    ev = R.evaluable_donors(flat).set_index("phenotype")
    assert ev.loc["B_cell", "n_donors"] == 20
    assert ev.loc["CD38_high_B", "n_donors"] == 13
    assert ev.loc["CD38_high_B", "variant_markers"] == "CD38"

    thin, dropped = R.reference_matrix(flat, panel.markers(), "zero", strict=False)
    assert dropped["CD38"], "the intersection must report what it cannot represent"
    assert (thin["B_cell"] == thin["CD38_high_B"]).all()


def test_marker_availability_is_a_property_of_the_donor():
    dm = panel.donor_markers()
    assert len(dm) == 20
    assert len(panel.donors_with(["CD38"])) == 13
    assert len(panel.donors_with(["b_Catenin1"])) == 4
    assert len(panel.donors_with(["CD38", "b_Catenin1"])) == 0, "no donor carries both"
    assert len(panel.donors_with(["CD3e"])) == 20
    for d in ("6523", "6566", "6591"):
        got = set(panel.donor_markers(d))
        assert not got & set(panel.PANEL_VARIANT_MARKERS), f"{d} carries neither"


def test_cd11b_is_in_the_panel_and_reaches_the_reference_matrix():
    """Excluded only for the retired RESTORE lane.  Every table so far constrains CD11b --
    the v4 Myeloid gate is CD68|CD11c|CD11b, the 59-plex table uses it on 20 of 36
    phenotypes -- so dropping it would silently remove that evidence."""
    assert "CD11b" in panel.markers()
    assert "IAPP" not in panel.markers(include_variant=True)
    if not R.RULES_CSV.exists():
        pytest.skip(f"no rule table at {R.RULES_CSV}")
    flat = R.flatten(R.load_rules())
    assert any("CD11b" in ref for ref in flat.values()), "no rule constrains CD11b"
    M, dropped = R.reference_matrix(flat, panel.markers(include_variant=True), "zero")
    assert "CD11b" not in dropped and "CD11b" in M.index


@needs_rules
def test_signed_encoding_resolves_phenotypes_that_zero_cannot():
    flat = R.flatten(R.load_rules())
    mk = panel.markers(include_variant=True)
    z, _ = R.reference_matrix(flat, mk, "zero", strict=False)
    s, _ = R.reference_matrix(flat, mk, "signed", strict=False)
    nz = len(z.columns) - sum(len(v) for v in R.distinct_references(z).values())
    ns = len(s.columns) - sum(len(v) for v in R.distinct_references(s).values())
    assert ns > nz, "signed must resolve at least the negation-defined phenotypes"
    assert (s.to_numpy() < 0).any() and not (z.to_numpy() < 0).any()


# ------------------------------------------------ the plan's verification items, asserted
def test_the_acquisition_group_statistic_reproduces_from_the_cohort_table_alone():
    """Plan verification item 4.

    CD57 must come back at ~35 and every other marker below the flag, derived from the donor
    IDs and the group table -- not read from a stored result.  This is the one measurement
    Step 2 acts on, so it is pinned rather than remembered.
    """
    out = panel.OUT / "marker_donor_usability.parquet"
    if not out.exists():
        pytest.skip("Step 1 has not been run on this machine")
    d = pd.read_parquet(out)
    grp = panel.cohort().set_index("donor_id")["acq_group"]
    r = {m: U.acquisition_group_ratio(g.set_index("donor_id")["p50"], grp)["acq_ratio"]
         for m, g in d.groupby("marker")}
    assert r["CD57"] == pytest.approx(35.1, abs=0.5)
    others = {m: v for m, v in r.items() if m != "CD57" and np.isfinite(v)}
    assert max(others.values()) < U.BARS["acq_ratio_flag"], \
        f"a second marker steps with the group: {max(others, key=others.get)}"


def test_step2_records_a_reference_donor_shift_and_provenance_for_every_marker():
    """Plan verification item 5, in the part that survived.

    The landmark half is void: refinement was measured and made every marker worse, so no
    landmark_annotations.csv exists and nothing is reproduced from one.  What must hold is
    that every shift is attributable.
    """
    f = panel.OUT / "shifts.json"
    if not f.exists():
        pytest.skip("Step 2 has not been run on this machine")
    meta = json.loads(f.read_text())
    assert meta["uniform_commit"]
    for m, v in meta["shifts"].items():
        assert v["reference_donor"] in meta["donors"], m
        assert set(v["scale"]) == set(meta["donors"]), m
        assert v["provenance"] in ("automatic_registration",
                                   "automatic_registration_over_bar",
                                   "refused_reagent_change"), (m, v["provenance"])
    refused = [m for m, v in meta["shifts"].items()
               if v["provenance"] == "refused_reagent_change"]
    assert refused == ["CD57"]
    assert all(s == 1.0 for s in meta["shifts"]["CD57"]["scale"].values()), \
        "a refused marker must be left unshifted, not shifted and labelled"


# ------------------------------------------------------------------- scoring (Step 3)
import score as SC                                    # noqa: E402


def test_cosine_is_invariant_to_scaling_the_whole_cell_vector():
    """The property this lane adopted cosine for: the dominant axis between cells is overall
    brightness, and rescaling a cell must not change which reference it matches."""
    rng = np.random.default_rng(0)
    X = rng.random((500, 8)).astype(np.float32)
    Rm = rng.random((8, 4)).astype(np.float32)
    a = SC.cosine_scores(X, Rm, normalise_ceiling=False)
    b = SC.cosine_scores(X * 37.0, Rm, normalise_ceiling=False)
    assert np.allclose(a, b, atol=1e-5)


def test_cosine_cannot_manufacture_a_pairwise_association():
    """Unlike dividing by a chosen marker, which turned a negative control from 0.1x to
    14.8x, rescaling the whole vector leaves every pairwise relationship intact."""
    rng = np.random.default_rng(1)
    X = rng.random((2000, 6)).astype(np.float32)
    order_before = np.argsort(X[:, 0])
    X2 = X * rng.random((2000, 1)).astype(np.float32)     # per-cell magnitude, not per-marker
    Rm = np.eye(6, dtype=np.float32)
    S1 = SC.cosine_scores(X, Rm, normalise_ceiling=False)
    S2 = SC.cosine_scores(X2, Rm, normalise_ceiling=False)
    assert np.allclose(S1, S2, atol=1e-5), "a per-cell rescaling changed the angles"
    assert order_before is not None


def test_a_reference_with_negatives_cannot_reach_one_and_the_ceiling_says_so():
    """Cell values are non-negative, so the best any cell can do against r is align with its
    positive part: max cos = ||r+|| / ||r||.  Measured on donor 6450 that ran 1.000 down to
    0.147 across the 37 references and predicted the call share at Spearman +0.48, so the
    argmax was partly choosing on reference geometry."""
    R_pos = np.array([[1.0], [1.0], [0.0]], dtype=np.float32)
    R_mix = np.array([[1.0], [-1.0], [-1.0]], dtype=np.float32)
    assert SC.reference_ceiling(R_pos)[0] == pytest.approx(1.0)
    assert SC.reference_ceiling(R_mix)[0] == pytest.approx(1 / np.sqrt(3), rel=1e-5)

    rng = np.random.default_rng(2)
    X = rng.random((3000, 3)).astype(np.float32)
    raw = SC.cosine_scores(X, np.hstack([R_pos, R_mix]), normalise_ceiling=False)
    assert np.nanmax(raw[:, 1]) < np.nanmax(raw[:, 0]), "the mixed reference must be capped lower"
    norm = SC.cosine_scores(X, np.hstack([R_pos, R_mix]), normalise_ceiling=True)
    assert np.nanmax(norm[:, 1]) == pytest.approx(1.0, abs=0.05), "ceiling did not equalise them"


def test_a_cell_with_no_signal_is_unscored_rather_than_called():
    """Zero norm means there is no angle.  Calling such a cell would be inventing evidence."""
    X = np.array([[0.0, 0.0], [1.0, 0.0]], dtype=np.float32)
    Rm = np.array([[1.0, 0.0], [0.0, 1.0]], dtype=np.float32)
    out = SC.call_from_scores(SC.cosine_scores(X, Rm), ["A", "B"])
    assert out.loc[0, "call"] == "Unscored" and np.isnan(out.loc[0, "score"])
    assert out.loc[1, "call"] == "A"


def test_the_margin_records_how_close_the_decision_was():
    X = np.array([[1.0, 1.0]], dtype=np.float32)          # exactly between two references
    Rm = np.eye(2, dtype=np.float32)
    out = SC.call_from_scores(SC.cosine_scores(X, Rm), ["A", "B"])
    assert out.loc[0, "margin"] == pytest.approx(0.0, abs=1e-6)


def test_scaling_clips_at_the_measured_floor_and_a_robust_top():
    scale = pd.DataFrame({"floor": [10.0], "top": [110.0], "span": [100.0]}, index=["M"])
    X = np.array([[0.0], [10.0], [60.0], [110.0], [1e6]], dtype=np.float32)
    out = SC.apply_scale(X, ["M"], scale)
    assert out.min() == 0.0 and out.max() == 1.0
    assert out[2, 0] == pytest.approx(0.5)


def test_the_dapi_ratio_is_applied_only_to_the_raw_nucleus_group():
    """DAPI is the denominator, never a numerator, and no REDSEA whole-cell marker is
    touched -- dividing one of those by DAPI compares two different measurements."""
    df = pd.DataFrame({"DAPI": [2.0, 4.0], "TOX": [1.0, 8.0], "CD3e": [3.0, 6.0]})
    out, changed = SC.nuclear_ratio(df.copy(), ["TOX", "CD3e"], "DAPI",
                                    panel.NUCLEUS_MEASURED)
    assert changed == ["TOX"]
    assert list(out["TOX"]) == [0.5, 2.0]
    assert list(out["CD3e"]) == [3.0, 6.0], "a whole-cell marker must not be divided"
    assert list(out["DAPI"]) == [2.0, 4.0], "the denominator must not divide itself"


def test_the_dapi_ratio_survives_a_zero_or_missing_denominator():
    df = pd.DataFrame({"DAPI": [0.0, np.nan, 4.0], "TOX": [1.0, 1.0, 8.0]})
    out, changed = SC.nuclear_ratio(df.copy(), ["TOX"], "DAPI", panel.NUCLEUS_MEASURED)
    assert changed == ["TOX"]
    assert np.isnan(out["TOX"].iloc[0]) and np.isnan(out["TOX"].iloc[1])
    assert out["TOX"].iloc[2] == 2.0
    unchanged, changed = SC.nuclear_ratio(df.drop(columns=["DAPI"]).copy(), ["TOX"], "DAPI",
                                          panel.NUCLEUS_MEASURED)
    assert changed == [] and list(unchanged["TOX"]) == [1.0, 1.0, 8.0]


# ---------------------------------------------------- the requirement gate (Step 3)
def test_an_or_group_is_satisfied_by_a_child_committing_to_one_member():
    """Epi-endocrine declares "any of B3TUBB|GCG|INS|SST"; Endo-alpha commits GCG to allpos
    and INS/SST to allneg.  Read off the flattened tokens only B3TUBB still carries anypos,
    so the group collapses to a hard requirement for it -- which took alpha cells from
    thousands to 79.  Read off the declaring row, GCG satisfies it and the group is dropped."""
    t = pd.DataFrame({
        "parent": ["all", "Endocrine"], "phenotype": ["Endocrine", "Alpha"],
        "CD56": ["allpos", None], "B3TUBB": ["anypos", None],
        "GCG": ["anypos", "allpos"], "INS": ["anypos", "allneg"], "SST": ["anypos", "allneg"]})
    rq = R.requirements(t)
    assert rq["Endocrine"] == (["CD56"], [["B3TUBB", "GCG", "INS", "SST"]])
    assert rq["Alpha"] == (["CD56", "GCG"], []), rq["Alpha"]


def test_two_or_groups_on_one_row_stay_separate():
    """T-exhausted is "(CD4 or CD8) AND (VISTA or LAG3 or TOX)".  Merged into one OR it would
    be satisfied by TOX alone, which is the opposite of requiring a lineage co-receptor."""
    t = pd.DataFrame({
        "parent": ["all"], "phenotype": ["Tex"],
        "CD3e": ["allpos"], "CD4": ["anypos2"], "CD8": ["anypos2"],
        "TOX": ["anypos"], "LAG3": ["anypos"]})
    all_of, groups = R.requirements(t)["Tex"]
    assert all_of == ["CD3e"]
    assert sorted(groups) == [["CD4", "CD8"], ["LAG3", "TOX"]], groups


def test_an_inherited_group_is_dropped_once_the_child_meets_it():
    """Imm-T-cell commits CD3e, which already satisfies Immune's "is an immune cell" group,
    so a T subset must not also be asked for CD68 or CD20."""
    t = pd.DataFrame({
        "parent": ["all", "Immune"], "phenotype": ["Immune", "Tcell"],
        "CD3e": ["anypos", "allpos"], "CD68": ["anypos", None], "CD20": ["anypos", None]})
    all_of, groups = R.requirements(t)["Tcell"]
    assert all_of == ["CD3e"] and groups == []


def test_the_gate_sets_the_candidates_and_the_cosine_only_ranks_within_them():
    X = np.array([[1.0, 0.0], [0.0, 1.0]], dtype=np.float32)
    Rm = np.array([[1.0, 0.9], [0.0, 0.1]], dtype=np.float32)   # A would win both on angle
    reqs = {"A": (["m0"], []), "B": (["m1"], [])}
    allowed = SC.meets_requirements(X, ["m0", "m1"], reqs, ["A", "B"], 0.2)
    out = SC.call_from_scores(SC.cosine_scores(X, Rm), ["A", "B"], allowed=allowed)
    assert list(out["call"]) == ["A", "B"], "the gate did not override the angle"


def test_a_cell_meeting_no_rule_is_ambiguous_with_a_reason():
    X = np.array([[0.01, 0.01]], dtype=np.float32)
    reqs = {"A": (["m0"], []), "B": (["m1"], [])}
    allowed = SC.meets_requirements(X, ["m0", "m1"], reqs, ["A", "B"], 0.2)
    out = SC.call_from_scores(SC.cosine_scores(X, np.eye(2, dtype=np.float32)),
                              ["A", "B"], allowed=allowed)
    assert out.loc[0, "call"] == "Ambiguous"
    assert out.loc[0, "reason"] == "no_reference_requirements_met"
    assert out.loc[0, "n_candidates"] == 0


# ------------------------------------------------------------- tissue domains (Step 4)
import utag_adapter as UT                              # noqa: E402


def test_utag_checkout_is_the_pinned_commit():
    utag = UT._import_utag()
    import inspect
    params = inspect.signature(utag).parameters
    for k in ("adata", "slide_key", "max_dist", "clustering_method", "resolutions"):
        assert k in params, k


def test_the_domain_input_never_contains_the_geometry_it_is_checked_against():
    """`cell_region` is polygon containment and the domains are expression plus position.
    Feeding the first into the second would make the concordance check circular."""
    import panel as P
    a = UT.donor_adata(P.cohort()["donor_id"].iloc[0], ["CD3e", "CD20", "E_cadherin"],
                       max_cells=2000)
    assert "cell_region" not in list(a.var_names)
    assert "islet_num" not in list(a.var_names)
    assert {"cell_region", "islet_num"} <= set(a.obs.columns), "carried for the check, not as input"
    assert a.obsm["spatial"].shape == (a.n_obs, 2)


def test_domain_concordance_is_an_enrichment_over_the_donor_own_base_rate():
    obs = pd.DataFrame({
        "d": ["islet"] * 100 + ["rest"] * 900,
        "cell_region": ["core"] * 50 + ["tissue"] * 50 + ["tissue"] * 900,
    })
    t = UT.domain_islet_concordance(obs, "d").set_index("domain")
    assert t.loc["islet", "core_frac"] == pytest.approx(0.5)
    assert t.loc["islet", "base_rate"] == pytest.approx(0.05)
    assert t.loc["islet", "enrichment"] == pytest.approx(10.0)
    assert t.loc["rest", "enrichment"] == pytest.approx(0.0)
