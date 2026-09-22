"""The refinement pass: the rule, the per-donor x per-stratum bound, and the single-pass guard.

Run from the repository root with .venv/bin/python -m pytest <this file>.

Fixtures are synthetic measurements.  Nothing here asserts a biological proportion, a demotion
rate or a cell count from the cohort: the pass is a rule about measurements, and a test that
reproduced one of its outputs would be measuring the expectation rather than the build.  The
two tests that touch the repository read the RULE TABLE and ``shifts.json`` -- inputs, not
results -- to check that the marker set is derived and not restated.
"""

from __future__ import annotations

import importlib.util
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pyarrow as pa
import pyarrow.parquet as pq
import pytest

REPO = Path(__file__).resolve().parents[4]
MODULE = REPO / "core/phenocycler/scripts/typing/refine.py"
sys.path.insert(0, str(MODULE.parent))
sys.path.insert(0, str(MODULE.parent.parent / "celltype"))
spec = importlib.util.spec_from_file_location("refine", MODULE)
refine = importlib.util.module_from_spec(spec)
spec.loader.exec_module(refine)


# --------------------------------------------------------------------- fixtures

MARKERS = ["INS", "GCG", "SST", "CD31", "Caveolin", "CD34", "CD3e", "CD68", "Iba1",
           "CD20", "CD79a", "MPO", "CD66", "CD11c", "CD209", "Ker8_18", "E_cadherin",
           "Pan_Cytokeratin", "B3TUBB", "CD56", "Podoplanin"]


def matrix(spec: dict, seed: int = 0):
    """A synthetic Z: per class, n cells whose own markers sit at a stated level in SD.

    ``spec`` is ``{class: [(n, own_level, stratum), ...]}`` or, for a group that must sit
    among cells unlike it, ``(n, own_level, stratum, "Exocrine")`` -- which also raises that
    class's own markers, so the cell lands in ITS cloud in marker space.  That is what makes a
    dim cell scattered, and a scattered cell is excluded from its own class's bound, exactly
    as it is on the cohort.  Every other marker is drawn at the population background, so the
    only thing separating two classes is whether each carries its own marker.
    """
    rng = np.random.default_rng(seed)
    rows, labels, strata = [], [], []
    for cls, groups in spec.items():
        own = [MARKERS.index(m) for m in refine.OWN_MARKERS[cls] if m in MARKERS]
        for group in groups:
            n, level, stratum = group[:3]
            like = group[3] if len(group) > 3 else None
            block = rng.normal(0.0, 0.05, size=(n, len(MARKERS)))
            block[:, own] = level + rng.normal(0.0, 0.02, size=(n, len(own)))
            if like is not None:
                cols = [MARKERS.index(m) for m in refine.OWN_MARKERS[like] if m in MARKERS]
                block[:, cols] = 2.0 + rng.normal(0.0, 0.02, size=(n, len(cols)))
            rows.append(block)
            labels += [cls] * n
            strata += [stratum] * n
    Z = np.vstack(rows).astype("float32")
    return Z, np.array(labels, dtype=object), np.array(strata, dtype=object)


def cells_frame(rows):
    return pd.DataFrame(rows, columns=["object_id", "cell_class", "cell_region_20um",
                                       "islet_num", "qc_analysis_eligible"])


# ----------------------------------------------------- the constants are pre-registered


def test_constants_are_the_recovered_ones():
    assert refine.LOW == 0.40
    assert refine.K == 30
    assert refine.BOUND_PCT == 5.0
    assert refine.MIN_CLEAR == 50
    assert refine.UNRESOLVED == "Unresolved (no lineage evidence)"


def test_require_min_is_imported_not_restated():
    """0.10 is pre-registered in the celltype lane; the pass must not carry its own copy."""
    import run_step3

    assert refine.immune_leaf.__defaults__[-1] is run_step3.REQUIRE_MIN
    assert run_step3.REQUIRE_MIN == 0.10


def test_own_markers_are_the_recovered_table():
    assert refine.OWN_MARKERS == {
        "Nerve": ("B3TUBB", "CD56"),
        "Lymphatic": ("Podoplanin",),
        "Blood vessel": ("CD31", "Caveolin", "CD34"),
        "Beta": ("INS",),
        "Alpha": ("GCG",),
        "Delta": ("SST",),
        "Immune": ("CD3e", "CD68", "Iba1", "CD20", "CD79a", "MPO", "CD66", "CD11c", "CD209"),
        "Exocrine": ("Ker8_18", "E_cadherin", "Pan_Cytokeratin"),
    }


# ------------------------------------------------------------------ the marker set


@pytest.mark.skipif(not (REPO / "data/development/celltype/shifts.json").exists(),
                    reason="shifts.json is not in this checkout")
@pytest.mark.data
def test_marker_set_is_derived_and_the_drop_list_with_it():
    markers = refine.derive_markers()
    assert len(markers) > 20
    for dropped in ("CD38", "b_Catenin1", "CD57"):
        assert dropped not in markers
    # CD57 must fall out through provenance, never by name: a future reagent change has to be
    # dropped by the same route without anyone editing this module.
    assert "CD57" not in refine.DROP_BY_NAME
    assert refine.DROP_PROVENANCE == "refused_reagent_change"
    for cls, own in refine.OWN_MARKERS.items():
        assert [m for m in own if m in markers], f"{cls} has no own marker in the derived set"


# ------------------------------------------------------------------------ the rule


def test_a_cell_without_its_own_marker_is_demoted():
    Z, lab, st = matrix({
        "Blood vessel": [(400, 3.0, ""), (16, -1.0, "", "Exocrine")],
        "Exocrine": [(600, 2.0, "")],
    })
    got = refine.apply_once(Z, MARKERS, lab, st, np.zeros(len(lab), bool))
    dim = np.zeros(len(lab), bool)
    dim[400:416] = True
    assert got["demoted"][dim].all()
    # Every other class loses only the bottom of its own percentile -- the floor the rule has
    # by construction, tested on its own below -- and nothing more.
    assert got["demoted"][lab == "Exocrine"].mean() < 0.07
    clear = (lab == "Blood vessel") & ~dim
    assert got["demoted"][clear].mean() < 0.07


def test_a_demoted_region_cell_falls_back_to_exocrine_without_immune_evidence():
    Z, lab, st = matrix({
        "Blood vessel": [(400, 3.0, ""), (16, -1.0, "", "Exocrine")],
        "Exocrine": [(600, 2.0, "")],
    })
    got = refine.apply_once(Z, MARKERS, lab, st, np.zeros(len(lab), bool))
    assert set(got["labels"][400:416]) == {"Exocrine"}


def test_the_same_demoted_cell_falls_back_to_immune_when_it_satisfies_a_leaf():
    """The recovered gap: without ``immune_leaf`` every demoted region cell went to Exocrine."""
    Z, lab, st = matrix({
        "Blood vessel": [(400, 3.0, ""), (16, -1.0, "", "Exocrine")],
        "Exocrine": [(600, 2.0, "")],
    })
    leaf = np.zeros(len(lab), bool)
    leaf[400:416] = True
    got = refine.apply_once(Z, MARKERS, lab, st, leaf)
    assert set(got["labels"][400:416]) == {"Immune"}


def test_a_demoted_exocrine_cell_has_no_complement_and_becomes_unresolved():
    Z, lab, st = matrix({
        "Exocrine": [(600, 2.0, ""), (24, -2.0, "", "Immune")],
        "Immune": [(200, 3.0, "")],
    })
    leaf = np.ones(len(lab), bool)          # even with immune evidence it is not made Immune
    got = refine.apply_once(Z, MARKERS, lab, st, leaf)
    assert set(got["labels"][600:624]) == {refine.UNRESOLVED}


def test_a_percentile_bound_cuts_about_five_percent_whether_or_not_anything_is_wrong():
    """A property of the rule, encoded so it is never read as a finding.

    ``bound_X`` is the 5th percentile of the class's own concordant members, so a class with
    no low tail at all still loses about 5 % of them.  A measured demotion rate near 5 % is
    therefore the floor of the instrument; only a rate well above it says something about the
    cells.  Exocrine's measured 5.2-5.4 % on the cohort sits on this floor.
    """
    Z, lab, st = matrix({"Immune": [(2_000, 3.0, "")], "Exocrine": [(2_000, 2.0, "")]})
    got = refine.apply_once(Z, MARKERS, lab, st, np.zeros(len(lab), bool))
    rate = got["demoted"][lab == "Immune"].mean()
    assert 0.04 < rate < 0.07


def test_the_bound_comes_from_the_class_own_concordant_members():
    """Not a constant: two classes at different levels get different bounds."""
    Z, lab, st = matrix({
        "Beta": [(400, 4.0, ">10")],
        "Immune": [(400, 1.0, "")],
        "Exocrine": [(400, 2.0, "")],
    })
    conc = refine.concordance(refine.neighbours(Z), lab)
    bounds = refine.fit_bounds(Z, MARKERS, lab, st, conc)
    assert bounds[("Beta", ">10")]["bound_sd"] > bounds[("Exocrine", "")]["bound_sd"]
    assert bounds[("Exocrine", "")]["bound_sd"] > bounds[("Immune", "")]["bound_sd"]


# ------------------------------------------------- per donor AND per size stratum


def test_the_bound_is_fitted_within_a_stratum_not_across_them():
    """A stratum whose members sit lower gets a lower bound, so its cells are not cut by the
    other stratum's level.  This is the 2026-09-18 ruling: reference levels per donor AND per
    size stratum, and the bound is a reference level."""
    Z, lab, st = matrix({
        "Delta": [(400, 5.0, ">10"), (400, 2.0, "4-10")],
        "Exocrine": [(400, 1.0, "")],
    })
    conc = refine.concordance(refine.neighbours(Z), lab)
    per_stratum = refine.fit_bounds(Z, MARKERS, lab, st, conc)
    pooled = refine.fit_bounds(Z, MARKERS, lab, np.full(len(lab), "", dtype=object), conc)
    assert per_stratum[("Delta", "4-10")]["bound_sd"] < per_stratum[("Delta", ">10")]["bound_sd"]
    demoted_per_stratum, _ = refine._demote(Z, MARKERS, lab, st, per_stratum)
    demoted_pooled, _ = refine._demote(Z, MARKERS, lab,
                                      np.full(len(lab), "", dtype=object), pooled)
    # Pooling the two strata puts the whole low stratum under one level and cuts into it.
    assert demoted_pooled[400:800].sum() > demoted_per_stratum[400:800].sum()


def test_an_underpowered_stratum_gets_no_bound_and_nothing_is_demoted_in_it():
    Z, lab, st = matrix({
        "Alpha": [(400, 4.0, ">10"), (refine.MIN_CLEAR - 10, 0.0, "1")],
        "Exocrine": [(400, 1.0, "")],
    })
    conc = refine.concordance(refine.neighbours(Z), lab)
    bounds = refine.fit_bounds(Z, MARKERS, lab, st, conc)
    assert bounds[("Alpha", "1")]["bound_sd"] is None
    demoted, why = refine._demote(Z, MARKERS, lab, st, bounds)
    small = st == "1"
    assert not demoted[small].any()
    assert all("bound unavailable" in w for w in why[small])


def test_an_underpowered_stratum_never_borrows_another_strata_bound():
    """The cells in the unfittable stratum would every one of them fail the other stratum's
    bound.  Borrowing it would demote them; the rule refuses instead, and says so."""
    Z, lab, st = matrix({
        "Alpha": [(400, 4.0, ">10"), (refine.MIN_CLEAR - 10, 0.0, "1")],
        "Exocrine": [(400, 1.0, "")],
    })
    conc = refine.concordance(refine.neighbours(Z), lab)
    bounds = refine.fit_bounds(Z, MARKERS, lab, st, conc)
    other = bounds[("Alpha", ">10")]["bound_sd"]
    own_max = refine.own_marker_max(Z, MARKERS, "Alpha")
    assert (own_max[st == "1"] < other).all()
    demoted, _ = refine._demote(Z, MARKERS, lab, st, bounds)
    assert not demoted[st == "1"].any()


# ---------------------------------------------------------------- the single pass


def test_apply_once_refuses_to_run_on_its_own_output():
    Z, lab, st = matrix({
        "Exocrine": [(600, 2.0, ""), (24, -2.0, "", "Immune")],
        "Immune": [(200, 3.0, "")],
    })
    once = refine.apply_once(Z, MARKERS, lab, st, np.zeros(len(lab), bool))
    with pytest.raises(refine.RefinementAlreadyApplied):
        refine.apply_once(Z, MARKERS, once["labels"], st, np.zeros(len(lab), bool))


def test_refine_donor_refuses_a_frame_that_has_already_been_refined():
    frame = cells_frame([("a", "Exocrine", "tissue", None, True)])
    frame["refined_class"] = "Exocrine"
    with pytest.raises(refine.RefinementAlreadyApplied):
        refine._refuse_second_pass(frame)


def test_refine_donor_refuses_a_frame_carrying_the_unresolved_label():
    frame = cells_frame([("a", refine.UNRESOLVED, "tissue", None, True)])
    with pytest.raises(refine.RefinementAlreadyApplied):
        refine._refuse_second_pass(frame)


def test_nothing_in_the_module_takes_a_step_or_iteration_argument():
    """The recovered run took ``PREV STEP`` on the command line and that is how it got looped.
    Step 2 was measured and rejected; there is no argument here that could raise it."""
    import inspect

    for fn in (refine.refine_donor, refine.apply_once, refine.main):
        names = set(inspect.signature(fn).parameters)
        assert not names & {"step", "steps", "passes", "iterations", "n_iter", "prev"}


# ------------------------------------------------- measurements are never invented


def test_a_missing_measurement_is_not_turned_into_a_number():
    params = pd.DataFrame({"cofactor": [1.0], "centre": [0.0], "spread": [1.0]},
                          index=pd.Index(["INS"], name="marker"))
    with pytest.raises(ValueError):
        refine.to_Z(np.array([[np.nan]], dtype="float32"), ["INS"], params)


def test_a_marker_with_no_representation_parameters_is_refused():
    params = pd.DataFrame({"cofactor": [1.0], "centre": [0.0], "spread": [1.0]},
                          index=pd.Index(["INS"], name="marker"))
    with pytest.raises(KeyError):
        refine.to_Z(np.array([[1.0, 2.0]], dtype="float32"), ["INS", "GCG"], params)


def test_a_class_whose_own_markers_the_donor_lacks_is_never_demoted():
    Z, lab, st = matrix({"Nerve": [(200, 3.0, "")], "Exocrine": [(400, 1.0, "")]})
    without = [m for m in MARKERS if m not in refine.OWN_MARKERS["Nerve"]]
    keep = [MARKERS.index(m) for m in without]
    conc = refine.concordance(refine.neighbours(Z[:, keep]), lab)
    bounds = refine.fit_bounds(Z[:, keep], without, lab, st, conc)
    assert ("Nerve", "") not in bounds
    demoted, why = refine._demote(Z[:, keep], without, lab, st, bounds)
    assert not demoted[lab == "Nerve"].any()
    assert set(why[lab == "Nerve"]) == {"no own marker on this donor"}


# --------------------------------------------------------------------- mechanics


def test_concordance_is_the_share_of_neighbours_with_the_same_label():
    labels = np.array(["A", "A", "B", "B"], dtype=object)
    nbr = np.array([[1, 2], [0, 3], [3, 0], [2, 1]])
    assert refine.concordance(nbr, labels).tolist() == [0.5, 0.5, 0.5, 0.5]


def test_the_neighbour_graph_is_reproducible():
    rng = np.random.default_rng(1)
    Z = rng.normal(size=(2_000, 8)).astype("float32")
    assert (refine.neighbours(Z) == refine.neighbours(Z)).all()


def test_concordance_is_recomputed_on_the_new_labels():
    Z, lab, st = matrix({
        "Blood vessel": [(400, 3.0, ""), (16, -1.0, "", "Exocrine")],
        "Exocrine": [(600, 2.0, "")],
    })
    got = refine.apply_once(Z, MARKERS, lab, st, np.zeros(len(lab), bool))
    assert not np.allclose(got["concordance"], got["concordance_before"])
    hand = refine.concordance(got["neighbours"], got["labels"])
    assert np.allclose(got["concordance"], hand)


# --------------------------------------------------------- strata come from the typing


def test_size_stratum_uses_the_one_authority():
    import subpanel_io

    assert refine.S.size_stratum is subpanel_io.size_stratum


def test_islet_size_counts_only_the_cells_the_analysis_admits(tmp_path, monkeypatch):
    """Structural classes and Excluded cells are in the islet but do not count toward its
    size, exactly as ``islet_identity.load`` counts it."""
    rows = [(f"u{i}", "Beta", "core", "islet-1", True) for i in range(4)]
    rows += [("v0", "Blood vessel", "core", "islet-1", True),
             ("v1", "Excluded", "core", "islet-1", True),
             ("w0", "Exocrine", "tissue", None, True)]
    path = tmp_path / "donor_id=X"
    path.mkdir()
    pq.write_table(pa.Table.from_pandas(cells_frame(rows), preserve_index=False),
                   path / "data.parquet")
    monkeypatch.setattr(refine, "CELLS", tmp_path)
    got = refine.load_cells("X").set_index("object_id")
    assert got.loc["u0", "islet_cells"] == 4
    assert got.loc["u0", "size_stratum"] == "4-10"
    assert got.loc["v0", "islet_cells"] == 4          # in the islet, not counted by it
    assert got.loc["w0", "size_stratum"] == ""        # outside every islet: no stratum
    assert np.isnan(got.loc["w0", "islet_cells"])


def test_a_cell_outside_an_islet_gets_no_stratum_rather_than_the_smallest(tmp_path,
                                                                          monkeypatch):
    rows = [(f"u{i}", "Exocrine", "tissue", None, True) for i in range(3)]
    path = tmp_path / "donor_id=X"
    path.mkdir()
    pq.write_table(pa.Table.from_pandas(cells_frame(rows), preserve_index=False),
                   path / "data.parquet")
    monkeypatch.setattr(refine, "CELLS", tmp_path)
    got = refine.load_cells("X")
    assert set(got["size_stratum"]) == {""}


# ------------------------------------------------------------- the per-donor frame


def _stub_donor(tmp_path, monkeypatch, rows, values, leaf=()):
    path = tmp_path / "donor_id=X"
    path.mkdir(exist_ok=True)
    pq.write_table(pa.Table.from_pandas(cells_frame(rows), preserve_index=False),
                   path / "data.parquet")
    monkeypatch.setattr(refine, "CELLS", tmp_path)
    monkeypatch.setattr(refine, "donor_values", lambda d, m: (values, MARKERS))
    monkeypatch.setattr(refine, "immune_leaf",
                        lambda d, m=None, require_min=0.1: pd.DataFrame(
                            {"object_id": [r[0] for r in rows],
                             "immune_leaf": [r[0] in set(leaf) for r in rows]}))
    return pd.DataFrame({"cofactor": 1.0, "centre": 0.0, "spread": 1.0},
                        index=pd.Index(MARKERS, name="marker"))


def test_object_id_is_preserved_verbatim_and_every_cell_comes_back(tmp_path, monkeypatch):
    rng = np.random.default_rng(3)
    n_clear, n_dim = 400, 16
    rows = [(f"bv{i}", "Blood vessel", "tissue", None, True) for i in range(n_clear + n_dim)]
    rows += [(f"ex{i}", "Exocrine", "tissue", None, True) for i in range(400)]
    Z, _, _ = matrix({"Blood vessel": [(n_clear, 3.0, ""), (n_dim, -1.0, "", "Exocrine")],
                      "Exocrine": [(400, 2.0, "")]}, seed=3)
    values = pd.DataFrame(np.sinh(Z).astype("float64"), columns=MARKERS)
    values.insert(0, "object_id", [r[0] for r in rows])
    values["qc_analysis_eligible"] = True
    params = _stub_donor(tmp_path, monkeypatch, rows, values)

    out, log = refine.refine_donor("X", params, MARKERS)
    assert out["object_id"].tolist() == [r[0] for r in rows]
    assert len(out) == len(rows) == log["n_cells"]
    assert log["n_measured"] == len(rows)
    dim = out["object_id"].isin([f"bv{i}" for i in range(n_clear, n_clear + n_dim)])
    assert out.loc[dim, "demoted"].all()
    assert set(out.loc[dim, "refined_class"]) == {"Exocrine"}
    assert log["iterated"] is False


def test_an_excluded_cell_is_left_alone_and_is_nobodys_neighbour(tmp_path, monkeypatch):
    rows = [(f"bv{i}", "Blood vessel", "tissue", None, True) for i in range(300)]
    rows += [(f"ex{i}", "Exocrine", "tissue", None, True) for i in range(300)]
    rows += [(f"xx{i}", "Excluded", "core", "islet-1", True) for i in range(80)]
    Z, _, _ = matrix({"Blood vessel": [(300, 3.0, "")],
                      "Exocrine": [(300, 2.0, ""), (80, 2.0, "")]}, seed=4)
    values = pd.DataFrame(np.sinh(Z).astype("float64"), columns=MARKERS)
    values.insert(0, "object_id", [r[0] for r in rows])
    values["qc_analysis_eligible"] = True
    params = _stub_donor(tmp_path, monkeypatch, rows, values)

    out, log = refine.refine_donor("X", params, MARKERS)
    excluded = out["cell_class"] == "Excluded"
    assert excluded.sum() == 80
    assert (out.loc[excluded, "refined_class"] == "Excluded").all()
    assert not out.loc[excluded, "demoted"].any()
    assert out.loc[excluded, "concordance"].isna().all()
    assert out.loc[excluded, "refine_pass"].eq(0).all()
    assert log["n_measured"] == len(rows) - 80
    assert "Excluded" not in {r["cell_class"] for r in log["after"]}


def test_a_cell_with_no_aligned_value_keeps_its_class_and_says_why(tmp_path, monkeypatch):
    rows = [(f"bv{i}", "Blood vessel", "tissue", None, True) for i in range(300)]
    rows += [(f"ex{i}", "Exocrine", "tissue", None, True) for i in range(300)]
    Z, _, _ = matrix({"Blood vessel": [(300, 3.0, "")], "Exocrine": [(300, 2.0, "")]}, seed=5)
    values = pd.DataFrame(np.sinh(Z).astype("float64"), columns=MARKERS)
    values.insert(0, "object_id", [r[0] for r in rows])
    values["qc_analysis_eligible"] = True
    values.loc[values["object_id"] == "bv0", "qc_analysis_eligible"] = False
    params = _stub_donor(tmp_path, monkeypatch, rows, values)

    out, log = refine.refine_donor("X", params, MARKERS)
    one = out[out["object_id"] == "bv0"].iloc[0]
    assert one["refined_class"] == "Blood vessel"
    assert not one["demoted"]
    assert np.isnan(one["concordance"])
    assert "not measured" in one["refine_note"]
    assert log["n_not_measured"] == 1


def test_the_log_records_where_every_demoted_cell_went(tmp_path, monkeypatch):
    rows = [(f"bv{i}", "Blood vessel", "tissue", None, True) for i in range(416)]
    rows += [(f"ex{i}", "Exocrine", "tissue", None, True) for i in range(400)]
    Z, _, _ = matrix({"Blood vessel": [(400, 3.0, ""), (16, -1.0, "", "Exocrine")],
                      "Exocrine": [(400, 2.0, "")]}, seed=6)
    values = pd.DataFrame(np.sinh(Z).astype("float64"), columns=MARKERS)
    values.insert(0, "object_id", [r[0] for r in rows])
    values["qc_analysis_eligible"] = True
    leaf = [f"bv{i}" for i in range(400, 408)]
    params = _stub_donor(tmp_path, monkeypatch, rows, values, leaf=leaf)

    out, log = refine.refine_donor("X", params, MARKERS)
    routed = log["routed"]
    dim = out["object_id"].isin([f"bv{i}" for i in range(400, 416)])
    assert out.loc[dim, "demoted"].all()
    # Eight of the sixteen carry immune evidence and go to Immune; every other demoted cell is
    # accounted for, and a demoted Exocrine cell becomes Unresolved rather than Immune.
    assert routed["Immune"] == 8
    assert routed["Unresolved"] == int((out["demoted"] & (out["cell_class"] == "Exocrine")).sum())
    assert log["n_demoted"] == sum(routed.values()) == int(out["demoted"].sum())
    assert set(out.loc[dim & out["object_id"].isin(leaf), "refined_class"]) == {"Immune"}
    assert not (out["demoted"] & (out["cell_class"] == "Exocrine")
                & (out["refined_class"] == "Immune")).any()
