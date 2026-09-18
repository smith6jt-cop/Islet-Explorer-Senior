"""pLN immune populations: the method transfers from the pancreas, and no number does.

Almost everything here is a statement about *definitions*, because the risk this build carries is
not an arithmetic slip -- it is a marker quietly joining or leaving the set that decides what
"immune" means, or a level creeping back in.  So the tests assert the sets against their
authorities (`lineages`, `core/panels/data/identity_core`, `subpanel`) rather than against a
frozen copy.

THE LOAD-BEARING TEST IS :func:`test_no_pancreas_number_and_no_reference_level_reaches_a_pln_call`.
Two generations of this lane have been withdrawn:

  2026-09-17 (a)  each donor's PANCREAS p99 was carried onto their pLN slide.  The transferred
                  CD66 level left 1.9 % of the slide beneath it and on 12 donors sat below that
                  slide's own 5th percentile.
  2026-09-17 (b)  the level was then fitted on the pLN itself, as the p99 of `context ==
                  'tissue'` -- which on this export is the whole slide, so every marker called
                  1 % of the slide positive BY CONSTRUCTION.  There is no context column and no
                  reference population on a node.

That test reads the module's SYNTAX TREE rather than its text, so the prose may go on explaining
either withdrawal without weakening the assertion, and no comment can satisfy it.

The tests that touch the 5.1M-row artifact skip when it is not on disk.

Run with
    PYTHONDONTWRITEBYTECODE=1 .venv/bin/python -m pytest -p no:cacheprovider \
      core/phenocycler/tests/current/test_pln_immune_pheno.py -q
"""

from __future__ import annotations

import ast
import importlib.util
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

REPO = Path(__file__).resolve().parents[4]
SCRIPT = REPO / "core/phenocycler/scripts/typing/build_pln_immune_pheno.py"
sys.path.insert(0, str(REPO / "core/phenocycler/scripts/typing"))
SPEC = importlib.util.spec_from_file_location("build_pln_immune_pheno", SCRIPT)
builder = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(builder)

import lineages as LIN       # noqa: E402
import subpanel as SP        # noqa: E402

OUT = REPO / "data/development/typing/pln/pheno"
CLUSTERS = REPO / "data/development/typing/pln/clusters/phenocycler"
IDENTITY_CORE = REPO / "core/panels/data/identity_core"


def _code_symbols() -> tuple[set[str], set[str]]:
    """Every name the module's CODE uses, and every string constant that is not a docstring.

    Reading the syntax tree instead of the file text is what makes this assertion mean something:
    a docstring or a comment explaining the withdrawal cannot satisfy it, and a call cannot hide
    behind one.
    """
    tree = ast.parse(SCRIPT.read_text())
    docstrings = set()
    for node in ast.walk(tree):
        if isinstance(node, (ast.Module, ast.FunctionDef, ast.AsyncFunctionDef, ast.ClassDef)):
            body = getattr(node, "body", [])
            if body and isinstance(body[0], ast.Expr) and isinstance(body[0].value, ast.Constant) \
                    and isinstance(body[0].value.value, str):
                docstrings.add(id(body[0].value))
    names: set[str] = set()
    strings: set[str] = set()
    for node in ast.walk(tree):
        if isinstance(node, ast.Name):
            names.add(node.id)
        elif isinstance(node, ast.Attribute):
            names.add(node.attr)
        elif isinstance(node, ast.alias):
            names.add(node.name.split(".")[-1])
            names.add(node.name)
        elif isinstance(node, ast.ImportFrom) and node.module:
            names.add(node.module)
        elif isinstance(node, ast.Constant) and isinstance(node.value, str) \
                and id(node) not in docstrings:
            strings.add(node.value)
    return names, strings


# --------------------------------------------------------------------------------------------
# No pancreas number, and no reference level of any kind, reaches a pLN call.
# --------------------------------------------------------------------------------------------
def test_no_pancreas_number_and_no_reference_level_reaches_a_pln_call():
    """Four independent ways, because any one of them can be edited around.

    1. The module's code names no pancreas accessor, no pancreas tree and no level-fitting
       function -- and `subpanel_io`, which owns every one of them, is not imported at all.
    2. The withdrawn functions and constants are gone from the module, not merely unused.
    3. No string constant in the code points at a pancreas tree.
    4. The published root holds none of the withdrawn outputs.
    """
    names, strings = _code_symbols()

    # (1) the lane that owns the pancreas levels is not reachable from here.
    assert "subpanel_io" not in names, "subpanel_io owns the pancreas aligned tree and its scale"
    for forbidden in ("fit_levels",          # fits a level from a reference population
                      "reference_level",     # the p99 of that population
                      "REFERENCE_QUANTILE",
                      "IMMUNE_REQUIRE_MIN",  # an absolute floor defined against such a level
                      "fires",               # compares values against that floor
                      "cell_source",         # the pancreas geometry accessor
                      "ALIGNED", "SCALE_CSV", "load_scale", "apply_scale", "RUN"):
        assert forbidden not in names, \
            f"{forbidden} is part of the withdrawn reference-level lane and must not run here"

    # `context` was the pancreas's partition of the section.  A node has none, and this lane must
    # not read one, synthesise one, or fit anything against one.
    assert "context" not in names
    assert "pln_context" not in names
    for column in ("context", "qupath_classification", "parent_raw"):
        assert column not in strings, \
            f"the export column {column!r} is read nowhere here: tissue does not mean anything"

    # (2) the withdrawn functions are gone from the module object.
    for gone in ("pancreas_levels", "transfer_check", "acquisition_factors", "fit_pln_levels",
                 "in_fires_units", "pln_context", "stage1_decide", "QUANTILE", "POSITIVE",
                 "NO_PANCREAS_LEVEL", "K_LOW", "K_HIGH", "FLOOR_QUANTILE"):
        assert not hasattr(builder, gone), f"{gone} belongs to a withdrawn generation"

    # (3) no string points at a pancreas tree.
    for fragment in ("typing/cells", "data/runs", "qupath_manifest", "aligned", "cells_redsea",
                     "restore"):
        assert not any(fragment in s for s in strings), f"{fragment} is a pancreas input"

    # (4) nothing of either withdrawn generation is presented as a current output.
    for stale in ("pancreas_reference_levels.parquet", "level_transfer_check.parquet",
                  "channel_acquisition.parquet"):
        assert not list(OUT.rglob(stale)), f"{stale} is an output of the withdrawn transfer"
    if (OUT / "manifest.json").exists():
        current = set(json.loads((OUT / "_reference/SUPERSEDED_2026-09-17.json").read_text())
                      ["current_outputs"])
        assert not any("level" in name for name in current)


def test_the_source_of_every_number_is_the_pln_section_itself():
    """The only cell tree this lane reads is `cells_pln`, and the only labels are the clusters."""
    _, strings = _code_symbols()
    roots = [s for s in strings if "data/" in s and s.endswith(("cells_pln", "phenocycler"))]
    assert "core/phenocycler/data/integration/cells_pln" in str(builder.PLN_CELLS)
    assert str(builder.CLUSTERS).endswith("data/development/typing/pln/clusters/phenocycler")
    assert roots, "the roots must be named in code, not assembled from fragments"


def test_no_donor_exclusion_list_is_restated_in_this_lane():
    """The cohort authority and the tracked pLN manifest decide; this file must not repeat them."""
    text = SCRIPT.read_text()
    for donor in ("6457", "6579"):
        assert f'"{donor}"' not in text and f"'{donor}'" not in text


# --------------------------------------------------------------------------------------------
# The statistic: a cluster against the other cells of its own section.
# --------------------------------------------------------------------------------------------
def _standout(values, codes, k, markers, offset_over=None):
    return builder.standout(pd.DataFrame(values), np.asarray(codes), k, markers,
                            offset_over=offset_over)


def test_auc_is_the_probability_a_cluster_cell_reads_higher_than_a_cell_elsewhere():
    """A cluster that holds the section's top half reads AUC 1.0 against the bottom half."""
    values = {"CD3e": np.arange(100, dtype=float), "SMA": np.arange(100, dtype=float)[::-1]}
    codes = np.array([0] * 50 + [1] * 50)
    st = _standout(values, codes, 2, ["CD3e", "SMA"])
    auc = st["auc"]
    assert auc[1, 0] == pytest.approx(1.0)      # cluster 1 holds CD3e 50..99
    assert auc[0, 0] == pytest.approx(0.0)
    assert auc[0, 1] == pytest.approx(1.0)      # and the bottom of SMA, which is its top


def test_a_cluster_indistinguishable_from_the_rest_of_its_section_sits_at_the_null():
    values = {"CD3e": np.tile([1.0, 1.0, 2.0, 2.0], 50)}
    codes = np.tile([0, 1], 100)          # each cluster then holds the same multiset
    st = _standout(values, codes, 2, ["CD3e"])
    assert st["auc"][0, 0] == pytest.approx(builder.AUC_NULL, abs=1e-9)
    assert not st["enriched"].any(), "an indistinguishable cluster is not enriched"


def test_no_absolute_threshold_survives_a_change_of_units():
    """The statistic reads ranks, so multiplying a channel by 1000 cannot move a single call.

    This is what "no level that has to mean the same thing on two slides" buys: two sections
    acquired at different exposures give the same answer for the same cell ordering.
    """
    values = {"CD3e": np.arange(200, dtype=float), "SMA": np.random.default_rng(0).random(200)}
    codes = np.array([0] * 100 + [1] * 100)
    a = _standout(values, codes, 2, ["CD3e", "SMA"])
    scaled = {"CD3e": values["CD3e"] * 1000.0 + 5.0, "SMA": values["SMA"] * 7.0}
    b = _standout(scaled, codes, 2, ["CD3e", "SMA"])
    assert np.allclose(a["auc"], b["auc"], equal_nan=True)
    assert (a["enriched"] == b["enriched"]).all()


def test_a_missing_measurement_is_unavailable_not_zero_and_not_the_bottom_rank():
    """A cell without a measurement is ranked among nobody; it does not become the dimmest cell."""
    v = np.arange(100, dtype=float)
    v[:50] = np.nan                      # cluster 0 is unmeasured, cluster 1 carries 50..99
    st = _standout({"CD3e": v}, np.array([0] * 50 + [1] * 50), 2, ["CD3e"])
    assert np.isnan(st["auc"][0, 0]), "an unmeasured cluster has no AUC, not AUC 0"
    assert st["n_finite"][0, 0] == 0
    assert not st["enriched"][0, 0]
    assert np.isnan(st["unit_rank"]["CD3e"][:50]).all()


def test_the_panel_offset_is_a_median_over_every_channel_and_cancels_in_a_comparison():
    """S removes a cluster's panel-wide offset -- and, stated because it would be easy to claim
    otherwise, it cannot change which of two markers of the SAME cluster is the larger."""
    rng = np.random.default_rng(1)
    n = 400
    codes = np.array([0] * 200 + [1] * 200)
    values = {m: rng.random(n) for m in ("CD3e", "SMA", "CD68", "CD20")}
    for m in values:                     # cluster 1 reads high on EVERYTHING
        values[m][200:] += 3.0
    st = _standout(values, codes, 2, list(values), offset_over=list(values))
    assert st["offset"][1] > 0.9 and st["offset"][0] < 0.1
    order_auc = np.argsort(st["auc"][1])
    order_rel = np.argsort(st["rel"][1])
    assert (order_auc == order_rel).all(), "a per-cluster constant must cancel in a comparison"
    assert np.allclose(st["rel"][1], st["auc"][1] - st["offset"][1])


def test_enrichment_needs_all_three_conjuncts():
    """Direction, an internal null, and being above the cluster's own panel-wide level."""
    rng = np.random.default_rng(2)
    n = 600
    codes = np.array([0] * 300 + [1] * 300)
    values = {m: rng.random(n) for m in ("CD3e", "SMA", "CD68", "CD20", "CD31")}
    for m in values:
        values[m][300:] += 5.0           # cluster 1 is uniformly bright: nothing stands out FOR it
    st = _standout(values, codes, 2, list(values), offset_over=list(values))
    assert (st["auc"][1] > builder.AUC_NULL).all()          # direction holds for every channel
    assert (st["p"][1] < st["alpha"]).all()                 # and so does the internal null
    assert st["enriched"][1].sum() <= len(values) // 2, \
        "a cluster that reads high on everything must not be enriched for everything"


# --------------------------------------------------------------------------------------------
# Stage 1: what "immune" means here, and where it came from.
# --------------------------------------------------------------------------------------------
def test_primaries_are_lineages_minus_supporting_only_and_nothing_else():
    """No marker is promoted by hand.  `lineages` is the authority and this is a projection of it."""
    for lineage, markers in builder.PRIMARY_IMMUNE.items():
        expected = tuple(m for m in LIN.LINEAGE_MARKERS[lineage] if m not in LIN.SUPPORTING_ONLY)
        assert markers == expected, f"{lineage} diverged from lineages.LINEAGE_MARKERS"
    assert set(builder.IMMUNE_LINEAGES) <= set(LIN.LINEAGE_MARKERS)


def test_supporting_only_markers_never_carry_a_call_alone():
    primaries = {m for markers in builder.PRIMARY_IMMUNE.values() for m in markers}
    assert not (primaries & LIN.SUPPORTING_ONLY)
    assert set(builder.SUPPORTING_IMMUNE) <= LIN.SUPPORTING_ONLY
    assert not (set(builder.STRUCTURAL) & LIN.SUPPORTING_ONLY)


def test_there_is_no_pan_leukocyte_gate_to_fall_back_on():
    """CD45/PTPRC is not on this panel, which is why stage 1 is a union and not one gate."""
    every = (set(builder.STRUCTURAL) | set(builder.SUPPORTING_IMMUNE)
             | {m for v in builder.PRIMARY_IMMUNE.values() for m in v})
    assert not ({"CD45", "PTPRC"} & every)


def test_structural_set_excludes_immune_primaries_and_names_what_it_dropped():
    primaries = {m for markers in builder.PRIMARY_IMMUNE.values() for m in markers}
    assert not (set(builder.STRUCTURAL) & primaries)
    assert "Vimentin" in builder.STRUCTURAL_DROPPED
    assert "leukocyte" in builder.STRUCTURAL_DROPPED["Vimentin"]
    assert "Vimentin" not in builder.STRUCTURAL


def test_a_channel_that_names_no_lineage_is_refused_on_the_taxonomy_and_not_on_a_threshold():
    named = set(LIN.identity_markers())
    for marker in builder.NOT_A_NAMING_MARKER:
        assert marker not in named, f"{marker} does name a lineage in `lineages`"
    assert {"CD14", "CD21", "HLA_DR", "CD57"} <= set(builder.NOT_A_NAMING_MARKER)


# --------------------------------------------------------------------------------------------
# Stage 1 and stage 2: the decision itself, on clusters.
# --------------------------------------------------------------------------------------------
def _call(rows, markers, members=None, magnitude=None):
    """Drive `call_clusters` from a hand-written (cluster x marker) standout table.

    A magnitude has to be supplied because the arbitration reads it.  Unless a test writes one
    itself, it is `1 + S`: a magnitude is a ratio of medians so it is positive, and `1` is the
    value it takes for a cluster that reads exactly the section's median, which is what S 0 says
    too -- so the hand-written tables order the same way under both statistics and a test that
    means to exercise the ORDER is not quietly exercising the change of statistic."""
    rel = np.array([[r.get(m, np.nan) for m in markers] for r in rows], dtype=float)
    enr = np.isfinite(rel) & (rel > 0)
    mag = (1.0 + rel) if magnitude is None else np.array(
        [[r.get(m, np.nan) for m in markers] for r in magnitude], dtype=float)
    st = {"markers": markers, "auc": rel + 0.5, "rel": rel, "enriched": enr,
          "magnitude": mag, "offset": np.zeros(len(rows))}
    members = members or builder.subpanel_membership(markers)
    ids = [f"c{i}" for i in range(len(rows))]
    return builder.call_clusters("T", ids, np.full(len(rows), 100), st, members)


MARKERS = ["CD3e", "CD4", "CD8", "CD20", "CD79a", "CD68", "Iba1", "CD163", "CD206", "CD11c",
           "CD209", "MPO", "CD66", "CD56", "CD57", "Granzyme_B", "FOXP3", "M2Gal3", "CD11b",
           "CD31", "CD34", "SMA", "Caveolin", "Collagen_IV", "E_cadherin", "EpCAM",
           "Pan_Cytokeratin", "Podoplanin", "B3TUBB"]


def test_a_cluster_is_immune_when_its_canonical_stands_out_over_the_structural_evidence():
    calls = _call([{"CD3e": 0.4, "SMA": 0.1}], MARKERS)
    assert bool(calls.loc[0, "is_immune"]) is True
    assert calls.loc[0, "stage1_rule"] == "immune_canonical_carries_more_magnitude_than_structural"
    assert calls.loc[0, "immune_class"] == "T_cell"


def test_a_cluster_whose_structural_evidence_wins_is_not_immune_and_gets_no_lineage():
    calls = _call([{"CD3e": 0.1, "Caveolin": 0.5}], MARKERS)
    assert bool(calls.loc[0, "is_immune"]) is False
    assert calls.loc[0, "stage1_rule"] == "structural_carries_more_magnitude_than_immune"
    # This lane names no non-immune lineage at all: the structural markers arbitrate and never
    # become a label.  A node has no exocrine and the panel could not split one anyway.
    assert calls.loc[0, "immune_class"] == ""


def test_a_cluster_where_nothing_stands_out_is_uncalled_and_is_never_folded_into_a_class():
    calls = _call([{"CD3e": -0.2, "SMA": -0.3}], MARKERS)
    assert calls.loc[0, "is_immune"] is None
    assert calls.loc[0, "stage1_rule"] == "no_lineage_marker_stands_out_on_this_section"
    assert calls.loc[0, "immune_class"] == ""


# --------------------------------------------------------------------------------------------
# The arbitration.  Corrected 2026-09-18: the rank statistic admits the evidence, the MAGNITUDE
# decides between two channels of one cluster.  These tests are about the mechanism; none of them
# passes by reproducing a count, a share or any previous result.
# --------------------------------------------------------------------------------------------
def _section(seed=7, n=600):
    """Two clusters on one synthetic section, built so that the two statistics DISAGREE.

    `CD3e` separates cluster 0 from cluster 1 by a factor of three, with the two distributions
    overlapping.  `EpCAM` is background on both -- 3 % apart -- but it is background with no
    overlap at all, so its rank separation is the CLEANER of the two.  That is the whole of the
    defect in one section: a channel with no population outranks a channel with one, because a
    rank cannot see that 1.03x and 3x are different things.
    """
    rng = np.random.default_rng(seed)
    codes = np.array([0] * n + [1] * n)
    values = {
        "CD3e": np.concatenate([rng.normal(30.0, 9.0, n), rng.normal(10.0, 3.0, n)]),
        "EpCAM": np.concatenate([rng.uniform(10.3, 10.4, n), rng.uniform(10.0, 10.1, n)]),
    }
    for filler in ("CD68", "CD20", "SMA", "CD31", "Caveolin", "CD34"):
        values[filler] = rng.uniform(1.0, 2.0, 2 * n)
    return values, codes, list(values)


def test_a_rank_and_a_magnitude_disagree_and_the_arbitration_reads_the_magnitude():
    values, codes, markers = _section()
    st = _standout(values, codes, 2, markers, offset_over=markers)
    st["magnitude"] = builder.magnitude(pd.DataFrame(values), codes, 2, markers)
    j = {m: i for i, m in enumerate(markers)}
    # the disagreement is REAL on this section and is asserted rather than assumed
    assert st["auc"][0, j["EpCAM"]] > st["auc"][0, j["CD3e"]], "the rank must favour the background"
    assert st["magnitude"][0, j["CD3e"]] > st["magnitude"][0, j["EpCAM"]]
    assert st["enriched"][0, j["CD3e"]] and st["enriched"][0, j["EpCAM"]]

    calls = builder.call_clusters("T", ["c0", "c1"], np.array([len(codes) // 2] * 2), st,
                                  builder.subpanel_membership(markers))
    assert bool(calls.loc[0, "is_immune"]) is True
    assert calls.loc[0, "top_immune_marker"] == "CD3e"
    assert calls.loc[0, "stage1_rule"] == "immune_canonical_carries_more_magnitude_than_structural"
    # and the statistic that was read is on the row, beside the one that was not
    assert calls.loc[0, "magnitude_margin"] > 0 > calls.loc[0, "separation_margin"]


def test_the_panel_offset_never_protected_the_arbitration():
    """Stated because the generation before this one believed it did.  S is AUC minus a
    per-cluster constant, so it cannot change which of two channels of THAT cluster is larger --
    the arbitration was on raw AUC all along, and subtracting the offset was never going to fix
    it."""
    values, codes, markers = _section()
    st = _standout(values, codes, 2, markers, offset_over=markers)
    j = {m: i for i, m in enumerate(markers)}
    for row in (0, 1):
        assert (np.argsort(st["auc"][row]) == np.argsort(st["rel"][row])).all()
    assert st["rel"][0, j["EpCAM"]] > st["rel"][0, j["CD3e"]]


def test_the_magnitude_is_a_ratio_inside_one_channel_and_survives_a_change_of_scale():
    """Two sections acquired at different exposures must give the same answer, which is the one
    thing a level could not do.  Each channel is scaled by a DIFFERENT constant, so nothing is
    preserved by accident."""
    values, codes, markers = _section()
    mag = builder.magnitude(pd.DataFrame(values), codes, 2, markers)
    scaled = {m: v * (3.0 + 2.0 * i) for i, (m, v) in enumerate(values.items())}
    mag2 = builder.magnitude(pd.DataFrame(scaled), codes, 2, markers)
    assert np.allclose(mag, mag2, equal_nan=True)

    st = _standout(values, codes, 2, markers, offset_over=markers)
    st["magnitude"] = mag
    members = builder.subpanel_membership(markers)
    a = builder.call_clusters("T", ["c0", "c1"], np.array([600, 600]), st, members)
    st2 = _standout(scaled, codes, 2, markers, offset_over=markers)
    st2["magnitude"] = mag2
    b = builder.call_clusters("T", ["c0", "c1"], np.array([600, 600]), st2, members)
    for col in ("is_immune", "stage1_rule", "immune_class", "top_immune_marker",
                "top_structural_marker"):
        assert list(a[col]) == list(b[col]), col


def test_a_channel_with_no_positive_section_median_has_no_magnitude_and_cannot_arbitrate():
    """`Ker8_18` and `Keratin_5` are on the structural list and absent from this export.  A
    missing magnitude is unavailable -- never zero, which would make the channel the quietest on
    the section, and never large, which would let it decide."""
    values, codes, markers = _section()
    values["Ker8_18"] = np.zeros(len(codes))
    markers = list(values)
    mag = builder.magnitude(pd.DataFrame(values), codes, 2, markers)
    assert np.isnan(mag[:, markers.index("Ker8_18")]).all()
    st = _standout(values, codes, 2, markers, offset_over=markers)
    st["magnitude"] = mag
    calls = builder.call_clusters("T", ["c0", "c1"], np.array([600, 600]), st,
                                  builder.subpanel_membership(markers))
    assert "Ker8_18" not in set(calls["top_structural_marker"])


def test_dynamic_range_is_the_per_channel_maximum_of_the_magnitude():
    """One definition backs the number the veto accounting reports and the number the arbitration
    reads, so the two can never drift apart."""
    values, codes, markers = _section()
    mag = builder.magnitude(pd.DataFrame(values), codes, 2, markers)
    dyn = builder.dynamic_range(pd.DataFrame(values), codes, 2, markers)
    for i, m in enumerate(markers):
        assert dyn[m] == pytest.approx(float(np.nanmax(mag[:, i])))


# --------------------------------------------------------------------------------------------
# Stage 2: the winner is held to its own rule-table requirement.
# --------------------------------------------------------------------------------------------
def test_the_requirement_is_the_all_of_clause_of_the_lane_s_own_rule_table():
    """Nothing is written by hand: each subpanel is mapped to a phenotype by name and the markers
    come from `subpanel.immune_rules()`.  `Imm-macrophage` requires nothing, and an empty
    requirement stays empty rather than being filled in."""
    table = SP.immune_rules()
    need = builder.subpanel_requirements()
    for name, rule in builder.SUBPANEL_RULE.items():
        assert need[name] == tuple(table[rule][0])
    assert need["T_cell"] == ("CD3e",)
    assert need["Myeloid"] == ()
    assert "NK_ILC" not in builder.SUBPANEL_RULE, \
        "NK reaches stage 2 through nk_licence(), which reads the whole rule and not all_of alone"


def test_no_cluster_is_named_t_on_cd4_or_cd8_without_cd3e():
    """The rule table says `Imm-T-cell` is `all_of ['CD3e']` and every T subtype in it repeats
    CD3e.  A cluster whose T evidence is CD4 or CD8 alone is UNCALLED with the missing marker
    named, and it is never handed to the runner-up."""
    calls = _call([{"CD4": 0.5, "CD8": 0.4, "CD68": 0.1, "SMA": 0.05}], MARKERS)
    assert bool(calls.loc[0, "is_immune"]) is True
    assert calls.loc[0, "immune_class"] == "Uncalled"
    assert "CD3e" in calls.loc[0, "stage2_reason"]
    assert "not reassigned" in calls.loc[0, "stage2_reason"]
    # and with CD3e enriched the same cluster is named
    named = _call([{"CD3e": 0.2, "CD4": 0.5, "CD8": 0.4, "SMA": 0.05}], MARKERS)
    assert named.loc[0, "immune_class"] == "T_cell"


def test_the_epithelial_negatives_of_the_rule_table_are_not_carried_over_as_enrichment():
    """Every immune phenotype in that table lists the epithelial channels as `none_of`, evaluated
    there against a LEVEL.  This lane has no level, and "above this cluster's own median channel"
    is a different statement -- carrying it across is what put 429,278 cells of a lymph node
    behind an epithelial veto.  The requirement is carried; the negatives are not."""
    for rule in builder.SUBPANEL_RULE.values():
        assert "EpCAM" in SP.immune_rules()[rule][2]
    calls = _call([{"CD3e": 0.4, "EpCAM": 0.9, "SMA": 0.05}],
                  MARKERS, magnitude=[{"CD3e": 3.0, "EpCAM": 1.05, "SMA": 1.0}])
    assert bool(calls.loc[0, "is_immune"]) is True
    assert calls.loc[0, "immune_class"] == "T_cell"


def test_the_per_cell_evidence_reads_the_magnitude_and_not_the_rank():
    """A cross-channel comparison INSIDE one cell has the same defect one level down: a cell's
    place in a channel that is background everywhere is still a place."""
    values, codes, markers = _section()
    st = _standout(values, codes, 2, markers, offset_over=markers)
    ev = builder.cell_evidence(st, builder.subpanel_membership(markers), len(codes))
    assert set(ev) == {"immune_magnitude", "structural_magnitude", "top_immune_marker"}
    # cluster 0 reads 3x on CD3e and 1.03x on EpCAM, so its cells lean immune on the magnitude
    assert np.nanmean(ev["immune_magnitude"][:600] > ev["structural_magnitude"][:600]) > 0.9
    # on the rank the same cells would not, which is the measurement this replaces
    rank_imm = np.vstack([st["unit_rank"][m] for m in ("CD3e", "CD20", "CD68")]).max(axis=0)
    rank_str = np.vstack([st["unit_rank"][m] for m in ("EpCAM", "SMA", "CD31")]).max(axis=0)
    assert np.nanmean(rank_imm[:600] > rank_str[:600]) < 0.9


def test_the_winning_subpanel_must_carry_an_enriched_canonical_of_its_own_lineage():
    """The Xenium draft failed exactly here: 55.4 % of its Granulocyte calls detected none of the
    granulocyte canonicals.  A winner with no enriched canonical is UNCALLED, and it is NEVER
    reassigned to the runner-up -- reassignment is how a larger subpanel absorbs a cluster that
    belongs to nobody."""
    # CD209 is a stage-1 myeloid canonical but is not in 06_immune_myeloid.csv, so it carries no
    # stage-2 evidence; nothing else stands out, so the cluster is immune and unresolved.
    calls = _call([{"CD209": 0.5, "SMA": 0.1, "CD3e": -0.1, "CD68": -0.2}], MARKERS)
    assert bool(calls.loc[0, "is_immune"]) is True
    assert calls.loc[0, "immune_class"] == "Uncalled"
    assert "no enriched canonical" in calls.loc[0, "stage2_reason"]
    assert "not reassigned" in calls.loc[0, "stage2_reason"]


def test_a_subpanel_never_wins_on_a_marker_that_is_not_its_own():
    """CD68 is myeloid; it cannot carry a T-cell call however loud it is."""
    calls = _call([{"CD68": 0.6, "CD3e": 0.2, "SMA": 0.0}], MARKERS)
    assert calls.loc[0, "immune_class"] == "Myeloid"
    assert calls.loc[0, "marker_Myeloid"] == "CD68"


def test_nk_is_licensed_by_the_rule_table_and_not_by_a_hand_picked_partner():
    """CD56 is `SUPPORTING_ONLY` and is the NK subpanel's only member on this panel, so NK_ILC has
    no naming marker.  The licence is `subpanel.immune_rules()['Imm-NK-cell']` -- CD56 with (CD57
    or Granzyme B) and CD3e negative -- which is the lane's own statement, and it runs without any
    fact about the pancreas because both partners are on the pLN panel."""
    members = builder.subpanel_membership(MARKERS)
    assert set(members["NK_ILC"]) <= LIN.SUPPORTING_ONLY
    assert not builder.naming_markers(members)["NK_ILC"]
    all_of, any_of, none_of = builder.nk_licence()
    rule = SP.immune_rules()[builder.SUPPORTING_ONLY_SUBPANEL["NK_ILC"]]
    assert all_of == list(rule[0]) == ["CD56"]
    assert any_of == [["CD57", "Granzyme_B"]]
    assert "CD3e" in none_of

    alone = _call([{"CD56": 0.6, "SMA": 0.0}], MARKERS)
    assert not bool(alone.loc[0, "nk_licensed_by_rule_table"]), "CD56 must not speak alone"
    partnered = _call([{"CD68": 0.15, "CD56": 0.6, "Granzyme_B": 0.3, "SMA": 0.0}], MARKERS)
    assert bool(partnered.loc[0, "nk_licensed_by_rule_table"])
    assert partnered.loc[0, "immune_class"] == "NK_ILC"
    vetoed = _call([{"CD68": 0.15, "CD56": 0.6, "Granzyme_B": 0.3, "CD3e": 0.1, "SMA": 0.0}],
                   MARKERS)
    assert not bool(vetoed.loc[0, "nk_licensed_by_rule_table"]), "CD3e vetoes the NK licence"
    assert vetoed.loc[0, "immune_class"] != "NK_ILC"


def test_a_zero_nk_count_is_a_property_of_the_panel_and_never_a_measurement():
    """CD56 is `SUPPORTING_ONLY` -- shared with the Neural lineage -- so it is not a stage-1
    primary, and a cluster standing out on CD56 ALONE never enters the immune compartment at all.
    `lineages` says it directly: "A lineage whose evidence rests only on these is Ambiguous, not
    named."  NK_ILC is therefore reachable only through a cluster admitted on another lineage's
    canonical, and the manifest must say so, because a reader who finds class_NK_ILC = 0 across
    fourteen donors would otherwise read it as a measurement that a lymph node holds no NK cell."""
    assert "CD56" in LIN.SUPPORTING_ONLY
    assert "CD56" not in {m for v in builder.PRIMARY_IMMUNE.values() for m in v}
    only_cd56 = _call([{"CD56": 0.9, "Granzyme_B": 0.5, "SMA": -0.1}], MARKERS)
    assert only_cd56.loc[0, "is_immune"] is None, "CD56 alone must not admit a cluster"
    assert only_cd56.loc[0, "immune_class"] == ""
    if (OUT / "manifest.json").exists():
        note = json.loads((OUT / "manifest.json").read_text())["nk_ilc_reachability"]
        assert "NOT a measurement" in note


# --------------------------------------------------------------------------------------------
# Stage 2: membership comes from core/panels and nowhere else.
# --------------------------------------------------------------------------------------------
def test_subpanel_membership_is_backed_by_an_identity_core_csv_or_a_declared_surrogate():
    from core.phenocycler.phenocycler.integration.vocab import (
        PROTEIN_GENE_CANDIDATES, SURROGATE_PANEL)
    members = builder.subpanel_membership(MARKERS)
    assert set(members) == set(builder.SUBPANEL_FILES)
    for name, hints in members.items():
        stem = builder.SUBPANEL_FILES[name][:-4]
        genes = set(pd.read_csv(IDENTITY_CORE / f"{stem}.csv")["gene"])
        for protein in hints:
            routed = SURROGATE_PANEL.get(protein) == stem
            direct = any(g in genes for g in PROTEIN_GENE_CANDIDATES.get(protein, []))
            assert routed or direct, f"{protein} is in {name} on no authority"


def test_a_canonical_the_csv_does_not_list_is_not_patched_in_by_hand():
    """CD209 is a myeloid canonical in `lineages` and CD66 is a granulocyte one, and neither gene
    is in its subpanel's CSV.  core/panels is the authority, so neither is added here."""
    members = builder.subpanel_membership(MARKERS)
    assert "CD209" not in members["Myeloid"]
    assert "CD66" not in members["Granulocyte"]
    assert builder.naming_markers(members)["Granulocyte"] == ["MPO"]


def test_state_markers_carry_no_stage_two_evidence():
    members = builder.subpanel_membership(
        MARKERS + ["PD_1", "LAG3", "TCF_1", "ICOS", "VISTA"])
    named = {m for hints in members.values() for m in hints}
    assert not (named & {"PD_1", "LAG3", "TCF_1", "ICOS", "VISTA"})


def test_no_acinar_or_ductal_label_can_be_produced():
    """The 60-plex carries no amylase, trypsin, CPA1, elastase, CK19, CA19-9, SOX9, HNF1B, CFTR or
    MUC1.  This lane names no non-immune lineage at all, so the labels it can emit are the five
    immune subpanels, `Uncalled` and the empty string."""
    emittable = set(builder.SUBPANEL_FILES) | {"Uncalled", ""}
    assert not any(word in label.lower() for label in emittable
                   for word in ("acinar", "ductal", "exocrine", "islet", "beta", "alpha"))
    if (OUT / "cluster_calls.csv").exists():
        seen = set(pd.read_csv(OUT / "cluster_calls.csv")["immune_class"].fillna(""))
        assert seen <= emittable


# --------------------------------------------------------------------------------------------
# Heterogeneity, which must be measured rather than assumed away.
# --------------------------------------------------------------------------------------------
def test_the_per_cell_subpanel_argmax_uses_one_marker_per_subpanel():
    """A max over a subpanel's whole set is a max of k, and k differs -- five naming markers for
    myeloid against one for granulocyte.  Measured: correcting this moved the median granulocyte
    agreement from 0.660 to 0.773 and the median T-cell agreement from 0.611 to 0.473, so the
    uncorrected statistic was reporting the shape of the panel, not the mixedness of a cluster."""
    unit = {"CD3e": np.array([0.9, 0.1], dtype=np.float32),
            "CD68": np.array([0.8, 0.2], dtype=np.float32),
            "Iba1": np.array([0.7, 0.95], dtype=np.float32)}
    calls = pd.DataFrame([{"cluster_id": "c0", "marker_T_cell": "CD3e", "marker_Myeloid": "CD68",
                           "marker_B_Plasma": "", "marker_NK_ILC": "", "marker_Granulocyte": ""}])
    got = builder.per_cell_subpanel(calls, ["c0"], np.array([0, 0]), {"unit_magnitude": unit})
    # Iba1 is the louder myeloid channel on cell 1 but it is not the marker that spoke for the
    # myeloid subpanel on this cluster, so it does not enter the comparison.
    assert list(got) == ["T_cell", "Myeloid"]


def test_every_called_cluster_reports_how_concentrated_its_evidence_is():
    if not (OUT / "cluster_heterogeneity.csv").exists():
        pytest.skip("pLN immune pheno has not been built in this checkout")
    het = pd.read_csv(OUT / "cluster_heterogeneity.csv")
    calls = pd.read_csv(OUT / "cluster_calls.csv")
    assert len(het) == len(calls)
    named = het[het["immune_class"].isin(set(builder.SUBPANEL_FILES))]
    assert named["cell_frac_agrees_class"].notna().all()
    assert named["cell_argmax_entropy"].between(0.0, 1.0).all()
    # A mixed cluster must be VISIBLE, not smoothed away.  If every cluster agreed with itself the
    # statistic would be measuring nothing, so its own spread is asserted.
    assert named["cell_frac_agrees_class"].min() < 0.5, \
        "no cluster is mixed, which means the heterogeneity statistic is not measuring anything"


# --------------------------------------------------------------------------------------------
# The artifact on disk, when it is there.
# --------------------------------------------------------------------------------------------
def _manifest():
    path = OUT / "manifest.json"
    if not path.exists():
        pytest.skip("pLN immune pheno has not been built in this checkout")
    return json.loads(path.read_text())


def test_the_manifest_declares_what_it_does_not_read():
    manifest = _manifest()
    for key in ("no_pancreas_input", "no_reference_level", "no_absolute_floor"):
        assert manifest[key].strip()
    assert "clusters" in manifest["clustering_consumed_from"]
    assert manifest["enriched_at"]["auc_above"] == builder.AUC_NULL
    assert "reference_quantile" not in manifest and "positive_at" not in manifest


def test_fourteen_typed_four_unavailable_one_withheld():
    states = pd.Series([s["state"] for s in _manifest()["donors"].values()]).value_counts()
    assert states.get("typed", 0) == 14
    assert states.get("unavailable", 0) == 4
    assert states.get("withheld", 0) == 1


def test_a_donor_with_no_intensity_measurement_is_unavailable_and_never_zero_immune_cells():
    manifest = _manifest()
    for donor in ("6551", "6566", "6593", "6623"):
        state = manifest["donors"][donor]
        assert state["state"] == "unavailable"
        assert "never 0 immune cells" in state["reason"]
        assert state.get("n_immune", "") == ""


def test_every_donor_that_is_not_typed_names_a_reason_on_disk():
    manifest = _manifest()
    for donor, state in manifest["donors"].items():
        if state["state"] == "typed":
            continue
        assert state.get("reason", "").strip(), f"{donor} is withheld with no reason"
        assert (OUT / f"donor_id={donor}/roi=pln_1/UNAVAILABLE.json").exists()


def test_a_withheld_donor_never_sits_beside_a_call_table():
    """Two generations side by side with nothing marking which is which is how this tree goes
    wrong.  A partition carries calls or a refusal, never both."""
    _manifest()
    for part in sorted(OUT.glob("donor_id=*/roi=*")):
        assert not ((part / "data.parquet").exists() and (part / "UNAVAILABLE.json").exists()), \
            f"{part} holds both a call table and a refusal"


def test_the_superseded_generation_is_named_so_this_tree_does_not_stack_two():
    _manifest()
    marker = OUT / "_reference/SUPERSEDED_2026-09-17.json"
    assert marker.exists()
    superseded = json.loads(marker.read_text())
    for name in superseded["superseded_files"]:
        assert "level" in name or "reference_population" in name
    assert "cluster_calls.csv" in superseded["current_outputs"]


def test_object_id_is_preserved_verbatim_and_unique_within_a_donor():
    _manifest()
    import pyarrow.parquet as pq
    tables = sorted(OUT.glob("donor_id=*/roi=*/data.parquet"))
    assert tables, "nothing was typed"
    for path in tables:
        donor = path.parents[1].name.split("=")[1]
        frame = pq.read_table(path, columns=["object_id"]).to_pandas()
        assert frame["object_id"].is_unique
        assert not pd.api.types.is_numeric_dtype(frame["object_id"]), \
            "object_id must stay the QuPath UUID string, never a row number"
        source = next(iter(builder.pln_source(donor).glob("*.parquet")))
        theirs = pq.read_table(source, columns=["object_id"]).to_pandas()["object_id"]
        assert set(frame["object_id"]) <= set(theirs), "object_id was not carried verbatim"


def test_a_cell_with_no_nucleus_is_excluded_from_every_calculation():
    """`qc_analysis_eligible` is the ingest's verdict and `missing_nucleus_geometry` is one of the
    reasons it withholds.  The call table must hold exactly the eligible, clustered cells."""
    _manifest()
    import pyarrow.parquet as pq
    counts = pd.read_csv(OUT / "per_donor_counts.csv", dtype={"donor_id": str})
    for _, row in counts.iterrows():
        donor = str(row["donor_id"])
        source = next(iter(builder.pln_source(donor).glob("*.parquet")))
        eligible = pq.read_table(source, columns=["qc_analysis_eligible"]).to_pandas()
        assert int(eligible["qc_analysis_eligible"].sum()) == row["n_eligible"]
        assert row["n_clustered"] == row["n_immune"] + row["n_non_immune"] + row["n_uncalled"]


def test_an_object_count_is_never_reported_as_a_cell_count():
    """8,961,168 rows in the root, 5,122,269 analysis-eligible.  A prior agent reported the
    10,192,986 QuPath objects as "cells" and inflated the layer by a third."""
    counts = pd.read_csv(OUT / "per_donor_counts.csv") if (OUT / "per_donor_counts.csv").exists() \
        else pytest.skip("not built")
    assert counts["n_eligible"].sum() == 5_122_269
    assert (counts["n_eligible"] < counts["n_rows"]).all()
    assert counts["n_clustered"].sum() == counts["n_eligible"].sum()


def test_every_cell_inherits_exactly_its_own_cluster_s_call():
    """Cells inherit; they are not typed one by one and then grouped."""
    _manifest()
    import pyarrow.parquet as pq
    calls = pd.read_csv(OUT / "cluster_calls.csv", dtype={"donor_id": str})
    for donor in ("115", "6534"):
        path = OUT / f"donor_id={donor}/roi=pln_1/data.parquet"
        cells = pq.read_table(path, columns=["cluster_id", "immune_class", "is_immune"]).to_pandas()
        mine = calls[calls["donor_id"] == donor].set_index("cluster_id")
        per = cells.groupby("cluster_id", observed=True)["immune_class"].nunique()
        assert (per == 1).all(), "a cluster's cells carry more than one class"
        got = cells.groupby("cluster_id", observed=True)["immune_class"].first()
        assert (got == mine["immune_class"].fillna("").reindex(got.index)).all()


def test_uncalled_stays_a_real_category_with_no_default_sink():
    counts = pd.read_csv(OUT / "per_donor_counts.csv") if (OUT / "per_donor_counts.csv").exists() \
        else pytest.skip("not built")
    assert counts["n_uncalled"].sum() > 0, "a node has no default sink; uncalled must survive"
    assert counts["class_Uncalled"].sum() > 0, \
        "no immune cluster failed stage 2, which would mean the requirement never binds"


def test_the_separation_each_donor_actually_got_is_reported():
    path = OUT / "per_donor_separation.csv"
    if not path.exists():
        pytest.skip("pLN immune pheno has not been built in this checkout")
    sep = pd.read_csv(path)
    assert len(sep) == 14
    # A donor is typed only where its clusters separate; the numbers that establish it are on
    # disk, so a reader judges rather than takes the calls on trust.
    assert (sep["clusters_with_an_enriched_immune_canonical"] > 0).all()
    assert (sep["max_immune_standout"] > 0).all()
    assert (sep["panel_offset_spread"] > 0).all()


def test_the_veto_accounting_carries_both_statistics_and_the_channel_s_span():
    """The table on which the defect was visible for a generation.  It must carry the number the
    arbitration READ (the magnitude), the number that admitted the channel (the rank), and the
    deciding channel's whole span on that section -- all three, so that an exclusion can be
    judged from the row instead of taken on trust."""
    calls = pd.DataFrame([
        {"cluster_id": "c0", "n_cells": 100, "is_immune": False,
         "top_structural_marker": "Caveolin", "structural_magnitude": 8.0,
         "immune_magnitude": 1.1, "structural_standout": 0.3, "immune_standout": 0.1},
        {"cluster_id": "c1", "n_cells": 40, "is_immune": True,
         "top_structural_marker": "EpCAM", "structural_magnitude": 1.3,
         "immune_magnitude": 2.4, "structural_standout": 0.4, "immune_standout": 0.2},
    ])
    rows = builder.structural_veto_rows("T", calls, {"Caveolin": 20.3, "EpCAM": 1.3})
    assert len(rows) == 1, "only an exclusion is accounted for here"
    row = rows[0]
    assert row["structural_marker"] == "Caveolin" and row["n_cells_excluded"] == 100
    assert row["median_structural_magnitude"] == 8.0
    assert row["median_immune_magnitude_of_those_clusters"] == 1.1
    assert row["median_structural_standout"] == 0.3
    assert row["marker_dynamic_range_on_this_section"] == 20.3


def test_every_exclusion_a_structural_marker_caused_is_accounted_for_with_its_span():
    """Every excluded cluster appears exactly once, and the cell counts reconcile with the calls.

    Before the 2026-09-18 correction this table said, on the cohort it was built from, that
    429,278 of 1,669,194 excluded cells were decided by an epithelial channel in a tissue with no
    epithelium -- EpCAM at a median span of 1.60x against 20.3x for Caveolin."""
    path = OUT / "structural_veto_accounting.csv"
    if not path.exists():
        pytest.skip("pLN immune pheno has not been built in this checkout")
    veto = pd.read_csv(path, dtype={"donor_id": str})
    calls = pd.read_csv(OUT / "cluster_calls.csv", dtype={"donor_id": str})
    excluded = calls[calls["is_immune"] == False]                                # noqa: E712
    assert veto["n_clusters_excluded"].sum() == len(excluded)
    assert veto["n_cells_excluded"].sum() == excluded["n_cells"].sum()
    assert veto["marker_dynamic_range_on_this_section"].notna().all()


# --------------------------------------------------------------------------------------------
# Coherence: internal contradictions only.  No expectation is ever a check.
# --------------------------------------------------------------------------------------------
def test_cd4_and_cd8_never_claim_more_cells_than_exist():
    """The UNION is the arithmetic impossibility: CD4-positive or CD8-positive cannot exceed 100 %
    of the T compartment, whatever a node contains.

    The SUM may exceed 100 % when many cells are positive for both, and a loose definition of
    positive produces that on its own -- so the sum is reported across a grid of section-internal
    definitions with the double-positive share beside it, and the q at which it happens is
    recorded rather than the grid being trimmed until it does not.  No call reads q: the calls are
    made on cluster enrichment, so this is a diagnostic of the diagnostic, not of the calls.
    """
    path = OUT / "_reference/coherence.csv"
    if not path.exists():
        pytest.skip("pLN immune pheno has not been built in this checkout")
    table = pd.read_csv(path, dtype={"donor_id": str})
    assert len(table) == 14
    for q in builder.COHERENCE_Q:
        tag = f"{q:g}".replace("0.", "p")
        union = table[f"cd4_or_cd8_pct_of_T_at_{tag}"]
        assert (union <= 100.0 + 1e-9).all(), \
            table.loc[union > 100.0, ["donor_id", f"cd4_or_cd8_pct_of_T_at_{tag}"]].to_string()
        assert np.allclose(table[f"cd4_plus_cd8_pct_of_T_at_{tag}"]
                           - table[f"cd4_and_cd8_pct_of_T_at_{tag}"], union, atol=1e-6)
    assert (table["union_over_100_at"].fillna("") == "").all()


def test_where_the_sum_passes_100_the_q_that_did_it_is_named():
    path = OUT / "_reference/coherence.csv"
    if not path.exists():
        pytest.skip("pLN immune pheno has not been built in this checkout")
    table = pd.read_csv(path, dtype={"donor_id": str})
    flagged = table[table["contradiction_cd4_cd8"]]
    assert (flagged["contradiction_cd4_cd8_at"].fillna("").str.len() > 0).all()
    # It must never be the whole grid: if the sum passed 100 % at the tightest definition too,
    # the calls themselves would be double-counting rather than the definition being loose.
    tightest = f"{min(builder.COHERENCE_Q):g}".replace("0.", "p")
    assert (table[f"cd4_plus_cd8_pct_of_T_at_{tightest}"] <= 100.0).all()


# --------------------------------------------------------------------------------------------
# Stage 3: the T_NK compartment, split by NEGATIVE GATING on CD3e.
#
# The withdrawn generation of this lane wrote class_NK_ILC = 0 on all fourteen donors and argued
# that NK could not be named because no single channel on this panel is NK-specific.  None has to
# be.  NK is CD3e-NEGATIVE with CD8 or CD56 and T cells are CD3e-POSITIVE, so the panel resolves
# the two by COMBINATION -- absence carrying as much as presence -- and these tests hold the split
# to that and to nothing else.  They are written against a SYNTHETIC section whose answer is known
# by construction, so nothing here reproduces a result from the cohort.
# --------------------------------------------------------------------------------------------
FILLERS = [f"Fill_{k}" for k in range(10)]
STAGE3_MARKERS = ["CD3e", "CD4", "CD8", "CD56", "FOXP3"] + FILLERS


def _synthetic_section(seed: int = 11, n_nk: int = 150, scale: dict | None = None):
    """A section of two clusters whose composition is known before anything is measured.

    Cluster c0 is a T zone: 900 CD4 T, 450 CD8 T, `n_nk` NK.  Cluster c1 is flat background.
    Every cell carries its own BRIGHTNESS, multiplying every channel it has, which is the
    confound the per-cell offset exists for: without it the dim cells of a cluster read low on
    CD3e, CD8 and CD56 alike and no gate can tell an NK cell from a dim T cell.
    """
    rng = np.random.default_rng(seed)
    scale = scale or {}
    n_t4, n_t8, n_bg = 900, 450, 1500
    n = n_t4 + n_t8 + n_nk + n_bg
    bright = rng.uniform(0.5, 2.0, n)
    v = {m: np.ones(n) for m in STAGE3_MARKERS}
    t4 = slice(0, n_t4)
    t8 = slice(n_t4, n_t4 + n_t8)
    nk = slice(n_t4 + n_t8, n_t4 + n_t8 + n_nk)
    for s in (t4, t8):
        v["CD3e"][s] = 10.0
    v["CD4"][t4], v["CD8"][t8] = 6.0, 6.0
    v["CD8"][nk], v["CD56"][nk] = 4.0, 9.0                 # CD3e stays at background for NK
    frame = pd.DataFrame({m: v[m] * bright * scale.get(m, 1.0) * rng.uniform(0.95, 1.05, n)
                          for m in STAGE3_MARKERS})
    codes = np.zeros(n, dtype=np.int64)
    codes[n_t4 + n_t8 + n_nk:] = 1
    truth = np.array(["T_CD4"] * n_t4 + ["T_CD8"] * n_t8 + ["NK_ILC"] * n_nk + [""] * n_bg,
                     dtype=object)
    return frame, codes, truth


def _resolve(frame, codes, compartment=("c0",)):
    ids = ["c0", "c1"]
    st = builder.standout(frame, codes, 2, list(frame.columns), offset_over=list(frame.columns))
    calls = pd.DataFrame({"cluster_id": ids,
                          "immune_class": [builder.T_NK_CLASS if c in compartment else "Myeloid"
                                           for c in ids]})
    return builder.resolve_t_nk("D", calls, ids, codes, st, list(frame.columns))


def test_nk_is_named_by_the_absence_of_cd3e_together_with_cd8_or_cd56():
    """The whole point.  These cells carry no NK-specific channel -- CD56 is shared with the
    Neural lineage and CD8 with T cells -- and they are named anyway, from the pattern."""
    frame, codes, truth = _synthetic_section()
    res, cells = _resolve(frame, codes)
    called = cells["t_nk_class"]
    assert res.loc[0, "split_accepted"], res.loc[0, "resolution_reason"]
    assert res.loc[0, "resolution_mode"] == "sub_cluster_gate_on_CD3e"
    named = (called == "NK_ILC")
    assert named.sum() > 0, "NK must be named, not reported unnameable"
    # Every NK cell of this section is found, and the cells wrongly swept in are the boundary
    # ones a threshold-free per-cell rule is bound to take: measured 3 of 153 here.  The test
    # asserts the recall exactly and the purity with the slack that boundary leaves.
    assert (truth[named] == "NK_ILC").sum() == (truth == "NK_ILC").sum()
    assert (truth[named] == "NK_ILC").mean() >= 0.95


def test_cd3e_negative_is_the_load_bearing_term_and_nothing_is_named_without_it():
    """Take the absence away -- give the NK cells the same CD3e as the T cells -- and the gate
    must stop naming them.  A gate that fires on CD8/CD56 alone would not be a negative gate."""
    frame, codes, truth = _synthetic_section()
    frame = frame.copy()
    # Raise their CD3e to the T level WITHOUT touching their brightness -- they keep the CD56
    # that made them findable and lose only the absence.  (Assigning a constant instead would
    # dim the bright cells' CD3e relative to their own panel and put the absence back.)
    frame.loc[truth == "NK_ILC", "CD3e"] *= 10.0
    res, cells = _resolve(frame, codes)
    assert int((cells["t_nk_class"] == "NK_ILC").sum()) == 0
    assert res.loc[0, "resolution_mode"] != "sub_cluster_gate_on_CD3e"


def test_the_gate_reads_no_threshold_and_survives_a_change_of_units():
    """Multiplying a channel by 1000 must not move one call.  The gate reads within-section ranks
    and a per-cell median of them, so no intensity, level, quantile or scale factor can enter."""
    frame, codes, _ = _synthetic_section()
    base = _resolve(frame, codes)[1]["t_nk_class"]
    for marker, factor in (("CD56", 1000.0), ("CD3e", 0.001), ("Fill_3", 1e6)):
        scaled = frame.copy()
        scaled[marker] = scaled[marker] * factor
        assert (_resolve(scaled, codes)[1]["t_nk_class"] == base).all(), marker


def test_the_per_cell_offset_is_a_median_over_every_channel_and_cancels_inside_one_cell():
    """`cell_offset` is the cell's own median place over the WHOLE panel, so it cannot be steered
    by the markers this lane cares about -- and because it is one number per cell it cancels in
    any comparison between two markers OF THAT CELL, which is why the CD4/CD8/FOXP3 argmax is
    identical with and without it.  Only comparisons BETWEEN cells move."""
    frame, codes, _ = _synthetic_section()
    markers = list(frame.columns)
    st = builder.standout(frame, codes, 2, markers, offset_over=markers)
    block = np.vstack([st["unit_rank"][m] for m in markers])
    assert np.allclose(st["cell_offset"], np.nanmedian(block, axis=0), atol=1e-6)
    raw = np.vstack([st["unit_rank"][m] for m in ("CD4", "CD8", "FOXP3")])
    corrected = raw - st["cell_offset"]
    assert (raw.argmax(axis=0) == corrected.argmax(axis=0)).all()


def test_a_cell_with_no_cd3e_measurement_is_unresolved_and_never_swept_onto_either_side():
    """An unavailable measurement is not a negative one.  A cell that carries no CD3e cannot be
    gated on CD3e, so it keeps the coarse label instead of being handed to the T side."""
    frame, codes, truth = _synthetic_section()
    frame = frame.copy()
    blind = np.zeros(len(frame), dtype=bool)
    blind[np.flatnonzero(truth == "T_CD4")[:50]] = True
    frame.loc[blind, "CD3e"] = np.nan
    res, cells = _resolve(frame, codes)
    called = cells["t_nk_class"]
    assert (called[blind] == builder.T_NK_RESIDUAL).all()
    assert int(res.loc[0, "n_cd3e_unmeasured"]) == 50


def test_the_subsets_are_mutually_exclusive_and_sum_to_the_compartment():
    """Internal coherence, which is legitimate because it is arithmetic: a cell of the
    compartment leaves stage 3 with exactly one label, so the four subsets plus the residual are
    the compartment and CD4 % + CD8 % among the CD3e-positive cells cannot pass 100 %."""
    frame, codes, _ = _synthetic_section()
    res, cells = _resolve(frame, codes)
    called = cells["t_nk_class"]
    n_comp = int((codes == 0).sum())
    assert int(sum((called == k).sum() for k in builder.T_NK_CLASSES)) == n_comp
    assert (called[codes == 1] == "").all(), "stage 3 must not label a cell outside its compartment"
    summary = builder.t_nk_summary("D", res, n_comp)
    assert summary["subsets_sum_to_the_compartment"]
    assert summary["pct_CD4_plus_CD8_of_CD3e_positive"] <= 100.0 + 1e-9


def test_a_cluster_that_does_not_split_keeps_the_coarse_label_and_says_why():
    """A refusal is a real answer and it has to carry its numbers.  A pure T cluster with no NK
    in it must come back with no NK named and a reason that states what was measured."""
    frame, codes, _ = _synthetic_section(n_nk=0)
    res, cells = _resolve(frame, codes)
    assert int((cells["t_nk_class"] == "NK_ILC").sum()) == 0
    assert not res.loc[0, "split_accepted"]
    reason = res.loc[0, "resolution_reason"]
    assert "AUC" in reason and "alpha" in reason
    for column in ("n_cd3e_negative", "nk_side_auc_cd8", "nk_side_auc_cd56", "nk_side_p",
                   "within_cluster_corr_cd3e_cd56", "eta2_at_the_zero_cut"):
        assert column in res.columns and pd.notna(res.loc[0, column]), column


def test_nk_ilc_is_never_written_as_zero_when_the_compartment_did_not_resolve():
    """THE ACCEPTANCE TEST.  0 is a measurement only when the gate actually ran.  When every
    cluster of a non-empty compartment came back unresolved -- or when the gate could not run at
    all -- the count is written EMPTY with the reason, because an unavailable measurement is not
    zero and a reader must not be able to read one as the other."""
    empty = pd.DataFrame(columns=["donor_id", "cluster_id", "n_cells", "resolution_mode",
                                  "gate_could_run"])
    out = builder.t_nk_summary("D", empty, compartment_cells=0)
    assert out["NK_ILC"] == "" and "unavailable" in out["nk_ilc_state"]

    stuck = pd.DataFrame([{"donor_id": "D", "cluster_id": "c0", "n_cells": 5000,
                           "resolution_mode": "did_not_resolve", "gate_could_run": True,
                           **{f"n_{k}": 0 for k in builder.T_NK_CLASSES}}])
    stuck.loc[0, f"n_{builder.T_NK_RESIDUAL}"] = 5000
    out = builder.t_nk_summary("D", stuck, compartment_cells=5000)
    assert out["NK_ILC"] == "", "an unresolved compartment must not report 0 NK cells"
    assert out["NK_ILC"] != 0 and out["nk_ilc_state"].startswith("unavailable")

    blind = stuck.copy()
    blind.loc[0, "gate_could_run"] = False
    out = builder.t_nk_summary("D", blind, compartment_cells=5000)
    assert out["NK_ILC"] == "" and "NOT zero NK cells" in out["nk_ilc_state"]

    ran = pd.DataFrame([{"donor_id": "D", "cluster_id": "c0", "n_cells": 5000,
                         "resolution_mode": "whole_cluster_CD3e_enriched", "gate_could_run": True,
                         **{f"n_{k}": 0 for k in builder.T_NK_CLASSES}}])
    ran.loc[0, "n_T_CD4"] = 5000
    out = builder.t_nk_summary("D", ran, compartment_cells=5000)
    assert out["NK_ILC"] == 0 and out["nk_ilc_state"] == "measured"


def test_the_gate_is_corrected_over_its_own_family_of_tests_and_not_stage_ones():
    """Stage 1 runs ~1,000 cluster x marker tests on a section; stage 3 runs one per compartment
    cluster per NK term.  Correcting the gate over stage 1's family would demand an alpha two to
    three orders of magnitude smaller than the family being run, and the count is fixed before
    any test so no result can change it."""
    frame, codes, _ = _synthetic_section()
    res, _ = _resolve(frame, codes)
    assert int(res.loc[0, "n_gate_tests_on_this_section"]) == 1 * len(builder.NK_POSITIVE)
    assert res.loc[0, "gate_bonferroni_alpha"] == pytest.approx(
        builder.ALPHA / int(res.loc[0, "n_gate_tests_on_this_section"]))
    assert res.loc[0, "gate_bonferroni_alpha"] > res.loc[0, "stage1_section_bonferroni_alpha"]


# --------------------------------------------------------------------------------------------
# D1.  No immune canonical may be listed as a non-immune structural marker.
# --------------------------------------------------------------------------------------------
def test_no_immune_canonical_is_listed_as_a_non_immune_structural_marker():
    """CD66 sat in `cluster_pln.PHENO_NON_IMMUNE` while `lineages` names it a Granulocyte
    canonical, so a granulocyte cluster was scored as structural by its own marker.  `lineages`
    is the authority on what an antibody means, and both lists are checked against it here rather
    than against a frozen copy, so the marker cannot come back under either name."""
    canonical = {m for lin in builder.IMMUNE_LINEAGES for m in LIN.LINEAGE_MARKERS[lin]}
    assert not (canonical & set(builder.STRUCTURAL)), "an immune canonical cannot exclude a cluster"

    source = (REPO / "core/phenocycler/scripts/pln/cluster_pln.py").read_text()
    tree = ast.parse(source)
    listed = None
    for node in ast.walk(tree):
        if (isinstance(node, ast.Assign) and len(node.targets) == 1
                and getattr(node.targets[0], "id", "") == "PHENO_NON_IMMUNE"):
            listed = {e.value for e in node.value.elts}
    assert listed is not None, "cluster_pln must still declare PHENO_NON_IMMUNE"
    assert not (canonical & listed), sorted(canonical & listed)
    assert "CD66" in LIN.LINEAGE_MARKERS["Granulocyte"] and "CD66" not in listed


def test_cd66_still_admits_a_cluster_as_granulocyte_evidence():
    """The other half of D1: removing CD66 from a structural list must not remove it from the
    evidence.  It is a stage-1 primary, and `lineages` records that an earlier lane dropping it
    "folded granulocytes into the myeloid label"."""
    assert "CD66" in builder.PRIMARY_IMMUNE["Granulocyte"]
    calls = _call([{"CD66": 0.4, "SMA": 0.1}], MARKERS)
    assert bool(calls.loc[0, "is_immune"]) is True
    assert calls.loc[0, "top_immune_marker"] == "CD66"


# --------------------------------------------------------------------------------------------
# D4.  A partial run never overwrites a cohort-wide artifact.
# --------------------------------------------------------------------------------------------
def test_a_partial_run_writes_its_own_scope_and_leaves_the_cohort_tables_alone():
    """`--donors X` and `--stage standout` used to write straight over `cluster_calls.csv` and
    `per_donor_separation.csv`, taking 293 cluster rows to 19 and 14 donor rows to 1.  A run that
    does not cover every clustered donor AND every stage is now scoped to `_partial/<scope>/`."""
    tree = ast.parse(SCRIPT.read_text())
    names = {n.id for n in ast.walk(tree) if isinstance(n, ast.Name)}
    assert {"COHORT", "partial"} <= names
    source = SCRIPT.read_text()
    for table in ("cluster_calls.csv", "per_donor_separation.csv", "per_donor_counts.csv",
                  "per_donor_status.csv", "cluster_heterogeneity.csv", "manifest.json",
                  "t_nk_resolution.csv", "per_donor_t_nk.csv"):
        assert f'OUT / "{table}"' not in source, f"{table} is still written to the cohort root"
        assert f'COHORT / "{table}"' in source, table


# --------------------------------------------------------------------------------------------
# Stage 3 against the artifact on disk.
# --------------------------------------------------------------------------------------------
def _t_nk_tables():
    counts = OUT / "per_donor_counts.csv"
    donors = OUT / "per_donor_t_nk.csv"
    clusters = OUT / "t_nk_resolution.csv"
    if not (counts.exists() and donors.exists() and clusters.exists()):
        pytest.skip("the pLN stage-3 artifact is not on disk")
    return pd.read_csv(counts), pd.read_csv(donors), pd.read_csv(clusters)


def test_on_disk_the_subsets_sum_to_the_t_nk_compartment_on_every_donor():
    counts, donors, _ = _t_nk_tables()
    subsets = [f"tnk_{k}" for k in builder.T_NK_CLASSES]
    total = counts[subsets].fillna(0).sum(axis=1).astype(int)
    assert (total == counts["class_T_cell"]).all(), counts.loc[total != counts["class_T_cell"]]
    assert bool(donors["subsets_sum_to_the_compartment"].all())


def test_on_disk_cd4_and_cd8_never_claim_more_than_the_cd3e_positive_cells():
    _, donors, _ = _t_nk_tables()
    over = donors.loc[donors["pct_CD4_plus_CD8_of_CD3e_positive"] > 100.0 + 1e-9]
    assert over.empty, over[["donor_id", "pct_CD4_plus_CD8_of_CD3e_positive"]]


def test_on_disk_every_cluster_of_the_compartment_says_which_way_it_went_and_why():
    _, _, clusters = _t_nk_tables()
    assert len(clusters) > 0
    assert clusters["resolution_mode"].isin(
        ["sub_cluster_gate_on_CD3e", "whole_cluster_CD3e_at_or_below_its_panel_level",
         "whole_cluster_CD3e_enriched", "did_not_resolve"]).all()
    assert clusters["resolution_reason"].str.len().gt(40).all()
    # every refusal carries the numbers behind it, so it can be read rather than taken on trust
    refused = clusters.loc[~clusters["split_accepted"].astype(bool)]
    for column in ("n_cd3e_negative", "nk_side_auc_negative_over_positive", "nk_side_p",
                   "gate_bonferroni_alpha", "n_cd3e_negative_and_cd8_or_cd56_positive"):
        assert refused[column].notna().all(), column


def test_on_disk_a_zero_nk_count_only_stands_where_the_gate_actually_ran():
    counts, donors, _ = _t_nk_tables()
    zeros = donors.loc[donors["NK_ILC"].astype(str).isin(["0", "0.0"])]
    assert (zeros["nk_ilc_state"] == "measured").all()
    assert (zeros["n_cells_where_the_CD3e_gate_could_run"] > 0).all()
    unavailable = donors.loc[donors["nk_ilc_state"] != "measured"]
    assert unavailable["NK_ILC"].isna().all() or (unavailable["NK_ILC"].astype(str) == "").all()
