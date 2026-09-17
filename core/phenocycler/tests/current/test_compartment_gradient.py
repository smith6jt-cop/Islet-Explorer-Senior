"""Contract tests for the compartment gradient.

Run from the repository root with .venv/bin/python -m pytest <this file>.

Fixtures are synthetic measurements, deliberately chosen so that no assertion
can pass by reproducing a biological expectation.  The immune fractions below
are not plausible, on purpose: what is checked is that the compartment cut, the
denominators, the stratifications and the tissue boundary behave as written, not
that a number came out the size someone hoped.
"""

from __future__ import annotations

import importlib.util
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

REPO = Path(__file__).resolve().parents[4]
BUILDER = REPO / "scripts/senior/build_compartment_gradient.py"
spec = importlib.util.spec_from_file_location("build_compartment_gradient", BUILDER)
builder = importlib.util.module_from_spec(spec)
sys.modules["build_compartment_gradient"] = builder
spec.loader.exec_module(builder)

OUT = REPO / "data" / "development" / "integration" / "compartments"


# --------------------------------------------------------------------------
# the compartment cut


def test_compartment_comes_from_the_signed_distance_and_the_expansion_is_20um():
    dist = pd.Series([-40.0, -0.5, 0.0, 0.001, 19.9, 20.0, 20.0001, 500.0])
    got = builder.compartment_from_signed_distance(dist).tolist()
    assert got == ["core", "core", "core", "peri", "peri", "peri",
                   "exocrine", "exocrine"]
    assert builder.EXPANSION_UM == 20.0


def test_a_cell_without_a_signed_distance_has_no_compartment_rather_than_a_default():
    dist = pd.Series([np.nan, 5.0])
    got = builder.compartment_from_signed_distance(dist)
    assert pd.isna(got.iloc[0]), "a missing distance is unavailable, never exocrine"
    assert got.iloc[1] == "peri"


# --------------------------------------------------------------------------
# the tissue boundary


def test_a_distance_is_refused_outside_the_pancreas():
    builder.refuse_cross_tissue("pancreas", "signed_dist_islet_um")  # does not raise
    with pytest.raises(builder.CrossTissueDistance) as err:
        builder.refuse_cross_tissue("pancreatic_lymph_node", "signed_dist_islet_um")
    assert "no islet" in str(err.value)


def test_the_pending_attachment_classification_is_named_and_not_inferred():
    # The maintainer's pLN-attachment dataset has not landed. The builder must
    # name the column it will honour and must not synthesise it.
    assert builder.LYMPHOID_ATTACHMENT_COLUMN == "attached_to_pancreas_section"
    source = BUILDER.read_text()
    assert "is not inferred" in source or "not inferred here" in source


# --------------------------------------------------------------------------
# denominators and internal coherence


def _fixture_frame():
    return pd.DataFrame({
        "is_immune": [True, True, True, False, False],
        "common": ["T_NK", "Myeloid", None, None, None],
    })


def test_subtypes_sum_to_the_immune_count_and_never_exceed_the_denominator():
    rows = builder._composition_rows("D", "pancreas", "phenocycler", "core",
                                     _fixture_frame(), "all_islets")
    total = next(r for r in rows if r["lineage"] == "Immune_total")
    lineages = [r for r in rows if r["lineage"] != "Immune_total"]
    assert total["n_cells"] == 5
    assert total["n_lineage"] == 3
    # One immune cell carries no resolved lineage. It is still immune: the
    # lineage rows sum to at most the immune count, never past it.
    assert sum(r["n_lineage"] for r in lineages) == 2 <= total["n_immune"]
    assert all(r["n_cells"] == total["n_cells"] for r in rows)


def test_an_empty_compartment_is_unavailable_with_a_reason_never_zero():
    empty = _fixture_frame().iloc[0:0]
    rows = builder._composition_rows("D", "pancreas", "phenocycler", "core",
                                     empty, "all_islets")
    assert all(r["state"] == "unavailable" for r in rows)
    assert all(r["reason"] for r in rows)
    assert all(pd.isna(r["fraction_of_compartment"]) for r in rows)


# --------------------------------------------------------------------------
# batch


def test_eta_squared_is_reported_with_a_permutation_null():
    rng = np.random.default_rng(7)
    values = rng.normal(size=24)
    groups = np.array(["a", "b", "c"] * 8, dtype=object)
    stat = builder.eta_with_null(values, groups, n_perm=200, seed=1)
    assert 0.0 <= stat["eta_squared"] <= 1.0
    assert 0.0 <= stat["null_median"] <= 1.0
    assert 0.0 < stat["p_permutation"] <= 1.0
    assert stat["n_perm"] == 200


def test_eta_squared_is_undefined_rather_than_zero_when_it_cannot_be_computed():
    assert np.isnan(builder.eta_squared(np.array([1.0, 1.0, 1.0]),
                                        np.array(["a", "a", "b"], dtype=object)))
    assert np.isnan(builder.eta_squared(np.array([1.0, 2.0]),
                                        np.array(["a", "a"], dtype=object)))


# --------------------------------------------------------------------------
# exclusions


def test_the_excluded_donor_list_is_derived_not_restated():
    source = BUILDER.read_text()
    assert "from phenocycler.cohort import DONOR_EXCLUSIONS" in source
    assert "6457" not in source, "the exclusion is central; do not restate it here"
    assert "6579" not in source


def test_the_centrally_excluded_donor_is_not_built():
    donors = builder.eligible_donors()
    assert "6457" not in donors
    if (builder.CELLS / "donor_id=6579").exists():
        assert "6579" in donors, "6579's exclusion was lifted 2026-09-17"


# --------------------------------------------------------------------------
# the built artifact, when it is on disk


requires_artifact = pytest.mark.skipif(
    not (OUT / "compartment_composition.csv").exists(),
    reason="compartment gradient not built in this checkout",
)


@requires_artifact
def test_built_pln_rows_carry_no_pancreas_geometry():
    comp = pd.read_csv(OUT / "compartment_composition.csv", dtype={"donor_id": str})
    pln = comp[comp["tissue"] != "pancreas"]
    assert len(pln)
    assert set(pln["compartment"]) == {"pln"}
    assert set(pln["stratum"]) == {"not_applicable_no_islets"}
    assert not {"distance_bin_um", "signed_dist_islet_um", "frame"} & set(comp.columns)


@requires_artifact
def test_built_counts_are_internally_coherent():
    comp = pd.read_csv(OUT / "compartment_composition.csv", dtype={"donor_id": str})
    m = comp[comp["state"] == "measured"]
    assert (m["n_lineage"] <= m["n_cells"]).all()
    assert (m["n_immune"] <= m["n_cells"]).all()
    key = ["donor_id", "tissue", "modality", "compartment", "stratum"]
    lineages = m[m["lineage"] != "Immune_total"].groupby(key)["n_lineage"].sum()
    totals = m[m["lineage"] == "Immune_total"].set_index(key)["n_immune"]
    joined = pd.concat([lineages, totals], axis=1).dropna()
    assert (joined["n_lineage"] <= joined["n_immune"]).all()


@requires_artifact
def test_built_manifest_records_the_domain_and_the_reasons():
    manifest = json.loads((OUT / "manifest.json").read_text())
    assert manifest["source_domain_root"].endswith("typing/app")
    assert manifest["parent_column"].startswith("never read")
    assert "donor_id" in manifest["cross_tissue"]
    assert manifest["pending_maintainer_classification"]["column"] == \
        builder.LYMPHOID_ATTACHMENT_COLUMN
    # Every pLN donor that produced nothing says why, by name.
    for note in manifest["pln_phenocycler"]:
        if note.get("state") != "measured":
            assert note.get("reason"), note


@requires_artifact
def test_both_stratifications_are_present_and_applied():
    strata = pd.read_csv(OUT / "islet_size_summary.csv", dtype={"donor_id": str})
    assert (strata["n_islets_ge_10"] <= strata["n_islets"]).all()
    # The fragmentation flag is measured, and the flagged donors are still built.
    flagged = set(strata.loc[strata["flag"].fillna("") != "", "donor_id"])
    comp = pd.read_csv(OUT / "compartment_composition.csv", dtype={"donor_id": str})
    assert flagged <= set(comp["donor_id"]), "a flagged donor is reported, never dropped"
    st = pd.read_csv(OUT / "islet_size_strata.csv", dtype={"donor_id": str})
    assert set(st["stratum"]) == {"islets_ge_10_core_cells"}
    assert set(st["compartment"]) == {"core", "peri"}


@requires_artifact
def test_batch_eta_squared_is_exported_on_the_finer_variables_with_a_null():
    eta = pd.read_csv(OUT / "batch_eta_squared.csv")
    assert {"stain_group", "acq_run"} <= set(eta["grouping"])
    # Never the Fusion version, and never the coarse batch, which is aliased
    # with disease status.
    assert not {"software", "coarse_batch"} & set(eta["grouping"])
    live = eta[eta["eta_squared"].notna()]
    assert len(live)
    assert live["null_median"].notna().all()
    assert live["p_permutation"].notna().all()
