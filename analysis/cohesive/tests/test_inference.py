import numpy as np
import pandas as pd

from analysis.cohesive.inference import bh, contrast, h6_association


def donors(values, groups, status=None):
    return pd.DataFrame({"donor_id": [str(i) for i in range(len(values))], "donor_effect": values,
                         "disease_status": groups, "status": status or ["available"] * len(values)})


def test_exact_enumeration_matches_hand_count():
    table = donors([3.0, 4.0, 0.0, 1.0, 2.0], ["T1D", "T1D", "ND", "ND", "ND"])
    result = contrast(table, 99999, 200, 1)
    assert result["method"] == "exact enumeration of disease labels"
    assert result["effect"] == 3.5 - 1.0
    # Of the C(5,2)=10 labellings the observed one and its mirror (values 0 and 1) reach |2.5|.
    assert result["p"] == 0.2


def test_unavailable_donors_and_missing_groups_never_become_zero():
    table = donors([1.0, np.nan, 2.0], ["T1D", "ND", "ND"], ["available", "unavailable", "available"])
    assert contrast(table, 999, 50, 1)["n_ND"] == 1
    only_t1d = donors([1.0, 2.0], ["T1D", "T1D"])
    result = contrast(only_t1d, 999, 50, 1)
    assert result["status"] == "unavailable" and np.isnan(result["effect"])


def test_stratified_permutation_keeps_labels_within_strata():
    values = [10.0, 11.0, 0.0, 1.0, 10.5, 0.5]
    groups = ["T1D", "T1D", "ND", "ND", "T1D", "ND"]
    table = donors(values, groups)
    strata = pd.Series(["a", "a", "b", "b", "a", "b"], index=table.index)
    # Disease is perfectly aliased with the stratum, so no within-stratum relabelling changes the effect.
    result = contrast(table, 2000, 50, 1, strata=strata)
    assert result["p"] == 1.0


def test_bh_matches_reference_values():
    q = bh(np.array([0.01, 0.04, 0.03, 0.2, 1.0, 1.0]))
    np.testing.assert_allclose(q, [0.06, 0.08, 0.08, 0.3, 1.0, 1.0])


def test_h6_needs_more_pairs_than_rank_and_reports_reason():
    paired = pd.DataFrame({"donor_id": ["1", "2"], "disease_status": ["ND", "T1D"], "paired_status": "available",
                           "H3": [0.1, 0.2], "architecture_excess": [0.3, 0.1], "mean_memory_state": [1.0, 2.0],
                           "mean_cytotoxic_state": [0.5, 0.4]})
    result, pairs = h6_association(paired, 99, 10, 1)
    assert result["status"] == "unavailable" and pairs is None and result["reason"]
