"""Synthetic tests for the Xenium typing lane: record shape, the cluster-award rule, the release label space."""
from __future__ import annotations

import pathlib
import sys

import numpy as np
import pandas as pd
import pytest

HERE = pathlib.Path(__file__).resolve().parents[1]
sys.path.insert(0, str(HERE))
import gate  # noqa: E402
import pins  # noqa: E402


def test_pins_serial_and_run_maps_cover_the_cohort():
    donors = pins.donors()
    assert len(donors) == 20
    for d in ("6436", "6442"):
        assert pins.serial(d) in ("0059907", "0059947") and pins.xenium_run(d).startswith("MM04")
    assert len(set(pins.xenium_run(d) for d in donors)) == 11


def test_validate_calls_refuses_a_missing_reason_and_an_unknown_label():
    df = pd.DataFrame({"cell_id": ["a", "b"], "lineage_native": ["Endocrine", gate.UNAVAILABLE], "celltype_fine": ["Beta", ""],
                       "confidence": [0.9, np.nan], "reason": ["", ""], "cycling": [False, False]})
    with pytest.raises(ValueError):
        gate.validate_calls(df)
    df.loc[1, "reason"] = "qc"
    gate.validate_calls(df)
    bad = df.copy(); bad.loc[0, "lineage_native"] = "Neural_X"
    with pytest.raises(ValueError):
        gate.validate_calls(bad)


def test_release_label_space_maps_every_lineage_to_the_common_vocabulary():
    import release
    from phenocycler.integration.vocab import COMMON_LINEAGES
    space = release.label_space("joint-run-aware")
    mapped = {r["label"]: r.get("common_lineage") for r in space["labels"] if "common_lineage" in r}
    assert set(mapped) == set(pins.LINEAGES)
    assert set(mapped.values()) <= set(COMMON_LINEAGES) - {"Unassigned"}
    resolutions = {r["label"]: r.get("common_resolution") for r in space["labels"] if "common_resolution" in r}
    assert resolutions == {gate.AMBIGUOUS: "ambiguous", gate.UNAVAILABLE: "unavailable"}


def test_baseline_cluster_award_is_exactly_one_lineage(monkeypatch):
    import baseline_analyzer as ba
    long = pd.DataFrame({
        "group": ["c1"] * 3 + ["c2"] * 3 + ["c3"] * 3,
        "lineage": ["Endocrine", "Exocrine", "Myeloid"] * 3,
        "clears": [True, False, False, True, True, False, False, False, False],
        "judgeable": [True] * 9,
    })
    monkeypatch.setattr(gate, "nullgate_table", lambda *a, **k: (long, None))
    lab, rsn, _ = ba.annotate_clusters(None, None, np.array(["c1", "c2", "c3"]), {})
    assert lab["c1"] == "Endocrine" and rsn["c1"] == ""
    assert lab["c2"] == gate.AMBIGUOUS and rsn["c2"] == "multi_lineage"
    assert lab["c3"] == gate.UNAVAILABLE and rsn["c3"] == "no_lineage_beats_null"
