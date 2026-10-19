"""Immune definitions veto epithelium only with markers every panel acquired (ruled 2026-10-18).

The 70-channel same-slide panel (0112, 6347, 6533, 6580) never acquired Ker8_18. While every immune
definition required Ker8_18 measured-negative, no definition could fire on those donors at all.
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd

REPO = Path(__file__).resolve().parents[4]
sys.path.insert(0, str(REPO / "core/phenocycler/scripts/typing"))
sys.path.insert(0, str(REPO / "core/phenocycler/scripts/celltype"))

import immune_subtype  # noqa: E402
import refine  # noqa: E402
import subpanel  # noqa: E402

EPITHELIAL = ("E_cadherin", "EpCAM", "Pan_Cytokeratin")


def t_cell(**values) -> pd.DataFrame:
    base = {"CD3e": 0.6, "E_cadherin": 0.0, "EpCAM": 0.0, "Pan_Cytokeratin": 0.0, "Ker8_18": np.nan}
    return pd.DataFrame([{**base, **values}])


def test_no_immune_definition_vetoes_with_an_unshared_marker():
    for rules in (subpanel.immune_rules(), immune_subtype.immune_rules()[1], refine.immune_leaves()[0]):
        vetoes = {m for _, _, none_of in rules.values() for m in none_of}
        assert "Ker8_18" not in vetoes
        assert set(EPITHELIAL) <= vetoes


def test_a_panel_without_ker8_18_can_still_call_an_immune_cell():
    assert subpanel.fires(t_cell(), subpanel.immune_rules()["Imm-T-cell"]).all()


def test_ker8_18_alone_does_not_veto():
    assert subpanel.fires(t_cell(Ker8_18=0.9), subpanel.immune_rules()["Imm-T-cell"]).all()


def test_each_shared_epithelial_marker_still_vetoes():
    rule = subpanel.immune_rules()["Imm-T-cell"]
    for marker in EPITHELIAL:
        assert not subpanel.fires(t_cell(**{marker: 0.9}), rule).any()


def test_a_cell_missing_an_acquired_negative_stays_undecided():
    satisfied, unavailable = subpanel.rule_evidence(t_cell(EpCAM=np.nan), subpanel.immune_rules()["Imm-T-cell"])
    assert not satisfied.any() and unavailable.all()
