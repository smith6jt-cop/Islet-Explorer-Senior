"""Composition metadata must join numeric IDs across zero-padded partitions."""

from __future__ import annotations

import sys
from pathlib import Path

HERE = Path(__file__).resolve()
LANE = HERE.parents[1]
sys.path.insert(0, str(LANE))

import build_endocrine_object_composition as composition  # noqa: E402
import pandas as pd


def test_inventory_status_matches_zero_padded_typing_partition():
    status_by_donor = composition.disease_status_map([
        {"donor_id": "112", "disease_status": "ND"},
        {"donor_id": "115", "disease_status": "ND"},
    ])

    assert status_by_donor[composition.donor_id_key("0112")] == "ND"
    assert status_by_donor[composition.donor_id_key("0115")] == "ND"


def test_inventory_status_rejects_conflicting_records_after_normalization():
    try:
        composition.disease_status_map([
            {"donor_id": "112", "disease_status": "ND"},
            {"donor_id": "0112", "disease_status": "T1D"},
        ])
    except ValueError as exc:
        assert "inconsistent disease status" in str(exc)
    else:
        raise AssertionError("conflicting status records were silently accepted")


def test_candidate_identity_never_restores_current_cell_analysis_eligibility():
    cells = pd.DataFrame({
        'object_id': ['valid', 'reviewed', 'single', 'explicit', 'failed_QC', 'current_excluded'],
        'cell_class': ['Alpha', 'Excluded', 'Endocrine (ambiguous subtype)', 'Alpha', 'Immune', 'Immune'],
        'qc_analysis_eligible': [True, True, True, True, False, True],
        'analysis_eligible': [True, False, False, True, False, False],
        'reviewed_identity': [None, 'Immune', None, None, None, None],
    })
    explicit = pd.DataFrame({'donor_id': ['112'], 'object_id': ['explicit']})
    mask = composition.composition_cell_analysis_mask(cells, '0112', explicit)
    assert cells.loc[mask, 'object_id'].tolist() == ['valid']
