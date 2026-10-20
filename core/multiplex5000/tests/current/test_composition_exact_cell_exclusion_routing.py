"""Exclude exact cells without discarding supported members of a Single annotation."""
from pathlib import Path
import sys

import numpy as np
import pandas as pd
import pytest

ROOT = Path(__file__).resolve().parents[4]
sys.path.insert(0, str(ROOT / "core/multiplex5000/scripts/typing"))
import build_endocrine_object_composition as composition


def two_cell_single():
    object_key = "e3541ebb-de44-4f40-b3e5-f0d916f45989"
    excluded_key = "f25dbb5c-a594-4817-b55a-af3c83137f0a"
    alpha_key = "d0d5cadf-1d67-4db8-ac40-b631d44b0f2b"
    cells = pd.DataFrame(dict(object_id=[excluded_key, alpha_key], endocrine_object_id=[object_key] * 2,
        endocrine_object_class=["EndocrineSingle"] * 2, cell_class=["Excluded", "Alpha"],
        refined_class=["Excluded", "Alpha"], is_endocrine=[False, True], qc_analysis_eligible=[True] * 2,
        resolved_by=["excluded_single_hormone", "single_hormone"],
        INS=[np.nan, .01], GCG=[np.nan, .8], SST=[np.nan, .02],
        immune_class=[None] * 2, immune_subtype=[None] * 2))
    row = dict(donor_id="115", endocrine_object_id=object_key, original_object_class="EndocrineSingle",
               n_cells_in_object=2, endocrine_denominator=1, identity_unresolved_count=0,
               primary_count__alpha=1, refined_count__alpha=1, count_endocrine__alpha=1,
               stromal_identity_available_in_donor=False, available__immune__all_immune_per_endocrine=True,
               object_area_um2=123.)
    for group in composition.OUTLINE_STRUCTURES:
        row[f"available__{group}__cell_bodies_per_endocrine"] = True
        row[f"feature__{group}__outline_coverage_fraction"] = .1
    for label in composition.IMMUNE_LABELS:
        row["count_immune__" + label] = 0
    for label in composition.IMMUNE_SUBTYPES:
        row["count_immune_subtype__" + label] = np.nan
    for hormone in ("INS", "GCG", "SST"):
        row["supplementary_available__" + hormone] = True
    exclusions = pd.DataFrame(dict(donor_id=["0115"], object_id=[excluded_key],
        endocrine_object_id=[object_key], exclusion_reason=["endocrine_single_no_usable_hormone"]))
    return pd.DataFrame([row]), cells, exclusions


def test_multicell_single_excludes_only_exact_unsupported_member_and_retains_alpha():
    objects, cells, exclusions = two_cell_single()
    original_objects, original_cells = objects.copy(deep=True), cells.copy(deep=True)
    out = composition.apply_and_rebuild_composition_cell_exclusions(objects, cells, "115", exclusions).iloc[0]
    assert out.n_original_cells_in_object == 2
    assert out.n_analysis_excluded_cells_in_object == 1
    assert out.n_analysis_admitted_cells_in_object == 1
    assert out.excluded_cell_ids == exclusions.object_id.iloc[0]
    assert out.composition_analysis_eligible
    assert out.endocrine_denominator == 1 and out.count_endocrine__alpha == 1
    assert out.original_endocrine_denominator == 1 and out.original_primary_count__alpha == 1
    assert out.endocrine_subtype_composition_complete
    assert out.feature__endocrine__alpha_fraction == 1
    assert out.feature__endocrine__beta_fraction == 0
    assert out.object_area_um2 == 123 and out.feature__vascular__outline_coverage_fraction == .1
    pd.testing.assert_frame_equal(objects, original_objects)
    pd.testing.assert_frame_equal(cells, original_cells)


def test_summary_only_helper_still_rejects_multicell_single():
    objects, _, exclusions = two_cell_single()
    with pytest.raises(ValueError, match="sole-cell Single"):
        composition.apply_composition_cell_exclusions(objects, exclusions)


@pytest.mark.parametrize("failure", ["wrong_object", "missing_cell", "duplicate_cell", "missing_object_key"])
def test_cell_level_routing_rejects_invalid_exact_exclusion_membership(failure):
    objects, cells, exclusions = two_cell_single()
    if failure == "wrong_object":
        exclusions.loc[0, "endocrine_object_id"] = "other-object"
    elif failure == "missing_cell":
        exclusions.loc[0, "object_id"] = "missing-cell"
    elif failure == "duplicate_cell":
        exclusions = pd.concat([exclusions, exclusions], ignore_index=True)
    else:
        exclusions.loc[0, "endocrine_object_id"] = None
    with pytest.raises(ValueError, match="membership|repeat or omit"):
        composition.apply_and_rebuild_composition_cell_exclusions(objects, cells, "115", exclusions)


def test_actual_sole_cell_single_keeps_summary_policy_and_empty_denominator():
    objects, cells, exclusions = two_cell_single()
    cells = cells.iloc[:1].copy()
    objects.loc[0, "n_cells_in_object"] = 1
    objects.loc[0, "endocrine_denominator"] = 0
    out = composition.apply_and_rebuild_composition_cell_exclusions(objects, cells, "115", exclusions).iloc[0]
    assert out.n_original_cells_in_object == out.n_analysis_excluded_cells_in_object == 1
    assert out.original_endocrine_denominator == 1
    assert out.endocrine_denominator == 0 and out.n_analysis_admitted_cells_in_object == 0
    assert not out.composition_analysis_eligible and np.isnan(out.feature__endocrine__alpha_fraction)


def test_other_donor_exclusion_does_not_remove_valid_alpha():
    objects, cells, exclusions = two_cell_single()
    exclusions.loc[0, "donor_id"] = "112"
    exclusions.loc[0, "object_id"] = cells.object_id.iloc[1]
    out = composition.apply_and_rebuild_composition_cell_exclusions(objects, cells, "115", exclusions).iloc[0]
    assert out.endocrine_denominator == out.n_analysis_admitted_cells_in_object == 1
    assert out.count_endocrine__alpha == 1
