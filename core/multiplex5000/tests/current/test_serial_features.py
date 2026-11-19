"""Analytic geometry and missingness checks, independent of historical outcomes."""
import numpy as np
import pytest
import shapely
from multiplex5000.integration.serial_features import (
    overlap_metrics, fractional_links, coverage_fidelity, paired_feature_summary, shift_sensitivity,
    footprint_coverages)


def test_exact_overlap_direction_and_holes():
    a=shapely.box(0,0,10,10);b=shapely.box(5,0,15,10)
    m=overlap_metrics(a,b)
    assert m['intersection_um2']==50 and m['dice']==.5 and m['iou']==pytest.approx(1/3)
    hole=a.difference(shapely.box(2,2,8,8))
    h=overlap_metrics(hole,a)
    assert h['source_coverage']==1 and h['target_coverage']==.64
    assert overlap_metrics(a,shapely.box(20,0,30,10))['dice']==0
    assert np.isnan(overlap_metrics(None,a)['dice'])


def test_split_assignment_conserves_area_without_renormalizing():
    src={'a':shapely.box(0,0,10,10)}
    targets={'left':shapely.box(0,0,3,10),'right':shapely.box(7,0,10,10)}
    edges,summary=fractional_links(src,targets,source_frame='scan:um',target_frame='scan:um')
    assert edges.exclusive_source_fraction.sum()==pytest.approx(.6)
    assert summary.unassigned_fraction.iloc[0]==pytest.approx(.4)
    assert summary.ambiguous_fraction.iloc[0]==0


def test_three_overlapping_targets_do_not_duplicate_source_area():
    src={'a':shapely.box(0,0,10,10)}
    target={'a':shapely.box(0,0,7,10),'b':shapely.box(3,0,10,10),'c':shapely.box(4,0,6,10)}
    e,s=fractional_links(src,target,source_frame='s:um',target_frame='s:um')
    assert s.ambiguous_fraction.iloc[0]==pytest.approx(.4)
    assert e.exclusive_source_fraction.sum()+s.ambiguous_fraction.iloc[0]+s.unassigned_fraction.iloc[0]==pytest.approx(1)
    reverse,_=fractional_links(src,dict(reversed(list(target.items()))),source_frame='s:um',target_frame='s:um')
    assert e.equals(reverse)


def test_missing_target_is_not_observed_absence():
    src={'a':shapely.box(0,0,10,10)}
    _,missing=fractional_links(src,{},source_frame='s:um',target_frame='s:um',target_available=False)
    _,empty=fractional_links(src,{},source_frame='s:um',target_frame='s:um')
    assert np.isnan(missing.unassigned_fraction.iloc[0])
    assert empty.unassigned_fraction.iloc[0]==1
    e,partial=fractional_links(src,{'known':shapely.box(0,0,3,10),'bad':None},source_frame='s:um',target_frame='s:um')
    assert partial.status.iloc[0]=='partial_target_geometry'
    assert partial.covered_fraction.iloc[0]==pytest.approx(.3)
    assert np.isnan(partial.exclusive_fraction.iloc[0])
    assert np.isnan(partial.unassigned_fraction.iloc[0])
    assert np.isnan(partial.ambiguous_fraction.iloc[0])
    assert np.isnan(e.exclusive_source_fraction.iloc[0])


def test_frame_and_invalid_geometry_fail_explicitly():
    with pytest.raises(ValueError,match='frame'):
        fractional_links({}, {},source_frame='rna:um',target_frame='protein:um')
    with pytest.raises(ValueError,match='positive-area'):
        overlap_metrics(shapely.LineString([(0,0),(1,1)]),shapely.box(0,0,1,1))


def test_coverage_fidelity_matches_independently_known_areas():
    row=coverage_fidelity([.8,.2,np.nan],[100,200,50],thresholds=[.5]).iloc[0]
    assert row.true_positive_area_um2==80
    assert row.false_positive_area_um2==pytest.approx(20)
    assert row.false_negative_area_um2==40
    assert row.precision==pytest.approx(.8) and row.recall==pytest.approx(2/3)
    assert row.dice==pytest.approx(160/220) and row.iou==pytest.approx(80/140)
    assert row.n_unavailable==1 and row.available_area_um2==300
    assert np.isnan(coverage_fidelity([np.nan],100).dice).all()


def test_zero_is_measured_and_missing_values_share_denominator():
    row=paired_feature_summary([0,10,np.nan,100],[4,6,7,np.nan],[.5,.5,1,1])
    assert row['n_supported']==2 and row['source_weighted_mean']==5
    assert row['target_weighted_mean']==5 and row['source_minus_target_weighted_mean']==0
    assert np.isnan(paired_feature_summary([1],[2],[0])['source_weighted_mean'])


def test_native_comparison_defaults_to_equal_weights_without_overlap_input():
    row=paired_feature_summary([0,100,np.nan],[4,0,7])
    assert row['n_supported']==2
    assert row['weight_sum']==2
    assert row['source_weighted_mean']==50
    assert row['target_weighted_mean']==2
    assert row['source_minus_target_weighted_mean']==48


def test_shift_bounds_include_baseline_and_expose_small_object_fragility():
    a=shapely.box(0,0,10,10)
    r=shift_sensitivity(a,a,radii_um=[20])[0]
    assert r['dice_max']==1 and r['dice_min']==0
    assert np.isnan(shift_sensitivity(None,a)[0]['dice_min'])


def test_fractional_sampling_unions_overlaps_and_preserves_holes():
    square=shapely.box(0,0,10,10)
    objects=[shapely.box(0,0,6,10),shapely.box(4,0,10,10)]
    assert footprint_coverages([square],objects)[0]==1
    assert footprint_coverages([square],[square.difference(shapely.box(2,2,8,8))])[0]==.64
    assert footprint_coverages([square],[])[0]==0
    assert np.isnan(footprint_coverages([square],None)[0])


@pytest.mark.parametrize('bad',[[-.1],[1.1],[np.inf]])
def test_bad_fraction_is_rejected(bad):
    with pytest.raises(ValueError):coverage_fidelity(bad,100)
    with pytest.raises(ValueError):paired_feature_summary([1],[2],bad)
