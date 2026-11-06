import numpy as np
import pandas as pd
import pytest
from shapely.geometry import box
from shapely.affinity import rotate,translate
from multiplex5000.integration.serial_context import (
    geometry_descriptor,neighbor_vectors,compare_neighbor_vectors,assign_geometry,
    assign_multicell_geometry,multicell_objects)


def test_exact_moments_are_translation_and_rotation_invariant():
    rectangle=box(-10,-5,10,5)
    a=geometry_descriptor(rectangle)
    b=geometry_descriptor(translate(rotate(rectangle,37),100000,200000))
    assert a['eccentricity']==pytest.approx(np.sqrt(.75))
    for key in ['area_um2','eccentricity','solidity','circularity']:
        assert a[key]==pytest.approx(b[key],rel=1e-10,abs=1e-10)
    assert geometry_descriptor(box(-10,-10,10,10).difference(box(-5,-5,5,5)))['eccentricity']==pytest.approx(0)


def test_focal_exclusion_preserves_context_with_zero_focal_overlap():
    a=np.array([[0,0],[200,0],[0,400],[-300,50],[100,-250]],float)
    b=a+[1000,500]
    ai,av=neighbor_vectors(a,0);bi,bv=neighbor_vectors(b,0)
    assert 0 not in ai and 0 not in bi
    assert compare_neighbor_vectors(av,bv)['score']==1
    assert box(-10,-10,10,10).intersection(box(990,490,1010,510)).area==0
    wrong=compare_neighbor_vectors(av,bv[:,::-1])
    assert wrong['score']<1


def test_neighbors_use_one_to_one_links_and_report_missing_evidence():
    a=np.array([[100,0],[100,1],[100,2]])
    b=np.array([[100,0],[1000,1000],[-1000,-1000]])
    m=compare_neighbor_vectors(a,b)
    assert m['n_agreeing']==1 and m['score']==pytest.approx(1/3)
    assert np.isnan(compare_neighbor_vectors(a,b[:2])['score'])


def test_geometry_weight_can_prefer_similar_shape_over_nearer_alternative():
    a=pd.DataFrame([dict(structure_id='a',**geometry_descriptor(box(-10,-5,10,5)))])
    b=pd.DataFrame([dict(structure_id='near_different',**geometry_descriptor(box(12,-10,28,10))),
                    dict(structure_id='far_similar',**geometry_descriptor(box(70,-5,90,5)))])
    assert assign_geometry(a,b,position_weight=.75).fixed_id.iloc[0]=='near_different'
    assert assign_geometry(a,b,position_weight=.25).fixed_id.iloc[0]=='far_similar'
    pd.testing.assert_frame_equal(assign_geometry(a,b,position_weight=.25),
                                  assign_geometry(a,b.iloc[::-1],position_weight=.25))


def test_invalid_geometry_and_out_of_range_candidates_are_not_forced():
    with pytest.raises(ValueError):geometry_descriptor(None)
    a=pd.DataFrame([dict(structure_id='a',**geometry_descriptor(box(0,0,10,10)))])
    b=pd.DataFrame([dict(structure_id='b',**geometry_descriptor(box(1000,0,1010,10)))])
    assert assign_geometry(a,b).empty
    with pytest.raises(ValueError):assign_geometry(pd.concat([a,a]),b)


def test_convex_object_solidity_stays_bounded_after_affine_roundoff():
    from shapely.geometry import MultiPoint
    from shapely.affinity import affine_transform
    rng=np.random.default_rng(7)
    for _ in range(30):
        g=affine_transform(MultiPoint(rng.normal(size=(20,2))*100).convex_hull,
                           [.935,.345,-.333,.925,1700,9200])
        descriptor=geometry_descriptor(g)
        assert 0<=descriptor['solidity']<=1
        assert descriptor['solidity']==pytest.approx(1)


@pytest.mark.parametrize('single_side', ['source', 'target', 'both'])
def test_singles_cannot_compete_with_multicell_assignment(single_side):
    def row(key, cls, x):
        return dict(structure_id=key, object_class=cls,
                    n_whole_cells=1 if cls=='EndocrineSingle' else 2,
                    **geometry_descriptor(box(x,0,x+10,10)))
    source=pd.DataFrame([row('rna_cluster','EndocrineCluster',0)])
    target=pd.DataFrame([row('protein_islet','Islet',50)])
    if single_side in ('source','both'):
        source=pd.concat([source,pd.DataFrame([row('rna_single','EndocrineSingle',50)])])
    if single_side in ('target','both'):
        target=pd.concat([target,pd.DataFrame([row('protein_single','EndocrineSingle',0)])])
    # Filtering the unrestricted result afterward loses the multicell pair:
    # each Single is closer and has identical geometry to the opposite object.
    unrestricted=assign_geometry(source,target)
    assert not ((unrestricted.moving_id=='rna_cluster') & (unrestricted.fixed_id=='protein_islet')).any()
    pairs=assign_multicell_geometry(source,target)
    assert list(zip(pairs.moving_id,pairs.fixed_id))==[('rna_cluster','protein_islet')]


def test_single_only_sections_have_no_correspondence_candidates():
    row=dict(structure_id='single',object_class='EndocrineSingle',n_whole_cells=1,
             **geometry_descriptor(box(0,0,10,10)))
    singles=pd.DataFrame([row])
    assert assign_multicell_geometry(singles,singles).empty
    mislabeled=singles.assign(object_class='EndocrineCluster')
    with pytest.raises(ValueError,match='at least two'):multicell_objects(mislabeled)
