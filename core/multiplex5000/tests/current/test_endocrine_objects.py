import json
import numpy as np
import pandas as pd
import pytest
from shapely.geometry import box, Polygon, mapping
from multiplex5000.integration.endocrine_objects import (
    rna_proposals, vertex_polygons, designate, mean_caliper_um,
    native_uuid, validate_constructed, execute_constructor, digest,
)


def cells():
    return pd.DataFrame(dict(cell_id=['one','two','three','no_nucleus','no_depth'],
        x_um=[5,105,115,205,305],y_um=[5]*5,nucleus_count=[1,1,1,0,1],
        depth=[20,30,40,30,0],seed_evidence=[True]*5))


def test_sparse_groups_and_actual_single_geometry():
    c=cells();geom=pd.DataFrame(dict(cell_id=c.cell_id,geometry=[box(x-5,0,x+5,10) for x in c.x_um]))
    proposals,seeds=rna_proposals(c,geom,acquisition_id='section-A',eps_um=20)
    assert sorted(proposals.n_seed_cells)==[1,2]
    assert set(seeds.cell_id)=={'one','two','three'}
    single=proposals.loc[proposals.n_seed_cells.eq(1)].iloc[0]
    assert single.geometry.equals_exact(geom.geometry.iloc[0],0)
    p2,_=rna_proposals(c.sample(frac=1,random_state=4),geom,acquisition_id='section-A',eps_um=20)
    assert set(proposals.proposal_id)==set(p2.proposal_id)
    p3,_=rna_proposals(c,geom,acquisition_id='section-B',eps_um=20)
    assert not set(p3.proposal_id)&set(proposals.proposal_id)


def test_missing_cell_geometry_does_not_become_a_circle():
    c=cells();g=pd.DataFrame(dict(cell_id=['one'],geometry=[box(0,0,10,10)]))
    with pytest.raises(ValueError,match='segmentation missing'):
        rna_proposals(c,g,acquisition_id='a')


def test_pubapp_category_boundary_and_whole_cell_count():
    assert designate(1,100)=='EndocrineSingle'
    assert designate(2,34.999)=='EndocrineCluster'
    assert designate(2,35)=='Islet'
    with pytest.raises(ValueError):designate(0,10)
    # Rectangle: diameter of minimum circle is diagonal, min caliper is height.
    assert mean_caliper_um(box(0,0,30,40))==pytest.approx(40)


def test_vendor_order_and_dangling_line_preserve_area():
    v=pd.DataFrame(dict(cell_id=['a']*7,vertex_x=[0,1,1,0,0,0,0],vertex_y=[0,0,1,1,0,-1,0]))
    p=vertex_polygons(v)
    assert p.geometry.iloc[0].area==pytest.approx(1)
    bad=v.copy();bad.loc[3,'cell_id']='b'
    with pytest.raises(ValueError,match='contiguous'):vertex_polygons(bad)


def test_single_geometry_and_unique_member_checks():
    g=box(0,0,10,10)
    feature=dict(id='object',geometry=mapping(g),properties=dict(classification=dict(name='EndocrineSingle'),measurements={'Whole cells':1}))
    native=pd.DataFrame(dict(cell_id=['native'],geometry=[g]))
    ob,links,qc=validate_constructed([feature],native,pixel_size_um=1)
    assert ob.object_class.tolist()==['EndocrineSingle']
    assert links.cell_id.tolist()==['native']
    assert qc['n_single_geometry_disagreements']==0
    second={**feature,'id':'duplicate'}
    _,_,qc=validate_constructed([feature,second],native,pixel_size_um=1)
    assert qc['n_multiply_owned_cells']==1


def test_native_uuid_binds_acquisition_and_preserves_native_sidecar():
    assert native_uuid('a','same')!=native_uuid('b','same')
    assert native_uuid('a','same')==native_uuid('a','same')


def test_constructor_refuses_tampered_source_before_execution(tmp_path):
    f=tmp_path/'source';f.write_text('original');h=digest(f);f.write_text('changed')
    config=tmp_path/'config.json';config.write_text('{}')
    with pytest.raises(ValueError,match='Pinned source changed'):
        execute_constructor(config_path=config,project_path=tmp_path/'project.qpproj',qupath='/no/engine',
                            source_hashes={str(f):h},log_path=tmp_path/'log')


def test_precision_bound_does_not_hide_single_boundary_change():
    g=box(1000,1000,1010,1010)
    altered=box(1000.1,1000,1010,1010)
    f=dict(id='a',geometry=mapping(altered),properties=dict(classification=dict(name='EndocrineSingle'),measurements={'Whole cells':1}))
    _,_,qc=validate_constructed([f],pd.DataFrame(dict(cell_id=['native'],geometry=[g])),pixel_size_um=1)
    assert qc['n_single_geometry_disagreements']==1


def test_replay_entry_rejects_unreviewed_constructor_before_writing(tmp_path):
    from multiplex5000.integration.endocrine_workflow import run_manifest
    manifest=tmp_path/'input.json'
    manifest.write_text(json.dumps(dict(schema='sequential_whole_endocrine_v1',sections=[],constructor={'sha256':'different'})))
    with pytest.raises(ValueError,match='Unreviewed PubApp constructor'):
        run_manifest(manifest,tmp_path/'result',qupath='/no/engine')
    assert not (tmp_path/'result').exists()


def test_lineage_keeps_competing_parents_and_does_not_count_boundary_touch():
    from multiplex5000.integration.endocrine_objects import proposal_lineage
    o=pd.DataFrame(dict(structure_id=['split'],geometry=[box(0,0,10,10)]))
    p=pd.DataFrame(dict(proposal_id=['first','second','touch'],geometry=[box(0,0,5,10),box(4,0,7,10),box(10,0,20,10)]))
    links=proposal_lineage(o,p)
    assert set(links.original_proposal_id)=={'first','second'}
