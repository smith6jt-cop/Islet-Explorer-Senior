import importlib.util
import json
from pathlib import Path
import sys

import numpy as np
import pandas as pd
import pytest
from scipy import sparse
from scipy.stats import poisson

HERE = Path(__file__).resolve().parents[1]
sys.path.insert(0,str(HERE)); sys.path.insert(0,str(HERE.parent/'xenium_celltype'))
import evidence as e
import gdna
import manifest
import workflow


def test_aliases_and_section_isolation():
    assert {manifest.canonical_donor(x) for x in ('HDL112','0112','112')} == {'112'}
    m = {'sections':[{'section_id':'112/panc/0081366'},{'section_id':'112/pln_1/0081366'}]}
    assert manifest.get_section(m,'112/panc/0081366')['section_id'].split('/')[1]=='panc'
    with pytest.raises(ValueError):manifest.get_section(m,'112')


def test_prior_zero_controls_and_ties():
    k=np.array([0,0,1,2]); lam,raw,_,a=gdna.background_rates(k,np.ones(4),21)
    assert a==21 and (lam[:2]>0).all() and (raw[:2]==0).all()
    reverse=gdna.background_rates(k[::-1],np.ones(4),21)[0][::-1]
    np.testing.assert_array_equal(lam,reverse)
    assert np.isfinite(gdna.background_rates([0,1],[np.nan,2],21)[0]).all()
    with pytest.raises(gdna.GdnaError):gdna.background_rates(k,np.zeros(4),21)
    with pytest.raises(gdna.GdnaError):gdna.background_rates(k,np.ones(4),0)
    with pytest.raises(gdna.GdnaError):gdna.background_rates([-1], [1],21)


def test_vendor_rejects_missing_invalid_and_disagreement():
    for metric in ({}, {gdna.VENDOR_FP_KEY:float('nan')},{gdna.VENDOR_FP_KEY:100}):
        with pytest.raises(gdna.GdnaError):gdna.verify_against_vendor(.01,5001,metric)
    assert gdna.verify_against_vendor(.01,5001,{gdna.VENDOR_FP_KEY:50.01})['relative_error']==0


def test_detection_extreme_counts_and_preserves_raw():
    x=sparse.csr_matrix([[0,1,2,1000000],[1,2,100,1000000]],dtype=np.int64)
    original=x.copy(); rates=np.array([.1,10000.])
    d=e.detected(x,rates)
    np.testing.assert_array_equal(d.toarray(),poisson.sf(x.toarray()-1,rates[:,None])<=.003)
    assert (x!=original).nnz==0 and x.dtype==np.int64
    with pytest.raises(ValueError):e.detected(sparse.csr_matrix([[1.5]]),[.1])
    with pytest.raises(ValueError):e.detected(x,[-1,0])


def test_independent_background_simulation_core_sizes():
    rng=np.random.default_rng(418)
    for size in (100,1000,10000):
        lam=np.linspace(.005,.3,size)
        x=sparse.csr_matrix(rng.poisson(lam[:,None],size=(size,100)))
        rate=e.detected(x,lam).nnz/(size*100)
        assert rate < .003 + 4*np.sqrt(.003/(size*100))


def test_nucleus_qc_counts_and_exclusions():
    x=sparse.csr_matrix(np.full((4,25),4,dtype=np.int32))
    cells=pd.DataFrame({'cell_area':[100,100,100,0],'nucleus_area':[20,0,20,20],
                        'total_counts':[110,110,110,110]})
    q,meta=e.qc_table(x,cells,np.array([.001,.001,.1,.1]),5001)
    assert q.qc_eligible.tolist()==[True,False,False,False]
    assert 'no_nucleus' in q.iloc[1].exclusion_reasons
    assert q.iloc[2].background_dominated and q.iloc[2].adjusted_qc_total<0
    assert q.vendor_total_counts.tolist()==[110]*4 and x.dtype==np.int32


def test_sex_controls_applicability_and_missing():
    x=sparse.csr_matrix(np.ones((100,3),int)*10)
    cells=pd.DataFrame({'raw_gene_counts':np.arange(100),'gdna_lambda':np.ones(100)*.01})
    d=e.detected(x,cells.gdna_lambda)
    female=e.control_report(x,d,e.Y_CLEAN,cells,'F','pancreas')
    male=e.control_report(x,d,e.Y_CLEAN,cells,'M','pancreatic_lymph_node')
    assert 'fail' in set(female.status) and 'fail' not in set(male.status)
    assert 'female_Y_negative_not_applicable' in set(male.status)
    assert (female[female.probe=='DDX3Y'].status=='missing_probe').all()
    assert 'NLRP7' not in set(male.probe)
    assert set(e.Y_CLEAN).issubset(e.SEX_GENES)


def test_partial_stage_is_not_success_and_output_tampering(tmp_path):
    p=tmp_path/'result';p.write_text('complete')
    rec=tmp_path/'current.json'
    workflow.write_json(rec,{'status':'running'})
    assert not workflow.valid_receipt(rec)
    workflow.write_json(rec,{'status':'validated_complete','fingerprint':'a',
         'outputs':[workflow.identity(p)],'validation':{'ok':True}})
    assert workflow.valid_receipt(rec,'a') and not workflow.valid_receipt(rec,'b')
    p.write_text('changed!');assert not workflow.valid_receipt(rec,'a')


def test_signed_excess_not_clipped():
    x=sparse.csr_matrix([[0,10],[0,10]])
    s=e.group_summary(x,e.detected(x,[1,1]),[1,1],['low','high'])
    assert s.iloc[0].signed_excess==-2 and s.iloc[0].below_section_background
    assert s.iloc[1].signed_excess==18


def test_overlapping_tissue_objects_no_double_count():
    from shapely.geometry import box
    from shapely import union_all
    from prepare_protein import partition_points
    xy=np.array([[1,1],[3,3],[5,5],[20,20]])
    tissue=box(0,0,10,10);nodes=union_all([box(2,2,4,4),box(3,3,6,6)])
    assert partition_points(xy,tissue,nodes).tolist()==['panc','pln_1','pln_1','outside']


def test_finite_measurement_transfer_weights():
    from build_joined import aggregate_features
    edges=pd.DataFrame({'cell_id':['x','x','y'],'object_id':['a','b','b'],'overlap_area_um2':[1.,3.,2.]})
    ph=pd.DataFrame({'object_id':['a','b'],'feat__Cell__CD3':[10.,np.nan],'feat__Cell__CD20':[4.,8.]})
    out,_=aggregate_features(['x','y'],edges,ph)
    assert out['feat__Cell__CD3'].iloc[0]==10 and np.isnan(out['feat__Cell__CD3'].iloc[1])
    assert out['feat__Cell__CD20'].tolist()==[7.,8.]


def test_reseed_footprint_keeps_only_seeded_xenium_area():
    import reseed_four
    rng=np.random.default_rng(0)
    moving=rng.uniform(0,1000,(5000,2)); fixed=rng.uniform(-2000,3000,(50000,2))
    keep=reseed_four.footprint_mask(moving,fixed,100.0)
    inside=(fixed>-100).all(1)&(fixed<1100).all(1); far=((fixed<-300)|(fixed>1300)).any(1)
    assert keep[inside].mean()>0.95 and not keep[far].any()


def test_join_dominant_partner_breaks_ties_by_uuid():
    from join_four import dominant_partner
    edges=pd.DataFrame({'cell_id':['x','x','y','y'],'object_id':['b','a','c','d'],'overlap_area_um2':[2.,2.,1.,3.]})
    ids,share=dominant_partner(pd.Index(['x','y','z']),edges)
    assert ids[:2].tolist()==['a','d'] and pd.isna(ids[2])
    assert share[:2].tolist()==[.5,.75] and np.isnan(share[2])
