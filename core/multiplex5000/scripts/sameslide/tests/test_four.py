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


def test_delivery_checksums_and_donor_spelling(tmp_path):
    import import_alignment as IA
    (tmp_path/'a.json').write_text('{}')
    digest=workflow.sha(tmp_path/'a.json')
    (tmp_path/'SHA256SUMS').write_text(f'{digest}  ./a.json\n{"0"*64}  ./b.json\n')
    sums=IA.checksums(tmp_path)
    assert IA.verified(tmp_path,'a.json',sums)['sha256']==digest
    (tmp_path/'b.json').write_text('{}')
    with pytest.raises(ValueError):IA.verified(tmp_path,'b.json',sums)
    assert IA.delivery_donor(['0112','115','6347'],'112')=='0112'
    with pytest.raises(ValueError):IA.delivery_donor(['115'],'112')


def test_alignment_frames_never_share_a_column():
    import import_alignment as IA
    lane=pd.DataFrame({'object_id':['p','n','gone'],'qc_analysis_eligible':[True,True,True],
                       'nucleus_area_um2':[10.,10.,10.]})
    cohort=pd.DataFrame({'object_id':['p','n'],'qc_analysis_eligible':[True,True],'cell_region':['tissue','pLN'],
                         'islet_num':['',''],'in_alignment_block':[True,False],'CD3e':[5.,np.nan]})
    node=pd.DataFrame({'object_id':['n'],'qc_analysis_eligible_with_nucleus':[True],'annotation_uuid':['u'],'CD3e':[7.]})
    t=IA.frame_table(lane,cohort,node,['CD3e'],['CD3e']).set_index('object_id')
    assert t.alignment_frame.tolist()==['pancreas_cohort','pln_tissue','']
    assert t.loc['p','CD3e']==5 and np.isnan(t.loc['p','pln__CD3e'])
    assert np.isnan(t.loc['n','CD3e']) and t.loc['n','pln__CD3e']==7
    # A prepared cell the delivery lacks is unavailable, never zero, and never eligible.
    assert t.loc['gone',['CD3e','pln__CD3e']].isna().all() and not t.loc['gone','analysis_eligible']
    both=node.assign(object_id='p')
    with pytest.raises(ValueError):IA.frame_table(lane,cohort,both,['CD3e'],['CD3e'])


def test_join_records_contributing_frames():
    from join_four import frame_of_cells
    xi=np.array([0,0,1,2,2]);frames=np.array(['pancreas_cohort','pancreas_cohort','pln_tissue','pln_tissue','pancreas_cohort'])
    assert frame_of_cells(xi,frames,4).tolist()==['pancreas_cohort','pln_tissue','mixed','']
    # The same rule names the tissue a cell's contributors lie in: a node a pancreas region images is pLN.
    tissue=np.array(['pancreas','pancreas','pLN','pLN','pancreas'])
    assert frame_of_cells(xi,tissue,4,('pancreas','pLN')).tolist()==['pancreas','pLN','mixed','']


def test_delivered_calls_read_by_donor_spelling_and_refuse_duplicates(tmp_path):
    import calls_four
    pd.DataFrame({'donor_id':['0112','0112','6580'],'object_id':['a','b','c'],'cell_class':['Beta','Excluded','Immune']}
                 ).to_parquet(tmp_path/'calls.parquet')
    calls,records=calls_four.delivered_calls(tmp_path,'112')
    assert calls.object_id.tolist()==['a','b'] and len(records)==1
    with pytest.raises(ValueError):calls_four.delivered_calls(tmp_path,'6347')
    pd.DataFrame({'donor_id':['0112'],'object_id':['a'],'cell_class':['Alpha']}).to_parquet(tmp_path/'again.parquet')
    with pytest.raises(ValueError):calls_four.delivered_calls(tmp_path,'112')


def test_one_class_and_its_levels_only_inside_their_group(tmp_path):
    import calls_four
    t=pd.DataFrame({'object_id':list('abcd'),'refined_class':['Immune','Exocrine','Immune',None],
                    'immune_class':['Imm-T-cell','Imm-macrophage',None,None]})
    out,record=calls_four.one_class(t,'refined_class',{'Immune':['immune_class']})
    # One class column; a subclass outside its group is dropped, a group without one stays unresolved.
    assert out.call.tolist()[:3]==['Immune','Exocrine','Immune'] and pd.isna(out.call.iloc[3])
    assert out.immune_class.tolist()[0]=='Imm-T-cell' and out.immune_class.iloc[1:].isna().all()
    assert record['immune_class']=={'group':'Immune','set_outside_group':1,'missing_inside_group':1}
    for kind,rows in (('pancreas',['a']),('nodes',['b'])):
        (tmp_path/kind).mkdir();pd.DataFrame({'donor_id':['0112'],'object_id':rows}).to_parquet(tmp_path/kind/'data.parquet')
    assert calls_four.delivered_calls(tmp_path,'112','pancreas')[0].object_id.tolist()==['a']


def test_each_targeted_level_belongs_to_one_table():
    import calls_four
    ok={'tables':{'pancreas':{'tissue':'pancreas','class_column':'refined_class','within':{'Immune':['immune_class']},'carry':['producer']},
                  'nodes':{'tissue':'pLN','class_column':'immune_class','within':{'T_cell':['t_nk_class']},'carry':['producer']}}}
    assert set(calls_four.table_specs(ok))=={'pancreas','nodes'}
    with pytest.raises(ValueError):calls_four.table_specs({'tables':{'a':{'within':{'X':['level']}},'b':{'within':{'Y':['level']}}}})
    with pytest.raises(ValueError):calls_four.table_specs({'tables':{'a':{'within':{'X':['level']}},'b':{'carry':['level']}}})


def test_node_cells_without_a_class_are_counted_by_the_producers_reason():
    import nodes_four
    out=nodes_four.why(['T_cell','','','','Uncalled',None],[True,False,pd.NA,True,True,None],{'','Uncalled'})
    assert out[0]=='' and out[1].startswith('not immune') and out[2]=='no lineage marker named'
    assert out[3]=='immune, no label' and out[4].startswith('Uncalled') and out[5]=='protein call unavailable'


def test_node_depth_is_positive_inside_and_negative_outside():
    import shapely
    import nodes_four
    d=nodes_four.depth(np.array([[50.,50],[10,50],[-5,50]]),shapely.box(0,0,100,100))
    assert np.allclose(d,[50,10,-5],atol=nodes_four.EDGE_STEP_UM/2)


def test_node_summaries_apply_to_a_pancreas_section_that_imaged_a_node(tmp_path):
    import completion_four as C
    panc={'section_id':'1/panc/9','donor_id':'1','roi':'panc','xenium_serial':'9','tissue':'pancreas','design':'same_slide'}
    assert C.applies(tmp_path,panc,'islets') and not C.applies(tmp_path,panc,'node_summaries')
    calls=C.receipt_path(tmp_path,panc,'calls');calls.parent.mkdir(parents=True)
    calls.write_text(json.dumps({'summary':{'status_by_winning_tissue':{'pancreas|awarded':5}}}))
    assert not C.applies(tmp_path,panc,'node_summaries')
    calls.write_text(json.dumps({'summary':{'status_by_winning_tissue':{'pancreas|awarded':5,'pLN|awarded':2}}}))
    assert C.applies(tmp_path,panc,'node_summaries')
    pln={**panc,'roi':'pln_1','tissue':'pancreatic_lymph_node'}
    assert C.applies(tmp_path,pln,'node_summaries') and not C.applies(tmp_path,pln,'islets')


def test_figure_call_status_covers_every_cell_and_names_the_gate_windows():
    import figures_four as F
    summary={'n_xenium_cells':100,'status_by_winning_tissue':{'pancreas|awarded':50,'pLN|awarded':20,'none|unavailable':4,
             'pLN|unavailable':1,'pancreas|not_a_class':2,'pLN|not_a_class':23},
             'not_a_class_by_label':{'pancreas|Excluded|nan':2,'nodes|""|False':15,'nodes|""|<NA>':3,'nodes|Uncalled|True':5}}
    s=F.call_status(summary)
    assert abs(sum(s.values())-1)<1e-12 and s['node: not immune']==.15 and s['no protein partner']==.04
    with pytest.raises(ValueError):F.call_status({**summary,'n_xenium_cells':101})
    assert F.short('only 7/14 islets (50%) were measurable; a verdict needs 60%')=='DAPI gate: 7 of 14 windows measurable, a verdict needs 60%'
    assert F.short('join: registration: no overlap; registration_recheck: x')=='no overlap'


def test_figure_writes_a_part_too_thin_to_see_beneath_its_bar():
    import matplotlib.pyplot as plt
    import figures_four as F
    assert [F.share_text(w) for w in (.0014,.01,1/13034,.45)]==['0.14 %','1 %','0.0077 %','45 %']
    fig,ax=plt.subplots()
    F.stacked(ax,[{'Exocrine':.9986,'Lymphatic':.0014},{'Exocrine':.995,'Nerve':.005,'Lymphatic':0.},{'Exocrine':.98,'Nerve':.02},None],
              ['Exocrine','Lymphatic','Nerve'],{'Exocrine':'k','Lymphatic':'m','Nerve':'y'},names={'Nerve':'nerve'})
    assert [t.get_text() for t in ax.texts]==['under 1 %: Lymphatic (0.14 %)','under 1 %: nerve (0.5 %)']
    plt.close(fig)


def test_islet_regions_come_from_the_signed_distance():
    import shapely
    import islets_four as I
    polys=np.array([shapely.box(0,0,100,100),shapely.box(300,0,310,10)])
    xy=np.array([[50.,50],[99,50],[110,50],[120.5,50],[305,5],[200,50]])
    d,idx=I.signed_distance(xy,polys)
    assert np.allclose(d,[-50,-1,10,20.5,-5,100]) and idx.tolist()==[0,0,0,0,1,0]
    # Core or tissue, as the approved pipeline: no expansion, however close to the outline.
    assert I.region_of(np.array([-3.,0,.01,20,20.01])).tolist()==['core','core','tissue','tissue','tissue']
    # Coverage reads the outline alone: a fully imaged object reads 1, one imaged on one side only less.
    grid=np.array([[x,y] for x in range(0,101,5) for y in range(0,101,5)],float)
    assert I.coverage(polys[:1],grid)[0]==1.0
    assert I.coverage(polys[:1],grid[grid[:,0]<50])[0]<.6
    cells=pd.DataFrame({'region':['core','tissue'],'object_class':['Islet','Islet'],'call':['Beta','Exocrine']})
    assert I.in_object(cells).object_class.isna().tolist()==[False,True]


def test_association_auroc_shift_and_quantile_ties():
    from scipy.stats import rankdata
    import associate_four as AS
    v=np.array([1.,2,3,4]);assert AS.auroc(rankdata(v),[0,0,1,1])==1 and AS.auroc(rankdata(v),[1,1,0,0])==0
    assert AS.auroc(rankdata(np.ones(4)),[0,1,0,1])==.5 and np.isnan(AS.auroc(rankdata(v),[1,1,1,1]))
    xy=np.array([[0.,0],[25,0],[50,0],[0,40]])
    assert AS.shifted_partner(xy).tolist()==[1,2,-1,-1]
    rows=AS.conditional(np.r_[np.zeros(50),np.arange(1.,51)],np.r_[np.zeros(50,bool),np.ones(50,bool)])
    # Fifty equal values are one bin, never split by row order.
    assert rows[0]['n']==50 and rows[0]['n_detected']==0 and sum(r['n'] for r in rows)==100


def test_comparison_ratios_need_detection_in_both_sections(tmp_path):
    import compare_sections as CS
    for design,rates in (('same_slide',[.2,.002,.0]),('xenium_only',[.1,.2,.0])):
        d=tmp_path/design;d.mkdir()
        pd.DataFrame({'gene':['A','B','C'],'n_cells':100,'observed_counts':[40,5,0],'signed_excess':[30,1,-1],
                      'corrected_detection_rate':rates,'raw_detection_rate':rates,'below_section_background':False}
                     ).to_parquet(d/'gene_summary.parquet')
    g,one=CS.gene_table({'same_slide':tmp_path/'same_slide','xenium_only':tmp_path/'xenium_only'},{'g':['A','B']})
    g=g.set_index('gene');assert one=={'same_slide':[],'xenium_only':[]}
    assert g.loc['A','log2_rate_ratio']==1 and g[['log2_rate_ratio']].loc[['B','C']].isna().all().all()
    assert g.loc['A','identity_groups']=='g' and g.loc['C','identity_groups']==''
    groups=CS.identity_groups('pancreas')
    assert 'pancreas_endocrine' in groups and 'pancreas_exocrine' in groups and 'immune_T_cell' in groups
    assert 'pancreas_endocrine' not in CS.identity_groups('pancreatic_lymph_node')


def test_recheck_reads_the_gate_that_judged_the_transform(tmp_path):
    import recheck_nuclei as RN
    (tmp_path/'alignment_four_refined_odd.json').write_text(json.dumps({'alignment':{'threshold_pass_um':1.65}}))
    (tmp_path/'alignment_other.json').write_text(json.dumps({'alignment':{'threshold_pass_um':9.0}}))
    assert RN.gate_threshold(tmp_path,'four_refined')[0]==1.65
    with pytest.raises(ValueError):RN.gate_threshold(tmp_path,'missing')


def test_completion_staleness_is_data_not_code_and_never_self(tmp_path):
    import completion_four as C
    data=tmp_path/'upstream.json';data.write_text('{}');code=tmp_path/'stage.py';code.write_text('x=1')
    own=tmp_path/'current.json';own.write_text('{}')
    r={'inputs':{'up':workflow.identity(data),'code':workflow.identity(code),'source':workflow.identity(own)}}
    own.write_text('{"rewritten": true}')
    assert C.changed(r,own=own)==([],[])
    code.write_text('x=20');assert C.changed(r,own=own)==([],['changed '+tmp_path.name+'/stage.py'])
    data.write_text('{"new": 1}');assert C.changed(r,own=own)[0]==['changed '+tmp_path.name+'/upstream.json']


def test_completion_marks_a_stage_of_the_other_tissue_not_applicable():
    import completion_four as C
    assert C.TISSUE_ONLY['islets']=='pancreas' and C.TISSUE_ONLY['node_summaries']=='pancreatic_lymph_node'
    assert set(C.TISSUE_ONLY)<=set(C.SAME_SLIDE) and 'islets' in C.SECTION_STAGES


def test_a_reprepared_segmentation_is_answered_by_the_recheck_not_staleness(tmp_path):
    import completion_four as C
    prep=tmp_path/'preparation';prep.mkdir();cur=prep/'current.json';cur.write_text('{}')
    reg=tmp_path/'registration';reg.mkdir()
    (reg/'current.json').write_text(json.dumps({'status':'validated_complete','inputs':{'preparation':workflow.identity(cur)}}))
    cur.write_text('{"re-prepared": true}')
    assert C.cell(reg/'current.json')['status']=='validated_complete'
