from pathlib import Path
import json
import pytest
import pyarrow as pa
import pyarrow.parquet as pq
from analysis.shared import typing_release_view as view


def selection(tmp_path, table, **updates):
    path=tmp_path/'cells.parquet';pq.write_table(table,path)
    scope={'scope_id':'pLN_RNA|serial|section-a','lane':'pLN_RNA','donor_id':'donor',
           'path':str(path),'sha256':view.sha256(path),'cell_key_column':'cell_id',
           'claim_column':'lineage','claim_status':'native_pending_review',
           'QC_column':'cluster_id','QC_kind':'nonnull','filters':{'section_key':'section-a'},
           'expected_source_rows':2,'expected_QC_rows':2}
    scope.update(updates)
    registry={'schema_version':1,'selection_status':'provisional','scopes':[scope],
              'across_scope_counts_are_additive':False,'manifest_bindings':[]}
    rp=tmp_path/'registry.json';rp.write_text(json.dumps(registry))
    return rp,scope,path


def test_pooled_sections_with_repeated_ids_are_filtered_before_counts(tmp_path):
    t=pa.table({'cell_id':['same','other','same','other'],
                'section_key':['section-a','section-a','section-b','section-b'],
                'lineage':['Unresolved','T','B','B'],'cluster_id':['a0','a1','b0','b0']})
    rp,s,_=selection(tmp_path,t)
    view.load_registry(rp)
    c=view.scope_counts(s)
    assert c['source_rows']==2
    assert c['QC_claim_counts']==[{'claim':'T','n':1},{'claim':'Unresolved','n':1}]
    assert c['counts_across_scopes_additive'] is False


def test_pooled_pLN_registry_cannot_omit_section_filter(tmp_path):
    t=pa.table({'cell_id':['same'],'section_key':['section-a'],'lineage':['T'],'cluster_id':['a0']})
    rp,_,_=selection(tmp_path,t,filters={})
    with pytest.raises(ValueError,match='exact section filter'):
        view.load_registry(rp)


def test_unknown_identity_is_retained_in_native_QC_denominator(tmp_path):
    t=pa.table({'cell_id':['a','b','c'],'lineage':['T','Unresolved','QC_failed'],
                'cluster_id':['a0','a1',None],'known_lineage_admitted':[True,False,False]})
    _,s,_=selection(tmp_path,t,lane='pancreas_RNA',filters={},expected_source_rows=3,
                    additional_admission_column='known_lineage_admitted',expected_additional_admission_rows=1,
                    additional_admission_semantics='Known lineage only, distinct from native QC')
    c=view.scope_counts(s)
    assert c['source_rows']==3 and c['QC_rows']==2 and c['additional_admission_rows']==1
    assert {'claim':'Unresolved','n':1} in c['QC_claim_counts']


def test_changed_selected_table_is_rejected(tmp_path):
    t=pa.table({'cell_id':['a','b'],'section_key':['section-a','section-a'],
                'lineage':['T','Unresolved'],'cluster_id':['a0','a1']})
    _,s,path=selection(tmp_path,t)
    pq.write_table(t.set_column(2,'lineage',pa.array(['B','B'])),path)
    with pytest.raises(ValueError,match='Selected cell table changed'):
        view.scope_counts(s)


def test_additional_admission_cannot_bypass_native_QC(tmp_path):
    t=pa.table({'cell_id':['a','b'],'lineage':['T','T'],'cluster_id':['a0',None],
                'known_lineage_admitted':[True,True]})
    _,s,_=selection(tmp_path,t,lane='pancreas_RNA',filters={},expected_QC_rows=1,
                    additional_admission_column='known_lineage_admitted',expected_additional_admission_rows=2)
    with pytest.raises(ValueError,match='cannot admit native QC-ineligible'):
        view.scope_counts(s)


def test_section_candidate_requires_matching_native_scope_and_qualified_keys(tmp_path):
    t=pa.table({'cell_id':['same','other'], 'serial_qualified_cell_key':['pLN|section-a|001|same','pLN|section-a|001|other'],
                'section_id':['section-a','section-a'],'lineage':['T','Unresolved'],'cluster_id':['a0','a1']})
    rp,s,_=selection(tmp_path,t,storage_scope='section_specific_candidate',section_id='section-a',source_serial='001',
                     cell_key_column='serial_qualified_cell_key',filters={'section_id':'section-a'},
                     native_source_filters={'section_key':'section-a'})
    view.load_registry(rp)
    assert view.scope_counts(s)['QC_rows']==2
    registry=json.loads(rp.read_text())
    registry['scopes'][0]['native_source_filters']={'section_key':'section-b'}
    rp.write_text(json.dumps(registry))
    with pytest.raises(ValueError,match='exact native scope'):
        view.load_registry(rp)
