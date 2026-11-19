#!/usr/bin/env python3
"""Synthetic contract/geometry/source-preservation tests; no real scientific data."""
import contextlib
import io
import json
from pathlib import Path
import tempfile
from types import SimpleNamespace
import unittest

import h5py
import numpy as np
import pandas as pd
import pyarrow as pa
import pyarrow.parquet as pq
from scipy import sparse

from multiplex5000.gammadelta.extract_corrected import POPULATION, extract, load_contract, sha, signature
from multiplex5000.gammadelta.validate_coordinates import independent_roi, validate

HERE=Path(tempfile.gettempdir())

def write(path,obj):
    path.write_text(json.dumps(obj,indent=2)+'\n')

def fixture(root,empty=False,duplicate=False,metadata_csv=False,nonfinite_assigned=False,fractional_centroids=False):
    bundle=root/'native';bundle.mkdir()
    sel=root/'selector';sel.mkdir();(sel/'validation').mkdir()
    names=['TRDC','TRGC2','CD3E','CD3G','NegativeControl']
    # The two candidates each express only ONE receptor and ONE CD3 gene.
    # The old TRDC AND TRGC2 gate would lose both; no such gate is allowed here.
    native=sparse.csc_matrix(np.array([[2,0,0],[0,2,0],[2,0,0],[0,2,0],[1,0,0]],dtype=np.int32))
    with h5py.File(bundle/'cell_feature_matrix.h5','w') as h:
        m=h.create_group('matrix')
        for k,v in [('shape',native.shape),('data',native.data),('indices',native.indices),('indptr',native.indptr)]:m.create_dataset(k,data=v)
        m.create_dataset('barcodes',data=np.asarray([b'a',b'b',b'c']))
        f=m.create_group('features');f.create_dataset('name',data=np.asarray([n.encode() for n in names]));f.create_dataset('id',data=np.asarray([n.encode() for n in names]));f.create_dataset('feature_type',data=np.asarray([b'Gene Expression']*4+[b'Negative Control Probe']))
    centers=[0.12345678901234568,1000.1234567890123] if fractional_centroids else [0.,1000.]
    metadata=pd.DataFrame({'cell_id':['a','b','c'],'x_centroid':centers+[5000.],'y_centroid':centers+[5000.],'nucleus_count':[2,0,1]})
    meta_name='cells.csv.gz' if metadata_csv else 'cells.parquet'
    if metadata_csv:metadata.to_csv(bundle/meta_name,index=False,compression='gzip')
    else:pq.write_table(pa.Table.from_pandas(metadata,preserve_index=False),bundle/meta_name)
    # Rows cover inclusive corners, just-outside ROI, low QV, unknown background
    # features, flagged assigned rows far outside ROI, native null/extra columns.
    source=pa.table({
        'transcript_id':pa.array([1,2,3,4,5,6,7,8,9,10,11,12,1 if duplicate else 13],type=pa.uint64()),
        'cell_id':['a','a','a','a','a','a','b','b','b','b','UNASSIGNED','UNASSIGNED','UNASSIGNED'],
        'feature_name':['TRDC','TRDC','CD3E','CD3E','NegativeControl','TRDC','TRGC2','TRGC2','CD3G','CD3G','unmeasured_background','TRDC','TRDC'],
        'x_location':pa.array([np.nan if nonfinite_assigned else 0.,30.,100.,-100.,-30.,90.,1000.,1030.,2000.,2000.,-30.,30.0001,0.],type=pa.float32()),
        'y_location':pa.array([0.,30.,100.,-100.,-30.,90.,1000.,1030.,2000.,2000.,30.,0.,0.],type=pa.float32()),
        'z_location':pa.array([1.,2.,3.,4.,5.,6.,7.,8.,9.,10.,11.,12.,-0.],type=pa.float32()),
        'qv':pa.array([20.,40.,30.,30.,21.,19.,20.,40.,30.,30.,8.,40.,None],type=pa.float32()),
        'overlaps_nucleus':[True,False,True,False,True,False,True,False,True,False,False,False,None],
        'native_extra':['original',None,'α','x','x','x','x','x','x','x','x','x','x'],
        'native_integer':pa.array(range(13),type=pa.int16()),
    }).replace_schema_metadata({b'native_fixture':b'preserve-this-exact-metadata'})
    pq.write_table(source,bundle/'transcripts.parquet',row_group_size=4)
    sources={}
    for name in ['cell_feature_matrix.h5',meta_name,'transcripts.parquet']:
        f=bundle/name;sources[name]={'path':str(f),'sha256':sha(f),'bytes':f.stat().st_size,'mtime_ns':f.stat().st_mtime_ns}
    sources['transcripts.parquet']['source_rows']=len(source)
    inventory=root/'inventory.json';write(inventory,{'runs':[{'run_id':'fixture','bundle':str(bundle),'sources':sources}]})
    config=root/'config.json';write(config,{'status':'USER_APPROVED_FOR_EXECUTION','criteria_sha256':'f'*64,'source_inventory_sha256':sha(inventory),'extraction_population':POPULATION,'preserve_flagged':True,'apply_nucleus_filter':False,'apply_tissue_filter':False,'roi_half_width_um':30.,'qv_reconciliation_minimum':20.,'preserve_native_transcript_columns':True})
    cells=pd.DataFrame({'run_id':['fixture','fixture'],'cell_id':['a','b'],'native_matrix_row':[0,1],'native_cell_key':['fixture:a','fixture:b'],'x_centroid':centers,'y_centroid':centers,'foreign_flag_union':[True,False],'nucleus_count':[2,0],'tissue_status':['unknown','unknown'],'TRDC':[2,0],'TRGC2':[0,2],'CD3E':[2,0],'CD3G':[0,2]})
    if empty:cells=cells.iloc[:0]
    cells.to_csv(sel/'eligible_raw_cells.csv',index=False)
    sparse.save_npz(sel/'eligible_native_counts.npz',native[:,cells.native_matrix_row.tolist()].T.tocsr())
    pd.DataFrame({'feature_index':range(len(names)),'feature_name':names}).to_csv(sel/'native_features.csv',index=False)
    manifest={'run_id':'fixture','extraction_population':POPULATION,'config_sha256':sha(config),'criteria_sha256':'f'*64,'source_inventory_sha256':sha(inventory),'n_eligible_cells':len(cells),'n_native_features':len(names),'h5_source':sources['cell_feature_matrix.h5'],'native_metadata_source':sources[meta_name],'output_sha256':{p.name:sha(p) for p in sel.iterdir() if p.is_file()}}
    write(sel/'selector_manifest.json',manifest)
    receipt={'status':'PASS','all_checks_pass':True,'run_id':'fixture','selector_manifest_sha256':sha(sel/'selector_manifest.json'),'eligible_cells_sha256':sha(sel/'eligible_raw_cells.csv'),'config_sha256':sha(config),'criteria_sha256':'f'*64,'source_inventory_sha256':sha(inventory)}
    write(sel/'validation/selection_validation.json',receipt)
    args=SimpleNamespace(config=config,inventory=inventory,selection_dir=sel,selection_validation=None,run_id='fixture',output_dir=root/'extract',batch_size=3)
    return args,source,cells

def v_args(a):
    return SimpleNamespace(**{**vars(a),'extraction_dir':a.output_dir,'output_dir':a.output_dir.parent/'validate','batch_size':5})

def repin_output(a,mutator):
    p=a.output_dir/'transcripts_roi_or_assigned.parquet'
    table=pq.read_table(p);table=mutator(table);pq.write_table(table,p,row_group_size=2)
    mpath=a.output_dir/'extraction_manifest.json';m=json.loads(mpath.read_text())
    m['outputs'][p.name]={'bytes':p.stat().st_size,'sha256':sha(p)}
    m['statistics']['exported_rows']=len(table)
    write(mpath,m)

class Coordinates(unittest.TestCase):
    def setUp(self):
        self.tmp=tempfile.TemporaryDirectory(prefix='coordinate_test_',dir=HERE)
        self.root=Path(self.tmp.name)
        self.quiet=contextlib.redirect_stdout(io.StringIO());self.quiet.__enter__()
    def tearDown(self):
        self.quiet.__exit__(None,None,None);self.tmp.cleanup()
    def test_all_flagged_candidates_full_source_columns_counts_and_closed_edges(self):
        a,source,cells=fixture(self.root);r=extract(a);v=validate(v_args(a))
        self.assertEqual(v['status'],'PASS');self.assertEqual(r['n_foreign_flagged_cells_included'],1)
        self.assertEqual(r['n_eligible_cells'],2);self.assertEqual(r['statistics']['exported_rows'],12)
        self.assertEqual(r['statistics']['assigned_rows_outside_all_rois'],5)
        got=pq.read_table(a.output_dir/'transcripts_roi_or_assigned.parquet')
        self.assertEqual(got['transcript_id'].to_pylist(),[x for x in range(1,14) if x!=12])
        self.assertIn('native_extra',got.column_names)
        self.assertEqual(v['n_cell_feature_pairs_tested_including_zeros'],10)
    def test_empty_selection(self):
        a,_,_=fixture(self.root,empty=True);extract(a);v=validate(v_args(a))
        self.assertEqual(v['n_eligible_cells'],0);self.assertEqual(v['statistics'].get('exported_rows',0),0)
        self.assertEqual(pq.ParquetFile(a.output_dir/'transcripts_roi_or_assigned.parquet').metadata.num_rows,0)
    def test_csv_metadata_pin(self):
        a,_,_=fixture(self.root,metadata_csv=True);extract(a);self.assertEqual(validate(v_args(a))['status'],'PASS')
    def test_source_pin_change_fails(self):
        a,_,_=fixture(self.root)
        with (self.root/'native/transcripts.parquet').open('ab') as f:f.write(b'x')
        with self.assertRaisesRegex(ValueError,'Pinned size'):extract(a)
        self.assertFalse(a.output_dir.exists())
    def test_stale_eligible_table_fails_before_extraction(self):
        a,_,cells=fixture(self.root);cells.iloc[1:].to_csv(a.selection_dir/'eligible_raw_cells.csv',index=False)
        with self.assertRaisesRegex(ValueError,'Eligible table differs'):extract(a)
    def test_native_coordinate_tamper_detected_after_output_repin(self):
        a,_,_=fixture(self.root);extract(a)
        def mutate(t):
            v=t['x_location'].to_pylist();v[0]+=1
            return t.set_column(t.schema.get_field_index('x_location'),t.schema.field('x_location'),pa.array(v,type=pa.float32()))
        repin_output(a,mutate)
        with self.assertRaisesRegex(ValueError,'Native floating values changed'):validate(v_args(a))
    def test_missing_lowq_background_detected_after_output_repin(self):
        a,_,_=fixture(self.root);extract(a)
        repin_output(a,lambda t:t.filter(pa.array([x!=11 for x in t['transcript_id'].to_pylist()])))
        with self.assertRaisesRegex(ValueError,'Native values changed|missing source-selected rows'):validate(v_args(a))
    def test_duplicate_transcript_ids_fail_independent_validation(self):
        a,_,_=fixture(self.root,duplicate=True);extract(a)
        with self.assertRaisesRegex(ValueError,'Duplicate transcript IDs'):validate(v_args(a))
    def test_derived_geometry_flag_tamper_detected(self):
        a,_,_=fixture(self.root);extract(a)
        def mutate(t):
            values=t['in_any_60um_roi'].to_pylist();values[0]=not values[0]
            return t.set_column(t.schema.get_field_index('in_any_60um_roi'),'in_any_60um_roi',pa.array(values))
        repin_output(a,mutate)
        with self.assertRaisesRegex(ValueError,'Derived flag mismatch'):validate(v_args(a))
    def test_no_overwrite(self):
        a,_,_=fixture(self.root);extract(a)
        with self.assertRaises(FileExistsError):extract(a)
        validate(v_args(a))
        with self.assertRaises(FileExistsError):validate(v_args(a))
    def test_independent_roi_nonfinite_and_closed_corner(self):
        xy=np.array([[30.,30.],[-30.,-30.],[30.000001,0.],[np.nan,0.],[np.inf,0.],[100.,100.]])
        self.assertEqual(independent_roi(xy,np.array([[0.,0.]])).tolist(),[True,True,False,False,False,False])
    def test_nonfinite_assigned_native_coordinate_retained(self):
        a,_,_=fixture(self.root,nonfinite_assigned=True);extract(a)
        self.assertEqual(validate(v_args(a))['status'],'PASS')
        self.assertTrue(np.isnan(pq.read_table(a.output_dir/'transcripts_roi_or_assigned.parquet')['x_location'][0].as_py()))
    def test_independent_roi_randomized_against_explicit_squares(self):
        rng=np.random.default_rng(42);centers=rng.normal(size=(30,2))*100
        xy=np.concatenate([rng.normal(size=(4000,2))*100,centers+30.,centers-30.])
        direct=np.zeros(len(xy),dtype=bool)
        for center in centers:direct|=(np.abs(xy-center)<=30.).all(axis=1)
        self.assertTrue(np.array_equal(independent_roi(xy,centers),direct))
    def test_native_double_centroids_survive_eligible_csv_exactly(self):
        a,_,_=fixture(self.root,fractional_centroids=True)
        source=pq.read_table(self.root/'native/cells.parquet').select(['x_centroid','y_centroid']).to_pandas().to_numpy()[:2]
        self.assertTrue(np.array_equal(load_contract(a)['xy'],source))
        lower_precision=pd.read_csv(a.selection_dir/'eligible_raw_cells.csv', dtype={'x_centroid': 'float32', 'y_centroid': 'float32'})[['x_centroid','y_centroid']].to_numpy()
        self.assertFalse(np.array_equal(lower_precision,source),'Fixture must distinguish native doubles from rounded coordinates')
        extract(a);self.assertEqual(validate(v_args(a))['status'],'PASS')

if __name__=='__main__':
    unittest.main(verbosity=2)
