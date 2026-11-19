#!/usr/bin/env python3
"""Synthetic-only selection tests. Temporary fixtures are removed after each test."""
import argparse, hashlib, itertools, json, tempfile, unittest
from pathlib import Path
import h5py
import numpy as np
import pandas as pd
from scipy import sparse
from multiplex5000.gammadelta import select_corrected as selector
from multiplex5000.gammadelta import validate_selection as validator

CRITERIA=Path(selector.__file__).with_name('reviewed_criteria.json')
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def write(p,o):Path(p).write_text(json.dumps(o,indent=2)+'\n')
def terms(node):
    if node['op']=='count':return [{node['gene']}]
    children=[terms(x) for x in node['args']]
    if node['op']=='any':return [s for c in children for s in c]
    return [set.union(*x) for x in itertools.product(*children)]

class SelectionTests(unittest.TestCase):
    def setUp(self):
        self.tmp=tempfile.TemporaryDirectory(prefix='gd-selection-synthetic-');self.root=Path(self.tmp.name);self.criteria=json.loads(CRITERIA.read_text());self.region='0000001'
    def tearDown(self):self.tmp.cleanup()
    def fixture(self,empty=False,parquet=False,missing=False,antibody=False,duplicate=False):
        names=list(reversed(validator.GENES))+['NEG_CONTROL'];nf=len(names);n=25
        x=np.zeros((n,nf),dtype=np.int64);cols={g:i for i,g in enumerate(names)}
        def setrow(i,gs):
            for g in gs:x[i,cols[g]]=2
        if not empty:
            setrow(0,['CD3E','TRDC']);setrow(1,['CD3G','TRGC2']);setrow(2,['TRDC','TRGC2']);x[2,cols['CD3E']]=x[2,cols['CD3G']]=1
            setrow(3,['CD3G']);setrow(4,['TRGC2'])
            setrow(5,['CD3E','TRDC']+list(terms(self.criteria['foreign_programs'][0]['predicate'])[0]|terms(self.criteria['foreign_programs'][1]['predicate'])[0]))
            for i,r in enumerate(self.criteria['foreign_programs']):setrow(6+i,['CD3E','TRDC']+list(terms(r['predicate'])[0]))
            setrow(21,['CD3E','TRDC','SDC1']);setrow(22,['CD3E','TRDC','TRAC','TRBC2'])
        x[0,cols['NEG_CONTROL']]=7
        b=self.root/'output-SYNTH__0000001__Region_1__20260101__000000';b.mkdir()
        ids=[f'cell{i:03}' for i in range(n)];types=['Gene Expression']*(nf-1)+['Negative Control Probe']
        if missing:names[cols['CD3E']]='UNKNOWN_GENE'
        if antibody:types[cols['CD3E']]='Antibody Capture'
        if duplicate:names[-1]='CD3E';types[-1]='Gene Expression'
        h=b/'cell_feature_matrix.h5';mat=sparse.csc_matrix(x.T)
        with h5py.File(h,'w') as f:
            m=f.create_group('matrix');m['shape']=np.array([nf,n]);m['data']=mat.data;m['indices']=mat.indices;m['indptr']=mat.indptr;m['barcodes']=np.asarray(ids,dtype='S')
            g=m.create_group('features');g['id']=np.asarray([f'G{i}' for i in range(nf)],dtype='S');g['name']=np.asarray(names,dtype='S');g['feature_type']=np.asarray(types,dtype='S')
        rng=np.random.default_rng(5101)
        meta=pd.DataFrame({'cell_id':ids,'x_centroid':rng.uniform(0,20000,n),'y_centroid':rng.uniform(0,20000,n),'nucleus_count':[0,2]+[1]*(n-2),'cell_area':rng.uniform(0,100,n),'nucleus_area':rng.uniform(0,100,n),'total_counts':x.sum(axis=1),'transcript_counts':x[:,:-1].sum(axis=1),'segmentation_method':'synthetic'})
        mp=b/('cells.parquet' if parquet else 'cells.csv.gz')
        if parquet:meta.to_parquet(mp,index=False)
        else:meta.to_csv(mp,index=False)
        pins=lambda p:dict(path=str(p),sha256=sha(p),bytes=p.stat().st_size,mtime_ns=p.stat().st_mtime_ns)
        inv={'runs':[dict(run_id=self.region,bundle=str(b),sources={'cell_feature_matrix.h5':pins(h),mp.name:pins(mp)},native_unit={'region_name':'Region_1'})]};ip=self.root/'inventory.json';write(ip,inv)
        actual_genes=list(reversed(validator.GENES))+['NEG_CONTROL'];c={g:i for i,g in enumerate(actual_genes)}
        t=(x[:,c['CD3E']]>=2)|(x[:,c['CD3G']]>=2);gamma=(x[:,c['TRDC']]>=2)|(x[:,c['TRGC2']]>=2);base=t&gamma
        old=(x[:,c['TRDC']]>=2)&(x[:,c['TRGC2']]>=2)&((x[:,c['CD3E']]+x[:,c['CD3G']])>=2)
        andp=self.root/'old_and.csv';orp=self.root/'old_or2.csv'
        for p,mask in [(andp,old),(orp,base)]:pd.DataFrame({'run_id':[self.region]*int(mask.sum()),'cell_id':np.array(ids)[mask]}).to_csv(p,index=False)
        cp=self.root/'config.json';write(cp,dict(status='USER_APPROVED_FOR_EXECUTION',criteria_sha256=sha(CRITERIA),source_inventory_sha256=sha(ip),historical_and_sha256_by_region={self.region:sha(andp)},historical_or2_sha256=sha(orp),preserve_flagged=True,apply_nucleus_filter=False,apply_tissue_filter=False,extraction_population='all_base_candidates_including_flagged',per_gene_minimum=2))
        self.args=argparse.Namespace(inventory=ip,config=cp,criteria=CRITERIA,region=self.region,out_dir=self.root/'selection',historical_and_cells=andp,historical_or2_cells=orp,chunk_cells=7)
        self.expected_t=t;self.expected_base=base;return self.args
    def v_args(self):
        d=vars(self.args).copy();d['selection_dir']=self.args.out_dir;d['out_dir']=self.args.out_dir/'validation';d['chunk_cells']=3;return argparse.Namespace(**d)
    def test_independent_rules_against_ast_and_boundaries(self):
        col={g:i for i,g in enumerate(validator.GENES)};vectors=[]
        for rule in self.criteria['foreign_programs']:
            for term in terms(rule['predicate']):
                a=np.zeros(67,dtype=int)
                for g in term:a[col[g]]=2
                vectors.append(a.copy())
                for g in term:
                    b=a.copy();b[col[g]]=1;vectors.append(b)
        for core in itertools.product(range(3),repeat=4):
            a=np.zeros(67,dtype=int);a[:4]=core;vectors.append(a)
        vectors.extend(np.random.default_rng(5000).integers(0,4,size=(1000,67)))
        x=np.asarray(vectors);t,g,f,ab=validator.independent_masks(x)
        expected=np.column_stack([selector.evaluate(r['predicate'],x,col) for r in self.criteria['foreign_programs']])
        np.testing.assert_array_equal(f,expected)
        np.testing.assert_array_equal(t&g,selector.evaluate(self.criteria['base_predicate'],x,col))
        self.assertEqual(sha(CRITERIA),selector.REVIEW_SHA)
    def test_complete_preserved_and_separate_nuclei(self):
        a=self.fixture();m=selector.run(a);v=validator.validate(self.v_args());self.assertTrue(v['all_checks_pass'])
        d=pd.read_csv(a.out_dir/'eligible_raw_cells.csv');self.assertEqual(len(d),int(self.expected_base.sum()));self.assertTrue({'cell000','cell001'}<=set(d.cell_id))
        self.assertEqual(m['n_t_gate_cells'],int(self.expected_t.sum()));self.assertTrue(d.foreign_flag_union.any());self.assertTrue((~d.foreign_flag_union).any())
        h=json.loads((a.out_dir/'historical_membership_summary.json').read_text());self.assertTrue(h['historical_or2']['exact_membership_match']);self.assertEqual(h['historical_and']['removed'],1)
    def test_empty_t_gate_outputs(self):
        a=self.fixture(empty=True);selector.run(a);self.assertEqual(validator.validate(self.v_args())['gate_counts']['base_candidates'],0)
    def test_native_parquet_metadata(self):
        a=self.fixture(parquet=True);selector.run(a);self.assertEqual(validator.validate(self.v_args())['status'],'PASS')
    def test_native_float32_metadata_value_preserved(self):
        a=self.fixture(parquet=True);i=json.loads(a.inventory.read_text());pin=i['runs'][0]['sources']['cells.parquet'];p=Path(pin['path']);d=pd.read_parquet(p)
        for c in ['x_centroid','y_centroid','cell_area','nucleus_area']:d[c]=d[c].astype(np.float32)
        d.to_parquet(p,index=False);pin.update(sha256=sha(p),bytes=p.stat().st_size,mtime_ns=p.stat().st_mtime_ns);write(a.inventory,i);c=json.loads(a.config.read_text());c['source_inventory_sha256']=sha(a.inventory);write(a.config,c)
        selector.run(a);self.assertEqual(validator.validate(self.v_args())['status'],'PASS')
    def test_missing_required_gene(self):
        with self.assertRaisesRegex(ValueError,'Missing/duplicate/non-RNA'):selector.run(self.fixture(missing=True))
    def test_antibody_feature_not_rna(self):
        with self.assertRaisesRegex(ValueError,'Missing/duplicate/non-RNA'):selector.run(self.fixture(antibody=True))
    def test_duplicate_required_gene(self):
        with self.assertRaisesRegex(ValueError,'Missing/duplicate/non-RNA'):selector.run(self.fixture(duplicate=True))
    def test_wrong_native_hash(self):
        a=self.fixture();i=json.loads(a.inventory.read_text());i['runs'][0]['sources']['cell_feature_matrix.h5']['sha256']='0'*64;write(a.inventory,i);c=json.loads(a.config.read_text());c['source_inventory_sha256']=sha(a.inventory);write(a.config,c)
        with self.assertRaisesRegex(ValueError,'Native source'):selector.run(a)
    def test_wrong_native_mtime(self):
        a=self.fixture();i=json.loads(a.inventory.read_text());i['runs'][0]['sources']['cell_feature_matrix.h5']['mtime_ns']-=1;write(a.inventory,i);c=json.loads(a.config.read_text());c['source_inventory_sha256']=sha(a.inventory);write(a.config,c)
        with self.assertRaisesRegex(ValueError,'Native source'):selector.run(a)
    def test_unapproved_or_filtering_config(self):
        a=self.fixture();c=json.loads(a.config.read_text());c['apply_nucleus_filter']=True;write(a.config,c)
        with self.assertRaisesRegex(ValueError,'Preserve-all'):selector.run(a)
        c['apply_nucleus_filter']=False;c['status']='PENDING';write(a.config,c)
        with self.assertRaisesRegex(ValueError,'approval'):selector.run(a)
    def test_stale_output_refused(self):
        a=self.fixture();a.out_dir.mkdir();(a.out_dir/'existing.txt').write_text('retain')
        with self.assertRaises(FileExistsError):selector.run(a)
    def test_self_consistently_rehashed_mask_corruption_rejected(self):
        a=self.fixture();selector.run(a);p=a.out_dir/'whole_raw_predicate_counts.npz'
        with np.load(p,allow_pickle=False) as z:values={k:z[k] for k in z.files}
        values['rule_flags'][0,0]=not values['rule_flags'][0,0];np.savez_compressed(p,**values)
        m=json.loads((a.out_dir/'selector_manifest.json').read_text());m['output_sha256'][p.name]=sha(p);write(a.out_dir/'selector_manifest.json',m)
        with self.assertRaisesRegex(AssertionError,'Independent mask'):validator.validate(self.v_args())
    def test_manifest_eligible_count_corruption_rejected(self):
        a=self.fixture();selector.run(a);p=a.out_dir/'selector_manifest.json';m=json.loads(p.read_text());m['n_eligible_cells']+=1;write(p,m)
        with self.assertRaisesRegex(AssertionError,'Manifest aggregate n_eligible_cells'):validator.validate(self.v_args())

if __name__=='__main__':unittest.main(verbosity=2)
