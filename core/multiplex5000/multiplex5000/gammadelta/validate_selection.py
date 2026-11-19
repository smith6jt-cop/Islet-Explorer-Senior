#!/usr/bin/env python3
"""Independent native-H5 validation; never import the selector/evaluate its AST."""
import argparse, csv, hashlib, json, os, sys
from datetime import datetime, timezone
from pathlib import Path
import h5py
import numpy as np
import pandas as pd
import pyarrow.parquet as pq
from scipy import sparse

REVIEW_SHA='b060488f17175d346277f112683c3e48805a663d3663ad16cbc81ce9d56eb1c3'
CORE=['TRDC','TRGC2','CD3E','CD3G']
RULE_IDS=['B_program','myeloid_program','ductal_epithelial_program','acinar_program','endocrine_program','endothelial_program','plasma_program','mural_stromal_program','fibroblast_program','mast_basophil_program','neutrophil_program','neural_program','cDC1_program','cDC2_program','pDC_program']
FOREIGN='CD79A CD79B MS4A1 CD19 CSF1R SPI1 CD14 CD163 FCGR1A EPCAM KRT19 CFTR MUC1 SOX9 CUZD1 CTRL AQP8 PDIA2 INSM1 NEUROD1 ISL1 CHGA SCG2 SCG5 PCSK1 CDH5 CLDN5 KDR PLVAP ADGRL4 MZB1 SDC1 TNFRSF17 PDGFRB RGS5 CSPG4 ACTA2 PDGFRA DPT COL5A1 COL5A2 FAP MS4A2 HDC KIT GATA2 CSF3R FCGR3B CXCR2 SOX10 MPZ PMP22 CLEC9A XCR1 BATF3 CD1C CLEC10A FCER1A CLEC4C TCF4 IL3RA'.split()
GENES=CORE+sorted(set(FOREIGN+['TRAC','TRBC2'])-set(CORE))

def sha(p):
    h=hashlib.sha256()
    with Path(p).open('rb') as f:
        for b in iter(lambda:f.read(8<<20),b''):h.update(b)
    return h.hexdigest()
def stat(p):
    s=Path(p).stat();return (s.st_size,s.st_mtime_ns,s.st_ino,s.st_dev)
def decode(a):return np.array([x.decode() if isinstance(x,bytes) else str(x) for x in a],dtype=object)
def require(cond,message):
    if not cond:raise AssertionError(message)
def independent_masks(x):
    # Exact fifteen reviewed definitions, authored independently from the JSON evaluator.
    p={g:x[:,i]>=2 for i,g in enumerate(GENES)}
    f=[
        p['CD79A']&p['CD79B']&(p['MS4A1']|p['CD19']),
        p['CSF1R']&p['SPI1']&(p['CD14']|p['CD163']|p['FCGR1A']),
        p['EPCAM']&p['KRT19']&(p['CFTR']|p['MUC1']|p['SOX9']),
        p['CUZD1']&p['CTRL']&(p['AQP8']|p['PDIA2']),
        (p['INSM1']|p['NEUROD1']|p['ISL1'])&(sum(p[g].astype(np.int8) for g in ['CHGA','SCG2','SCG5','PCSK1'])>=2),
        p['CDH5']&p['CLDN5']&(p['KDR']|p['PLVAP']|p['ADGRL4']),
        p['MZB1']&p['SDC1']&p['TNFRSF17'],
        p['PDGFRB']&p['RGS5']&(p['CSPG4']|p['ACTA2']),
        p['PDGFRA']&p['DPT']&(p['COL5A1']|p['COL5A2']|p['FAP']),
        p['MS4A2']&p['HDC']&(p['KIT']|p['GATA2']),
        p['CSF3R']&p['FCGR3B']&p['CXCR2'],
        p['SOX10']&p['MPZ']&p['PMP22'],
        p['CLEC9A']&p['XCR1']&p['BATF3'],
        p['CD1C']&p['CLEC10A']&p['FCER1A'],
        p['CLEC4C']&p['TCF4']&p['IL3RA']]
    t=p['CD3E']|p['CD3G'];g=p['TRDC']|p['TRGC2'];flags=np.column_stack(f)
    return t,g,flags,p['TRAC']&p['TRBC2']
def bools(series):
    values=[]
    for x in series:
        if x is True or isinstance(x,np.bool_) and x:values.append(True)
        elif x is False or isinstance(x,np.bool_) and not x:values.append(False)
        elif isinstance(x,str) and x in {'True','False'}:values.append(x=='True')
        else:raise AssertionError('Non-Boolean value in Boolean output column')
    return np.asarray(values,dtype=bool)
def get_sources(inv,region):
    rows=[r for r in inv.get('runs',[]) if r['run_id']==region]
    require(len(rows)==1,'Canonical inventory region uniqueness');r=rows[0]
    meta=[x for x in r['sources'] if x in {'cells.csv.gz','cells.csv','cells.parquet'}]
    require(len(meta)==1,'Exactly one metadata source');return r,r['sources']['cell_feature_matrix.h5'],r['sources'][meta[0]]
def hist(path,region,expected):
    require(sha(path)==expected,'Historical hash')
    d=pd.read_csv(path,dtype={'cell_id':str,'run_id':str},float_precision='round_trip')
    if 'run_id' in d:d=d[d.run_id==region]
    require(not d.cell_id.duplicated().any(),'Historical native IDs duplicated');return set(d.cell_id)
def validate(a):
    code_sha=sha(__file__);config=json.loads(Path(a.config).read_text());criteria_sha=sha(a.criteria);config_sha=sha(a.config);inventory_sha=sha(a.inventory)
    require(criteria_sha==REVIEW_SHA==config['criteria_sha256'],'Reviewed criteria digest')
    require(config['status']=='USER_APPROVED_FOR_EXECUTION' and config['preserve_flagged'] is True and config['apply_nucleus_filter'] is False and config['apply_tissue_filter'] is False and config['extraction_population']=='all_base_candidates_including_flagged','Approval/preservation contract')
    require(config['source_inventory_sha256']==inventory_sha,'Inventory bound by approval')
    require(config.get('per_gene_minimum',2)==2,'Approved per-gene minimum')
    criteria=json.loads(Path(a.criteria).read_text());require([r['id'] for r in criteria['foreign_programs']]==RULE_IDS,'Exact rule ordering')
    pin,hp,mp=get_sources(json.loads(Path(a.inventory).read_text()),a.region);h5=Path(hp['path']);meta_path=Path(mp['path'])
    require(h5.parent==meta_path.parent==Path(pin['bundle']),'Exact pinned source bundle')
    hb,mb=stat(h5),stat(meta_path)
    require(hb[0]==hp['bytes'] and hb[1]==hp['mtime_ns'] and sha(h5)==hp['sha256'],'Full native H5 hash/size/mtime')
    require(mb[0]==mp['bytes'] and mb[1]==mp['mtime_ns'] and sha(meta_path)==mp['sha256'],'Full native metadata hash/size/mtime')
    selection=Path(a.selection_dir);sm_path=selection/'selector_manifest.json';sm_hash=sha(sm_path);sm=json.loads(sm_path.read_text())
    require(sm['config_sha256']==config_sha and sm['criteria_sha256']==criteria_sha and sm['source_inventory_sha256']==inventory_sha,'Selector pins')
    require(sm['run_id']==a.region and sm['bundle']==pin['bundle'] and sm['extraction_population']=='all_base_candidates_including_flagged','Selector population identity')
    require(sm['nucleus_filter_applied'] is False and sm['tissue_filter_applied'] is False,'No nucleus/tissue filter')
    for name,digest in sm['output_sha256'].items():require(sha(selection/name)==digest,'Selector output hash: '+name)
    require(sm['required_genes']==GENES and len(GENES)==67,'Exact independently declared gene order')
    with np.load(selection/'whole_raw_predicate_counts.npz',allow_pickle=False) as archive:
        cache={key:archive[key] for key in archive.files}  # Decompress each array exactly once.
    require(cache['genes'].tolist()==GENES and cache['rule_ids'].tolist()==RULE_IDS,'Cache axes')
    ttable=pd.read_csv(selection/'t_gate_cells.csv.gz',dtype={'run_id':str,'cell_id':str,'source_slide_id':str},float_precision='round_trip')
    table=pd.read_csv(selection/'eligible_raw_cells.csv',dtype={'run_id':str,'cell_id':str,'source_slide_id':str},float_precision='round_trip')
    tm=sparse.load_npz(selection/'t_gate_native_counts.npz').tocsr();em=sparse.load_npz(selection/'eligible_native_counts.npz').tocsr()
    t_all=[];g_all=[];f_all=[];a_all=[];t_offset=0
    with h5py.File(h5,'r') as f:
        m=f['matrix'];nf,nc=map(int,m['shape'][:]);ids=decode(m['barcodes'][:]);ptr=m['indptr'][:]
        require(len(ids)==nc and len(set(ids))==nc and len(ptr)==nc+1 and ptr[0]==0 and np.all(np.diff(ptr)>=0) and ptr[-1]==len(m['data']),'Native CSC/ID structure')
        require(np.array_equal(cache['cell_ids'],ids),'Cache all-native barcode order')
        names=decode(m['features/name'][:]);types=decode(m['features/feature_type'][:]);lookup=np.full(nf,-1,dtype=np.int64)
        feature_csv=pd.read_csv(selection/'native_features.csv',keep_default_na=False)
        require(np.array_equal(feature_csv.feature_index.to_numpy(),np.arange(nf)),'Feature indices')
        for name,node in m['features'].items():
            if isinstance(node,h5py.Dataset) and node.ndim==1 and len(node)==nf:require(np.array_equal(feature_csv['feature_'+name].to_numpy(str),decode(node[:])),'Native feature column '+name)
        for j,gene in enumerate(GENES):
            hits=np.flatnonzero(names==gene);require(len(hits)==1 and types[hits[0]]=='Gene Expression','Required gene absent/duplicate/antibody-type: '+gene);lookup[hits[0]]=j
        require(tm.shape[1]==em.shape[1]==nf,'All native feature columns preserved')
        for start in range(0,nc,a.chunk_cells):
            stop=min(nc,start+a.chunk_cells);lo,hi=int(ptr[start]),int(ptr[stop]);v=m['data'][lo:hi];ix=m['indices'][lo:hi]
            require(np.isfinite(v).all() and np.all(v>=0) and np.all(v==np.floor(v)) and np.all(ix>=0) and np.all(ix<nf),'Native finite integer counts')
            # Independent CSC interpretation: entry-to-cell repeats + indexed accumulation.
            cell=np.repeat(np.arange(stop-start),np.diff(ptr[start:stop+1]));wanted=lookup[ix]>=0
            x=np.zeros((stop-start,len(GENES)),dtype=np.int64);np.add.at(x,(cell[wanted],lookup[ix[wanted]]),v[wanted].astype(np.int64))
            require(np.array_equal(x,cache['counts'][start:stop]),'All-native 67-gene counts differ')
            t,g,flags,ab=independent_masks(x);t_all.append(t);g_all.append(g);f_all.append(flags);a_all.append(ab)
            full=sparse.coo_matrix((v.astype(np.int64),(cell,ix)),shape=(stop-start,nf)).tocsr()[t]
            require((full!=tm[t_offset:t_offset+int(t.sum())]).nnz==0,'Preserved full T-gate count matrix differs')
            t_offset+=int(t.sum())
    t=np.concatenate(t_all);g=np.concatenate(g_all);flags=np.concatenate(f_all);ab=np.concatenate(a_all);base=t&g;union=flags.any(axis=1)
    for key,value in {'n_native_cells':nc,'n_native_features':nf,'n_t_gate_cells':int(t.sum()),'n_eligible_cells':int(base.sum()),'n_base_flagged':int((base&union).sum()),'n_base_unflagged':int((base&~union).sum())}.items():require(sm[key]==value,'Manifest aggregate '+key)
    require(sm['rule_ids']==RULE_IDS,'Manifest rule identities')
    for key,expected in [('h5_source',hp),('native_metadata_source',mp)]:
        require(all(sm[key][k]==expected[k] for k in ['path','sha256','bytes','mtime_ns']),'Manifest native source binding '+key)
    for key,value in {'t_gate':t,'receptor_or':g,'rule_flags':flags,'base_candidate':base,'foreign_flag_union':union,'alpha_beta_flag':ab}.items():require(np.array_equal(cache[key],value),'Independent mask differs: '+key)
    require(t_offset==len(ttable)==tm.shape[0],'All T-gate rows preserved')
    for label,d,mask in [('T',ttable,t),('base',table,base)]:
        ix=np.flatnonzero(mask);require(d.cell_id.tolist()==ids[mask].tolist() and np.array_equal(d.native_matrix_row,ix),'Exact membership/order '+label)
        require(d.run_id.eq(a.region).all() and d.source_bundle.eq(pin['bundle']).all() and d.native_cell_key.tolist()==[pin['bundle']+'|'+s for s in ids[mask]],'Exact source keys '+label)
        require(np.array_equal(d[GENES].to_numpy(np.int64),cache['counts'][mask]),'Stored marker counts '+label)
        expected={'t_gate':t[mask],'receptor_or':g[mask],'base_candidate':base[mask],'foreign_flag_union':union[mask],'base_flagged':(base&union)[mask],'base_unflagged':(base&~union)[mask],'alpha_beta_flag':ab[mask]}
        expected.update({'flag__'+name:flags[mask,j] for j,name in enumerate(RULE_IDS)})
        for key,val in expected.items():require(np.array_equal(bools(d[key]),val),'Per-cell flag '+label+'/'+key)
        require(not bools(d.nucleus_used_for_selection).any() and not bools(d.tissue_used_for_selection).any() and d.tissue_eligibility.eq('not_evaluated').all(),'Separate technical/tissue eligibility')
    require((em!=tm[g[t]]).nnz==0 and em.shape[0]==len(table),'Base all-feature matrix is exact T subset')
    pd.testing.assert_frame_equal(ttable.loc[g[t]].reset_index(drop=True),table,check_dtype=False,check_exact=True)
    require(np.array_equal(np.asarray(tm.sum(axis=1)).ravel(),ttable.native_all_feature_count.to_numpy(np.int64)),'Native all-feature totals')
    # Metadata is independently streamed and matched by native barcode, with no eligibility filtering.
    row_index={s:i for i,s in enumerate(ttable.cell_id)};seen=set();meta_columns=[]
    batches=(b.to_pandas() for b in pq.ParquetFile(meta_path).iter_batches(batch_size=77777)) if meta_path.suffix=='.parquet' else pd.read_csv(meta_path,dtype={'cell_id':str},chunksize=77777,float_precision='round_trip')
    for batch in batches:
        meta_columns=[c for c in batch if c in ttable and c not in GENES]
        for row in batch[batch.cell_id.isin(row_index)].itertuples(index=False):
            r=row._asdict();key=str(r['cell_id']);require(key not in seen,'Duplicate metadata selected barcode');seen.add(key);stored=ttable.iloc[row_index[key]]
            for k in meta_columns:
                v=r[k];w=stored[k]
                require((pd.isna(v) and pd.isna(w)) or (str(v)==str(w)) or (isinstance(v,(int,float,np.number)) and isinstance(w,(int,float,np.number)) and v==w),'Native metadata mismatch '+key+'/'+k)
    require(seen==set(row_index),'Native metadata membership complete')
    expected_nucleus=np.select([ttable.nucleus_count.eq(0),ttable.nucleus_count.eq(1),ttable.nucleus_count.gt(1)],['zero','single','multiple'],default='unknown')
    require(np.array_equal(ttable.nucleus_status.to_numpy(),expected_nucleus),'Nucleus status is descriptive metadata')
    for name,mask in [('base_unflagged.csv',base&~union),('base_flagged.csv',base&union)]:
        d=pd.read_csv(selection/name,dtype={'run_id':str,'cell_id':str,'source_slide_id':str},float_precision='round_trip');require(d.cell_id.tolist()==ids[mask].tolist(),'Flagged/unflagged view exact membership')
        pd.testing.assert_frame_equal(table.loc[bools(table.base_flagged) if name=='base_flagged.csv' else bools(table.base_unflagged)].reset_index(drop=True),d,check_dtype=False,check_exact=True)
    audit=json.loads((selection/'gate_counts.json').read_text())
    totals={'all_raw_cells':nc,'t_gate':int(t.sum()),'t_gate_foreign_flagged':int((t&union).sum()),'t_gate_unflagged':int((t&~union).sum()),'receptor_or_all_raw':int(g.sum()),'base_candidates':int(base.sum()),'base_flagged':int((base&union).sum()),'base_unflagged':int((base&~union).sum()),'base_alpha_beta_flag':int((base&ab).sum()),'extraction_eligible':int(base.sum())}
    for key,val in totals.items():require(audit[key]==val,'Gate aggregate '+key)
    prior=np.zeros(nc,dtype=bool)
    for i,r in enumerate(audit['named_rules']):
        f=flags[:,i];expected={'all_raw':int(f.sum()),'t_gate':int((f&t).sum()),'base':int((f&base).sum()),'ordered_newly_flagged_t':int((f&t&~prior).sum()),'ordered_newly_flagged_base':int((f&base&~prior).sum())};prior|=f
        expected.update(ordered_t_unflagged_remaining=int((t&~prior).sum()),ordered_base_unflagged_remaining=int((base&~prior).sum()))
        require(r['id']==RULE_IDS[i] and all(r[k]==v for k,v in expected.items()),'Rule accounting '+RULE_IDS[i])
    require(len(audit['named_rules'])==15 and len(audit['pairwise_overlaps'])==105,'Complete rule/overlap accounting')
    for r in audit['pairwise_overlaps']:
        both=flags[:,RULE_IDS.index(r['rule_a'])]&flags[:,RULE_IDS.index(r['rule_b'])]
        require(r['all_raw']==int(both.sum()) and r['t_gate']==int((both&t).sum()) and r['base']==int((both&base).sum()),'Pairwise overlap')
    patterns=pd.read_csv(selection/'flag_overlap_patterns.csv',keep_default_na=False)
    bitmask=(flags.astype(np.uint32)*(1<<np.arange(15,dtype=np.uint32))).sum(axis=1)
    for label,mask in [('all_raw',np.ones(nc,dtype=bool)),('t_gate',t),('base',base)]:
        u,n=np.unique(bitmask[mask],return_counts=True);saved=patterns[patterns.universe==label]
        require(saved.bitmask.astype(int).tolist()==u.astype(int).tolist() and saved.cells.astype(int).tolist()==n.astype(int).tolist(),'Exact multi-rule overlap patterns')
        require(saved.rule_ids.tolist()==[';'.join(RULE_IDS[j] for j in range(15) if int(b)&(1<<j)) for b in u],'Overlap program identities')
    old_and=hist(a.historical_and_cells,a.region,config['historical_and_sha256_by_region'][a.region]);old_or=hist(a.historical_or2_cells,a.region,config['historical_or2_sha256']);new=set(ids[base]);summary=json.loads((selection/'historical_membership_summary.json').read_text())
    histrows=pd.read_csv(selection/'historical_membership_changes.csv',dtype={'cell_id':str,'run_id':str},float_precision='round_trip');expected_keys=set(old_and)|set(old_or)|new
    require(set(histrows.cell_id)==expected_keys and not histrows.cell_id.duplicated().any(),'Historical comparison union')
    for name,old in [('historical_and',old_and),('historical_or2',old_or)]:
        expected=dict(old_cells=len(old),new_base_cells=len(new),shared=len(old&new),added=len(new-old),removed=len(old-new),exact_membership_match=old==new);require(summary[name]==expected,'Historical summary '+name)
        require(np.array_equal(bools(histrows[name]),histrows.cell_id.isin(old).to_numpy()),'Historical membership '+name)
    require(np.array_equal(bools(histrows.new_base),histrows.cell_id.isin(new).to_numpy()),'New membership changes')
    lookup={c:i for i,c in enumerate(ids)};histix=np.array([lookup[c] for c in histrows.cell_id],dtype=int)
    require(np.array_equal(histrows.native_matrix_row.to_numpy(),histix) and histrows.native_cell_key.tolist()==[pin['bundle']+'|'+c for c in histrows.cell_id],'Historical exact native row/keys')
    require(np.array_equal(histrows[CORE].to_numpy(np.int64),cache['counts'][histix,:4]),'Historical four-marker evidence')
    require(np.array_equal(bools(histrows.new_base_flagged),(base&union)[histix]) and np.array_equal(bools(histrows.new_base_unflagged),(base&~union)[histix]),'Historical flagged/unflagged membership')
    require(stat(h5)==hb and stat(meta_path)==mb,'Sources stable during validation')
    require(sha(sm_path)==sm_hash and sha(a.config)==config_sha and sha(a.criteria)==criteria_sha and sha(a.inventory)==inventory_sha and sha(__file__)==code_sha,'Witnesses/code stable during validation')
    out=Path(a.out_dir)
    if out.exists() and any(out.iterdir()):raise FileExistsError('Refusing nonempty validation output')
    out.mkdir(parents=True,exist_ok=True)
    result=dict(schema='gammadelta-selection-validation/1',status='PASS',all_checks_pass=True,created_utc=datetime.now(timezone.utc).isoformat(),run_id=a.region,selector_manifest_sha256=sm_hash,eligible_cells_sha256=sha(selection/'eligible_raw_cells.csv'),config_sha256=config_sha,criteria_sha256=criteria_sha,source_inventory_sha256=inventory_sha,script_sha256=code_sha,independence='No selector import or predicate evaluator; hardcoded fifteen approved rules; native CSC entry-to-cell accumulation and exact all-raw count/mask comparison.',extraction_population='all_base_candidates_including_flagged',gate_counts=totals,historical_membership=summary,nucleus_filter_applied=False,tissue_filter_applied=False,biological_identity_accepted=False)
    (out/'selection_validation.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps({'status':'PASS','run_id':a.region,**totals}));return result

def parser():
    p=argparse.ArgumentParser(description=__doc__)
    for n in ['inventory','config','criteria','selection-dir','historical-and-cells','historical-or2-cells']:p.add_argument('--'+n,required=True,type=Path)
    p.add_argument('--region','--run-id',dest='region',required=True);p.add_argument('--out-dir','--output-dir',dest='out_dir',required=True,type=Path);p.add_argument('--chunk-cells',type=int,default=7777);return p
if __name__=='__main__':
    a=parser().parse_args()
    if a.chunk_cells<1:raise ValueError('chunk-cells must be positive')
    validate(a)
