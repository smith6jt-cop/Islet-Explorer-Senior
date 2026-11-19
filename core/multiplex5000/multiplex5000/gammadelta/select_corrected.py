#!/usr/bin/env python3
"""Pinned native RNA selection; preserve every T-gate and base gamma candidate."""
import argparse, hashlib, json, os, platform, socket, sys
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
META=['cell_id','x_centroid','y_centroid','transcript_counts','control_probe_counts','genomic_control_counts','control_codeword_counts','unassigned_codeword_counts','deprecated_codeword_counts','total_counts','cell_area','nucleus_area','nucleus_count','segmentation_method']

def sha(p):
    h=hashlib.sha256()
    with Path(p).open('rb') as f:
        for b in iter(lambda:f.read(8<<20),b''):h.update(b)
    return h.hexdigest()
def stat(p):
    s=Path(p).stat();return dict(size=s.st_size,mtime_ns=s.st_mtime_ns,inode=s.st_ino,device=s.st_dev)
def strings(a):return np.asarray([x.decode() if isinstance(x,bytes) else str(x) for x in a],dtype=object)
def dump(p,o):Path(p).write_text(json.dumps(o,indent=2)+'\n')
def node_genes(n):
    if n.get('op')=='count':
        if set(n)!={'op','gene','ge'} or type(n['ge']) is not int or n['ge']!=2:raise ValueError('Only individual Gene Expression count >=2 is approved')
        return {n['gene']}
    if n.get('op') not in {'all','any'} or set(n)!={'op','args'} or not n['args']:raise ValueError('Unsupported predicate')
    return set().union(*(node_genes(c) for c in n['args']))
def evaluate(n,x,col):
    if n['op']=='count':return x[:,col[n['gene']]]>=2
    vals=[evaluate(c,x,col) for c in n['args']]
    return np.logical_and.reduce(vals) if n['op']=='all' else np.logical_or.reduce(vals)
def load_contract(config_path,criteria_path,inventory_path):
    config=json.loads(Path(config_path).read_text());criteria=json.loads(Path(criteria_path).read_text())
    if config.get('status')!='USER_APPROVED_FOR_EXECUTION':raise ValueError('Runtime approval absent')
    if config.get('preserve_flagged') is not True or config.get('apply_nucleus_filter') is not False or config.get('apply_tissue_filter') is not False or config.get('extraction_population')!='all_base_candidates_including_flagged':raise ValueError('Preserve-all/separate-eligibility contract mismatch')
    if sha(criteria_path)!=REVIEW_SHA or config.get('criteria_sha256')!=REVIEW_SHA:raise ValueError('Reviewed criteria hash mismatch')
    if config.get('source_inventory_sha256')!=sha(inventory_path):raise ValueError('Source inventory hash mismatch')
    if config.get('per_gene_minimum',2)!=2:raise ValueError('Runtime threshold differs from approved per-gene >=2')
    rules=criteria['foreign_programs']
    if [r['id'] for r in rules]!=RULE_IDS:raise ValueError('Fifteen reviewed rule identities/order changed')
    required=node_genes(criteria['base_predicate'])|{'TRAC','TRBC2'}
    for r in rules:required|=node_genes(r['predicate'])
    genes=CORE+sorted(required-set(CORE))
    if len(genes)!=67:raise ValueError('Expected exactly 67 approved/audit genes')
    if 'required_genes' in config and set(config['required_genes'])!=set(genes):raise ValueError('Required gene set mismatch')
    return config,criteria,genes
def get_pin(path,region):
    inventory=json.loads(Path(path).read_text())
    pins=[p for p in inventory.get('runs',inventory.get('sources',[])) if p['run_id']==region]
    if len(pins)!=1:raise ValueError('Region must occur exactly once in source inventory')
    pin=pins[0]
    if 'sources' not in pin:return pin
    src=pin['sources'];meta=[k for k in src if k in {'cells.csv.gz','cells.csv','cells.parquet'}]
    if len(meta)!=1:raise ValueError('Canonical inventory needs exactly one native metadata source')
    def oldform(x):return dict(path=x['path'],sha256=x['sha256'],before={'size':x['bytes'],'mtime_ns':x['mtime_ns']})
    return dict(run_id=pin['run_id'],bundle=pin['bundle'],native_h5_source=oldform(src['cell_feature_matrix.h5']),native_metadata_source=oldform(src[meta[0]]),native_unit=pin.get('native_unit',{}))
def check_source(path,pin):
    before=stat(path)
    if str(Path(path))!=pin['path'] or before['size']!=pin['before']['size'] or before['mtime_ns']!=pin['before']['mtime_ns'] or sha(path)!=pin['sha256']:raise ValueError('Native source path/size/mtime/SHA mismatch: '+str(path))
    return before
def read_metadata(path,wanted):
    wanted=set(wanted);parts=[]
    if Path(path).suffix=='.parquet':
        f=pq.ParquetFile(path);cols=[c for c in META if c in f.schema_arrow.names]
        batches=(b.to_pandas() for b in f.iter_batches(columns=cols,batch_size=100000))
    else:batches=pd.read_csv(path,dtype={'cell_id':str},usecols=lambda c:c in META,chunksize=100000,float_precision='round_trip')
    for d in batches:
        d['cell_id']=strings(d.cell_id);parts.append(d[d.cell_id.isin(wanted)])
    if not parts:raise ValueError('Native metadata has no schema/batches')
    d=pd.concat(parts,ignore_index=True)
    # CSV round trips must preserve the exact numeric value used for centroids.
    # Widen any native float32 field before serialization instead of shortening it.
    for c in d.select_dtypes(include=['float32']).columns:d[c]=d[c].astype(np.float64)
    if d.cell_id.duplicated().any() or set(d.cell_id)!=wanted:raise ValueError('Native metadata missing/duplicate selected IDs')
    if 'nucleus_count' not in d:d['nucleus_count']=np.nan
    for k in ['x_centroid','y_centroid']:
        if k not in d or not np.isfinite(d[k].to_numpy(float)).all():raise ValueError('Native centroid absent/nonfinite; no cells silently dropped')
    return d
def historical(path,region,expected_hash):
    if sha(path)!=expected_hash:raise ValueError('Historical membership SHA mismatch')
    d=pd.read_csv(path,dtype={'run_id':str,'cell_id':str},float_precision='round_trip')
    if 'run_id' in d:d=d[d.run_id==region]
    if d.cell_id.duplicated().any():raise ValueError('Historical IDs duplicated within region')
    return set(d.cell_id)
def audit_counts(t,gamma,flags,alpha):
    base=t&gamma;union=flags.any(axis=1);remaining_t=t.copy();remaining_b=base.copy()
    rules=[];pairs=[]
    for i,name in enumerate(RULE_IDS):
        f=flags[:,i]
        rules.append(dict(id=name,all_raw=int(f.sum()),t_gate=int((t&f).sum()),base=int((base&f).sum()),ordered_newly_flagged_t=int((remaining_t&f).sum()),ordered_newly_flagged_base=int((remaining_b&f).sum())))
        remaining_t&=~f;remaining_b&=~f
        rules[-1].update(ordered_t_unflagged_remaining=int(remaining_t.sum()),ordered_base_unflagged_remaining=int(remaining_b.sum()))
        for j in range(i):pairs.append(dict(rule_a=RULE_IDS[j],rule_b=name,all_raw=int((f&flags[:,j]).sum()),t_gate=int((t&f&flags[:,j]).sum()),base=int((base&f&flags[:,j]).sum())))
    return dict(all_raw_cells=len(t),t_gate=int(t.sum()),t_gate_foreign_flagged=int((t&union).sum()),t_gate_unflagged=int((t&~union).sum()),receptor_or_all_raw=int(gamma.sum()),base_candidates=int(base.sum()),base_flagged=int((base&union).sum()),base_unflagged=int((base&~union).sum()),base_alpha_beta_flag=int((base&alpha).sum()),extraction_eligible=int(base.sum()),named_rules=rules,pairwise_overlaps=pairs,ordered_counts_are_noncausal_accounting_only=True,nucleus_filter_applied=False,tissue_filter_applied=False)

def run(a):
    if 'Islet-Explorer-Senior-v3' in sys.executable:raise ValueError('Use isolated technical runtime, not v3')
    config,criteria,genes=load_contract(a.config,a.criteria,a.inventory);pin=get_pin(a.inventory,a.region)
    code_sha=sha(__file__);config_sha=sha(a.config);inventory_sha=sha(a.inventory)
    h5=Path(pin['native_h5_source']['path']);mp=Path(pin['native_metadata_source']['path'])
    if h5.parent!=Path(pin['bundle']) or mp.parent!=h5.parent:raise ValueError('Sources must belong to exact pinned bundle')
    hb=check_source(h5,pin['native_h5_source']);mb=check_source(mp,pin['native_metadata_source'])
    old_and=historical(a.historical_and_cells,a.region,config['historical_and_sha256_by_region'][a.region])
    old_or=historical(a.historical_or2_cells,a.region,config['historical_or2_sha256'])
    out=Path(a.out_dir)
    if out.exists() and any(out.iterdir()):raise FileExistsError('Refusing nonempty selector output')
    out.mkdir(parents=True,exist_ok=True)
    all_counts=[];all_flags=[];all_t=[];all_gamma=[];all_alpha=[];tmat=[];trows=[];t_totals=[]
    with h5py.File(h5,'r') as f:
        m=f['matrix'];nf,nc=map(int,m['shape'][:]);ids=strings(m['barcodes'][:]);ptr=m['indptr'][:]
        if len(ids)!=nc or len(set(ids))!=nc or len(ptr)!=nc+1 or ptr[0]!=0 or np.any(np.diff(ptr)<0) or ptr[-1]!=len(m['data']):raise ValueError('Malformed native barcode/CSC structure')
        fd={'feature_index':np.arange(nf)}
        for name,node in m['features'].items():
            if isinstance(node,h5py.Dataset) and node.ndim==1 and len(node)==nf:fd['feature_'+name]=strings(node[:])
        features=pd.DataFrame(fd)
        if 'feature_id' not in features or features.feature_id.duplicated().any():raise ValueError('Native feature IDs absent/duplicated')
        indices=[]
        for gene in genes:
            hit=np.flatnonzero(features.feature_name.to_numpy()==gene)
            if len(hit)!=1 or features.iloc[hit[0]].feature_feature_type!='Gene Expression':raise ValueError('Missing/duplicate/non-RNA feature: '+gene)
            indices.append(int(hit[0]))
        col={g:i for i,g in enumerate(genes)}
        for start in range(0,nc,a.chunk_cells):
            stop=min(nc,start+a.chunk_cells);lo,hi=int(ptr[start]),int(ptr[stop]);vals=m['data'][lo:hi];ix=m['indices'][lo:hi]
            if not np.isfinite(vals).all() or np.any(vals<0) or np.any(vals!=np.floor(vals)) or np.any(ix<0) or np.any(ix>=nf):raise ValueError('Invalid native count/index')
            mat=sparse.csr_matrix((vals.astype(np.int64),ix,ptr[start:stop+1]-lo),shape=(stop-start,nf));mat.sum_duplicates()
            counts=mat[:,indices].toarray();t=(counts[:,col['CD3E']]>=2)|(counts[:,col['CD3G']]>=2);gamma=(counts[:,0]>=2)|(counts[:,1]>=2)
            flags=np.column_stack([evaluate(r['predicate'],counts,col) for r in criteria['foreign_programs']]);alpha=(counts[:,col['TRAC']]>=2)&(counts[:,col['TRBC2']]>=2)
            all_counts.append(counts);all_flags.append(flags);all_t.append(t);all_gamma.append(gamma);all_alpha.append(alpha)
            if t.any():tmat.append(mat[t]);trows.append(start+np.flatnonzero(t));t_totals.append(np.asarray(mat[t].sum(axis=1)).ravel())
    if stat(h5)!=hb:raise RuntimeError('Native H5 changed during selection')
    counts=np.concatenate(all_counts);flags=np.concatenate(all_flags);t=np.concatenate(all_t);gamma=np.concatenate(all_gamma);alpha=np.concatenate(all_alpha)
    base=t&gamma;union=flags.any(axis=1);ix_t=np.concatenate(trows) if trows else np.array([],dtype=np.int64)
    tm=sparse.vstack(tmat,format='csr') if tmat else sparse.csr_matrix((0,nf),dtype=np.int64)
    frame=pd.DataFrame(counts[t],columns=genes);frame.insert(0,'native_matrix_row',ix_t);frame.insert(0,'cell_id',ids[t]);frame.insert(0,'run_id',a.region)
    frame['native_cell_key']=[pin['bundle']+'|'+c for c in frame.cell_id]
    frame['source_bundle']=pin['bundle'];frame['source_slide_id']=a.region;frame['source_region_name']=Path(pin['bundle']).name.split('__')[2]
    frame['CD3E_plus_CD3G']=frame.CD3E+frame.CD3G
    frame['native_all_feature_count']=np.concatenate(t_totals) if t_totals else np.array([],dtype=np.int64)
    frame['t_gate']=True;frame['receptor_or']=gamma[t];frame['base_candidate']=base[t];frame['foreign_flag_union']=union[t];frame['base_flagged']=(base&union)[t];frame['base_unflagged']=(base&~union)[t];frame['alpha_beta_flag']=alpha[t]
    for j,name in enumerate(RULE_IDS):frame['flag__'+name]=flags[t,j]
    meta=read_metadata(mp,frame.cell_id);frame=frame.merge(meta,on='cell_id',how='left',validate='one_to_one',sort=False)
    if stat(mp)!=mb:raise RuntimeError('Native metadata changed during read')
    frame['nucleus_status']=np.select([frame.nucleus_count.eq(0),frame.nucleus_count.eq(1),frame.nucleus_count.gt(1)],['zero','single','multiple'],default='unknown')
    frame['tissue_eligibility']='not_evaluated';frame['nucleus_used_for_selection']=False;frame['tissue_used_for_selection']=False
    eligible=frame[frame.base_candidate].copy().reset_index(drop=True)
    frame.to_csv(out/'t_gate_cells.csv.gz',index=False);eligible.to_csv(out/'eligible_raw_cells.csv',index=False)
    eligible[eligible.base_unflagged].to_csv(out/'base_unflagged.csv',index=False);eligible[eligible.base_flagged].to_csv(out/'base_flagged.csv',index=False)
    sparse.save_npz(out/'t_gate_native_counts.npz',tm);sparse.save_npz(out/'eligible_native_counts.npz',tm[frame.base_candidate.to_numpy()]);features.to_csv(out/'native_features.csv',index=False)
    np.savez_compressed(out/'whole_raw_predicate_counts.npz',counts=counts,cell_ids=np.asarray(ids,dtype=str),genes=np.asarray(genes,dtype=str),rule_ids=np.asarray(RULE_IDS,dtype=str),rule_flags=flags,t_gate=t,receptor_or=gamma,base_candidate=base,foreign_flag_union=union,alpha_beta_flag=alpha)
    np.savez_compressed(out/'whole_raw_gate_counts.npz',counts=counts[:,:4],cell_ids=np.asarray(ids,dtype=str),genes=np.asarray(CORE,dtype=str))
    audit=audit_counts(t,gamma,flags,alpha);dump(out/'gate_counts.json',audit)
    # Every exact overlap pattern is preserved; no hierarchy is used for assignment.
    overlap=[]
    bits=(flags.astype(np.uint32)*(1<<np.arange(15,dtype=np.uint32))).sum(axis=1)
    for label,mask in [('all_raw',np.ones(nc,dtype=bool)),('t_gate',t),('base',base)]:
        patterns,nums=np.unique(bits[mask],return_counts=True)
        overlap += [dict(universe=label,bitmask=int(b),rule_ids=';'.join(RULE_IDS[j] for j in range(15) if int(b)&(1<<j)),cells=int(n)) for b,n in zip(patterns,nums)]
    pd.DataFrame(overlap).to_csv(out/'flag_overlap_patterns.csv',index=False)
    idrow={c:i for i,c in enumerate(ids)}
    if not (old_and|old_or)<=set(ids):raise ValueError('Historical native IDs absent from pinned H5')
    keys=sorted(old_and|old_or|set(ids[base]),key=idrow.get);changes=[]
    for c in keys:
        i=idrow[c];changes.append(dict(run_id=a.region,cell_id=c,native_matrix_row=i,native_cell_key=pin['bundle']+'|'+c,historical_and=c in old_and,historical_or2=c in old_or,new_base=bool(base[i]),new_base_flagged=bool(base[i]&union[i]),new_base_unflagged=bool(base[i]&~union[i]),**{g:int(counts[i,j]) for j,g in enumerate(CORE)}))
    pd.DataFrame(changes,columns=['run_id','cell_id','native_matrix_row','native_cell_key','historical_and','historical_or2','new_base','new_base_flagged','new_base_unflagged']+CORE).to_csv(out/'historical_membership_changes.csv',index=False)
    diffs={}
    for name,old in [('historical_and',old_and),('historical_or2',old_or)]:
        new=set(ids[base]);diffs[name]=dict(old_cells=len(old),new_base_cells=len(new),shared=len(old&new),added=len(new-old),removed=len(old-new),exact_membership_match=old==new)
    dump(out/'historical_membership_summary.json',diffs)
    (out/'approved_runtime_config.json').write_bytes(Path(a.config).read_bytes());(out/'reviewed_criteria.json').write_bytes(Path(a.criteria).read_bytes())
    if sha(__file__)!=code_sha or sha(a.config)!=config_sha or sha(a.inventory)!=inventory_sha or sha(a.criteria)!=REVIEW_SHA:raise RuntimeError('Code/config/pins changed during run')
    if stat(h5)!=hb or stat(mp)!=mb:raise RuntimeError('Native sources changed before publication')
    if sha(a.historical_and_cells)!=config['historical_and_sha256_by_region'][a.region] or sha(a.historical_or2_cells)!=config['historical_or2_sha256']:raise RuntimeError('Historical membership inputs changed during run')
    outputs={p.name:sha(p) for p in out.iterdir() if p.is_file()}
    def source_doc(p,pin,before):return dict(path=str(p),sha256=pin['sha256'],bytes=before['size'],mtime_ns=before['mtime_ns'],before=before,after=stat(p))
    manifest=dict(schema='gammadelta-corrected-selector/1',status='COMPLETE_AWAITING_INDEPENDENT_VALIDATION',created_utc=datetime.now(timezone.utc).isoformat(),run_id=a.region,bundle=pin['bundle'],config_sha256=config_sha,predicate_config_sha256=config_sha,criteria_sha256=REVIEW_SHA,source_inventory_sha256=inventory_sha,h5_source=source_doc(h5,pin['native_h5_source'],hb),native_metadata_source=source_doc(mp,pin['native_metadata_source'],mb),historical_inputs={'and':{'path':str(Path(a.historical_and_cells).resolve()),'sha256':sha(a.historical_and_cells)},'or2':{'path':str(Path(a.historical_or2_cells).resolve()),'sha256':sha(a.historical_or2_cells)}},extraction_population='all_base_candidates_including_flagged',n_native_cells=nc,n_native_features=nf,n_t_gate_cells=len(frame),n_eligible_cells=len(eligible),n_base_flagged=int((base&union).sum()),n_base_unflagged=int((base&~union).sum()),required_genes=genes,rule_ids=RULE_IDS,nucleus_filter_applied=False,tissue_filter_applied=False,output_sha256=outputs,script_sha256=code_sha,runtime={'python':sys.version,'executable':sys.executable,'hostname':socket.gethostname(),'slurm_job_id':os.getenv('SLURM_JOB_ID')},biological_identity_accepted=False)
    dump(out/'selector_manifest.json',manifest)
    print(json.dumps({k:manifest[k] for k in ['run_id','n_native_cells','n_t_gate_cells','n_eligible_cells','n_base_flagged','n_base_unflagged']}))
    return manifest

def parser():
    p=argparse.ArgumentParser(description=__doc__)
    for name in ['inventory','config','criteria','historical-and-cells','historical-or2-cells']:p.add_argument('--'+name,required=True,type=Path)
    p.add_argument('--region','--run-id',dest='region',required=True);p.add_argument('--out-dir','--output-dir',dest='out_dir',required=True,type=Path);p.add_argument('--chunk-cells',type=int,default=20000)
    return p
if __name__=='__main__':
    args=parser().parse_args()
    if args.chunk_cells<1:raise ValueError('chunk-cells must be positive')
    run(args)
