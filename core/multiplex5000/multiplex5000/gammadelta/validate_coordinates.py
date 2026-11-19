#!/usr/bin/env python3
"""Independent source-to-export coordinate/union and native-H5 count validation.

Only artifact binding helpers are shared with the extractor. ROI membership uses
independent sorted-X interval geometry (not its KD tree), the complete original
Parquet stream is compared in order across every native field, and counts are
re-read directly from original H5. Biological eligibility comes exclusively from
the separately validated corrected selection; no old hardcoded AND gate exists.
"""
from __future__ import annotations
import argparse
from collections import Counter
from datetime import datetime, timezone
import json
from pathlib import Path
import sqlite3

import h5py
import numpy as np
import pandas as pd
import pyarrow as pa
import pyarrow.compute as pc
import pyarrow.parquet as pq
from scipy import sparse

from .extract_corrected import FLAGS, POPULATION, check_pin, dump_new, load_contract, require, sha, signature

def decode(values):
    return [x.decode('utf8') if isinstance(x,bytes) else str(x) for x in values]

def independent_roi(xy, centers, half=30.):
    """Closed squares evaluated by sorted-X intervals and explicit Y comparisons."""
    result = np.zeros(len(xy),dtype=bool)
    finite = np.flatnonzero(np.isfinite(xy).all(axis=1))
    order = finite[np.argsort(xy[finite,0],kind='stable')]
    x = xy[order,0]
    for cx,cy in centers:
        lo = np.searchsorted(x,np.nextafter(cx-half,-np.inf),side='left')
        hi = np.searchsorted(x,np.nextafter(cx+half,np.inf),side='right')
        candidates = order[lo:hi]
        result[candidates] |= (np.abs(xy[candidates,0]-cx) <= half) & (np.abs(xy[candidates,1]-cy) <= half)
    return result

class Cursor:
    """Read an exact number of output rows irrespective of Parquet batch boundaries."""
    def __init__(self,path,batch_size):
        self.pf = pq.ParquetFile(path)
        self.batches = iter(self.pf.iter_batches(batch_size=batch_size))
        self.batch = None
        self.offset = 0
        self.consumed = 0
    def take(self,n):
        pieces = []
        remaining = n
        while remaining:
            if self.batch is None or self.offset == self.batch.num_rows:
                self.batch = next(self.batches,None)
                self.offset = 0
                require(self.batch is not None,'Export is missing source-selected rows')
            count = min(remaining,self.batch.num_rows-self.offset)
            pieces.append(self.batch.slice(self.offset,count))
            self.offset += count
            self.consumed += count
            remaining -= count
        return pa.Table.from_batches(pieces,schema=self.pf.schema_arrow)

def exact_columns(observed,expected):
    require(observed.schema.equals(expected.schema,check_metadata=True),'Native output schema/metadata changed')
    for name in expected.column_names:
        a,b = observed[name].combine_chunks(),expected[name].combine_chunks()
        require(a.is_null().equals(b.is_null()),f'Native null mask differs: {name}')
        if pa.types.is_floating(a.type):
            valid = a.is_valid().to_numpy(zero_copy_only=False)
            av = a.to_numpy(zero_copy_only=False)[valid]
            bv = b.to_numpy(zero_copy_only=False)[valid]
            # Compare exact floating values including signed zero/NaN bit patterns.
            require(av.tobytes() == bv.tobytes(),f'Native floating values changed: {name}')
        else:
            require(a.equals(b),f'Native values changed: {name}')

def native_h5_selected(pin,cells):
    """Read selected columns directly from native features×cells CSC H5."""
    before = check_pin(pin,full_hash=True)
    indices = cells.native_matrix_row.to_numpy(np.int64)
    require(np.array_equal(indices,cells.native_matrix_row.to_numpy()),'Noninteger native cell indices')
    rr,cc,vv = [],[],[]
    with h5py.File(pin['path'],'r') as h:
        m = h['matrix']; shape = tuple(int(x) for x in m['shape'][:])
        require(len(shape)==2,'Native H5 matrix is not 2D')
        require(((indices>=0)&(indices<shape[1])).all(),'Native cell indices out of bounds')
        barcodes = decode(m['barcodes'][:])
        require(len(set(barcodes))==len(barcodes)==shape[1],'Native barcodes not unique/shape-aligned')
        require([barcodes[i] for i in indices] == cells.cell_id.tolist(),'Eligible selected IDs/indices differ from native H5')
        feature_names = decode(m['features']['name'][:])
        require(len(feature_names)==shape[0] and len(set(feature_names))==len(feature_names),'Ambiguous native feature names')
        indptr = m['indptr'][:]
        require(len(indptr)==shape[1]+1 and indptr[0]==0 and np.all(np.diff(indptr)>=0),'Malformed native CSC offsets')
        require(indptr[-1]==len(m['data'])==len(m['indices']),'Native CSC data length mismatch')
        for output_row,native_row in enumerate(indices):
            start,end = int(indptr[native_row]),int(indptr[native_row+1])
            fi = m['indices'][start:end].astype(np.int64)
            values = m['data'][start:end].astype(np.int64)
            require(((fi>=0)&(fi<shape[0])).all() and (values>=0).all(),'Invalid native feature count/index')
            rr.extend([output_row]*len(fi)); cc.extend(fi.tolist()); vv.extend(values.tolist())
    require(signature(pin['path'])==before,'Native H5 changed during validation')
    return sparse.csr_matrix((vv,(rr,cc)),shape=(len(cells),shape[0]),dtype=np.int64),feature_names

def validate(args):
    require(args.batch_size>0,'batch-size must be positive')
    c = load_contract(args)
    pins = c['inventory_record']['sources']
    before = {n:check_pin(pins[n]) for n in ('transcripts.parquet','cell_feature_matrix.h5',c['metadata_source_name'])}
    extraction_path = args.extraction_dir/'extraction_manifest.json'
    initial_extraction_sha = sha(extraction_path)
    extraction = json.loads(extraction_path.read_text())
    require(extraction['status']=='EXTRACTION_COMPLETE_PENDING_COORDINATE_VALIDATION','Unexpected extraction status')
    for key,value in c['bindings'].items():
        require(extraction[key]==value,f'Extraction contract mismatch: {key}')
    require(str(extraction['run_id'])==str(args.run_id),'Extraction run mismatch')
    require(extraction['extraction_population']==POPULATION,'Extraction dropped a candidate stratum')
    require(extraction['n_eligible_cells']==len(c['cells']),'Extraction candidate count mismatch')
    require(extraction['n_foreign_flagged_cells_included']==c['n_foreign_flagged'],'Flagged candidate inclusion mismatch')
    require(extraction['roi_half_width_um']==30.,'ROI changed')
    require(extraction['foreign_flag_filter'] is False and extraction['nucleus_filter'] is False and extraction['tissue_filter'] is False,'Unauthorized population filter')
    require(extraction['gene_filter_at_extraction'] is None and extraction['quality_filter_at_extraction'] is None,'Unauthorized transcript filter')
    require(extraction['source_sha256']==pins['transcripts.parquet']['sha256'],'Original transcript hash binding mismatch')
    require(Path(extraction['source_path']).resolve()==Path(pins['transcripts.parquet']['path']).resolve(),'Original transcript path binding mismatch')
    require(extraction['source_hash_status']=='full_file_sha256_verified_against_pinned_original','Transcript hash was not verified by extraction')
    for n,sig in before.items():
        require(sig==extraction['source_after'][n],f'Native source changed after extraction: {n}')
    for name,record in extraction['outputs'].items():
        require(Path(name).name==name,'Invalid extraction output basename')
        path=args.extraction_dir/name
        require(path.stat().st_size==record['bytes'] and sha(path)==record['sha256'],f'Extraction output hash mismatch: {name}')
    source_pf = pq.ParquetFile(pins['transcripts.parquet']['path'])
    output = Cursor(args.extraction_dir/'transcripts_roi_or_assigned.parquet',args.batch_size)
    native_columns = source_pf.schema_arrow.names
    expected_schema = pa.schema(list(source_pf.schema_arrow)+[pa.field(n,pa.bool_()) for n in FLAGS],metadata=source_pf.schema_arrow.metadata)
    require(output.pf.schema_arrow.equals(expected_schema,check_metadata=True),'Output does not retain exact native schema plus three flags')
    require(extraction['native_columns_retained']==native_columns,'Native columns manifest mismatch')
    require(extraction['source_rows']==source_pf.metadata.num_rows,'Native source row count mismatch')
    require(output.pf.metadata.num_rows==extraction['statistics'].get('exported_rows',0),'Output row count/manifest mismatch')
    native,names = native_h5_selected(pins['cell_feature_matrix.h5'],c['cells'])
    saved = sparse.load_npz(args.selection_dir/'eligible_native_counts.npz').tocsr()
    require(saved.shape==native.shape and (saved-native).nnz==0,'Selector selected matrix differs from original H5')
    saved_features = pd.read_csv(args.selection_dir/'native_features.csv',dtype=str,keep_default_na=False)
    require(saved_features.feature_name.tolist()==names,'Saved feature order differs from original H5')
    require(native.shape[1]==c['selector']['n_native_features'],'Native feature count differs from selector')
    cells = c['cells']; ci={s:i for i,s in enumerate(cells.cell_id)}; fi={s:i for i,s in enumerate(names)}
    args.output_dir.mkdir(parents=True,exist_ok=False)
    dbpath=args.output_dir/'transcript_id_uniqueness.sqlite'
    db=sqlite3.connect(dbpath)
    db.execute('CREATE TABLE seen (id TEXT PRIMARY KEY) WITHOUT ROWID')
    totals=Counter(); count20=Counter(); countall=Counter(); unknown=Counter()
    try:
        # A mathematically empty union needs no source row iteration.
        batches=source_pf.iter_batches(batch_size=args.batch_size,use_threads=True) if len(cells) else []
        for number,batch in enumerate(batches,1):
            totals['source_rows_scanned']+=batch.num_rows
            xy=np.column_stack([batch[n].to_numpy(zero_copy_only=False) for n in ('x_location','y_location')]).astype(float)
            roi=independent_roi(xy,c['xy'])
            source_ids=batch['cell_id'].to_pylist()
            assigned=np.fromiter((x in ci for x in source_ids),dtype=bool,count=batch.num_rows)
            keep=roi|assigned
            if not keep.any():
                continue
            expected=pa.Table.from_batches([batch.filter(pa.array(keep))])
            observed=output.take(expected.num_rows)
            exact_columns(observed.select(native_columns),expected)
            # Actual output flags are independently recomputed from untouched source fields.
            qv=np.asarray(expected['qv'].to_pylist(),dtype=float)
            highq=qv>=20
            expected_flags=[roi[keep],assigned[keep],highq]
            for flag,value in zip(FLAGS,expected_flags):
                require(observed[flag].null_count==0 and np.array_equal(observed[flag].to_numpy(zero_copy_only=False),value),f'Derived flag mismatch: {flag}')
            d=observed.select(['transcript_id','cell_id','feature_name']).to_pydict()
            tid=d['transcript_id']
            require(all(x is not None and str(x)!='' for x in tid),'Missing transcript ID')
            try:
                db.executemany('INSERT INTO seen(id) VALUES (?)',[(str(x),) for x in tid])
            except sqlite3.IntegrityError as exc:
                raise ValueError('Duplicate transcript IDs in extraction union') from exc
            db.commit()
            totals['exported_rows']+=observed.num_rows
            totals['rows_in_rois']+=int(expected_flags[0].sum())
            totals['rows_assigned_to_eligible_cells']+=int(expected_flags[1].sum())
            totals['assigned_rows_outside_all_rois']+=int((expected_flags[1]&~expected_flags[0]).sum())
            totals['rows_qv_ge20']+=int(highq.sum())
            totals['rows_qv_lt20_or_missing']+=int((~highq).sum())
            for j in np.flatnonzero(expected_flags[1]):
                cell_index=ci[d['cell_id'][j]]; feature=d['feature_name'][j]
                if feature not in fi:
                    unknown[(str(feature),bool(highq[j]))]+=1
                    continue
                key=(cell_index,fi[feature]);countall[key]+=1
                if highq[j]:
                    count20[key]+=1
            if number%100==0:
                print(json.dumps({'run_id':args.run_id,'stage':'independent_coordinate_validation','batches':number,**dict(totals)}),flush=True)
    finally:
        db.close()
    require(totals['source_rows_scanned']==(source_pf.metadata.num_rows if len(cells) else 0),'Independent original-source scan incomplete')
    require(output.consumed==output.pf.metadata.num_rows,'Export contains extra rows not in exact source union')
    for key in ['exported_rows','rows_in_rois','rows_assigned_to_eligible_cells','assigned_rows_outside_all_rois','rows_qv_ge20','rows_qv_lt20_or_missing']:
        require(totals[key]==extraction['statistics'].get(key,0),f'Extraction statistic mismatch: {key}')
    def to_sparse(counter):
        keys=list(counter)
        return sparse.csr_matrix(([counter[k] for k in keys],([k[0] for k in keys],[k[1] for k in keys])),shape=native.shape,dtype=np.int64)
    observed20=to_sparse(count20);observedall=to_sparse(countall);delta=(observed20-native).tocoo()
    sparse.save_npz(args.output_dir/'reconciled_assigned_q20_counts.npz',observed20)
    sparse.save_npz(args.output_dir/'reconciled_assigned_all_qv_counts.npz',observedall)
    mismatches=pd.DataFrame({'cell_id':cells.cell_id.to_numpy()[delta.row],'feature_name':np.asarray(names)[delta.col],'observed_minus_native':delta.data})
    mismatches.to_csv(args.output_dir/'count_mismatches.csv',index=False)
    pd.DataFrame([{'feature_name':g,'qv_ge20':q,'assigned_rows':n} for (g,q),n in sorted(unknown.items())],columns=['feature_name','qv_ge20','assigned_rows']).to_csv(args.output_dir/'features_absent_from_native_matrix.csv',index=False)
    percell=cells[['run_id','cell_id','native_matrix_row','native_cell_key','foreign_flag_union']].copy()
    percell['native_all_feature_count']=np.asarray(native.sum(axis=1)).ravel()
    percell['assigned_q20_all_feature_count']=np.asarray(observed20.sum(axis=1)).ravel()
    percell['assigned_all_qv_feature_count']=np.asarray(observedall.sum(axis=1)).ravel()
    percell['n_mismatched_features']=np.bincount(delta.row,minlength=len(cells))
    percell.to_csv(args.output_dir/'per_cell_count_reconciliation.csv',index=False)
    require(delta.nnz==0,f'Native H5 reconciliation failed: {delta.nnz} cell-feature pairs')
    require(not any(n for (g,q),n in unknown.items() if q),'Assigned QV>=20 features absent from native H5')
    counts_csv=pd.read_csv(args.extraction_dir/'qv_ge20_assigned_counts.csv',dtype={'cell_id':str,'feature_name':str},keep_default_na=False)
    require(not counts_csv.duplicated(['cell_id','feature_name']).any(),'Duplicate count summary keys')
    expected_count_map={(cells.cell_id.iloc[i],names[j]):v for (i,j),v in count20.items()}
    actual_count_map={(r.cell_id,r.feature_name):int(r.qv_ge20_assigned_count) for r in counts_csv.itertuples()}
    require(actual_count_map==expected_count_map,'Extraction count summary differs from direct native transcript tally')
    for name,sig in before.items():
        require(signature(pins[name]['path'])==sig,f'Original source changed during validation: {name}')
    require(load_contract(args)['bindings']==c['bindings'],'Selection contract changed during validation')
    require(sha(extraction_path)==initial_extraction_sha,'Extraction manifest changed')
    for name,record in extraction['outputs'].items():
        require(sha(args.extraction_dir/name)==record['sha256'],f'Extraction output changed during validation: {name}')
    report={'schema':'corrected-native-coordinate-validation-v1','status':'PASS','all_checks_pass':True,
        'completed_utc':datetime.now(timezone.utc).isoformat(),'run_id':str(args.run_id),**c['bindings'],
        'native_unit':c['inventory_record'].get('native_unit',{}),
        'extraction_manifest_sha256':sha(extraction_path),'extraction_population':POPULATION,
        'n_eligible_cells':len(cells),'n_foreign_flagged_cells_included':c['n_foreign_flagged'],
        'n_native_features':len(names),'n_cell_feature_pairs_tested_including_zeros':int(native.shape[0]*native.shape[1]),
        'n_cell_feature_mismatches':0,'statistics':dict(totals),
        'checks':{'independently_verified_eligible_ids_bound':True,'flagged_candidates_preserved':True,'all_native_columns_exact_source_values':True,'exact_complete_ROI_or_assignment_union':True,'closed_square_geometry_independently_recomputed':True,'no_duplicate_or_empty_transcript_ids':True,'all_feature_assigned_QV20_counts_match_original_H5':True,'source_and_output_hash_bindings':True,'no_nucleus_or_tissue_filter':True},
        'source_verification':{'transcripts':'Full source SHA verified at extraction; unchanged file signature checked, then full source union re-streamed with every selected native field compared.','h5':'Full pinned SHA independently rechecked; selected CSC columns read directly.','cells':'Pinned source/selector validation receipt and unchanged signature; centroids independently checked by upstream selection validator.'},
        'script_sha256':sha(__file__),'binding_helper_sha256':sha(Path(__file__).with_name('extract_corrected.py')),
        'output_sha256':{p.name:sha(p) for p in args.output_dir.iterdir() if p.is_file()},
        'limits':['Technical equality does not validate biological identity, segmentation purity or transcript origin.','Native coordinate units follow the original bundle provenance; no coordinate transform is applied.','Nonfinite source values, if present in assigned rows, are preserved exactly rather than silently dropped.']}
    dump_new(args.output_dir/'coordinate_validation.json',report)
    print(json.dumps({'run_id':args.run_id,'status':'PASS','candidates':len(cells),'flagged_included':c['n_foreign_flagged'],'rows':totals['exported_rows']}),flush=True)
    return report

def parse_args():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--config',type=Path,required=True);p.add_argument('--inventory',type=Path,required=True)
    p.add_argument('--selection-dir',type=Path,required=True);p.add_argument('--selection-validation',type=Path)
    p.add_argument('--run-id','--region',dest='run_id',required=True)
    p.add_argument('--extraction-dir',type=Path,required=True)
    p.add_argument('--output-dir','--out-dir',dest='output_dir',type=Path,required=True)
    p.add_argument('--batch-size',type=int,default=262144)
    return p.parse_args()

if __name__=='__main__':
    a=parse_args();a._extraction_sha=sha(a.extraction_dir/'extraction_manifest.json')
    validate(a)
