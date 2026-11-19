#!/usr/bin/env python3
"""Export ALL native columns/QVs in 60um squares OR assigned to ALL base candidates.

The corrected selector's independent PASS receipt is mandatory. No biological
predicate, foreign-program flag, nucleus or tissue field is an extraction filter.
"""
from __future__ import annotations
import argparse
from collections import Counter
from datetime import datetime, timezone
import hashlib
import json
import os
from pathlib import Path
import platform
import time

import numpy as np
import pandas as pd
import pyarrow as pa
import pyarrow.compute as pc
import pyarrow.parquet as pq
from scipy.spatial import cKDTree

POPULATION = 'all_base_candidates_including_flagged'
FLAGS = ['in_any_60um_roi', 'is_assigned_to_eligible_cell', 'qv_ge20']
REQUIRED = ['transcript_id','cell_id','feature_name','x_location','y_location','z_location','qv','overlaps_nucleus']

def require(ok, message):
    if not ok:
        raise ValueError(message)

def sha(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as handle:
        for block in iter(lambda: handle.read(8 << 20), b''):
            h.update(block)
    return h.hexdigest()

def signature(path):
    s = Path(path).stat()
    return {'bytes':s.st_size, 'mtime_ns':s.st_mtime_ns, 'inode':s.st_ino, 'device':s.st_dev}

def dump_new(path, data):
    with Path(path).open('x') as f:
        json.dump(data, f, indent=2, allow_nan=False)
        f.write('\n')

def check_pin(record, full_hash=False):
    path = Path(record['path'])
    sig = signature(path)
    require(sig['bytes'] == record['bytes'], f'Pinned size differs: {path}')
    require(sig['mtime_ns'] == record['mtime_ns'], f'Pinned mtime differs: {path}')
    require(len(record['sha256']) == 64, f'Missing full source hash: {path}')
    if full_hash:
        require(sha(path) == record['sha256'], f'Pinned SHA256 differs: {path}')
        require(signature(path) == sig, f'Source changed during hashing: {path}')
    return sig

def load_contract(args):
    """Shared artifact binding only; selection logic is verified upstream independently."""
    config_path, inv_path, sel = args.config, args.inventory, args.selection_dir
    config = json.loads(config_path.read_text())
    require(config.get('status') == 'USER_APPROVED_FOR_EXECUTION', 'Runtime config not USER_APPROVED_FOR_EXECUTION')
    require(config.get('extraction_population') == POPULATION and config.get('preserve_flagged') is True, 'Approved config must preserve ALL base candidates')
    require(config.get('apply_nucleus_filter') is False and config.get('apply_tissue_filter') is False, 'Approved config must not filter nucleus/tissue annotations')
    require(config.get('roi_half_width_um') == 30. and config.get('qv_reconciliation_minimum') == 20., 'Unexpected approved ROI or reconciliation threshold')
    require(config.get('preserve_native_transcript_columns') is True, 'Approved config must retain native columns')
    require(config.get('source_inventory_sha256') == sha(inv_path), 'Approved config bound to different source inventory')
    inventory = json.loads(inv_path.read_text())
    matches = [x for x in inventory['runs'] if str(x['run_id']) == str(args.run_id)]
    require(len(matches) == 1, 'Inventory must identify exactly one native run')
    record = matches[0]
    require({'transcripts.parquet','cell_feature_matrix.h5'} <= set(record['sources']), 'Missing required native source pins')
    metadata_names = [n for n in record['sources'] if n in ('cells.parquet','cells.csv.gz','cells.csv')]
    require(len(metadata_names) == 1, 'Inventory must pin exactly one native metadata source')
    metadata_name = metadata_names[0]
    manifest_path = sel / 'selector_manifest.json'
    manifest = json.loads(manifest_path.read_text())
    receipt_path = args.selection_validation or sel / 'validation' / 'selection_validation.json'
    receipt = json.loads(receipt_path.read_text())
    eligible = sel / 'eligible_raw_cells.csv'
    require(receipt.get('status') == 'PASS' and receipt.get('all_checks_pass') is True, 'Independent selection validation did not pass')
    require(manifest.get('extraction_population') == POPULATION, 'Selector population is not ALL base candidates')
    require(str(manifest['run_id']) == str(args.run_id) == str(receipt['run_id']), 'Run ID disagreement')
    for key, actual in [('config_sha256',sha(config_path)), ('source_inventory_sha256',sha(inv_path))]:
        require(manifest[key] == actual == receipt[key], f'{key} binding mismatch')
    require(manifest['criteria_sha256'] == receipt['criteria_sha256'] == config['criteria_sha256'], 'Approved criteria binding mismatch')
    require(receipt['selector_manifest_sha256'] == sha(manifest_path), 'Selection receipt is stale')
    require(receipt['eligible_cells_sha256'] == sha(eligible), 'Eligible table differs from independently verified table')
    for name, digest in manifest['output_sha256'].items():
        require(Path(name).name == name, 'Selector output name must be a basename')
        require(sha(sel / name) == digest, f'Selector output hash differs: {name}')
    require(manifest['output_sha256'].get(eligible.name) == sha(eligible), 'Eligible CSV missing from selector output pins')
    for field, name in [('h5_source','cell_feature_matrix.h5'),('native_metadata_source',metadata_name)]:
        source = manifest[field]
        pin = record['sources'][name]
        require(Path(source['path']).resolve() == Path(pin['path']).resolve(), f'{field} source path disagreement')
        require(source['sha256'] == pin['sha256'], f'{field} source hash disagreement')
    cells = pd.read_csv(eligible, dtype={'run_id':str,'cell_id':str,'native_cell_key':str}, keep_default_na=False, float_precision='round_trip')
    require({'run_id','cell_id','native_matrix_row','native_cell_key','x_centroid','y_centroid','foreign_flag_union'} <= set(cells), 'Missing eligible columns')
    require(len(cells) == manifest['n_eligible_cells'], 'Eligible row count differs from selector manifest')
    require(not cells.cell_id.duplicated().any() and not cells.native_matrix_row.duplicated().any(), 'Duplicate eligible cell ID or native row')
    require(not cells.native_cell_key.duplicated().any(), 'Duplicate full native cell key')
    require(cells.run_id.eq(str(args.run_id)).all(), 'Eligible table contains another run')
    require(cells.cell_id.ne('').all(), 'Empty eligible cell ID')
    xy = cells[['x_centroid','y_centroid']].to_numpy(float)
    require(np.isfinite(xy).all(), 'Nonfinite centroid: review source instead of silently filtering a candidate')
    foreign = cells.foreign_flag_union.astype(str).str.lower()
    require(foreign.isin(['true','false','1','0']).all(), 'Malformed foreign_flag_union')
    return {'config':config,'inventory_record':record,'metadata_source_name':metadata_name,'selector':manifest,'receipt':receipt,
            'receipt_path':receipt_path,'manifest_path':manifest_path,'eligible_path':eligible,
            'cells':cells,'xy':xy,'n_foreign_flagged':int(foreign.isin(['true','1']).sum()),
            'bindings':{'config_sha256':sha(config_path),'criteria_sha256':manifest['criteria_sha256'],
                        'source_inventory_sha256':sha(inv_path),'selector_manifest_sha256':sha(manifest_path),
                        'selection_validation_sha256':sha(receipt_path),'eligible_cells_sha256':sha(eligible)}}

def parse_args():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--config',type=Path,required=True)
    p.add_argument('--inventory',type=Path,required=True)
    p.add_argument('--selection-dir',type=Path,required=True)
    p.add_argument('--selection-validation',type=Path)
    p.add_argument('--run-id','--region',dest='run_id',required=True)
    p.add_argument('--output-dir','--out-dir',dest='output_dir',type=Path,required=True)
    p.add_argument('--batch-size',type=int,default=262144)
    return p.parse_args()

def extract(args):
    require(args.batch_size > 0, 'batch-size must be positive')
    contract = load_contract(args)
    pins = contract['inventory_record']['sources']
    before = {name:check_pin(pins[name]) for name in ('transcripts.parquet','cell_feature_matrix.h5',contract['metadata_source_name'])}
    # The H5 and cells full hashes are verified by the mandatory upstream receipt.
    # Independently hash the complete original transcript source in this stage.
    transcript = pins['transcripts.parquet']
    check_pin(transcript, full_hash=True)
    source = Path(transcript['path'])
    pf = pq.ParquetFile(source)
    require(pf.metadata.num_rows == transcript['source_rows'], 'Native transcript row count differs from inventory')
    require(set(REQUIRED) <= set(pf.schema_arrow.names), 'Native transcript schema lacks required columns')
    require(not set(FLAGS) & set(pf.schema_arrow.names), 'Derived names collide with native columns')
    require(pa.types.is_string(pf.schema_arrow.field('cell_id').type) or pa.types.is_large_string(pf.schema_arrow.field('cell_id').type), 'Native cell_id must be string encoded')
    out = args.output_dir
    out.mkdir(parents=True, exist_ok=False)
    start = time.monotonic()
    cells, xy = contract['cells'], contract['xy']
    tree = cKDTree(xy) if len(cells) else None
    ids = pa.array(cells.cell_id.tolist(), type=pf.schema_arrow.field('cell_id').type)
    schema = pa.schema(list(pf.schema_arrow)+[pa.field(n,pa.bool_()) for n in FLAGS], metadata=pf.schema_arrow.metadata)
    target = out / 'transcripts_roi_or_assigned.parquet'
    stats = Counter()
    tally = Counter()
    pa.set_cpu_count(max(1,int(os.getenv('SLURM_CPUS_PER_TASK','2'))))
    pa.set_io_thread_count(2)
    # Empty eligibility deliberately produces a schema-correct empty table.
    batches = pf.iter_batches(batch_size=args.batch_size,use_threads=True) if len(cells) else []
    with pq.ParquetWriter(target,schema,compression='zstd') as writer:
        for number,batch in enumerate(batches,1):
            stats['source_rows_scanned'] += batch.num_rows
            coords = np.column_stack([batch[n].to_numpy(zero_copy_only=False) for n in ('x_location','y_location')])
            finite = np.isfinite(coords).all(axis=1)
            roi = np.zeros(batch.num_rows,dtype=bool)
            if finite.any():
                distance,_ = tree.query(coords[finite],k=1,p=np.inf,distance_upper_bound=np.nextafter(30.,np.inf),workers=1)
                roi[finite] = distance <= 30.
            assigned = pc.fill_null(pc.is_in(batch['cell_id'],value_set=ids),False).to_numpy(zero_copy_only=False)
            keep = roi | assigned
            stats['source_rows_nonfinite_xy'] += int((~finite).sum())
            if keep.any():
                chosen = pa.Table.from_batches([batch.filter(pa.array(keep))])
                selected_roi, selected_assigned = roi[keep], assigned[keep]
                highq = pc.fill_null(pc.greater_equal(chosen['qv'],20.),False)
                chosen = chosen.append_column(FLAGS[0],pa.array(selected_roi))
                chosen = chosen.append_column(FLAGS[1],pa.array(selected_assigned))
                chosen = chosen.append_column(FLAGS[2],highq)
                writer.write_table(chosen,row_group_size=131072)
                stats['exported_rows'] += chosen.num_rows
                stats['rows_in_rois'] += int(selected_roi.sum())
                stats['rows_assigned_to_eligible_cells'] += int(selected_assigned.sum())
                stats['assigned_rows_outside_all_rois'] += int((selected_assigned & ~selected_roi).sum())
                stats['rows_qv_ge20'] += int(pc.sum(pc.cast(highq,pa.int64())).as_py() or 0)
                stats['rows_qv_lt20_or_missing'] += chosen.num_rows - int(pc.sum(pc.cast(highq,pa.int64())).as_py() or 0)
                counted = chosen.filter(pc.and_(highq,pa.array(selected_assigned)))
                if counted.num_rows:
                    # Count every row, including control features; null IDs fail independent QA.
                    groups = counted.select(['cell_id','feature_name']).group_by(['cell_id','feature_name']).aggregate([('cell_id','count')]).to_pydict()
                    tally.update({(c,g):int(n) for c,g,n in zip(groups['cell_id'],groups['feature_name'],groups['cell_id_count'])})
            if number % 100 == 0:
                print(json.dumps({'run_id':args.run_id,'stage':'extraction','batches':number,**dict(stats)}),flush=True)
    require(stats['source_rows_scanned'] == (pf.metadata.num_rows if len(cells) else 0),'Incomplete native transcript scan')
    for name,sig in before.items():
        require(signature(pins[name]['path']) == sig, f'Original source changed during extraction: {name}')
    # Reject changed selection/config/validation artifacts after the long scan.
    end_contract = load_contract(args)
    require(end_contract['bindings'] == contract['bindings'],'Selection contract changed during extraction')
    count_path = out / 'qv_ge20_assigned_counts.csv'
    pd.DataFrame([{'cell_id':c,'feature_name':g,'qv_ge20_assigned_count':n} for (c,g),n in sorted(tally.items())],columns=['cell_id','feature_name','qv_ge20_assigned_count']).to_csv(count_path,index=False)
    report = {
        'schema':'corrected-native-coordinate-extraction-v1','status':'EXTRACTION_COMPLETE_PENDING_COORDINATE_VALIDATION',
        'completed_utc':datetime.now(timezone.utc).isoformat(),'run_id':str(args.run_id),
        'native_unit':contract['inventory_record'].get('native_unit',{}),
        'extraction_population':POPULATION,'n_eligible_cells':len(cells),'n_foreign_flagged_cells_included':contract['n_foreign_flagged'],
        **contract['bindings'],'selection_validation_path':str(contract['receipt_path'].resolve()),
        'source_path':str(source.resolve()),'source_sha256':transcript['sha256'],'source_hash_status':'full_file_sha256_verified_against_pinned_original',
        'source_before':before,'source_after':{name:signature(pins[name]['path']) for name in before},
        'source_rows':pf.metadata.num_rows,'source_row_groups':pf.metadata.num_row_groups,
        'native_columns_retained':pf.schema_arrow.names,'native_schema':str(pf.schema_arrow),'added_columns':FLAGS,
        'roi_half_width_um':30.,'extraction_rule':'Every native column and every QV within any centroid ±30 micrometre square OR assigned to ANY corrected base candidate, including all foreign-program-flagged candidates.',
        'coordinates_modified':False,'cell_assignments_modified':False,'gene_filter_at_extraction':None,'quality_filter_at_extraction':None,
        'foreign_flag_filter':False,'nucleus_filter':False,'tissue_filter':False,'counts_quality_threshold':20,
        'statistics':dict(stats),'batch_size':args.batch_size,'elapsed_seconds':time.monotonic()-start,
        'script_sha256':sha(__file__),'hostname':platform.node(),'pyarrow_version':pa.__version__,
        'slurm_job_id':os.getenv('SLURM_JOB_ID'),'outputs':{p.name:{'bytes':p.stat().st_size,'sha256':sha(p)} for p in [target,count_path]},
        'limitation':'Technical extraction only. ROI is not segmentation; foreign-program flags and technical metadata are preserved upstream, never population filters.'}
    dump_new(out/'extraction_manifest.json',report)
    print(json.dumps({'run_id':args.run_id,'status':report['status'],'candidates':len(cells),'flagged_included':contract['n_foreign_flagged'],'rows':stats['exported_rows']}),flush=True)
    return report

if __name__ == '__main__':
    extract(parse_args())
