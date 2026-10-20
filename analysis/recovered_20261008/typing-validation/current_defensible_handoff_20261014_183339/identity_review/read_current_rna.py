#!/usr/bin/env python3
"""Read pinned RNA annotation sidecars with exact native nucleus evidence.

No labels, QC, native measurements, keys or geometry are changed. Biological
filtering is explicit, separate from native QC and named-lineage admission.
"""
from pathlib import Path
import hashlib
import json
import pyarrow as pa
import pyarrow.compute as pc
import pyarrow.parquet as pq

INVENTORY = Path(__file__).with_name('selected_current_rna_inventory.json')
APPENDED_COLUMNS = (
    'handoff_exact_native_cell_id', 'handoff_native_section_key',
    'handoff_native_nucleus_count', 'handoff_native_nucleus_measured',
    'handoff_positive_nucleus', 'handoff_native_QC_eligible',
    'handoff_named_lineage_admission', 'handoff_biological_use_eligible',
    'handoff_named_lineage_use_eligible',
)

def sha256(path):
    h=hashlib.sha256()
    with Path(path).open('rb') as f:
        for block in iter(lambda:f.read(8<<20),b''):h.update(block)
    return h.hexdigest()

def require(value,message):
    if not value:raise ValueError(message)

def load_inventory(path=INVENTORY):
    d=json.loads(Path(path).read_text())
    require(d.get('status')=='CURRENT_BYTES_AND_EXACT_NATIVE_SECTION_NUCLEUS_BINDINGS_VERIFIED',
            'A verified current RNA sidecar/native nucleus inventory is required')
    require(d.get('across_scope_counts_additive') is False,'Scopes must retain nonadditivity')
    scopes=d['scopes']
    require(len({s['scope_id'] for s in scopes})==len(scopes),'Duplicate RNA scope IDs')
    return {s['scope_id']:s for s in scopes}

def read_current_rna(scope, *, biological_only=True, named_lineage_only=False,
                     columns=None, verify_hashes=True):
    """Return one Arrow table; default excludes missing/nonpositive nuclei.

    ``scope`` is one entry from the verified inventory. ``biological_only=False``
    is audit/provenance use and retains every native source row. For named-lineage
    estimates use ``named_lineage_only=True``; this additionally enforces the
    existing admission flag. Original columns/values and serialized keys remain.
    The full scope is validated before any filtered table is returned.
    """
    require(scope['lane'] in ('pancreas_RNA','pLN_RNA'),'RNA scope required')
    require(scope['cell_key_column']=='serial_qualified_cell_key','Exact qualified key required')
    require(isinstance(biological_only,bool) and isinstance(named_lineage_only,bool),'Boolean use modes required')
    candidate=Path(scope['path']);native=Path(scope['native_source_path'])
    if verify_hashes:
        require(sha256(candidate)==scope['sha256'],'Candidate sidecar SHA256 mismatch')
        require(sha256(native)==scope['native_source_sha256'],'Native call SHA256 mismatch')
    pf=pq.ParquetFile(candidate);nf=pq.ParquetFile(native)
    required={scope['cell_key_column'],scope['claim_column'],scope['QC_column'],
              scope['additional_admission_column'],'cell_id'}
    filters=scope.get('candidate_filters',{})
    required.update(filters)
    require(required.issubset(pf.schema_arrow.names),'Missing contracted annotation fields')
    require(not set(APPENDED_COLUMNS).intersection(pf.schema_arrow.names),'Derived column collision')
    requested=None if columns is None else sorted(required|set(columns))
    a=pf.read(columns=requested).combine_chunks()
    for col,value in filters.items():
        require(pc.all(pc.equal(a[col],value)).as_py() is True,'Candidate contains rows outside exact scope')
    require(len(a)==scope['expected_source_rows'],'Candidate source-row denominator mismatch')
    for col in ('cell_id',scope['cell_key_column']):
        require(a[col].null_count==0 and len(pc.unique(a[col]))==len(a),
                'Candidate exact keys are missing or duplicated')
    native_columns=['cell_id','section_key','nucleus_count','nucleus_measured']
    require(set(native_columns).issubset(nf.schema_arrow.names),'Explicit native nucleus evidence unavailable')
    b=nf.read(columns=native_columns).combine_chunks()
    section=scope['section_id']
    native_filters=scope.get('native_source_filters',{'section_key':section})
    require(native_filters=={'section_key':section},'Exact native section filter required')
    b=b.filter(pc.equal(b['section_key'],section)).combine_chunks()
    require(len(b)==len(a),'Candidate/native exact section row-count mismatch')
    require(b['cell_id'].null_count==0 and len(pc.unique(b['cell_id']))==len(b),
            'Native exact section IDs missing or duplicated')
    index=pc.index_in(a['cell_id'],value_set=b['cell_id'])
    require(index.null_count==0,'Candidate exact cell ID missing from native section')
    require(pa.types.is_boolean(b['nucleus_measured'].type),'Native nucleus_measured must be boolean')
    counts=pc.take(b['nucleus_count'],index)
    measured=pc.take(b['nucleus_measured'],index)
    numeric=pc.cast(counts,pa.float64())
    positive=pc.fill_null(pc.and_kleene(pc.equal(measured,True),
        pc.and_kleene(pc.is_finite(numeric),pc.greater(numeric,0))),False)
    qcfield=a[scope['QC_column']]
    if scope['QC_kind']=='nonnull':qc=pc.is_valid(qcfield)
    else:
        require(scope['QC_kind']=='boolean' and pa.types.is_boolean(qcfield.type),'Invalid native QC contract')
        qc=pc.fill_null(qcfield,False)
    admission=a[scope['additional_admission_column']]
    require(pa.types.is_boolean(admission.type),'Named-lineage admission must be boolean')
    adm=pc.fill_null(admission,False)
    def count(v):return pc.sum(pc.cast(v,pa.int64())).as_py() or 0
    require(count(qc)==scope['expected_QC_rows'],'Native QC denominator mismatch')
    require(count(adm)==scope['expected_additional_admission_rows'],'Named-lineage admission denominator mismatch')
    require(count(pc.and_(adm,pc.invert(qc)))==0,'Named admission includes native-QC ineligible rows')
    bio=pc.and_(qc,positive)
    named=pc.and_(bio,adm)
    derived=(pc.take(b['cell_id'],index),pc.take(b['section_key'],index),counts,measured,
             positive,qc,admission,bio,named)
    for name,values in zip(APPENDED_COLUMNS,derived,strict=True):a=a.append_column(name,values)
    if named_lineage_only:return a.filter(named)
    if biological_only:return a.filter(bio)
    return a

def iter_current_rna(scope, *, biological_only=True, named_lineage_only=False,
                     columns=None, verify_hashes=True, batch_size=131072):
    """Yield Arrow tables after validating the full exact scope and source hashes."""
    require(isinstance(batch_size,int) and batch_size>0,'Positive batch_size required')
    table=read_current_rna(scope,biological_only=biological_only,
        named_lineage_only=named_lineage_only,columns=columns,verify_hashes=verify_hashes)
    for batch in table.to_batches(max_chunksize=batch_size):
        yield pa.Table.from_batches([batch])

if __name__=='__main__':
    import argparse
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--inventory',type=Path,default=INVENTORY)
    p.add_argument('--scope-id',required=True)
    p.add_argument('--output',type=Path,required=True)
    p.add_argument('--audit-all-source-rows',action='store_true')
    p.add_argument('--named-lineage-only',action='store_true')
    args=p.parse_args()
    scope=load_inventory(args.inventory)[args.scope_id]
    table=read_current_rna(scope,biological_only=not args.audit_all_source_rows,
                           named_lineage_only=args.named_lineage_only)
    args.output.parent.mkdir(parents=True,exist_ok=True)
    pq.write_table(table,args.output,compression='zstd')
    print(json.dumps({'scope_id':args.scope_id,'rows':len(table),'output':str(args.output),
        'biological_only':not args.audit_all_source_rows,'named_lineage_only':args.named_lineage_only,
        'labels_or_source_data_changed':False,'across_scope_counts_additive':False}))
