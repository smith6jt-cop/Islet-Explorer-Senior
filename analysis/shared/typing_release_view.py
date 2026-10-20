#!/usr/bin/env python3
"""Read explicit, source-bound typing selections without modifying generations.

Scopes are analysis views, not a promise of distinct physical cells. In particular,
embedded pLN views can overlap the pancreas source. Counts must retain their scope.
"""
from pathlib import Path
import argparse
import hashlib
import json
from collections import Counter
import pyarrow as pa
import pyarrow.compute as pc
import pyarrow.parquet as pq


def sha256(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as f:
        for b in iter(lambda: f.read(8 << 20), b''):
            h.update(b)
    return h.hexdigest()


def load_registry(path, *, verify_bindings=True):
    registry = json.loads(Path(path).read_text())
    if registry.get('schema_version') != 1:
        raise ValueError('Unsupported typing registry schema')
    scopes = registry['scopes']
    ids = [s['scope_id'] for s in scopes]
    if len(ids) != len(set(ids)):
        raise ValueError('Duplicate analysis scope identifier')
    if registry.get('across_scope_counts_are_additive') is not False:
        raise ValueError('Registry must explicitly retain overlapping analysis views')
    if verify_bindings:
        for binding in registry['manifest_bindings']:
            if sha256(binding['path']) != binding['sha256']:
                raise ValueError('Selected source manifest changed: ' + binding['path'])
    for s in scopes:
        if not s['cell_key_column'] or not s['claim_column'] or not s['QC_column']:
            raise ValueError('A selected scope lacks an explicit consumer contract')
        if s['lane'] == 'pLN_RNA':
            if s.get('storage_scope') == 'section_specific_candidate':
                section=s.get('section_id')
                if (not section or s.get('filters', {}).get('section_id') != section
                        or s.get('native_source_filters', {}).get('section_key') != section
                        or not s.get('source_serial')
                        or s['cell_key_column'] != 'serial_qualified_cell_key'):
                    raise ValueError('Section-specific pLN candidate requires source-qualified identity and exact native scope')
            elif not s.get('filters', {}).get('section_key'):
                raise ValueError('Donor-pooled pLN RNA requires an exact section filter')
    return registry


def iter_scope_batches(scope, *, columns=None, verify_file=True, batch_size=131072):
    path = Path(scope['path'])
    if verify_file and sha256(path) != scope['sha256']:
        raise ValueError('Selected cell table changed: ' + str(path))
    needed = {scope['cell_key_column'], scope['claim_column'], scope['QC_column']}
    if scope.get('additional_admission_column'):
        needed.add(scope['additional_admission_column'])
    if columns is not None:
        needed.update(columns)
    needed.update(scope.get('filters', {}))
    pf = pq.ParquetFile(path)
    if not needed.issubset(pf.schema_arrow.names):
        raise ValueError('Selected cell table lacks contracted columns')
    requested = None if columns is None else sorted(needed)
    for batch in pf.iter_batches(batch_size=batch_size, columns=requested):
        table = pa.Table.from_batches([batch])
        for column, value in scope.get('filters', {}).items():
            table = table.filter(pc.equal(table[column], value))
        if table.num_rows:
            yield table


def scope_counts(scope, *, verify_file=True):
    total = eligible = admitted = 0
    all_claims, QC_claims = Counter(), Counter()
    # A null partition ID encodes native QC in the pLN RNA release. Other lanes
    # expose a preserved boolean flag; no identity call is used as a QC proxy.
    QC_kind = scope['QC_kind']
    if QC_kind not in ('boolean', 'nonnull'):
        raise ValueError('Unsupported QC contract')
    cols = [scope['claim_column'], scope['QC_column']]
    for table in iter_scope_batches(scope, columns=cols, verify_file=verify_file):
        claim = table[scope['claim_column']]
        q = table[scope['QC_column']]
        keep = pc.is_valid(q) if QC_kind == 'nonnull' else pc.fill_null(q, False)
        if QC_kind == 'boolean' and not pa.types.is_boolean(q.type):
            raise ValueError('QC flag must be boolean')
        total += table.num_rows
        eligible += pc.sum(pc.cast(keep, pa.int64())).as_py()
        if scope.get('additional_admission_column'):
            admission = table[scope['additional_admission_column']]
            if not pa.types.is_boolean(admission.type):
                raise ValueError('Additional admission flag must be boolean')
            admission = pc.fill_null(admission, False)
            if pc.any(pc.and_(admission, pc.invert(keep))).as_py():
                raise ValueError('Additional admission cannot admit native QC-ineligible rows')
            admitted += pc.sum(pc.cast(admission, pa.int64())).as_py()
        for record in pc.value_counts(claim).to_pylist():
            all_claims[record['values']] += record['counts']
        for record in pc.value_counts(pc.filter(claim, keep)).to_pylist():
            QC_claims[record['values']] += record['counts']
    if total != scope['expected_source_rows'] or eligible != scope['expected_QC_rows']:
        raise ValueError(f'Scope denominator mismatch for {scope["scope_id"]}: {total}/{eligible}')
    if scope.get('additional_admission_column') and admitted != scope['expected_additional_admission_rows']:
        raise ValueError('Additional admission denominator mismatch for ' + scope['scope_id'])
    return {'scope_id': scope['scope_id'], 'lane': scope['lane'],
            'claim_status': scope['claim_status'], 'source_rows': total, 'QC_rows': eligible,
            'all_row_claim_counts': [{'claim': k, 'n': v} for k,v in sorted(all_claims.items(), key=lambda x:str(x[0]))],
            'QC_claim_counts': [{'claim': k, 'n': v} for k,v in sorted(QC_claims.items(), key=lambda x:str(x[0]))],
            'additional_admission_rows':admitted if scope.get('additional_admission_column') else None,
            'additional_admission_semantics':scope.get('additional_admission_semantics'),
            'counts_across_scopes_additive': False,
            'biological_identity_accuracy_validated': False}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--registry', type=Path, required=True)
    sub = parser.add_subparsers(dest='command', required=True)
    sub.add_parser('inventory')
    count = sub.add_parser('counts')
    count.add_argument('--scope-id', required=True)
    audit = sub.add_parser('audit-counts')
    audit.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    registry = load_registry(args.registry)
    if args.command == 'inventory':
        print(json.dumps({'selection_status':registry['selection_status'],
                          'scopes':[{k:s[k] for k in ('scope_id','lane','donor_id','claim_status','expected_source_rows','expected_QC_rows')} for s in registry['scopes']],
                          'across_scope_counts_are_additive':False}, indent=2))
    elif args.command == 'counts':
        match = [s for s in registry['scopes'] if s['scope_id'] == args.scope_id]
        if len(match) != 1:
            raise ValueError('Scope must identify exactly one analysis view')
        print(json.dumps(scope_counts(match[0]), indent=2))
    else:
        receipts, verified = [], set()
        for scope in registry['scopes']:
            # Shared donor-pooled files retain independent section filters, while
            # their actual immutable bytes only need hashing once per invocation.
            binding = (scope['path'],scope['sha256'])
            receipts.append(scope_counts(scope,verify_file=binding not in verified))
            verified.add(binding)
            print(scope['scope_id'], 'consumer denominator readback passed', flush=True)
        result = {'status':'all_selected_scope_consumer_counts_verified',
                  'registry_path':str(args.registry.resolve()),'registry_sha256':sha256(args.registry),
                  'consumer_sha256':sha256(Path(__file__)), 'n_scopes':len(receipts),
                  'scopes':receipts,'across_scope_counts_are_additive':False,
                  'scientific_goal_complete':False,'biological_identity_accuracy_validated':False}
        args.output.parent.mkdir(parents=True,exist_ok=True)
        args.output.write_text(json.dumps(result,indent=2)+'\n')


if __name__ == '__main__':
    main()
