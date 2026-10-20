#!/usr/bin/env python3
"""Explicitly abstain where the current producer assigned identity by complement/state.

This additive candidate does not infer new identities and uses no new marker gate.
It preserves every current assignment key, original call, QC, and admission flag.
"""
from pathlib import Path
from collections import Counter
import hashlib
import json
import resource
import time
import numpy as np
import pandas as pd
import pyarrow as pa
import pyarrow.parquet as pq

HERE = Path(__file__).resolve().parent
GEN = Path('/data/io60-calibration/own_compartment_frozen_reference_repair_20261008')
SOURCE = GEN / 'pancreas_analysis_assignments_current_generation'
OUT = HERE / 'conservative_assignment_candidate'
COLS = ['cell_id', 'object_id', 'donor_id', 'cell_class', 'refined_class',
        'native_cell_class', 'native_refined_class', 'resolved_by', 'immune_leaf',
        'qc_analysis_eligible', 'analysis_eligible', 'analysis_exclusion_reason',
        'mandatory_analysis_exclusion', 'is_endocrine', 'endocrine_subtype_status']


def sha(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as f:
        for b in iter(lambda: f.read(8 << 20), b''):
            h.update(b)
    return h.hexdigest()


def derive(d):
    qc = d.qc_analysis_eligible.fillna(False).to_numpy(bool)
    admitted = d.analysis_eligible.fillna(False).to_numpy(bool)
    active = qc & admitted
    native = d.cell_class.fillna('Unavailable').to_numpy(object)
    refined = d.refined_class.fillna('Unavailable').to_numpy(object)
    route = d.resolved_by.fillna('').to_numpy(object)
    leaf = d.immune_leaf.fillna(False).to_numpy(bool)
    candidate = native.copy()
    reason = np.full(len(d), 'affirmative_native_producer_route_retained', dtype=object)
    exo = active & (native == 'Exocrine')
    state_only = active & (native == 'Immune') & (route == 'supporting') & ~leaf
    candidate[exo | state_only] = 'Identity unresolved'
    reason[exo] = 'non-endocrine fallthrough; affirmative epithelial identity not established'
    reason[state_only] = 'supporting-marker-only immune; affirmative immune leaf absent'
    bad_route = active & (native == 'Immune') & ~np.isin(route, ['definition', 'supporting']) & ~leaf
    candidate[bad_route] = 'Identity unresolved'
    reason[bad_route] = 'immune route lacks recorded affirmative lineage evidence'
    reason[active & (native == 'Immune') & leaf] = 'existing_affirmative_immune_leaf_retained'
    candidate[~qc] = 'Unavailable'
    reason[~qc] = 'native_QC_ineligible'
    existing_exclusion = (native == 'Excluded') | (qc & ~admitted)
    candidate[existing_exclusion] = 'Excluded'
    reason[existing_exclusion] = 'native_analysis_exclusion_preserved'
    # Refinement can remove a claim without demonstrating the identity of its complement.
    refined_candidate = refined.copy()
    refined_reason = np.full(len(d), 'affirmative_native_producer_route_retained', dtype=object)
    refined_exo = active & (refined == 'Exocrine')
    refined_state_only = active & (refined == 'Immune') & ~leaf & (route != 'definition')
    refined_candidate[refined_exo | refined_state_only] = 'Identity unresolved'
    refined_reason[refined_exo & (native == 'Exocrine')] = 'non-endocrine fallthrough; affirmative epithelial identity not established'
    refined_reason[refined_exo & (native != 'Exocrine')] = 'demotion to exocrine complement; affirmative epithelial identity not established'
    refined_reason[refined_state_only] = 'supporting-marker-only immune; affirmative immune leaf absent'
    refined_candidate[~qc] = 'Unavailable'
    refined_reason[~qc] = 'native_QC_ineligible'
    refined_candidate[existing_exclusion] = 'Excluded'
    refined_reason[existing_exclusion] = 'native_analysis_exclusion_preserved'
    out = d.copy()
    out['candidate_cell_class'] = candidate
    out['candidate_identity_reason'] = reason
    out['candidate_refined_class'] = refined_candidate
    out['candidate_refined_identity_reason'] = refined_reason
    out['candidate_primary_changed'] = candidate != native
    out['candidate_refined_changed'] = refined_candidate != refined
    # A retained native route is affirmative producer evidence, not external truth.
    out['independent_identity_validation'] = False
    return out


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    manifest_path = SOURCE / 'manifest.json'
    original_manifest_sha = sha(manifest_path)
    manifest = json.loads(manifest_path.read_text())
    inputs = []
    groups = Counter()
    changed = Counter()
    n_all = 0
    started = time.time()
    status = {'status': 'running', 'created_UTC': time.strftime('%Y-%m-%dT%H:%M:%SZ', time.gmtime()),
              'source_manifest': str(manifest_path), 'source_manifest_sha256': original_manifest_sha,
              'sources': [], 'method': 'conservative claim correction from exhaustive current producer routes',
              'labels_in_original_generation_modified': False, 'production_alias_promoted': False,
              'new_marker_gate': False, 'independent_identity_accuracy_validated': False,
              'initial_candidate_exclusion_precedence_snapshot': str(HERE/'candidate_exclusion_precedence_initial_snapshot'),
              'exclusion_precedence': 'Existing Excluded claims, including QC-ineligible mandatory keys, remain Excluded; native QC/admission flags remain verbatim.'}
    (OUT / 'manifest.json').write_text(json.dumps(status, indent=2) + '\n')
    for directory in sorted(SOURCE.glob('donor_id=*')):
        source = directory / 'data.parquet'
        if not source.exists():
            continue
        donor = directory.name.split('=', 1)[1]
        target_dir = OUT / directory.name
        target_dir.mkdir(exist_ok=True)
        target = target_dir / 'data.parquet'
        affected = target_dir / 'changed_exact_keys.parquet'
        source_sha = sha(source)
        writer = affected_writer = None
        n = n_changed = 0
        for batch in pq.ParquetFile(source).iter_batches(batch_size=131072, columns=COLS):
            d = batch.to_pandas()
            if not (d.donor_id.astype(str) == donor).all():
                raise ValueError(f'donor namespace mismatch {donor}')
            result = derive(d)
            table = pa.Table.from_pandas(result, preserve_index=False)
            if writer is None:
                writer = pq.ParquetWriter(target, table.schema, compression='zstd', use_dictionary=True)
            writer.write_table(table)
            selected = result.loc[result.candidate_primary_changed | result.candidate_refined_changed,
                                  ['cell_id', 'object_id', 'donor_id', 'candidate_primary_changed',
                                   'candidate_refined_changed', 'candidate_identity_reason',
                                   'candidate_refined_identity_reason']]
            at = pa.Table.from_pandas(selected, preserve_index=False)
            if affected_writer is None:
                affected_writer = pq.ParquetWriter(affected, at.schema, compression='zstd', use_dictionary=True)
            affected_writer.write_table(at)
            for key, count in result.groupby(['cell_class', 'refined_class', 'resolved_by',
                                              'candidate_cell_class', 'candidate_identity_reason',
                                              'candidate_refined_class', 'candidate_refined_identity_reason',
                                              'analysis_eligible'], dropna=False).size().items():
                groups[(donor,) + tuple('' if pd.isna(x) else x for x in key)] += int(count)
            n += len(d)
            n_changed += len(selected)
            changed['primary'] += int(result.candidate_primary_changed.sum())
            changed['refined'] += int(result.candidate_refined_changed.sum())
        writer.close()
        affected_writer.close()
        # Independent readback: exact source column values and candidate decisions on every row.
        source_batches = pq.ParquetFile(source).iter_batches(batch_size=131072, columns=COLS)
        output_batches = pq.ParquetFile(target).iter_batches(batch_size=131072)
        n_readback = 0
        for a, b in zip(source_batches, output_batches, strict=True):
            aa = a.to_pandas()
            bb = b.to_pandas()
            pd.testing.assert_frame_equal(aa, bb[COLS], check_dtype=False)
            expected = derive(aa)
            pd.testing.assert_frame_equal(expected, bb, check_dtype=False)
            n_readback += len(bb)
        if n_readback != n or sha(source) != source_sha:
            raise ValueError('incomplete readback or changed immutable input')
        record = {'donor_id': donor, 'source': str(source), 'source_sha256': source_sha,
                  'output': str(target), 'output_sha256': sha(target), 'rows': n,
                  'affected_keys': str(affected), 'affected_keys_sha256': sha(affected),
                  'rows_with_any_corrected_claim': n_changed, 'all_rows_read_back': n_readback,
                  'exact_source_fields_preserved': True}
        inputs.append(record)
        n_all += n
        status['sources'] = inputs
        status['n_rows_completed'] = n_all
        status['peak_RSS_GiB'] = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / (1024 ** 2)
        (OUT / 'manifest.json').write_text(json.dumps(status, indent=2) + '\n')
        print(f'{donor}: {n:,} rows; corrected claims in {n_changed:,}; full readback passed', flush=True)
    columns = ['donor_id', 'source_cell_class', 'source_refined_class', 'source_resolved_by',
               'candidate_cell_class', 'candidate_identity_reason', 'candidate_refined_class',
               'candidate_refined_identity_reason', 'analysis_eligible', 'n_cells']
    census = pd.DataFrame([list(key) + [value] for key, value in groups.items()], columns=columns)
    census.to_csv(OUT / 'exhaustive_route_census.csv', index=False)
    donor_counts = census.groupby(['donor_id', 'candidate_cell_class', 'candidate_identity_reason',
                                  'analysis_eligible'], dropna=False).n_cells.sum().reset_index()
    donor_counts.to_csv(OUT / 'candidate_identity_counts_by_donor.csv', index=False)
    if sha(manifest_path) != original_manifest_sha or len(inputs) != 25 or n_all != 37642209:
        raise ValueError('wrong full pancreas source scope')
    status.update(status='complete_additive_conservative_candidate', n_donors=len(inputs),
                  n_rows=n_all, n_primary_claims_corrected=changed['primary'],
                  n_refined_claims_corrected=changed['refined'], elapsed_seconds=round(time.time()-started, 2),
                  script_sha256=sha(__file__), exhaustive_route_census_sha256=sha(OUT/'exhaustive_route_census.csv'),
                  scope='all current pancreas assignments including QC-ineligible and excluded native rows',
                  limitations=['Unresolved candidate labels abstain from unsupported identities; they do not prove the cells lack epithelial or immune biology.',
                               'Current affirmative producer rules are retained as operational support; independent accuracy remains unestablished.',
                               'Dedicated and embedded pLN candidates are not generated by this pancreas-only route correction.'])
    (OUT / 'manifest.json').write_text(json.dumps(status, indent=2) + '\n')
    print(json.dumps({k: status[k] for k in ['status','n_rows','n_primary_claims_corrected','n_refined_claims_corrected','peak_RSS_GiB','elapsed_seconds']}), flush=True)


if __name__ == '__main__':
    main()
