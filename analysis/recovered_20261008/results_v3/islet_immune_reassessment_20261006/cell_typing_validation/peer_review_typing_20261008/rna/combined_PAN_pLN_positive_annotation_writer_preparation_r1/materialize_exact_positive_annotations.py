#!/usr/bin/env python3
"""Materialize an explicitly enumerated, source-bound RNA review contract.

This writer has no marker, depth, reference-label or abundance decision rule.
It applies only root-adjudicated exact native-QC memberships. It preserves all
other values, source-relative flag semantics and the original Arrow schemas.
Selection, source measurements, notebooks and figures are outside its writes.
"""
from __future__ import annotations
if not __debug__:
    raise RuntimeError('This materializer requires enabled assertions; Python -O is unsupported.')

import argparse
from collections import Counter, defaultdict
from datetime import datetime, timezone
import gc
import hashlib
import json
from pathlib import Path
import shutil
import traceback

import numpy as np
import pandas as pd
import pyarrow as pa
import pyarrow.parquet as pq

HERE = Path(__file__).resolve().parent
SCHEMA = 'exact_reviewed_shared_broad_positive_RNA_candidate_materialization_v2'
# Existing shared-taxonomy pairs only. A root-accepted exact contract chooses
# the label; this writer neither scores the evidence nor invents subtypes.
BROAD_PAIRS = {
    'Exocrine': 'Exocrine', 'Endocrine': 'Endocrine', 'Vascular': 'Vascular',
    'Myeloid': 'Myeloid', 'Granulocyte': 'Granulocyte', 'Neural': 'Neural',
    'Stromal_Fibroblast_Stellate': 'Stromal', 'Stromal_Pericyte_SMC': 'Stromal',
    'Immune_T': 'T_NK', 'Immune_B_Plasma': 'B_Plasma',
}
PLN_COMMON = {'Immune_T': 'T_NK', 'Immune_B_Plasma': 'B_Plasma',
              'Myeloid': 'Myeloid', 'Granulocyte': 'Granulocyte', 'Non_immune': ''}


def pLN_operations(label, reason):
    """Exact root/A-reviewed broad-pLN field semantics; no subtype algorithm."""
    assert label in PLN_COMMON and isinstance(reason, str) and reason.strip()
    named = label != 'Non_immune'
    values = {
        'candidate_cluster_lineage': label, 'candidate_lineage': label,
        'candidate_common_lineage': PLN_COMMON[label],
        'candidate_common_lineage_analysis_eligible': named,
        'candidate_reason': reason,
        'candidate_subtype_descriptor': 'qualified_broad_cluster_programme_annotation_no_per_cell_subtype_gate' if named else 'qualified_coarse_Non_immune_structural_programme_no_named_admission',
        'candidate_exact_unit_decision': 'qualified_broad_cluster_annotation_original_QC_members' if named else 'qualified_coarse_Non_immune_structural_annotation_no_named_admission',
        'reference_evidence_status': 'literal_complete_profile_qualified_broad_annotation_reference_diagnostic_only',
    }
    operations = {k: {'operation': 'set_exact_value', 'value': v} for k, v in values.items()}
    for flag, source, candidate in [
        ('candidate_stage2_assignment_changed', 'cluster_lineage', 'candidate_cluster_lineage'),
        ('candidate_assignment_changed', 'lineage', 'candidate_lineage'),
        ('candidate_common_lineage_changed', 'common_lineage', 'candidate_common_lineage'),
    ]:
        operations[flag] = {'operation': 'recompute_source_relative_change_flag',
                            'source_column': source, 'candidate_column': candidate}
    return operations
KEYS = {'cell_id', 'serial_qualified_cell_key', 'subcluster_id', 'cluster_id',
        'section_id', 'source_serial'}
FLAG_PAIRS = {
    'candidate_assignment_changed': ('lineage', 'candidate_lineage'),
    'candidate_stage1_assignment_changed': ('stage1', 'candidate_stage1'),
    'candidate_stage2_assignment_changed': ('cluster_lineage', 'candidate_cluster_lineage'),
    'candidate_common_lineage_changed': ('common_lineage', 'candidate_common_lineage'),
}


def sha(path):
    with Path(path).open('rb') as handle:
        return hashlib.file_digest(handle, 'sha256').hexdigest()


def bind(path):
    path = Path(path).resolve()
    return {'path': str(path), 'sha256': sha(path)}


def verify(binding):
    assert sha(binding['path']) == binding['sha256'], binding['path']


def save(path, data):
    with Path(path).open('x') as handle:
        json.dump(data, handle, indent=2, allow_nan=False)
        handle.write('\n')


def same(a, b):
    return (a.eq(b) | (a.isna() & b.isna())).fillna(False).to_numpy(dtype=bool)


def exact(a, b):
    pd.testing.assert_frame_equal(a.reset_index(drop=True), b.reset_index(drop=True), check_exact=True)


def key_digest(values):
    return hashlib.sha256(('\n'.join(sorted(values.astype(str))) + '\n').encode()).hexdigest()


def distribution(series):
    if pd.api.types.is_bool_dtype(series.dtype):
        assert series.notna().all()
        return Counter({'true': int(series.sum()), 'false': int((~series).sum())})
    return Counter({str(k): int(v) for k, v in series.fillna('__NULL__').astype(str).value_counts().items()})


def normalized(values):
    return {k: int(v) for k, v in values.items() if v}


def distributions(table, qc):
    return {c: distribution(table.loc[qc, c]) for c in table if c not in KEYS}


def update_counts(target, values, factor=1):
    for col, counts in values.items():
        for value, count in counts.items():
            target[col][value] += factor * int(count)
            assert target[col][value] >= 0, (col, value)


def eligible(tissue):
    return 'candidate_lineage_analysis_eligible' if tissue == 'pancreas' else 'candidate_common_lineage_analysis_eligible'


def claim(tissue):
    return 'candidate_lineage' if tissue == 'pancreas' else 'candidate_cluster_lineage'


BASELINE_SPEC = {'count_field': 'independent_r12_r9_candidate_sidecar_nonkey_field_counts_json', 'generation_field': 'independent_r12_r9_candidate_generation', 'path_field': 'independent_r12_r9_candidate_sidecar_path', 'sha256_field': 'independent_r12_r9_candidate_sidecar_sha256', 'QC_count_field': 'r12_r9_n_original_QC_cells'}


def load_bound_baseline_ledger(path, spec, parents):
    """Read the explicitly bound current generation, preserving historical fields.

    The r12/r9 candidate ledger retains historical selected-generation fields.
    Only its explicitly independently read r12/r9 distributions supply this
    materialization's baseline; the candidate need not have been selected.
    """
    assert spec == BASELINE_SPEC, 'Baseline fields must name the explicit candidate-generation evidence.'
    fixed = ['partition_key', 'tissue', 'section_id', 'source_serial', 'n_cells']
    ledger = pd.read_parquet(path, columns=fixed + list(spec.values()))
    assert ledger.partition_key.is_unique
    sections = {}
    for tissue, parent in parents.items():
        for section in parent['sections']:
            key = (tissue, section['section_id'], section['source_serial'])
            assert key not in sections
            sections[key] = section
    for row in ledger.to_dict('records'):
        tissue = row['tissue']
        section = sections[(tissue, row['section_id'], row['source_serial'])]
        assert row[spec['generation_field']] == parents[tissue]['candidate_generation']
        assert row[spec['path_field']] == section['candidate_sidecar_path']
        assert row[spec['sha256_field']] == section['candidate_sidecar_sha256']
        assert row[spec['QC_count_field']] == row['n_cells']
        values = json.loads(row[spec['count_field']])
        assert values and all(sum(v.values()) == row['n_cells'] for v in values.values())
    return ledger[fixed + [spec['count_field']]].rename(columns={
        spec['count_field']: 'materialization_baseline_field_counts_json'})


def read_source(section, tissue):
    columns = ['cell_id', 'section_key', 'lineage', 'common_lineage']
    columns += ['subcluster_id', 'source_serial', 'serial_qualified_cell_key'] if tissue == 'pancreas' else ['cluster_id', 'stage1', 'cluster_lineage']
    source = pd.read_parquet(section['source_calls_path'], columns=columns)
    if tissue == 'pLN':
        source = source.loc[source.section_key.eq(section['section_id'])].reset_index(drop=True)
    return source


def materialize_section(section, tissue, cases, destination, baseline, global_after):
    verify({'path': section['source_calls_path'], 'sha256': section['source_calls_sha256']})
    pf = pq.ParquetFile(section['candidate_sidecar_path'])
    before = pf.read().to_pandas()
    source = read_source(section, tissue)
    cc = 'subcluster_id' if tissue == 'pancreas' else 'cluster_id'
    assert len(before) == len(source) == section.get('n_source_cells', section.get('n_source_rows'))
    assert before.cell_id.is_unique and before.serial_qualified_cell_key.is_unique
    assert before.cell_id.eq(source.cell_id).all() and same(before[cc], source[cc]).all()
    assert source.section_key.eq(section['section_id']).all()
    prefix = section['source_serial'] + '|' + section['section_id'] + '|' if tissue == 'pancreas' else 'pLN|' + section['section_id'] + '|' + section['source_serial'] + '|'
    expected_keys = prefix + source.cell_id
    assert before.serial_qualified_cell_key.eq(expected_keys).all()
    if tissue == 'pancreas':
        assert source.source_serial.eq(section['source_serial']).all() and source.serial_qualified_cell_key.eq(expected_keys).all()
    else:
        assert before.section_id.eq(section['section_id']).all() and before.source_serial.eq(section['source_serial']).all()
    qc = source[cc].notna().to_numpy()
    assert int(qc.sum()) == section['n_native_QC_cells']
    before_counts = distributions(before, qc)
    assert set(before_counts) == set(baseline)
    assert all(normalized(before_counts[c]) == normalized(baseline[c]) for c in before_counts)
    revised = before.copy()
    allowed = {c: np.zeros(len(before), dtype=bool) for c in before}
    masks = {}
    for case in cases:
        tissue_key, sid, serial, cluster = case['partition_key'].split('|')
        assert (tissue_key, sid, serial) == (tissue, section['section_id'], section['source_serial'])
        mask = before[cc].eq(cluster).to_numpy()
        assert qc[mask].all() and int(mask.sum()) == case['n_original_QC_cells']
        assert key_digest(before.loc[mask, 'serial_qualified_cell_key']) == case['original_QC_key_digest']
        assert before.loc[mask, claim(tissue)].eq(case['current_partition_claim']).all()
        assert int(before.loc[mask, eligible(tissue)].sum()) == case['n_prior_admissions']
        operations = case['field_operations']
        assert operations and not (set(operations) & KEYS)
        for col, op in operations.items():
            assert col in before
            allowed[col] |= mask
            if op['operation'] == 'set_exact_value':
                assert col not in FLAG_PAIRS, 'Source-relative flags must be recomputed explicitly.'
                if pa.types.is_boolean(pf.schema_arrow.field(col).type):
                    assert type(op['value']) is bool
                else:
                    assert isinstance(op['value'], str)
                revised.loc[mask, col] = op['value']
            else:
                assert op['operation'] == 'recompute_source_relative_change_flag' and col in FLAG_PAIRS
                assert (op['source_column'], op['candidate_column']) == FLAG_PAIRS[col]
        # Derived flags use the fully applied candidate values, independent of
        # operation ordering and without a constant-true assumption.
        for col, op in operations.items():
            if op['operation'] == 'recompute_source_relative_change_flag':
                sc, candidate = FLAG_PAIRS[col]
                revised.loc[mask, col] = ~same(source.loc[mask, sc], revised.loc[mask, candidate])
        assert int(revised.loc[mask, eligible(tissue)].sum()) == case['n_new_admissions']
        masks[case['partition_key']] = mask
    path = destination / 'calls' / ('donor_id=' + section['donor_id']) / ('section=' + section['section_id']) / 'candidate_sidecar.parquet'
    path.parent.mkdir(parents=True, exist_ok=False)
    table = pa.Table.from_pandas(revised, schema=pf.schema_arrow, preserve_index=False).replace_schema_metadata(pf.schema_arrow.metadata)
    pq.write_table(table, path, compression='zstd')
    del table, revised
    out_pf = pq.ParquetFile(path)
    assert out_pf.schema_arrow.equals(pf.schema_arrow, check_metadata=True)
    actual = out_pf.read().to_pandas()
    assert list(actual) == list(before)
    case_receipts = []
    for case in cases:
        mask = masks[case['partition_key']]
        for col, op in case['field_operations'].items():
            if op['operation'] == 'set_exact_value':
                assert actual.loc[mask, col].eq(op['value']).all()
            else:
                sc, candidate = FLAG_PAIRS[col]
                assert np.array_equal(actual.loc[mask, col].to_numpy(), ~same(source.loc[mask, sc], actual.loc[mask, candidate]))
        fixed = [c for c in before if c not in case['field_operations']]
        exact(before.loc[mask, fixed], actual.loc[mask, fixed])
        case_receipts.append({'partition_key': case['partition_key'], 'original_QC_cells': int(mask.sum()),
            'prior_admissions': int(before.loc[mask, eligible(tissue)].sum()), 'new_admissions': int(actual.loc[mask, eligible(tissue)].sum()),
            'original_QC_key_digest': key_digest(actual.loc[mask, 'serial_qualified_cell_key']),
            'changed_values_by_field': {c: int((~same(before.loc[mask, c], actual.loc[mask, c])).sum()) for c in before}})
    deltas = []
    for col in before:
        changed = ~same(before[col], actual[col])
        assert not changed[~allowed[col]].any(), (section['section_id'], col)
        deltas.append({'tissue': tissue, 'section_id': section['section_id'], 'source_serial': section['source_serial'],
            'field': col, 'n_values_changed': int(changed.sum()), 'changes_outside_accepted_operations': 0})
    for flag, (sc, candidate) in FLAG_PAIRS.items():
        if flag in actual:
            assert np.array_equal(actual[flag].to_numpy(), ~same(source[sc], actual[candidate])), (section['section_id'], flag)
    assert not actual.loc[~qc, eligible(tissue)].any()
    actual_counts = distributions(actual, qc)
    update_counts(global_after, before_counts, -1)
    update_counts(global_after, actual_counts, 1)
    all_bool = {c: int(actual[c].sum()) for c in actual if pa.types.is_boolean(pf.schema_arrow.field(c).type)}
    receipt = {'tissue': tissue, 'section_id': section['section_id'], 'source_serial': section['source_serial'],
        'prior_sidecar': bind(section['candidate_sidecar_path']), 'new_sidecar': bind(path), 'source_calls': bind(section['source_calls_path']),
        'n_source_rows': len(actual), 'n_original_QC_rows': int(qc.sum()), 'n_fields': len(actual.columns),
        'n_actual_output_values_checked': len(actual) * len(actual.columns), 'all_field_deltas': deltas,
        'all_original_arrow_schema_metadata_preserved': True, 'all_native_keys_order_QC_and_source_relative_flags_verified': True,
        'every_value_outside_exact_accepted_operations_preserved': True, 'cases': case_receipts,
        'QC_field_counts': {c: normalized(v) for c, v in actual_counts.items()}, 'all_source_boolean_counts': all_bool}
    stable = ['donor_id', 'section_id', 'source_serial', 'source_calls_path', 'source_calls_sha256', 'source_matrix_path', 'source_matrix_sha256', 'n_source_cells', 'n_source_rows', 'n_native_QC_cells']
    state = {k: section[k] for k in stable if k in section}
    state.update({'candidate_sidecar_path': str(path), 'candidate_sidecar_sha256': sha(path),
        'sidecar_storage': 'fresh_exact_root_adjudicated_operations', 'parent_sidecar': bind(section['candidate_sidecar_path']),
        'n_candidate_common_lineage_analysis_eligible_QC_cells': int(actual.loc[qc, eligible(tissue)].sum()),
        'all_source_boolean_counts': all_bool, 'all_rows_all_fields_producer_readback': True})
    if tissue == 'pancreas':
        state.update({'n_candidate_lineage_analysis_eligible_QC_cells': state['n_candidate_common_lineage_analysis_eligible_QC_cells'],
            'n_candidate_known_lineage_QC_admissions': state['n_candidate_common_lineage_analysis_eligible_QC_cells'],
            'n_candidate_assignments_changed': all_bool['candidate_assignment_changed'],
            'n_candidate_assignments_withheld': int((actual.candidate_assignment_changed & ~actual.candidate_lineage_analysis_eligible & qc).sum())})
    else:
        state['n_candidate_common_lineage_analysis_eligible'] = state['n_candidate_common_lineage_analysis_eligible_QC_cells']
        for name, flag in [('stage1', 'candidate_stage1_assignment_changed'), ('stage2', 'candidate_stage2_assignment_changed'), ('stage3', 'candidate_assignment_changed')]:
            state['n_' + name + '_assignments_changed'] = all_bool[flag]
    print('fresh_section_PASS', tissue, section['section_id'], len(actual), len(cases), flush=True)
    del before, actual, source
    gc.collect()
    return state, receipt


def preflight_decision_counts(parent, tissue, cases, expected_uncertainty,
                              expected_annotations, expected_superseded):
    """Reject a decision-history mismatch before any candidate namespace exists."""
    assert type(expected_uncertainty) is int and expected_uncertainty >= 0
    assert type(expected_annotations) is int and expected_annotations >= 0
    # This reviewed positive extension has no accepted operative supersessions.
    assert type(expected_superseded) is int and expected_superseded == 0
    stem = 'uncertainty_decisions' if tissue == 'pancreas' else 'decisions'
    verify({'path': parent[stem + '_path'], 'sha256': parent[stem + '_sha256']})
    prior = pd.read_csv(parent[stem + '_path'], dtype=str, keep_default_na=False)
    assert 'partition_key' in prior and prior.partition_key.is_unique
    accepted = {c['partition_key'] for c in cases}
    assert len(accepted) == len(cases)
    assert not set(prior.partition_key) & accepted, 'This batch cannot supersede an operative uncertainty decision.'
    assert len(prior) == expected_uncertainty
    prior_positive = parent.get('positive_annotation_decisions')
    n_prior_positive = 0
    if prior_positive:
        verify(prior_positive)
        positive = pd.read_csv(prior_positive['path'], dtype=str, keep_default_na=False)
        assert 'partition_key' in positive and positive.partition_key.is_unique
        assert not set(positive.partition_key) & accepted
        n_prior_positive = len(positive)
    assert n_prior_positive + len(cases) == expected_annotations


def decisions(parent, tissue, cases, destination, acceptance):
    stem = 'uncertainty_decisions' if tissue == 'pancreas' else 'decisions'
    verify({'path': parent[stem + '_path'], 'sha256': parent[stem + '_sha256']})
    prior = pd.read_csv(parent[stem + '_path'], dtype=str, keep_default_na=False)
    assert 'partition_key' in prior, 'The r12/r9 parent must carry exact operative keys.'
    assert prior.partition_key.is_unique
    accepted = {c['partition_key'] for c in cases}
    prior_positive_binding = parent.get('positive_annotation_decisions')
    if prior_positive_binding:
        verify(prior_positive_binding)
        prior_positive = pd.read_csv(prior_positive_binding['path'], dtype=str, keep_default_na=False)
        assert prior_positive.partition_key.is_unique
        # Reannotation of an existing positive operation requires a separately
        # designed contract; this positive-only extension never supersedes one.
        assert not set(prior_positive.partition_key) & accepted
        shutil.copy2(prior_positive_binding['path'], destination / 'historical_predecessor_positive_annotation_decisions.csv')
    else:
        prior_positive = None
    fixed = ~prior.partition_key.isin(accepted)
    combined = prior.loc[fixed].reset_index(drop=True)
    # A superseded withholding is historical evidence, not an operative
    # uncertainty decision in the new positively annotated generation.
    shutil.copy2(parent[stem + '_path'], destination / 'historical_predecessor_uncertainty_decisions.csv')
    path = destination / 'serial_qualified_uncertainty_decisions.csv'
    if fixed.all():
        # In particular, all fifteen pLN operative uncertainty rows are
        # disjoint from accepted scopes and remain byte-identical.
        shutil.copy2(parent[stem + '_path'], path)
    else:
        combined.to_csv(path, index=False)
    exact(combined, pd.read_csv(path, dtype=str, keep_default_na=False))
    additions = []
    for case in cases:
        key = case['partition_key']
        _, section, serial, cluster = key.split('|')
        old = prior.loc[prior.partition_key.eq(key)].to_dict('records')
        additions.append({'partition_key': key, 'section_id': section, 'source_serial': serial,
            'cluster_id': cluster, 'n_cells': str(case['n_original_QC_cells']),
            'prior_lineage': case['current_partition_claim'],
            'candidate_lineage': case['proposed_candidate_lineage'] if tissue == 'pancreas' else case['qualified_broad_annotation'],
            'reason': case['scientific_reason'], 'accepted_scope': case['accepted_scope'],
            'root_acceptance_path': acceptance['path'], 'root_acceptance_sha256': acceptance['sha256'],
            'n_prior_admissions': str(case['n_prior_admissions']),
            'n_new_admissions': str(case['n_new_admissions']),
            'original_QC_key_digest': case['original_QC_key_digest'],
            'field_operations_json': json.dumps(case['field_operations']),
            'superseded_withholding_records_json': json.dumps(old)})
    annotation_path = destination / 'serial_qualified_positive_annotation_decisions.csv'
    if additions or prior_positive is not None:
        annotation = pd.DataFrame(additions)
        if prior_positive is not None:
            if additions:
                assert list(annotation) == list(prior_positive), 'Prior positive schema must be retained.'
                annotation = pd.concat([prior_positive, annotation], ignore_index=True)
            else:
                annotation = prior_positive.copy()
        assert annotation.partition_key.is_unique
        annotation.to_csv(annotation_path, index=False)
        exact(annotation, pd.read_csv(annotation_path, dtype=str, keep_default_na=False))
        if prior_positive is not None:
            exact(prior_positive, annotation.iloc[:len(prior_positive)])
    else:
        annotation_path = None
        annotation = pd.DataFrame()
    return (path, len(combined), int(combined.n_cells.astype(int).sum()), annotation_path,
            len(prior) - len(combined), 0 if prior_positive is None else len(prior_positive),
            len(annotation), 0 if annotation.empty else int(annotation.n_cells.astype(int).sum()))


def validate_case(case):
    parts = case['partition_key'].split('|')
    assert len(parts) == 4 and all(parts) and parts[0] in {'pancreas', 'pLN'}
    assert type(case['n_original_QC_cells']) is int and case['n_original_QC_cells'] > 0
    assert type(case['n_prior_admissions']) is int and case['n_prior_admissions'] == 0
    assert len(case['original_QC_key_digest']) == 64 and all(c in '0123456789abcdef' for c in case['original_QC_key_digest'])
    assert isinstance(case['scientific_reason'], str) and case['scientific_reason'].strip()
    if parts[0] == 'pLN':
        assert case['case_kind'] == 'exact_qualified_broad_pLN_partition_annotation'
        label = case['qualified_broad_annotation']
        assert label in PLN_COMMON
        assert type(case['n_new_admissions']) is int
        assert case['n_new_admissions'] == (0 if label == 'Non_immune' else case['n_original_QC_cells'])
        assert case['accepted_scope'] == 'qualified_broad_native_partition_annotation_not_exclusive_individual_cell_truth'
        assert case['field_operations'] == pLN_operations(label, case['scientific_reason'])
        return
    assert case['case_kind'] == 'exact_qualified_broad_positive_annotation'
    assert type(case['n_new_admissions']) is int and case['n_new_admissions'] == case['n_original_QC_cells']
    assert case['scientific_reason'].strip() and case['replacement_identity_assigned'] is True
    lineage = case['proposed_candidate_lineage']
    assert lineage in BROAD_PAIRS
    assert case['proposed_candidate_common_lineage'] == BROAD_PAIRS[lineage]
    assert case['accepted_scope'] == 'qualified_broad_' + lineage + '_cluster_annotation_original_QC_members'
    ops = case['field_operations']
    expected = {
        'candidate_lineage': {'operation': 'set_exact_value', 'value': lineage},
        'candidate_common_lineage': {'operation': 'set_exact_value', 'value': BROAD_PAIRS[lineage]},
        'candidate_reason': {'operation': 'set_exact_value', 'value': case['scientific_reason']},
        'candidate_lineage_analysis_eligible': {'operation': 'set_exact_value', 'value': True},
        'candidate_assignment_changed': {'operation': 'recompute_source_relative_change_flag',
            'source_column': 'lineage', 'candidate_column': 'candidate_lineage'},
        'reference_evidence_status': {'operation': 'set_exact_value',
            'value': 'root_adjudicated_qualified_broad_' + lineage + '_annotation_reference_diagnostic_only'},
    }
    assert ops == expected


def main(acceptance_path, digest):
    acceptance = {'path': str(acceptance_path), 'sha256': digest}
    verify(acceptance)
    contract = json.loads(acceptance_path.read_text())
    assert contract['schema_version'] == SCHEMA and contract['candidate_materialization_authorized'] is True
    assert contract['selected_registry_change_authorized_by_this_document'] is False
    for b in [contract['selected_registry'], contract['baseline_7052_ledger'], *contract['parent_releases'].values(),
              contract['parent_independent_all52_readback'], *contract['scientific_bindings']]:
        verify(b)
    cases = contract['cases']
    assert cases and len({c['partition_key'] for c in cases}) == len(cases)
    for case in cases:
        validate_case(case)
    parents = {t: json.loads(Path(b['path']).read_text()) for t, b in contract['parent_releases'].items()}
    destinations = {t: Path(d) for t, d in contract['destination_directories'].items()}
    run = Path(contract['materialization_run_directory'])
    assert set(parents) == set(destinations) == {'pancreas', 'pLN'}
    for case in cases:
        tissue, sid, serial, _ = case['partition_key'].split('|')
        assert sum(s['section_id'] == sid and s['source_serial'] == serial for s in parents[tissue]['sections']) == 1, case['partition_key']
    assert not run.exists() and all(not d.exists() for d in destinations.values())
    ledger = load_bound_baseline_ledger(contract['baseline_7052_ledger']['path'],
        contract['ledger_baseline_field_semantics'], parents)
    assert len(ledger) == 7052 and ledger.partition_key.is_unique and int(ledger.n_cells.sum()) == 25255757
    global_before = {t: defaultdict(Counter) for t in parents}
    section_before = defaultdict(lambda: defaultdict(Counter))
    for row in ledger.itertuples(index=False):
        values = json.loads(row.materialization_baseline_field_counts_json)
        assert all(sum(v.values()) == row.n_cells for v in values.values())
        update_counts(global_before[row.tissue], values)
        update_counts(section_before[(row.tissue, row.section_id, row.source_serial)], values)
    global_after = {t: defaultdict(Counter, {c: Counter(v) for c, v in global_before[t].items()}) for t in parents}
    protected_sidecars = []
    protected_decisions = []
    for tissue, parent in parents.items():
        lane_cases = [c for c in cases if c['partition_key'].split('|')[0] == tissue]
        preflight_decision_counts(parent, tissue, lane_cases,
            contract['expected_lane_operative_uncertainty_decisions'][tissue],
            contract['expected_lane_operative_annotation_decisions'][tissue],
            contract['expected_lane_superseded_uncertainty_decisions'][tissue])
        stem = 'uncertainty_decisions' if tissue == 'pancreas' else 'decisions'
        decision_binding = {'path': parent[stem + '_path'], 'sha256': parent[stem + '_sha256']}
        verify(decision_binding)
        protected_decisions.append(decision_binding)
        if parent.get('positive_annotation_decisions'):
            verify(parent['positive_annotation_decisions'])
            protected_decisions.append(parent['positive_annotation_decisions'])
        for section in parent['sections']:
            b = {'path': section['candidate_sidecar_path'], 'sha256': section['candidate_sidecar_sha256']}
            verify(b)
            protected_sidecars.append(b)
    run.mkdir(parents=True)
    shutil.copy2(Path(__file__), run / 'executed_producer.py')
    save(run / 'frozen_acceptance_binding.json', acceptance)
    all_receipts = []
    releases = {}
    for tissue, parent in parents.items():
        destination = destinations[tissue]
        destination.mkdir(parents=True)
        lane_cases = [c for c in cases if c['partition_key'].split('|')[0] == tissue]
        states = []
        fresh = []
        for section in parent['sections']:
            chosen = [c for c in lane_cases if c['partition_key'].split('|')[1:3] == [section['section_id'], section['source_serial']]]
            if chosen:
                state, receipt = materialize_section(section, tissue, chosen, destination,
                    section_before[(tissue, section['section_id'], section['source_serial'])], global_after[tissue])
                fresh.append(receipt)
                all_receipts.append(receipt)
            else:
                state = dict(section)
                state['sidecar_storage'] = 'inherited_exact_immutable_parent_file'
                state['all_rows_all_fields_fresh_readback_by_this_producer'] = False
            states.append(state)
        assert sum(len(r['cases']) for r in fresh) == len(lane_cases)
        (decision_path, n_decisions, n_decision_cells, annotation_path, n_superseded,
         n_prior_positive, n_operative_positive, n_positive_QC) = decisions(parent, tissue, lane_cases, destination, acceptance)
        assert n_decisions == contract['expected_lane_operative_uncertainty_decisions'][tissue]
        assert n_operative_positive == contract['expected_lane_operative_annotation_decisions'][tissue]
        assert n_superseded == contract['expected_lane_superseded_uncertainty_decisions'][tissue] == 0
        after = {c: normalized(v) for c, v in global_after[tissue].items()}
        admissions = int(global_after[tissue][eligible(tissue)].get('true', 0))
        prior_admissions = int(global_before[tissue][eligible(tissue)].get('true', 0))
        assert prior_admissions - admissions == sum(c['n_prior_admissions'] - c['n_new_admissions'] for c in lane_cases)
        assert admissions == contract['expected_lane_admissions'][tissue]
        assert sum(after[claim(tissue)].values()) == parent['n_native_QC_cells']
        release = {'release_status': 'UNSELECTED_EXACT_SCIENTIFIC_CORRECTIONS_PRODUCER_VERIFIED_PENDING_INDEPENDENT_READBACK',
            'candidate_generation': contract['candidate_generations'][tissue], 'tissue': tissue,
            'script_path': str(Path(__file__)), 'script_sha256': sha(__file__),
            'parent_release': contract['parent_releases'][tissue], 'parent_independent_all52_readback': contract['parent_independent_all52_readback'],
            'count_baseline_ledger': contract['baseline_7052_ledger'],
            'count_baseline_field_semantics': contract['ledger_baseline_field_semantics'],
            'n_sections': len(states), 'n_source_cells': parent['n_source_cells'], 'n_native_QC_cells': parent['n_native_QC_cells'],
            'n_exact_native_partitions': int(ledger.tissue.eq(tissue).sum()),
            'n_native_partitions': int(ledger.tissue.eq(tissue).sum()),
            'n_new_changed_section_sidecars': len(fresh), 'n_inherited_unchanged_section_sidecars': len(states) - len(fresh),
            'n_candidate_common_lineage_analysis_eligible_QC_cells': admissions,
            'candidate_QC_lineage_counts': after[claim(tissue)], 'candidate_QC_all_nonkey_field_counts': after,
            'candidate_QC_all_boolean_field_totals': {c: int(v.get('true', 0)) for c, v in after.items() if set(v) <= {'true', 'false'}},
            'n_operative_exact_uncertainty_decisions': n_decisions, 'n_QC_cells_in_operative_exact_decisions': n_decision_cells,
            'n_withheld_partitions': n_decisions, 'n_withheld_cells': n_decision_cells,
            'n_new_exact_unit_uncertainty_decisions': 0,
            'n_new_named_admissions_withheld': 0,
            'n_new_qualified_broad_annotation_units': len(lane_cases),
            'n_new_qualified_broad_annotation_QC_cells': sum(c['n_original_QC_cells'] for c in lane_cases),
            'n_new_qualified_named_annotation_units': sum(c['n_new_admissions'] > 0 for c in lane_cases),
            'n_new_qualified_named_annotation_QC_admissions': sum(c['n_new_admissions'] for c in lane_cases),
            'n_new_qualified_coarse_nonimmune_annotation_units': sum(c.get('qualified_broad_annotation') == 'Non_immune' for c in lane_cases),
            'n_new_qualified_coarse_nonimmune_annotation_QC_cells': sum(c['n_original_QC_cells'] for c in lane_cases if c.get('qualified_broad_annotation') == 'Non_immune'),
            'n_superseded_operative_withholding_decisions': n_superseded,
            'positive_annotation_decisions': bind(annotation_path) if annotation_path else None,
            'n_inherited_operative_positive_annotation_decisions': n_prior_positive,
            'n_operative_positive_annotation_decisions': n_operative_positive,
            'n_QC_cells_in_operative_positive_annotation_decisions': n_positive_QC,
            'n_existing_reason_only_corrections': 0,
            'count_semantics': 'Peer exact-scope decisions, identity-uncertainty flags and all unadmitted QC observations are distinct. Lane/view counts are nonadditive.',
            'root_scientific_scope_acceptance': acceptance, 'source_evidence': contract['scientific_bindings'],
            'sections': states, 'decisions_path': str(decision_path), 'decisions_sha256': sha(decision_path),
            'uncertainty_decisions_path': str(decision_path), 'uncertainty_decisions_sha256': sha(decision_path),
            'source_QC_membership_measurements_geometry_native_records_or_parent_files_modified': False,
            'selected_registry_notebooks_figures_or_production_alias_updated': False,
            'replacement_identity_assigned': any(c['n_new_admissions'] > 0 for c in lane_cases),
            'qualified_annotation_fields_changed': bool(lane_cases), 'independent_materialization_verification_required_before_selection': True,
            'scientific_identity_accuracy_purity_function_or_molecule_origin_established': False,
            'limits': contract['limits']}
        if tissue == 'pancreas':
            release.update({'n_candidate_known_lineage_QC_admissions': admissions,
                'n_candidate_lineage_analysis_eligible_QC_cells': admissions,
                'baseline_release_path': parent['baseline_release_path'], 'baseline_release_sha256': parent['baseline_release_sha256']})
        else:
            release.update({'native_release_path': parent['native_release_path'], 'native_release_sha256': parent['native_release_sha256'],
                'candidate_QC_stage2_partition_claim_counts': after['candidate_cluster_lineage'],
                'candidate_QC_stage1_context_counts': after['candidate_stage1'],
                'candidate_specific_NK_lineage_admitted_cells': after['candidate_specific_NK_lineage_analysis_eligible'].get('true', 0),
                'candidate_stable_Treg_lineage_admitted_cells': after['candidate_stable_Treg_lineage_analysis_eligible'].get('true', 0)})
        save(destination / 'release.json', release)
        releases[tissue] = bind(destination / 'release.json')
    consumed = [c['partition_key'] for r in all_receipts for c in r['cases']]
    assert len(consumed) == len(cases) and set(consumed) == {c['partition_key'] for c in cases}
    save(run / 'all_changed_sections_producer_full_row_field_readbacks.json', {'status': 'PRODUCER_ALL_CHANGED_ROWS_AND_FIELDS_PASS', 'sections': all_receipts})
    for b in [acceptance, contract['selected_registry'], contract['baseline_7052_ledger'], *contract['parent_releases'].values(),
              contract['parent_independent_all52_readback'], *contract['scientific_bindings'], *protected_sidecars, *protected_decisions]:
        verify(b)
    save(run / 'terminal_unselected_candidate_materialization.json', {
        'status': 'TERMINAL_PRODUCER_EXACT_REVIEWED_CANDIDATE_MATERIALIZATION_PASS_INDEPENDENT_REVIEW_REQUIRED',
        'completed_UTC': datetime.now(timezone.utc).isoformat(), 'acceptance': acceptance,
        'producer': bind(__file__), 'releases': releases,
        'full_row_field_readbacks': bind(run / 'all_changed_sections_producer_full_row_field_readbacks.json'),
        'n_fresh_changed_sections': len(all_receipts), 'n_fresh_source_rows_checked': sum(r['n_source_rows'] for r in all_receipts),
        'n_fresh_output_values_checked': sum(r['n_actual_output_values_checked'] for r in all_receipts),
        'n_exact_operation_units': len(cases), 'all52_parent_sidecars_rehashed_and_unchanged': True,
        'all_parent_decision_CSVs_rehashed_and_unchanged': True, 'all_contract_cases_consumed_exactly_once': True,
        'all_source_native_QC_keys_geometry_and_measurements_unchanged': True,
        'selected_registry_modified': False, 'biological_identity_accuracy_validated': False})
    print(json.dumps({'releases': releases, 'terminal': bind(run / 'terminal_unselected_candidate_materialization.json')}), flush=True)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--root-acceptance', type=Path, required=True)
    parser.add_argument('--root-acceptance-sha256', required=True)
    args = parser.parse_args()
    try:
        main(args.root_acceptance.resolve(), args.root_acceptance_sha256)
    except Exception:
        # The failed namespace is preserved. There is no implicit cleanup or
        # permission to overwrite/restart an existing materialization.
        print(traceback.format_exc(), flush=True)
        raise
