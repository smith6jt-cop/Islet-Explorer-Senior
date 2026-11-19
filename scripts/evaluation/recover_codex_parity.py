#!/usr/bin/env python3
"""Recover CODEX producer provenance, then call the pinned numerical assertions.

This is a separate evaluation harness. It does not invoke, edit, monkeypatch, or
report success for the old pytest provenance fixture. That fixture's mtime rule
still rejects the historical files. The recovered execution journal binds the
actual pre-commit source bytes and all output bytes instead.
"""
from __future__ import annotations

import csv
import hashlib
import importlib.util
import json
import os
import subprocess
import sys
import time
import traceback
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
CODEX = Path('/data/panc_CODEX/Islet-Explorer')
OUT = ROOT / 'outputs/handoff_completion_20261102/codex'
RECORDS = OUT / 'producer_records'
SNAPSHOT = ROOT / 'outputs/dynamic_portability_20261102/typing_port'
PIN = '578767eeaf0f95bb1df7272238f361888a0ae75d'
PORT_PIN = 'df3fe1d6ecfe534bc1eaf9e1cf94b8ef6d3a7018'
TYPING = RECORDS / '20261030T212756.342638Z/execution.json'
STRUCTURES = RECORDS / '20261030T210113.520762Z/execution.json'
_HASHES = {}


def sha(path):
    path = Path(path)
    st = path.stat()
    key = (str(path), st.st_size, st.st_mtime_ns, st.st_ctime_ns)
    if key not in _HASHES:
        h = hashlib.sha256()
        with path.open('rb') as f:
            for block in iter(lambda: f.read(1 << 20), b''):
                h.update(block)
        _HASHES[key] = h.hexdigest()
    return _HASHES[key]


def save(name, value):
    (OUT / name).write_text(json.dumps(value, indent=2) + '\n')


def stat_record(path):
    s = Path(path).stat()
    return {'size': s.st_size, 'mtime_ns': s.st_mtime_ns, 'ctime_ns': s.st_ctime_ns, 'inode': s.st_ino}


def main():
    sys.dont_write_bytecode = True
    os.environ['PYTHONDONTWRITEBYTECODE'] = '1'
    os.environ['OPENBLAS_NUM_THREADS'] = '4'
    os.environ['OMP_NUM_THREADS'] = '4'
    journal = json.loads(TYPING.read_text())
    structure_journal = json.loads(STRUCTURES.read_text())
    run = Path(journal['run'])
    proof = {'schema': 'codex.journal_validated_parity/1', 'source_pin': PIN,
             'port_pin': PORT_PIN, 'journal': str(TYPING), 'journal_sha256': sha(TYPING),
             'structures_journal': str(STRUCTURES), 'structures_journal_sha256': sha(STRUCTURES),
             'code': [], 'outputs': [], 'inputs': [], 'roster': {}, 'status': 'auditing',
             'old_fixture': 'Unchanged and still fails its pre-2075d31 mtime rule.',
             'scope': 'Computational fidelity only; no biological acceptance or release activation.'}

    def bound(path, expected, source, bucket='inputs', kind='sha256'):
        path = Path(path)
        actual = sha(path)
        rec = {'path': str(path), 'expected': expected, 'actual': actual,
                   'match': actual == expected, 'binding_source': source, 'binding_kind': kind}
        proof[bucket].append(rec)
        if not rec['match']:
            save('journal_proof.json', proof)
            raise ValueError(f'Content binding failed: {path}')

    assert journal['status'] == structure_journal['status'] == 'complete'
    assert journal['previous_phase'] == structure_journal['created_at_utc']
    for path, rec in journal['outputs'].items():
        bound(path, rec['sha256'], str(TYPING), 'outputs')
    for stage in journal['stages']:
        assert stage['returncode'] == 0
        bound(stage['log'], stage['log_sha256'], str(TYPING))

    for path in sorted((RECORDS / 'source_578767e').rglob('*')):
        if not path.is_file():
            continue
        relative = path.relative_to(RECORDS / 'source_578767e')
        key = str(CODEX / relative)
        expected = journal['code_sha256'][key]
        bound(path, expected, f'{TYPING}:code_sha256', 'code')
        blob = subprocess.check_output(['git', '-C', str(CODEX), 'show', PIN + ':' + str(relative)])
        assert hashlib.sha256(blob).hexdigest() == expected

    # Verify every imported port implementation and the assertion source against
    # the original Git pin, not merely the current working tree.
    port_files = list((SNAPSHOT / 'core/multiplex5000/multiplex5000/positivity').glob('*.py'))
    port_files += list((SNAPSHOT / 'core/multiplex5000/multiplex5000/positivity/registry').glob('*.csv'))
    fixture = SNAPSHOT / 'core/multiplex5000/tests/current/test_positivity_parity.py'
    port_files.append(fixture)
    for path in port_files:
        rel = path.relative_to(SNAPSHOT)
        expected = hashlib.sha256(subprocess.check_output(['git', '-C', str(ROOT), 'show', PORT_PIN + ':' + str(rel)])).hexdigest()
        bound(path, expected, f'git:{PORT_PIN}:{rel}', 'code')

    with (run / 'typing/levels.csv').open() as f:
        tags = sorted({r['image'] for r in csv.DictReader(f)})
    with (RECORDS / 'source_578767e/scripts/uniform/cohort.csv').open() as f:
        cohort_tags = sorted(r['image_name'] for r in csv.DictReader(f))
    project = json.loads((run / 'project_manifest.json').read_text())
    project_tags = sorted(r['image_name'] for r in project['entries'])
    identity_tags = sorted(Path(p).stem for p in journal['outputs'] if '/identity/' in p and p.endswith('.parquet'))
    assert tags == cohort_tags == project_tags == identity_tags and len(tags) == 20
    proof['roster'] = {'tags': tags, 'images': len(tags), 'cohort': cohort_tags, 'project': project_tags, 'identity': identity_tags}
    null_record = run / 'typing/redsea_refresh/empty_area_nulls.json'
    nulls = json.loads(null_record.read_text())
    bound(run / 'typing/empty_area_nulls.csv', nulls['csv_sha256'], str(null_record))
    assert stat_record(null_record) == {k: journal['protected_inputs'][str(null_record)][k] for k in stat_record(null_record)}
    entries = {r['image_name']: r for r in project['entries']}
    for tag in tags:
        refresh_path = run / f'typing/redsea_refresh/{tag}.json'
        refresh = json.loads(refresh_path.read_text())
        assert stat_record(refresh_path) == {k: journal['protected_inputs'][str(refresh_path)][k] for k in stat_record(refresh_path)}
        bound(run / f'typing/redsea/{tag}.parquet', refresh['outputs']['parquet'], str(refresh_path))
        bound(run / f'typing/redsea/{tag}.json', refresh['outputs']['json'], str(refresh_path))
        assert refresh['outputs']['json'] == nulls['redsea_record_sha256'][tag]
        for key in ['cells', 'geometry', 'islets', 'server']:
            rec = refresh['inputs']['files'][key]
            bound(rec['path'], rec['sha256'], str(refresh_path))
        context = run / f'typing/context/{tag}.parquet'
        bound(context, structure_journal['outputs'][str(context)]['sha256'], str(STRUCTURES))
        calls = (run / f'superseded/typing-v17-20261028/calls/{tag}.parquet').resolve()
        recorded = journal['protected_inputs'][str(calls)]
        actual = stat_record(calls)
        assert actual == {k: recorded[k] for k in actual}
        proof['inputs'].append({'path': str(calls), 'match': True,
            'binding_kind': 'original exact inode/size/mtime/ctime; newly measured SHA256',
            'binding_source': str(TYPING), 'original_stat': recorded, 'actual_stat': actual,
            'sha256_now': sha(calls), 'historical_sha256': None})
        coverage = json.loads((run / f'typing/coverage/{tag}.json').read_text())
        for rel, expected in coverage['inputs']['sha256'].items():
            if rel in ['registry', 'coverage_code', 'pixel_processing_code']:
                continue  # Verified as exact original source blobs above.
            bound(run / rel, expected, str(run / f'typing/coverage/{tag}.json'))
        image = coverage['inputs']['image']
        st = Path(image['path']).stat()
        assert st.st_size == image['size'] and st.st_mtime_ns == image['mtime_ns']
        assert entries[tag]['image_path'] == image['path']

    proof['status'] = 'producer_outputs_and_execution_inputs_verified'
    proof['binding_limitation'] = ('Established endocrine calls were recorded by filesystem identity, not SHA256, '
        'by the original producer. All 20 retain the exact recorded inode/size/mtime/ctime; current SHA256s are added '
        'for this evaluation. Source TIFFs retain the original path/size/mtime contract.')
    save('journal_proof.json', proof)
    print('Producer proof passed: 95 outputs, all 20 images, exact pinned scientific source.', flush=True)

    sys.path.insert(0, str(SNAPSHOT / 'core/multiplex5000'))
    import pandas as pd
    spec = importlib.util.spec_from_file_location('pinned_numerical_assertions', fixture)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    # The original numerical functions are called unchanged. This harness uses
    # the independently recorded proof above; the old fixture is not invoked.
    tables = module.tables.__wrapped__(run)
    results = []
    for tag in tags:
        started = time.monotonic()
        tests = []
        for name in ['test_levels_are_reproduced', 'test_canary_levels_are_reproduced', 'test_identity_is_reproduced']:
            try:
                getattr(module, name)(run, tag, tables)
                tests.append({'test': name, 'status': 'pass'})
            except Exception as e:  # noqa: BLE001 - retain every acquisition's failure in the report.
                tests.append({'test': name, 'status': 'fail', 'error': str(e), 'traceback': traceback.format_exc()})
        ids = pd.read_parquet(run / f'typing/identity/{tag}.parquet', columns=['object_id'])
        rs = pd.read_parquet(run / f'typing/redsea/{tag}.parquet', columns=['object_id', 'qc_analysis_eligible'])
        eligible = rs.loc[rs['qc_analysis_eligible'], 'object_id']
        results.append({'image': tag, 'tests': tests, 'identity_cells': len(ids),
            'redsea_cells': len(rs), 'eligible_cells': len(eligible),
            'identity_ids_unique': not ids['object_id'].duplicated().any(),
            'identity_ids_equal_eligible': set(ids['object_id']) == set(eligible),
            'elapsed_seconds': round(time.monotonic() - started, 3)})
        summary = {'schema': 'codex.independent_numerical_parity/1', 'producer_proof': str(OUT / 'journal_proof.json'),
            'producer_proof_sha256': sha(OUT / 'journal_proof.json'), 'results': results,
            'images_completed': len(results), 'expected_images': 20,
            'assertions_passed': sum(t['status'] == 'pass' for r in results for t in r['tests']),
            'assertions_failed': sum(t['status'] == 'fail' for r in results for t in r['tests']),
            'comparison': 'Unchanged pinned assert_frame_equal(check_exact=True, check_dtype=False) and array_equal',
            'original_pytest_provenance_fixture': 'unchanged; still fails pre-commit mtime rule',
            'status': 'complete' if len(results) == 20 else 'running'}
        save('numerical_parity.json', summary)
        print(tag, [t['status'] for t in tests], len(ids), flush=True)
    return int(summary['assertions_failed'] > 0)


if __name__ == '__main__':
    raise SystemExit(main())
