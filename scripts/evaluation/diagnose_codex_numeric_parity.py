#!/usr/bin/env python3
"""Separate exact source-to-port comparison and CSV precision diagnostics.

No comparison tolerance is added. The pinned artifact assertions keep their
failures; direct source-to-port assertions use the same exact criterion.
"""
from __future__ import annotations

import importlib
import itertools
import json
import os
import sys
import traceback
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / 'outputs/handoff_completion_20261102/codex'
SNAPSHOT = ROOT / 'outputs/dynamic_portability_20261102/typing_port'
ORIGINAL = OUT / 'producer_records/source_578767e/scripts/uniform/typing'


def main():
    sys.dont_write_bytecode = True
    threads = os.environ.get('CODEX_PARITY_THREADS', '4')
    os.environ['OPENBLAS_NUM_THREADS'] = threads
    os.environ['OMP_NUM_THREADS'] = threads
    os.environ['PYTHONDONTWRITEBYTECODE'] = '1'
    import numpy as np
    import pandas as pd

    report_name = os.environ.get('CODEX_PARITY_DIAGNOSTIC_NAME', 'source_to_port_diagnostics')
    if Path(report_name).name != report_name:
        raise ValueError('Diagnostic name must be a plain filename stem')
    tables_only = os.environ.get('CODEX_PARITY_TABLES_ONLY') == '1'

    proof = json.loads((OUT / 'journal_proof.json').read_text())
    assert proof['status'] == 'producer_outputs_and_execution_inputs_verified'
    journal = json.loads(Path(proof['journal']).read_text())
    run = Path(journal['run'])
    sys.path.insert(0, str(ORIGINAL))
    old_lv, old_cn, old_id = [importlib.import_module(n) for n in ['levels', 'canary', 'identity']]
    sys.path.insert(0, str(SNAPSHOT / 'core/multiplex5000'))
    from multiplex5000.positivity import canary as new_cn
    from multiplex5000.positivity import identity as new_id
    from multiplex5000.positivity import levels as new_lv
    from multiplex5000.positivity import registry, sources

    names = {'levels': 'levels.csv', 'levels_differentiation': 'levels_differentiation.csv',
             'selection': 'levels_selection.csv', 'differentiation_selection': 'levels_differentiation_selection.csv',
             'canary_levels': 'canary_levels.csv'}
    tables = {key: pd.read_csv(run / 'typing' / name) for key, name in names.items()}
    roundtrip = {key: pd.read_csv(run / 'typing' / name, float_precision='round_trip') for key, name in names.items()}
    diff = tables['differentiation_selection']
    selected = tables['selection']
    chosen = {r.marker: (r.measurement, r.method) for r in itertools.chain(
        selected[~selected['marker'].isin(diff['marker'])].itertuples(), diff.itertuples())}
    positivity = {r.marker: (r.measurement, r.method) for r in selected[selected['marker'].isin(diff['marker'])].itertuples()}
    reg = registry.load('codex-uniform-32', 'pancreas')
    source = sources.CodexRunSource(run)
    original_nulls = old_lv.empty_nulls(run)
    port_nulls = new_lv.empty_nulls(source, reg)
    assert original_nulls == port_nulls

    def same(got, expected, key):
        got = got[list(expected.columns)].sort_values(key).reset_index(drop=True)
        expected = expected.sort_values(key).reset_index(drop=True)
        try:
            pd.testing.assert_frame_equal(got, expected, check_exact=True, check_dtype=False)
            return {'status': 'pass'}
        except AssertionError as e:
            return {'status': 'fail', 'error': str(e)}

    def difference(got, reference, key):
        got = got[list(reference.columns)].sort_values(key).reset_index(drop=True)
        reference = reference.sort_values(key).reset_index(drop=True)
        columns = []
        for col in reference:
            a, b = got[col], reference[col]
            missing_disagreements = int((a.isna() != b.isna()).sum())
            if pd.api.types.is_numeric_dtype(a) and pd.api.types.is_numeric_dtype(b):
                av, bv = a.to_numpy(dtype=float), b.to_numpy(dtype=float)
                valid = np.isfinite(av) & np.isfinite(bv)
                delta = np.abs(av[valid] - bv[valid])
                n = int((delta != 0).sum())
                if n or missing_disagreements:
                    columns.append({'column': col, 'unequal_finite': n,
                        'missing_disagreements': missing_disagreements,
                        'max_absolute_difference': float(delta.max()) if len(delta) else None,
                        'max_relative_difference': float(np.max(delta / np.maximum(np.abs(bv[valid]), np.finfo(float).tiny))) if len(delta) else None})
            else:
                n = int((~((a == b) | (a.isna() & b.isna()))).sum())
                if n or missing_disagreements:
                    columns.append({'column': col, 'unequal': n, 'missing_disagreements': missing_disagreements})
        return columns

    results = []
    for tag in proof['roster']['tags']:
        record = {'image': tag, 'direct_source_to_port': [], 'artifact_diagnostics': []}
        try:
            original = pd.DataFrame(old_lv.image_pass(run, tag, ('codex-uniform-32', 'pancreas'),
                'final', chosen, original_nulls, positivity))
            port = pd.DataFrame(new_lv.image_pass(source, tag, reg.spec, 'final', chosen, port_nulls, positivity))
            record['direct_source_to_port'].append({'computation': 'levels', **same(port, original, ['reading', 'marker'])})
            for reading, name in [('positivity', 'levels'), ('differentiation', 'levels_differentiation')]:
                got = port[port['reading'] == reading]
                native = original[original['reading'] == reading]
                ref = tables[name][tables[name]['image'] == tag]
                rt = roundtrip[name][roundtrip[name]['image'] == tag]
                record['artifact_diagnostics'].append({'table': name,
                    'port_vs_default_csv': difference(got, ref, ['marker']),
                    'original_vs_default_csv': difference(native, ref, ['marker']),
                    'port_vs_roundtrip_csv': difference(got, rt, ['marker'])})
            original_canary = pd.DataFrame(old_cn.image_levels(run, tag))
            port_canary = pd.DataFrame(new_cn.image_levels(source, tag, reg.spec, tables))
            record['direct_source_to_port'].append({'computation': 'canary', **same(port_canary, original_canary, ['marker'])})
            ref = tables['canary_levels'][tables['canary_levels']['image'] == tag]
            rt = roundtrip['canary_levels'][roundtrip['canary_levels']['image'] == tag]
            record['artifact_diagnostics'].append({'table': 'canary_levels',
                'port_vs_default_csv': difference(port_canary, ref, ['marker']),
                'original_vs_default_csv': difference(original_canary, ref, ['marker']),
                'port_vs_roundtrip_csv': difference(port_canary, rt, ['marker'])})
            if not tables_only:
                original_identity, _ = old_id.classify_image(run, tag)
                port_identity, _ = new_id.classify_image(source, tag, reg.spec, tables)
                record['direct_source_to_port'].append({'computation': 'identity', **same(port_identity, original_identity, ['object_id'])})
        except Exception as e:  # noqa: BLE001 - keep every acquisition's diagnostic failure visible.
            record['error'] = str(e)
            record['traceback'] = traceback.format_exc()
        results.append(record)
        summary = {'schema': 'codex.exact_source_port_and_csv_diagnostics/1', 'results': results,
            'images_completed': len(results), 'expected_images': 20,
            'criterion': 'Exact equality; no tolerance added or criterion relaxed.',
            'computations': ['levels', 'canary'] if tables_only else ['levels', 'canary', 'identity'],
            'threads': threads,
            'status': 'complete' if len(results) == 20 else 'running',
            'note': 'Round-trip CSV parsing is a diagnostic only. The original fixture was not edited.'}
        (OUT / f'{report_name}.json').write_text(json.dumps(summary, indent=2) + '\n')
        print(tag, [(t['computation'], t['status']) for t in record['direct_source_to_port']], record.get('error', ''), flush=True)
    return int(any(r.get('error') or any(t['status'] != 'pass' for t in r['direct_source_to_port']) for r in results))


if __name__ == '__main__':
    raise SystemExit(main())
