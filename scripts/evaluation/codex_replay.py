#!/usr/bin/env python3
"""Content-bound CODEX replay with separate producer inputs and exact references.

The pinned port remains external. Its numerical assertions and scientific code
are unchanged. Only the validation fixture's provenance and reference loading
are repaired; source files and prior evaluation artifacts are never written.
"""
from __future__ import annotations

import argparse
import hashlib
import importlib.metadata
import importlib.util
import json
import os
import platform
import sys
import time
import traceback
from pathlib import Path


class ContractError(ValueError):
    """The declared replay inputs, implementation, or runtime are incompatible."""


def sha256(path: Path) -> str:
    with Path(path).open('rb') as stream:
        return hashlib.file_digest(stream, 'sha256').hexdigest()


def file_identity(path: Path) -> dict:
    st = Path(path).stat()
    return {'size': st.st_size, 'mtime_ns': st.st_mtime_ns, 'ctime_ns': st.st_ctime_ns, 'inode': st.st_ino}


def require_hash(path: Path, expected: str):
    actual = sha256(path)
    if actual != expected:
        raise ContractError(f'Content binding mismatch: {path}: {actual} != {expected}')
    return actual


def verify_provenance(proof_path: Path, contract: dict) -> tuple[dict, dict]:
    """Verify the anchored producer proof and every consumed content binding."""
    require_hash(proof_path, contract['proof_sha256'])
    proof = json.loads(proof_path.read_text())
    if proof['status'] != 'producer_outputs_and_execution_inputs_verified':
        raise ContractError('Producer proof is incomplete')
    expected = contract['images']
    if len(expected) != len(set(expected)):
        raise ContractError('Duplicate acquisition in replay contract')
    for key in ['tags', 'cohort', 'project', 'identity']:
        if proof['roster'][key] != expected:
            raise ContractError(f'Incomplete or inconsistent acquisition roster: {key}')
    if proof['source_pin'] != contract['source_pin'] or proof['port_pin'] != contract['port_pin']:
        raise ContractError('Producer or port pin disagrees with the replay contract')
    require_hash(Path(proof['journal']), proof['journal_sha256'])
    require_hash(Path(proof['structures_journal']), proof['structures_journal_sha256'])
    identity = {}
    verified = []
    cache = {}
    for group in ['code', 'outputs', 'inputs']:
        for row in proof[group]:
            path = Path(row['path'])
            if row['binding_kind'] == 'sha256':
                expected_hash = row['expected']
            elif row.get('historical_sha256') is None and 'original_stat' in row:
                current = file_identity(path)
                if current != {k: row['original_stat'][k] for k in current}:
                    raise ContractError(f'Original established-call file identity changed: {path}')
                expected_hash = row['sha256_now']
            else:
                raise ContractError(f'Unsupported provenance binding: {row["binding_kind"]}')
            if str(path) not in cache:
                cache[str(path)] = require_hash(path, expected_hash)
            elif cache[str(path)] != expected_hash:
                raise ContractError(f'Conflicting content bindings: {path}')
            identity[str(path)] = file_identity(path)
            verified.append({'path': str(path), 'sha256': expected_hash, 'group': group})
    if len(proof['outputs']) != contract['producer_outputs']:
        raise ContractError('Incomplete producer output closure')
    return proof, {'verified': verified, 'file_identities': identity,
                   'proof_sha256': contract['proof_sha256'],
                   'historical_limit': proof['binding_limitation']}


def check_unchanged(identities: dict):
    changed = [p for p, expected in identities.items() if file_identity(Path(p)) != expected]
    if changed:
        raise ContractError(f'Inputs changed during replay: {changed}')


def verify_runtime(contract: dict) -> dict:
    expected = contract['runtime']
    if platform.python_version() != expected['python']:
        raise ContractError(f"Python {platform.python_version()} != required {expected['python']}")
    packages = {name: importlib.metadata.version(name) for name in expected['packages']}
    if packages != expected['packages']:
        wrong = {name: value for name, value in packages.items() if value != expected['packages'][name]}
        raise ContractError(f'Runtime package versions differ: {wrong}')
    for name, value in expected['thread_environment'].items():
        if os.environ.get(name) != value:
            raise ContractError(f'{name} must be {value} before importing numerical libraries')
    # Load the libraries before inspecting actual pools; an environment variable
    # alone cannot constrain a library already loaded with another thread count.
    import numpy  # noqa: F401
    import scipy.linalg  # noqa: F401
    import sklearn  # noqa: F401
    from threadpoolctl import threadpool_info
    pools = threadpool_info()
    if any(pool['num_threads'] != 1 for pool in pools):
        raise ContractError(f'Numerical thread pools must all use one thread: {pools}')
    return {'python': platform.python_version(), 'executable': sys.executable,
            'packages': packages, 'thread_environment': expected['thread_environment'], 'thread_pools': pools,
            'scope': 'Verified replay runtime, not a claim that the historical producer supplied a package lock.'}


def load_table_views(run: Path) -> tuple[dict, dict]:
    """Preserve producer parsing for inputs; decode expected CSV floats losslessly."""
    import pandas as pd
    names = {'levels': 'levels.csv', 'levels_differentiation': 'levels_differentiation.csv',
             'selection': 'levels_selection.csv', 'differentiation_selection': 'levels_differentiation_selection.csv',
             'canary_levels': 'canary_levels.csv'}
    inputs = {key: pd.read_csv(run / 'typing' / name) for key, name in names.items()}
    references = {key: pd.read_csv(run / 'typing' / name, float_precision='round_trip')
                  for key, name in names.items()}
    return inputs, references


def assertion_views(inputs: dict, references: dict) -> dict:
    # Levels uses selections and empty-null source inputs, while canary consumes
    # producer-parsed levels. A single shared DataFrame cannot serve both the
    # producer's input-parsing contract and a lossless expected-output contract.
    return {'test_levels_are_reproduced': {**inputs, 'levels': references['levels'],
                'levels_differentiation': references['levels_differentiation']},
            'test_canary_levels_are_reproduced': {**inputs, 'canary_levels': references['canary_levels']},
            'test_identity_is_reproduced': inputs}


def import_assertions(snapshot: Path, proof: dict):
    fixture = snapshot / 'core/multiplex5000/tests/current/test_positivity_parity.py'
    bound = {Path(row['path']).resolve(): row['expected'] for row in proof['code']}
    if fixture.resolve() not in bound:
        raise ContractError('Selected port snapshot is not the content-verified snapshot')
    require_hash(fixture, bound[fixture.resolve()])
    sys.path.insert(0, str(snapshot / 'core/multiplex5000'))
    spec = importlib.util.spec_from_file_location('pinned_codex_replay_assertions', fixture)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    for imported in [module.LV, module.CN, module.ID, module.RG, module.sources]:
        path = Path(imported.__file__).resolve()
        if path not in bound:
            raise ContractError(f'An unpinned scientific module was imported: {path}')
        require_hash(path, bound[path])
    return module


def replay(contract_path: Path, evidence: Path, snapshot: Path, out: Path) -> int:
    sys.dont_write_bytecode = True
    contract = json.loads(contract_path.read_text())
    # This entry point imports only stdlib before enforcing the declared pool settings.
    for key, value in contract['runtime']['thread_environment'].items():
        os.environ[key] = value
    runtime = verify_runtime(contract)
    out.mkdir(parents=True, exist_ok=False)
    record = {'schema': 'codex.replay/1', 'status': 'validating_provenance', 'contract': str(contract_path.resolve()),
              'contract_sha256': sha256(contract_path), 'runtime': runtime,
              'original_fixture_unchanged': True, 'criterion': 'check_exact=True; check_dtype=False; array_equal',
              'expected_images': len(contract['images']), 'images': []}

    def save():
        (out / 'replay.json').write_text(json.dumps(record, indent=2) + '\n')

    save()
    try:
        proof, verified = verify_provenance(evidence / 'journal_proof.json', contract)
        (out / 'verified_provenance.json').write_text(json.dumps(verified, indent=2) + '\n')
        journal = json.loads(Path(proof['journal']).read_text())
        run = Path(journal['run'])
        assertions = import_assertions(snapshot, proof)
        inputs, references = load_table_views(run)
        views = assertion_views(inputs, references)
        record['run'] = str(run)
        record['status'] = 'running'
        import pandas as pd
        for tag in contract['images']:
            start = time.monotonic()
            checks = []
            for name, tables in views.items():
                try:
                    getattr(assertions, name)(run, tag, tables)
                    checks.append({'check': name, 'status': 'pass'})
                except Exception as exc:  # noqa: BLE001 - report every acquisition's actual failed check.
                    checks.append({'check': name, 'status': 'fail', 'error': str(exc), 'traceback': traceback.format_exc()})
            ids = pd.read_parquet(run / f'typing/identity/{tag}.parquet', columns=['object_id'])
            rs = pd.read_parquet(run / f'typing/redsea/{tag}.parquet', columns=['object_id', 'qc_analysis_eligible'])
            eligible = rs.loc[rs['qc_analysis_eligible'], 'object_id']
            if ids['object_id'].duplicated().any() or set(ids['object_id']) != set(eligible):
                raise ContractError(f'Cell identity/count reconciliation failed: {tag}')
            record['images'].append({'image': tag, 'checks': checks, 'identity_cells': len(ids),
                'measured_objects': len(rs), 'eligible_cells': len(eligible),
                'elapsed_seconds': round(time.monotonic() - start, 3)})
            save()
            print(tag, [r['status'] for r in checks], flush=True)
        check_unchanged(verified['file_identities'])
        record.update(status='complete', passed=sum(c['status'] == 'pass' for r in record['images'] for c in r['checks']),
                      failed=sum(c['status'] == 'fail' for r in record['images'] for c in r['checks']),
                      unchanged_scientific_inputs_and_baselines=True)
        save()
        return int(record['failed'] != 0)
    except Exception as exc:
        record.update(status='failed', error=str(exc), traceback=traceback.format_exc())
        save()
        raise


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--contract', type=Path, default=Path(__file__).with_name('codex_replay_contract.json'))
    parser.add_argument('--evidence', type=Path, required=True)
    parser.add_argument('--snapshot', type=Path, required=True)
    parser.add_argument('--out', type=Path, required=True)
    args = parser.parse_args()
    raise SystemExit(replay(args.contract, args.evidence, args.snapshot, args.out))
