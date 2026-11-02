"""Validation-fixture failures must remain loud without changing scientific code."""
from __future__ import annotations

import hashlib
import importlib.util
import json
import os
from pathlib import Path

import pytest

SPEC = importlib.util.spec_from_file_location('codex_replay', Path(__file__).parents[1] / 'codex_replay.py')
REPLAY = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(REPLAY)


def fixture_proof(tmp_path):
    baseline = tmp_path / 'baseline.csv'
    baseline.write_text('image,value\nS1,1.0\n')
    journal = tmp_path / 'producer.json'
    journal.write_text('{"status":"complete"}\n')
    digest = hashlib.sha256(baseline.read_bytes()).hexdigest()
    proof = {'status': 'producer_outputs_and_execution_inputs_verified',
        'source_pin': 'original', 'port_pin': 'port', 'journal': str(journal),
        'journal_sha256': REPLAY.sha256(journal), 'structures_journal': str(journal),
        'structures_journal_sha256': REPLAY.sha256(journal),
        'roster': {key: ['S1'] for key in ['tags', 'cohort', 'project', 'identity']},
        'code': [], 'inputs': [],
        'outputs': [{'path': str(baseline), 'expected': digest, 'binding_kind': 'sha256'}],
        'binding_limitation': 'Synthetic fixture'}
    path = tmp_path / 'proof.json'
    path.write_text(json.dumps(proof))
    contract = {'proof_sha256': REPLAY.sha256(path), 'images': ['S1'], 'producer_outputs': 1,
                'source_pin': 'original', 'port_pin': 'port'}
    return path, contract, baseline


def test_old_timestamp_cannot_override_valid_producer_content(tmp_path):
    proof, contract, baseline = fixture_proof(tmp_path)
    os.utime(baseline, (1, 1))
    _, verified = REPLAY.verify_provenance(proof, contract)
    assert len(verified['verified']) == 1


def test_mutated_baseline_cannot_pass_with_new_timestamp(tmp_path):
    proof, contract, baseline = fixture_proof(tmp_path)
    baseline.write_text('image,value\nS1,2.0\n')
    with pytest.raises(REPLAY.ContractError, match='Content binding mismatch'):
        REPLAY.verify_provenance(proof, contract)


def test_proof_itself_must_match_reviewed_anchor(tmp_path):
    proof, contract, _ = fixture_proof(tmp_path)
    proof.write_text(proof.read_text() + ' ')
    with pytest.raises(REPLAY.ContractError, match='Content binding mismatch'):
        REPLAY.verify_provenance(proof, contract)


def test_incomplete_acquisition_roster_is_rejected(tmp_path):
    proof, contract, _ = fixture_proof(tmp_path)
    contract['images'].append('S2')
    with pytest.raises(REPLAY.ContractError, match='acquisition roster'):
        REPLAY.verify_provenance(proof, contract)


def test_duplicate_acquisition_roster_is_rejected(tmp_path):
    proof, contract, _ = fixture_proof(tmp_path)
    contract['images'].append('S1')
    with pytest.raises(REPLAY.ContractError, match='Duplicate acquisition'):
        REPLAY.verify_provenance(proof, contract)


def test_wrong_runtime_refused_before_computation():
    with pytest.raises(REPLAY.ContractError, match='Python'):
        REPLAY.verify_runtime({'runtime': {'python': '0.0.0'}})


def test_unsafe_thread_settings_are_refused(monkeypatch):
    monkeypatch.setenv('OPENBLAS_NUM_THREADS', '4')
    with pytest.raises(REPLAY.ContractError, match='OPENBLAS_NUM_THREADS must be 1'):
        REPLAY.verify_runtime({'runtime': {'python': REPLAY.platform.python_version(), 'packages': {},
            'thread_environment': {'OPENBLAS_NUM_THREADS': '1'}}})


def test_lossless_references_do_not_change_producer_input_parsing(tmp_path):
    typing = tmp_path / 'typing'
    typing.mkdir()
    for name in ['levels', 'levels_differentiation', 'levels_selection',
                 'levels_differentiation_selection', 'canary_levels']:
        (typing / f'{name}.csv').write_text('image,marker,value\nS1,M,0.123456789123456789\n')
    inputs, references = REPLAY.load_table_views(tmp_path)
    assert inputs['levels'].iloc[0]['value'] != references['levels'].iloc[0]['value']
    assert references['levels'].iloc[0]['value'] == float('0.123456789123456789')
    views = REPLAY.assertion_views(inputs, references)
    assert views['test_levels_are_reproduced']['levels'] is references['levels']
    assert views['test_canary_levels_are_reproduced']['levels'] is inputs['levels']
    assert views['test_canary_levels_are_reproduced']['canary_levels'] is references['canary_levels']
    assert views['test_identity_is_reproduced'] is inputs


def test_inputs_changed_during_replay_fail(tmp_path):
    proof, contract, baseline = fixture_proof(tmp_path)
    _, verified = REPLAY.verify_provenance(proof, contract)
    baseline.write_text('changed')
    with pytest.raises(REPLAY.ContractError, match='Inputs changed during replay'):
        REPLAY.check_unchanged(verified['file_identities'])


def test_real_contract_declares_the_complete_twenty_image_cohort():
    contract = json.loads((Path(__file__).parents[1] / 'codex_replay_contract.json').read_text())
    assert len(contract['images']) == len(set(contract['images'])) == 20
    assert contract['producer_outputs'] == 95
    assert set(contract['runtime']['thread_environment'].values()) == {'1'}
