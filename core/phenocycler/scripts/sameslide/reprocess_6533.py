#!/usr/bin/env python3
"""Snapshot donor 6533 and submit the recorded October InstanSeg recipe.

The upstream scripts are read without modifying them. A local job-script copy adds
the strict duplicate override recorded in the October exports. Submission is
explicit and idempotent: an existing submission receipt must never be overwritten.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import os
from pathlib import Path
import re
import shutil
import subprocess
import sys
import fcntl
from datetime import datetime, timezone

REPO = Path(__file__).resolve().parents[4]
UPSTREAM = Path('/home/smith6jt/IO60panc2nd')
EXPORT = Path('/blue/maigan/smith6jt/IO60panc2nd_v3/results/v3/export')
IMAGE = '0080870_6533_pLNpanc_Scan3.er.qptiff - resolution #1'
REFERENCE = EXPORT / 'provenance/0080502_6580_pLNpanc_Scan1.er.qptiff_-_resolution_1/hpg'
COMMIT = '6d4dbad0379dbc7a10a271c8f675cb1df79d6f0a'


def digest(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as f:
        for chunk in iter(lambda: f.read(8 * 1024 * 1024), b''):
            h.update(chunk)
    return h.hexdigest()


def write_json(path, value):
    tmp = Path(str(path) + '.tmp')
    tmp.write_text(json.dumps(value, indent=2) + '\n')
    tmp.replace(path)


def strict_job_copy(text):
    needle = '--steps "$STEPS" --images "$IMAGE"'
    if text.count(needle) != 1:
        raise ValueError('Upstream job interface changed: expected exactly one run_v3 invocation')
    return text.replace(needle, '--steps "$STEPS" --override "0:DUP_RESOLVE_MAX_FRACTION=0" --images "$IMAGE"')


def prepare(out):
    out.mkdir(parents=True, exist_ok=True)
    project = UPSTREAM / 'project.qpproj'
    entry = [e for e in json.loads(project.read_text())['images'] if e['imageName'] == IMAGE]
    if len(entry) != 1 or str(entry[0]['entryID']) != '37':
        raise ValueError('6533 must resolve uniquely to the recorded project entry 37')
    commit = subprocess.check_output(['git', '-C', str(UPSTREAM / 'scripts'), 'rev-parse', 'HEAD'], text=True).strip()
    if commit != COMMIT:
        raise ValueError(f'Upstream commit differs from October recipe: {commit}')
    sources = {}
    for step in ('preflight', '0', '2', 'C_cells', 'C_geojson'):
        rec = json.loads((REFERENCE / f'{step}.json').read_text())
        source = Path(rec['script']['source_path'])
        actual = digest(source)
        if actual != rec['script']['source_sha256']:
            raise ValueError(f'{step}: source differs from successful October recipe: {source}')
        sources[step] = {'path': str(source), 'sha256': actual,
                         'overrides': rec['script'].get('overrides', {})}
    snapshot = out / 'snapshot'
    snapshot.mkdir(exist_ok=True)
    snapshot_files = {}
    for source in (project, *(UPSTREAM / 'data/37').iterdir()):
        if not source.is_file():
            continue
        target = snapshot / source.name
        if not target.exists():
            temp = Path(str(target) + '.copying')
            before = source.stat()
            shutil.copy2(source, temp)
            after = source.stat()
            if (before.st_size, before.st_mtime_ns) != (after.st_size, after.st_mtime_ns):
                raise ValueError(f'Source changed while snapshotting: {source}')
            if digest(temp) != digest(source):
                raise ValueError(f'Snapshot did not verify: {source}')
            temp.replace(target)
        snapshot_files[source.name] = {'source': str(source), 'snapshot': str(target),
                                       'bytes': target.stat().st_size, 'sha256': digest(target)}
    for kind in ('gpu', 'cpu'):
        source = UPSTREAM / f'scripts/v3/hpg/cells_{kind}.sbatch'
        (out / f'6533_{kind}.sbatch').write_text(strict_job_copy(source.read_text()))
    record = {'created_utc': datetime.now(timezone.utc).isoformat(), 'image': IMAGE,
              'entry_id': '37', 'source_commit': commit, 'reference': str(REFERENCE),
              'sources': sources, 'snapshot': snapshot_files,
              'segmentation_overrides': {'DUP_RESOLVE_MAX_FRACTION': '0'},
              'measurement_threads': 16, 'measurement_heap_gb': 160}
    write_json(out / 'preparation.json', record)
    return record


def submit(out, account='clive'):
    if account not in ('clive', 'maigan'):
        raise ValueError('account must be clive or maigan')
    receipt = out / 'submission.json'
    if not (out / 'preparation.json').is_file():
        raise ValueError('Run preparation before submitting')
    # Prevent a second caller from launching duplicate image jobs.
    lock = out / 'submission.lock'
    with lock.open('a') as handle:
        fcntl.flock(handle, fcntl.LOCK_EX | fcntl.LOCK_NB)
        record = json.loads(receipt.read_text()) if receipt.exists() else {
            'image': IMAGE, 'account': account,
            'submitted_utc': datetime.now(timezone.utc).isoformat()}
        if record.get('account', account) != account:
            raise ValueError('partial submission belongs to another account')
        if record.get('cpu_job'):
            return record
        if record.get('submission_in_progress'):
            raise ValueError('ambiguous interrupted submission; reconcile Slurm before retrying')
        dry = subprocess.check_output(['bash', str(UPSTREAM / 'scripts/v3/hpg/submit_cells.sh'),
             '--account', account, '--image', IMAGE, '--out', str(EXPORT),
             '--cpu-heap', '160', '--cpu-cpus', '20', '--dry-run'], text=True)
        match = re.search(r'-> (hpg-[\w-]+) gpu:([\w_]+)', dry)
        if not match:
            raise ValueError(f'Unable to resolve upstream GPU selection: {dry}')
        partition, gpu = match.groups()
        common = ['sbatch', '--parsable', f'--account={account}', f'--qos={account}']
        env = f'ALL,IMAGE={IMAGE},V3_OUT={EXPORT},GPU_TYPE={gpu},STAGE_TO_TMPDIR=0'
        if not record.get('gpu_job'):
            record['submission_in_progress'] = 'gpu'
            write_json(receipt, record)
            try:
                gpu_id = subprocess.check_output(common + [f'--partition={partition}', f'--gres=gpu:{gpu}:1',
                    '--cpus-per-task=16', '--mem=128gb', '--time=02:15:00', f'--export={env}',
                    str(out / '6533_gpu.sbatch')], text=True).strip().split(';')[0]
            except subprocess.CalledProcessError:
                record.pop('submission_in_progress')
                write_json(receipt, record)
                raise
            record.update(gpu_job=gpu_id, partition=partition, gpu_type=gpu)
            record.pop('submission_in_progress')
            write_json(receipt, record)
        gpu_id = record['gpu_job']
        record['submission_in_progress'] = 'cpu'
        write_json(receipt, record)
        try:
            cpu_id = subprocess.check_output(common + ['--partition=hpg-turin,hpg-milan',
                 '--cpus-per-task=20', '--mem=184gb', '--time=04:45:00',
                 f'--dependency=afterok:{gpu_id}',
                 f'--export=ALL,IMAGE={IMAGE},V3_OUT={EXPORT},STEP2_HEAP_GB=160',
                 str(out / '6533_cpu.sbatch')], text=True).strip().split(';')[0]
        except subprocess.CalledProcessError:
            record.pop('submission_in_progress')
            write_json(receipt, record)
            raise
        record['cpu_job'] = cpu_id
        record.pop('submission_in_progress')
        write_json(receipt, record)
        return record


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--out', type=Path, default=REPO / 'data/development/sameslide_four/upstream_6533')
    ap.add_argument('--submit', action='store_true')
    ap.add_argument('--account', choices=('clive', 'maigan'), default='clive')
    args = ap.parse_args(argv)
    prepared = (json.loads((args.out / 'preparation.json').read_text())
                if (args.out / 'preparation.json').exists() else prepare(args.out))
    print(json.dumps(submit(args.out, args.account) if args.submit else prepared, indent=2))
    return 0


if __name__ == '__main__':
    sys.exit(main())
