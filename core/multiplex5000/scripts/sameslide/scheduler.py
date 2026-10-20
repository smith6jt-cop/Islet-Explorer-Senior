"""Restart-safe Slurm submissions shared by the four-slide development stages."""
import fcntl
import hashlib
import json
import os
from pathlib import Path
import shlex
import subprocess

from workflow import write_json,now,identity,REPO,DEFAULT_OUT

PYTHON=Path('/blue/maigan/smith6jt/miniforge3/envs/xenium_analysis/bin/python')
ACTIVE={'PENDING','RUNNING','CONFIGURING','COMPLETING','SUSPENDED','REQUEUED','RESIZING'}


def job_state(job):
    p=subprocess.run(['sacct','-j',str(job),'--noheader','--parsable2','--format=JobIDRaw,State'],capture_output=True,text=True,check=True)
    rows=[l.split('|') for l in p.stdout.splitlines() if l.split('|')[0]==str(job)]
    if not rows:raise RuntimeError(f'Slurm accounting has no record for {job}; reconciliation required')
    return rows[0][1].split()[0].rstrip('+')


def submit(key,script,args,*,account,mem_gb,cpus=4,hours=4,dependencies=(),root=DEFAULT_OUT):
    if account not in ('maigan','clive'):raise ValueError('only the two selected investment accounts are allowed')
    dest=Path(root)/'jobs'/key;dest.mkdir(parents=True,exist_ok=True)
    script=Path(script).resolve();cmd=[str(PYTHON),str(script),*map(str,args)]
    spec={'command':cmd,'source':identity(script),'account':account,'mem_gb':mem_gb,'cpus':cpus}
    fingerprint=hashlib.sha256(json.dumps(spec,sort_keys=True).encode()).hexdigest()
    with (dest/'submission.lock').open('a') as lock:
        fcntl.flock(lock,fcntl.LOCK_EX)
        path=dest/'current.json';previous=json.loads(path.read_text()) if path.exists() else None
        if previous and previous.get('job_id'):
            state=job_state(previous['job_id'])
            if state in ACTIVE:
                if previous['fingerprint']!=fingerprint:raise RuntimeError(f'{key}: different code is already running as {previous["job_id"]}')
                return previous['job_id']
            if state=='COMPLETED' and previous['fingerprint']==fingerprint:
                return previous['job_id']  # Worker output still requires its own validated receipt.
        elif previous and previous.get('status')=='submitting':
            raise RuntimeError(f'{key}: interrupted submission; reconcile the recorded unique job name before retry')
        attempt=dest/(fingerprint[:12]+'-'+now().replace(':','').replace('+','_'))
        attempt.mkdir()
        jobname='ss4-'+hashlib.sha256(str(attempt).encode()).hexdigest()[:12]
        dep=[str(j) for j in dependencies if job_state(j) in ACTIVE]
        body=['#!/bin/bash',f'#SBATCH --job-name={jobname}',f'#SBATCH --account={account}',f'#SBATCH --qos={account}',
              '#SBATCH --partition=hpg-default',f'#SBATCH --cpus-per-task={cpus}',f'#SBATCH --mem={mem_gb}G',
              f'#SBATCH --time={hours}:00:00',f'#SBATCH --output={attempt}/job_%j.log','set -euo pipefail',
              'cd '+shlex.quote(str(REPO)),f'export OMP_NUM_THREADS={cpus}',f'export OPENBLAS_NUM_THREADS={cpus}',
              f'export MKL_NUM_THREADS={cpus}',
              'export PYTHONPATH='+shlex.quote(str(Path(root)/'python_deps')),
              'exec '+shlex.join(cmd)]
        sbatch=attempt/'job.sbatch';sbatch.write_text('\n'.join(body)+'\n')
        rec={'status':'submitting','specification':spec,'fingerprint':fingerprint,'job_name':jobname,
             'attempt':str(attempt),'dependencies':dep,'started_utc':now()}
        write_json(path,rec);write_json(attempt/'submission.json',rec)
        command=['sbatch','--parsable']
        if dep:command+=['--dependency=afterok:'+':'.join(dep)]
        result=subprocess.run(command+[str(sbatch)],capture_output=True,text=True)
        rec['stdout']=result.stdout;rec['stderr']=result.stderr
        if result.returncode:
            rec['status']='rejected';write_json(path,rec);write_json(attempt/'submission.json',rec)
            raise RuntimeError(f'{key}: sbatch rejected: {result.stderr}')
        token=result.stdout.strip().split(';')[0]
        if not token.isdigit():raise RuntimeError(f'{key}: ambiguous sbatch response; reconcile {jobname}')
        rec.update(status='submitted',job_id=int(token),submitted_utc=now())
        write_json(path,rec);write_json(attempt/'submission.json',rec)
        return int(token)
