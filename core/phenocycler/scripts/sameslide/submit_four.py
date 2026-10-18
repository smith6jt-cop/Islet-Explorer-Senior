#!/usr/bin/env python3
"""Submit frozen section workers across the clive and maigan investment accounts."""
import argparse
import fcntl
import hashlib
import json
from pathlib import Path
import shlex
import shutil
import subprocess
import sys

from manifest import load_manifest
from workflow import DEFAULT_OUT, HERE, write_json, now

PYTHON = '/blue/maigan/smith6jt/miniforge3/envs/xenium_analysis/bin/python'


def freeze(root):
    sources = [HERE/n for n in ('workflow.py','evidence.py','manifest.py')]
    sources += [HERE.parent/'xenium_celltype/gdna.py']
    key = hashlib.sha256(b''.join(p.name.encode()+p.read_bytes() for p in sources)).hexdigest()[:16]
    dest = root/'code'/key
    for p in sources:
        q = dest/'core/phenocycler/scripts'/p.parent.name/p.name
        if q.exists():
            if q.read_bytes()!=p.read_bytes():raise ValueError('frozen code changed')
        else:
            q.parent.mkdir(parents=True,exist_ok=True); shutil.copy2(p,q)
    return key,dest/'core/phenocycler/scripts/sameslide/workflow.py'


def submit(manifest_path, root, donors, design, accounts=('maigan','clive'), memory_gb=24):
    root=Path(root).resolve();m=load_manifest(manifest_path)
    key,worker=freeze(root)
    selected=[s for s in m['sections'] if s['donor_id'] in donors and (design=='all' or s['design']==design)]
    jobs=root/'jobs'/key;jobs.mkdir(parents=True,exist_ok=True)
    config=jobs/'manifest.json'
    if config.exists() and json.loads(config.read_text())!=m:raise ValueError('frozen manifest changed')
    if not config.exists():write_json(config,m)
    # Each account runs a bounded wave, with afterany allowing a failed section to
    # leave independent sections runnable. These are memory limits, not cohort caps.
    slots={'maigan':[None]*4,'clive':[None]*2}
    per_account={a:0 for a in accounts}
    for i,s in enumerate(selected):
        account=accounts[i%len(accounts)]
        jobkey=hashlib.sha256((key+s['section_id']).encode()).hexdigest()[:12]
        receipt=jobs/(jobkey+'.json'); script=jobs/(jobkey+'.sbatch')
        with (jobs/(jobkey+'.lock')).open('a') as lock:
            fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
            if receipt.exists():
                r=json.loads(receipt.read_text())
                if not r.get('job_id'):raise ValueError(f'Interrupted sbatch requires reconciliation: {receipt}')
                account=r['account']
            else:
                cmd=[PYTHON,str(worker),'--manifest',str(config),'--section',s['section_id'],'--out',str(root)]
                script.write_text('#!/bin/bash\nset -euo pipefail\nexport OMP_NUM_THREADS=4 OPENBLAS_NUM_THREADS=4 MKL_NUM_THREADS=4 NUMBA_NUM_THREADS=4\n'+shlex.join(cmd)+'\n')
                slot=per_account[account]%len(slots[account]);dep=slots[account][slot]
                args=['sbatch','--parsable',f'--account={account}',f'--qos={account}',
                      '--partition=hpg-turin,hpg-milan','--cpus-per-task=4',f'--mem={memory_gb}G','--time=02:00:00',
                      f'--job-name=ss4-{jobkey}',f'--output={jobs}/{jobkey}_%j.log',
                      f'--chdir={root}']
                if dep:args += [f'--dependency=afterany:{dep}']
                r={'section_id':s['section_id'],'account':account,'code':key,'command':cmd,
                   'resources':{'cpus':4,'memory_gb':memory_gb},'submitted_utc':now(),'status':'submitting'}
                write_json(receipt,r)
                try:
                    jid=subprocess.check_output(args+[str(script)],text=True).strip().split(';')[0]
                except subprocess.CalledProcessError:
                    # A confirmed rejection can be retried; ambiguous interruptions remain recorded.
                    r['status']='submission_rejected';write_json(receipt,r);raise
                r.update(job_id=jid,status='submitted');write_json(receipt,r)
            slot=per_account[account]%len(slots[account]);slots[account][slot]=r['job_id'];per_account[account]+=1
            print(json.dumps(r),flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--manifest',required=True,type=Path);p.add_argument('--out',type=Path,default=DEFAULT_OUT)
    p.add_argument('--donors',nargs='+',required=True);p.add_argument('--design',choices=['all','same_slide','xenium_only'],default='all')
    p.add_argument('--accounts',nargs='+',choices=['clive','maigan'],default=['maigan','clive'])
    p.add_argument('--memory-gb',type=int,default=24)
    a=p.parse_args();submit(a.manifest,a.out,a.donors,a.design,a.accounts,a.memory_gb)
