#!/usr/bin/env python3
"""Protein-guided RNA context with held-out evaluation, one section at a time.

``context.evaluate`` on the section's RNA base and its transferred protein calls: two spatial
folds, each cell's class prior fitted on the OTHER fold's RNA only, raw counts never replaced.
Its output is protein-GUIDED support and is reported apart from ``associate_four.py``'s
independent protein-RNA association, whose RNA side no protein value enters.
"""
import argparse
import hashlib
import json
import os
from pathlib import Path
import traceback
from uuid import uuid4

import context
from join_four import current
from manifest import load_manifest,get_section
from workflow import DEFAULT_OUT,section_root,identity,write_json,now,valid_receipt


def run(manifest_path,sid,root):
    section=get_section(load_manifest(manifest_path),sid);base=section_root(root,section)
    paths={'calls':base/'calls/current.json','rna':base/'rna/current.json'}
    calls_rec=current(paths['calls'],'calls')
    for p in ((paths['rna'],) if calls_rec['status']=='scientifically_unavailable' else paths.values()):
        if not valid_receipt(p):raise ValueError(f'{p}: not validated')
    params={'min_cells':context.MIN_CLASS_CELLS,'min_purity':context.MIN_CALL_PURITY,'min_overlap':context.MIN_OVERLAP}
    inputs={**{k:identity(p) for k,p in paths.items()},'code':[identity(__file__),identity(context.__file__)],'parameters':params}
    fp=hashlib.sha256(json.dumps(inputs,sort_keys=True).encode()).hexdigest()
    stage=base/'context';stage.mkdir(parents=True,exist_ok=True)
    if valid_receipt(stage/'current.json',fp):print(f'{sid}: context reused',flush=True);return
    out=stage/'attempts'/(fp[:16]+'-'+uuid4().hex[:8]);out.mkdir(parents=True)
    rec={'status':'running','section_id':sid,'fingerprint':fp,'inputs':inputs,'started_utc':now(),
         'attempt':str(out),'job_id':os.environ.get('SLURM_JOB_ID')}
    write_json(out/'receipt.json',rec)
    if calls_rec['status']=='scientifically_unavailable':
        # No class to guide the RNA: the section's current context says so rather than standing empty.
        write_json(out/'unavailable.json',{'calls':calls_rec.get('reason')})
        rec.update(status='scientifically_unavailable',reason=f"calls: {calls_rec.get('reason','')}",finished_utc=now(),
                   validation={'calls_available':False},outputs=[identity(out/'unavailable.json')])
        write_json(out/'receipt.json',rec);write_json(stage/'current.json',rec)
        print(f"{sid}: context unavailable: {rec['reason']}",flush=True);return
    try:
        rna=Path(json.loads(paths['rna'].read_text())['attempt'])/'rna.h5ad'
        calls=Path(json.loads(paths['calls'].read_text())['attempt'])/'protein_calls.parquet'
        meta=context.evaluate(rna,calls,out,**params)
        rec.update(status='validated_complete',finished_utc=now(),
                   summary={k:meta[k] for k in ('n_cells','n_genes','n_supported_entries','fallback_counts','split_axis','split_um')},
                   validation={'held_out_spatial_folds':2,'raw_counts_unchanged':True,'posterior_calibrated':meta['posterior_calibrated'],
                               'reported_apart_from_independent_association':True},
                   outputs=[identity(p) for p in out.iterdir() if p.is_file() and p.name!='receipt.json'])
        write_json(out/'receipt.json',rec);write_json(stage/'current.json',rec)
        print(f"{sid}: context support {meta['n_supported_entries']:,} entries; fallbacks {meta['fallback_counts']}",flush=True)
    except Exception as exc:
        rec.update(status='failed',error=str(exc),traceback=traceback.format_exc(),finished_utc=now())
        write_json(out/'receipt.json',rec);write_json(stage/'last_failure.json',rec);raise


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__,formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('--manifest',required=True,type=Path);p.add_argument('--section',required=True)
    p.add_argument('--out',type=Path,default=DEFAULT_OUT)
    a=p.parse_args();run(a.manifest,a.section,a.out)
