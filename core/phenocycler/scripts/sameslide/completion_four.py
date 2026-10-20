#!/usr/bin/env python3
"""The completion matrix: what every stage holds for every section, read from receipts only.

A cell is the stage's current receipt as it stands -- validated_complete, scientifically_unavailable
with its reason, failed, or not run -- and STALE when a DATA file it recorded as an input (an upstream
receipt, a delivered table, a bundle file) has changed size or mtime since. A changed code file is
reported beside it, not as staleness: an edit need not change what the stage wrote, and rerunning is
the way to find out. A receipt's record of its own current.json (refine_four's source) is not an
input. Nothing is recomputed and no result is judged here.
"""
import argparse
import hashlib
import json
import os
from pathlib import Path
from uuid import uuid4

import pandas as pd

from manifest import load_manifest
from workflow import DEFAULT_OUT,section_root,identity,write_json,now

SECTION_STAGES=('rna','rna_description','registration','registration_recheck','joined','association','calls','context','islets',
                'node_summaries')
SAME_SLIDE=('rna','rna_description','preparation','alignment','registration','registration_recheck','joined',
            'association','calls','context','comparison','islets','node_summaries')
XENIUM_ONLY=('rna','rna_description','comparison')
#: A transform registers images. A re-prepared segmentation is answered by registration_recheck --
#: the transform scored, unchanged, against the new nuclei, which the join requires -- not by
#: registering again (maintainer decision 2026-10-19), so it does not make the registration stale.
RECHECKED={'registration':'preparation/current.json'}
#: A stage that belongs to one tissue only; elsewhere the cell is n/a, not missing...
TISSUE_ONLY={'islets':'pancreas','node_summaries':'pancreatic_lymph_node'}
#: ...unless a section of the other tissue imaged it: 6347's and 6580's pancreas regions also cover
#: their same-block node, which their calls receipts record as Xenium cells won by pLN cells.
IMAGED_ELSEWHERE={'node_summaries':'pLN|'}


def recorded(value):
    """Every identity record (path, bytes, mtime_ns, sha256) anywhere inside a receipt's inputs."""
    if isinstance(value,dict):
        if {'path','bytes','mtime_ns','sha256'}<=set(value):yield value
        else:
            for v in value.values():yield from recorded(v)
    elif isinstance(value,list):
        for v in value:yield from recorded(v)


def changed(receipt,own=None):
    """(data, code) inputs whose file is gone or differs in size or mtime from the record."""
    data,code=[],[]
    for r in recorded(receipt.get('inputs',{})):
        p=Path(r['path'])
        if own is not None and p.resolve()==Path(own).resolve():continue
        bucket=code if p.suffix=='.py' else data
        if not p.exists():bucket.append(f'missing {p.name}')
        else:
            s=p.stat()
            if (s.st_size,s.st_mtime_ns)!=(r['bytes'],r['mtime_ns']):bucket.append(f'changed {p.parent.name}/{p.name}')
    return data,code


def cell(path):
    stage=Path(path).parent
    if not Path(path).exists():
        failure=stage/'last_failure.json'
        if failure.exists():return {'status':'failed','reason':json.loads(failure.read_text()).get('error','')[:200]}
        return {'status':'not_run','reason':''}
    r=json.loads(Path(path).read_text())
    stale,code=changed(r,own=path)
    if stage.name in RECHECKED:
        covered=[x for x in stale if x.endswith(RECHECKED[stage.name])]
        stale=[x for x in stale if x not in covered]
    status=r.get('status','')
    if stale and status in ('validated_complete','scientifically_unavailable'):status='stale'
    return {'status':status,'reason':('; '.join(stale) if status=='stale' else r.get('reason') or '')[:300],
            'attempt':r.get('attempt',''),'stale_inputs':'; '.join(stale),'code_changed_since':'; '.join(code)}


def receipt_path(root,section,stage):
    d=section['donor_id']
    if stage in SECTION_STAGES:return section_root(root,section)/stage/'current.json'
    if stage in ('preparation','alignment'):return Path(root)/'protein'/f'donor_id={d}'/stage/'current.json'
    if stage=='comparison':return Path(root)/'comparisons'/f'donor_id={d}'/f"roi={section['roi']}"/'current.json'
    raise ValueError(stage)


def applies(root,section,stage):
    """Whether a stage belongs to this section, from the manifest and, across tissues, its calls receipt."""
    if TISSUE_ONLY.get(stage,section['tissue'])==section['tissue']:return True
    p=receipt_path(root,section,'calls')
    if stage not in IMAGED_ELSEWHERE or not p.exists():return False
    won=(json.loads(p.read_text()).get('summary') or {}).get('status_by_winning_tissue',{})
    return any(k.startswith(IMAGED_ELSEWHERE[stage]) for k in won)


def run(manifest_path,root):
    m=load_manifest(manifest_path);rows=[]
    for s in m['sections']:
        for stage in (SAME_SLIDE if s['design']=='same_slide' else XENIUM_ONLY):
            c=cell(receipt_path(root,s,stage)) if applies(root,s,stage) else {'status':'not_applicable','reason':''}
            rows.append({'section_id':s['section_id'],'design':s['design'],'tissue':s['tissue'],'stage':stage,**c})
    table=pd.DataFrame(rows)
    fp=hashlib.sha256(table.to_csv(index=False).encode()).hexdigest()
    stage=Path(root)/'completion';out=stage/'attempts'/(fp[:16]+'-'+uuid4().hex[:8]);out.mkdir(parents=True)
    table.to_csv(out/'completion.csv',index=False)
    short={'validated_complete':'done','scientifically_unavailable':'unavailable','not_run':'-',
           'failed':'FAILED','stale':'STALE','not_applicable':'n/a'}
    lines=[]
    for design,stages in (('same_slide',SAME_SLIDE),('xenium_only',XENIUM_ONLY)):
        t=table[table.design==design]
        piv=t.pivot(index='section_id',columns='stage',values='status').reindex(columns=list(stages)).fillna('').replace(short)
        lines+=[f'## {design.replace("_","-")} sections','','| section | '+' | '.join(stages)+' |','|'+'---|'*(len(stages)+1)]
        lines+=['| '+sid+' | '+' | '.join(piv.loc[sid])+' |' for sid in piv.index]+['']
    notes=table[table.status.isin(['scientifically_unavailable','failed','stale'])]
    lines+=['## Reasons','']+[f'- **{r.section_id} {r.stage}** ({short[r.status]}): {r.reason}' for r in notes.itertuples()]
    edited=table[table.get('code_changed_since',pd.Series('',index=table.index)).fillna('').astype(bool)]
    if len(edited):
        lines+=['','## Code edited since the receipt (rerun to confirm the output)','']
        lines+=[f'- {r.section_id} {r.stage}: {r.code_changed_since}' for r in edited.itertuples()]
    (out/'completion.md').write_text('\n'.join(lines)+'\n')
    rec={'status':'validated_complete','fingerprint':fp,'finished_utc':now(),'attempt':str(out),'job_id':os.environ.get('SLURM_JOB_ID'),
         'inputs':{'manifest':identity(manifest_path),'code':identity(__file__)},
         'summary':table.status.value_counts().to_dict(),'validation':{'receipts_only':True},
         'outputs':[identity(p) for p in out.iterdir() if p.is_file()]}
    write_json(out/'receipt.json',rec);write_json(stage/'current.json',rec)
    print((out/'completion.md').read_text())


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--manifest',required=True,type=Path)
    p.add_argument('--out',type=Path,default=DEFAULT_OUT);a=p.parse_args();run(a.manifest,a.out)
