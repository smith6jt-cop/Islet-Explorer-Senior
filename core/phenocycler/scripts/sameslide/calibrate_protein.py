#!/usr/bin/env python3
"""Tissue-specific UniFORM fits for the four acquisitions. Superseded by the full-dataset UniFORM; see HANDOFF_2026-10-18.md."""
import argparse
import contextlib
import hashlib
import io
import json
import os
from pathlib import Path
import sys
import traceback
from uuid import uuid4

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(HERE.parent/'celltype'))
import uniform_adapter as U
import score as S
import panel
from manifest import load_manifest
from workflow import DEFAULT_OUT,identity,write_json,now,valid_receipt


def run(manifest_path,roi,root):
    sections=[r for r in load_manifest(manifest_path)['sections'] if r['design']=='same_slide' and r['roi']==roi]
    if len(sections)!=4:raise ValueError('calibration requires all four tissue-matched acquisitions')
    sources={}
    for r in sections:
        p=Path(root)/'protein'/f"donor_id={r['donor_id']}"/'redsea/current.json'
        if not valid_receipt(p):raise ValueError(f'{p}: REDSEA not validated')
        sources[r['donor_id']]=p
    inputs={'redsea':{d:identity(p) for d,p in sources.items()},'code':identity(__file__),
'uniform_commit':U.UNIFORM_COMMIT,'roi':roi,
            'scale_code':identity(S.__file__),'adapter_code':identity(U.__file__)}
    fingerprint=hashlib.sha256(json.dumps(inputs,sort_keys=True).encode()).hexdigest()
    base=Path(root)/'protein'/f'calibration_{roi}';base.mkdir(parents=True,exist_ok=True)
    if valid_receipt(base/'current.json',fingerprint):return
    out=base/'attempts'/(fingerprint[:16]+'-'+uuid4().hex[:8]);out.mkdir(parents=True)
    rec={'status':'running','fingerprint':fingerprint,'inputs':inputs,'started_utc':now(),'attempt':str(out),'job_id':os.environ.get('SLURM_JOB_ID')}
    write_json(out/'receipt.json',rec)
    try:
        reg,norm,_=U._import_uniform()
        frames={d:pd.read_parquet(Path(json.loads(p.read_text())['attempt'])/roi/'expression.parquet') for d,p in sources.items()}
        markers=sorted(set.union(*(set(panel.markers(include_variant=True))&set(f.columns) for f in frames.values())))
        shifts={};usability=[]
        for marker in markers:
            vals={d:f.loc[f.qc_analysis_eligible,marker].to_numpy(float) for d,f in frames.items() if marker in f}
            available={d:v for d,v in vals.items() if (np.isfinite(v)&(v>=1)).sum()>=1000 and np.nanstd(v)>0}
            if len(available)<2:
                shifts[marker]={'status':'unavailable','reason':'fewer_than_two_measurable_acquisitions'};continue
            counts,ranges=U.marker_histograms(available);donors=list(available)
            hist={marker:{'counts':[counts[d] for d in donors],'bin_edges':[None]*len(donors)}}
            with contextlib.redirect_stdout(io.StringIO()):
                shift_map,refs,_=reg.automatic_registration(hist,[marker],[marker],donors,None,None,U.NUM_BINS)
            plt.close('all')
            bins=dict(zip(donors,shift_map[marker]));before=U.background_mode_spread(counts,ranges)
            after=U.background_mode_spread(counts,ranges,bins)
            logs=dict(zip(donors,norm.calculate_shift_in_log_pixels({marker:ranges},[marker],U.NUM_BINS,shift_map)[marker]))
            scales={d:float(np.exp(-logs[d])) for d in donors}
            shifts[marker]={'status':'aligned','reference_donor':refs[0],'scale':scales,'mode_sd_before':float(before['sd_log']),
                            'mode_sd_after':float(after['sd_log']),'cross_donor_cutoff_supported':bool(after['sd_log']<=.10)}
            for d in frames:
                if marker not in frames[d]:continue
                if d in scales:frames[d][marker]*=scales[d]
                else:frames[d][marker]=np.nan
            for d,v in vals.items():
                finite=np.isfinite(v);usability.append({'donor_id':d,'marker':marker,'n_finite':int(finite.sum()),
                    'missing_fraction':float(1-finite.mean()),'scale':scales.get(d,np.nan),
                    'mode_sd_after':float(after['sd_log']),'cross_donor_cutoff_supported':after['sd_log']<=.10})
            print(f'{roi}/{marker}: UniFORM mode SD {before["sd_log"]:.3f} -> {after["sd_log"]:.3f}',flush=True)
        write_json(out/'shifts.json',shifts)
        pd.DataFrame(usability).to_csv(out/'marker_usability.csv',index=False)
        rng=np.random.default_rng(0);pool={m:[] for m in markers}
        for d,f in frames.items():
            sub=f[f.qc_analysis_eligible]
            if len(sub)>200000:sub=sub.iloc[rng.choice(len(sub),200000,replace=False)]
            for marker in markers:
                if marker in sub:pool[marker].append(sub[marker].to_numpy(float))
        pool={m:np.concatenate(v) for m,v in pool.items() if v};floors={}
        for marker,v in pool.items():
            v=v[np.isfinite(v)&(v>0)]
            if len(v)<1000:continue
            h,b=np.histogram(np.log1p(v),bins=512)
            floors[marker]=float(np.expm1(b[1:][int(np.argmax(h[1:]))+1]))
        scale=S.marker_scale(pool,floors);scale.to_csv(out/'scale.csv')
        # No typing here: it belongs to a per-donor stage on the full-dataset UniFORM (HANDOFF_2026-10-18.md).
        for d,f in frames.items():
            dest=out/f'donor_id={d}';dest.mkdir()
            f.to_parquet(dest/'aligned.parquet',index=False)
        rec.update(status='validated_complete',finished_utc=now(),validation={'four_acquisitions':True,'separate_tissue_fit':True,
                   'missing_measurements_not_zero':True},
                   outputs=[identity(p) for p in out.rglob('*') if p.is_file() and p.name!='receipt.json'])
        write_json(out/'receipt.json',rec);write_json(base/'current.json',rec)
    except Exception as exc:
        rec.update(status='failed',error=str(exc),traceback=traceback.format_exc(),finished_utc=now())
        write_json(out/'receipt.json',rec);write_json(base/'last_failure.json',rec);raise


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--manifest',required=True,type=Path)
    p.add_argument('--roi',required=True,choices=['panc','pln_1']);p.add_argument('--out',type=Path,default=DEFAULT_OUT)
    a=p.parse_args();run(a.manifest,a.roi,a.out)
