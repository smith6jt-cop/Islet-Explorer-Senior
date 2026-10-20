#!/usr/bin/env python3
"""The validated transform, unchanged, against the CURRENT preparation's nuclei.

A transform registers images, not segmentations: the DAPI window gate reads only the Xenium
nuclei and the PhenoCycler DAPI image, so another segmentation of the same image leaves the
transform and its gate verdict as they were (maintainer decision 2026-10-19, for the workstation
segmentation of 6347 and 6533). What a segmentation can still get wrong is its own frame -- pixel
size, origin, a flip. So the registration's transform is scored, with no refit, against the
current preparation's eligible nuclei, and the median mutual-pair distance is held to the gate's
own pass threshold. The join requires this receipt to name the preparation it reads.
"""
import argparse
import hashlib
import json
import os
from pathlib import Path
import sys
import traceback
from uuid import uuid4

import pandas as pd

HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(HERE.parents[1]));sys.path.insert(0,str(HERE.parent/'restain'))
import nuclei
import refine
from phenocycler.integration.transform import Transform
from manifest import load_manifest,get_section
from workflow import DEFAULT_OUT,section_root,identity,write_json,now,valid_receipt


def gate_threshold(directory,stem):
    """threshold_pass_um of the gate that judged this transform, from that gate's own record."""
    hits=sorted(Path(directory).glob(f'alignment_{stem}*.json'))
    if len(hits)!=1:raise ValueError(f'{directory}: expected one gate record for {stem}, found {[h.name for h in hits]}')
    a=json.loads(hits[0].read_text())['alignment']
    return float(a['threshold_pass_um']),hits[0]


def run(manifest_path,sid,root):
    section=get_section(load_manifest(manifest_path),sid);donor=section['donor_id'];base=section_root(root,section)
    paths={'registration':base/'registration/current.json',
           'preparation':Path(root)/'protein'/f'donor_id={donor}'/'preparation/current.json'}
    reg=json.loads(paths['registration'].read_text())
    if not valid_receipt(paths['preparation']):raise ValueError('protein preparation is not validated')
    inputs={**{k:identity(p) for k,p in paths.items()},'code':identity(__file__)}
    fp=hashlib.sha256(json.dumps(inputs,sort_keys=True).encode()).hexdigest()
    stage=base/'registration_recheck';stage.mkdir(parents=True,exist_ok=True)
    if valid_receipt(stage/'current.json',fp):print(f'{sid}: recheck reused',flush=True);return
    out=stage/'attempts'/(fp[:16]+'-'+uuid4().hex[:8]);out.mkdir(parents=True)
    rec={'status':'running','section_id':sid,'fingerprint':fp,'inputs':inputs,'started_utc':now(),
         'attempt':str(out),'job_id':os.environ.get('SLURM_JOB_ID')}
    write_json(out/'receipt.json',rec)
    try:
        if reg.get('status')!='validated_complete':
            write_json(out/'unavailable.json',{'registration':reg.get('status'),'reason':reg.get('reason')})
            rec.update(status='scientifically_unavailable',reason=f"registration {reg.get('status')}: {reg.get('reason')}",
                       finished_utc=now(),validation={'registration_validated':False},outputs=[identity(out/'unavailable.json')])
            write_json(out/'receipt.json',rec);write_json(stage/'current.json',rec);return
        tf=Transform.load(Path(reg['transform_directory']),stem=reg['transform_stem'])
        threshold,gate=gate_threshold(reg['transform_directory'],reg['transform_stem'])
        xy,_=nuclei.xenium_nucleus_centroids(Path(section['bundle']))
        prep=Path(json.loads(paths['preparation'].read_text())['attempt'])
        fixed=pd.read_parquet(prep/section['roi']/'cells.parquet',columns=['nuc_x_um','nuc_y_um','qc_analysis_eligible'])
        fixed=fixed[fixed.qc_analysis_eligible][['nuc_x_um','nuc_y_um']].dropna().to_numpy()
        score=refine.score_transform(tf,xy,fixed)
        recorded=Path(reg['transform_directory'])/'nucleus_residuals.json'
        p50=score['pair_distance_um'].get('p50')
        ok=p50 is not None and p50<=threshold
        write_json(out/'recheck.json',{'score':score,'gate_threshold_pass_um':threshold,'gate_record':str(gate),
                   'registration_recorded_score':json.loads(recorded.read_text())['score'] if recorded.exists() else None,
                   'rule':'median mutual nucleus-pair distance of the unchanged transform <= the gate pass threshold'})
        rec.update(status='validated_complete' if ok else 'scientifically_unavailable',finished_utc=now(),
                   reason=f'median nucleus pair {p50:.2f} um {"<=" if ok else ">"} gate pass {threshold:.2f} um' if p50 is not None else 'no nucleus pairs',
                   summary={'pair_p50_um':p50,'paired_fraction':score['paired_fraction'],'threshold_pass_um':threshold,
                            'matched_within_2um':score['matched_within'].get('2um'),'matched_within_2um_null':score['matched_within_null'].get('2um')},
                   validation={'transform_unchanged':True,'no_refit':True,'within_gate_threshold':bool(ok)},
                   outputs=[identity(out/'recheck.json')])
        write_json(out/'receipt.json',rec);write_json(stage/'current.json',rec)
        print(f"{sid}: recheck {rec['status']}: {rec['reason']}",flush=True)
    except Exception as exc:
        rec.update(status='failed',error=str(exc),traceback=traceback.format_exc(),finished_utc=now())
        write_json(out/'receipt.json',rec);write_json(stage/'last_failure.json',rec);raise


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__,formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('--manifest',required=True,type=Path);p.add_argument('--section',required=True)
    p.add_argument('--out',type=Path,default=DEFAULT_OUT)
    a=p.parse_args();run(a.manifest,a.section,a.out)
