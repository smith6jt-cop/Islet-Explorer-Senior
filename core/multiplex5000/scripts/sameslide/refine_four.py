#!/usr/bin/env python3
"""One DAPI residual correction fitted on even windows, gated on withheld odd windows."""
import argparse
import hashlib
import json
import os
from pathlib import Path
import sys
from uuid import uuid4

import register_four as registration
import run_intensity_refine as IR
from phenocycler.integration.transform import Transform
from manifest import load_manifest,get_section
from workflow import DEFAULT_OUT,section_root,identity,write_json,now,sha


def run(manifest_path,sid,root):
    section=get_section(load_manifest(manifest_path),sid)
    base=section_root(root,section)/'registration';source_path=base/'current.json'
    source=json.loads(source_path.read_text())
    if source['status']=='validated_complete':return
    if source.get('refinement_attempted'):return
    if source['status']!='scientifically_unavailable':raise ValueError('source registration is incomplete')
    for v in source['outputs']:
        if sha(v['path'])!=v['sha256']:raise ValueError('source registration was modified')
    donor=section['donor_id'];old=Path(source['transform_directory'])
    gate=json.loads((old/'alignment_four.json').read_text())
    inputs={'source':identity(source_path),'code':identity(__file__),'refiner':identity(IR.__file__)}
    fp=hashlib.sha256(json.dumps(inputs,sort_keys=True).encode()).hexdigest()
    out=base/'attempts'/(fp[:16]+'-'+uuid4().hex[:8]);dest=out/'valis'/f'donor_id={donor}';dest.mkdir(parents=True)
    rec={'status':'running','fingerprint':fp,'inputs':inputs,'attempt':str(out),'started_utc':now(),
         'section_id':sid,'job_id':os.environ.get('SLURM_JOB_ID'),'refinement_attempted':True}
    write_json(out/'receipt.json',rec)
    try:
        train=IR._training_windows(gate,'even')
        fitted=IR._fit_correction(train);fitted.update(training_fold='even',withheld_fold='odd')
        tf=IR._corrected_transform(Transform.load(old,stem=source['transform_stem']),fitted,gate['alignment']['windows'])
        tf.save(dest,stem='four_refined')
        write_json(dest/'refinement.json',{'method':'existing bounded robust affine residual correction',
             'training_fold':'even','withheld_gate_fold':'odd','n_training':len(train),
             'training_median_before_um':fitted['baseline_median_um'],'training_median_after_um':fitted['residual_median_um']})
        ph=section['pheno'];facts=registration.RI.qptiff_registration_facts(Path(ph['image']),allow_dtypes=('uint8','uint16'),pixel_um_override=ph['pixel_size_um'])
        r=registration.RI.RegistrationInput(donor_id=donor,serial=section['xenium_serial'],roi=section['roi'],qptiff=Path(ph['image']),
             bundle=Path(section['bundle']),geojson=None,provenance_source=str(source_path),
             **{k:facts[k] for k in ('pixel_um','dapi_page','channel_names','pyramid_shapes','pixel_size_source')})
        result=registration.run_align.align_one(r,mode='nucleus',windows='occupied',window_fold='odd',transform='four_refined',out_dir=out)
        result['registration_input']['pixel_size_source']=facts['pixel_size_source']
        result['alignment']['certifies_mapping']['source']=facts['pixel_size_source']
        write_json(dest/'alignment_four_refined_odd.json',result)
        verdict=result['alignment']['verdict']
        rec.update(status='validated_complete' if verdict=='PASS' else 'scientifically_unavailable',verdict=verdict,
           reason=result['alignment']['verdict_why'],transform_directory=str(dest),transform_stem='four_refined',
           validation={'DAPI_image_gate':verdict,'withheld_gate_fold':'odd','training_fold':'even'})
    except IR.IntensityRefineError as exc:
        write_json(out/'refinement_unavailable.json',{'reason':str(exc),'source_registration':str(old)})
        rec.update(status='scientifically_unavailable',reason=str(exc),verdict=source['verdict'],
           transform_directory=str(old),transform_stem=source['transform_stem'],validation=source['validation'])
    rec.update(finished_utc=now(),outputs=[identity(p) for p in out.rglob('*') if p.is_file() and p.name!='receipt.json'])
    write_json(out/'receipt.json',rec);write_json(base/'current.json',rec)
    print(f"{sid}: {rec['status']}: {rec['reason']}",flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--manifest',required=True,type=Path)
    p.add_argument('--section',required=True);p.add_argument('--out',type=Path,default=DEFAULT_OUT)
    a=p.parse_args();run(a.manifest,a.section,a.out)
