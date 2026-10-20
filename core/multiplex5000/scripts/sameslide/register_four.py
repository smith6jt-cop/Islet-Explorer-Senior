#!/usr/bin/env python3
"""Register each physical section and apply the established DAPI window gate."""
import argparse
import hashlib
import json
import os
from pathlib import Path
import sys
import traceback
from uuid import uuid4

import numpy as np
import pandas as pd

HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(HERE.parents[1]));sys.path.insert(0,str(HERE.parent/'restain'))
import export_dapi
import registration_inputs as RI
import run_image_seed
import run_align
import nuclei
import refine
from manifest import load_manifest,get_section
from workflow import DEFAULT_OUT,section_root,identity,write_json,now,valid_receipt
from register_points import fit,held_out


def run(manifest_path,sid,root):
    section=get_section(load_manifest(manifest_path),sid);donor=section['donor_id']
    prep_path=Path(root)/'protein'/f'donor_id={donor}'/'preparation/current.json'
    if not valid_receipt(prep_path):raise ValueError('protein geometry preparation is not validated')
    prep=json.loads(prep_path.read_text());prepared=Path(prep['attempt'])
    inputs={'section':section,'preparation':identity(prep_path),'code':identity(__file__)}
    fingerprint=hashlib.sha256(json.dumps(inputs,sort_keys=True).encode()).hexdigest()
    base=section_root(root,section)/'registration';base.mkdir(parents=True,exist_ok=True)
    if valid_receipt(base/'current.json',fingerprint):
        print(f'{sid}: registration reused',flush=True);return
    out=base/'attempts'/(fingerprint[:16]+'-'+uuid4().hex[:8]);out.mkdir(parents=True)
    report={'status':'running','section_id':sid,'fingerprint':fingerprint,'inputs':inputs,
            'started_utc':now(),'attempt':str(out),'job_id':os.environ.get('SLURM_JOB_ID')}
    write_json(out/'receipt.json',report)
    try:
        ph=section['pheno'];facts=RI.qptiff_registration_facts(Path(ph['image']),allow_dtypes=('uint8','uint16'),pixel_um_override=ph['pixel_size_um'])
        rec=RI.RegistrationInput(donor_id=donor,serial=section['xenium_serial'],roi=section['roi'],qptiff=Path(ph['image']),
             bundle=Path(section['bundle']),geojson=None,provenance_source=str(prepared/'geometry_qc.json'),
             **{k:facts[k] for k in ('pixel_um','dapi_page','channel_names','pyramid_shapes','pixel_size_source')})
        export_dapi.export_one(rec,out_root=out/'valis')
        dest=out/'valis'/f'donor_id={donor}'
        xy,xids=nuclei.xenium_nucleus_centroids(rec.bundle)
        np.save(dest/'xenium_nuclei_um.npy',xy)
        fixed=pd.read_parquet(prepared/section['roi']/'cells.parquet',columns=['nuc_x_um','nuc_y_um','qc_analysis_eligible'])
        fixed=fixed[fixed.qc_analysis_eligible][['nuc_x_um','nuc_y_um']].dropna().to_numpy()
        try:
            tf,image_report=run_image_seed.fit_image_seed(dest)
            write_json(dest/'image_seed_report.json',image_report)
            seed_method='DAPI_SIFT_ECC'
        except run_image_seed.ImageSeedError as exc:
            # A distinct predeclared geometry seed may rescue failed image feature matching.
            # Both candidates still face the same independent final image gate.
            write_json(dest/'image_seed_failure.json',{'error':str(exc)})
            tf,point_report=fit(xy,fixed)
            write_json(dest/'point_seed.json',point_report)
            seed_method='nucleus_density_rotation_search_ICP'
        # Refinement uses nuclei, never protein identity or RNA expression.
        tf=refine.refine_icp(tf,xy,fixed)
        tf.save(dest,stem='four')
        score=refine.score_transform(tf,xy,fixed)
        split=held_out(xy,fixed,tf)
        write_json(dest/'nucleus_residuals.json',{'score':score,'spatial_split':split,'seed_method':seed_method})
        gate=run_align.align_one(rec,mode='nucleus',windows='occupied',transform='four',out_dir=out)
        gate['registration_input']['pixel_size_source']=facts['pixel_size_source']
        write_json(dest/'alignment_four.json',gate)
        verdict=gate['alignment']['verdict']
        report.update(status='validated_complete' if verdict=='PASS' else 'scientifically_unavailable',
                      verdict=verdict,reason=gate['alignment']['verdict_why'],
                      transform_directory=str(dest),transform_stem='four',finished_utc=now(),
                      validation={'DAPI_image_gate':verdict,'nucleus_residuals_recorded':True,'spatial_holdout_recorded':True})
        report['outputs']=[identity(p) for p in out.rglob('*') if p.is_file() and p.name!='receipt.json']
        write_json(out/'receipt.json',report);write_json(base/'current.json',report)
        print(f'{sid}: registration {verdict}: {report["reason"]}',flush=True)
    except Exception as exc:
        report.update(status='failed',error=str(exc),traceback=traceback.format_exc(),finished_utc=now())
        write_json(out/'receipt.json',report);write_json(base/'last_failure.json',report);raise


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--manifest',required=True,type=Path)
    p.add_argument('--section',required=True);p.add_argument('--out',type=Path,default=DEFAULT_OUT)
    a=p.parse_args();run(a.manifest,a.section,a.out)
