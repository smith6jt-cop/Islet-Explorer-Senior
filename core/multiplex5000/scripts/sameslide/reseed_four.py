#!/usr/bin/env python3
"""Re-register a partial-ROI section from a declared seed, against its Xenium footprint only.

The 6533 sections were damaged after the PhenoCycler run, so Xenium was acquired on undamaged
partial ROIs (maintainer, 2026-10-17). The PhenoCycler section is therefore larger than the
Xenium one, and a seed searched against the whole PhenoCycler section can land on the wrong
tissue. Here the seed is an image-derived affine named on the command line (its displacement
field, if any, is dropped), PhenoCycler nuclei are restricted to the seeded Xenium footprint
plus a margin -- geometry only, no residuals, protein or RNA -- and the unchanged ICP, nucleus
hold-out and DAPI window gate follow. refine_four.py may then apply its even/odd correction.
"""
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
from scipy import ndimage

HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(HERE.parents[1]));sys.path.insert(0,str(HERE.parent/'restain'))
import export_dapi
import registration_inputs as RI
import run_align
import nuclei
import refine
from phenocycler.integration.transform import Transform
from manifest import load_manifest,get_section
from workflow import DEFAULT_OUT,section_root,identity,write_json,now
from register_points import held_out

GRID_UM=50.0


def footprint_mask(warped,fixed,margin_um):
    """Fixed points within margin_um of the occupied GRID_UM cells of the warped moving nuclei."""
    lo=np.minimum(warped.min(0),fixed.min(0))-margin_um
    shape=tuple((np.ceil((np.maximum(warped.max(0),fixed.max(0))+margin_um-lo)/GRID_UM)+1).astype(int)[::-1])
    occ=np.zeros(shape,bool);ij=((warped-lo)//GRID_UM).astype(int);occ[ij[:,1],ij[:,0]]=True
    occ=ndimage.binary_closing(occ,iterations=2)
    occ=ndimage.binary_dilation(occ,iterations=int(np.ceil(margin_um/GRID_UM)))
    fj=((fixed-lo)//GRID_UM).astype(int)
    return occ[fj[:,1],fj[:,0]]


def run(manifest_path,sid,seed_path,root,margin_um):
    section=get_section(load_manifest(manifest_path),sid);donor=section['donor_id']
    prep_path=Path(root)/'protein'/f'donor_id={donor}'/'preparation/current.json'
    prep=json.loads(prep_path.read_text())
    if prep.get('status')!='validated_complete':raise ValueError('protein geometry preparation is not validated')
    prepared=Path(prep['attempt']);seed_path=Path(seed_path).resolve()
    inputs={'section':section,'preparation':identity(prep_path),'code':identity(__file__),
            'seed':identity(seed_path),'footprint_margin_um':margin_um}
    fingerprint=hashlib.sha256(json.dumps(inputs,sort_keys=True).encode()).hexdigest()
    base=section_root(root,section)/'registration'
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
        xy,_ids=nuclei.xenium_nucleus_centroids(rec.bundle)
        np.save(dest/'xenium_nuclei_um.npy',xy)
        fixed=pd.read_parquet(prepared/section['roi']/'cells.parquet',columns=['nuc_x_um','nuc_y_um','qc_analysis_eligible'])
        fixed=fixed[fixed.qc_analysis_eligible][['nuc_x_um','nuc_y_um']].dropna().to_numpy()
        seed_json=json.loads(seed_path.read_text())
        seed=Transform(affine=np.array(seed_json['affine']),meta={'seed_file':str(seed_path),'displacement_dropped':bool(seed_json.get('has_displacement'))})
        if seed.is_mirrored:raise ValueError('declared seed is mirrored; a same-slide section cannot be')
        keep=footprint_mask(seed.apply(xy),fixed,margin_um);cropped=fixed[keep]
        write_json(dest/'footprint.json',{'method':'seeded Xenium nucleus occupancy, closed and dilated',
            'grid_um':GRID_UM,'margin_um':margin_um,'n_fixed_all':len(fixed),'n_fixed_in_footprint':int(keep.sum()),
            'n_moving':len(xy),'seed':identity(seed_path),'seed_rotation_deg':seed.rotation_deg,
            'seed_translation_um':seed.translation_um,'displacement_dropped':seed.meta['displacement_dropped']})
        tf=refine.refine_icp(seed,xy,cropped)
        tf.save(dest,stem='four')
        score=refine.score_transform(tf,xy,cropped)
        split=held_out(xy,cropped,tf)
        write_json(dest/'nucleus_residuals.json',{'score':score,'spatial_split':split,'seed_method':'declared_image_affine_footprint_ICP'})
        gate=run_align.align_one(rec,mode='nucleus',windows='occupied',transform='four',out_dir=out)
        gate['registration_input']['pixel_size_source']=facts['pixel_size_source']
        write_json(dest/'alignment_four.json',gate)
        verdict=gate['alignment']['verdict']
        report.update(status='validated_complete' if verdict=='PASS' else 'scientifically_unavailable',
                      verdict=verdict,reason=gate['alignment']['verdict_why'],seed_method='declared_image_affine_footprint_ICP',
                      transform_directory=str(dest),transform_stem='four',finished_utc=now(),
                      validation={'DAPI_image_gate':verdict,'nucleus_residuals_recorded':True,'spatial_holdout_recorded':True})
        report['outputs']=[identity(p) for p in out.rglob('*') if p.is_file() and p.name!='receipt.json']
        write_json(out/'receipt.json',report);write_json(base/'current.json',report)
        print(f'{sid}: reseeded registration {verdict}: {report["reason"]}',flush=True)
    except Exception as exc:
        report.update(status='failed',error=str(exc),traceback=traceback.format_exc(),finished_utc=now())
        write_json(out/'receipt.json',report);write_json(base/'last_failure.json',report);raise


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--manifest',required=True,type=Path)
    p.add_argument('--section',required=True);p.add_argument('--seed',required=True,type=Path)
    p.add_argument('--margin-um',type=float,default=100.0);p.add_argument('--out',type=Path,default=DEFAULT_OUT)
    a=p.parse_args();run(a.manifest,a.section,a.seed,a.out,a.margin_um)
