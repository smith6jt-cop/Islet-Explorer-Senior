#!/usr/bin/env python3
"""Full-resolution REDSEA selected expression, isolated from the main cohort."""
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
from scipy import sparse
from skimage.segmentation import find_boundaries

HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(HERE.parents[1]));sys.path.insert(0,str(HERE.parent/'celltype'))
from phenocycler import redsea as R
from phenocycler.marker_registry import load_registry
import panel
import usability
from manifest import load_manifest
from workflow import DEFAULT_OUT,identity,write_json,now,valid_receipt
from read_pheno_tsv import sanitize


def run(manifest_path,donor,root):
    section=next(r for r in load_manifest(manifest_path)['sections'] if r['donor_id']==donor and r['design']=='same_slide')
    ph=section['pheno'];base=Path(root)/'protein'/f'donor_id={donor}'
    prep_path=base/'preparation/current.json'
    if not valid_receipt(prep_path):raise ValueError('geometry preparation is not validated')
    prep=json.loads(prep_path.read_text());prepared=Path(prep['attempt'])
    inputs={'geometry':identity(ph['geometry']),'image':identity(ph['image']),'preparation':identity(prep_path),
            'script':identity(__file__),'redsea_code':identity(R.__file__),
            'registry':identity(HERE.parents[1]/'phenocycler/marker_registry.json')}
    fingerprint=hashlib.sha256(json.dumps(inputs,sort_keys=True).encode()).hexdigest()
    base=base/'redsea';base.mkdir(exist_ok=True)
    if valid_receipt(base/'current.json',fingerprint):print(f'{donor}: REDSEA reused',flush=True);return
    out=base/'attempts'/(fingerprint[:16]+'-'+uuid4().hex[:8]);out.mkdir(parents=True)
    record={'status':'running','fingerprint':fingerprint,'inputs':inputs,'attempt':str(out),
            'started_utc':now(),'job_id':os.environ.get('SLURM_JOB_ID')};write_json(out/'receipt.json',record)
    try:
        frames={roi:pd.read_parquet(prepared/roi/'cells.parquet') for roi in ('panc','pln_1')}
        cells=pd.concat(frames.values(),ignore_index=True).set_index('object_id',drop=False)
        tf,pages,names,shape,_=R.qptiff_channels(Path(ph['image']),0)
        markers=[R.SANITIZE(n) for n in names]
        if len(set(markers))!=len(markers):raise ValueError('duplicate marker names')
        policies={x.name:x for x in load_registry().markers}
        if set(markers)-set(policies):raise ValueError(f'unknown panel markers {set(markers)-set(policies)}')
        mask,oids=R.rasterize_mask(Path(ph['geometry']),shape,1.0)
        oids=np.asarray(oids,dtype=str);n=len(oids)
        nuc=cells.nucleus_area_um2.reindex(oids).fillna(0).to_numpy()>0
        allowed=np.r_[False,nuc]
        # No-nucleus objects cannot influence compensation or any biological calculation.
        mask[~allowed[mask]]=0
        sizes=np.bincount(mask.ravel(),minlength=n+1)[1:].astype(float)
        band=find_boundaries(mask,mode='inner',connectivity=2,background=0)&(mask>0)
        band_idx=np.flatnonzero(band.ravel());del band
        flat=mask.ravel();mask_band=flat[band_idx]
        sums=np.zeros((n,len(markers)));edges=np.zeros_like(sums)
        print(f'{donor}: raster {shape}; {n:,} cells; {len(band_idx):,} boundary pixels',flush=True)
        for k,marker in enumerate(markers):
            img=R.read_channel(pages,k).ravel()
            sums[:,k]=np.bincount(flat,weights=img,minlength=n+1)[1:]
            edges[:,k]=np.bincount(mask_band,weights=img[band_idx],minlength=n+1)[1:]
            del img
            if (k+1)%10==0:print(f'{donor}: measured {k+1}/{len(markers)} channels',flush=True)
        tf.close();del flat,mask_band,band_idx
        contact=R.contact_matrix_gapbridge(mask,n,1);del mask
        raw=np.divide(sums,sizes[:,None],out=np.full_like(sums,np.nan),where=sizes[:,None]>0)
        alignment=[]
        for marker in ('DAPI','INS','Pan_Cytokeratin','CD3e','Vimentin','CD31'):
            if marker not in markers:continue
            k=markers.index(marker);col='feat__Cell__'+sanitize(names[k])
            if col not in cells:continue
            a=cells[col].reindex(oids).to_numpy(float);b=raw[:,k];ok=np.isfinite(a)&np.isfinite(b)
            corr=float(np.corrcoef(a[ok],b[ok])[0,1]);positive=ok&(a>0)
            ratio=float(np.median(b[positive]/a[positive]))
            alignment.append({'marker':marker,'n':int(ok.sum()),'pearson_r':corr,'median_ratio':ratio})
            if ok.sum()<10 or not np.isfinite(corr) or corr<.8 or not .25<=ratio<=4:
                raise ValueError(f'pixel/QuPath alignment failed: {alignment[-1]}')
        if len(alignment)<2:raise ValueError('insufficient channel alignment probes')
        alpha,mode,_=R.apply_no_spillover_mask(R.RedseaParams(),markers,1.,1)
        corrected,clamp,isolated=R.compensate(sums,edges,sizes,contact,mode,alpha,norm_form='donor')
        sparse.save_npz(out/'contact_matrix.npz',contact)
        pd.DataFrame({'object_id':oids,'raster_pixels':sizes}).to_parquet(out/'raster_ids.parquet',index=False)
        expr=pd.DataFrame(index=pd.Index(oids,name='object_id'))
        for k,marker in enumerate(markers):
            nuclear=marker in panel.NUCLEUS_MEASURED
            col='feat__'+('Nucleus' if nuclear else 'Cell')+'__'+sanitize(names[k])
            original=cells[col].reindex(oids).to_numpy(float) if col in cells else np.full(n,np.nan)
            values=original.copy() if nuclear else corrected[:,k].copy()
            values[~np.isfinite(original)|~nuc]=np.nan
            expr[marker]=values.astype(np.float32)
        for roi,frame in frames.items():
            dest=out/roi;dest.mkdir()
            meta_cols=['object_id','donor_id','serial','roi','x_um','y_um','nuc_x_um','nuc_y_um','area_um2',
                       'nucleus_area_um2','qc_analysis_eligible','exclusion_reasons']
            result=frame[meta_cols].merge(expr.reset_index(),on='object_id',validate='one_to_one')
            # Zero raster area is unavailable even when the vector geometry looked valid.
            pixels=pd.Series(sizes,index=oids).reindex(result.object_id).to_numpy()
            result.loc[pixels<=0,'qc_analysis_eligible']=False
            result.loc[pixels<=0,'exclusion_reasons']+='raster_empty;'
            result.to_parquet(dest/'expression.parquet',index=False)
            rows=[]
            for marker in markers:
                stat=usability.shape_stats(result.loc[result.qc_analysis_eligible,marker].to_numpy())
                rows.append({'marker':marker,'measurement':'raw_nucleus' if marker in panel.NUCLEUS_MEASURED else 'redsea_cell',**stat})
            pd.DataFrame(rows).to_csv(dest/'marker_usability.csv',index=False)
        write_json(out/'redsea.json',{'parameters':{'downsample':1,'edge_radius':0,'gap_bridge':1,'alpha':1.,'comp_mode':1,'norm_form':'donor'},
            'alignment':alignment,'markers':markers,'n_isolated':isolated,'clamp_fraction':dict(zip(markers,clamp.tolist())),
            'nuclear_representation':'unscaled raw QuPath nuclear means; no DAPI ratio',
            'raw_compartments':str(prepared),'unaligned_expression_preserved':True})
        record.update(status='validated_complete',finished_utc=now(),
                      validation={'channel_alignment':True,'UUID_join':True,'no_nucleus_excluded':True},
                      outputs=[identity(p) for p in out.rglob('*') if p.is_file() and p.name!='receipt.json'])
        write_json(out/'receipt.json',record);write_json(base/'current.json',record)
        print(f'{donor}: REDSEA complete',flush=True)
    except Exception as exc:
        record.update(status='failed',finished_utc=now(),error=str(exc),traceback=traceback.format_exc())
        write_json(out/'receipt.json',record);write_json(base/'last_failure.json',record);raise


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--manifest',required=True,type=Path)
    p.add_argument('--donor',required=True);p.add_argument('--out',type=Path,default=DEFAULT_OUT)
    a=p.parse_args();run(a.manifest,a.donor,a.out)
