#!/usr/bin/env python3
"""Read-only QuPath annotation export and UUID-joined tissue/geometry preparation."""
from __future__ import annotations
import argparse
import json
import os
from pathlib import Path
import re
import shutil
import subprocess
import sys

import numpy as np
import pandas as pd
import pyarrow as pa
import pyarrow.parquet as pq
import shapely
from shapely.affinity import scale
from shapely.geometry import shape

HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(HERE.parents[1]))
from phenocycler.artifacts import iter_qupath_geojson_features
from manifest import load_manifest, canonical_donor
from workflow import DEFAULT_OUT, identity, write_json, now, valid_receipt, sha
from read_pheno_tsv import plan_columns


def annotations(image,out):
    project=Path('/home/smith6jt/IO60panc2nd/project.qpproj')
    spec=json.loads(project.read_text())
    hits=[r for r in spec['images'] if r['imageName']==image]
    if len(hits)!=1:raise ValueError('ambiguous QuPath entry')
    entry=str(hits[0]['entryID']);tag=re.sub(r'[^A-Za-z0-9._-]+','_',image)
    sentinel=out/'manual_annotations'/f'{tag}.json'
    source=Path('/home/smith6jt/IO60panc2nd')/'data'/entry/'data.qpdata'
    before=source.stat()
    work=out/'project';work.mkdir(parents=True,exist_ok=True)
    shutil.copy2(project,work/'project.qpproj')
    for name in ('data','classifiers','scripts','resources'):
        p=work/name
        if not p.exists():p.symlink_to(project.parent/name,target_is_directory=True)
    env=os.environ.copy()
    env['JAVA_TOOL_OPTIONS']='-Xmx48g -Djava.util.prefs.userRoot=/home/smith6jt/QuPath/pipeline_prefs'
    env['LD_LIBRARY_PATH']='/blue/maigan/smith6jt/IO60panc2nd_v3/lib_stdcxx'+(':'+env['LD_LIBRARY_PATH'] if env.get('LD_LIBRARY_PATH') else '')
    script=project.parent/'scripts/v3/export_annotations.groovy'
    with (out/'annotation_export.log').open('w') as log:
        subprocess.run(['/home/smith6jt/QuPath/bin/QuPath','script','-p',str(work/'project.qpproj'),
              '-i',image,'-a',str(out),str(script)],env=env,stdout=log,stderr=subprocess.STDOUT,check=True)
    after=source.stat()
    if (before.st_size,before.st_mtime_ns)!=(after.st_size,after.st_mtime_ns):
        raise ValueError('QuPath data changed during read-only annotation export')
    record=json.loads(sentinel.read_text())
    if str(record['entry_id'])!=entry or record['image']!=image:raise ValueError('annotation export identity mismatch')
    return Path(record['all_geojson']),{'entry_id':entry,'qpdata':identity(source),'sentinel':record,'script':identity(script)}


def partition_points(xy, tissue, pln):
    """pLN union has priority inside overlapping Tissue; every UUID occurs once."""
    pts=shapely.points(xy)
    lymph=shapely.covers(pln,pts)
    panc=shapely.covers(tissue,pts)&~lymph
    return np.where(lymph,'pln_1',np.where(panc,'panc','outside'))


def prepare(manifest_path,donor,root):
    m=load_manifest(manifest_path)
    entries=[r for r in m['sections'] if r['design']=='same_slide' and r['donor_id']==canonical_donor(donor)]
    if len(entries)!=2:raise ValueError('requires one pancreas and one pLN record')
    ph=entries[0]['pheno'];measurement=Path(ph['measurements']);geo=Path(ph['geometry'])
    inputs={k:identity(ph[k]) for k in ('measurements','geometry')}
    if ph.get('annotations'):inputs['annotations']=identity(ph['annotations'])
    inputs['source']=identity(__file__)
    receipts={}
    for step in ('0','2','C_cells','C_geojson'):
        p=Path(ph['provenance'])/(step+'.json');r=json.loads(p.read_text())
        if r['status']!='pass' or not r['postcondition']['pass']:raise ValueError(f'upstream {step} did not pass')
        receipts[step]={'receipt':identity(p),'record':r}
    import hashlib
    fingerprint=hashlib.sha256(json.dumps(inputs,sort_keys=True).encode()).hexdigest()
    base=Path(root)/'protein'/f'donor_id={donor}'/'preparation';base.mkdir(parents=True,exist_ok=True)
    if valid_receipt(base/'current.json',fingerprint):
        print(f'{donor}: validated protein preparation reused',flush=True);return
    from uuid import uuid4
    out=base/'attempts'/(fingerprint[:16]+'-'+uuid4().hex[:8]);out.mkdir(parents=True)
    record={'status':'running','fingerprint':fingerprint,'inputs':inputs,'started_utc':now(),
            'upstream':receipts,'job_id':os.environ.get('SLURM_JOB_ID'),'attempt':str(out)}
    write_json(out/'receipt.json',record)
    try:
        if ph.get('annotations'):
            # A segmentation made at another site comes with that project's own Tissue and pLN
            # objects; they are read as delivered, never re-exported from this site's project.
            ann_path=Path(ph['annotations']);ann_meta={'source':'delivered','file':identity(ann_path)}
        else:
            image=Path(ph['image']).name+' - resolution #1'
            ann_path,ann_meta=annotations(image,out)
        px=float(ph['pixel_size_um'])
        anns=json.loads(ann_path.read_text())['features']
        geoms={}
        for f in anns:
            cls=(f['properties'].get('classification') or {}).get('name','')
            g=shapely.make_valid(scale(shape(f['geometry']),xfact=px,yfact=px,origin=(0,0)))
            geoms.setdefault(cls,[]).append(g)
        if not geoms.get('Tissue') or not geoms.get('pLN'):raise ValueError('explicit Tissue and pLN objects are required')
        tissue=shapely.union_all(geoms['Tissue']);pln=shapely.union_all(geoms['pLN'])
        header=pq.read_schema(measurement).names
        meta,features,_=plan_columns(header)
        extra=[c for c in header if c=='Cell: Solidity']
        df=pd.read_parquet(measurement,columns=list(meta)+list(features)+extra).rename(columns={**meta,**features,'Cell: Solidity':'cell_solidity'})
        df['object_id']=df.cell_id.astype(str)
        if df.object_id.duplicated().any():raise ValueError('duplicate QuPath UUIDs')
        table_index=pd.Index(df.object_id)
        seen=set();batches=[];writer=None;sample_px=[];sample_um=[]
        n_nuc=0
        try:
            for n,f in enumerate(iter_qupath_geojson_features(geo)):
                oid=str(f['id'])
                if oid in seen:raise ValueError('duplicate geometry UUID')
                seen.add(oid)
                cell=shapely.make_valid(shape(f['geometry']));nuclear=f.get('nucleusGeometry')
                nuc=shapely.make_valid(shape(nuclear)) if nuclear else None
                cent=cell.centroid
                if len(sample_px)<10000:
                    ix=table_index.get_indexer([oid])[0]
                    if ix<0:raise ValueError('unknown geometry UUID')
                    sample_px.append([cent.x,cent.y]);sample_um.append(df.iloc[ix][['x_um','y_um']].to_numpy(float))
                cell=scale(cell,xfact=px,yfact=px,origin=(0,0))
                if nuc is not None and not nuc.is_empty:
                    nuc=scale(nuc,xfact=px,yfact=px,origin=(0,0));nc=nuc.centroid;n_nuc+=1
                    nx,ny,na=nc.x,nc.y,nuc.area
                else:nx=ny=np.nan;na=0.
                batches.append({'object_id':oid,'wkb':cell.wkb,'nucleus_wkb':nuc.wkb if nuc is not None else None,
                                'nuc_x_um':nx,'nuc_y_um':ny,'nucleus_area_um2':na,'geometry_area_um2':cell.area})
                if len(batches)==10000:
                    table=pa.Table.from_pylist(batches)
                    if writer is None:writer=pq.ParquetWriter(out/'geometry.parquet',table.schema)
                    writer.write_table(table);batches=[]
                if (n+1)%100000==0:print(f'{donor}: {n+1:,} geometries',flush=True)
            if batches:
                table=pa.Table.from_pylist(batches)
                if writer is None:writer=pq.ParquetWriter(out/'geometry.parquet',table.schema)
                writer.write_table(table)
        finally:
            if writer is not None:writer.close()
        if seen!=set(table_index):raise ValueError('geometry and measurement UUID populations disagree')
        residual=np.linalg.norm(np.asarray(sample_px)*px-np.asarray(sample_um),axis=1)
        if np.median(residual)>5:raise ValueError('pixel calibration does not reproduce measurement centroids')
        gd=pd.read_parquet(out/'geometry.parquet',columns=['object_id','nuc_x_um','nuc_y_um','nucleus_area_um2','geometry_area_um2'])
        df=df.merge(gd,on='object_id',validate='one_to_one')
        df['roi']=partition_points(df[['x_um','y_um']].to_numpy(),tissue,pln)
        df['qc_analysis_eligible']=False;df['exclusion_reasons']=''
        for roi in ('panc','pln_1'):
            mask=df.roi==roi;sub=df[mask].copy()
            nuc_ok=np.isfinite(sub.nucleus_area_um2)&(sub.nucleus_area_um2>0)
            area=sub.area_um2.to_numpy(float);median=float(np.median(area[nuc_ok]))
            reasons=np.full(len(sub),'',object)
            rules={'no_nucleus':~nuc_ok.to_numpy(),'area_below_20':area<20,'area_above_4x_median':area>4*median,
                   'low_or_invalid_solidity':~np.isfinite(sub.cell_solidity)|(sub.cell_solidity<.7),
                   'invalid_geometry_area':~np.isfinite(sub.geometry_area_um2)|(sub.geometry_area_um2<=0)}
            for name,bad in rules.items():reasons[bad]=np.char.add(reasons[bad].astype(str),name+';')
            sub['qc_analysis_eligible']=reasons=='';sub['exclusion_reasons']=reasons
            sub['donor_id']=donor;sub['serial']=entries[0]['xenium_serial']
            dest=out/roi;dest.mkdir();sub.to_parquet(dest/'cells.parquet',index=False)
            sub.loc[~sub.qc_analysis_eligible,['object_id','exclusion_reasons']].to_parquet(dest/'exclusions.parquet',index=False)
            print(f'{donor}/{roi}: {len(sub):,} cells; {int(sub.qc_analysis_eligible.sum()):,} geometry eligible',flush=True)
        summary={'n_cells':len(df),'n_nuclei':n_nuc,'n_outside_tissue_objects':int((df.roi=='outside').sum()),
            'pixel_size_um':px,'median_calibration_residual_um':float(np.median(residual)),
            'tissue_routing':'pLN union first; Tissue minus pLN; one row per UUID',
            'pLN_object_count':len(geoms['pLN']),'tissue_pLN_overlap_um2':float(tissue.intersection(pln).area),
            'annotation_export':ann_meta,'marker_columns':features}
        write_json(out/'geometry_qc.json',summary)
        outputs=[identity(p) for p in out.rglob('*') if p.is_file() and not p.is_symlink()
                 and 'project' not in p.relative_to(out).parts and p.name!='receipt.json']
        record.update(status='validated_complete',finished_utc=now(),outputs=outputs,
                      validation={'uuid_population_exact':True,'nucleus_exclusion':True,'calibration':True,'tissue_partition_disjoint':True})
        write_json(out/'receipt.json',record);write_json(base/'current.json',record)
    except Exception as exc:
        import traceback
        record.update(status='failed',error=str(exc),traceback=traceback.format_exc(),finished_utc=now())
        write_json(out/'receipt.json',record);write_json(base/'last_failure.json',record);raise


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--manifest',required=True,type=Path)
    p.add_argument('--donor',required=True);p.add_argument('--out',type=Path,default=DEFAULT_OUT)
    a=p.parse_args();prepare(a.manifest,a.donor,a.out)
