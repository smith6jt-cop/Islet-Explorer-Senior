#!/usr/bin/env python3
"""Cell-level protein -> RNA transfer for one registered four-slide section.

Rows are exactly the section's eligible RNA cells (the rna stage's h5ad obs), so every
downstream RNA stage sees the same population. Protein comes from the eligible PhenoCycler
cells of the same tissue (``qc_analysis_eligible`` with a nucleus outline) by area-weighted
polygon overlap in the validated registration frame -- the operator and edge set are
``build_joined.py``'s, so continuous measurements and any later categorical transfer rest on
the same contributing UUIDs and overlap weights.

Both the UniFORM-aligned measurement (``feat__<marker>``) and the unaligned REDSEA
measurement (``feat__redsea__<marker>``) are carried; the calibration's per-marker alignment
status travels in join.json. A Xenium cell with no overlap has no protein measurement:
NaN, never zero. Calls are not transferred here -- the protein typing stage reuses
``overlap_edges.parquet`` with ``lift_calls.transfer_calls``.

Refuses any section whose registration receipt is not validated_complete: failed
registration permits section-level QC only.
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
import shapely

HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(HERE.parents[1]))
from phenocycler.integration.transform import Transform
from manifest import load_manifest,get_section
from workflow import DEFAULT_OUT,section_root,identity,write_json,now
from build_joined import xenium_polygons,overlap_edges,aggregate_features

META=('object_id','donor_id','serial','roi','x_um','y_um','nuc_x_um','nuc_y_um','area_um2',
      'nucleus_area_um2','qc_analysis_eligible','exclusion_reasons')


def current(path,what):
    r=json.loads(Path(path).read_text())
    if r.get('status')!='validated_complete':
        raise ValueError(f'{what} is {r.get("status")}; cell-level transfer requires validated_complete')
    return r


def eligible_rna_ids(h5ad):
    import anndata as ad
    a=ad.read_h5ad(h5ad,backed='r')
    obs=a.obs
    if not obs.qc_eligible.all() or (obs.nucleus_area<=0).any():raise ValueError('RNA rows are not nucleus-aware QC eligible')
    ids=pd.Index(a.obs_names.astype(str))
    a.file.close()
    if not ids.is_unique:raise ValueError('nonunique RNA cell ids')
    return ids


def protein_table(aligned_path,redsea_path):
    al=pd.read_parquet(aligned_path);rs=pd.read_parquet(redsea_path)
    markers=[c for c in al.columns if c not in META]
    if [c for c in rs.columns if c not in META]!=markers:raise ValueError('aligned and REDSEA marker sets differ')
    if not al.object_id.astype(str).equals(rs.object_id.astype(str)):raise ValueError('aligned and REDSEA rows differ')
    keep=al.qc_analysis_eligible.to_numpy(bool)&(al.nucleus_area_um2.to_numpy(float)>0)
    t=pd.DataFrame({'object_id':al.object_id.astype(str).to_numpy()[keep]})
    for m in markers:
        t['feat__'+m]=al[m].to_numpy(float)[keep]
        t['feat__redsea__'+m]=rs[m].to_numpy(float)[keep]
    return t,markers,int((~keep).sum())


def dominant_partner(xids,edges):
    """Largest shared area per Xenium cell, ties broken by UUID; share of overlapped area."""
    total=edges.groupby('cell_id',sort=False).overlap_area_um2.sum()
    w=edges.sort_values(['overlap_area_um2','object_id'],ascending=[False,True],kind='stable').drop_duplicates('cell_id').set_index('cell_id')
    w=w.reindex(xids)
    return w.object_id.to_numpy(),(w.overlap_area_um2/total.reindex(xids)).to_numpy(float)


def run(manifest_path,sid,root):
    section=get_section(load_manifest(manifest_path),sid);donor=section['donor_id'];roi=section['roi']
    base=section_root(root,section);root=Path(root)
    receipts={'registration':base/'registration/current.json','rna':base/'rna/current.json',
              'preparation':root/'protein'/f'donor_id={donor}'/'preparation/current.json',
              'redsea':root/'protein'/f'donor_id={donor}'/'redsea/current.json',
              'calibration':root/'protein'/f'calibration_{roi}'/'current.json'}
    rec={k:current(p,k) for k,p in receipts.items()}
    inputs={'section':section,'receipts':{k:identity(p) for k,p in receipts.items()},
            'code':[identity(HERE/'join_four.py'),identity(HERE/'build_joined.py')]}
    fingerprint=hashlib.sha256(json.dumps(inputs,sort_keys=True).encode()).hexdigest()
    out=base/'joined'/'attempts'/(fingerprint[:16]+'-'+uuid4().hex[:8]);out.mkdir(parents=True)
    report={'status':'running','section_id':sid,'fingerprint':fingerprint,'inputs':inputs,
            'started_utc':now(),'attempt':str(out),'job_id':os.environ.get('SLURM_JOB_ID')}
    write_json(out/'receipt.json',report)
    try:
        reg=rec['registration'];tf=Transform.load(Path(reg['transform_directory']),stem=reg['transform_stem'])
        rna_ids=eligible_rna_ids(Path(rec['rna']['attempt'])/'rna.h5ad')
        xids,xpolys=xenium_polygons(Path(section['bundle']))
        xids=pd.Index(np.asarray(xids).astype(str))
        missing=rna_ids.difference(xids)
        if len(missing):raise ValueError(f'{len(missing)} eligible RNA cells have no Xenium boundary')
        pos=xids.get_indexer(rna_ids);xpolys=xpolys[pos];xids=rna_ids
        cal=Path(rec['calibration']['attempt']);red=Path(rec['redsea']['attempt'])
        prot,markers,n_protein_excluded=protein_table(cal/f'donor_id={donor}'/'aligned.parquet',red/roi/'expression.parquet')
        geo=pd.read_parquet(Path(rec['preparation']['attempt'])/'geometry.parquet',columns=['object_id','wkb','nucleus_wkb'])
        geo=geo.assign(object_id=geo.object_id.astype(str))
        geo=geo[geo.object_id.isin(prot.object_id)&geo.nucleus_wkb.notna()][['object_id','wkb']]
        prot=prot[prot.object_id.isin(geo.object_id)]
        edges,xarea=overlap_edges(xids.to_numpy(),xpolys,geo,tf)
        if edges.duplicated(['cell_id','object_id']).any():raise ValueError('duplicate overlap edges')
        if not edges.object_id.isin(prot.object_id).all():raise ValueError('edge to an ineligible protein cell')
        df,_=aggregate_features(xids.to_numpy(),edges,prot)
        xi=xids.get_indexer(edges.cell_id);n=len(xids)
        covered=np.bincount(xi,weights=edges.overlap_area_um2.to_numpy(float),minlength=n)
        cent=shapely.centroid(xpolys)
        df['x_um']=shapely.get_x(cent);df['y_um']=shapely.get_y(cent)
        df['area_um2']=shapely.area(xpolys);df['registered_area_um2']=xarea
        df['n_pheno_cells']=np.bincount(xi,minlength=n)
        df['overlap_fraction']=np.divide(covered,xarea,out=np.full(n,np.nan),where=xarea>0)
        df['dominant_object_id'],df['dominant_share']=dominant_partner(xids,edges)
        df['protein_available']=df.n_pheno_cells>0
        df.insert(1,'donor_id',donor);df.insert(2,'roi',roi);df.insert(3,'serial',section['xenium_serial'])
        df.to_parquet(out/'joined_cells.parquet',index=False)
        edges.to_parquet(out/'overlap_edges.parquet',index=False)
        # Area-weighted means must lie inside their contributors' range; check one marker exactly.
        m='feat__'+markers[0];v=prot.set_index('object_id')[m].reindex(edges.object_id).to_numpy()
        lo=pd.Series(v).groupby(xi).min();hi=pd.Series(v).groupby(xi).max();got=df[m].to_numpy()[lo.index]
        bounded=bool(np.all((got>=lo.to_numpy()-1e-4*np.abs(lo.to_numpy())-1e-6)&(got<=hi.to_numpy()+1e-4*np.abs(hi.to_numpy())+1e-6)))
        no_overlap_nan=bool(df.loc[~df.protein_available,[c for c in df if c.startswith('feat__')]].isna().all().all())
        touched=edges.object_id.nunique()
        shifts=json.loads((cal/'shifts.json').read_text())
        summary={'n_rna_eligible_cells':n,'n_protein_eligible_cells':len(prot),'n_protein_cells_excluded':n_protein_excluded,
                 'n_edges':len(edges),'cells_with_protein':int(df.protein_available.sum()),
                 'fraction_cells_with_protein':float(df.protein_available.mean()),
                 'protein_cells_touched':int(touched),
                 'overlap_fraction_median':float(np.nanmedian(df.overlap_fraction)),
                 'overlap_fraction_gt_1_01':int((df.overlap_fraction>1.01).sum()),
                 'dominant_share_median':float(np.nanmedian(df.dominant_share)),
                 'n_pheno_cells_median':float(np.median(df.n_pheno_cells[df.protein_available])) if df.protein_available.any() else None,
                 'markers':markers,'marker_alignment':{k:{'status':v.get('status'),'cross_donor_cutoff_supported':v.get('cross_donor_cutoff_supported')} for k,v in shifts.items()},
                 'transform':{'directory':reg['transform_directory'],'stem':reg['transform_stem'],'verdict':reg.get('verdict'),'reason':reg.get('reason')}}
        write_json(out/'join.json',summary)
        checks={'rows_equal_rna_eligible':len(df)==n and df.cell_id.astype(str).equals(pd.Series(xids)),
                'edges_unique_positive':bool(np.isfinite(edges.overlap_area_um2).all() and (edges.overlap_area_um2>0).all()),
                'edges_eligible_protein_only':True,'no_overlap_is_missing_not_zero':no_overlap_nan,
                'weighted_mean_within_contributors':bounded,'registration_validated':True}
        ok=all(checks.values()) and summary['cells_with_protein']>0
        report.update(status='validated_complete' if ok else 'failed',validation=checks,summary={k:v for k,v in summary.items() if k not in('markers','marker_alignment')},
                      finished_utc=now(),outputs=[identity(p) for p in out.iterdir() if p.is_file() and p.name!='receipt.json'])
        write_json(out/'receipt.json',report)
        if not ok:write_json(base/'joined/last_failure.json',report);raise ValueError(f'validation failed: {checks}')
        write_json(base/'joined/current.json',report)
        print(f"{sid}: joined {summary['cells_with_protein']:,}/{n:,} cells with protein; median overlap {summary['overlap_fraction_median']:.2f}",flush=True)
    except Exception as exc:
        if report['status']!='failed':
            report.update(status='failed',error=str(exc),traceback=traceback.format_exc(),finished_utc=now())
            write_json(out/'receipt.json',report);write_json(base/'joined/last_failure.json',report)
        raise


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__,formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('--manifest',required=True,type=Path);p.add_argument('--section',required=True)
    p.add_argument('--out',type=Path,default=DEFAULT_OUT)
    a=p.parse_args();run(a.manifest,a.section,a.out)
