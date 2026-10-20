#!/usr/bin/env python3
"""Cell-level protein -> RNA transfer for one registered four-slide section.

Rows are exactly the section's eligible RNA cells (the rna stage's h5ad obs), so every
downstream RNA stage sees the same population. Protein comes from the analysis-eligible
PhenoCycler cells of the whole image (the lane's nucleus-aware QC and the delivery's own, with a
nucleus outline) by area-weighted overlap in the validated registration frame -- the operator and
edge set are ``build_joined.py``'s, so continuous measurements and any later categorical transfer
rest on the same contributing UUIDs and overlap weights. The Xenium region decides which RNA cells
are rows; the overlap decides which protein cells contribute. A Xenium region is not a tissue:
6580's pancreas region also images its same-block lymph node, whose 61k cells were left without
protein while the join read the pancreas tissue alone. ``protein_tissue`` records which annotation
(Tissue -> ``pancreas``, pLN -> ``pLN``) a cell's contributors lie in, so no stage mixes the two.

The values are the full-dataset UniFORM alignment imported by ``import_alignment.py`` (maintainer
ruling 2026-10-18), in two unit systems that never share a column: ``feat__<marker>`` in the
25-donor pancreas alignment's units and ``feat__pln__<marker>`` in the pLN-tissue alignment's.
``protein_frame`` says which a Xenium cell's contributors carried. The lane's own REDSEA values
are a different measurement version and are not carried. A Xenium cell with no overlap has no
protein measurement: NaN, never zero. Calls are not transferred here -- the typing delivery
reuses ``overlap_edges.parquet`` with ``lift_calls.transfer_calls``.

Cell-level transfer needs every upstream receipt validated, including ``registration_recheck``
made against the very preparation read here. When registration, its recheck or the alignment is
scientifically unavailable, the join is recorded scientifically unavailable as the section's
current state -- section-level QC only -- so an earlier join cannot be read in its place.
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
import pyarrow.parquet as pq
import shapely

HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(HERE.parents[1]))
from phenocycler.integration.transform import Transform
from manifest import load_manifest,get_section
from workflow import DEFAULT_OUT,section_root,identity,write_json,now
from build_joined import xenium_polygons,overlap_edges,aggregate_features

#: Alignment frame -> the column prefix its values travel under. One frame per column, always.
FRAME_PREFIX={'pancreas_cohort':'feat__','pln_tissue':'feat__pln__'}
FRAME_SOURCE={'pancreas_cohort':('pancreas_markers',''),'pln_tissue':('pln_markers','pln__')}
#: The preparation's routing (pLN union first, then Tissue) -> the tissue a protein cell lies in.
TISSUE={'panc':'pancreas','pln_1':'pLN'}


def current(path,what):
    """An upstream receipt that is finished: validated, or scientifically unavailable with a reason."""
    r=json.loads(Path(path).read_text())
    if r.get('status') not in ('validated_complete','scientifically_unavailable'):
        raise ValueError(f'{what} is {r.get("status")}; cell-level transfer requires a finished upstream stage')
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


def protein_table(alignment):
    """Analysis-eligible protein cells of both tissues; only the frames present among them get columns."""
    meta=json.loads((alignment/'alignment.json').read_text())
    al=pd.concat([pq.ParquetFile(alignment/roi/'aligned.parquet').read().to_pandas().assign(protein_tissue=tissue)
                  for roi,tissue in TISSUE.items()],ignore_index=True)
    if al.object_id.astype(str).duplicated().any():raise ValueError('a protein cell is routed to both tissues')
    keep=al.analysis_eligible.to_numpy(bool)
    t=pd.DataFrame({'object_id':al.object_id.astype(str).to_numpy()[keep],
                    'alignment_frame':al.alignment_frame.to_numpy()[keep],
                    'protein_tissue':al.protein_tissue.to_numpy()[keep]})
    present=[f for f in FRAME_PREFIX if (t.alignment_frame==f).any()]
    columns={}
    for frame in present:
        key,source=FRAME_SOURCE[frame]
        for m in meta[key]:
            columns[FRAME_PREFIX[frame]+m]=al[source+m].to_numpy(float)[keep]
    t=pd.concat([t,pd.DataFrame(columns,index=t.index)],axis=1)
    markers={f:list(meta[FRAME_SOURCE[f][0]]) for f in present}
    return t,markers,int((~keep).sum())


def dominant_partner(xids,edges):
    """Largest shared area per Xenium cell, ties broken by UUID; share of overlapped area."""
    total=edges.groupby('cell_id',sort=False).overlap_area_um2.sum()
    w=edges.sort_values(['overlap_area_um2','object_id'],ascending=[False,True],kind='stable').drop_duplicates('cell_id').set_index('cell_id')
    w=w.reindex(xids)
    return w.object_id.to_numpy(),(w.overlap_area_um2/total.reindex(xids)).to_numpy(float)


def frame_of_cells(xi,edge_frames,n,labels=tuple(FRAME_PREFIX)):
    """Which frame (or tissue) contributed to each Xenium cell: the one, 'mixed', or '' without protein."""
    seen=np.stack([np.bincount(xi,weights=(edge_frames==f).astype(float),minlength=n)>0 for f in labels])
    k=seen.sum(axis=0);out=np.full(n,'',dtype=object)
    out[k>1]='mixed';one=k==1
    out[one]=np.asarray(labels,dtype=object)[seen[:,one].argmax(axis=0)]
    return out


def bounded(df,edges,prot,column,xi):
    """Area-weighted means lie inside their contributors' range, for one marker of one frame."""
    v=prot.set_index('object_id')[column].reindex(edges.object_id).to_numpy()
    ok=np.isfinite(v)
    if not ok.any():return True
    lo=pd.Series(v[ok]).groupby(xi[ok]).min();hi=pd.Series(v[ok]).groupby(xi[ok]).max();got=df[column].to_numpy()[lo.index]
    return bool(np.all((got>=lo.to_numpy()-1e-4*np.abs(lo.to_numpy())-1e-6)&(got<=hi.to_numpy()+1e-4*np.abs(hi.to_numpy())+1e-6)))


def run(manifest_path,sid,root):
    section=get_section(load_manifest(manifest_path),sid);donor=section['donor_id'];roi=section['roi']
    base=section_root(root,section);root=Path(root)
    receipts={'registration':base/'registration/current.json','rna':base/'rna/current.json',
              'preparation':root/'protein'/f'donor_id={donor}'/'preparation/current.json',
              'alignment':root/'protein'/f'donor_id={donor}'/'alignment/current.json',
              'registration_recheck':base/'registration_recheck/current.json'}
    rec={k:current(p,k) for k,p in receipts.items()}
    # The transform's nuclei check must be about the preparation this join reads, not an older one.
    if rec['registration_recheck']['inputs']['preparation']['sha256']!=identity(receipts['preparation'])['sha256']:
        raise ValueError('registration_recheck was made against another preparation; rerun recheck_nuclei.py')
    inputs={'section':section,'receipts':{k:identity(p) for k,p in receipts.items()},
            'code':[identity(HERE/'join_four.py'),identity(HERE/'build_joined.py')]}
    fingerprint=hashlib.sha256(json.dumps(inputs,sort_keys=True).encode()).hexdigest()
    out=base/'joined'/'attempts'/(fingerprint[:16]+'-'+uuid4().hex[:8]);out.mkdir(parents=True)
    report={'status':'running','section_id':sid,'fingerprint':fingerprint,'inputs':inputs,
            'started_utc':now(),'attempt':str(out),'job_id':os.environ.get('SLURM_JOB_ID')}
    write_json(out/'receipt.json',report)
    unavailable={k:r.get('reason','') for k,r in rec.items() if r['status']=='scientifically_unavailable'}
    if unavailable:
        # The section keeps section-level QC only. Recorded as current so no earlier join -- one
        # built on a superseded measurement among them -- can still be read as this section's.
        write_json(out/'unavailable.json',unavailable)
        report.update(status='scientifically_unavailable',reason='; '.join(f'{k}: {v}' for k,v in unavailable.items()),
                      finished_utc=now(),validation={'upstream_unavailable':sorted(unavailable)},
                      outputs=[identity(out/'unavailable.json')])
        write_json(out/'receipt.json',report);write_json(base/'joined/current.json',report)
        print(f'{sid}: join unavailable: {report["reason"]}',flush=True);return
    try:
        reg=rec['registration'];tf=Transform.load(Path(reg['transform_directory']),stem=reg['transform_stem'])
        rna_ids=eligible_rna_ids(Path(rec['rna']['attempt'])/'rna.h5ad')
        xids,xpolys=xenium_polygons(Path(section['bundle']))
        xids=pd.Index(np.asarray(xids).astype(str))
        missing=rna_ids.difference(xids)
        if len(missing):raise ValueError(f'{len(missing)} eligible RNA cells have no Xenium boundary')
        pos=xids.get_indexer(rna_ids);xpolys=xpolys[pos];xids=rna_ids
        alignment=Path(rec['alignment']['attempt'])
        prot,markers,n_protein_excluded=protein_table(alignment)
        geo=pd.read_parquet(Path(rec['preparation']['attempt'])/'geometry.parquet',columns=['object_id','wkb','nucleus_wkb'])
        geo=geo.assign(object_id=geo.object_id.astype(str))
        geo=geo[geo.object_id.isin(prot.object_id)&geo.nucleus_wkb.notna()][['object_id','wkb']]
        prot=prot[prot.object_id.isin(geo.object_id)]
        edges,xarea=overlap_edges(xids.to_numpy(),xpolys,geo,tf)
        if edges.duplicated(['cell_id','object_id']).any():raise ValueError('duplicate overlap edges')
        if not edges.object_id.isin(prot.object_id).all():raise ValueError('edge to an ineligible protein cell')
        by_id=prot.set_index('object_id')
        edges['alignment_frame']=edges.object_id.map(by_id.alignment_frame).to_numpy()
        edges['protein_tissue']=edges.object_id.map(by_id.protein_tissue).to_numpy()
        df,_=aggregate_features(xids.to_numpy(),edges,prot.drop(columns=['alignment_frame','protein_tissue']))
        xi=xids.get_indexer(edges.cell_id);n=len(xids)
        covered=np.bincount(xi,weights=edges.overlap_area_um2.to_numpy(float),minlength=n)
        cent=shapely.centroid(xpolys)
        df['x_um']=shapely.get_x(cent);df['y_um']=shapely.get_y(cent)
        df['area_um2']=shapely.area(xpolys);df['registered_area_um2']=xarea
        df['n_pheno_cells']=np.bincount(xi,minlength=n)
        df['overlap_fraction']=np.divide(covered,xarea,out=np.full(n,np.nan),where=xarea>0)
        df['dominant_object_id'],df['dominant_share']=dominant_partner(xids,edges)
        df['protein_frame']=frame_of_cells(xi,edges.alignment_frame.to_numpy(),n)
        df['protein_tissue']=frame_of_cells(xi,edges.protein_tissue.to_numpy(),n,tuple(TISSUE.values()))
        df['protein_available']=df.n_pheno_cells>0
        df.insert(1,'donor_id',donor);df.insert(2,'roi',roi);df.insert(3,'serial',section['xenium_serial'])
        df.to_parquet(out/'joined_cells.parquet',index=False)
        edges.to_parquet(out/'overlap_edges.parquet',index=False)
        features=[c for c in df if c.startswith('feat__')]
        within={f:bounded(df,edges,prot,FRAME_PREFIX[f]+ms[0],xi) for f,ms in markers.items() if ms}
        no_overlap_nan=bool(df.loc[~df.protein_available,features].isna().all().all())
        frames=json.loads((alignment/'frames.json').read_text())
        summary={'n_rna_eligible_cells':n,'n_protein_eligible_cells':len(prot),'n_protein_cells_excluded':n_protein_excluded,
                 'n_edges':len(edges),'cells_with_protein':int(df.protein_available.sum()),
                 'fraction_cells_with_protein':float(df.protein_available.mean()),
                 'protein_cells_touched':int(edges.object_id.nunique()),
                 'cells_by_protein_frame':pd.Series(df.protein_frame).replace('','none').value_counts().to_dict(),
                 'cells_by_protein_tissue':pd.Series(df.protein_tissue).replace('','none').value_counts().to_dict(),
                 'cells_by_frame_and_tissue':{f'{f or "none"}|{t or "none"}':int(c) for (f,t),c in
                                              df.groupby(['protein_frame','protein_tissue']).size().items()},
                 'overlap_fraction_median':float(np.nanmedian(df.overlap_fraction)),
                 'overlap_fraction_gt_1_01':int((df.overlap_fraction>1.01).sum()),
                 'dominant_share_median':float(np.nanmedian(df.dominant_share)),
                 'n_pheno_cells_median':float(np.median(df.n_pheno_cells[df.protein_available])) if df.protein_available.any() else None,
                 'markers_by_frame':markers,'frame_prefix':{f:FRAME_PREFIX[f] for f in markers},
                 'alignment_notes':frames['notes'],
                 'transform':{'directory':reg['transform_directory'],'stem':reg['transform_stem'],'verdict':reg.get('verdict'),'reason':reg.get('reason')}}
        write_json(out/'join.json',summary)
        checks={'rows_equal_rna_eligible':len(df)==n and df.cell_id.astype(str).equals(pd.Series(xids)),
                'edges_unique_positive':bool(np.isfinite(edges.overlap_area_um2).all() and (edges.overlap_area_um2>0).all()),
                'edges_eligible_protein_only':True,'no_overlap_is_missing_not_zero':no_overlap_nan,
                'weighted_mean_within_contributors':all(within.values()) and bool(within),
                'one_frame_per_column':True,'registration_validated':True}
        ok=all(checks.values()) and summary['cells_with_protein']>0
        report.update(status='validated_complete' if ok else 'failed',validation=checks,
                      summary={k:v for k,v in summary.items() if k not in('markers_by_frame','alignment_notes')},
                      finished_utc=now(),outputs=[identity(p) for p in out.iterdir() if p.is_file() and p.name!='receipt.json'])
        write_json(out/'receipt.json',report)
        if not ok:write_json(base/'joined/last_failure.json',report);raise ValueError(f'validation failed: {checks}')
        write_json(base/'joined/current.json',report)
        print(f"{sid}: joined {summary['cells_with_protein']:,}/{n:,} cells with protein {summary['cells_by_frame_and_tissue']}; "
              f"median overlap {summary['overlap_fraction_median']:.2f}",flush=True)
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
