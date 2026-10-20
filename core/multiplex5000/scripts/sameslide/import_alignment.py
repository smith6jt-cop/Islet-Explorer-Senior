#!/usr/bin/env python3
"""Full-dataset UniFORM values on this lane's own objects, one donor at a time.

Maintainer ruling 2026-10-18: protein here is the cohort's own-compartment REDSEA measurement,
aligned by the full-dataset UniFORM run on the workstation. It replaces the lane's own UniFORM fit
over these four acquisitions alone (calibrate_protein.py, removed; git has it). Two deliveries
(manifest ``protein_alignment``), two unit systems:

* ``cohort``    -- the 25-donor pancreas alignment. Same-block tissue only: a different-block
  node is NaN there with ``in_alignment_block`` False.
* ``pln_nodes`` -- the different-block nodes, aligned among pLN tissue on their shared markers.

They are never put on one scale. Pancreas-frame markers keep their names, pLN-frame markers are
written ``pln__<marker>``, and a cell carries finite values in its own frame only.

Every delivered file read is checked against its delivery's SHA256SUMS. The delivered objects must
be THIS lane's objects -- every delivered UUID a prepared UUID -- or the donor was aligned on
another segmentation of the same image and is recorded scientifically_unavailable; a centroid
match is not a substitute. Joined by object_id only. A prepared cell the delivery lacks has no
aligned value: NaN, never zero.
"""
import argparse
import hashlib
import json
import os
from pathlib import Path
import traceback
from uuid import uuid4

import numpy as np
import pandas as pd
import pyarrow.parquet as pq

from manifest import load_manifest,canonical_donor
from workflow import DEFAULT_OUT,identity,write_json,now,valid_receipt

META=('object_id','donor_id','serial','roi','x_um','y_um','nuc_x_um','nuc_y_um','area_um2',
      'nucleus_area_um2','qc_analysis_eligible','exclusion_reasons')
COHORT_META=('donor_id','object_id','qc_analysis_eligible','X_centroid','Y_centroid','cell_region',
             'islet_num','in_alignment_block')
NODE_META=('object_id','X_centroid','Y_centroid','qc_analysis_eligible_with_nucleus','sample_id',
           'donor_id','source','annotation_uuid')
PLN_PREFIX='pln__'
COHORT_FILES=('generation_config.json','shifts.json','scale.parquet','scale_provenance.json',
              'marker_donor_usability.parquet','unavailable_alignment_histograms.json')
NODE_FILES=('provenance.json','shifts.json','scale.parquet','marker_alignment_QC.csv','uniform_pln_nodes.py')


def read(path,columns=None):
    """One parquet file, never a hive dataset: a donor_id=<d> path segment is not a column."""
    return pq.ParquetFile(path).read(columns=columns).to_pandas()


def checksums(root):
    """A delivery's SHA256SUMS as relative path -> digest."""
    out={}
    for line in (Path(root)/'SHA256SUMS').read_text().splitlines():
        if line.strip():
            digest,name=line.split(None,1)
            out[name.strip().lstrip('*').removeprefix('./')]=digest
    return out


def verified(root,rel,sums):
    """identity() of a delivered file, refused unless it matches its delivery checksum."""
    record=identity(Path(root)/rel)
    if sums.get(rel)!=record['sha256']:
        raise ValueError(f'{Path(root)/rel} does not match its delivery checksum')
    return record


def delivery_donor(donors,donor):
    """The delivery's spelling of a canonical donor id ('0112' for '112'), found, never assumed."""
    hits=[d for d in donors if canonical_donor(d)==donor]
    if len(hits)!=1:raise ValueError(f'{donor}: expected one delivery donor, found {hits}')
    return hits[0]


def frame_table(lane,cohort,node,pan_markers,pln_markers):
    """Lane rows joined to both frames by object_id; finite values in the cell's own frame only."""
    c=cohort.rename(columns={'qc_analysis_eligible':'cohort_qc_eligible'})[
        ['object_id','cohort_qc_eligible','cell_region','islet_num','in_alignment_block',*pan_markers]]
    n=node.rename(columns={'qc_analysis_eligible_with_nucleus':'node_qc_eligible',
                           **{m:PLN_PREFIX+m for m in pln_markers}})[
        ['object_id','node_qc_eligible','annotation_uuid',*[PLN_PREFIX+m for m in pln_markers]]]
    t=lane.merge(c,on='object_id',how='left',validate='one_to_one').merge(n,on='object_id',how='left',validate='one_to_one')
    in_pan=t.in_alignment_block.astype('boolean').fillna(False).to_numpy(bool)
    in_pln=t.annotation_uuid.notna().to_numpy()
    if (in_pan&in_pln).any():raise ValueError('a cell is in both alignment frames')
    t['alignment_frame']=np.where(in_pan,'pancreas_cohort',np.where(in_pln,'pln_tissue',''))
    # The lane's nucleus-aware QC AND the delivery's own, in the cell's frame.
    frame_ok=np.where(in_pan,t.cohort_qc_eligible.astype('boolean').fillna(False).to_numpy(bool),
                      np.where(in_pln,t.node_qc_eligible.astype('boolean').fillna(False).to_numpy(bool),False))
    t['analysis_eligible']=(t.qc_analysis_eligible.to_numpy(bool)&(t.nucleus_area_um2.to_numpy(float)>0)&frame_ok)
    for m in [*pan_markers,*[PLN_PREFIX+m for m in pln_markers]]:
        t[m]=t[m].astype(np.float32)
    return t


def finite(v):
    """A JSON-safe number. Unavailable -- NaN, or the delivery's own 'unavailable' -- is null."""
    try:
        v=float(v)
    except (TypeError,ValueError):
        return None
    return v if np.isfinite(v) else None


def marker_record(shifts,scale,key,markers):
    """Per marker: this sample's scale, the alignment's verdict and the cohort-wide floor/top/span."""
    out={}
    for m in markers:
        s=shifts[m];row=scale.loc[m] if m in scale.index else None
        out[m]={'scale':finite(s.get('scale',{}).get(key)),'provenance':s.get('provenance'),
                'reference':s.get('reference_donor',s.get('reference_sample')),'over_bar':bool(s.get('over_bar')),
                'mode_sd_after':finite(s.get('mode_sd_after')),'acq_ratio':finite(s.get('acq_ratio')),
                **{k:None if row is None else finite(row[k]) for k in ('floor','top','span')}}
    return out


def run(manifest_path,donor,root):
    m=load_manifest(manifest_path);donor=canonical_donor(donor)
    sections=[s for s in m['sections'] if s['design']=='same_slide' and s['donor_id']==donor]
    if sorted(s['roi'] for s in sections)!=['panc','pln_1']:raise ValueError('requires one pancreas and one pLN record')
    src=m['protein_alignment'];cohort,nodes=Path(src['cohort']),Path(src['pln_nodes'])
    prep_path=Path(root)/'protein'/f'donor_id={donor}'/'preparation/current.json'
    if not valid_receipt(prep_path):raise ValueError('protein preparation is not validated')
    prep=Path(json.loads(prep_path.read_text())['attempt'])
    inputs={'preparation':identity(prep_path),'code':identity(__file__),'ruling':src['ruling'],
            'cohort':{'root':str(cohort),'readme':identity(cohort/'README.md'),'checksums':identity(cohort/'SHA256SUMS')},
            'pln_nodes':{'root':str(nodes),'readme':identity(nodes/'README.md'),'checksums':identity(nodes/'SHA256SUMS')}}
    fingerprint=hashlib.sha256(json.dumps(inputs,sort_keys=True).encode()).hexdigest()
    base=Path(root)/'protein'/f'donor_id={donor}'/'alignment';base.mkdir(parents=True,exist_ok=True)
    if valid_receipt(base/'current.json',fingerprint):print(f'{donor}: alignment import reused',flush=True);return
    out=base/'attempts'/(fingerprint[:16]+'-'+uuid4().hex[:8]);out.mkdir(parents=True)
    rec={'status':'running','donor_id':donor,'fingerprint':fingerprint,'inputs':inputs,'started_utc':now(),
         'attempt':str(out),'job_id':os.environ.get('SLURM_JOB_ID')}
    write_json(out/'receipt.json',rec)
    try:
        csum,nsum=checksums(cohort),checksums(nodes)
        cshift=json.loads((cohort/'shifts.json').read_text());nshift=json.loads((nodes/'shifts.json').read_text())
        dd=delivery_donor(cshift['donors'],donor)
        node_rel=sorted(str(p.relative_to(nodes)) for p in (nodes/'aligned').glob(f'node__{dd}__*.parquet'))
        files={'cohort':{r:verified(cohort,r,csum) for r in (*COHORT_FILES,f'uniform/donor_id={dd}/data.parquet',
                                                              f'availability/donor_id={dd}/data.parquet')},
               'pln_nodes':{r:verified(nodes,r,nsum) for r in (*NODE_FILES,*node_rel)}}
        rec['delivery_files']=files;rec['delivery_donor_id']=dd
        coh=read(cohort/f'uniform/donor_id={dd}/data.parquet')
        if set(coh.donor_id.astype(str))!={dd}:raise ValueError('cohort partition holds another donor')
        node=(pd.concat([read(nodes/r) for r in node_rel],ignore_index=True) if node_rel
              else pd.DataFrame({c:pd.Series(dtype=object) for c in NODE_META}))
        if len(node) and set(node.donor_id.astype(str))!={dd}:raise ValueError('node file holds another donor')
        for name,t in (('cohort',coh),('pln_nodes',node)):
            if t.object_id.isna().any() or t.object_id.astype(str).duplicated().any():raise ValueError(f'{name}: null or duplicate UUIDs')
            t['object_id']=t.object_id.astype(str)
        prepared=set(read(prep/'geometry.parquet',['object_id']).object_id.astype(str))
        delivered=set(coh.object_id)|set(node.object_id)
        relation={'delivery_donor_id':dd,'prepared_objects':len(prepared),'cohort_rows':len(coh),'node_rows':len(node),
                  'delivered_in_prepared':len(delivered&prepared),'delivered_not_prepared':len(delivered-prepared),
                  'prepared_not_delivered':len(prepared-delivered)}
        write_json(out/'uuid_relation.json',relation)
        if relation['delivered_not_prepared']:
            reason=(f"{relation['delivered_not_prepared']:,} of {len(delivered):,} delivered UUIDs are not this lane's "
                    f"objects ({relation['delivered_in_prepared']:,} shared): the alignment was made on another "
                    f"segmentation of this image. Values attach by object_id only.")
            rec.update(status='scientifically_unavailable',reason=reason,finished_utc=now(),
                       validation={'delivery_checksums_verified':True,'delivered_objects_are_prepared_objects':False},
                       outputs=[identity(out/'uuid_relation.json')])
            write_json(out/'receipt.json',rec);write_json(base/'current.json',rec)
            print(f'{donor}: alignment unavailable: {reason}',flush=True);return
        pan_markers=[c for c in coh.columns if c not in COHORT_META]
        pln_markers=[c for c in node.columns if c not in NODE_META]
        block=coh.in_alignment_block.astype(bool)
        # The delivery's own statement, checked rather than trusted: a node cell is an unaligned
        # (NaN, out-of-block) cohort row, and no out-of-block row carries a value.
        nodec=coh[coh.object_id.isin(set(node.object_id))]
        node_unaligned=bool(len(nodec)==len(node) and not nodec.in_alignment_block.any())
        out_of_block_empty=bool(coh.loc[~block,pan_markers].isna().all().all())
        if not (node_unaligned and out_of_block_empty):raise ValueError('frames overlap in the delivery')
        frames={'pancreas_cohort':{'source':str(cohort),'units':'cohort own-compartment REDSEA x per-donor UniFORM scale',
                                   'column_prefix':'','markers':marker_record(cshift['shifts'],read(cohort/'scale.parquet').set_index('marker'),dd,pan_markers)},
                'pln_tissue':{'source':str(nodes),'units':'raw QuPath means x per-sample UniFORM scale among pLN tissue',
                              'column_prefix':PLN_PREFIX,'samples':sorted(node.sample_id.astype(str).unique().tolist()) if len(node) else [],
                              'markers':{}}}
        for sample in frames['pln_tissue']['samples']:
            frames['pln_tissue']['markers'][sample]=marker_record(nshift,read(nodes/'scale.parquet').set_index('marker'),sample,pln_markers)
        frames['notes']={'CD8':'aligned CD8 carries UniFORM\'s shift, the method for all 25 donors (maintainer ruling 2026-10-19). '
                                'Not ruled: a within-donor step. The factors are per marker, so on these four donors CD8 is '
                                'multiplied 3.7-19x as much as CD4, and a call that sets CD8 against CD4 in one cell is exposed; '
                                'a within-section rank is not',
                         'pln_over_bar':'pLN-frame markers over the 0.10 mode-SD bar are not shown comparable across samples by this alignment',
                         'never_pooled':'pancreas_cohort and pln_tissue values are different units and never share a scale'}
        write_json(out/'frames.json',frames)
        summary={}
        for roi in ('panc','pln_1'):
            lane=read(prep/roi/'cells.parquet',list(META));lane['object_id']=lane.object_id.astype(str)
            t=frame_table(lane,coh,node,pan_markers,pln_markers)
            region=t.cell_region
            known=region.notna().to_numpy()
            disagree=int((known&((region.astype(str)=='pLN').to_numpy()!=(roi=='pln_1'))).sum())
            if disagree:raise ValueError(f'{roi}: {disagree:,} cells routed against the delivery\'s Parent')
            missing=~t.alignment_frame.astype(bool).to_numpy()
            values=[*pan_markers,*[PLN_PREFIX+m for m in pln_markers]]
            if t.loc[missing,values].notna().any().any():raise ValueError('a cell without a delivered row carries a value')
            dest=out/roi;dest.mkdir();t.to_parquet(dest/'aligned.parquet',index=False)
            summary[roi]={'rows':len(t),'lane_eligible':int(t.qc_analysis_eligible.sum()),
                          'analysis_eligible':int(t.analysis_eligible.sum()),
                          'by_frame':t.alignment_frame.replace('','not_delivered').value_counts().to_dict(),
                          'eligible_by_frame':t.loc[t.analysis_eligible,'alignment_frame'].value_counts().to_dict(),
                          'lane_eligible_frame_ineligible':int((t.qc_analysis_eligible&~t.analysis_eligible&~missing).sum()),
                          'lane_eligible_not_delivered':int((t.qc_analysis_eligible&missing).sum()),
                          'cell_region':region.fillna('not_delivered').value_counts().to_dict()}
            print(f"{donor}/{roi}: {summary[roi]['analysis_eligible']:,} of {len(t):,} cells analysis-eligible; frames {summary[roi]['eligible_by_frame']}",flush=True)
        write_json(out/'alignment.json',{'relation':relation,'rois':summary,'pancreas_markers':pan_markers,'pln_markers':pln_markers})
        rec.update(status='validated_complete',finished_utc=now(),summary=summary,
                   validation={'delivery_checksums_verified':True,'delivered_objects_are_prepared_objects':True,
                               'frames_disjoint':True,'node_cells_unaligned_in_cohort':node_unaligned,
                               'tissue_routing_consistent':True,'missing_not_zero':True},
                   outputs=[identity(p) for p in out.rglob('*') if p.is_file() and p.name!='receipt.json'])
        write_json(out/'receipt.json',rec);write_json(base/'current.json',rec)
    except Exception as exc:
        rec.update(status='failed',error=str(exc),traceback=traceback.format_exc(),finished_utc=now())
        write_json(out/'receipt.json',rec);write_json(base/'last_failure.json',rec);raise


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__,formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('--manifest',required=True,type=Path);p.add_argument('--donor',required=True)
    p.add_argument('--out',type=Path,default=DEFAULT_OUT)
    a=p.parse_args();run(a.manifest,a.donor,a.out)
