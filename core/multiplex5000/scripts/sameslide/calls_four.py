#!/usr/bin/env python3
"""Carry the workstation's protein class onto one joined section's Xenium cells.

Maintainer decisions 2026-10-19. The four donors are typed on the workstation by each tissue's
current producer and the calls are delivered keyed by object_id (manifest ``protein_calls``), one
table per tissue: ``pancreas`` from the subpanel chain under the 2026-10-18 rulings, ``nodes`` so far
from ``build_pln_immune_pheno.py``, which types the lymph node's immune cells with the lineage marker
sets (``typing/lineages.py``) and the shared identity panels (``core/panels``). Every cell is to be
typed, a cell type the node shares with the pancreas by the pancreas methods (maintainer,
2026-10-19); this stage carries whatever the nodes table delivers. There is one class.
The typing puts cells in large groups by identifying features that are mutually exclusive by a
priori knowledge, then uses the full marker set in targeted ways inside each group; the delivered
designations are steps to that one result. Each table's ``class_column`` names the column that
holds it -- ``refined_class`` for the pancreas, the lineage ``immune_class`` for the nodes -- and
no other class column is carried. A group's targeted level (``within``: ``immune_class`` and
``immune_subtype`` inside the pancreas's Immune, ``t_nk_class`` inside the nodes' T_cell) is
carried only where the cell's class is that group.

This lane decides nothing about identity. It checks that the calls are about THIS lane's objects,
that no object is in two tables, and that each protein cell lies in the tissue its table types (the
join's ``protein_tissue``), then transfers them over the join's own overlap edges with
``lift_calls.transfer_calls``: largest shared area wins, and ``call_share`` and ``call_purity``
record what the transfer is worth. Every protein cell the join used competes, called or not, so a
Xenium cell lying mostly over a node is not handed a pancreas class by a sliver. A label in a
table's ``not_a_class`` is carried but never awarded, so it can never become a context class: the
pancreas's ``Excluded`` (ruled out by the producer, with its reason), and the node producer's empty
label (its cluster read structural or named no lineage; the producer names immune lineages only)
and ``Uncalled`` (immune, but no subpanel met its own requirement). A protein cell without a
delivered class is unavailable, never a class.
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

HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(HERE.parents[1]))
from import_alignment import checksums,verified,read
from join_four import current
from lift_calls import transfer_calls
from manifest import load_manifest,get_section,canonical_donor
from workflow import DEFAULT_OUT,section_root,identity,write_json,now,valid_receipt


def delivered_calls(root,donor,table=None):
    """Every delivered row for one donor, from each parquet under root (or root/table); checksums when listed."""
    root=Path(root);sums=checksums(root) if (root/'SHA256SUMS').exists() else {}
    files=sorted((root/table if table else root).rglob('*.parquet'))
    if not files:raise ValueError(f'{root}: no delivered call tables')
    records,frames=[],[]
    for f in files:
        rel=str(f.relative_to(root))
        records.append(verified(root,rel,sums) if sums else identity(f))
        t=read(f)
        if 'donor_id' not in t or 'object_id' not in t:raise ValueError(f'{f}: calls need donor_id and object_id')
        spelled=sorted(set(t.donor_id.astype(str)))
        keep=[d for d in spelled if canonical_donor(d)==donor]
        if keep:frames.append(t[t.donor_id.astype(str).isin(keep)])
    if not frames:raise ValueError(f'{donor}: no delivered calls')
    calls=pd.concat(frames,ignore_index=True)
    if calls.object_id.isna().any() or calls.object_id.astype(str).duplicated().any():raise ValueError('null or duplicate call UUIDs')
    return calls.assign(object_id=calls.object_id.astype(str)),records


def one_class(calls,column,within):
    """The class and each group's targeted levels, a level kept only inside its own group."""
    out=calls[['object_id']].assign(call=calls[column].to_numpy())
    record={}
    for group,levels in within.items():
        inside=out.call.eq(group).to_numpy()
        for level in levels:
            have=calls[level].notna().to_numpy()&calls[level].astype(str).ne('').to_numpy()
            record[level]={'group':group,'set_outside_group':int((have&~inside).sum()),'missing_inside_group':int((~have&inside).sum())}
            out[level]=calls[level].where(inside).to_numpy()
    return out,record


def unavailable(rec,out,stage,reason,why):
    """The section's current state is unavailable, with its reason, so no earlier calls stand in."""
    write_json(out/'unavailable.json',why)
    rec.update(status='scientifically_unavailable',reason=reason,finished_utc=now(),validation={'calls_awarded':False},
               outputs=[identity(out/'unavailable.json')])
    write_json(out/'receipt.json',rec);write_json(stage/'current.json',rec)
    print(f"{rec['section_id']}: calls unavailable: {reason}",flush=True)


def table_specs(spec):
    """Each delivered table's spec. A targeted level belongs to one table, and no other table carries it."""
    tables=spec['tables'];owner={}
    for name,t in tables.items():
        for level in (l for v in t.get('within',{}).values() for l in v):
            if level in owner:raise ValueError(f'{level} is a targeted level of both {owner[level]} and {name}')
            owner[level]=name
    for name,t in tables.items():
        clash=[c for c in t.get('carry',[]) if owner.get(c,name)!=name]
        if clash:raise ValueError(f'{name} carries {clash}, a targeted level of another table')
    return tables


def run(manifest_path,sid,root):
    m=load_manifest(manifest_path);section=get_section(m,sid);donor=section['donor_id']
    spec=m.get('protein_calls')
    if not spec:raise ValueError('no protein_calls declared in the manifest: the workstation typing has not been delivered')
    tables=table_specs(spec)
    base=section_root(root,section)
    paths={'joined':base/'joined/current.json','preparation':Path(root)/'protein'/f'donor_id={donor}'/'preparation/current.json'}
    joined_rec=current(paths['joined'],'joined')
    if not valid_receipt(paths['preparation']):raise ValueError(f"{paths['preparation']}: not validated")
    inputs={**{k:identity(p) for k,p in paths.items()},'code':[identity(__file__),identity(HERE/'lift_calls.py')],
            'calls':{'root':spec['root'],'tables':tables}}
    delivered_root=Path(spec['root'])
    if (delivered_root/'SHA256SUMS').exists():inputs['calls']['checksums']=identity(delivered_root/'SHA256SUMS')
    fp=hashlib.sha256(json.dumps(inputs,sort_keys=True).encode()).hexdigest()
    stage=base/'calls';stage.mkdir(parents=True,exist_ok=True)
    if valid_receipt(stage/'current.json',fp):print(f'{sid}: calls reused',flush=True);return
    out=stage/'attempts'/(fp[:16]+'-'+uuid4().hex[:8]);out.mkdir(parents=True)
    rec={'status':'running','section_id':sid,'fingerprint':fp,'inputs':inputs,'started_utc':now(),
         'attempt':str(out),'job_id':os.environ.get('SLURM_JOB_ID')}
    write_json(out/'receipt.json',rec)
    try:
        if joined_rec['status']=='scientifically_unavailable':
            unavailable(rec,out,stage,f"join: {joined_rec.get('reason','')}",{'joined':joined_rec.get('reason')});return
        frames,records={},{}
        for name,t in tables.items():
            calls,records[name]=delivered_calls(delivered_root,donor,name)
            absent=[c for c in (t['class_column'],*t.get('carry',[]),*(l for v in t.get('within',{}).values() for l in v)) if c not in calls]
            if absent:raise ValueError(f'the delivered {name} calls carry no {absent} column')
            frames[name]=calls
        rec['delivered_files']=records
        owner=pd.concat([c[['object_id']].assign(call_table=n) for n,c in frames.items()],ignore_index=True)
        if owner.object_id.duplicated().any():raise ValueError(f'{int(owner.object_id.duplicated().sum()):,} objects are in more than one delivered table')
        prep=Path(json.loads(paths['preparation'].read_text())['attempt'])
        prepared=set(read(prep/'geometry.parquet',['object_id']).object_id.astype(str))
        foreign=int((~owner.object_id.isin(prepared)).sum())
        if foreign:raise ValueError(f'{foreign:,} delivered call UUIDs are not this lane\'s objects')
        joined=Path(joined_rec['attempt'])
        edges=read(joined/'overlap_edges.parquet')
        cells=read(joined/'joined_cells.parquet',['cell_id','registered_area_um2'])
        xids=cells.cell_id.astype(str).to_numpy()
        objects=edges.drop_duplicates('object_id')[['object_id','protein_tissue']].reset_index(drop=True)
        levels,classes={},[]
        for name,t in tables.items():
            c,levels[name]=one_class(frames[name],t['class_column'],t.get('within',{}))
            carried=frames[name][['object_id',*t.get('carry',[])]]
            classes.append(c.merge(carried,on='object_id',how='left',validate='one_to_one').assign(call_table=name))
        classes=pd.concat(classes,ignore_index=True)
        # Every protein cell the join used competes; one no delivered table covers has no class.
        lifted=objects.merge(classes,on='object_id',how='left',validate='one_to_one')
        # A table types one tissue; a protein cell the join placed in the other tissue is a routing error.
        typed=lifted.call_table.notna()
        misplaced=typed&lifted.protein_tissue.ne(lifted.call_table.map({n:t['tissue'] for n,t in tables.items()}))
        if misplaced.any():
            raise ValueError(f'{int(misplaced.sum()):,} protein cells lie in one tissue by the join and in the other tissue\'s table: '
                             +str(lifted[misplaced].groupby(['protein_tissue','call_table']).size().to_dict()))
        refused={n:set(t.get('not_a_class',[])) for n,t in tables.items()}
        is_refused=np.array([c in refused[t] if isinstance(t,str) else False for c,t in zip(lifted.call,lifted.call_table)],dtype=bool)
        coverage={'protein_cells_in_edges':len(objects),'in_delivered_tables':int(typed.sum()),
                  'with_class':int((lifted.call.notna()&~is_refused).sum()),'not_a_class':int(is_refused.sum()),
                  'by_join_tissue_and_table':{f'{a}|{b}':int(n) for (a,b),n in
                                              lifted.fillna({'call_table':'none'}).groupby(['protein_tissue','call_table']).size().items()}}
        out_cells=transfer_calls(xids,edges[['cell_id','object_id','overlap_area_um2']],lifted[['object_id','call']],cells.registered_area_um2.to_numpy(float))
        by_object=lifted.set_index('object_id')
        out_cells['winning_tissue']=out_cells.winning_object_id.map(by_object.protein_tissue).to_numpy()
        out_cells['call_table']=out_cells.winning_object_id.map(by_object.call_table).to_numpy()
        # Carried, never awarded: a producer's "ruled out" or "no label" stays visible and cannot become a class.
        for name,labels in refused.items():
            out_cells.loc[out_cells.call_table.eq(name)&out_cells.call.isin(labels),'status']='not_a_class'
        out_cells.loc[out_cells.winning_object_id.notna()&out_cells.call_table.isna(),'status']='no_class_delivered'
        extra=[c for c in classes.columns if c not in ('object_id','call','call_table')]
        carried=lifted[['object_id',*extra]].rename(columns={c:'protein__'+c for c in extra})
        out_cells=out_cells.merge(carried,left_on='winning_object_id',right_on='object_id',how='left',validate='many_to_one').drop(columns=['object_id'])
        if not np.array_equal(out_cells.cell_id.astype(str).to_numpy(),xids):raise ValueError('transfer changed the Xenium row order')
        awarded=out_cells.status.eq('awarded')
        refused_cells=out_cells[out_cells.status.eq('not_a_class')]
        refused_cells=refused_cells.assign(call=refused_cells.call.replace('','""'))
        why=['call_table','call',*(['protein__is_immune'] if 'protein__is_immune' in refused_cells else [])]
        summary={'n_xenium_cells':len(out_cells),'awarded':int(awarded.sum()),
                 'class_column':{n:t['class_column'] for n,t in tables.items()},
                 'status':out_cells.status.value_counts().to_dict(),
                 'status_by_winning_tissue':{f'{t or "none"}|{s}':int(n) for (t,s),n in
                                             out_cells.fillna({'winning_tissue':''}).groupby(['winning_tissue','status']).size().items()},
                 'calls':out_cells.loc[awarded,'call'].value_counts().to_dict(),
                 'calls_by_table':{f'{t}|{c}':int(n) for (t,c),n in out_cells[awarded].groupby(['call_table','call']).size().items()},
                 'not_a_class_by_label':{'|'.join(map(str,k)):int(n) for k,n in
                                         refused_cells.astype({c:str for c in why}).groupby(why).size().items()},
                 'call_share_median':float(np.nanmedian(out_cells.loc[awarded,'call_share'])) if awarded.any() else None,
                 'call_purity_median':float(np.nanmedian(out_cells.loc[awarded,'call_purity'])) if awarded.any() else None,
                 'call_purity_ge_0_8':float((out_cells.loc[awarded,'call_purity']>=.8).mean()) if awarded.any() else None,
                 'coverage':coverage,'targeted_levels':levels}
        write_json(out/'calls.json',summary)
        if not awarded.any():
            unavailable(rec,out,stage,'no Xenium cell\'s largest protein partner carries a delivered class '
                        f"(the delivered tables cover {coverage['in_delivered_tables']:,} of the {len(objects):,} protein cells this join used)",summary);return
        out_cells.to_parquet(out/'protein_calls.parquet',index=False)
        rec.update(status='validated_complete',finished_utc=now(),summary=summary,
                   validation={'delivered_objects_are_prepared_objects':True,'tables_disjoint':True,'join_tissue_matches_table':True,
                               'same_edges_as_join':True,'rows_equal_join':True,'not_a_class_never_awarded':True,'one_class_column':True},
                   outputs=[identity(p) for p in out.iterdir() if p.is_file() and p.name!='receipt.json'])
        write_json(out/'receipt.json',rec);write_json(stage/'current.json',rec)
        print(f"{sid}: {summary['awarded']:,} of {len(out_cells):,} Xenium cells awarded a class; "
              f"purity median {summary['call_purity_median']}",flush=True)
    except Exception as exc:
        rec.update(status='failed',error=str(exc),traceback=traceback.format_exc(),finished_utc=now())
        write_json(out/'receipt.json',rec);write_json(stage/'last_failure.json',rec);raise


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__,formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('--manifest',required=True,type=Path);p.add_argument('--section',required=True)
    p.add_argument('--out',type=Path,default=DEFAULT_OUT)
    a=p.parse_args();run(a.manifest,a.section,a.out)
