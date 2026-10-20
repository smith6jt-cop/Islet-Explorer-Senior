#!/usr/bin/env python3
"""Islet summaries for one same-slide pancreas section.

The islet objects are the workstation's own (manifest ``pheno.islets``), outlines eroded to whole
cells on 2026-10-02: ``Islet`` (two or more whole cells, at least 35 µm across),
``EndocrineCluster`` (two or more, under 35 µm) and ``EndocrineSingle`` (one). A ``Hormone
Fragment`` holds no whole cell and is never an islet, so it is not read. Every Xenium cell is
placed by its own centroid, mapped through the validated transform into the PhenoCycler frame the
objects were drawn in, and its signed distance to the nearest object is measured there (negative
inside). Where a Xenium cell sits does not depend on which PhenoCycler cell it overlaps.

The approved pipeline has no peri-islet expansion: every cell carries its distance to the nearest
islet, so a cell is in an object's ``core`` or in ``tissue`` and the distance stays a continuous
measurement, with no band drawn in its place (the delivered pancreas calls: ``cell_region`` and
``signed_dist_islet_um``, measured to these same three object classes). This stage does the same.

Composition is the transferred class (``calls_four.py``: one class, the typing's result) of the
Xenium cells in pancreas tissue. Per object it is set beside the object's own PhenoCycler cells --
the workstation's class of the cells it names (``endocrine_object_id``) -- which says what the
transfer kept. ``coverage`` is the share of an object's outline that Xenium imaged, so an object
cut by the edge of the Xenium region is not read as a small one. RNA is summarised per region and
class as the gDNA-corrected detection of each identity group of the panel taxonomy.

Descriptive only: one section per donor, run and batch confounded, and no test between donors or
disease groups (maintainer, 2026-10-19: no hypotheses until the data are properly processed).
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
from shapely.affinity import scale
from shapely.geometry import shape
from scipy.spatial import cKDTree
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(HERE.parents[1]))
from phenocycler.integration.transform import Transform
import analysis_input
from calls_four import delivered_calls
from compare_sections import identity_groups
from import_alignment import read
from join_four import current
from manifest import load_manifest,get_section
from workflow import DEFAULT_OUT,section_root,identity,write_json,now,valid_receipt

WHOLE_CELL=('Islet','EndocrineCluster','EndocrineSingle')
REGIONS=('core','tissue')
#: Coverage is read on a grid this fine; a grid point counts as imaged when a Xenium centroid lies this close.
COVER_STEP_UM,COVER_RADIUS_UM=5.0,10.0
PROFILE_EDGES=np.arange(-30.,101.,5.)
ENDOCRINE=('Alpha','Beta','Delta','Endocrine (ambiguous subtype)')
INK,MUTED,GRID,SURFACE=('#0b0b0b','#52514e','#e4e3df','#fcfcfb')
SERIES=('#2a78d6','#eb6834','#1baf7a')


def islet_objects(path,px):
    """Whole-cell islet objects in µm; a Hormone Fragment is never an islet."""
    rows=[]
    for f in json.loads(Path(path).read_text())['features']:
        cls=(f['properties'].get('classification') or {}).get('name','')
        if cls in WHOLE_CELL:
            g=shapely.make_valid(scale(shape(f['geometry']),xfact=px,yfact=px,origin=(0,0)))
            rows.append({'object_id':str(f['id']),'object_class':cls,'name':f['properties'].get('name'),'geometry':g})
    t=pd.DataFrame(rows)
    if t.empty or t.object_id.duplicated().any():raise ValueError(f'{path}: no islet objects, or duplicate UUIDs')
    return t


def signed_distance(xy,polys):
    """Signed distance to the nearest object (negative inside, the deepest object wins) and its row."""
    pts=shapely.points(xy);tree=shapely.STRtree(polys)
    near,dist=tree.query_nearest(pts,return_distance=True,all_matches=False)
    d=np.full(len(xy),np.inf);idx=np.full(len(xy),-1)
    d[near[0]]=dist;idx[near[0]]=near[1]
    pi,oi=tree.query(pts,predicate='within')
    if len(pi):
        depth=-shapely.distance(pts[pi],shapely.boundary(polys[oi]))
        order=np.lexsort((depth,pi))
        first=np.r_[True,pi[order][1:]!=pi[order][:-1]]
        d[pi[order][first]]=depth[order][first];idx[pi[order][first]]=oi[order][first]
    return d,idx


def region_of(d):
    """Inside or on an object's outline is its core; every other cell is tissue. No band."""
    return np.where(d<=0,'core','tissue')


def coverage(polys,xy,step=COVER_STEP_UM,radius=COVER_RADIUS_UM):
    """Share of each object's outline with a Xenium centroid within radius of a grid point."""
    tree=cKDTree(xy);out=np.zeros(len(polys))
    for k,p in enumerate(polys):
        x0,y0,x1,y1=p.bounds
        gx,gy=np.meshgrid(np.arange(x0,x1+step,step),np.arange(y0,y1+step,step))
        g=np.c_[gx.ravel(),gy.ravel()];g=g[shapely.contains_xy(p,g[:,0],g[:,1])]
        if not len(g):g=np.array([[p.centroid.x,p.centroid.y]])
        out[k]=float(np.mean(np.isfinite(tree.query(g,k=1,distance_upper_bound=radius)[0])))
    return out


def counts_by(frame,keys,column='call'):
    return frame.groupby([*keys,column],dropna=False).size().rename('n').reset_index()


def in_object(frame):
    """A core cell is counted under its object's class; a tissue cell is in no object, whatever is nearest."""
    return frame.assign(object_class=frame.object_class.where(frame.region.eq('core')))


def figure(comp,profile,objects,title,path):
    parts={'Islet core':comp[comp.region.eq('core')&comp.object_class.eq('Islet')],'tissue':comp[comp.region.eq('tissue')]}
    classes=list(pd.concat(parts.values()).groupby('call').n.sum().sort_values(ascending=False).index)
    fig=plt.figure(figsize=(15,5.4),facecolor=SURFACE,constrained_layout=True)
    grid=fig.add_gridspec(1,3)
    ax=fig.add_subplot(grid[0,0]);ax.set_facecolor(SURFACE)
    y=np.arange(len(classes))[::-1]
    for k,((label,part),colour) in enumerate(zip(parts.items(),SERIES)):
        r=part.set_index('call').reindex(classes)
        ax.scatter(100*r.fraction.fillna(0),y+(.5-k)*.3,s=46,color=colour,edgecolor=SURFACE,linewidth=1.5,zorder=3,label=label)
    ax.set_xscale('symlog',linthresh=.1);ax.set_xlim(0,100);ax.set_yticks(y);ax.set_yticklabels(classes,fontsize=8.5,color=INK)
    ax.xaxis.set_major_formatter(matplotlib.ticker.FuncFormatter(lambda v,_:f'{v:g}%'))
    ax.set_xlabel('share of the region\'s classed Xenium cells',color=MUTED,fontsize=9)
    ax.set_title('Composition of the Islet core and of tissue',fontsize=10,color=INK,loc='left')
    ax.legend(frameon=False,fontsize=8,loc='lower right');ax.grid(color=GRID,lw=.8)
    ax=fig.add_subplot(grid[0,1]);ax.set_facecolor(SURFACE)
    for (label,members),colour in zip((('endocrine',ENDOCRINE),('Immune',('Immune',)),('Blood vessel',('Blood vessel',))),SERIES):
        p=profile[profile.call.isin(members)].groupby('bin_mid').n.sum()
        tot=profile.groupby('bin_mid').n.sum()
        share=(p.reindex(tot.index).fillna(0)/tot).sort_index()
        ax.plot(share.index,100*share.to_numpy(),color=colour,lw=2,label=label,solid_capstyle='round')
    ax.axvline(0,color=MUTED,lw=1,ls=(0,(3,3)))
    ax.set_xlabel('signed distance to the nearest Islet object (µm; negative inside)',color=MUTED,fontsize=9)
    ax.set_ylabel('share of classed Xenium cells',color=MUTED,fontsize=9)
    ax.yaxis.set_major_formatter(matplotlib.ticker.PercentFormatter(100,decimals=0))
    ax.set_title('Class share by distance from the islet edge',fontsize=10,color=INK,loc='left')
    ax.legend(frameon=False,fontsize=8);ax.grid(color=GRID,lw=.8)
    ax=fig.add_subplot(grid[0,2]);ax.set_facecolor(SURFACE)
    q=objects[objects.object_class.eq('Islet')&(objects.coverage>=.95)&(objects.core_awarded>=5)&(objects.native_classed>=5)]
    ax.plot([0,1],[0,1],color=MUTED,lw=1,zorder=1)
    ax.scatter(q.native_endocrine_fraction,q.core_endocrine_fraction,s=30,color=SERIES[0],edgecolor=SURFACE,linewidth=1.2,zorder=3)
    ax.set_xlim(-.02,1.02);ax.set_ylim(-.02,1.02)
    for a in (ax.xaxis,ax.yaxis):a.set_major_formatter(matplotlib.ticker.PercentFormatter(1.0,decimals=0))
    ax.set_xlabel('endocrine share, the object\'s PhenoCycler cells',color=MUTED,fontsize=9)
    ax.set_ylabel('endocrine share, Xenium cells in its core',color=MUTED,fontsize=9)
    ax.set_title(f'What the transfer kept, per Islet object (n = {len(q)})',fontsize=10,color=INK,loc='left')
    ax.grid(color=GRID,lw=.8)
    for a in fig.axes:
        a.tick_params(labelsize=8,colors=MUTED)
        for s in a.spines.values():s.set_color(GRID)
    fig.suptitle(f'{title} · descriptive; one section per donor',x=.01,ha='left',fontsize=10,color=INK)
    fig.savefig(path.with_suffix('.png'),dpi=160,facecolor=SURFACE);fig.savefig(path.with_suffix('.pdf'),facecolor=SURFACE);plt.close(fig)


def run(manifest_path,sid,root):
    m=load_manifest(manifest_path);section=get_section(m,sid);base=section_root(root,section)
    if section['tissue']!='pancreas' or section['design']!='same_slide':raise ValueError(f'{sid}: islet summaries are for same-slide pancreas sections')
    ph=section['pheno'];spec=m['protein_calls'];pspec=spec['tables']['pancreas']
    paths={k:base/k/'current.json' for k in ('registration','joined','calls','rna')}
    rec_in={k:current(p,k) for k,p in paths.items()}
    inputs={**{k:identity(p) for k,p in paths.items()},'islets':identity(ph['islets']),'islet_table':identity(ph['islet_table']),
            'native_calls':identity(Path(spec['root'])/'SHA256SUMS'),
            'code':[identity(__file__),identity(HERE/'compare_sections.py')],
            'parameters':{'regions':REGIONS,'whole_cell':WHOLE_CELL,'cover_step_um':COVER_STEP_UM,'cover_radius_um':COVER_RADIUS_UM}}
    fp=hashlib.sha256(json.dumps(inputs,sort_keys=True).encode()).hexdigest()
    stage=base/'islets';stage.mkdir(parents=True,exist_ok=True)
    if valid_receipt(stage/'current.json',fp):print(f'{sid}: islets reused',flush=True);return
    out=stage/'attempts'/(fp[:16]+'-'+uuid4().hex[:8]);out.mkdir(parents=True)
    rec={'status':'running','section_id':sid,'fingerprint':fp,'inputs':inputs,'started_utc':now(),
         'attempt':str(out),'job_id':os.environ.get('SLURM_JOB_ID')}
    write_json(out/'receipt.json',rec)
    try:
        gone={k:r.get('reason') for k,r in rec_in.items() if r['status']=='scientifically_unavailable'}
        if gone:
            write_json(out/'unavailable.json',gone)
            rec.update(status='scientifically_unavailable',reason='; '.join(f'{k}: {v}' for k,v in gone.items()),finished_utc=now(),
                       validation={'upstream_unavailable':sorted(gone)},outputs=[identity(out/'unavailable.json')])
            write_json(out/'receipt.json',rec);write_json(stage/'current.json',rec);print(f'{sid}: islets unavailable',flush=True);return
        for k,p in paths.items():
            if not valid_receipt(p):raise ValueError(f'{p}: not validated')
        reg=rec_in['registration'];tf=Transform.load(Path(reg['transform_directory']),stem=reg['transform_stem'])
        joined=read(Path(rec_in['joined']['attempt'])/'joined_cells.parquet',['cell_id','x_um','y_um','protein_tissue'])
        calls=read(Path(rec_in['calls']['attempt'])/'protein_calls.parquet',['cell_id','call','status','winning_object_id'])
        if not np.array_equal(joined.cell_id.astype(str).to_numpy(),calls.cell_id.astype(str).to_numpy()):raise ValueError('calls and join rows differ')
        objs=islet_objects(ph['islets'],float(ph['pixel_size_um']))
        polys=objs.geometry.to_numpy()
        xy=tf.apply(joined[['x_um','y_um']].to_numpy(float))
        d,idx=signed_distance(xy,polys)
        # Every cell keeps its nearest object and its distance to it, the approved pipeline's spelling.
        cells=pd.DataFrame({'cell_id':joined.cell_id.astype(str),'x_pheno_um':xy[:,0],'y_pheno_um':xy[:,1],'signed_dist_islet_um':d,
                            'object_id':objs.object_id.to_numpy()[idx],'object_class':objs.object_class.to_numpy()[idx],
                            'region':region_of(d),'tissue':joined.protein_tissue.to_numpy(),'status':calls.status.to_numpy(),
                            'call':calls.call.where(calls.status.eq('awarded')).to_numpy()})
        cells.to_parquet(out/'cells.parquet',index=False)
        panc=cells[cells.tissue.eq('pancreas')]
        classed=panc[panc.call.notna()]
        comp=counts_by(in_object(classed),['region','object_class'])
        comp['fraction']=comp.n/comp.groupby(['region','object_class'],dropna=False).n.transform('sum')
        comp.to_parquet(out/'composition.parquet',index=False)
        status=counts_by(in_object(panc),['region','object_class'],'status')
        status.to_parquet(out/'status.parquet',index=False)
        near=classed[classed.object_class.eq('Islet')]
        mids=(PROFILE_EDGES[:-1]+PROFILE_EDGES[1:])/2
        profile=near.assign(bin_mid=pd.cut(near.signed_dist_islet_um,PROFILE_EDGES,labels=mids).astype(float)).dropna(subset=['bin_mid'])
        profile=counts_by(profile,['bin_mid']);profile.to_parquet(out/'profile.parquet',index=False)
        # Per object: its Xenium core, its coverage, and its own PhenoCycler cells.
        table=pd.read_csv(ph['islet_table'],sep='\t').rename(columns={'Object ID':'object_id','Whole cells':'whole_cells',
                          'Mean caliper diameter µm':'caliper_um','Num Detections':'detections'})
        objects=objs.drop(columns=['geometry']).assign(area_um2=shapely.area(polys),coverage=coverage(polys,xy))
        objects=objects.merge(table[['object_id','whole_cells','caliper_um','detections']],on='object_id',how='left',validate='one_to_one')
        core=panc[panc.region.eq('core')]
        objects['core_cells']=objects.object_id.map(core.object_id.value_counts()).fillna(0).astype(int)
        objects['core_awarded']=objects.object_id.map(core[core.call.notna()].object_id.value_counts()).fillna(0).astype(int)
        wide=core[core.call.notna()].groupby(['object_id','call']).size().unstack(fill_value=0).add_prefix('core__')
        # A count is a count: an object with no cell of a class has zero of it.
        objects=objects.merge(wide,left_on='object_id',right_index=True,how='left').fillna({c:0 for c in wide.columns})
        core_end=sum((objects[f'core__{c}'] for c in ENDOCRINE if f'core__{c}' in objects),pd.Series(0,index=objects.index))
        objects['core_endocrine_fraction']=np.divide(core_end.to_numpy(float),objects.core_awarded.to_numpy(float),
                                                     out=np.full(len(objects),np.nan),where=objects.core_awarded.to_numpy()>0)
        native,_=delivered_calls(spec['root'],section['donor_id'],'pancreas')
        native=native[native.endocrine_object_id.notna()&native[pspec['class_column']].notna()
                      &~native[pspec['class_column']].isin(pspec.get('not_a_class',[]))]
        native=native[native.endocrine_object_id.astype(str).isin(set(objects.object_id))]
        nw=native.groupby([native.endocrine_object_id.astype(str),pspec['class_column']]).size().unstack(fill_value=0)
        objects['native_classed']=objects.object_id.map(nw.sum(axis=1)).fillna(0).astype(int)
        nend=nw[[c for c in ENDOCRINE if c in nw]].sum(axis=1)
        objects['native_endocrine_fraction']=objects.object_id.map(nend/nw.sum(axis=1))
        objects=objects.merge(nw.add_prefix('native__'),left_on='object_id',right_index=True,how='left').fillna({f'native__{c}':0 for c in nw.columns})
        objects.to_parquet(out/'objects.parquet',index=False)
        # RNA: share of each identity group's genes detected (gDNA-corrected), per region and class.
        rna=analysis_input.load_base(Path(rec_in['rna']['attempt'])/'rna.h5ad',donor=section['donor_id'],roi=section['roi'],serial=section['xenium_serial'])
        if not np.array_equal(rna.obs_names.astype(str).to_numpy(),cells.cell_id.to_numpy()):raise ValueError('RNA and join rows differ')
        det=rna.layers['detected_gdna'].tocsc() if hasattr(rna.layers['detected_gdna'],'tocsc') else rna.layers['detected_gdna']
        keep=(cells.tissue.eq('pancreas')&cells.call.notna()).to_numpy()
        rows=[]
        for group,genes in identity_groups('pancreas').items():
            j=[rna.var_names.get_loc(g) for g in genes if g in rna.var_names]
            if not j:continue
            share=np.asarray(det[:,j].astype(np.int32).sum(axis=1)).ravel()/len(j)
            frame=pd.DataFrame({'region':cells.region.to_numpy()[keep],'call':cells.call.to_numpy()[keep],'share':share[keep]})
            for (region,call),g in frame.groupby(['region','call']):
                rows.append({'group':group,'genes_on_panel':len(j),'region':region,'call':call,'n_cells':len(g),
                             'mean_share_detected':float(g.share.mean()),'cells_detecting_any':float((g.share>0).mean())})
        pd.DataFrame(rows).to_parquet(out/'rna.parquet',index=False)
        label='HDL112' if section['donor_id']=='112' else section['donor_id']
        figure(comp,profile,objects,f"{label} · pancreas · Xenium {section['xenium_serial']}",out/'islets')
        imaged=objects[objects.core_cells>0]
        summary={'objects':objects.object_class.value_counts().to_dict(),
                 'objects_with_xenium_core_cells':imaged.object_class.value_counts().to_dict(),
                 'objects_fully_imaged':objects[objects.coverage>=.95].object_class.value_counts().to_dict(),
                 'xenium_cells_by_region':panc.region.value_counts().to_dict(),
                 'classed_by_region':classed.region.value_counts().to_dict(),
                 'cells_not_in_pancreas_tissue':int((~cells.tissue.eq('pancreas')).sum()),
                 'islet_core_composition':comp[comp.object_class.eq('Islet')&comp.region.eq('core')].set_index('call').fraction.round(4).to_dict(),
                 'tissue_composition':comp[comp.region.eq('tissue')].set_index('call').fraction.round(4).to_dict(),
                 'transform_scale':tf.scale}
        write_json(out/'islets.json',{**summary,'parameters':inputs['parameters'],
                   'interpretation':'descriptive; one section per donor; no test between donors or disease groups'})
        rec.update(status='validated_complete',finished_utc=now(),summary=summary,
                   validation={'core_or_tissue_from_signed_distance':bool(set(cells.region)<=set(REGIONS)),
                               'hormone_fragments_never_islets':True,'rows_equal_join':True,'pancreas_tissue_only':True},
                   outputs=[identity(p) for p in out.iterdir() if p.is_file() and p.name!='receipt.json'])
        write_json(out/'receipt.json',rec);write_json(stage/'current.json',rec)
        print(f"{sid}: {summary['objects_with_xenium_core_cells']} objects with Xenium core cells; regions {summary['xenium_cells_by_region']}",flush=True)
    except Exception as exc:
        rec.update(status='failed',error=str(exc),traceback=traceback.format_exc(),finished_utc=now())
        write_json(out/'receipt.json',rec);write_json(stage/'last_failure.json',rec);raise


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__,formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('--manifest',required=True,type=Path);p.add_argument('--section',required=True)
    p.add_argument('--out',type=Path,default=DEFAULT_OUT)
    a=p.parse_args();run(a.manifest,a.section,a.out)
