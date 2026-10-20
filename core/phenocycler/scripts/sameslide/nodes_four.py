#!/usr/bin/env python3
"""Lymph-node composition and depth summaries for one same-slide section that imaged a node.

The nodes are the image's pLN annotations (manifest ``pheno.annotations``; a feature's id is the
QuPath UUID that the node table's ``pln_annotation_id`` names). A Xenium cell belongs to the node
its largest protein partner lies in (``calls_four``: ``protein__pln_annotation_id``), so a node is
read in whichever section imaged it -- a different-block node in its pLN section, the same-block
node of 6347 and 6580 in their pancreas sections. Every Xenium cell is placed by its own centroid,
mapped through the validated transform into the PhenoCycler frame the outlines were drawn in, and
its depth is its distance inside that node's outline (negative outside it).

Composition is the transferred class (``calls_four.py``: one class, the lineage, with the T/NK
subset carried inside T_cell) of each node's Xenium cells, set beside the same node's PhenoCycler
cells that this section's Xenium cells overlap -- the workstation's class of those very cells --
which says what the transfer kept. A cell without a class is counted by why, never as a class: the
producer names immune lineages only, so a cell its cluster read structural has none. ``coverage``
is the share of the node's outline that Xenium imaged. RNA is summarised per class and subset as
the gDNA-corrected detection of each shared identity group of the panel taxonomy.

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
from islets_four import coverage,COVER_STEP_UM,COVER_RADIUS_UM
from join_four import current
from manifest import load_manifest,get_section
from workflow import DEFAULT_OUT,section_root,identity,write_json,now,valid_receipt

NODE_CLASS='pLN'
#: Depth is measured to the outline densified at this spacing, so it is exact to half of it.
EDGE_STEP_UM=1.0
DEPTH_STEP_UM=50.0
SUBSET_OF='T_cell'
INK,MUTED,GRID,SURFACE=('#0b0b0b','#52514e','#e4e3df','#fcfcfb')
SERIES=('#2a78d6','#eb6834','#1baf7a','#8a5cc2','#c9a227','#d6457a')


def node_outlines(path,px):
    """The image's lymph-node objects in µm, by QuPath UUID."""
    rows=[]
    for f in json.loads(Path(path).read_text())['features']:
        if (f['properties'].get('classification') or {}).get('name','')==NODE_CLASS:
            rows.append({'node_id':str(f['id']),'geometry':shapely.make_valid(scale(shape(f['geometry']),xfact=px,yfact=px,origin=(0,0)))})
    t=pd.DataFrame(rows)
    if t.empty or t.node_id.duplicated().any():raise ValueError(f'{path}: no {NODE_CLASS} objects, or duplicate UUIDs')
    return t


def depth(xy,poly,step=EDGE_STEP_UM):
    """Distance inside the outline, negative outside it, to the outline densified at step."""
    edge=shapely.get_coordinates(shapely.segmentize(shapely.boundary(poly),step))
    d=cKDTree(edge).query(xy,k=1)[0]
    return np.where(shapely.contains_xy(poly,xy[:,0],xy[:,1]),d,-d)


def why(call,is_immune,refused):
    """Why a cell carries no class, from the producer's own label; '' where it carries one."""
    call=pd.Series(call,dtype=object).reset_index(drop=True)
    imm=pd.Series(is_immune,dtype=object).reset_index(drop=True)
    imm=imm.map(lambda v:None if v is None or v is pd.NA or (isinstance(v,float) and np.isnan(v)) else bool(v))
    out=np.full(len(call),'',dtype=object)
    out[call.isna().to_numpy()]='protein call unavailable'
    empty=call.eq('').to_numpy()
    out[empty&imm.eq(False).to_numpy()]='not immune: the producer names immune lineages only'
    out[empty&imm.isna().to_numpy()]='no lineage marker named'
    out[empty&imm.eq(True).to_numpy()]='immune, no label'
    other=call.isin(refused).to_numpy()&~empty
    out[other]=call[other].map(lambda v:f'{v}: immune, no subpanel named a lineage' if v=='Uncalled' else f'{v}: not a class').to_numpy()
    return out


def shares(frame,keys,value,denominator):
    """Counts of value within keys, and each count's share of the denominator rows of the same keys."""
    n=frame.groupby([*keys,value]).size().rename('n').reset_index()
    total=denominator.groupby(keys).size().rename('total').reset_index()
    n=n.merge(total,on=keys,how='left')
    return n.assign(share=n.n/n.total).drop(columns='total')


def figure(comp,profile,rna,objects,title,path):
    fig=plt.figure(figsize=(15,9.6),facecolor=SURFACE,constrained_layout=True)
    grid=fig.add_gridspec(2,2)
    classes=list(comp[comp.level.eq('class')].groupby('value').n.sum().sort_values(ascending=False).index)
    colour={c:SERIES[k%len(SERIES)] for k,c in enumerate(classes)}
    for k,(level,label) in enumerate((('class','class, among classed cells'),('subset',f'{SUBSET_OF} subset, among {SUBSET_OF} cells'))):
        ax=fig.add_subplot(grid[0,k]);ax.set_facecolor(SURFACE)
        part=comp[comp.level.eq(level)];values=list(part.groupby('value').n.sum().sort_values(ascending=False).index)
        pal={v:colour.get(v,SERIES[j%len(SERIES)]) for j,v in enumerate(values)}
        bars=[(n,s) for n in objects.node_id for s in ('xenium','native')]
        for y,(node,source) in enumerate(bars[::-1]):
            r=part[part.node_id.eq(node)&part.source.eq(source)].set_index('value').share.reindex(values).fillna(0)
            left=0.
            for v,w in r.items():
                ax.barh(y,w,left=left,color=pal[v],edgecolor=SURFACE,linewidth=.8,height=.72)
                left+=w
        ax.set_yticks(range(len(bars)))
        ax.set_yticklabels([f"{n[:8]} · {'Xenium cells' if s=='xenium' else 'its PhenoCycler cells'}" for n,s in bars[::-1]],fontsize=8.5,color=INK)
        ax.set_xlim(0,1);ax.xaxis.set_major_formatter(matplotlib.ticker.PercentFormatter(1.0,decimals=0))
        ax.set_title(f'Composition by {label}',fontsize=10,color=INK,loc='left')
        ax.legend(handles=[matplotlib.patches.Patch(color=pal[v],label=v) for v in values],frameon=False,fontsize=8,
                  loc='upper left',bbox_to_anchor=(0,-.07),ncol=min(len(values),5))
        ax.grid(axis='x',color=GRID,lw=.8)
    ax=fig.add_subplot(grid[1,0]);ax.set_facecolor(SURFACE)
    big=objects.sort_values('xenium_cells',ascending=False).node_id.iloc[0]
    p=profile[profile.node_id.eq(big)];tot=p.groupby('bin_mid').n.sum()
    # The bin's classed cells, behind the shares, so a share read from a few cells is seen as one.
    tw=ax.twinx();tw.bar(tot.index,tot.to_numpy(),width=.9*DEPTH_STEP_UM,color=GRID,zorder=0)
    tw.set_ylabel('classed Xenium cells in the bin',color=MUTED,fontsize=9);ax.set_zorder(tw.get_zorder()+1);ax.patch.set_visible(False)
    for c in classes:
        s=(p[p.call.eq(c)].groupby('bin_mid').n.sum().reindex(tot.index).fillna(0)/tot).sort_index()
        ax.plot(s.index,100*s.to_numpy(),color=colour[c],lw=2,label=c,solid_capstyle='round')
    ax.axvline(0,color=MUTED,lw=1,ls=(0,(3,3)))
    ax.set_xlabel(f'depth inside the node outline (µm; {big[:8]}; negative outside)',color=MUTED,fontsize=9)
    ax.set_ylabel('share of classed Xenium cells',color=MUTED,fontsize=9)
    ax.yaxis.set_major_formatter(matplotlib.ticker.PercentFormatter(100,decimals=0))
    ax.set_title('Class share by depth from the node edge',fontsize=10,color=INK,loc='left')
    ax.legend(frameon=False,fontsize=8);ax.grid(color=GRID,lw=.8)
    ax=fig.add_subplot(grid[1,1]);ax.set_facecolor(SURFACE)
    h=rna[rna.node_id.eq(big)&rna.subset.eq('')].pivot_table(index='group',columns='call',values='mean_share_detected')
    h=h.reindex(columns=[c for c in classes if c in h.columns])
    im=ax.imshow(h.to_numpy(float),cmap='Blues',aspect='auto',vmin=0)
    ax.set_xticks(range(h.shape[1]));ax.set_xticklabels(h.columns,fontsize=8.5,color=INK,rotation=20,ha='right')
    ax.set_yticks(range(h.shape[0]));ax.set_yticklabels(h.index,fontsize=8.5,color=INK)
    cb=fig.colorbar(im,ax=ax,shrink=.8);cb.ax.tick_params(labelsize=8,colors=MUTED)
    cb.set_label('mean share of the group\'s genes detected (gDNA-corrected)',color=MUTED,fontsize=8.5)
    ax.set_title(f'RNA identity groups by transferred class ({big[:8]})',fontsize=10,color=INK,loc='left')
    for a in fig.axes:
        a.tick_params(labelsize=8,colors=MUTED)
        for s in a.spines.values():s.set_color(GRID)
    fig.suptitle(f'{title} · descriptive; one section per donor',x=.01,ha='left',fontsize=10,color=INK)
    fig.savefig(path.with_suffix('.png'),dpi=160,facecolor=SURFACE);fig.savefig(path.with_suffix('.pdf'),facecolor=SURFACE);plt.close(fig)


def run(manifest_path,sid,root):
    m=load_manifest(manifest_path);section=get_section(m,sid);base=section_root(root,section)
    if section['design']!='same_slide':raise ValueError(f'{sid}: node summaries are for same-slide sections')
    ph=section['pheno'];spec=m['protein_calls'];nspec=spec['tables']['nodes']
    paths={k:base/k/'current.json' for k in ('registration','joined','calls','rna')}
    rec_in={k:current(p,k) for k,p in paths.items()}
    inputs={**{k:identity(p) for k,p in paths.items()},'annotations':identity(ph['annotations']),
            'native_calls':identity(Path(spec['root'])/'SHA256SUMS'),
            'code':[identity(__file__),identity(HERE/'islets_four.py'),identity(HERE/'compare_sections.py')],
            'parameters':{'node_class':NODE_CLASS,'edge_step_um':EDGE_STEP_UM,'depth_step_um':DEPTH_STEP_UM,'subset_of':SUBSET_OF,
                          'cover_step_um':COVER_STEP_UM,'cover_radius_um':COVER_RADIUS_UM}}
    fp=hashlib.sha256(json.dumps(inputs,sort_keys=True).encode()).hexdigest()
    stage=base/'node_summaries';stage.mkdir(parents=True,exist_ok=True)
    if valid_receipt(stage/'current.json',fp):print(f'{sid}: node summaries reused',flush=True);return
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
            write_json(out/'receipt.json',rec);write_json(stage/'current.json',rec);print(f'{sid}: node summaries unavailable',flush=True);return
        for k,p in paths.items():
            if not valid_receipt(p):raise ValueError(f'{p}: not validated')
        reg=rec_in['registration'];tf=Transform.load(Path(reg['transform_directory']),stem=reg['transform_stem'])
        joined=read(Path(rec_in['joined']['attempt'])/'joined_cells.parquet',['cell_id','x_um','y_um'])
        calls=read(Path(rec_in['calls']['attempt'])/'protein_calls.parquet',
                   ['cell_id','call','status','call_table','protein__pln_annotation_id','protein__block_relationship',
                    'protein__t_nk_class','protein__is_immune'])
        if not np.array_equal(joined.cell_id.astype(str).to_numpy(),calls.cell_id.astype(str).to_numpy()):raise ValueError('calls and join rows differ')
        node=calls.call_table.eq('nodes').to_numpy()
        if not node.any():raise ValueError(f"{sid}: no Xenium cell's largest protein partner lies in a lymph node; node summaries do not apply")
        outlines=node_outlines(ph['annotations'],float(ph['pixel_size_um'])).set_index('node_id')
        xy=tf.apply(joined[['x_um','y_um']].to_numpy(float))
        awarded=calls.status.eq('awarded').to_numpy()
        refused=set(nspec.get('not_a_class',[]))
        cells=pd.DataFrame({'cell_id':joined.cell_id.astype(str).to_numpy(),'x_pheno_um':xy[:,0],'y_pheno_um':xy[:,1],
                            'node_id':calls.protein__pln_annotation_id.astype(object).to_numpy(),
                            'block':calls.protein__block_relationship.astype(object).to_numpy(),'status':calls.status.to_numpy(),
                            'call':calls.call.where(awarded).to_numpy(),
                            'subset':calls.protein__t_nk_class.where(awarded&calls.call.eq(SUBSET_OF).to_numpy()).to_numpy(),
                            'reason':np.where(awarded,'',why(calls.call.to_numpy(),calls.protein__is_immune.to_numpy(),refused))})[node]
        unknown=sorted(set(cells.node_id.dropna())-set(outlines.index))
        if unknown or cells.node_id.isna().any():raise ValueError(f'node cells without a {NODE_CLASS} outline in the annotations: {unknown[:5]}')
        cells['depth_um']=np.nan
        for nid,g in cells.groupby('node_id'):
            cells.loc[g.index,'depth_um']=depth(g[['x_pheno_um','y_pheno_um']].to_numpy(float),outlines.geometry[nid])
        cells.to_parquet(out/'cells.parquet',index=False)
        # The node's own PhenoCycler cells that this section's Xenium cells overlap, classed by the same rule.
        used=set(read(Path(rec_in['joined']['attempt'])/'overlap_edges.parquet',['object_id']).object_id.astype(str))
        native,_=delivered_calls(spec['root'],section['donor_id'],'nodes')
        native=native[native.object_id.isin(used)&native.pln_annotation_id.astype(str).isin(set(cells.node_id))]
        nc=native[nspec['class_column']].astype(object)
        named=nc.notna()&~nc.isin(refused)
        native=pd.DataFrame({'node_id':native.pln_annotation_id.astype(str).to_numpy(),'call':nc.where(named).to_numpy(),
                             'subset':native.t_nk_class.astype(object).where(named&nc.eq(SUBSET_OF)).to_numpy(),
                             'reason':np.where(named,'',why(nc.to_numpy(),native.is_immune.to_numpy(),refused))})
        frames={'xenium':cells,'native':native}
        comp=[]
        for source,f in frames.items():
            f=f.assign(source=source)
            classed=f[f.call.notna()]
            comp+=[shares(classed,['node_id','source'],'call',classed).rename(columns={'call':'value'}).assign(level='class'),
                   shares(classed[classed.subset.notna()],['node_id','source'],'subset',classed[classed.call.eq(SUBSET_OF)])
                   .rename(columns={'subset':'value'}).assign(level='subset'),
                   shares(f[f.reason.ne('')],['node_id','source'],'reason',f).rename(columns={'reason':'value'}).assign(level='no_class')]
        comp=pd.concat(comp,ignore_index=True)[['node_id','source','level','value','n','share']]
        comp.to_parquet(out/'composition.parquet',index=False)
        classed=cells[cells.call.notna()]
        top=float(np.ceil(max(classed.depth_um.max(),DEPTH_STEP_UM)/DEPTH_STEP_UM)*DEPTH_STEP_UM)
        low=float(np.floor(min(classed.depth_um.min(),0.)/DEPTH_STEP_UM)*DEPTH_STEP_UM)
        edges=np.arange(low,top+DEPTH_STEP_UM,DEPTH_STEP_UM);mids=(edges[:-1]+edges[1:])/2
        profile=classed.assign(bin_mid=pd.cut(classed.depth_um,edges,labels=mids,include_lowest=True).astype(float))
        profile=profile.groupby(['node_id','bin_mid','call']).size().rename('n').reset_index()
        profile.to_parquet(out/'profile.parquet',index=False)
        objects=[]
        for nid,g in cells.groupby('node_id'):
            poly=outlines.geometry[nid];nat=native[native.node_id.eq(nid)]
            objects.append({'node_id':nid,'block':str(g.block.dropna().iloc[0]) if g.block.notna().any() else '',
                            'area_um2':float(shapely.area(poly)),'coverage':float(coverage(np.array([poly]),xy)[0]),
                            'xenium_cells':len(g),'xenium_classed':int(g.call.notna().sum()),
                            'native_cells_overlapped':len(nat),'native_classed':int(nat.call.notna().sum()),
                            'depth_max_um':float(g.depth_um.max()),'outside_outline':int((g.depth_um<0).sum())})
        objects=pd.DataFrame(objects);objects.to_parquet(out/'objects.parquet',index=False)
        # RNA: share of each shared identity group's genes detected (gDNA-corrected), per node, class and subset.
        rna=analysis_input.load_base(Path(rec_in['rna']['attempt'])/'rna.h5ad',donor=section['donor_id'],roi=section['roi'],serial=section['xenium_serial'])
        if not np.array_equal(rna.obs_names.astype(str).to_numpy(),joined.cell_id.astype(str).to_numpy()):raise ValueError('RNA and join rows differ')
        det=rna.layers['detected_gdna']
        rows_=np.flatnonzero(node)[cells.call.notna().to_numpy()]
        det=det.tocsr()[rows_].tocsc() if hasattr(det,'tocsr') else det[rows_]
        keyed=classed.assign(subset=classed.subset.fillna(''))
        rows=[]
        for group,genes in identity_groups('pancreatic_lymph_node').items():
            j=[rna.var_names.get_loc(g) for g in genes if g in rna.var_names]
            if not j:continue
            share=np.asarray(det[:,j].astype(np.int32).sum(axis=1)).ravel()/len(j)
            frame=keyed[['node_id','call','subset']].assign(share=share)
            for by in (['node_id','call'],['node_id','call','subset']):
                for k,g in frame.groupby(by):
                    if len(by)==3 and k[2]=='':continue
                    rows.append({'group':group,'genes_on_panel':len(j),'node_id':k[0],'call':k[1],'subset':k[2] if len(by)==3 else '',
                                 'n_cells':len(g),'mean_share_detected':float(g.share.mean()),'cells_detecting_any':float((g.share>0).mean())})
        rna_t=pd.DataFrame(rows);rna_t.to_parquet(out/'rna.parquet',index=False)
        label='HDL112' if section['donor_id']=='112' else section['donor_id']
        tissue={'pancreatic_lymph_node':'pLN'}.get(section['tissue'],section['tissue'])
        figure(comp,profile,rna_t,objects,f"{label} · {tissue} section · Xenium {section['xenium_serial']} · lymph node",out/'nodes')
        view=lambda source,level:{n:g.set_index('value').share.round(4).to_dict() for n,g in comp[comp.source.eq(source)&comp.level.eq(level)].groupby('node_id')}
        summary={'nodes':objects.set_index('node_id').round(4).to_dict('index'),
                 'xenium_class':view('xenium','class'),'native_class':view('native','class'),
                 'xenium_subset':view('xenium','subset'),'native_subset':view('native','subset'),
                 'xenium_no_class':view('xenium','no_class'),'transform_scale':tf.scale}
        write_json(out/'nodes.json',{**summary,'parameters':inputs['parameters'],
                   'interpretation':'descriptive; one section per donor; no test between donors or disease groups'})
        rec.update(status='validated_complete',finished_utc=now(),summary=summary,
                   validation={'depth_from_the_outline':True,'node_cells_only':True,'rows_equal_join':True,
                               'cells_without_a_class_counted_by_reason':True},
                   outputs=[identity(p) for p in out.iterdir() if p.is_file() and p.name!='receipt.json'])
        write_json(out/'receipt.json',rec);write_json(stage/'current.json',rec)
        print(f"{sid}: {len(objects)} node(s); {len(cells):,} Xenium node cells, {len(classed):,} classed",flush=True)
    except Exception as exc:
        rec.update(status='failed',error=str(exc),traceback=traceback.format_exc(),finished_utc=now())
        write_json(out/'receipt.json',rec);write_json(stage/'last_failure.json',rec);raise


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__,formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('--manifest',required=True,type=Path);p.add_argument('--section',required=True)
    p.add_argument('--out',type=Path,default=DEFAULT_OUT)
    a=p.parse_args();run(a.manifest,a.section,a.out)
