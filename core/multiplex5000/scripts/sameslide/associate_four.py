#!/usr/bin/env python3
"""Independent protein-RNA association and conditional detection, one joined section at a time.

The RNA side is the rna stage's protein-independent evidence: raw integer counts and the
gDNA-background Poisson detection (``layers['detected_gdna']``). No protein value, call or label
enters it. The protein side is the joined full-dataset UniFORM value, compared only WITHIN this
section, within one alignment frame and within one tissue: every statistic is a rank or a
within-stratum quantile, so no unit is compared across sections or frames, and a lymph node a
pancreas region happens to image (6580's, 6347's same-block nodes) is its own stratum rather than
a pool of bright immune cells inside the pancreas ranks.

Per panel pair (``vocab.protein_gene_pairs`` on this section's genes):

* ``auroc``          P(protein of an RNA-detected cell > protein of an undetected one), ties half.
* ``auroc_shifted``  the same with each cell's protein read from the cell nearest a point SHIFT_UM
                     away -- what tissue-scale co-variation alone gives. Only the excess over it
                     is evidence at the scale of a cell.
* ``spearman``       protein against raw count.
* conditional detection: the RNA detection rate in each within-section protein decile, Wilson 95 %.

Population, per stratum (frame x tissue): Xenium cells whose protein came from that single frame
and that single tissue, with overlap_fraction at least ``context.MIN_OVERLAP``. Descriptive and
conditional on the registration: it validates neither modality, and nothing here is a target for
any other stage.
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
from scipy.spatial import cKDTree
from scipy.stats import rankdata
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(HERE.parents[1]))
from phenocycler.integration.vocab import protein_gene_pairs
import analysis_input
import evidence
from context import MIN_OVERLAP
from join_four import FRAME_PREFIX,TISSUE
from manifest import load_manifest,get_section
from workflow import DEFAULT_OUT,section_root,identity,write_json,now,valid_receipt

#: The displaced null the registration gate already uses: a partner this far away shares the
#: tissue, not the cell.
SHIFT_UM=25.0
#: A shifted point with no Xenium cell centroid this close has no null partner.
MATCH_UM=5.0
N_BINS=10
INK,MUTED,GRID,SURFACE=('#0b0b0b','#52514e','#e4e3df','#fcfcfb')
SERIES=('#2a78d6','#eb6834')


def auroc(ranks,label):
    """Mann-Whitney AUROC from precomputed average ranks; NaN without both classes."""
    y=np.asarray(label,bool);n1=int(y.sum());n0=len(y)-n1
    if not n1 or not n0:return np.nan
    return float((ranks[y].sum()-n1*(n1+1)/2)/(n1*n0))


def shifted_partner(xy,shift_um=SHIFT_UM,match_um=MATCH_UM):
    """Row of the cell nearest each cell's +x shifted position, or -1 when none is within match_um."""
    d,i=cKDTree(xy).query(xy+np.array([shift_um,0.0]),k=1,distance_upper_bound=match_um)
    return np.where(np.isfinite(d),i,-1)


def spearman(a,b):
    ra,rb=rankdata(a),rankdata(b)
    if np.ptp(ra)==0 or np.ptp(rb)==0:return np.nan
    return float(np.corrcoef(ra,rb)[0,1])


def conditional(value,detected,bins=N_BINS):
    """Detection per within-section quantile bin of the protein value; equal values share a bin."""
    q=pd.qcut(value,bins,labels=False,duplicates='drop')
    rows=[]
    for b in np.unique(q):
        m=q==b;n=int(m.sum());k=int(detected[m].sum());lo,hi=evidence.wilson(k,n)
        rows.append({'bin':int(b)+1,'n':n,'n_detected':k,'rate':k/n,'ci95_low':lo,'ci95_high':hi,
                     'protein_median':float(np.median(value[m]))})
    return rows


def pair_statistics(value,counts,detected,partner):
    """Every statistic for one (marker, gene) over one population; inputs share row order."""
    ok=np.isfinite(value);v=value[ok];c=counts[ok];d=detected[ok]
    out={'n_cells':int(ok.sum()),'n_detected':int(d.sum()),'detection_rate':float(d.mean()) if ok.any() else np.nan,
         'raw_detected_cells':int((c>0).sum())}
    out['auroc']=auroc(rankdata(v),d) if ok.any() else np.nan
    p=partner[ok];have=p>=0
    shifted=np.where(have,value[np.where(have,p,0)],np.nan)
    sh=np.isfinite(shifted)
    out['n_shifted']=int(sh.sum())
    out['auroc_shifted']=auroc(rankdata(shifted[sh]),d[sh]) if sh.any() else np.nan
    out['auroc_excess']=out['auroc']-out['auroc_shifted']
    out['spearman_count']=spearman(v,c) if ok.sum()>2 else np.nan
    return out,(conditional(v,d) if ok.sum()>=N_BINS else [])


def scatter(ax,pairs,heading):
    """Cell against shifted-partner AUROC for one stratum's pairs."""
    ax.set_facecolor(SURFACE)
    ax.plot([0,1],[0,1],color=MUTED,lw=1,zorder=1)
    for colour,(label,mask) in zip(SERIES,(('identity (TYPE)',pairs.role.eq('TYPE')),('state / other',~pairs.role.eq('TYPE')))):
        ax.scatter(pairs.loc[mask,'auroc_shifted'],pairs.loc[mask,'auroc'],s=40,color=colour,edgecolor=SURFACE,linewidth=1.5,zorder=3,label=label)
    placed=[]
    for _,r in pairs.sort_values('auroc_excess',ascending=False).head(8).iterrows():
        # Greedy, in data units: a label whose extent would cross a placed one is dropped, never
        # stacked. ~0.0085 units per character at this size on this axis.
        text=f"{r.marker}/{r.gene}";w=.0085*len(text)
        for side,lo,hi in (('left',r.auroc_shifted+.008,r.auroc_shifted+.008+w),('right',r.auroc_shifted-.008-w,r.auroc_shifted-.008)):
            if all(abs(r.auroc-py)>.03 or hi<plo or lo>phi for plo,phi,py in placed):
                ax.annotate(text,(r.auroc_shifted,r.auroc),xytext=(lo if side=='left' else hi,r.auroc),ha=side,va='center',
                            fontsize=8,color=INK)
                placed.append((lo,hi,r.auroc));break
    ax.set_xlim(.3,1);ax.set_ylim(.3,1);ax.set_xlabel(f'AUROC, protein from the cell {SHIFT_UM:g} µm away',color=MUTED)
    ax.set_ylabel('AUROC, protein of the cell itself',color=MUTED)
    ax.set_title(heading,fontsize=10,color=INK,loc='left')
    ax.grid(color=GRID,lw=.8);ax.legend(frameon=False,fontsize=8,loc='lower right')
    for s in ax.spines.values():s.set_color(GRID)


def figure(pairs,conditional_rows,title,path):
    """One scatter per stratum (largest first, at most four); deciles for the first stratum's top identity pairs."""
    strata=pairs.groupby(['frame','tissue'],sort=False).n_population.first().sort_values(ascending=False).head(4)
    frame,tissue=strata.index[0]
    main=pairs[(pairs.frame==frame)&(pairs.tissue==tissue)]
    panels=main[main.role.eq('TYPE')].sort_values('auroc_excess',ascending=False).head(12)
    ncol=4;nrow=1+int(np.ceil(len(panels)/ncol));span=2 if len(strata)<=2 else 1
    fig=plt.figure(figsize=(14,3.2*nrow+(.8 if len(strata)>1 else 0)),facecolor=SURFACE,constrained_layout=True)
    grid=fig.add_gridspec(nrow,ncol)
    for k,((f,t),n) in enumerate(strata.items()):
        scatter(fig.add_subplot(grid[0,k*span:(k+1)*span]),pairs[(pairs.frame==f)&(pairs.tissue==t)],
                f'{t} cells, {f} frame (n = {n:,})')
    note=(f"{title}\n{len(main)} marker-gene pairs on this panel; strata are frame x tissue, overlap ≥ {MIN_OVERLAP:g}\n"
          "RNA detection: gDNA-background Poisson test (protein-independent)\n"
          "points above the diagonal carry cell-level association")
    if len(strata)==1:
        ax=fig.add_subplot(grid[0,2:]);ax.axis('off');ax.text(0,.95,note,va='top',fontsize=9,color=INK)
    else:
        fig.suptitle(note,x=.01,ha='left',fontsize=9,color=INK)
    cond=pd.DataFrame(conditional_rows)
    for k,(_,r) in enumerate(panels.iterrows()):
        a=fig.add_subplot(grid[1+k//ncol,k%ncol]);a.set_facecolor(SURFACE)
        c=cond[(cond.frame==r.frame)&(cond.tissue==r.tissue)&(cond.marker==r.marker)&(cond.gene==r.gene)].sort_values('bin')
        a.fill_between(c.bin,c.ci95_low,c.ci95_high,color=SERIES[0],alpha=.1,lw=0)
        a.plot(c.bin,c.rate,color=SERIES[0],lw=2,solid_capstyle='round')
        a.set_title(f"{r.marker} → {r.gene}   AUROC {r.auroc:.2f} (shifted {r.auroc_shifted:.2f})",fontsize=8.5,color=INK,loc='left')
        a.set_xticks(range(1,N_BINS+1));a.set_xlabel(f'protein decile within {tissue} cells',fontsize=8,color=MUTED)
        a.yaxis.set_major_formatter(matplotlib.ticker.PercentFormatter(1.0,decimals=1))
        a.set_ylabel('cells with RNA detected',fontsize=8,color=MUTED);a.grid(color=GRID,lw=.8);a.tick_params(labelsize=7,colors=MUTED)
        for s in a.spines.values():s.set_color(GRID)
    fig.savefig(path.with_suffix('.png'),dpi=160,facecolor=SURFACE);fig.savefig(path.with_suffix('.pdf'),facecolor=SURFACE);plt.close(fig)


def run(manifest_path,sid,root):
    section=get_section(load_manifest(manifest_path),sid);base=section_root(root,section)
    joined_path,rna_path=base/'joined/current.json',base/'rna/current.json'
    for p in (joined_path,rna_path):
        if not valid_receipt(p):raise ValueError(f'{p}: not validated')
    joined=json.loads(joined_path.read_text());rna=json.loads(rna_path.read_text())
    inputs={'joined':identity(joined_path),'rna':identity(rna_path),'code':identity(__file__),
            'parameters':{'shift_um':SHIFT_UM,'match_um':MATCH_UM,'bins':N_BINS,'min_overlap':MIN_OVERLAP}}
    fp=hashlib.sha256(json.dumps(inputs,sort_keys=True).encode()).hexdigest()
    stage=base/'association';stage.mkdir(parents=True,exist_ok=True)
    if valid_receipt(stage/'current.json',fp):print(f'{sid}: association reused',flush=True);return
    out=stage/'attempts'/(fp[:16]+'-'+uuid4().hex[:8]);out.mkdir(parents=True)
    rec={'status':'running','section_id':sid,'fingerprint':fp,'inputs':inputs,'started_utc':now(),
         'attempt':str(out),'job_id':os.environ.get('SLURM_JOB_ID')}
    write_json(out/'receipt.json',rec)
    try:
        a=analysis_input.load_base(Path(rna['attempt'])/'rna.h5ad',donor=section['donor_id'],roi=section['roi'],serial=section['xenium_serial'])
        df=pd.read_parquet(Path(joined['attempt'])/'joined_cells.parquet')
        if not np.array_equal(df.cell_id.astype(str).to_numpy(),a.obs_names.astype(str).to_numpy()):
            raise ValueError('joined rows and RNA rows differ')
        alignment=json.loads(Path(joined['inputs']['receipts']['alignment']['path']).read_text())
        frames=json.loads((Path(alignment['attempt'])/'frames.json').read_text())
        pairs=[p for p in protein_gene_pairs(list(a.var_names)) if p.genes]
        genes=sorted({g for p in pairs for g in p.genes})
        j=[a.var_names.get_loc(g) for g in genes]
        counts=a.layers['counts'][:,j].toarray();detected=a.layers['detected_gdna'][:,j].toarray().astype(bool)
        col={g:k for k,g in enumerate(genes)}
        xy=df[['x_um','y_um']].to_numpy(float)
        rows,cond=[],[]
        for frame,prefix in FRAME_PREFIX.items():
            spec=frames[frame]['markers']
            qc=spec if frame=='pancreas_cohort' else next(iter(spec.values()),{})
            for tissue in TISSUE.values():
                pop=(df.protein_frame.eq(frame)&df.protein_tissue.eq(tissue)&(df.overlap_fraction>=MIN_OVERLAP)).to_numpy()
                if not pop.any():continue
                idx=np.flatnonzero(pop);partner=shifted_partner(xy[idx])
                for p in pairs:
                    column=prefix+p.protein
                    if column not in df:continue
                    value=df[column].to_numpy(float)[idx]
                    for g in p.genes:
                        stats,curve=pair_statistics(value,counts[idx,col[g]],detected[idx,col[g]],partner)
                        q=qc.get(p.protein,{})
                        rows.append({'frame':frame,'tissue':tissue,'marker':p.protein,'gene':g,'role':p.role,'n_population':len(idx),
                                     'marker_over_bar':q.get('over_bar'),'marker_provenance':q.get('provenance'),**stats})
                        cond+= [{'frame':frame,'tissue':tissue,'marker':p.protein,'gene':g,**c} for c in curve]
        if not rows:raise ValueError('no pair could be evaluated')
        table=pd.DataFrame(rows);table.to_parquet(out/'pairs.parquet',index=False);table.to_csv(out/'pairs.csv',index=False)
        pd.DataFrame(cond).to_parquet(out/'conditional.parquet',index=False)
        label='HDL112' if section['donor_id']=='112' else section['donor_id']
        figure(table,cond,f"{label} · {section['tissue']} · Xenium {section['xenium_serial']}",out/'association')
        strata=[]
        for (f,t),g in table.groupby(['frame','tissue'],sort=False):
            ty=g[g.role.eq('TYPE')]
            strata.append({'frame':f,'tissue':t,'n_population':int(g.n_population.iloc[0]),'n_pairs':len(g),'type_pairs':len(ty),
                           'type_pairs_auroc_median':float(ty.auroc.median()),'type_pairs_auroc_shifted_median':float(ty.auroc_shifted.median()),
                           'type_pairs_excess_positive':int((ty.auroc_excess>0).sum())})
        strata.sort(key=lambda s:-s['n_population'])
        # The headline is the section's own tissue; another tissue it images is its own stratum.
        own='pancreas' if section['tissue']=='pancreas' else 'pLN'
        head=next((s for s in strata if s['tissue']==own),strata[0])
        summary={'n_pairs':len(table),'frames':sorted(table.frame.unique()),'strata':strata,
                 'headline_stratum':{'frame':head['frame'],'tissue':head['tissue']},
                 **{k:head[k] for k in ('type_pairs_auroc_median','type_pairs_auroc_shifted_median','type_pairs_excess_positive','type_pairs')}}
        write_json(out/'association.json',{**summary,'parameters':inputs['parameters'],
                   'interpretation':'descriptive, within section, frame and tissue; conditional on registration; validates neither modality'})
        rec.update(status='validated_complete',finished_utc=now(),summary=summary,
                   validation={'rows_equal_rna_eligible':True,'rna_evidence_protein_independent':True,
                               'within_section_ranks_only':True,'pairs_evaluated':len(table)},
                   outputs=[identity(p) for p in out.iterdir() if p.is_file() and p.name!='receipt.json'])
        write_json(out/'receipt.json',rec);write_json(stage/'current.json',rec)
        print(f"{sid}: {len(table)} pairs; TYPE AUROC median {summary['type_pairs_auroc_median']:.3f} "
              f"(shifted {summary['type_pairs_auroc_shifted_median']:.3f})",flush=True)
    except Exception as exc:
        rec.update(status='failed',error=str(exc),traceback=traceback.format_exc(),finished_utc=now())
        write_json(out/'receipt.json',rec);write_json(stage/'last_failure.json',rec);raise


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__,formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('--manifest',required=True,type=Path);p.add_argument('--section',required=True)
    p.add_argument('--out',type=Path,default=DEFAULT_OUT)
    a=p.parse_args();run(a.manifest,a.section,a.out)
