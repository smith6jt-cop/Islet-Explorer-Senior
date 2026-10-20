#!/usr/bin/env python3
"""Summary figures across the same-slide sections, read from current receipts only.

One figure per result, each gathering every section that holds it:

  0 status        the completion matrix
  1 transfer      registration recheck, protein coverage, call status, call purity, context use
  2 association   identity-pair AUROC of the independent protein-RNA association, per section and tissue
  3 islets        composition of the Islet core and of tissue; class share by signed distance to the
                  nearest Islet object (no expansion); per Islet object, what the transfer kept
  4 nodes         per node, the transferred class and T/NK subset beside the node's own PhenoCycler
                  cells; node cells without a class, by the producer's reason; class share by depth
  5 rna_identity  RNA identity-group detection (gDNA-corrected) by transferred class
  6 rna_compare   same-slide against Xenium-only sections, RNA only

Nothing a stage wrote is recomputed, and a section whose stage is unavailable is drawn as
unavailable with its reason, never left out. Descriptive only: one section per donor, run and batch
confounded, and no test between donors or disease groups (maintainer, 2026-10-19: no hypotheses
until the data are properly processed).
"""
import argparse
import hashlib
import json
import os
from pathlib import Path
import re
import sys
from uuid import uuid4

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.colors import PowerNorm,TwoSlopeNorm
from matplotlib.patches import Patch

HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(HERE.parents[1]))
from context import MIN_CLASS_CELLS
from import_alignment import read
from islets_four import ENDOCRINE
from manifest import load_manifest
from workflow import DEFAULT_OUT,section_root,identity,write_json,now,valid_receipt

STAGES=('registration_recheck','joined','association','calls','context','islets','node_summaries')
INK,MUTED,GRID,SURFACE=('#0b0b0b','#52514e','#e4e3df','#fcfcfb')
DONOR={'112':'#0072B2','6347':'#E69F00','6533':'#009E73','6580':'#CC79A7'}
PANC_CLASS={'Beta':'#0072B2','Alpha':'#56B4E9','Delta':'#009E73','Endocrine (ambiguous subtype)':'#9fd3c7',
            'Immune':'#E69F00','Blood vessel':'#D55E00','Lymphatic':'#CC79A7','Nerve':'#c9a227','Exocrine':'#c9c6be'}
NODE_CLASS={'T_cell':'#0072B2','Myeloid':'#E69F00','B_Plasma':'#009E73','Granulocyte':'#CC79A7'}
SUBSET={'T_CD4':'#0072B2','T_CD8':'#56B4E9','T_reg':'#009E73','NK_ILC':'#c9a227','T_NK_unresolved':'#c9c6be'}
#: Why a node cell carries no class, in nodes_four's words, shortened for a legend.
REASON={'not immune':'#7d7a73','no lineage marker named':'#a8a59d','Uncalled':'#cfccc4','immune, no label':'#b9b5ab',
        'protein call unavailable':'#e6e3dc'}
STATUS={'pancreas class':'#0072B2','node class':'#009E73','pancreas: Excluded':'#D55E00','node: not immune':'#7d7a73',
        'node: no lineage named':'#a8a59d','node: Uncalled':'#cfccc4','node: immune, no label':'#b9b5ab',
        'call unavailable':'#e6e3dc','no delivered call':'#f2c94c','no protein partner':'#f3f1ec'}
COVER={'pancreas':'#0072B2','pLN':'#009E73','mixed':'#CC79A7','none':'#f3f1ec'}
PURITY=(('≥ 0.95',.95,1.01,'#08519c'),('0.8–0.95',.8,.95,'#6baed6'),('0.5–0.8',.5,.8,'#c6dbef'),('< 0.5',-.01,.5,'#a8a59d'))
CONTEXT={'':('used as context','#0072B2'),'low_overlap':('overlap < 0.5','#a8a59d'),'low_call_purity':('purity < 0.8','#cfccc4'),
         'unresolved_protein_call':('no class','#7d7a73'),'class_absent_or_small_in_training_fold':('class too small','#e6e3dc')}
CAVEAT='descriptive; one section per donor, run and batch confounded; no test between donors'
#: The DAPI gate (restain/align.py) calls its evaluation windows islets, in a lymph node too.
GATE=re.compile(r'only (\d+)/(\d+) islets \((\d+)%\) were measurable; a verdict needs (\d+)%')


def name(donor):
    return 'HDL112' if donor=='112' else donor


def label(s):
    return f"{name(s['donor_id'])} {'pancreas' if s['roi']=='panc' else 'pLN'}"


def current(path):
    """A stage's current receipt and, when it holds a result, its attempt folder."""
    p=Path(path)
    if not p.exists():return {'status':'not_run'},None
    r=json.loads(p.read_text())
    return r,(Path(r['attempt']) if r.get('status')=='validated_complete' else None)


def short(reason,n=70):
    """The cause of an unavailable stage: upstream prefixes dropped, the DAPI gate's windows named as such."""
    text=str(reason or '')
    g=GATE.search(text)
    if g:return f'DAPI gate: {g[1]} of {g[2]} windows measurable, a verdict needs {g[4]}%'
    head=text.split(';')[0].rsplit(': ',1)[-1]
    return head if len(head)<=n else head[:n-1]+'…'


def why_missing(base,stage):
    """Why a section has no result for a stage: its own reason, or, never run, the join's."""
    r,_=current(base/stage/'current.json')
    if r.get('status')=='not_run':
        j,_=current(base/'joined/current.json')
        return 'not run; '+short(j.get('reason')) if j.get('status')=='scientifically_unavailable' else 'not run'
    return short(r.get('reason'))


def style(ax,grid='x'):
    ax.set_facecolor(SURFACE)
    ax.tick_params(labelsize=8,colors=MUTED)
    for s in ax.spines.values():s.set_color(GRID)
    if grid:ax.grid(axis=grid,color=GRID,lw=.8,zorder=0)
    ax.set_axisbelow(True)


def title(ax,text):
    ax.set_title(text,fontsize=10,color=INK,loc='left')


def share_text(w):
    """A share as a percentage to two significant digits, never in scientific notation."""
    return np.format_float_positional(100*w,precision=2,fractional=False,trim='-')+' %'


def stacked(ax,rows,parts,colours,labels=None,height=.66,y=None,thin=.01,names=None):
    """Horizontal stacked shares, one bar per row; rows is a list of dicts part -> share.

    A part above zero but under ``thin`` is too thin to see, so it is also written beneath its bar with
    its share (``names`` maps a part to its legend label); a legend part neither drawn nor written has
    no cells in that bar."""
    y=np.arange(len(rows))[::-1] if y is None else np.asarray(y)
    names=names or {};written=False
    for yi,r in zip(y,rows):
        left=0.
        if r is None:continue
        for p in parts:
            w=float(r.get(p,0.))
            if w>0:ax.barh(yi,w,left=left,color=colours[p],edgecolor=SURFACE,linewidth=.6,height=height,zorder=2)
            left+=w
        small=[f'{names.get(p,p)} ({share_text(r[p])})' for p in parts if 0<float(r.get(p,0.))<thin]
        if small:
            ax.text(0,yi-height/2-.03,f'under {share_text(thin)}: '+' · '.join(small),va='top',ha='left',
                    fontsize=7,color=MUTED,zorder=3);written=True
    if written:
        lo,hi=ax.get_ylim();ax.set_ylim(min(lo,float(np.min(y))-height/2-.36),hi)
    ax.set_yticks(y)
    ax.set_yticklabels(labels or ['']*len(rows),fontsize=8.5,color=INK)
    ax.set_xlim(0,1);ax.xaxis.set_major_formatter(matplotlib.ticker.PercentFormatter(1.0,decimals=0))
    return y


def legend(ax,colours,parts=None,ncol=3,anchor=(0,-.1)):
    parts=parts or list(colours)
    ax.legend(handles=[Patch(color=colours[p],label=p) for p in parts],frameon=False,fontsize=8,
              loc='upper left',bbox_to_anchor=anchor,ncol=ncol)


def unavailable_text(ax,y,text):
    ax.text(.01,y,f'unavailable: {text}',va='center',ha='left',fontsize=7.5,color=MUTED,style='italic',
            transform=ax.get_yaxis_transform(),zorder=3)


def call_status(summary):
    """Every Xenium cell of a section in one call-status category, from its calls summary."""
    out={k:0 for k in STATUS}
    for key,n in summary['status_by_winning_tissue'].items():
        tissue,status=key.split('|',1)
        if status=='awarded':out['pancreas class' if tissue=='pancreas' else 'node class']+=n
        elif status=='unavailable':out['no protein partner' if tissue=='none' else 'call unavailable']+=n
        elif status=='no_class_delivered':out['no delivered call']+=n
    for key,n in summary.get('not_a_class_by_label',{}).items():
        table,lab,imm=(key.split('|')+['',''])[:3]
        if table=='pancreas':out['pancreas: Excluded' if lab=='Excluded' else 'call unavailable']+=n
        elif lab=='Uncalled':out['node: Uncalled']+=n
        elif imm=='False':out['node: not immune']+=n
        elif imm=='True':out['node: immune, no label']+=n
        else:out['node: no lineage named']+=n
    total=sum(out.values())
    if total!=summary['n_xenium_cells']:raise ValueError(f"call status covers {total:,} of {summary['n_xenium_cells']:,} Xenium cells")
    return {k:v/total for k,v in out.items()}


def reason_key(text):
    """A node cell's reason for carrying no class, as nodes_four writes it, to its legend key."""
    return next((k for k in REASON if str(text).startswith(k)),'protein call unavailable')


def fig_status(table,sections):
    order=[s['section_id'] for s in sections]
    stages=[c for c in dict.fromkeys(table.stage)]
    piv=table.pivot(index='section_id',columns='stage',values='status').reindex(index=order,columns=stages)
    colour={'validated_complete':'#9fd3c7','scientifically_unavailable':'#f2c94c','not_applicable':'#f3f1ec',
            'not_run':'#ffffff','failed':'#D55E00','stale':'#E69F00'}
    text={'validated_complete':'done','scientifically_unavailable':'unavailable','not_applicable':'n/a','not_run':'–',
          'failed':'FAILED','stale':'STALE'}
    fig,ax=plt.subplots(figsize=(15,6.2),facecolor=SURFACE,constrained_layout=True)
    style(ax,grid=None)
    for i,sid in enumerate(order):
        for j,st in enumerate(stages):
            v=piv.loc[sid,st]
            v='not_applicable' if pd.isna(v) else v
            ax.add_patch(plt.Rectangle((j,len(order)-1-i),.96,.9,color=colour.get(v,'#ffffff'),ec=GRID,lw=.6))
            ax.text(j+.48,len(order)-1-i+.45,text.get(v,v),ha='center',va='center',fontsize=7.5,
                    color=INK if v!='not_applicable' else MUTED)
    lab={s['section_id']:f"{label(s)} · {'same slide' if s['design']=='same_slide' else 'Xenium only'} · {s['xenium_serial']}" for s in sections}
    ax.set_xlim(0,len(stages));ax.set_ylim(0,len(order))
    ax.set_xticks(np.arange(len(stages))+.48);ax.set_xticklabels([s.replace('_',' ') for s in stages],rotation=30,ha='right',fontsize=8.5,color=INK)
    ax.set_yticks(np.arange(len(order))+.45);ax.set_yticklabels([lab[s] for s in order][::-1],fontsize=8.5,color=INK)
    ax.tick_params(length=0)
    title(ax,'What every stage holds, from its current receipt')
    fig.suptitle('Same-slide four · completion',x=.01,ha='left',fontsize=11,color=INK)
    return fig


def fig_transfer(sections,root):
    rows={k:[] for k in ('recheck','cover','status','purity','context')}
    notes={}
    for s in sections:
        base=section_root(root,s)
        rc,_=current(base/'registration_recheck/current.json')
        rows['recheck'].append(rc.get('summary') if rc.get('status')=='validated_complete' else None)
        if rc.get('status')!='validated_complete':notes[label(s)]=short(rc.get('reason'))
        j,ja=current(base/'joined/current.json')
        if ja:
            t=json.loads((ja/'join.json').read_text())
            n=t['n_rna_eligible_cells'];rows['cover'].append({k:v/n for k,v in t['cells_by_protein_tissue'].items()})
        else:rows['cover'].append(None)
        c,ca=current(base/'calls/current.json')
        if ca:
            rows['status'].append(call_status(c['summary']))
            pc=read(ca/'protein_calls.parquet',['status','call_purity'])
            p=pc.loc[pc.status.eq('awarded'),'call_purity'].to_numpy(float)
            rows['purity'].append({b:float(((p>=lo)&(p<hi)).mean()) for b,lo,hi,_ in PURITY})
        else:rows['status'].append(None);rows['purity'].append(None)
        cx,_=current(base/'context/current.json')
        if cx.get('status')=='validated_complete':
            fb=cx['summary']['fallback_counts'];n=cx['summary']['n_cells']
            rows['context'].append({CONTEXT.get(k,(k,))[0]:v/n for k,v in fb.items()})
        else:rows['context'].append(None)
    labels=[label(s) for s in sections]
    fig=plt.figure(figsize=(19,10.5),facecolor=SURFACE,constrained_layout=True)
    g=fig.add_gridspec(2,3)
    # (a) registration recheck
    ax=fig.add_subplot(g[0,0]);style(ax)
    y=np.arange(len(sections))[::-1]
    for yi,r,lab in zip(y,rows['recheck'],labels):
        if r is None:unavailable_text(ax,yi,notes.get(lab,''));continue
        ax.plot([0,r['threshold_pass_um']],[yi,yi],color=GRID,lw=6,solid_capstyle='butt',zorder=1)
        ax.plot(r['threshold_pass_um'],yi,marker='|',ms=14,mew=2,color=MUTED,zorder=2)
        ax.plot(r['pair_p50_um'],yi,'o',ms=7,color=DONOR[lab.split()[0].replace('HDL','')],mec=SURFACE,zorder=3)
    ax.set_yticks(y);ax.set_yticklabels(labels,fontsize=8.5,color=INK);ax.set_xlim(0,2.1)
    ax.set_xlabel('median nucleus pair distance after the validated transform (µm); bar = the gate\'s pass threshold',color=MUTED,fontsize=8.5)
    title(ax,'(a) Registration, rechecked on the current nuclei')
    # (b) nuclei matched within 2 µm, against the shifted null
    ax=fig.add_subplot(g[0,1]);style(ax)
    for yi,r,lab in zip(y,rows['recheck'],labels):
        if r is None:unavailable_text(ax,yi,'registration');continue
        ax.plot([r['matched_within_2um_null'],r['matched_within_2um']],[yi,yi],color=GRID,lw=2,zorder=1)
        ax.plot(r['matched_within_2um_null'],yi,'o',ms=7,mfc=SURFACE,mec=MUTED,zorder=2)
        ax.plot(r['matched_within_2um'],yi,'o',ms=7,color=DONOR[lab.split()[0].replace('HDL','')],mec=SURFACE,zorder=3)
    ax.set_yticks(y);ax.set_yticklabels(labels,fontsize=8.5,color=INK);ax.set_xlim(0,1)
    ax.xaxis.set_major_formatter(matplotlib.ticker.PercentFormatter(1.0,decimals=0))
    ax.set_xlabel('Xenium nuclei with a PhenoCycler nucleus within 2 µm (filled) and under the shifted null (open)',color=MUTED,fontsize=8.5)
    title(ax,'(b) Nuclei matched after registration')
    # (c) protein coverage by contributing tissue
    ax=fig.add_subplot(g[0,2]);style(ax)
    cover={'pancreas':'pancreas cells','pLN':'node cells','mixed':'both tissues','none':'no protein'}
    stacked(ax,rows['cover'],list(COVER),COVER,labels,names=cover)
    for yi,r in zip(y,rows['cover']):
        if r is None:unavailable_text(ax,yi,'no join')
    legend(ax,{cover[k]:c for k,c in COVER.items()},ncol=4)
    title(ax,'(c) Xenium cells by the tissue of their protein contributors')
    # (d) call status
    ax=fig.add_subplot(g[1,0]);style(ax)
    parts=[k for k in STATUS if any(r and r.get(k,0)>0 for r in rows['status'])]
    stacked(ax,rows['status'],parts,STATUS,labels)
    for yi,r in zip(y,rows['status']):
        if r is None:unavailable_text(ax,yi,'no join')
    legend(ax,STATUS,parts,ncol=3)
    title(ax,'(d) Xenium cells by transferred call')
    # (e) purity of awarded calls
    ax=fig.add_subplot(g[1,1]);style(ax)
    pc={b:c for b,_,_,c in PURITY}
    stacked(ax,rows['purity'],list(pc),pc,labels)
    for yi,r in zip(y,rows['purity']):
        if r is None:unavailable_text(ax,yi,'no join')
    legend(ax,pc,ncol=4)
    title(ax,'(e) Call purity of awarded cells (overlapped area agreeing with the call)')
    # (f) context use
    ax=fig.add_subplot(g[1,2]);style(ax)
    cc={v[0]:v[1] for v in CONTEXT.values()}
    cc.update({k:'#cccccc' for r in rows['context'] if r for k in r if k not in cc})
    stacked(ax,rows['context'],list(cc),cc,labels)
    for yi,r in zip(y,rows['context']):
        if r is None:unavailable_text(ax,yi,'no join')
    legend(ax,cc,ncol=3)
    title(ax,'(f) Protein-guided RNA context: whose class served as context')
    fig.suptitle(f'Same-slide four · protein to Xenium transfer · {CAVEAT}',x=.01,ha='left',fontsize=11,color=INK)
    return fig


def fig_association(sections,root):
    cols,frames,missing=[],[],[]
    for i,s in enumerate(sections):
        base=section_root(root,s)
        a,aa=current(base/'association/current.json')
        if not aa:
            missing.append(f'{label(s)}: no association ({why_missing(base,"association")})');continue
        p=pd.read_parquet(aa/'pairs.parquet');p=p[p.role.eq('TYPE')]
        for tissue,g in p.groupby('tissue',sort=False):
            where=('pancreas tissue' if tissue=='pancreas' else ('same-block node' if s['roi']=='panc' else 'node'))
            key=f"{name(s['donor_id'])} {'pancreas' if s['roi']=='panc' else 'pLN'} section · {where}"
            cols.append(((tissue!='pancreas',i),key));frames.append(g.assign(column=key))
    order=[k for _,k in sorted(cols)]
    n_panc=sum(1 for (node,_),_ in cols if not node)
    p=pd.concat(frames,ignore_index=True)
    p['pair']=p.marker.astype(str)+' → '+p.gene.astype(str)
    a=p.pivot_table(index='pair',columns='column',values='auroc').reindex(columns=order)
    a=a.loc[a.mean(axis=1).sort_values(ascending=False).index]
    fig=plt.figure(figsize=(17,11.5),facecolor=SURFACE,constrained_layout=True)
    g=fig.add_gridspec(1,2,width_ratios=[3.1,1.1])
    ax=fig.add_subplot(g[0,0]);style(ax,grid=None)
    norm=TwoSlopeNorm(vmin=.25,vcenter=.5,vmax=1.)
    im=ax.imshow(a.to_numpy(float),cmap='RdBu_r',norm=norm,aspect='auto')
    for i in range(a.shape[0]):
        for j in range(a.shape[1]):
            v=a.iat[i,j]
            if np.isfinite(v):ax.text(j,i,f'{v:.2f}',ha='center',va='center',fontsize=7,color=INK if .35<v<.85 else SURFACE)
    ax.axvline(n_panc-.5,color=SURFACE,lw=4)
    ax.set_xticks(range(a.shape[1]));ax.set_xticklabels(a.columns,rotation=30,ha='right',fontsize=8.5,color=INK)
    ax.set_yticks(range(a.shape[0]));ax.set_yticklabels(a.index,fontsize=8.5,color=INK)
    cb=fig.colorbar(im,ax=ax,shrink=.6,pad=.01);cb.ax.tick_params(labelsize=8,colors=MUTED)
    cb.set_label('AUROC: protein rank separates cells with the gene detected (0.5 = no separation)',color=MUTED,fontsize=8.5)
    title(ax,'(a) Identity pairs: does the protein rank the cells where its gene is detected above the rest? '
             '(blank: marker not in that frame\'s alignment)')
    ax=fig.add_subplot(g[0,1]);style(ax)
    strata=p.groupby('column').agg(auroc=('auroc','median'),shifted=('auroc_shifted','median'),n=('n_population','first')).reindex(order)
    y=np.arange(len(strata))[::-1]
    for yi,(k,r) in zip(y,strata.iterrows()):
        ax.plot([r.shifted,r.auroc],[yi,yi],color=GRID,lw=2,zorder=1)
        ax.plot(r.shifted,yi,'o',ms=7,mfc=SURFACE,mec=MUTED,zorder=2)
        ax.plot(r.auroc,yi,'o',ms=7,color=DONOR[k.split()[0].replace('HDL','')],mec=SURFACE,zorder=3)
        ax.text(.995,yi,f'{int(r.n):,} cells',transform=ax.get_yaxis_transform(),ha='right',va='center',fontsize=7.5,color=MUTED)
    ax.axvline(.5,color=MUTED,lw=1,ls=(0,(3,3)))
    ax.set_yticks(y);ax.set_yticklabels(strata.index,fontsize=8,color=INK);ax.set_xlim(.4,1.)
    ax.set_xlabel('median identity-pair AUROC (filled);\nshifted null, protein from the cell 25 µm away (open)',color=MUTED,fontsize=8.5)
    title(ax,'(b) Per section and tissue')
    note='; '.join(missing)
    fig.suptitle(f'Same-slide four · independent protein–RNA association (RNA side: gDNA-corrected detection, no protein value enters) · {CAVEAT}'
                 +(f'\n{note}' if note else ''),x=.01,ha='left',fontsize=10.5,color=INK)
    return fig


def islet_tables(sections,root):
    out={}
    for s in sections:
        if s['roi']!='panc':continue
        r,a=current(section_root(root,s)/'islets/current.json')
        if a:out[s['donor_id']]={k:pd.read_parquet(a/f'{k}.parquet') for k in ('composition','profile','objects','rna')}
    return out


def fig_islets(sections,root):
    t=islet_tables(sections,root);donors=list(t)
    fig=plt.figure(figsize=(19,16),facecolor=SURFACE,constrained_layout=True)
    g=fig.add_gridspec(3,4,height_ratios=[1.25,1,1])
    # (a) composition of the Islet core and of tissue: the approved pipeline's two regions, no expansion
    ax=fig.add_subplot(g[0,:2]);style(ax)
    rows,labels=[],[]
    for d in donors:
        c=t[d]['composition']
        for word,where in (('Islet core',c.region.eq('core')&c.object_class.eq('Islet')),('tissue',c.region.eq('tissue'))):
            rows.append(c[where].set_index('call').fraction.to_dict());labels.append(f'{name(d)} · {word}')
    parts=[k for k in PANC_CLASS if any(k in r for r in rows)]
    stacked(ax,rows,parts,PANC_CLASS,labels,height=.62)
    legend(ax,PANC_CLASS,parts,ncol=5,anchor=(0,-.06))
    title(ax,'(a) Transferred class of pancreas-tissue Xenium cells: Islet core and tissue')
    # (b) objects: all, with Xenium core cells, fully imaged
    ax=fig.add_subplot(g[0,2:]);style(ax)
    kinds=('Islet','EndocrineCluster','EndocrineSingle');ylab,yy,top=[],[],0
    k=0
    for d in donors:
        o=t[d]['objects']
        for kind in kinds:
            ok=o[o.object_class.eq(kind)]
            vals=(len(ok),int((ok.core_cells>0).sum()),int(((ok.coverage>=.95)&(ok.core_cells>0)).sum()))
            for v,c,h in zip(vals,('#e6e3dc','#9ecae1','#2171b5'),(.8,.55,.3)):
                ax.barh(-k,v,height=h,color=c,zorder=2)
            ax.text(vals[0]+8,-k,f'{vals[2]} / {vals[1]} / {vals[0]}',va='center',fontsize=7.5,color=MUTED)
            top=max(top,vals[0])
            ylab.append(f'{name(d)} · {kind}');yy.append(-k);k+=1
        k+=.6
    ax.set_yticks(yy);ax.set_yticklabels(ylab,fontsize=8.5,color=INK);ax.set_xlim(0,top*1.25)
    ax.legend(handles=[Patch(color='#e6e3dc',label='objects'),Patch(color='#9ecae1',label='with Xenium cells in the core'),
                       Patch(color='#2171b5',label='and ≥ 95 % of the outline imaged')],frameon=False,fontsize=8,loc='upper right')
    ax.set_xlabel('islet objects (eroded to whole cells); numbers: fully imaged / with Xenium core cells / all',color=MUTED,fontsize=8.5)
    title(ax,'(b) Islet objects the Xenium region imaged')
    # (c) class share by signed distance, one line per donor
    groups=(('endocrine (Alpha, Beta, Delta, ambiguous)',ENDOCRINE),('Immune',('Immune',)),('Blood vessel',('Blood vessel',)))
    for i,(lab,members) in enumerate(groups):
        ax=fig.add_subplot(g[1,i]);style(ax,grid='both')
        for d in donors:
            pr=t[d]['profile'];tot=pr.groupby('bin_mid').n.sum()
            sh=(pr[pr.call.isin(members)].groupby('bin_mid').n.sum().reindex(tot.index).fillna(0)/tot).sort_index()
            ax.plot(sh.index,100*sh.to_numpy(),color=DONOR[d],lw=2,label=name(d),solid_capstyle='round')
        ax.axvline(0,color=MUTED,lw=1,ls=(0,(3,3)))
        ax.set_xlabel('signed distance to the nearest Islet object (µm; negative inside)',color=MUTED,fontsize=8.5)
        ax.set_ylabel('share of classed Xenium cells',color=MUTED,fontsize=8.5)
        ax.yaxis.set_major_formatter(matplotlib.ticker.PercentFormatter(100,decimals=0))
        title(ax,f'(c{i+1}) {lab}')
        if i==0:ax.legend(frameon=False,fontsize=8)
    ax=fig.add_subplot(g[1,3]);ax.axis('off')
    ax.text(0,.95,'Reading (a): a class under 1 % of a bar is too thin\nto see, so it is written beneath the bar with its\nshare. A class in the legend that is neither drawn\nnor written has no cells in that bar.\n\nReading (c): every cell\'s signed distance to its\nnearest islet object, measured in the PhenoCycler\nframe; cells whose nearest object is an Islet. As\nin the approved pipeline there is no peri-islet\nexpansion: the distance stays continuous, and\nthe dashed line is the outline. Each Xenium cell\nis placed by its own centroid, so where it sits\ndoes not depend on which PhenoCycler cell it\noverlaps.\n\nReading (d): Islet objects with ≥ 95 % of the\noutline imaged, ≥ 5 Xenium cells in the core and\n≥ 5 named PhenoCycler cells. The line is equality.',
            va='top',fontsize=8.5,color=MUTED)
    # (d) what the transfer kept, per Islet object
    for i,d in enumerate(donors):
        ax=fig.add_subplot(g[2,i]);style(ax,grid='both')
        o=t[d]['objects']
        q=o[o.object_class.eq('Islet')&(o.coverage>=.95)&(o.core_awarded>=5)&(o.native_classed>=5)]
        ax.plot([0,1],[0,1],color=MUTED,lw=1,zorder=1)
        ax.scatter(q.native_endocrine_fraction,q.core_endocrine_fraction,s=22,color=DONOR[d],edgecolor=SURFACE,linewidth=.8,zorder=3)
        r=q[['native_endocrine_fraction','core_endocrine_fraction']].corr().iloc[0,1] if len(q)>2 else np.nan
        med=(q.native_endocrine_fraction-q.core_endocrine_fraction).abs().median()
        ax.text(.03,.97,f'n = {len(q)}\nr = {r:.2f}\nmedian |difference| = {med:.3f}',transform=ax.transAxes,va='top',fontsize=8,color=INK)
        ax.set_xlim(-.02,1.02);ax.set_ylim(-.02,1.02)
        for a_ in (ax.xaxis,ax.yaxis):a_.set_major_formatter(matplotlib.ticker.PercentFormatter(1.0,decimals=0))
        ax.set_xlabel('endocrine share, the object\'s own PhenoCycler cells',color=MUTED,fontsize=8.5)
        if i==0:ax.set_ylabel('endocrine share, Xenium cells in its core',color=MUTED,fontsize=8.5)
        title(ax,f'(d) {name(d)} · what the transfer kept')
    fig.suptitle(f'Same-slide four · islets · {CAVEAT}',x=.01,ha='left',fontsize=11,color=INK)
    return fig


def node_tables(sections,root):
    out=[]
    for s in sections:
        base=section_root(root,s)
        r,a=current(base/'node_summaries/current.json')
        if not a:
            if r.get('status')=='scientifically_unavailable':out.append({'section':s,'unavailable':short(r.get('reason'))})
            continue
        t={k:pd.read_parquet(a/f'{k}.parquet') for k in ('composition','profile','objects','rna')}
        for nid in t['objects'].node_id:
            o=t['objects'].set_index('node_id').loc[nid]
            where='pancreas section' if s['roi']=='panc' else 'pLN section'
            out.append({'section':s,'node':nid,'label':f"{name(s['donor_id'])} · {o.block} node {nid[:8]} · {where}",
                        'short':f"{name(s['donor_id'])} · {o.block} node",
                        **{k:v[v.node_id.eq(nid)] for k,v in t.items()},'object':o})
    return out


def fig_nodes(sections,root):
    nodes=node_tables(sections,root);ok=[n for n in nodes if 'node' in n]
    gone=[f"{label(n['section'])}: no node summary ({n['unavailable']})" for n in nodes if 'unavailable' in n]
    fig=plt.figure(figsize=(19,16),facecolor=SURFACE,constrained_layout=True)
    g=fig.add_gridspec(3,len(ok),height_ratios=[1.4,.95,1])
    top=g[0,:].subgridspec(1,2,width_ratios=[1.45,1])
    # (a) and (b): one pair of bars per node, Xenium cells above the node's own PhenoCycler cells
    y=[];labels=[];pos=0.
    for n in ok:
        for word in ('Xenium cells','its PhenoCycler cells'):
            y.append(-pos);labels.append(f"{n['short']} · {word}");pos+=1
        pos+=.5
    for k,(level,colours,head) in enumerate((('class',NODE_CLASS,'(a) Transferred class, among classed cells'),
                                            ('subset',SUBSET,'(b) T/NK subset, among T_cell cells'))):
        ax=fig.add_subplot(top[0,k]);style(ax)
        rows=[]
        for n in ok:
            c=n['composition'];c=c[c.level.eq(level)]
            rows+=[c[c.source.eq(src)].set_index('value').share.to_dict() for src in ('xenium','native')]
        parts=[p for p in colours if any(p in r for r in rows)]
        stacked(ax,rows,parts,colours,labels if k==0 else None,height=.62,y=y)
        legend(ax,colours,parts,ncol=5,anchor=(0,-.05));title(ax,head)
    # (c) every node cell: classed, or why not
    ax=fig.add_subplot(g[1,:]);style(ax)
    rows,labels=[],[]
    for n in ok:
        c=n['composition'];x=c[c.source.eq('xenium')]
        r={'classed':1-float(x[x.level.eq('no_class')].share.sum())}
        for _,v in x[x.level.eq('no_class')].iterrows():r[reason_key(v.value)]=r.get(reason_key(v.value),0.)+float(v.share)
        rows.append(r);o=n['object']
        labels.append(f"{n['short']} · {int(o.xenium_cells):,} cells · {o.coverage:.0%} of outline imaged")
    colours={'classed':'#0072B2',**REASON}
    parts=[p for p in colours if any(p in r for r in rows)]
    stacked(ax,rows,parts,colours,labels,height=.66)
    legend(ax,colours,parts,ncol=6,anchor=(0,-.1))
    title(ax,'(c) Every Xenium node cell: classed, or the producer\'s reason it is not. The node producer names immune '
             'lineages only; typing every node cell is asked of the workstation (handoff item D)')
    # (d) class share by depth, one panel per node; grey: classed cells per bin, relative to the busiest
    for i,n in enumerate(ok):
        ax=fig.add_subplot(g[2,i]);style(ax,grid='both')
        p=n['profile'];tot=p.groupby('bin_mid').n.sum().sort_index()
        ax.fill_between(tot.index,0,100*tot.to_numpy()/tot.max(),step='mid',color=GRID,zorder=0,lw=0)
        for c,col in NODE_CLASS.items():
            sh=(p[p.call.eq(c)].groupby('bin_mid').n.sum().reindex(tot.index).fillna(0)/tot)
            ax.plot(sh.index,100*sh.to_numpy(),color=col,lw=1.8,label=c)
        ax.axvline(0,color=MUTED,lw=1,ls=(0,(3,3)));ax.set_ylim(0,100)
        ax.yaxis.set_major_formatter(matplotlib.ticker.PercentFormatter(100,decimals=0))
        ax.set_xlabel('depth inside the node outline (µm)',color=MUTED,fontsize=8)
        if i==0:
            ax.set_ylabel('share of classed Xenium cells\n(grey: cells per 50 µm bin, of the busiest bin)',color=MUTED,fontsize=8)
            ax.legend(frameon=False,fontsize=7.5,loc='upper right')
        title(ax,f"(d) {n['short']}")
    fig.suptitle(f'Same-slide four · lymph nodes · {CAVEAT}'+(f"\n{'; '.join(gone)}" if gone else ''),x=.01,ha='left',fontsize=11,color=INK)
    return fig


def heat(ax,frame,columns,head,norm):
    """Detection by class: colour and number are the same value, the number in %."""
    h=frame.reindex(columns=columns)
    style(ax,grid=None)
    im=ax.imshow(h.to_numpy(float),cmap='Blues',norm=norm,aspect='auto')
    for i in range(h.shape[0]):
        for j in range(h.shape[1]):
            v=h.iat[i,j]
            if np.isfinite(v):ax.text(j,i,f'{100*v:.1f}',ha='center',va='center',fontsize=6.5,color=INK if norm(v)<.6 else SURFACE)
    ax.set_xticks(range(h.shape[1]));ax.set_xticklabels(h.columns,rotation=35,ha='right',fontsize=7.5,color=INK)
    ax.set_yticks(range(h.shape[0]));ax.set_yticklabels(h.index,fontsize=7.5,color=INK)
    title(ax,head)
    return im


def fig_rna_identity(sections,root):
    islets=islet_tables(sections,root);nodes=[n for n in node_tables(sections,root) if 'node' in n]
    width=max(len(islets),len(nodes))
    panels=[]
    for d,t in islets.items():
        r=t['rna'];r=r.assign(w=r.mean_share_detected*r.n_cells)
        agg=r.groupby(['group','call']).agg(w=('w','sum'),n=('n_cells','sum'))
        size=agg.reset_index().groupby('call').n.max()
        panels.append((0,(agg.w/agg.n).unstack('call'),[c for c in PANC_CLASS if c in size.index and size[c]>=MIN_CLASS_CELLS],
                       f'{name(d)} pancreas'))
    for n in nodes:
        r=n['rna'];r=r[r.subset.eq('')]
        size=r.groupby('call').n_cells.max()
        panels.append((1,r.pivot_table(index='group',columns='call',values='mean_share_detected'),
                       [c for c in NODE_CLASS if c in size.index and size[c]>=MIN_CLASS_CELLS],n['short']))
    fig=plt.figure(figsize=(4.1*width+1.5,12.5),facecolor=SURFACE,constrained_layout=True)
    g=fig.add_gridspec(2,width+1,width_ratios=[1]*width+[.08])
    for row in (0,1):
        mine=[(h,c,t_) for r_,h,c,t_ in panels if r_==row]
        if not mine:continue
        vmax=max(float(np.nanmax(h.reindex(columns=c).to_numpy(float))) for h,c,_ in mine)
        norm=PowerNorm(gamma=.5,vmin=0,vmax=vmax)
        for i,(h,c,t_) in enumerate(mine):im=heat(fig.add_subplot(g[row,i]),h,c,t_,norm)
        cb=fig.colorbar(im,cax=fig.add_subplot(g[row,width]));cb.ax.tick_params(labelsize=8,colors=MUTED)
        cb.ax.yaxis.set_major_formatter(matplotlib.ticker.FuncFormatter(lambda v,_:f'{100*v:.1f}%'))
        cb.set_label('genes of the group detected per cell, mean over the class',color=MUTED,fontsize=8.5)
    fig.suptitle('Same-slide four · RNA identity groups by transferred protein class: % of the group\'s genes detected per cell '
                 f'(gDNA-corrected), mean over the class; classes with ≥ {MIN_CLASS_CELLS} cells; one colour scale per tissue.\n'
                 f'Protein-guided description, reported apart from the independent association · {CAVEAT}',x=.01,ha='left',fontsize=10,color=INK)
    return fig


def fig_rna_compare(sections,root):
    rows=[]
    for s in sections:
        r,a=current(Path(root)/'comparisons'/f"donor_id={s['donor_id']}"/f"roi={s['roi']}"/'current.json')
        if a:rows.append({'label':label(s),'donor':s['donor_id'],'roi':s['roi'],**json.loads((a/'comparison.json').read_text())})
    _,sa=current(Path(root)/'comparisons/summary/current.json')
    qc=pd.read_csv(sa/'summary.csv',dtype={'donor_id':str})
    fig=plt.figure(figsize=(19,6.8),facecolor=SURFACE,constrained_layout=True)
    g=fig.add_gridspec(1,3)
    y=np.arange(len(rows))[::-1]
    ax=fig.add_subplot(g[0,0]);style(ax)
    for off,(k,lab,mk,col) in zip((.2,0,-.2),(('median_log2_count_ratio','raw counts per gene','o','#0072B2'),
                                              ('median_log2_raw_rate_ratio','raw detection rate','s','#56B4E9'),
                                              ('median_log2_rate_ratio','gDNA-corrected detection rate','D','#D55E00'))):
        ax.scatter([r[k] for r in rows],y+off,marker=mk,s=46,color=col,edgecolor=SURFACE,zorder=3,label=lab)
    ax.axvline(0,color=MUTED,lw=1)
    ax.set_yticks(y);ax.set_yticklabels([r['label'] for r in rows],fontsize=8.5,color=INK)
    ax.set_xlabel('median log2 ratio, same-slide over Xenium-only section, across shared genes',color=MUTED,fontsize=8.5)
    ax.legend(frameon=False,fontsize=8,loc='upper left',bbox_to_anchor=(0,-.1),ncol=3)
    title(ax,'(a) Raw counts and raw detection similar or higher; gDNA-corrected detection lower')
    for i,(col,head,xl,logx) in enumerate((('median_gdna_lambda','(b) Median gDNA background per gene (λ)',
                                           'expected background transcripts per gene per cell; line: λ = α = 0.003',True),
                                          ('fraction_critical_count_ge_2','(c) Cells where a gene needs ≥ 2 transcripts to clear background',
                                           'share of QC-eligible cells',False))):
        ax=fig.add_subplot(g[0,i+1]);style(ax)
        for yi,r in zip(y,rows):
            q=qc[qc.donor_id.eq(r['donor'])&qc.roi.eq(r['roi'])].set_index('design')[col]
            ax.plot([q.get('xenium_only'),q.get('same_slide')],[yi,yi],color=GRID,lw=2,zorder=1)
            ax.plot(q.get('xenium_only'),yi,'o',ms=7,mfc=SURFACE,mec=MUTED,zorder=2)
            ax.plot(q.get('same_slide'),yi,'o',ms=7,color=DONOR[r['donor']],mec=SURFACE,zorder=3)
        if logx:ax.set_xscale('log');ax.axvline(.003,color=MUTED,lw=1)
        else:ax.xaxis.set_major_formatter(matplotlib.ticker.PercentFormatter(1.0,decimals=0))
        ax.set_yticks(y);ax.set_yticklabels([r['label'] for r in rows],fontsize=8.5,color=INK)
        ax.set_xlabel(xl+' · filled: same slide, open: Xenium only',color=MUTED,fontsize=8.5);title(ax,head)
    fig.suptitle('Same-slide four · same-slide against Xenium-only sections, RNA only. Different sections of the same donor; run, age '
                 'and sex confounded, so no treatment-loss estimate',x=.01,ha='left',fontsize=10.5,color=INK)
    return fig


def run(manifest_path,root):
    m=load_manifest(manifest_path);root=Path(root)
    same=[s for s in m['sections'] if s['design']=='same_slide']
    same=sorted(same,key=lambda s:(s['roi']!='panc',int(s['donor_id'])))
    paths={f"{s['section_id']}|{st}":section_root(root,s)/st/'current.json' for s in same for st in STAGES}
    paths.update({f"comparison|{s['section_id']}":root/'comparisons'/f"donor_id={s['donor_id']}"/f"roi={s['roi']}"/'current.json' for s in same})
    paths.update({'comparison|summary':root/'comparisons/summary/current.json','completion':root/'completion/current.json'})
    if not valid_receipt(paths['completion']):raise ValueError('no validated completion matrix: run completion_four.py first')
    inputs={'receipts':{k:identity(p) for k,p in paths.items() if p.exists()},'manifest':identity(manifest_path),
            'code':[identity(__file__)]}
    fp=hashlib.sha256(json.dumps(inputs,sort_keys=True).encode()).hexdigest()
    stage=root/'figures';stage.mkdir(parents=True,exist_ok=True)
    if valid_receipt(stage/'current.json',fp):print('figures reused',flush=True);return
    out=stage/'attempts'/(fp[:16]+'-'+uuid4().hex[:8]);out.mkdir(parents=True)
    rec={'status':'running','fingerprint':fp,'inputs':inputs,'started_utc':now(),'attempt':str(out),'job_id':os.environ.get('SLURM_JOB_ID')}
    write_json(out/'receipt.json',rec)
    completion=pd.read_csv(Path(json.loads(paths['completion'].read_text())['attempt'])/'completion.csv',dtype=str)
    order=same+sorted([s for s in m['sections'] if s['design']!='same_slide'],key=lambda s:(s['roi']!='panc',int(s['donor_id'])))
    figures={'0_status':lambda:fig_status(completion,order),'1_transfer':lambda:fig_transfer(same,root),
             '2_association':lambda:fig_association(same,root),'3_islets':lambda:fig_islets(same,root),
             '4_nodes':lambda:fig_nodes(same,root),'5_rna_identity':lambda:fig_rna_identity(same,root),
             '6_rna_compare':lambda:fig_rna_compare(same,root)}
    with PdfPages(out/'same_slide_summary.pdf') as pdf:
        for key,draw in figures.items():
            fig=draw()
            fig.savefig(out/f'{key}.png',dpi=170,facecolor=SURFACE);fig.savefig(out/f'{key}.pdf',facecolor=SURFACE)
            pdf.savefig(fig,facecolor=SURFACE);plt.close(fig)
            print(f'drew {key}',flush=True)
    write_json(out/'figures.json',{'figures':list(figures),'interpretation':CAVEAT,'sections':[s['section_id'] for s in same]})
    rec.update(status='validated_complete',finished_utc=now(),summary={'figures':list(figures)},
               validation={'receipts_only':True,'unavailable_sections_drawn_with_reason':True},
               outputs=[identity(p) for p in out.iterdir() if p.is_file() and p.name!='receipt.json'])
    write_json(out/'receipt.json',rec);write_json(stage/'current.json',rec)
    print(out,flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__,formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('--manifest',required=True,type=Path);p.add_argument('--out',type=Path,default=DEFAULT_OUT)
    a=p.parse_args();run(a.manifest,a.out)
