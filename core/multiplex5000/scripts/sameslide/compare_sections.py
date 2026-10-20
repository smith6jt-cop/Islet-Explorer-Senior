#!/usr/bin/env python3
"""Same-slide Xenium sections against their Xenium-only comparison sections, RNA only.

For each donor and tissue the manifest pairs one same-slide section (PhenoCycler first, Xenium
after) with one Xenium-only section of the same donor and tissue. They are DIFFERENT physical
sections -- different cells, different area, and, as the manifest says, acquisition run, age and
sex confounded -- so nothing here is cell-matched and nothing is a treatment-loss estimate. What
it states is how the two sections' RNA evidence compares, read from each section's own rna and
rna_description stages: QC, the protein-independent detection per gene, the identity genes of
the panel taxonomy (``core/panels``, its one authority), and the controls.

    python compare_sections.py --manifest <json> --donor 112 --roi panc
    python compare_sections.py --manifest <json> --summary
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
from scipy.stats import poisson
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

import evidence
from manifest import load_manifest,canonical_donor
from workflow import DEFAULT_OUT,REPO,section_root,identity,write_json,now,valid_receipt

PANELS=REPO/'core/panels/data'
#: Identity groups per tissue, read from the taxonomy's own files. Subtype hints in them are not
#: read: Exocrine is one label, never acinar or ductal.
IDENTITY={'pancreas':[PANELS/'identity_core',PANELS/'tissues/pancreas/identity_core'],
          'pancreatic_lymph_node':[PANELS/'identity_core']}
VENDOR=('num_cells_detected','region_area','median_genes_per_cell','median_transcripts_per_cell',
        'fraction_transcripts_assigned','decoded_transcripts_per_100um2','fraction_transcripts_decoded_q20',
        'adjusted_genomic_control_probe_rate','estimated_number_of_false_positive_transcripts_per_cell_including_genomic_counts')
DESIGNS=('same_slide','xenium_only')
INK,MUTED,GRID,SURFACE=('#0b0b0b','#52514e','#e4e3df','#fcfcfb')
SERIES=('#2a78d6','#eb6834')


def identity_groups(tissue):
    """group -> genes, one group per taxonomy file; the file stem without its order prefix."""
    out={}
    for root in IDENTITY[tissue]:
        for f in sorted(Path(root).glob('*.csv')):
            out[f.stem.split('_',1)[1]]=sorted(set(pd.read_csv(f).gene.astype(str)))
    return out


def section_qc(attempt,description,bundle):
    """One section's QC row from its own stage outputs and the vendor's metrics."""
    qc=pd.read_parquet(Path(attempt)/'cell_qc.parquet')
    e=qc[qc.qc_eligible]
    norm=json.loads((Path(attempt)/'normalization.json').read_text())
    ctl=json.loads((Path(description)/'control_validation.json').read_text())
    vendor=pd.read_csv(Path(bundle)/'metrics_summary.csv').iloc[0]
    row={'n_cells':len(qc),'n_eligible':len(e),'eligible_fraction':len(e)/len(qc),
         **{f'excluded_{k}':v for k,v in norm['qc']['exclusion_counts'].items()},
         'median_counts':float(e.raw_gene_counts.median()),'median_genes':float(e.raw_detected_genes.median()),
         'median_gdna_lambda':float(e.gdna_lambda.median()),'median_count_density':float(e.count_density.median()),
         'median_cell_area_um2':float(e.cell_area.median()),'median_nucleus_area_um2':float(e.nucleus_area.median()),
         'background_dominated_fraction':float(e.background_dominated.mean()),
         # The fewest transcripts a gene needs in a cell to clear that cell's gDNA background at the
         # nominal alpha -- the same critical count context.py uses. It rises from 1 to 2 where the
         # background rate passes alpha, which moves every corrected rate without any RNA changing.
         'median_critical_count':float(np.median(poisson.isf(evidence.ALPHA,e.gdna_lambda.to_numpy(float))+1)),
         'fraction_critical_count_ge_2':float(np.mean(poisson.isf(evidence.ALPHA,e.gdna_lambda.to_numpy(float))+1>=2)),
         'female_Y_negative_test':ctl['female_Y_negative_test'],
         'held_out_technical_above_target':len(ctl['held_out_technical_probes_above_target']),
         'vendor_fp_relative_error':norm['vendor_check']['relative_error']}
    row.update({f'vendor_{k}':float(vendor[k]) for k in VENDOR if k in vendor and pd.notna(vendor[k])})
    return row


def gene_table(attempts,groups):
    """Per SHARED gene, both designs side by side, with the identity groups it belongs to.

    Only genes both sections measured are compared. The Xenium-only sections carry a custom
    add-on the same-slide sections do not; those genes are returned separately, by design.
    """
    frames=[]
    for design,attempt in attempts.items():
        g=pd.read_parquet(Path(attempt)/'gene_summary.parquet').set_index('gene')
        frames.append(pd.DataFrame({f'{design}__corrected_rate':g.corrected_detection_rate,
            f'{design}__raw_rate':g.raw_detection_rate,f'{design}__mean_count':g.observed_counts/g.n_cells,
            f'{design}__excess_per_cell':g.signed_excess/g.n_cells,f'{design}__below_background':g.below_section_background}))
    one_sided={d:sorted(set(f.index)-set.intersection(*(set(x.index) for x in frames))) for d,f in zip(attempts,frames)}
    t=pd.concat(frames,axis=1,join='inner')
    member={}
    for name,genes in groups.items():
        for gene in genes:member.setdefault(gene,[]).append(name)
    t['identity_groups']=[';'.join(member.get(g,[])) for g in t.index]
    # Ratios only where both sections detect the gene above the nominal background rate: a ratio
    # of two background-level numbers is noise, and a zero is not a level.
    both=(t['same_slide__corrected_rate']>evidence.ALPHA)&(t['xenium_only__corrected_rate']>evidence.ALPHA)
    t['log2_rate_ratio']=np.where(both,np.log2(t['same_slide__corrected_rate']/t['xenium_only__corrected_rate']),np.nan)
    t['log2_count_ratio']=np.where(both,np.log2(t['same_slide__mean_count']/t['xenium_only__mean_count']),np.nan)
    t['log2_raw_rate_ratio']=np.where(both,np.log2(t['same_slide__raw_rate']/t['xenium_only__raw_rate']),np.nan)
    return t.rename_axis('gene').reset_index(),one_sided


def group_table(genes,groups):
    rows=[]
    for name,members in groups.items():
        g=genes[genes.gene.isin(members)]
        rows.append({'group':name,'genes_declared':len(members),'genes_on_panel':len(g),
                     **{f'{d}__median_corrected_rate':float(g[f'{d}__corrected_rate'].median()) if len(g) else np.nan for d in DESIGNS},
                     'genes_with_ratio':int(g.log2_rate_ratio.notna().sum()),
                     'median_log2_rate_ratio':float(g.log2_rate_ratio.median()) if g.log2_rate_ratio.notna().any() else np.nan})
    return pd.DataFrame(rows)


def pair_figure(genes,groups,title,path):
    """Corrected and raw detection side by side: the gDNA background moves only the first."""
    fig,axes=plt.subplots(1,3,figsize=(17,5.4),facecolor=SURFACE,constrained_layout=True,gridspec_kw={'width_ratios':[1,1,1.25]})
    floor=1e-5
    for ax,kind,name in ((axes[0],'corrected','corrected detection rate (clears gDNA background)'),
                         (axes[1],'raw','raw detection rate (at least one transcript)')):
        ax.set_facecolor(SURFACE)
        x=genes[f'xenium_only__{kind}_rate'].clip(lower=floor);y=genes[f'same_slide__{kind}_rate'].clip(lower=floor)
        ax.plot([floor,1],[floor,1],color=MUTED,lw=1,zorder=1)
        ax.scatter(x,y,s=9,color=SERIES[0],alpha=.35,lw=0,zorder=2)
        ax.set_xscale('log');ax.set_yscale('log');ax.set_xlim(floor/2,1.5);ax.set_ylim(floor/2,1.5)
        ax.set_xlabel('Xenium-only section',color=MUTED);ax.set_ylabel('same-slide section',color=MUTED)
        ax.set_title(f'{name}\n{len(genes):,} shared genes; below {floor:g} drawn at {floor:g}',fontsize=9,color=INK,loc='left')
        ax.grid(color=GRID,lw=.8,zorder=0)
    ax=axes[2];ax.set_facecolor(SURFACE)
    g=groups.dropna(subset=['median_log2_rate_ratio']).sort_values('median_log2_rate_ratio').reset_index(drop=True)
    raw=[]
    for _,r in g.iterrows():
        members=genes[genes.identity_groups.str.split(';').apply(lambda s,n=r.group:n in s)]
        raw.append(members.log2_raw_rate_ratio.median())
    ax.axvline(0,color=MUTED,lw=1,zorder=1)
    ax.scatter(g.median_log2_rate_ratio,g.index,s=48,color=SERIES[0],edgecolor=SURFACE,linewidth=1.5,zorder=3,label='corrected detection')
    ax.scatter(raw,g.index,s=48,color=SERIES[1],edgecolor=SURFACE,linewidth=1.5,zorder=3,label='raw detection')
    ax.set_yticks(g.index);ax.set_yticklabels([f'{n}  ({k})' for n,k in zip(g.group,g.genes_with_ratio)],fontsize=8)
    ax.set_xlabel('median log2(same-slide / Xenium-only), genes above background in both',color=MUTED)
    ax.set_title('Identity groups of the panel taxonomy (genes compared)',fontsize=9,color=INK,loc='left')
    ax.grid(axis='x',color=GRID,lw=.8,zorder=0);ax.legend(frameon=False,fontsize=8,loc='lower right')
    for a in axes:
        a.tick_params(colors=MUTED,labelsize=8)
        for s in a.spines.values():s.set_color(GRID)
    fig.suptitle(title+'  ·  different physical sections; acquisition run, age and sex confounded',fontsize=10,color=INK)
    fig.savefig(path.with_suffix('.png'),dpi=160,facecolor=SURFACE);fig.savefig(path.with_suffix('.pdf'),facecolor=SURFACE);plt.close(fig)


def run_pair(manifest_path,donor,roi,root):
    m=load_manifest(manifest_path);donor=canonical_donor(donor)
    sections={s['design']:s for s in m['sections'] if s['donor_id']==donor and s['roi']==roi}
    if set(sections)!=set(DESIGNS):raise ValueError(f'{donor}/{roi}: requires one section of each design')
    receipts={}
    for d,s in sections.items():
        for stage in ('rna','rna_description'):
            p=section_root(root,s)/stage/'current.json'
            if not valid_receipt(p):raise ValueError(f'{p}: not validated')
            receipts[f'{d}__{stage}']=p
    tissue=sections['same_slide']['tissue']
    taxonomy=[f for r in IDENTITY[tissue] for f in sorted(Path(r).glob('*.csv'))]
    inputs={'receipts':{k:identity(p) for k,p in receipts.items()},'taxonomy':[identity(f) for f in taxonomy],
            'code':identity(__file__),'interpretation':m['interpretation']}
    fp=hashlib.sha256(json.dumps(inputs,sort_keys=True).encode()).hexdigest()
    stage=Path(root)/'comparisons'/f'donor_id={donor}'/f'roi={roi}';stage.mkdir(parents=True,exist_ok=True)
    if valid_receipt(stage/'current.json',fp):print(f'{donor}/{roi}: comparison reused',flush=True);return
    out=stage/'attempts'/(fp[:16]+'-'+uuid4().hex[:8]);out.mkdir(parents=True)
    rec={'status':'running','donor_id':donor,'roi':roi,'fingerprint':fp,'inputs':inputs,'started_utc':now(),
         'attempt':str(out),'job_id':os.environ.get('SLURM_JOB_ID')}
    write_json(out/'receipt.json',rec)
    try:
        attempts={d:json.loads(receipts[f'{d}__rna'].read_text())['attempt'] for d in DESIGNS}
        descriptions={d:json.loads(receipts[f'{d}__rna_description'].read_text())['attempt'] for d in DESIGNS}
        qc=pd.DataFrame([{'design':d,'section_id':sections[d]['section_id'],'acquisition_run':sections[d]['acquisition_run'],
                          **section_qc(attempts[d],descriptions[d],sections[d]['bundle'])} for d in DESIGNS])
        qc.to_csv(out/'qc.csv',index=False)
        groups=identity_groups(tissue)
        genes,one_sided=gene_table(attempts,groups)
        if len(genes)!=m['expected_shared_genes'] or genes[[f'{d}__corrected_rate' for d in DESIGNS]].isna().any().any():
            raise ValueError(f"{len(genes)} shared genes; the manifest declares {m['expected_shared_genes']}")
        genes.to_parquet(out/'genes.parquet',index=False)
        grouped=group_table(genes,groups);grouped.to_csv(out/'groups.csv',index=False)
        both=genes[['same_slide__corrected_rate','xenium_only__corrected_rate']].dropna()
        rho=float(both.rank().corr().iloc[0,1])
        label='HDL112' if donor=='112' else donor
        pair_figure(genes,grouped,f'{label} · {tissue}',out/'comparison')
        summary={'genes':len(genes),'genes_with_ratio':int(genes.log2_rate_ratio.notna().sum()),
                 'spearman_corrected_rate_across_genes':rho,
                 'median_log2_rate_ratio':float(genes.log2_rate_ratio.median()),
                 'median_log2_count_ratio':float(genes.log2_count_ratio.median()),
                 'median_log2_raw_rate_ratio':float(genes.log2_raw_rate_ratio.median()),
                 'median_critical_count':dict(zip(qc.design,qc.median_critical_count)),
                 'median_genes':dict(zip(qc.design,qc.median_genes)),'median_counts':dict(zip(qc.design,qc.median_counts)),
                 'eligible_fraction':dict(zip(qc.design,qc.eligible_fraction))}
        write_json(out/'comparison.json',{**summary,'interpretation':m['interpretation'],
                   'genes_in_one_design_only':{d:{'n':len(g),'genes':g} for d,g in one_sided.items()},
                   'ratio_rule':f'genes whose corrected detection rate exceeds the nominal alpha {evidence.ALPHA} in both sections'})
        rec.update(status='validated_complete',finished_utc=now(),summary=summary,
                   validation={'both_rna_stages_validated':True,'shared_genes_equal_declared':True,
                               'not_cell_matched':True,'no_treatment_effect_estimate':True},
                   outputs=[identity(p) for p in out.iterdir() if p.is_file() and p.name!='receipt.json'])
        write_json(out/'receipt.json',rec);write_json(stage/'current.json',rec)
        print(f"{donor}/{roi}: rho {rho:.3f}; median log2 rate ratio {summary['median_log2_rate_ratio']:+.3f}",flush=True)
    except Exception as exc:
        rec.update(status='failed',error=str(exc),traceback=traceback.format_exc(),finished_utc=now())
        write_json(out/'receipt.json',rec);write_json(stage/'last_failure.json',rec);raise


def run_summary(manifest_path,root):
    m=load_manifest(manifest_path)
    pairs=sorted({(s['donor_id'],s['roi']) for s in m['sections']})
    paths={p:Path(root)/'comparisons'/f'donor_id={p[0]}'/f'roi={p[1]}'/'current.json' for p in pairs}
    missing=[f'{d}/{r}' for (d,r),p in paths.items() if not valid_receipt(p)]
    if missing:raise ValueError(f'comparisons not validated: {missing}')
    inputs={'pairs':{f'{d}/{r}':identity(p) for (d,r),p in paths.items()},'code':identity(__file__)}
    fp=hashlib.sha256(json.dumps(inputs,sort_keys=True).encode()).hexdigest()
    stage=Path(root)/'comparisons'/'summary';stage.mkdir(parents=True,exist_ok=True)
    if valid_receipt(stage/'current.json',fp):return
    out=stage/'attempts'/(fp[:16]+'-'+uuid4().hex[:8]);out.mkdir(parents=True)
    rows=[]
    for (d,r),p in paths.items():
        a=Path(json.loads(p.read_text())['attempt'])
        qc=pd.read_csv(a/'qc.csv');c=json.loads((a/'comparison.json').read_text())
        for _,q in qc.iterrows():
            rows.append({'donor_id':d,'roi':r,**q.to_dict(),'spearman_corrected_rate_across_genes':c['spearman_corrected_rate_across_genes'],
                         'median_log2_rate_ratio':c['median_log2_rate_ratio'],'median_log2_count_ratio':c['median_log2_count_ratio'],
                         'median_log2_raw_rate_ratio':c['median_log2_raw_rate_ratio']})
    t=pd.DataFrame(rows);t.to_csv(out/'summary.csv',index=False)
    labels=[f"{'HDL112' if d=='112' else d} {r}" for d,r in pairs]
    fig,axes=plt.subplots(1,4,figsize=(19,4.8),facecolor=SURFACE,constrained_layout=True)
    for ax,(col,name) in zip(axes,(('median_genes','median detected genes per eligible cell'),
                                   ('median_counts','median gene counts per eligible cell'),
                                   ('eligible_fraction','fraction of cells QC-eligible'),
                                   ('median_gdna_lambda','median gDNA background per gene (λ), log'))):
        ax.set_facecolor(SURFACE)
        for colour,design,dy in zip(SERIES,DESIGNS,(-.12,.12)):
            v=[t[(t.donor_id==d)&(t.roi==r)&(t.design==design)][col].iloc[0] for d,r in pairs]
            ax.scatter(v,np.arange(len(pairs))+dy,s=48,color=colour,edgecolor=SURFACE,linewidth=1.5,label=design.replace('_',' '),zorder=3)
        ax.set_yticks(range(len(pairs)));ax.set_yticklabels(labels,fontsize=8);ax.invert_yaxis()
        ax.set_title(name,fontsize=9.5,color=INK,loc='left');ax.grid(axis='x',color=GRID,lw=.8);ax.tick_params(colors=MUTED)
        for s in ax.spines.values():s.set_color(GRID)
    ax=axes[3];ax.set_xscale('log');ax.axvline(evidence.ALPHA,color=MUTED,lw=1,zorder=1)
    ax.set_title(f'median gDNA background per gene (λ), log\nline: λ = α = {evidence.ALPHA:g}; right of it a gene needs 2 transcripts',
                 fontsize=9.5,color=INK,loc='left')
    axes[0].legend(frameon=False,fontsize=8,loc='lower right')
    fig.suptitle('Same-slide vs Xenium-only sections, RNA QC  ·  different sections; run, age and sex confounded',fontsize=10,color=INK)
    fig.savefig(out/'summary.png',dpi=160,facecolor=SURFACE);fig.savefig(out/'summary.pdf',facecolor=SURFACE);plt.close(fig)
    rec={'status':'validated_complete','fingerprint':fp,'inputs':inputs,'finished_utc':now(),'attempt':str(out),
         'job_id':os.environ.get('SLURM_JOB_ID'),'validation':{'all_pairs_validated':True,'pairs':len(pairs)},
         'outputs':[identity(p) for p in out.iterdir() if p.is_file() and p.name!='receipt.json']}
    write_json(out/'receipt.json',rec);write_json(stage/'current.json',rec)
    print(f'comparison summary over {len(pairs)} pairs -> {out}',flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__,formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('--manifest',required=True,type=Path);p.add_argument('--donor');p.add_argument('--roi',choices=['panc','pln_1'])
    p.add_argument('--summary',action='store_true');p.add_argument('--out',type=Path,default=DEFAULT_OUT)
    a=p.parse_args()
    if a.summary:run_summary(a.manifest,a.out)
    else:run_pair(a.manifest,a.donor,a.roi,a.out)
