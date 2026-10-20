#!/usr/bin/env python3
"""Independent RNA control audit, descriptive embedding and spatial summaries."""
import argparse
import hashlib
import json
import os
from pathlib import Path
from uuid import uuid4

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from sklearn.decomposition import TruncatedSVD

import evidence
import analysis_input
from manifest import load_manifest,get_section
from workflow import DEFAULT_OUT,section_root,identity,write_json,now,valid_receipt

SPATIAL_GENES=('INS','GCG','SST','PPY','IAPP','CHGA','CHGB','GHRL','PTPRC','CD3D','CD3E','CD4','CD8A','CD8B',
 'MS4A1','CD79A','CD19','NKG7','GNLY','LST1','CD68','CD163','FOXP3','TCF7','PDCD1','HLA-A','HLA-DRA','PECAM1',
 'VWF','PDPN','COL1A1','COL3A1','KRT8','KRT18','PRSS1','CPA1','KRT19','CFTR')


def run(manifest_path,sid,root):
    s=get_section(load_manifest(manifest_path),sid);base=section_root(root,s)
    source=base/'rna/current.json'
    if not valid_receipt(source):raise ValueError('RNA base is not validated')
    src=Path(json.loads(source.read_text())['attempt'])
    inputs={'rna':identity(source),'code':identity(__file__),'contract':identity(analysis_input.__file__)}
    fp=hashlib.sha256(json.dumps(inputs,sort_keys=True).encode()).hexdigest();base=base/'rna_description';base.mkdir(parents=True,exist_ok=True)
    if valid_receipt(base/'current.json',fp):return
    out=base/'attempts'/(fp[:16]+'-'+uuid4().hex[:8]);out.mkdir(parents=True)
    rec={'status':'running','fingerprint':fp,'inputs':inputs,'attempt':str(out),'started_utc':now(),'job_id':os.environ.get('SLURM_JOB_ID')}
    write_json(out/'receipt.json',rec)
    a=analysis_input.load_base(src/'rna.h5ad',donor=s['donor_id'],roi=s['roi'],serial=s['xenium_serial'])
    controls=pd.read_csv(src/'sex_and_curated_controls.csv')
    y=controls[controls['set'].isin(['Y_clean','Y_crosshyb']) & controls.probe.ne('ANY')]
    female=s['sex']=='F';failed=y[y.status.eq('fail')];missing=y[y.status.eq('missing_probe')]
    tech=pd.read_parquet(src/'technical_controls.parquet')
    intervals=[evidence.wilson(int(k),int(n)) for k,n in zip(tech.corrected_detected_cells,tech.n_cells)]
    tech['detection_ci95_low']=[v[0] for v in intervals];tech['detection_ci95_high']=[v[1] for v in intervals]
    tech['above_nominal_target']=tech.corrected_detection_rate>evidence.ALPHA
    tech.to_parquet(out/'technical_control_audit.parquet',index=False)
    audit={'recorded_sex':s['sex'],'assignment_source':'data/donor_metadata_panc.xlsx; verified by RNA receipt',
       'female_Y_negative_test':'not_applicable' if not female else 'failed' if len(failed) else 'unavailable_missing_probe' if len(missing) else 'passed',
       'nominal_per_gene_target':evidence.ALPHA,'panel_wide_FDR_guarantee':False,
       'failed_Y_probes_and_strata':failed[['set','probe','stratum','corrected_detection_rate','ci95_low','ci95_high']].to_dict('records'),
       'missing_Y_probes':missing.probe.unique().tolist(),
       'sex_consistency':'Y and XIST/TSIX measurements are diagnostics; no expression-based reassignment',
       'male_negative_control_policy':'XIST/TSIX are not validated male negatives',
       'genomic_control_role':'fitted controls; not independent held-out validation',
       'held_out_technical_probes_above_target':tech.loc[tech.evaluation_role.eq('held_out_technical') & tech.above_nominal_target,'gene'].tolist(),
       'claim_restrictions':'Female Y-negative validation fails in the explicitly listed strata; affected specificity claims are blocked.' if female and (len(failed) or len(missing)) else 'No female Y-negative restriction; technical and curated controls remain probe-specific diagnostics.'}
    write_json(out/'control_validation.json',audit)
    keep=a.var.identity_eligible.to_numpy(bool);x=a.X[:,keep]
    rng=np.random.default_rng(0);sample=np.sort(rng.choice(a.n_obs,min(100000,a.n_obs),replace=False))
    model=TruncatedSVD(n_components=20,n_iter=3,random_state=0)
    model.fit(x[sample]);embedding=np.vstack([model.transform(x[i:i+20000]) for i in range(0,a.n_obs,20000)]).astype('float32')
    pd.DataFrame({'cell_id':a.obs_names,**{f'LSI{i+1}':embedding[:,i] for i in range(20)}}).to_parquet(out/'embedding.parquet',index=False)
    np.savez_compressed(out/'embedding_model.npz',components=model.components_,genes=a.var_names[keep].to_numpy(str),
          explained_variance_ratio=model.explained_variance_ratio_,training_rows=sample)
    write_json(out/'embedding.json',{'method':'TruncatedSVD of standard median-depth log1p expression; no background subtraction',
          'fit_n_cells':len(sample),'n_components':20,'random_seed':0,'sex_and_curated_controls_excluded':True,
          'interpretation':'descriptive uncentered expression embedding; no protein labels or disease information used'})
    xy=a.obsm['spatial'];grid=np.floor(xy/100).astype(int)
    cells=pd.DataFrame({'cell_id':a.obs_names,'x_um':xy[:,0],'y_um':xy[:,1],'grid_x':grid[:,0],'grid_y':grid[:,1],
        'raw_counts':a.obs.raw_gene_counts.to_numpy(),'gdna_lambda':a.obs.gdna_lambda.to_numpy()})
    for gene in SPATIAL_GENES:
        if gene in a.var_names:
            j=a.var_names.get_loc(gene);cells[gene+'_observed']=a.layers['counts'][:,j].toarray().ravel()
            cells[gene+'_detected']=a.layers['detected_gdna'][:,j].toarray().ravel()
    grouped=cells.drop(columns=['cell_id']).groupby(['grid_x','grid_y']).sum()
    grouped['n_cells']=cells.groupby(['grid_x','grid_y']).size();grouped['x_um']/=grouped.n_cells;grouped['y_um']/=grouped.n_cells
    grouped.rename(columns={'gdna_lambda':'summed_expected_background_per_gene'},inplace=True)
    grouped.to_parquet(out/'spatial_100um.parquet')
    fig,axes=plt.subplots(2,2,figsize=(11,9),constrained_layout=True)
    for ax,coords,values,title in [(axes[0,0],xy,np.log1p(a.obs.raw_gene_counts),'Spatial / log1p observed gene count'),
        (axes[0,1],xy,a.obs.gdna_lambda,'Spatial / estimated background per gene'),
        (axes[1,0],embedding[:,:2],np.log1p(a.obs.raw_gene_counts),'Expression embedding / log1p count'),
        (axes[1,1],embedding[:,:2],a.obs.gdna_lambda,'Expression embedding / estimated background')]:
        im=ax.hexbin(coords[:,0],coords[:,1],C=values,gridsize=100,mincnt=1,reduce_C_function=np.mean,cmap='viridis')
        ax.set_title(title,fontsize=10);fig.colorbar(im,ax=ax,shrink=.7)
        if coords is xy:ax.set_aspect('equal');ax.invert_yaxis();ax.set_xlabel('µm');ax.set_ylabel('µm')
        else:ax.set_xlabel('SVD 1');ax.set_ylabel('SVD 2')
    label='HDL112' if s['donor_id']=='112' else s['donor_id']
    fig.suptitle(f"{label} · {s['tissue']} · {s['design']} · {a.n_obs:,} eligible cells")
    fig.savefig(out/'rna_qc.png',dpi=180);fig.savefig(out/'rna_qc.pdf');plt.close(fig)
    if not np.array_equal(pd.read_parquet(out/'embedding.parquet',columns=['cell_id']).cell_id.to_numpy(str),a.obs_names.to_numpy(str)):
        raise ValueError('eligible cell population changed')
    rec.update(status='validated_complete',finished_utc=now(),validation={'same_eligible_cell_ids':True,'control_genes_excluded_from_embedding_fit':True},
       outputs=[identity(p) for p in out.iterdir() if p.is_file() and p.name!='receipt.json'])
    write_json(out/'receipt.json',rec);write_json(base/'current.json',rec)
    print(f'{sid}: RNA description complete',flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--manifest',required=True,type=Path);p.add_argument('--section',required=True)
    p.add_argument('--out',type=Path,default=DEFAULT_OUT);a=p.parse_args();run(a.manifest,a.section,a.out)
