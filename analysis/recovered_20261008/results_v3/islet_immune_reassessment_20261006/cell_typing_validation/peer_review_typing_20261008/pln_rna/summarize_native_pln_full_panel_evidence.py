#!/usr/bin/env python3
"""Source-qualified pLN stage 1/2/3 review; no nearest-reference relabels."""
from pathlib import Path
import sys,json
import numpy as np,pandas as pd
ROOT=Path('/home/smith6jt/IO60panc2nd/Islet-Explorer-Senior-v3');sys.path.insert(0,str(ROOT))
from analysis.xenium.scripts.validate_rna_external_lineage_reference import broad,sha
from prepare_native_pln_full_panel_evidence import PROGRAMS
OUT=Path(__file__).resolve().parent;RNA=OUT.parent/'rna';IMMUNE={'Immune_T','Immune_NK_ILC','Immune_B_Plasma','Myeloid','Granulocyte'};ER={'XBP1','FKBP11','DERL3'}
def maptype(x):return None if pd.isna(x) or not str(x) else broad(str(x))
def decorate(meta,pred,fold=None,label='cluster_lineage',key='partition_key'):
 a=meta.set_index(key).copy()
 for mode,short in [('full_measured_panel','full'),('native_naming_basis_removed','naming_removed'),('shared_ER_removed','ER_removed')]:
  p=pred[pred.comparison.eq(mode)].set_index('key').loc[a.index];assert len(p)==len(a)
  a[short+'_reference_category']=p.source_type;a[short+'_reference_broad_context']=p.source_type.map(maptype);a[short+'_default_pruned']=p.pruned_source_type.isna();a[short+'_delta_next']=p.delta_next
 if fold is not None:
  for mode,short in [('full_measured_panel','full'),('native_naming_basis_removed','naming_removed')]:
   f=fold[fold.comparison.eq(mode)].copy();f['broad']=f.source_type.map(maptype);f['all_ref']=f.key.map(a[short+'_reference_broad_context']);f['matches']=f.broad.eq(f.all_ref)
   z=f.groupby('key').agg(n_donor_folds=('key','size'),n_distinct_broad_contexts=('broad','nunique'),n_all_ref_matches=('matches','sum'),n_pruned=('pruned_source_type',lambda x:int(x.isna().sum())))
   for c in z:a[short+'_holdout_'+c]=z[c].loc[a.index]
 a['same_reference_context_across_gene_sensitivities']=a.full_reference_broad_context.eq(a.naming_removed_reference_broad_context)&a.full_reference_broad_context.eq(a.ER_removed_reference_broad_context)
 a['independent_same_cell_accuracy_established']=False;a['reference_context_only']=True;a['native_call_modified']=False
 return a

def program_table(meta,z,key,label):
 rows=[];gs=z['genes'];gi={g:i for i,g in enumerate(gs)};mean=z['mean_lognorm'];pos=z['n_raw_positive'];raw=z['sum_raw_counts'];norm=z['n_normalized_positive'];assert np.array_equal(meta[key].to_numpy(),z['partition_keys'] if key=='partition_key' else z['stratum_keys'])
 for sid,idx in meta.groupby('section_id').indices.items():
  nn=meta.n_cells.to_numpy()[idx];totaln=int(nn.sum());tp=pos[idx].sum(axis=0);tn=norm[idx].sum(axis=0);ts=(mean[idx]*nn[:,None]).sum(axis=0)
  for i in idx:
   info=meta.iloc[i];rest=totaln-int(info.n_cells)
   for role,genes in PROGRAMS.items():
    for g in genes:
     j=gi.get(g);measured=j is not None and np.isfinite(mean[i,j]);r={'key':info[key],'section_id':sid,'donor_id':info.donor_id,'source_serial':info.source_serial,'cluster_id':info.cluster_id,'native_label':info[label],'n_cells':int(info.n_cells),'gene':g,'gene_role':role,'gene_measured':measured,'mean_lognorm':float(mean[i,j]) if measured else None,'n_raw_positive':int(pos[i,j]) if measured else None,'sum_raw_counts':int(raw[i,j]) if measured else None,'raw_positive_fraction':float(pos[i,j]/info.n_cells) if measured else None,'historical_normalized_positive_fraction':float(norm[i,j]/info.n_cells) if measured else None,'rest_of_section_raw_positive_fraction':float((tp[j]-pos[i,j])/rest) if measured and rest else None,'rest_of_section_historical_normalized_positive_fraction':float((tn[j]-norm[i,j])/rest) if measured and rest else None,'rest_of_section_mean_lognorm':float((ts[j]-mean[i,j]*info.n_cells)/rest) if measured and rest else None,'context_not_exclusive_identity':True}
     rows.append(r)
 return pd.DataFrame(rows)
def main():
 prep=json.load(open(OUT/'full_panel_preparation_manifest.json'));meta=pd.read_csv(OUT/'pln_native_partition_metadata.csv',dtype={'section_id':str,'source_serial':str,'donor_id':str});cov=pd.read_csv(OUT/'native_all27_partition_covariates.csv');pred=pd.read_csv(OUT/'SingleR_native_partition_predictions.csv');fold=pd.read_csv(OUT/'SingleR_native_reference_donor_perturbations.csv');a=decorate(meta,pred,fold)
 a['full_reference_stage1_immune_context']=a.full_reference_broad_context.isin(IMMUNE);a['native_stage2_named_immune']=a.cluster_lineage.isin(IMMUNE);a['full_reference_stage2_context_matches']=a.full_reference_broad_context.eq(a.cluster_lineage)&a.native_stage2_named_immune
 a['stage2_evidence_status']='native_ambiguous_retained'
 named=a.native_stage2_named_immune
 a.loc[named,'stage2_evidence_status']='reference_context_concordant_with_sensitivity'
 a.loc[named&~a.full_reference_stage2_context_matches,'stage2_evidence_status']='reference_context_discordant_requires_program_review'
 a.loc[named&a.full_default_pruned,'stage2_evidence_status']='reference_default_pruned_requires_program_review'
 stable=named&a.full_reference_stage2_context_matches&a.same_reference_context_across_gene_sensitivities&~a.full_default_pruned&~a.naming_removed_default_pruned&a.full_holdout_n_distinct_broad_contexts.eq(1)&a.naming_removed_holdout_n_distinct_broad_contexts.eq(1)
 a.loc[stable,'stage2_evidence_status']='reference_context_concordant_across_declared_sensitivities'
 non=a.stage1.eq('Non_immune');a.loc[non&~a.full_reference_stage1_immune_context,'stage2_evidence_status']='coarse_nonimmune_context_no_structural_subtype_claim'
 a.loc[non&a.full_reference_stage1_immune_context,'stage2_evidence_status']='stage1_nonimmune_reference_immune_context_program_review'
 a=a.join(cov.set_index('partition_key')[[c for c in cov if c not in a and c!='partition_key']]);a.reset_index().to_csv(OUT/'all27_stage1_stage2_partition_reference_evidence.csv',index=False);a.reset_index().to_parquet(OUT/'all27_stage1_stage2_partition_reference_evidence.parquet',index=False)
 a.groupby(['stage1','cluster_lineage','stage2_evidence_status'],dropna=False).agg(partitions=('n_cells','size'),QC_cells=('n_cells','sum'),donors=('donor_id','nunique'),median_partition_depth=('depth_q50','median'),median_partition_cell_area=('cell_area_q50','median'),median_multinucleus_fraction=('fraction_multinucleus','median')).reset_index().to_csv(OUT/'stage1_stage2_reference_evidence_summary.csv',index=False)
 a.groupby(['cluster_lineage','full_reference_broad_context']).agg(partitions=('n_cells','size'),QC_cells=('n_cells','sum')).reset_index().to_csv(OUT/'stage2_native_reference_context_cross_counts.csv',index=False)
 z=np.load(OUT/'pln_native_all_gene_profiles.npz');program=program_table(meta,z,'partition_key','cluster_lineage');program.to_parquet(OUT/'all27_partition_program_and_competitor_evidence.parquet',index=False)
 # Old native naming admissions, without pretending shared secretory state is
 # unique B/plasma identity. All source admission markers are retained.
 old=pd.read_parquet('/data/xenium_typing/pln_20261003/clusters/cluster_gene_stats.parquet');admissions=[]
 for core,file in {'Immune_T':'03_immune_T_cell.csv','Immune_B_Plasma':'04_immune_B_plasma.csv','Immune_NK_ILC':'05_immune_NK_ILC.csv','Myeloid':'06_immune_myeloid.csv','Granulocyte':'07_immune_granulocyte.csv'}.items():
  c=pd.read_csv(ROOT/'core/panels/data/identity_core'/file);names=set(c.loc[c.evidence.isin(['canonical','both']),'gene'])
  for _,r in meta[meta.cluster_lineage.eq(core)].iterrows():
   p=old[(old.section_key.astype(str)==r.section_id)&old.cluster_id.eq(r.cluster_id)&old.standout];ng=set(p.gene)&names;admissions.append({'partition_key':r.partition_key,'native_stage2_label':core,'n_cells':int(r.n_cells),'source_naming_standout_genes':';'.join(sorted(ng)),'n_naming_standouts':len(ng),'only_shared_ER_naming_standouts':bool(ng and not ng-ER),'non_ER_naming_standouts':';'.join(sorted(ng-ER))})
 pd.DataFrame(admissions).to_csv(OUT/'native_stage2_naming_admission_evidence.csv',index=False)
 sm=pd.read_csv(OUT/'pln_stage3_stratum_metadata.csv',dtype={'section_id':str,'source_serial':str,'donor_id':str});sp=pd.read_csv(OUT/'SingleR_native_stage3_stratum_predictions.csv');ss=sm[sm.stage2_lineage.isin(['Immune_T','Immune_NK_ILC'])].copy();sa=decorate(ss,sp,label='stage3_label',key='stratum_key')
 sa['native_stage3_broad_identity']=np.where(sa.stage3_label.str.startswith('T_')&~sa.stage3_label.eq('T_NK_unresolved'),'Immune_T',np.where(sa.stage3_label.eq('NK_ILC'),'Immune_NK_ILC','Unresolved_T_NK'))
 sa['full_reference_broad_context_matches_stage3_named_identity']=sa.full_reference_broad_context.eq(sa.native_stage3_broad_identity)
 sa['stage3_evidence_status']=np.where(sa.stage3_label.eq('T_NK_unresolved'),'native_unresolved_retained',np.where(sa.full_reference_broad_context_matches_stage3_named_identity,'aggregate_reference_broad_context_concordant_subtype_unvalidated','aggregate_reference_context_discordant_individual_cells_not_adjudicated'))
 sa['Treg_function_or_stable_lineage_established']=False;sa['NK_vs_ILC_subtype_established']=False
 sc=pd.read_csv(OUT/'native_all27_stage3_stratum_covariates.csv');sa=sa.join(sc.set_index('stratum_key')[[c for c in sc if c not in sa and c!='stratum_key']]);sa.reset_index().to_csv(OUT/'all27_stage3_stratum_reference_evidence.csv',index=False)
 sa.groupby(['stage3_label','stage3_evidence_status','full_reference_broad_context'],dropna=False).agg(strata=('n_cells','size'),native_stage3_cells=('n_cells','sum')).reset_index().to_csv(OUT/'stage3_reference_context_summary.csv',index=False)
 sz=np.load(OUT/'pln_stage3_stratified_all_gene_profiles.npz');program_table(sm,sz,'stratum_key','stage3_label').to_parquet(OUT/'all27_stage3_stratum_program_evidence.parquet',index=False)
 # Exact same-cell context distributions: no new classifier, gene-count gate,
 # truth labels, or redistribution of sparse/unresolved cells.
 controls=[];dist=[];cokeys=[]
 for integ in prep['section_integrity']:
  d=pd.read_parquet(integ['stage3_cell_evidence_path']);assert len(d)==integ['n_stage3_T_NK_cells'];d['stratum_key']='pLN|'+d.section_id+'|'+d.source_serial+'|'+d.cluster_id+'|'+d.lineage
  def raw(g):return d['raw__'+g].to_numpy() if 'raw__'+g in d else np.full(len(d),np.nan)
  def nd(g):return raw(g)*10000/np.maximum(d.depth.to_numpy(),1)>=1
  d['raw_TRAC_positive']=raw('TRAC')>0;d['norm_TRAC_positive']=nd('TRAC');d['norm_TCR_support_coexpressed']=nd('TRAC')&(nd('CD5')|nd('CD6')|nd('THEMIS'))
  d['raw_treg_support_any']=(raw('IL2RA')>0)|(raw('CTLA4')>0)|(raw('IKZF2')>0);d['norm_FOXP3_IL2RA_CTLA4']=nd('FOXP3')&nd('IL2RA')&nd('CTLA4');d['raw_FOXP3_IL2RA_CTLA4']=(raw('FOXP3')>0)&(raw('IL2RA')>0)&(raw('CTLA4')>0)
  d['norm_myeloid_alternative']=nd('CSF1R')&(nd('SPI1')|nd('CD68')|nd('CD14'));d['norm_B_alternative']=nd('PAX5')&(nd('CD79A')|nd('MS4A1'))
  d['norm_NK_receptor_any']=np.logical_or.reduce([nd(g) for g in PROGRAMS['NK_receptor_context'] if 'raw__'+g in d]);d['norm_shared_cytotoxic_any']=np.logical_or.reduce([nd(g) for g in PROGRAMS['shared_cytotoxic_state'] if 'raw__'+g in d]);d['norm_ILC_Th_shared_any']=np.logical_or.reduce([nd(g) for g in PROGRAMS['ILC_Th_shared_context'] if 'raw__'+g in d])
  flags=['raw_TRAC_positive','norm_TRAC_positive','norm_TCR_support_coexpressed','raw_treg_support_any','norm_FOXP3_IL2RA_CTLA4','raw_FOXP3_IL2RA_CTLA4','norm_myeloid_alternative','norm_B_alternative','norm_NK_receptor_any','norm_shared_cytotoxic_any','norm_ILC_Th_shared_any']
  for key,f in d.groupby('stratum_key',sort=True):
   r={'stratum_key':key,'section_id':integ['section_id'],'source_serial':integ['source_serial'],'cluster_id':f.cluster_id.iloc[0],'native_stage3_label':f.lineage.iloc[0],'n_cells':len(f),'median_depth':float(f.depth.median()),'median_cell_area':float(f.cell_area.median()),'fraction_multinucleus':float(f.nucleus_count.gt(1).mean()),'n_native_Treg_support_zero':int(f.n_treg_support_detected.eq(0).sum())}
   for flag in flags:r['n_'+flag]=int(f[flag].sum());r['fraction_'+flag]=float(f[flag].mean())
   controls.append(r)
  d['depth_band']=pd.cut(d.depth,bins=[0,25,50,100,200,500,1000,np.inf],right=False).astype(str)
  d['within_section_TNK_cell_area_quartile']=pd.qcut(d.cell_area,4,labels=False,duplicates='drop')
  for (lab,depth,area),f in d.groupby(['lineage','depth_band','within_section_TNK_cell_area_quartile'],observed=True):
   r={'section_id':integ['section_id'],'source_serial':integ['source_serial'],'native_stage3_label':lab,'depth_band':depth,'within_section_TNK_cell_area_quartile':int(area),'n_cells':len(f),'n_TRAC_normpositive':int(f.norm_TRAC_positive.sum()),'n_Treg_support_zero':int(f.n_treg_support_detected.eq(0).sum())};dist.append(r)
  # The exact compound flags can be independently read back by source-qualified
  # cell key, without changing any native identity or membership field.
  fpath=OUT/'stage3_cell_context_flags'/('section='+integ['section_id'])/'data.parquet';fpath.parent.mkdir(parents=True,exist_ok=True);d[['serial_qualified_cell_key','cell_id','section_id','source_serial','cluster_id','lineage']+flags].to_parquet(fpath,index=False);cokeys.append({'path':str(fpath),'sha256':sha(fpath),'n_cells':len(d)})
 cc=pd.DataFrame(controls);cc.to_csv(OUT/'stage3_same_cell_program_coexpression_context.csv',index=False);pd.DataFrame(dist).to_csv(OUT/'stage3_depth_morphology_context_distributions.csv',index=False)
 for f in cc.columns:
  if f.startswith('n_') and f!='n_cells':cc[f]=cc[f].astype(np.int64)
 nums=[f for f in cc.columns if f.startswith('n_')];cc.groupby('native_stage3_label')[nums].sum().reset_index().to_csv(OUT/'stage3_same_cell_context_cohort_counts.csv',index=False)
 # Independently observed reference positive fractions reveal that a gene can
 # be absent from a repository's T core yet expressed in source T populations.
 ref=pd.read_parquet(RNA/'independent_reference_all_gene_donor_profiles.parquet');nkgenes=json.load(open('/data/xenium_typing/pln_20261003/release.json'))['panel']['nk_exclusive_naming_genes'];rr=ref.loc[ref.gene.isin(nkgenes)&ref.native_panel_measured].copy();rr['reference_broad_context']=rr.source_cell_type.map(maptype);rr.to_csv(OUT/'native_NK_exclusive_genes_independent_reference_expression.csv',index=False)
 cv=pd.read_csv(OUT/'SingleR_held_out_reference_profile_predictions.csv');cv['source_broad']=cv.true_source_type.map(maptype);cv['predicted_broad']=cv.source_type.map(maptype);cv['broad_context_matches']=cv.source_broad.eq(cv.predicted_broad);cv.groupby(['comparison','source_broad']).agg(n_donor_category_profiles=('key','size'),n_donors=('held_out_reference_donor','nunique'),mean_profile_context_match=('broad_context_matches','mean'),n_reference_cells_represented=('true_reference_cells','sum')).reset_index().to_csv(OUT/'reference_donor_holdout_context_summary.csv',index=False)
 report={'status':'all27_pLN_stage1_stage2_stage3_full_panel_descriptive_evidence_complete','n_sections':27,'n_QC_cells':int(meta.n_cells.sum()),'n_partitions':len(meta),'n_native_TNK_cells':int(ss.n_cells.sum()),'n_native_stage3_TNK_strata':len(ss),'n_stage1_nonimmune_reference_immune_context_partitions':int((non&a.full_reference_stage1_immune_context).sum()),'n_stage1_nonimmune_reference_immune_context_cells':int(a.loc[non&a.full_reference_stage1_immune_context,'n_cells'].sum()),'stage2_context_summary':a.groupby(['cluster_lineage','stage2_evidence_status']).n_cells.agg(['size','sum']).reset_index().to_dict('records'),'n_stage2_B_admitted_only_shared_ER':int(sum(r['only_shared_ER_naming_standouts'] for r in admissions if r['native_stage2_label']=='Immune_B_Plasma')),'native_labels_modified':False,'independent_same_cell_accuracy_established':False,'stage3_cell_flags':cokeys,'script_path':str(Path(__file__)),'script_sha256':sha(Path(__file__)),'outputs':[{'path':str(p),'sha256':sha(p)} for p in OUT.iterdir() if p.is_file() and p.name not in ['pLN_full_panel_evidence_manifest.json'] and p.suffix not in ['.log','.R','.py']]};(OUT/'pLN_full_panel_evidence_manifest.json').write_text(json.dumps(report,indent=2,allow_nan=False)+'\n');print('done',len(meta),len(ss),int(ss.n_cells.sum()),flush=True);print(pd.read_csv(OUT/'stage3_same_cell_context_cohort_counts.csv').to_string(index=False))
if __name__=='__main__':main()
